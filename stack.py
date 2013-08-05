import numpy as np
from obspy.core import read
from obspy.core import stream
from obspy.core import trace
from obspy.core import  UTCDateTime
import os
import shutil

import matplotlib.pyplot as plt
import TOOLS.read_xml as rxml
import TOOLS.correlations as corr


def stack(xmlinput):

    """
    Here comes the help and descrition
    """

    #==============================================================================================
    # Initialize
    #==============================================================================================

    #- Read the input from xml file----------------------------------------------------------------
    inp1=rxml.read_xml(xmlinput)
    inp1=inp1[1]

    indir=inp1['directories']['indir']
    outdir=inp1['directories']['outdir']
    if os.path.isdir(outdir)==False:
        os.mkdir(outdir)
        
    verbose=bool(int(inp1['verbose']))
    plotting=bool(int(inp1['plotting']))
    
    save_lin=inp1['stacks']['save_linstack']
    save_pws=inp1['stacks']['save_pwstack']
    pws_nu=int(inp1['stacks']['pws_nu'])
    
    startday=UTCDateTime(inp1['timethings']['startdate'])
    endday=UTCDateTime(inp1['timethings']['enddate'])
    win_len=int(inp1['timethings']['winlen'])
    olap=int(inp1['timethings']['olap'])

    channels=inp1['channels']['channel_list'].split(' ')
    mix_channels=bool(int(inp1['channels']['mix_channels']))
    
    if verbose:
        print 'startday = ',  startday
        print 'endday = ',  endday
        print 'Window length = ', win_len, ' s'
        print 'Overlap = ', olap, ' s'
        print 'Channels: ', channels
        print 'Mix channels:', mix_channels
    
    #- Input parameters for correlation: Correlation type, maximum lag, nu parameter
    corr_type=inp1['correlations']['corr_type']
    #- Force maximum lag to be integer. In seconds:
    maxlag=int(inp1['correlations']['max_lag'])
    pcc_nu=int(inp1['correlations']['pcc_nu'])
    
    #- Get rid of annoyances
    if os.path.exists(indir+'/.DS_Store'):
        os.remove(indir+'/.DS_Store')
    
    #- copy the input xml to the output directory for documentation
    shutil.copy(xmlinput,outdir)

    #- list available data files (should be one per channel)
    record_list=os.listdir(indir+'/')

    #- Find out what are the relevant combinations. This returns a list of tuples with identifiers that are to be correlated (e. g. ('G.ECH.00.BHE','G.CAN.00.BHE'))
    corr_ch=find_pairs(record_list,channels,mix_channels)

    if verbose:
        print 'number of potential correlations: '+str(len(corr_ch))

    #==============================================================================================    
    #-Loop over station pairs======================================================================
    #==============================================================================================

    for chpair in corr_ch:

        #==========================================================================================
        #- compute pairwise correlations
        #==========================================================================================

        if verbose:
            print '================================================================================'
            print 'Stacking correlations for:'
            print chpair[0]
            print chpair[1]
            
        #- Open files. Assumes one trace per file, which should be done by the preprocessing
        try:
            dat1=read(indir+'/'+chpair[0])[0]
            dat2=read(indir+'/'+chpair[1])[0]
        except (IOError,TypeError):
            if verbose: print 'One or both files not found. Skipping this correlation.'
            continue
     
        #- Test if sampling rates are same. Should have been enforced by the preprocessing.
        if dat1.stats.sampling_rate!=dat2.stats.sampling_rate:
            if verbose: print 'Unequal sampling rates. Skipping this correlation.'
            continue


        #- Compute stacked correlations ===========================================================
        (correlation_stack, coherence_stack, windows, n, n_skip)=stack_windows(dat1, dat2, startday, endday, win_len, olap, corr_type, maxlag, pcc_nu, verbose)
        if verbose: 
            print 'Number of successfully stacked time windows: ', n
            print 'Number of skipped time windows: ', n_skip
        
        #==========================================================================================
        #- produce written and visual output provided that at least one time window could be used
        #==========================================================================================    

        if n>0:

            #- plot correlation function, if wanted ===============================================

            if plotting:
                t=np.linspace(-maxlag*dat1.stats.sampling_rate,maxlag*dat1.stats.sampling_rate, len(correlation_stack))
                plt.subplot(311)
                plt.plot(t,correlation_stack)
                plt.ylabel('linear stack')
                plt.subplot(312)
                plt.plot(t,np.abs(coherence_stack))
                plt.ylabel('phase coherence')
                plt.subplot(313)
                plt.plot(t,correlation_stack*np.abs(coherence_stack))
                plt.xlabel('t [s]')
                plt.ylabel('phase-weighted stack')
                plt.show()
            
            #- Write correlation function to a file ===============================================
        
            #- Create a trace object and fill in the basic information
            tr1=trace.Trace()
            tr1.stats.sampling_rate=dat1.stats.sampling_rate
            tr1.data=correlation_stack

            #- fill in all the necessary information
        
            #- open file and write correlation function
            fn1=dat1.stats.network+'.'+dat1.stats.station+'.'+dat1.stats.location+'.'+dat1.stats.channel
            fn2=dat2.stats.network+'.'+dat2.stats.station+'.'+dat2.stats.location+'.'+dat2.stats.channel
            fileid=outdir+'/'+fn1+'-'+fn2+'.lin_stack.MSEED'

            #- if file already exists, open the old one and add traces together
            if os.path.exists(fileid)==True:
                if verbose==True: "Correlation function already exists. Add to previous one."
                tr_old=read(fileid)
                tr1.data+=tr_old[0].data

            tr1.write(fileid, format="MSEED")

            #- Write time windows to a file for documentation =====================================
            filename=outdir+'/'+fn1+'-'+fn2+'.windows'
            fid=open(filename,'a')
            for window in windows:
               fid.write(str(window)+'\n')
            fid.close()

        #Create a trace object
#        tr2=trace.Trace()
#        tr2.stats.sampling_rate=Fs
#        tr2.data=xcorrstack
        
#        fileid=outdir+'/'+chpair[0]+'.'+chpair[1]+'.pw_stack'
        #append start and end date to fileid?
#        tr2.write(fileid, format="MSEED")
        


#==================================================================================================
# find pairs of recordings
#==================================================================================================

def find_pairs(record_list,channels,mix_channels):
    
    """
    Find pairs of recordings.

    ccpairs=find_pairs(record_list,channels,mix_channels):

    record_list:    list of seismogram files following the naming convention network.station.location.channel
    channels:       list of channels to be considered, e.g. ['BHZ','LHE']
    mix_channels:   boolean parameter determining if pairs are allowed to have different channels

    ccpairs:        list of seismogram pairs

    """

    ccpairs=[]
    
    for i in range(len(record_list)):
        for j in range(len(record_list)):
            if i<j: continue

            #- stations
            sta1=record_list[i].split('.')[1]
            sta2=record_list[j].split('.')[1]

            #- channels
            cha1=record_list[i].split('.')[3]
            cha2=record_list[j].split('.')[3]

            #- times
            t11=UTCDateTime(record_list[i].split('.')[4])
            t12=UTCDateTime(record_list[i].split('.')[5])
            t21=UTCDateTime(record_list[j].split('.')[4])
            t22=UTCDateTime(record_list[j].split('.')[5])

            chcomb=()
            make_pair=False

            if (cha1 in channels) & (cha2 in channels):
                #- if channels differ, make pair only when channel mixing is allowed
                if (cha1!=cha2 and mix_channels):
                    if (t12>t21) and (t11<t12) and (t21<t22):
                        make_pair=True
                #- make a pair when channels are identical
                elif (cha1==cha2):
                    if (t12>t21) and (t11<t12) and (t21<t22):
                        make_pair=True

            if make_pair:
                chcomb=(record_list[i],record_list[j])
                ccpairs.append(chcomb)
                    
    return ccpairs
    
    
#==================================================================================================
# Compute stecked correlation functions
#==================================================================================================

def stack_windows(dat1, dat2, startday, endday, win_len, olap, corr_type, maxlag, pcc_nu, verbose):

    """
    Compute stacked correlation functions.

    correlation_stack,coherence_stack,windows,n,n_skip=stack_windows(dat1, dat2, startday, endday, win_len, olap, corr_type, maxlag, pcc_nu, phase, verbose):

    dat1:       first time series
    dat2:       second time series
    startday:   UTC starting time of the first correlation window
    win_len:    length of the time windows to be correlated
    olap:       overlap of time windows
    endday:     UTC time when the last correlation window starts
    corr_type:  type of correlation functions
    maxlag:     maximum time lag in the correlation functions
    pcc_nu:     exponent in the phase cross-correlation
    verbose:    talk or not

    correlation_stack:  stacked correlations
    coherence_stack:    stacked phase coherences
    windows:            time windows used in the stack
    n:                  number of successfully stacked correlations
    n_skip:             number of discarded time windows

    """

    from scipy.signal import hilbert

    #- initialisations ----------------------------------------------------------------------------

    #- Prepare stack
    #stack=np.zeros((int(2.0*float(maxlag)*dat1.stats.sampling_rate+1.0), ))
       
    #- Counter for number of successfully correlated time windows
    n=0

    #- Counter for failed time windows
    n_skip=0

    #- Initial time window and initial time window pairs for documentation
    t1=startday
    t2=t1+win_len
    windows=[]

    if corr_type not in ['td_classic','fd_classic','pcc']:
        if verbose: 'Correlation type '+corr_type+' not supported'
        return([], [], [], 0, 0)

    #- Loop over time windows and update stack ----------------------------------------------------
    while t2<=endday:

        #- Get the portion of the trace that you want
        tr1=dat1.copy()
        tr2=dat2.copy()
            
        tr1.trim(starttime=t1, endtime=t2)
        tr2.trim(starttime=t1, endtime=t2)
    
        #- Perform a series of check on the time series
        correlate=True

        if len(tr1.data)!=len(tr2.data):
            #if verbose: print "Traces of unequal length (%d, %d) in time window %s to %s, skipped" % (len(tr1.data),len(tr2.data),str(t1),str(t2))
            correlate=False
        if len(tr1.data)==0:
            #if verbose: print "No data for station %s in time window %s to %s, skipped" % (tr1.stats.station,str(t1),str(t2))
            correlate=False
        if len(tr2.data)==0:
            #if verbose: print "No data for station %s in time window %s to %s, skipped" % (tr2.stats.station,str(t1),str(t2))
            correlate=False
        if (True in np.isnan(tr1.data)) or (True in np.isnan(tr2.data)):
            #if verbose: print 'Traces contain NaN in time window '+str(t1)+' to '+str(t2)+', skipped'
            correlate=False
        if (True in np.isinf(tr1.data)) or (True in np.isinf(tr2.data)):
            #if verbose: print 'Traces contain Inf in time window '+str(t1)+' to '+str(t2)+', skipped'
            correlate=False

        #- Compute correlations, provided that time series are okay
        if correlate==False:
            n_skip+=1
        else:
            if corr_type=='td_classic':
                correlation=corr.xcorrelation_td(tr1, tr2, maxlag)
            elif corr_type=='fd_classic':
                correlation=corr.xcorrelation_fd(tr1, tr2)
            elif corr_type=='pcc':
                correlation=corr.phase_xcorrelation(tr1, tr2, maxlag, pcc_nu)[0]

            #- update statistics
            n+=1

            #- linear stack =======================================================================
            if n==1:
                correlation_stack=correlation
            else:
                correlation_stack+=correlation

            #- phase coherence stack ==============================================================
            coherence=hilbert(correlation)
            tol=np.mean(np.abs(coherence))/10.0
            coherence=coherence/(np.abs(coherence)+tol)
            if n==1:
                coherence_stack=coherence
            else:
                coherence_stack+=coherence

            #- make time window pairs for documentation
            window=(t1,t2)
            windows.append(window)

            #- phase-weighted stack
            # *** to be implemented ***
            #- At this point one could simply compute the phase weighting function, based on "stack".

        #- go to the next time window
        t1=t2-olap
        t2=t1+win_len
        
    if n==0:
        return([], [], [], n, n_skip)
    else:
        return(correlation_stack, coherence_stack, windows, n, n_skip)
