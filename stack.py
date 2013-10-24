import numpy as np
from obspy.core import read
from obspy.core import stream
from obspy.core import trace
from obspy.core import  UTCDateTime
import os
import shutil
import re
import matplotlib.pyplot as plt
import TOOLS.read_xml as rxml
import TOOLS.correlations as corr
from math import log
from glob import glob

def stack(xmlinput):

    """
    A routine to calculate and stack cross-correlations/phase cross-correlations.
    xmlinput: Path to file that contains all input parameters.
    """

    #==============================================================================================
    # Initialize
    #==============================================================================================

    #- Read the input from xml file----------------------------------------------------------------
    inp1=rxml.read_xml(xmlinput)
    inp1=inp1[1]

    indir=inp1['directories']['indir']
    outdir=inp1['directories']['outdir']
    fformat=inp1['directories']['format']
    if os.path.isdir(outdir)==False:
        os.mkdir(outdir)
        
    verbose=bool(int(inp1['verbose']))
    plotting=bool(int(inp1['plotting']))
    
    
    startday=UTCDateTime(inp1['timethings']['startdate'])
    endday=UTCDateTime(inp1['timethings']['enddate'])
    win_len=int(inp1['timethings']['winlen'])
    Fs=float(inp1['timethings']['Fs'])
    olap=int(inp1['timethings']['olap'])
    
    win_len=wl_adjust(win_len, Fs, verbose)
    
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
    
    
    #- copy the input xml to the output directory for documentation
    shutil.copy(xmlinput,outdir)
    #copy also the preprocessing input to the output directory
    #This overwrites the first file if they are called the same! Dont call them the same!
    xmlfile=glob(indir+"/*.xml")
    shutil.copy(xmlfile, outdir)

    #- list of available data files (should be one per channel)
    record_list=os.listdir(indir+'/')

    #- Find out what are the relevant combinations. This returns a list of tuples with identifiers that are to be correlated (e. g. ('G.ECH.00.BHE','G.CAN.00.BHE'))
    corr_ch=find_pairs(record_list,fformat, channels,mix_channels,win_len, verbose)

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
     
        #- Test if sampling rates are the same. Should have been enforced by the preprocessing.
        if dat1.stats.sampling_rate!=dat2.stats.sampling_rate:
            if verbose: print 'Unequal sampling rates. Skipping this correlation.'
            continue
            
        #- Compute stacked correlations ===========================================================
        (correlation_stack, coherence_stack, windows, n, n_skip)=stack_windows(dat1, dat2, startday, endday, win_len, olap, corr_type, maxlag, pcc_nu,  verbose)
        if verbose: 
            print 'Number of successfully stacked time windows: ', n
            print 'Number of skipped time windows: ', n_skip

        #==========================================================================================
        #- write metadata
        #==========================================================================================

        fn1=dat1.stats.network+'.'+dat1.stats.station+'.'+dat1.stats.location+'.'+dat1.stats.channel
        fn2=dat2.stats.network+'.'+dat2.stats.station+'.'+dat2.stats.location+'.'+dat2.stats.channel
        mdfilename=outdir+'/'+fn1+'-'+fn2+corrtype+'.metadata'
        
        if os.path.exists(mdfilename)==False:
            fid=open(mdfilename,'w')
            fid.write('= station information ================================\n')
            fid.write('network_1: '+dat1.stats.network+'\n')
            fid.write('station_1: '+dat1.stats.station+'\n')
            fid.write('location_1: '+dat1.stats.location+'\n')
            fid.write('channel_1: '+dat1.stats.channel+'\n')
            fid.write('network_2: '+dat2.stats.network+'\n')
            fid.write('station_2: '+dat2.stats.station+'\n')
            fid.write('location_2: '+dat2.stats.location+'\n')
            fid.write('channel_2: '+dat2.stats.channel+'\n')
            fid.write('= contributing time windows ==========================\n')
            fid.close()

        
        #==========================================================================================
        #- produce written and visual output provided that at least one time window could be used
        #==========================================================================================    

        if n>0:

            #- plot correlation function, if wanted ===============================================

            if plotting:
                t=np.linspace(-maxlag*dat1.stats.sampling_rate,maxlag*dat1.stats.sampling_rate, len(correlation_stack))
                plt.subplot(311)
                plt.plot(correlation_stack)
                plt.ylabel('linear stack')
                plt.subplot(312)
                plt.plot(np.abs(coherence_stack))
                plt.ylabel('phase coherence')
                plt.subplot(313)
                plt.plot(correlation_stack*np.abs(coherence_stack))
                plt.xlabel('t [s]')
                plt.ylabel('phase-weighted stack')
                plt.show()
            
            #- Write correlation function to a file ===============================================
        
            #- Create a trace object and fill in the basic information
            tr_correlation_stack=trace.Trace()
            tr_coherence_stack_real=trace.Trace()
            tr_coherence_stack_imag=trace.Trace()

            tr_correlation_stack.stats.sampling_rate=dat1.stats.sampling_rate
            tr_correlation_stack.stats.starttime=UTCDateTime(2000,1,1,0,0)-dat1.stats.delta*float((len(correlation_stack)-1))/2.0
            tr_correlation_stack.data=correlation_stack

            tr_coherence_stack_real.stats.sampling_rate=dat1.stats.sampling_rate
            tr_coherence_stack_real.stats.starttime=UTCDateTime(2000,1,1,0,0)-dat1.stats.delta*float((len(correlation_stack)-1))/2.0
            tr_coherence_stack_real.data=np.real(coherence_stack)

            tr_coherence_stack_imag.stats.sampling_rate=dat1.stats.sampling_rate
            tr_coherence_stack_imag.stats.starttime=UTCDateTime(2000,1,1,0,0)-dat1.stats.delta*float((len(correlation_stack)-1))/2.0
            tr_coherence_stack_imag.data=np.imag(coherence_stack.imag)


            #- open file and write correlation function
            fileid_correlation_stack=outdir+'/'+fn1+'-'+fn2+'.'+corr_type+'_stack.MSEED'
            fileid_coherence_stack_real=outdir+'/'+fn1+'-'+fn2+'.coherence_stack_real.MSEED'
            fileid_coherence_stack_imag=outdir+'/'+fn1+'-'+fn2+'.coherence_stack_imag.MSEED'

            #- linear stack
            if os.path.exists(fileid_correlation_stack)==True:
                if verbose: 
                    print "Correlation stack already exists. Add to previous one."
                tr_old=read(fileid_correlation_stack)
                tr_correlation_stack.data=tr_correlation_stack.data+tr_old[0].data
                tr_correlation_stack.write(fileid_correlation_stack, format="MSEED")
            else:
                tr_correlation_stack.write(fileid_correlation_stack, format="MSEED")

            #- real part of coherence stack
            if os.path.exists(fileid_coherence_stack_real)==True:
                if verbose: 
                    print "Real part of coherence stack already exists. Add to previous one."
                tr_old=read(fileid_coherence_stack_real)
                tr_coherence_stack_real.data=tr_coherence_stack_real.data+tr_old[0].data
                tr_coherence_stack_real.write(fileid_coherence_stack_real, format="MSEED")
            else:
                tr_coherence_stack_real.write(fileid_coherence_stack_real, format="MSEED")

            #- imaginary part of coherence stack
            if os.path.exists(fileid_coherence_stack_imag)==True:
                if verbose:
                    print "Imaginary part of coherence stack already exists. Add to previous one."
                tr_old=read(fileid_coherence_stack_imag)
                tr_coherence_stack_imag.data=tr_coherence_stack_imag.data+tr_old[0].data
                tr_coherence_stack_imag.write(fileid_coherence_stack_imag, format="MSEED")
            else:
                tr_coherence_stack_imag.write(fileid_coherence_stack_imag, format="MSEED")

            #- Write time windows to a file for documentation =====================================
            #filename=outdir+'/'+fn1+'-'+fn2+'.'+corr_type+'.metadata'
            fid=open(mdfilename,'a')
            for window in windows:
               fid.write(str(window[0].year)+' '+str(window[0].month)+' '+str(window[0].day)+' '+str(window[0].hour)+' '+str(window[0].minute)+' '+str(window[0].second)+', '+str(window[1].year)+' '+str(window[1].month)+' '+str(window[1].day)+' '+str(window[1].hour)+' '+str(window[1].minute)+' '+str(window[1].second)+'\n')
            fid.close()


#==================================================================================================
# find pairs of recordings
#==================================================================================================

def find_pairs(record_list,format,channels,mix_channels,win_len,verbose):
    
    """
    Find pairs of recordings with overlapping time windows.

    ccpairs=find_pairs(record_list,channels,mix_channels):

    record_list:    list of seismogram files following the naming convention network.station.location.channel
    format:         string specifying the file format, mseed/sac
    channels:       list of channels to be considered, e.g. ['BHZ','LHE']
    mix_channels:   boolean parameter determining if pairs are allowed to have different channels
    win_len:        length of the time windows to be correlated in seconds, used to check minimum length of traces

    ccpairs:        list of seismogram pairs

    """

    ccpairs=[]
    format=re.compile(format, re.I)
    
    for i in range(len(record_list)):
        for j in range(len(record_list)):
            if i<j: continue
            if format.search(record_list[i]) is None: continue

            #- get station and channel names, as well as start and end times ----------------------

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

            #- perform tests to compile the station pair list -------------------------------------

            if (cha1 in channels) & (cha2 in channels) & (t12>=t11+win_len) & (t22>=t21+win_len):

                #- if channels differ, make pair only when channel mixing is allowed
                if (cha1!=cha2 and mix_channels):
                    if (t12>t21) and (t11<t12):
                        make_pair=True
                #- make a pair when channels are identical
                elif (cha1==cha2):
                    if (t12>t21) and (t11<t12):
                        make_pair=True

            #- add to the list of pairs -----------------------------------------------------------
    
            if make_pair:
                chcomb=(record_list[i],record_list[j])
                ccpairs.append(chcomb)
                    
    return ccpairs
    
    
#==================================================================================================
# Compute stacked correlation functions
#==================================================================================================

def stack_windows(dat1, dat2, startday, endday, win_len, olap, corr_type, maxlag, pcc_nu,verbose):

    """
    Compute stacked correlation functions.

    correlation_stack,coherence_stack,windows,n,n_skip=stack_windows(dat1, dat2, startday, endday, win_len, olap, corr_type, maxlag, pcc_nu, phase, verbose):

    dat1:       first time series
    dat2:       second time series
    startday:   UTC starting time of the first correlation window
    win_len:    length of the time windows to be correlated
    olap:       overlap of time windows
    endday:     UTC time when the last correlation window starts
    corr_type:  type of correlation functions, ccc or pcc
    maxlag:     maximum time lag in the correlation functions
    pcc_nu:     exponent in the phase cross-correlation
    verbose:    talk or not

    correlation_stack:  stacked correlations
    coherence_stack:    stacked phase coherences, this is the complex coherence before taking the absolute value
    windows:            time windows used in the stack
    n:                  number of successfully stacked correlations
    n_skip:             number of discarded time windows

    """
   
    from scipy.signal import hilbert

    #- initialisations ----------------------------------------------------------------------------
    if verbose: print 'Hi from stacking routine.'
   
    #- Get sampling frequency
    Fs=dat1.stats.sampling_rate
    
    #- Counter for number of successfully correlated time windows
    n=0

    #- Counter for failed time windows
    n_skip=0

    #- check how far the traces go!
    startday=max(startday, dat1.stats.starttime, dat2.stats.starttime)
    endday=min(endday, dat1.stats.endtime, dat2.stats.endtime)
    
    #- Initial time window and initial time window pairs for documentation
    t1=startday
    t2=t1+win_len
   
    windows=[]

    if corr_type not in ['ccc', 'fcc','pcc']:
        if verbose: 'Correlation type '+corr_type+' not supported'
        return([], [], [], 0, 0)
        
    

    #- Loop over time windows and update stack ----------------------------------------------------
    while t2<=endday:
        
        correlate=True

        #- Get the portion of the trace that you want
        tr1=dat1.copy()
        tr2=dat2.copy()
        
        if (dat1.stats.starttime<=t1) and (dat1.stats.endtime>=t2) and (dat2.stats.starttime<=t1) and (dat2.stats.endtime>=t2):
            tr1.trim(starttime=t1, endtime=t2-1/Fs)
            tr2.trim(starttime=t1, endtime=t2-1/Fs)
         
            #- Perform a series of checks on the time series
            if len(tr1.data)!=len(tr2.data):
                if verbose: print "Traces of unequal length (%d, %d) in time window %s to %s, skipped" % (len(tr1.data),len(tr2.data),str(t1),str(t2))
                correlate=False
            if len(tr1.data)==0:
                if verbose: print "No data for station %s in time window %s to %s, skipped" % (tr1.stats.station,str(t1),str(t2))
                correlate=False
            if len(tr2.data)==0:
                if verbose: print "No data for station %s in time window %s to %s, skipped" % (tr2.stats.station,str(t1),str(t2))
                correlate=False
            if (True in np.isnan(tr1.data)) or (True in np.isnan(tr2.data)):
                if verbose: print 'Traces contain NaN in time window '+str(t1)+' to '+str(t2)+', skipped'
                correlate=False
            if (True in np.isinf(tr1.data)) or (True in np.isinf(tr2.data)):
                if verbose: print 'Traces contain Inf in time window '+str(t1)+' to '+str(t2)+', skipped'
                correlate=False
        
            
            
        else:
            #if verbose: print "Traces do not cover the desired window. T1 is %s whereas windows start at %s  %s\n T2 is %s whereas windows end at %s  %s." % (str(t1), str(dat1.stats.starttime), str(dat2.stats.starttime), str(t2), str(dat1.stats.endtime),str(dat2.stats.endtime))
            correlate=False
       
       
        
        #- Compute correlations, provided that time series are okay
        if correlate==False:
            n_skip+=1
        else:
            if corr_type=='ccc':
                correlation=corr.xcorrelation_td(tr1, tr2, maxlag)
            elif corr_type=='fcc':
                correlation=corr.xcorrelation_fd(tr1, tr2)
            elif corr_type=='pcc':
                correlation=corr.phase_xcorrelation(tr1, tr2, maxlag, pcc_nu)
        
            #- update statistics
            n+=1

            #- linear stack =======================================================================
            if n==1:
                correlation_stack=correlation
            elif len(correlation_stack)==len(correlation):
                correlation_stack+=correlation

            #- phase coherence stack ==============================================================
            coherence=hilbert(correlation)
            tol=np.mean(np.abs(coherence))/10.0
            coherence=coherence/(np.abs(coherence)+tol)
            if n==1:
                coherence_stack=coherence
            elif len(coherence_stack)==len(coherence):
                coherence_stack+=coherence

            #- make time window pairs for documentation
            window=(t1,t2)
            windows.append(window)
           

        #- go to the next time window
        t1=t2-olap
        t2=t1+win_len
        
    if n==0:
        return([], [], [], n, n_skip)
    else:
        return(correlation_stack, coherence_stack, windows, n, n_skip)
        
        
#==================================================================================================
# Adjust window lenght to be power of 2
#==================================================================================================
"""Get the nearest power of two in terms of samples; then determine the corresponding window length in seconds. 
win_len: Integer, User-defined window length in seconds
Fs: Integer, Sampling rate """
from math import ceil,  log

def wl_adjust(win_len, Fs, verbose):
    
    #current window length
    cwl=Fs*win_len;
    nwl=int((2**ceil(log(cwl)/log(2)))/Fs)
    
    if cwl!=nwl:
        if verbose:
            print 'Window length adjusted to '+str(nwl)
    return nwl
