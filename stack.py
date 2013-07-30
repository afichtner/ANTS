import numpy as np
from obspy.core import read
from obspy.core import stream
from obspy.core import trace
from obspy.core import  UTCDateTime
import os

import matplotlib.pyplot as plt
import TOOLS.read_xml as rxml
import TOOLS.correlations as corr

#==================================================================================================
# Initialize
#==================================================================================================
def stack(xmlinput):
    #Read the input
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
    
    if verbose:
        print 'startday = ',  startday
        print 'endday = ',  endday
        print 'Window length = ', win_len, ' s'
        print 'Overlap = ', olap, ' s'
    
    #Input parameters for correlation: Correlation type, maximum lag, nu parameter
    corr_type=inp1['correlations']['corr_type']
    #Force maximum lag to be integer. In seconds:
    maxlag=int(inp1['correlations']['max_lag'])
    pcc_nu=int(inp1['correlations']['pcc_nu'])
    
    #Get rid of annoyances
    if os.path.exists(indir+'/.DS_Store'):
        os.remove(indir+'/.DS_Store')
    
    #list available data files (should be one per channel)
    listch=os.listdir(indir+'/')
    #Find out what are the relevant channel combinations. This returns a list of tuples with identifiers that are to be correlated (e. g. ('G.ECH.00.BHE','G.CAN.00.BHE'))
    corr_ch=matchchannels2(listch)
    
#==================================================================================================    
#-Loop over station pairs============================================================================
#==================================================================================================

    for chpair in corr_ch:
        if verbose:
            print '================================================================================'
            print 'Stacking correlations for:'
            print chpair[0]
            print chpair[1]
            
        #open files, get the sampling rate
        try:
            dat1=read(indir+'/'+chpair[0]+'*')[0]
            dat2=read(indir+'/'+chpair[1]+'*')[0]
        except IOError:
            if verbose: print 'One or both files not found. Skipping this correlation.'
            continue
        
        Fs=dat1.stats.sampling_rate
        
        #test if sampling rates are same
        if Fs!=dat2.stats.sampling_rate:
            if verbose: print 'Unequal sampling rates. Skipping this correlation.'
            continue
    
#-Step 1: Get the phase stack=======================================================================
        phase=True
        (cps, n, n_skip)=stack_windows(dat1, dat2, startday, win_len, olap, endday, corr_type, maxlag, pcc_nu,  phase, verbose)
        
#-Get the stacked correlation========================================================================
        phase=False
        (xcorrstack, n, n_skip)=stack_windows(dat1, dat2, startday, win_len, olap, endday, corr_type, maxlag, pcc_nu, phase, verbose)
        if verbose: 
            print 'Number of successfully stacked time windows: ', n
            print 'Number of skipped time windows: ', n_skip
        
        pwstack=xcorrstack*abs(cps)**pws_nu
        
        if plotting:
            x=np.linspace(-maxlag*dat1.stats.sampling_rate,maxlag*dat1.stats.sampling_rate, len(xcorrstack))
            plt.figure(1)
            plt.subplot(211)
            plt.plot(x, xcorrstack,linewidth=1)
            plt.plot(x, pwstack, linewidth=1.8)
            plt.subplot(212)
            plt.plot(x, abs(cps))
            plt.show()
            
#-Write it to a file=================================================================================
        #Create a trace object
        tr1=trace.Trace()
        tr1.stats.sampling_rate=Fs
        tr1.data=xcorrstack
        
        fileid=outdir+'/'+chpair[0]+'.'+chpair[1]+'.lin_stack'
        #append start and end date to fileid?
        tr1.write(fileid, format="MSEED")
        
        #Create a trace object
        tr2=trace.Trace()
        tr2.stats.sampling_rate=Fs
        tr2.data=xcorrstack
        
        fileid=outdir+'/'+chpair[0]+'.'+chpair[1]+'.pw_stack'
        #append start and end date to fileid?
        tr2.write(fileid, format="MSEED")
        






def matchchannels2(channels):
    
    ccpairs=[]
    
    for i in range(len(channels)):
        for j in range(len(channels)):
            if i<j: continue
            inf1=channels[i].split('.')
            inf2=channels[j].split('.')
    
            #station
            sta1=inf1[1]
            sta2=inf2[1]
            #channel
            cha1=inf1[3]
            cha2=inf2[3]
          
            chcomb=()
            if  sta1!=sta2 and cha1==cha2:
                name1=inf1[0]+'.'+inf1[1]+'.'+inf1[2]+'.'+inf1[3]
                name2=inf2[0]+'.'+inf2[1]+'.'+inf2[2]+'.'+inf2[3]
                chcomb=(name1,name2)
                ccpairs.append(chcomb)
    
    return ccpairs
    
    
def stack_windows(dat1, dat2, startday, win_len, olap, endday, corr_type, maxlag, pcc_nu, phase, verbose):
        #Prepare stack
        if phase:
            stack=np.zeros((2*maxlag*int(dat1.stats.sampling_rate)+1, ), dtype=complex)
        else:
            stack=np.zeros((2*maxlag*int(dat1.stats.sampling_rate)+1, ))
            
        
        #Counter for days
        n=0
        #Counter for failed time windows
        n_skip=0
        #Time window bounds
        t1=startday
        t2=startday+win_len-1
        
#-Loop over days, update stack======================================================================
        while t2<=endday:
            
            # Get the portion of the trace that you want
            tr1=dat1.copy()
            tr2=dat2.copy()
            
            tr1.trim(starttime=t1, endtime=t2)
            tr2.trim(starttime=t1, endtime=t2)
            
            #At this point, check for nodata (we dont want to correlate data gaps)
            if len(tr1.data)!=len(tr2.data) or len(tr1.data)==0:
                t1=t2-olap
                t2=t1+win_len-1
                n_skip+=1
                continue 
            
            #Obtain correlation
            if corr_type=='1':
                stack+=corr.xcorrelation_td(tr1, tr2, maxlag, phase)
            elif corr_type=='2':
                stack+=corr.xcorrelation_fd(tr1, tr2, phase)
            elif corr_type=='3':
                stack+=corr.phase_xcorrelation(tr1, tr2, maxlag, pcc_nu, phase)[0]
            else: 
                if verbose: print 'Invalid Correlation type in input file.'
                return()
                #how to properly deal with this error so that the program does not crash?
            n+=1
            
            #update time window
            t1=t2-olap
            t2=t1+win_len-1
            
        #Linear stack is easy, what about phase weighted stack? Complicated....
        stack/=n
        
        
        return(stack, n, n_skip)
        
