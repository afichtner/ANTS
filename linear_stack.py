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
def stack(stack_input):
    #Read the input
    inp1=rxml.read_xml(stack_input)
    inp1=inp1[1]

    indir=inp1['directories']['indir']
    outdir=inp1['directories']['outdir']
    if os.path.isdir(outdir)==False:
        os.mkdir(outdir)
        
    verbose=bool(int(inp1['verbose']))
    plotting=bool(int(inp1['plotting']))
    
    
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
    nu=int(inp1['correlations']['pcc_nu'])
    
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
            
        #test if the starting times are the same, otherwise apply a correction
        st1=dat1.stats.starttime
        st2=dat2.stats.starttime

        if st1!=st2:
            offset=(st1-st2)*int(Fs)
            if verbose: print 'Warning: Unequal start times.\nOffset of', offset,'samples'
    
    
#-Get the stacked correlation========================================================================

        (xcorrstack, n, n_skip)=stack_windows(dat1, dat2, startday, win_len, olap, endday, corr_type, maxlag, nu, verbose)
        if verbose: 
            print 'Number of successfully stacked time windows: ', n
            print 'Number of skipped time windows: ', n_skip
        
        
        if plotting:
            x=np.linspace(-maxlag*dat1.stats.sampling_rate,maxlag*dat1.stats.sampling_rate, len(xcorrstack))
            plt.plot(x, xcorrstack, linewidth=1.8)
            plt.show()
            
#-Write it to a file=================================================================================
        #Create a trace object
        tr=trace.Trace()
        tr.stats.sampling_rate=Fs
        tr.data=xcorrstack
        
        fileid=outdir+'/'+chpair[0]+'.'+chpair[1]
        #append start and end date to fileid?
        tr.write(fileid, format="MSEED")
        






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
    
    
def stack_windows(dat1, dat2, startday, win_len, olap, endday, corr_type, maxlag, nu, verbose):
        #Prepare stacking loop
        xcorrstack=np.zeros((2*maxlag*int(dat1.stats.sampling_rate)+1, ))
        
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
               
                if verbose: print 'Data gap found in this time window, skipping correlation.'
                t1=t2-olap
                t2=t1+win_len-1
                n_skip+=1
                continue 
            
            #Obtain correlation
            if corr_type=='1':
                xcorrstack+=corr.xcorrelation_td(tr1, tr2, maxlag)
            elif corr_type=='2':
                xcorrstack+=corr.xcorrelation_fd(tr1, tr2)
            elif corr_type=='3':
                xcorrstack+=corr.phase_xcorrelation(tr1, tr2, maxlag, nu)[0]
            else: 
                if verbose: print 'Invalid Correlation type in input file.'
                return()
                #how to properly deal with this error so that the program does not crash?
            n+=1
            
            #update time window
            t1=t2-olap
            t2=t1+win_len-1
            
        #Linear stack is easy, what about phase weighted stack? Complicated....
        xcorrstack/=n
        
        
        return(xcorrstack, n, n_skip)
        
