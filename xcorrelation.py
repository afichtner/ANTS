import numpy as np
from obspy.core import read
from obspy.core import stream
import os

from obspy.signal import cross_correlation
import matplotlib.pyplot as plt
import TOOLS.read_xml as rxml
import TOOLS.correlations as corr


#==================================================================================================
# Run different types of correlation on a number of files
#==================================================================================================
def av_xcorr(corr_input):
    
    #read input file
    inp1=rxml.read_xml(corr_input)
    inp1=inp1[1]

    indir=inp1['directories']['indir']
    outdir=inp1['directories']['outdir']
    verbose=bool(int(inp1['verbose']))
    
    
    #Create filelist
    cmd='ls '+indir+' > filelist.txt'
    os.system(cmd)
    #Open filelist
    filelist = open('filelist.txt', 'r')
    filelist=filelist.read().split('\n')[:-1]
    
    #Determine the maximum lag of correlations
    max_lag=float(inp1['correlations']['max_lag'])
    
    #Cross-correlation loop over files
    for i in range(len(filelist)):
       for j in range(len(filelist)):
            if i<j: continue
            
            filename1=indir+'/'+filelist[i]
            filename2=indir+'/'+filelist[j]
            dat1=read(filename1)[0]
            dat2=read(filename2)[0]        
           
            if  dat1.stats.station!=dat2.stats.station and dat1.stats.channel==dat2.stats.channel:
                
                if dat1.stats.sampling_rate!=dat2.stats.sampling_rate:
                    print 'Unequal sampling rates. Cannot crosscorrelate.'
                    continue
        
                if verbose:
                    print '---------------------------------------------------------------------'
                    print 'Crosscorrelation of...'
                    print filename1 
                    print filename2
                
                #- Classical cross-correlation time domain============================================================================
                if inp1['correlations']['cc_td']=='1':
                    if verbose:
                        print '* Time-domain classical cross-correlation'
                    cc_td=corr.xcorrelation_td(dat1, dat2, max_lag)
                    taxis=np.linspace(-max_lag, max_lag, len(cc_td), True)
                    plt.plot(taxis, cc_td, linewidth=2.0)
                
                
                #- Classical cross-correlation frequency domain=======================================================================
                if inp1['correlations']['cc_fd']=='1':
                    if verbose:
                        print '* Frequency-domain classical cross-correlation'
                    cc_fd=corr.xcorrelation_fd(dat1, dat2)
                    taxis=np.linspace(0,len(cc_fd)/dat1.stats.sampling_rate, len(cc_fd), True)
                    plt.plot(taxis, cc_fd, linewidth=2.0)
                    
                #- Classical cross-correlation frequency domain=======================================================================
                if inp1['correlations']['pcc']['doit']=='1':
                    if verbose:
                        print '* Phase cross-correlation'
                    nu=int(inp1['correlations']['pcc']['nu'])
                    pcc=corr.phase_xcorrelation(dat1, dat2, max_lag, nu)
                    taxis=np.linspace(-max_lag, max_lag, len(pcc), True)
                    plt.plot(taxis, pcc, linewidth=2.0)
                plt.xlabel('Lage time seconds')
                plt.ylabel('Cross-Correlation')
                plt.legend(('Time domain CC', 'Freq. domain CC', 'Phase CC'))
                plt.show()

    #Clean up
    os.system('rm filelist.txt')
