#function to calculate classical cross-correlations
def av_xcorr_td(filelist, max_lag):
    import numpy as np
    import obspy as obs
    import sys
    import os
    sys.path.append("./processing")
    sys.path.append("./general")
    from obspy.signal import cross_correlation
    import matplotlib.pyplot as plt
 
    #Open filelist
    filelist = open(filelist, 'r')
    filelist=filelist.read().split('\n')[:-1]
    
    #Cross-correlation loop over files
    for i in range(len(filelist)):
       for j in range(len(filelist)):
            if i<j: continue
            
            filename1=filelist[i]
            filename2=filelist[j]
            dat1=obs.read(filename1)[0]
            dat2=obs.read(filename2)[0]        
           
            if  dat1.stats.station!=dat2.stats.station and dat1.stats.channel==dat2.stats.channel:
                
                if dat1.stats.sampling_rate!=dat2.stats.sampling_rate:
                    print 'Unequal sampling rates. Cannot crosscorrelate.'
                    continue
                    
                else:
                    numsamples=int(max_lag*dat1.stats.sampling_rate)
                    #Create an option to specify maximum lag
                
                
                print '---------------------------------------------------------------------'
                print 'Crosscorrelation of...'
                print filename1 
                print filename2
                
                xcorr=cross_correlation.xcorr(dat1.data, dat2.data, numsamples,True)[2]
                            
                taxis=np.linspace(-max_lag, max_lag, len(xcorr))
                plt.plot(taxis,xcorr)
                plt.xlabel('Time (s)')
                plt.ylabel('Cross-correlation amplitude')
                plt.show()

