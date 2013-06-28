#script to calculate cross-correlations
import numpy as np
import obspy as obs
import sys
import os
sys.path.append("./processing")
sys.path.append("./general")
from fourier_tools import xcorr_fd
import matplotlib.pyplot as plt

#Paths to files with time series to be crosscorrelated should be listed in this file:
if len(sys.argv)<2: print 'No channel list specified. Exiting.'; exit()

print 'Channel list to be cross-correlated:'
print sys.argv[1]


#Open filelist
filelist = open(sys.argv[1], 'r')
filelist = filelist.readlines()

endl=os.linesep


#Cross-correlation loop over files
for i in range(len(filelist)):
   for j in range(len(filelist)):
        
        if i<j: continue
        
        filename1=filelist[i].strip(endl)
        filename2=filelist[j].strip(endl)
        dat1=obs.read(filename1)[0]
        dat2=obs.read(filename2)[0]
       
        if  dat1.stats.station!=dat2.stats.station and dat1.stats.channel==dat2.stats.channel:
            print '---------------------------------------------------------------------'
            print 'Crosscorrelation of...'
            print filename1 
            print filename2
            
            dt=dat1.stats.delta
            xcorr=xcorr_fd(dat1.data, dat2.data, dt)
                        
            taxis=np.linspace(3600,3600+100*dt, 100)
            plt.plot(taxis,xcorr[3600*1/dt:3600*1/dt+100])
            plt.xlabel('Time (s)')
            plt.ylabel('Cross-correlation amplitude')
            plt.show()

