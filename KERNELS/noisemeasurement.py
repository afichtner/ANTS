import numpy as np
from obspy import read
from math import log
import matplotlib.pyplot as plt

class Nmeasure(object):
    """
    A class to obtain noise measurements
    Gets a correlation trace and a group speed
    Holds the information: Where are the stations? (coordinates), what is their distance in m?, 
    what is the signal to noise ratio in the chosen windows w. r. t. the surrounding trace,
    what is the amplitude ratio causal to acausal?
    """
    
    def __init__(self,file,min_snr,min_win,groupspeed,whw,prefilter=None):
    
        self.corr = read(file)[0]
        self.min_snr = min_snr
        self.min_win = min_win
        self.groupspeed = groupspeed
        self.whw = whw
        self.filt = prefilter
        self.id = self.corr.id+'--'+self.corr.stats.sac['kuser0'].strip()+\
                    '.'+self.corr.stats.sac['kevnm'].strip()+'.'+\
                        self.corr.stats.sac['kuser1'].strip()+'.'+\
                            self.corr.stats.sac['kuser2'].strip()
        self.nw = self.corr.stats.sac['user0']
        if self.filt:
            self.corr.filter('bandpass',freqmin=self.filt[0],\
                                freqmax=self.filt[1],corners=self.filt[2],\
                                    zerophase=True)
    
    
    
    def plot(self,savefig=False):
        maxlag = self.corr.stats.sac['e']
        l = len(self.corr.data)
        lag = np.linspace(-maxlag,maxlag,l)
        win = self.getwin_fun()
        msr = self.take_measurement()
        nw = self.nw
        winlen = self.corr.stats.sac['user1']
        (snc,sna) = self.check_snr()[0:2]
        (x1,y1) = (lag[0]+100,np.max(self.corr.data)/2)
        (x2,y2) = (lag[:-1]-100,np.max(self.corr.data)/2)
        
        plt.plot()
        plt.plot(lag,self.corr.data,'k')
        plt.plot(lag,win*np.max(self.corr.data),'r--')
        plt.plot()           
        plt.title(self.id)
        plt.xlabel('Lag (sec)')
        plt.ylabel('Correlation')
        plt.annotate('ln(amplitude ratio): %5.4f\ncausal window s/n: %5.4f\
        \nacausal window s/n: %5.4f\nnr. of stacked windows: %g\n window length (s): %g' \
                            %(msr,snc,sna,nw,winlen),\
                                xy=(x1,y1),xytext=(x1,y1),\
                                    bbox=dict(boxstyle="round", fc="0.8"))
        plt.show()
      
      
        
    def geoinf(self):
        # returns lat1,lon1,lat2,lon2,dist in m
        return (self.corr.stats.sac['stla'],self.corr.stats.sac['stlo'],\
                    self.corr.stats.sac['evla'],self.corr.stats.sac['evlo'],\
                        self.corr.stats.sac['dist'])
                        
    
    def getwin_fun(self):
        # Get a window function; this can be extended to a Gaussian, a hanning...
        if self.nw < self.min_win:
            return None
        else:
            l = len(self.corr.data)
            win = np.zeros(l)
            
            (i0,i1,i2,i3) = self.getwin_ind()
            win[i0:i1] = 1.
            win[i2:i3] = 1.
            return win
                    
    def getwin_ind(self):
        # Determine the bounds of the causal window. (the acausal is just the one symmetric to 0)
        # The window is returned in samples, not seconds
        if self.corr.stats.sac['user0'] < self.min_win:
            return None
            
        else:
            Fs = self.corr.stats.sampling_rate
            dist = self.geoinf()[4]
            
            m = (len(self.corr.data)-1)/2
            t = int(dist/self.groupspeed*Fs)
            w = int(self.whw*Fs)
            win_ind = (m-t-w, m-t+w, m+t-w, m+t+w)            

            # Test: With a simple example, e. g. dist 9000 km and g=3000km/s and 
            # a sampling rate of 1 Hz, do we get the right indices?
            
            return win_ind
            
    def take_measurement(self):
        
        win_ind = self.getwin_ind()
        if win_ind is not None:
            (i0,i1,i2,i3) = win_ind
            
        sig_a = self.corr.data[i0:i1]
        sig_c = self.corr.data[i2:i3]
        
        msr = log(np.sum(np.power(sig_c,2))/np.sum(np.power(sig_a,2)))
        
        return msr
        # test: put in an exactly symmetric function or other known functions
        
    def check_snr(self):
        
        win_ind = self.getwin_ind()
        snrc = 0.0
        snra = 0.0
        
        if win_ind is not None:
            
            Fs = self.corr.stats.sampling_rate
            (i0,i1,i2,i3) = win_ind
            nw = int(self.whw*Fs)
            
            winsig_a = self.corr.data[i0:i1]
            #winnoi_a = np.array(self.corr.data[i0-nw:i0]) + \
            #            np.array(self.corr.data[i1:i1+nw])
            winnoi_a = np.array(self.corr.data[i0-2*nw:i0])
                        
            winsig_c = self.corr.data[i2:i3]
            #winnoi_c = list(self.corr.data[i2-nw:i2]) + \
            #            list(self.corr.data[i3:i3+nw])
            winnoi_c = np.array(self.corr.data[i3:i3+2*nw])
            
            # Test: the winsig and winnoi must have the same length
                        
            snrc = np.sum(np.power(winsig_c,2))/ \
                    np.sum(np.power(winnoi_c,2))

            snra = np.sum(np.power(winsig_a,2))/ \
                    np.sum(np.power(winnoi_a,2))

        return (snrc,snra,win_ind)
            
            
            
                        
    
        
