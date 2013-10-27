
from obspy.core import trace
import numpy as np
import time

# different types of cross-correlations and similar (e. g. phase cross-correlation) are collected here
# the idea is to compare the different correlation approaches - - both their timing and their results
# Only two traces of the same length may be used as input


#==================================================================================================
# Frequency domain cross-correlation
#==================================================================================================
def xcorrelation_fd(dat1, dat2):
    import scipy.fftpack as fft
    
    # Remove mean
    dat1.detrend('demean')
    dat2.detrend('demean')

    spec1=fft.fft(dat1.data)
    spec2=fft.fft(dat2.data)
    
    specx=spec1*np.conj(spec2)
    crosscorr=np.real(fft.ifft(specx))
    crosscorr/=max(crosscorr)
    
    return crosscorr


#==================================================================================================
# Time domain cross-correlation
#==================================================================================================
def xcorrelation_td(dat1, dat2, max_lag):
    from obspy.signal import cross_correlation
    
    numsamples=int(float(max_lag)*dat1.stats.sampling_rate)
    print 'obtained number of samples.'
    print time.time()
    
    xcorr=cross_correlation.xcorr(dat1.data, dat2.data, numsamples ,True)[2]
    print 'Cross correlation obtained.'
    print time.time()
    return xcorr
   


#==================================================================================================
# Time domain cross-correlation II - only implemented to check speed
#==================================================================================================
def xcorrelation_td2(dat1, dat2, max_lag):
#Max lag in sec --> convert to sample numbers
    Fs=dat1.stats.sampling_rate
    max_lag=int(max_lag*Fs)
    
    s1=dat1.data
    s2=dat2.data
    
    n=len(s1)
    
    xcorr_td2=np.zeros((2*max_lag+1, ))
    
    #- loop over positive time lags
    for k in range(0, max_lag+1):
        s=np.abs(s1[k:]*s2[:(n-k)])
        xcorr_td2[max_lag+k]=np.sum(s)
    #- loop over negative time lags
    for k in range(1,max_lag+1):
        s=np.abs(s2[k:]*s1[:(n-k)])
        xcorr_td2[max_lag-k]=np.sum(s)
    
    return xcorr_td2


#==================================================================================================
# Phase cross correlation II (Schimmel 1999, with fixed correlation window length)
#==================================================================================================    
from scipy.signal import hilbert
import time

def phase_xcorrelation(dat1, dat2, max_lag=10, nu=1):
    # Initialize arrays
    s1=np.zeros((len(dat1),),  dtype=np.float)
    s2=np.zeros((len(dat1),),  dtype=np.float)
    print 'Arrays initialized.'
    print time.time()
    
    #Obtain analytic signal
    s1=hilbert(dat1.data)
    s2=hilbert(dat2.data)
    print 'Hilbert transform done.'
    print time.time()
    
    #Normalization
    s1=s1/(np.abs(s1))
    s2=s2/(np.abs(s2))
    print 'Normalized.'
    print time.time()
   
    #Max lag in sec --> convert to sample numbers
    Fs=dat1.stats.sampling_rate
    max_lag=int(max_lag*Fs)
    print 'Max lag in samples determined.'
    print time.time()
   
    #Correlation window length (in samples)
    T=min(len(s1), len(s2))-2*max_lag
   
    if T<=0:
        print 'Not enough samples available to calculate correlation at maximum lag.'
        return ()
    
    print 'Window length checked.'
    print time.time()
    
    #Initialize pcc array
    pxc=np.zeros((2*max_lag+1,), dtype=float)
    print 'Array for pcc initialized.'
    print time.time()
    

    #And this is a la Schimmel
    for k in range(-max_lag, max_lag+1):
        
        i1=max_lag+k
        i2=len(s1)-max_lag+k
        
        s=np.abs(s1[i1:i2]+s2[i1-k:i2-k])**nu - np.abs(s1[i1:i2]-s2[i1-k:i2-k])**nu
        pxc[max_lag+k]=np.sum(s)
    
    print 'Correlation obtained.'
    print time.time()   

#Normalization???
    pxc/=(2*T)
    
    return pxc
    print 'Done.'
    print time.time()


