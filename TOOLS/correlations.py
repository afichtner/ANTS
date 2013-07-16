#different types of cross-correlations and similar (e. g. phase cross-correlation) are collected here
# the idea is to compare the different correlation approaches - - both their timing and their results

import numpy as np
from obspy.core import trace


#==================================================================================================
# Frequency domain cross-correlation
#==================================================================================================
def xcorrelation_fd(dat1, dat2):
    import scipy.fftpack as fft
 
#Remove mean
    dat1.detrend('demean')
    dat2.detrend('demean')
    spec1=fft.fft(dat1)
    spec2=fft.fft(dat2)
    
    if len(spec1)>len(spec2): 
        print 'Warning: Different length arrays. First array is longer by nr. of samples:'
        print len(spec1)-len(spec2)
        spec1=spec1[:len(spec2)]
        
    if len(spec1)<len(spec2): 
        print 'Warning: Different length arrays. First array is shorter by nr. of samples:'
        print len(spec1)-len(spec2)
        spec2=spec2[:len(spec1)]    
    
    specx=spec1*np.conj(spec2)
    crosscorr=np.real(fft.ifft(specx))
    crosscorr/=max(crosscorr)

    return crosscorr


#==================================================================================================
# Time domain cross-correlation
#==================================================================================================
def xcorrelation_td(dat1, dat2, max_lag):
    from obspy.signal import cross_correlation
    
    numsamples=int(max_lag*dat1.stats.sampling_rate)
    xcorr=cross_correlation.xcorr(dat1.data, dat2.data, numsamples,True)[2]

    return xcorr
    
    
#==================================================================================================
# Phase cross correlation (Schimmel 1999)
#==================================================================================================    

def phase_xcorrelation(dat1, dat2, max_lag=10.0, nu=1):
    from scipy.signal import hilbert
    # Initialize as complex arrays
    s1=np.zeros((len(dat1), 1),  dtype=np.complex)
    s2=np.zeros((len(dat1), 1),  dtype=np.complex)
    
    #Obtain analytic signal
    s1=hilbert(dat1.data)
    s2=hilbert(dat2.data)
    
    #Normalization
    s1=s1/abs(s1)
    s2=s2/abs(s2)
    
    #Max lag in sec --> convert to sample numbers
    Fs=dat1.stats.sampling_rate
    max_lag=int(max_lag*Fs)
    
    #Correlation window length (in samples)
    T=min(len(s1), len(s2))-2*max_lag

    if T<=0:
        print 'Not enough samples available to calculate correlation at maximum lag.'
        return ()
    
    #Array for the phase crosscorrelation
    pxc=np.zeros((2*max_lag+1, 1))

    #Summation
    ind=0

    for lag in range(1, max_lag+1):
        lag=abs(lag-(max_lag+1))
        for k in range(max_lag, len(s1)-max_lag):
            pxc[ind]+=abs(s1[k]+s2[k+lag])**nu-abs(s1[k]-s2[k+lag])**nu
        ind+=1

    for lag in range(0, max_lag+1):
        for k in range(max_lag, len(s1)-max_lag):
            pxc[ind]+=abs(s1[k+lag]+s2[k])**nu-abs(s1[k+lag]-s2[k])**nu
        ind+=1

    #Normalization
    pxc/=(2*T)
    
    return pxc





