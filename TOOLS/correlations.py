
from obspy.core import trace
import numpy as np

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

    xcorr=cross_correlation.xcorr(dat1.data, dat2.data, numsamples ,True)[2]

    return xcorr
   


#==================================================================================================
# Time domain cross-correlation II - only implemented to control my phase cross correlation script
#==================================================================================================
def xcorrelation_td2(dat1, dat2, max_lag, phase):
    
#Max lag in sec --> convert to sample numbers
    Fs=dat1.stats.sampling_rate
    max_lag=int(max_lag)*Fs
    
    s1=dat1.data
    s2=dat2.data
    
    #Correlation window length (in samples)
    T=min(len(s1), len(s2))-2*max_lag

    if T<=0:
        print 'Not enough samples available to calculate correlation at maximum lag.'
        return ()
    
    #Array for the phase crosscorrelation
    xcorr_td2=np.zeros((2*max_lag+1, 1))

    #Summation
    ind=0

    for lag in range(1, max_lag+1):
        lag=abs(lag-(max_lag+1))
        for k in range(max_lag, len(s1)-max_lag):
            xcorr_td2[ind]+=s1[k]*s2[k+lag]
        ind+=1

    for lag in range(0, max_lag+1):
        for k in range(max_lag, len(s1)-max_lag):
            xcorr_td2[ind]+=s1[k+lag]*s2[k]
        ind+=1

    #Normalization
    xcorr_td2/=(2*T)
    xcorr_td2/=np.std(s1)*np.std(s2)
    
    return xcorr_td2, max_lag

#==================================================================================================
# Phase cross correlation (Schimmel 1999)
#==================================================================================================    

def phase_xcorrelation(dat1, dat2, max_lag=10, nu=1):

    from scipy.signal import hilbert

    # Initialize phase correlation as complex array
    n=len(dat1.data)
    print n
    pxc=np.ones(n, dtype=np.complex)
    
    #- Obtain analytic signal
    s1=hilbert(dat1.data)
    s2=hilbert(dat2.data)
    
    #- Normalization
    tol1=np.mean(np.abs(s1))/10.0
    tol2=np.mean(np.abs(s2))/10.0
    s1=s1/(np.abs(s1)+tol1)
    s2=s2/(np.abs(s2)+tol2)

    #- loop over time lags
    for k in range(n):
        s1[k:]

    #- Max lag in sec --> convert to sample numbers
#    Fs=dat1.stats.sampling_rate
#    max_lag=int(max_lag*Fs)


   
    #Correlation window length (in samples)
#    T=min(len(s1), len(s2))-2*max_lag

#    if T<=0:
#        print 'Not enough samples available to calculate correlation at maximum lag.'
#        return ()
    
    #Array for the phase crosscorrelation
#    pxc=np.zeros((2*max_lag+1, ))

    #Summation
#    ind=0

#    for lag in range(1, max_lag):
#        lag=abs(lag-(max_lag+1))
#        for k in range(max_lag, len(s1)-max_lag-1):
#            pxc[ind]+=abs(s1[k]+s2[k+lag])**nu-abs(s1[k]-s2[k+lag])**nu
#        ind+=1
    
    #zero lag
#    for k in range(max_lag, len(s1)-max_lag-1):
#        pxc[ind]+=abs(s1[k]+s2[k])**nu-abs(s1[k]-s2[k])**nu
#    ind+=1

#    for lag in range(1, max_lag):
#        for k in range(max_lag, len(s1)-max_lag-1):
#            pxc[ind]+=abs(s1[k+lag]+s2[k])**nu-abs(s1[k+lag]-s2[k])**nu
#        ind+=1

    #Normalization
#    pxc/=(2*T)
    
#    if phase:
#        s1=np.ndarray(len(pxc), dtype=complex)
#        s1=hilbert(pxc)
#        s1=s1/abs(s1)
#        return s1
#    else:
    return np.abs(pxc)



