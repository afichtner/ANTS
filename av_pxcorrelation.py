# Phase cross correlation of signal 1 and 2
# Reference: Martin Schimmel (1999): Phase cross-correlations: design, comparisons and applications, BSSA 89.

def pcc(datarray1, datarray2, Fs, max_lag=10, nu=1):
    from scipy.signal import hilbert
    import numpy as np
    # Initialize as complex arrays
    s1=np.zeros((len(datarray1), 1),  dtype=np.complex)
    s2=np.zeros((len(datarray1), 1),  dtype=np.complex)
    
    #Obtain analytic signal
    s1=hilbert(datarray1)
    s2=hilbert(datarray2)
    
    #Normalization
    s1=s1/abs(s1)
    s2=s2/abs(s2)
    
    #Min and max lag are in seconds --> convert to sample numbers
    max_lag=int(max_lag*Fs)
    
    #Correlation window length (in samples)
    T=min(len(s1), len(s2))-2*max_lag

    if T<=0:
        print 'Not enough samples available to calculate correlation at maximum lag.'
        return ()
    
    pxc=np.zeros((2*max_lag+1, 1))
    
    #Summation
    ind=0
    for lag in range(-max_lag, 0):
        for k in range(max_lag, len(s1)-max_lag):
            pxc[ind]+=abs(s1[k]+s2[k+lag])**nu-abs(s1[k]-s2[k+lag])**nu
                   
        ind+=1
    
    for k in range(max_lag, len(s1)-max_lag+1):
        pxc[ind]+=abs(s1[k]+s2[k])**nu-abs(s1[k]-s2[k])**nu
    ind+=1
    
    for lag in range(1, max_lag+1):
        for k in range(max_lag, len(s1)-max_lag):
            pxc[ind]+=abs(s1[k+lag]+s2[k])**nu-abs(s1[k+lag]-s2[k])**nu
                   
        ind+=1
    
    
    pxc/=(2*T)
    
    return pxc
