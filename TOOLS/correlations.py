
from obspy.core import trace
from math import sqrt
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
    xcorr=cross_correlation.xcorr(dat1.data, dat2.data, numsamples ,True)[2]
    
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
    
    #- loop over time lags
    for k in range(0, max_lag+1):
        
        xcorr_td2[max_lag+k]=np.sum(s1[k:]*s2[:(n-k)])
        norm[max_lag+k]=sqrt(np.sum(s1[k:]))
        xcorr_td2[max_lag-k]=np.sum(s2[k:]*s1[:(n-k)])
        #norm[max_lag-k]
   
    return xcorr_td2


#==================================================================================================
# Phase cross correlation (Schimmel 1999)
#==================================================================================================    
from scipy.signal import hilbert
import time

def phase_xcorrelation(dat1, dat2, max_lag=10, nu=1, varwl=True):
    
    # Initialize arrays
    s1=np.zeros((len(dat1),),  dtype=np.float)
    s2=np.zeros((len(dat1),),  dtype=np.float)
    
    
    #Obtain analytic signal
    s1=hilbert(dat1.data)
    s2=hilbert(dat2.data)
    
    #Normalization
    s1=s1/(np.abs(s1))
    s2=s2/(np.abs(s2))
    
    
    #Max lag in sec --> convert to sample numbers
    Fs=dat1.stats.sampling_rate
    max_lag=int(max_lag*Fs)
   
    if varwl==True:
        return phase_xcorr(s1,s2,max_lag,nu)
   
    else:
         #Correlation window length (in samples)
         T=min(len(s1), len(s2))-2*max_lag
        
         if T<=0:
             print 'Not enough samples available to calculate correlation at maximum lag.'
             return ()
         else:
             return phase_xcorr_eqweight(s1,s2,max_lag,nu) 
         
    
#==================================================================================================
# Phase cross correlation (Schimmel 1999); this is obtained with a variable window length
#=================================================================================================  

def phase_xcorr(data1,data2,max_lag,nu=1):
    """
    
    data1, data2: Numpy arrays containing the analytic signal normalized sample by sample
    by their absolute value (ie containing only the instantaneous phase information)
    max_lag: maximum lag in number of samples, integer
    
    """
    
    #Initialize pcc array:
    pxc=np.zeros((2*max_lag+1,), dtype=float)
    
   
    for k in range(0,max_lag+1):
        i11=0
        i12=len(data1)-k
        i21=k
        i22=len(data1)
        
        
        pxc[max_lag+k]=1.0/float(2*len(data1)-k)*(np.sum(np.abs(data1[i11:i12]+data2[i21:i22])**nu) - np.sum(np.abs(data1[i11:i12]-data2[i21:i22])**nu))
        pxc[max_lag-k]=1.0/float(2*len(data1)-k)*(np.sum(np.abs(data1[i21:i22]+data2[i11:i12])**nu) - np.sum(np.abs(data1[i21:i22]-data2[i11:i12])**nu))
        
    return pxc
    

#==================================================================================================
# Phase cross correlation (Schimmel 1999); Here an equal window length is used for every lag. This comes at the expense that data from one of the traces are thrown away at each end (max_lag samples at each end)
#=================================================================================================   
    
def phase_xcorr_eqweight(data1,data2,max_lag,nu=1):
    """
    
    data1, data2: Numpy arrays containing the analytic signal normalized sample by sample
    by their absolute value (ie containing only the instantaneous phase information)
    max_lag: maximum lag in number of samples, integer
    
    """

    # Initialize pcc array
    pxc=np.zeros((2*max_lag+1,), dtype=float)
    
    #And this is a la Schimmel
    for k in range(-max_lag, max_lag+1):
        
        i1=max_lag+k
        i2=len(s1)-max_lag+k
        
        pxc[max_lag+k]=np.sum(np.abs(s1[i1:i2]+s2[i1-k:i2-k])**nu) - np.sum(np.abs(s1[i1:i2]-s2[i1-k:i2-k])**nu)
        
    #Normalization
    pxc/=float(2*(len(s1)-2*max_lag))
           
    return pxc
    
#==================================================================================================
# Our own simple analytic function subroutine
#=================================================================================================


def asig(data):
    from math import floor
    from scipy.fftpack import fft, ifft
    """
    Determine the analytic signal 'by hand'
    data: numpy array, containing the data to be transformed.
    """
    
    
    n=len(data)
    nhalf=int(floor(n/2))
    
    s=fft(data)
    s[1:nhalf+1]*=2
    s[nhalf:]=0
    
    return ifft(s)
    
    
    
    
    
#==================================================================================================
# Discarded parts
#=================================================================================================
    
  # fid=open('asig1.txt','w')
  # for i in range(len(s1)):
  #     fid.write(str(i))
  #     fid.write(str(s1[i]))
  #     fid.write("\n")
  # fid.close()
  #     
  # fid=open('asig2.txt','w')
  # for i in range(len(s2)):
  #     fid.write(str(i))
  #     fid.write(str(s2[i]))
  #     fid.write("\n")
  # fid.close()
