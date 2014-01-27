from __future__ import print_function
import os as os
import numpy as np
import scipy as sp
from obspy.core import read
from obspy.core import trace
import obspy.signal as obs
import matplotlib.pyplot as plt

#==================================================================================================
# ONE-BIT NORMALISATION
#==================================================================================================

def onebit(data,verbose,ofid):

    """
    Perform one-bit normalisation.
    """

    if verbose==True: print('* one-bit normalisation',file=ofid)
    data.data=np.copysign(1, data.data)
    return data


#==================================================================================================
# RMS CLIPPING
#==================================================================================================

def clip(data,verbose,ofid):

    """
    Clip data at their rms value.
    """

    if verbose==True: print('* rms clipping',file=ofid)
    rms_amp=np.sqrt(np.mean(data.data*data.data))
    data.data=np.clip(data.data, -rms_amp, rms_amp )
    return data
    

#==================================================================================================
# MOVING AVERAGE NORMALISATION
#==================================================================================================

def ram_normal(data,window_length,verbose,ofid):

    """
    Running average normalisation with window length in seconds. Divide data by running average.
    """

    if verbose==True: print('* running average with '+str(window_length)+' s window',file=ofid)

    #- number of samples in the window, as an odd nr of samples
    N=np.round(window_length/data.stats.delta)
    N=2*int(np.round(float(N)/2.0))+1
    
    #- make weighting array
    weightarray=np.convolve(np.abs(data.data),np.ones((N,))/N)
    weightarray=weightarray[(N/2-1):(N/2-1)+len(data.data)]

    #- normalise
    data.data=data.data/weightarray

    return data


#==================================================================================================
# MOVING AVERAGE NORMALISATION
#==================================================================================================

def waterlevel(data,level,verbose,ofid):

    """
    Waterlevel normalisation. Iteratively clip data at a multiple of their rms value. 
    Continues until all data points are below that multiple of the rms amplitude.
    """

    if verbose==True: print('* iterative clipping at '+str(level)+' times the rms value',file=ofid)

    rms_amp=np.sqrt(np.mean(data.data*data.data))
    while max(np.abs(data.data)) >= level*rms_amp:
        new_rms_amp=np.sqrt(np.mean(data.data*data.data))
        data.data=np.clip(data.data, -new_rms_amp, new_rms_amp )

    return data


#==================================================================================================
# SPECTRAL WHITENING
#==================================================================================================

def whiten(data,smoothing,verbose,ofid):

    """
    Perform spectral whitening. The amplitude spectrum can be smoothed.
    """

    if verbose==True: print('* spectral whitening',file=ofid)

    datarrayft=sp.fftpack.fft(data.data)
    datarrayft_smooth=obs.util.smooth(abs(datarrayft),smoothing)
    datarrayft=datarrayft/datarrayft_smooth
    data.data=np.real(sp.fftpack.ifft(datarrayft))

    return data


    


