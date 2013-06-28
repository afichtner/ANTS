#Preprocessing module
import numpy as np
import scipy as sp
import obspy.signal as obs
from fourier_tools import fourier_spectrum
import matplotlib.pylab as plt

def onebit(datarray):
    datarray=np.copysign(1, datarray)
    return datarray
    
    
def clip(datarray):
    #This should hopefully result in element-wise multiplication
    rms_amp=np.sqrt(np.mean(datarray*datarray))
    #for i in range(0,len(datarray)-1):
      #  if datarray(i)>=rms_amp:
        #    datarray(i)=rms_amp
    datarray=np.clip(datarray, -rms_amp, rms_amp )
    return datarray
        
        
        
def ram_normal(datarray, windowlength):
    #Window length has to be passed in in nr. of samples. CHECK that the smoothing is done correctly!
    weightarray=np.convolve(abs(datarray), np.ones((windowlength,))/windowlength)[(windowlength-1):]
    datarray=datarray/weightarray
    return datarray
    

def waterlevel(datarray, level):
    #Waterlevel has to be selected as integer multiple of daily rms amplitude
    #Please test with a simple bunch of numbers.
    rms_amp=np.sqrt(np.mean(datarray*datarray))
    while max(datarray) >= level*rms_amp:
        new_rms_amp=np.sqrt(np.mean(datarray*datarray))
        datarray=np.clip(datarray, -new_rms_amp, new_rms_amp )
    return datarray


def whiten(datarray,smoothing, samplerate):
#You can choose between inverse weighting the spectrum by its amplitude spectrum, or by a smoother version
    datarrayft=sp.fftpack.fft(datarray)
    #datarrayft=sp.fftpack.rfft(datarray)
    datarrayft_smooth=obs.util.smooth(abs(datarrayft),smoothing)
    datarrayft=datarrayft/datarrayft_smooth
    #datarray=sp.fftpack.irfft(datarrayft)
    datarray=sp.fftpack.ifft(datarrayft)
   
    return datarray
    

    


