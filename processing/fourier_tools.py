#a little script to fast fourier transform a time series and return the one sided frequency axis plus the spectrum. 
def fourier_spectrum(datarray,Fs):
    #Fs is sampling rate of time series
    from math import floor,  log
    import scipy.fftpack as fft
    import numpy as np
    
    #Optimal length (this speeds up calculation)
    N=2**(floor(log(len(datarray))/log(2))+1)    
    
    #determine spectrum
    spectrum=fft.fft(datarray,N)
    #frequency axis
    freqaxis=np.arange(N)
    T=N/Fs 
    freqaxis=freqaxis/T
    
#    #determine spectrum
#    spectrum=fft.fft(datarray)
#    freqaxis=np.arange(len(datarray))
#    #Time duration of window that is ffted
#    T=len(datarray)/Fs 
#    freqaxis=freqaxis/T
    
    
    return freqaxis,spectrum



def xcorr_fd(trace1, trace2, delta):
    import numpy as np
    import scipy.fftpack as fft
    import matplotlib.pyplot as plt
    
    spec1=fft.fft(trace1)
    spec2=fft.fft(trace2)
    
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
    taxis=np.arange(0, delta*len(trace1), delta)
    #plt.plot(taxis[0:3000], np.real(crosscorr[0:3000]))
    #plt.show()
    
    return crosscorr
