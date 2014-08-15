import numpy as np
from obspy.signal.tf_misfit import cwt
from scipy.signal import hilbert



def cwt(data, dt, fmin, fmax, wl='morlet', w0=12.0):
    
    """
    Returns time-frequency coherence of the input signal
    
    """
    
    data = hilbert(data)
    
    # Get the local phase
    tol = np.max(data)/10000
    data = data/(np.absolute(data)+tol)
    
    # Go to time-scale-domain
    cwt = cwt(data, dt, w0, fmin/10.0, fmax*10,nf=1000)
    
    return cwt
    
    
    
def t_pws(pcc):
    """
    Returns time-domain phase weight.
    """
    
    data = hilbert(data)
    
    # Get the local phase
    tol = np.max(data)/10000
    data = data/(np.absolute(data)+tol)
    
    return data
    

def cwtui(cwt, f_min, f_max,  wl='morlet', w0=12.0, verbose=False):
    
    """
    
    Rudimentary function to approximate a function by summing up it's 
    wavelet transform coefficients.
    
    """
    
    nfreq = cwt[:,0].size
    
    f = np.logspace(np.log10(f_min), np.log10(f_max), nfreq)
    
    
    if wl == 'morlet':
        scale = lambda f: w0 / (2 * np.pi * f)
        psi = lambda t: np.pi ** (-.25) * np.exp(1j * w0 * t) * \
            np.exp(-t ** 2 / 2.)
    else:
        raise ValueError('wavelet type "' + wl + '" not defined!')
        
    
    # Obtain the scaling constant?
    
    aax = np.zeros(nfreq, dtype=np.complex)
    
    for n, _f in enumerate(f):
        
        a = scale(_f)
        # Divide lines by a**1.5
        cwt[n,:] /= (a*np.sqrt(a))
        # axis needed for integration
        aax[n] = a
        
    # sum up the columns:
    cwt = np.trapz(cwt, aax, axis=0)
    cwt=np.real(cwt)
    cwt/=np.max(np.max(np.abs(cwt)))
    
    return cwt