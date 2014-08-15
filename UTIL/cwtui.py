import numpy as np

def cwtui(cwt, f_min, f_max, w0, wl='morlet', verbose=False):
    
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