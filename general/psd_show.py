#A short script to plot the power spectrum of a passed-in time series. Window length is in seconds.
def psdplot(timeseries, samplerate=80, winlength=100, savefig=0, title='Power spectral density'):
    from matplotlib.mlab import psd
    import matplotlib.pyplot as plt
    import numpy as np
    nfft=int(samplerate*winlength)
    (powerspec, freqs)=psd(timeseries,nfft,samplerate,  sides='onesided')
    
    
    h=plt.figure(1)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power spectral density')
    plt.plot(freqs, powerspec, 'k-',linewidth=2)
    plt.axis([0., 2., np.amin(powerspec), 1.05*np.amax(powerspec)])
    if savefig:
        plt.savefig('psd_plot.png')
    plt.show()
