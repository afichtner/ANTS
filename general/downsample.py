def downsample(datstr, channelname, Fs, Fsnew, format):
    import obspy.signal as sig
    from obspy import UTCDateTime
    import matplotlib.pyplot as plt
    import numpy as np
    #print datstr
    #datstr.plot()
    
    
    Fsnew=float(Fsnew)
    
    if Fs==Fsnew:
        print 'Current and new sampling rate are equal. No downsampling performed.'
        return datstr
    
    #Downsampling---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print 'Downsampling...'
    newstr=datstr.copy()
    
    #Plotting
    tr=datstr[0]
    t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
    
    
    #Determine integer decimation 
    dec=int(Fs/Fsnew)
    
    #downsample
    newstr.decimate(dec, no_filter=True)
    
    #Plotting
    tr_new=newstr[0]
    t_new = np.arange(0, tr_new.stats.npts / tr_new.stats.sampling_rate,
                      tr_new.stats.delta)
    plt.plot(t[1:20000], tr.data[1:20000], 'g', label='Raw', alpha=0.9)
    plt.plot(t_new[1:20000/dec], tr_new.data[1:20000/dec], 'r--', label='Downsampled', alpha=0.9)
    
    
    
    #lowpassfilter
    newstr.filter('lowpass', freq=0.4*Fsnew)
    tr_new=newstr[0]
    plt.plot(t_new[1:10000], tr_new.data[1:10000],  'b--', label='Low-passed/Downsampled', alpha=0.9)
    
    plt.xlabel('Time [s]')
    plt.suptitle(tr.stats.starttime)
    plt.legend()
    #plt.xlim(UTCDateTime(2013, 3, 3, 8, 00, 00), UTCDateTime(2013, 3,3, 8, 00, 10) )
    plt.show()
   
    #write to file
    newfilename=channelname+'_downsampled_'+str(int(Fsnew))
    newstr.write(newfilename, format)
    
    
    print 'Done Downsampling.'
    return newstr
    
    

