import os
import time
import matplotlib.pyplot as plt
import numpy as np
from TOOLS.tukey import tukeywin
from obspy import Trace, read, UTCDateTime
#from obspy.geodetics import gps2dist_azimuth
import INPUT.input_synthetics as inp

from UTIL.map_xyz import map_xyz
#from obspy.signal.util import next_pow_2
from obspy.signal.util import nextpow2
from obspy.core.util import gps2DistAzimuth
from math import log
from warnings import warn
from mpi4py import MPI

if __name__=='__main__':
    import synthesize_correlation3 as sc
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    
    print 'Rank '+str(rank)+' started at: '+\
    time.strftime('%d.%m.%Y., %H.%M')
    sc.synthesize_correlation(rank,size)
    print 'Rank '+str(rank)+' finished at: '+\
    time.strftime('%d.%m.%Y., %H.%M')
    
def synthesize_correlation(rank,size):
    
    # =============================================================
    # Preliminaries
    # =============================================================
    if rank == 0 and not os.path.exists(inp.outdir):
        os.mkdir(inp.outdir)
        os.system('cp '+inp.mask+' '+inp.outdir)
        os.system('cp synthesize_correlation2.py '+inp.outdir)
        os.system('cp INPUT/input_synthetics.py '+inp.outdir)
        
    
    # =============================================================
    # Match station pairs
    # =============================================================
    pairs = []
    stations = open(inp.stations,'r')
    stas = stations.read().split('\n')
    i = 0
    
    while i < len(stas):
        sta_0 = stas[i].strip()
        stas_other = stas[i+1:]
        i += 1
        if sta_0 == '': continue
        
        for sta in stas_other:
            if sta == '':
                continue
            else:
                pairs.append([sta_0,sta])
    
    # Got entries in the format
    # [['net1 sta1 lat1 lon1','net2 sta2 lat2 lon2'],[...]]
    
    # =============================================================
    # Read in the noise mask
    # Consider storing noise mask as numpy format!!
    # =============================================================
    mask = np.load(inp.mask)
    print 'Rank %g loaded noise source mask' %rank
    # =============================================================
    # plot the noise mask
    # =============================================================
    if inp.plot_mask == True and rank==0:
        map_xyz(mask[1,:],mask[2,:],mask[3,:])
        os.system('cp new_map.png '+inp.outdir)
    # =============================================================
    # Each entry pair: Extract seismograms on the fly, correlate
    # =============================================================
    #numpairs = int(len(pairs)/size)
    
    #mypairs = pairs[rank*numpairs:(rank+1)*numpairs]
    #try:
    #    mypairs.append(pairs[numpairs*size+rank])
        
    #except IndexError:
    #   pass
    
    mypairs = pairs[rank:len(pairs):size]
    
    for pair in mypairs:
        
        id1 = pair[0]
        id2 = pair[1]
        sta1 = id1.split()[0]+'.'+id1.split()[1]
        sta2 = id2.split()[0]+'.'+id2.split()[1]
    
        # Place sources in the two receiver locations
        lat1 = float(id1.split()[2])
        lon1 = float(id1.split()[3])
        lat2 = float(id2.split()[2])
        lon2 = float(id2.split()[3])
        
    # =============================================================
    # Each entry pair: Check if calculating was already done....
    # =============================================================    
        filename = inp.outdir +sta1+'--'+sta2+'.SAC'
        print filename
        
        if os.path.exists(filename):
            continue
        
        correlation, K, corr_test = get_synthetic_correlation(sta1,sta2,lat1,lat2,lon1,lon2,mask)
        
        
        
    # =============================================================
    # Each entry pair: Store resulting correlations
    # =============================================================
        
        correlation.write(filename=filename,format = 'SAC')
        corr_test.write(filename=filename+'_test',format='SAC')
        
        # Save the kernel
        filename_k = inp.outdir +sta1+'--'+sta2+'.pKern.npy'
        print filename_k
        np.save(filename_k,K) # might introduce to not allow pickling if this is a portability issue
    
def get_synthetic_correlation(sta1,sta2,lat1,lat2,lon1,lon2,mask):
                
    if inp.input == 'specfem_sac' or inp.input == 'specfem_mseed':
        
        mask_i = mask[0,:].reshape(-1)
        mask_x = mask[1,:].reshape(-1)
        mask_y = mask[2,:].reshape(-1)
        mask_z = mask[3,:].reshape(-1)
        
        #if mask_i.ndim > 1:
            
        #    mask_i_temp = []
        #    mask_z_temp = []
        #    mask_x_temp = []
        #    mask_y_temp = []
        #    for j in range(len(mask_i[0,:])):
        #        for k in range(len(mask_i[:,0])):
        #            mask_i_temp.append(mask_i[k,j])
        #            mask_z_temp.append(mask_z[k,j])
        #            mask_y_temp.append(mask_y[k,j])
        #            mask_x_temp.append(mask_x[k,j])
        #    mask_i = np.array(mask_i_temp)
        #    mask_z = np.array(mask_z_temp)
        #    mask_x = np.array(mask_y_temp)
        #    mask_y = np.array(mask_x_temp)
         
        if inp.input == 'specfem_sac':    
            test_seism = read(os.path.join(inp.database,sta1,\
            'OUTPUT_FILES/SRC.00000000.MXZ.sem.sac'))[0]
        else:
            tr1 = read(os.path.join(inp.database,sta1,\
            'OUTPUT_FILES',sta1+'.mseed'))
            tr2 = read(os.path.join(inp.database,sta2,\
            'OUTPUT_FILES',sta2+'.mseed'))
            test_seism = tr1[0]
        
        # Determine the length of the frequency axis for one-sided, real-input fft
        if len(test_seism.data)%2 == 0:
            n = len(test_seism.data)/2 + 1
        else:
            n = (len(test_seism.data)+1)/2
        
        
        # Determine the frequency axis. It is one-sided. Only positive frequencies and zero are kept.
        freq = np.fft.fftfreq(test_seism.stats.npts,d = 1./test_seism.stats.sampling_rate)
        freq = freq[0:n]
        if freq[-1] < 0:
            msg = 'rounding error leading to negative frequencies on freq. axis.\n\
            To avoid this, consider resampling.'
            warn(msg)
            
        correlation = np.zeros(n,dtype='complex')
        correlation_test = np.zeros(len(test_seism.data),dtype='complex')
        
        test_seism.stats.starttime = UTCDateTime(2000,01,01)
        
        
        # Get the source power spectrum shape from, um, somewhere.
        # Maybe actually just define a Gaussian, we can make it more complicated later on.
        if inp.source_spectrum is not None:
            if type(inp.source_spectrum) != tuple or\
            len(inp.source_spectrum) != 2:
                msg = 'Source spectrum must be None or tuple with (mean, half-bandwidth) of \
                Gaussian source spectrum'
                raise ValueError(msg)
            sigma2 = inp.source_spectrum[1] ** 2 / (2. * log(2.))
            print sigma2
            gauss = np.exp(-(freq-inp.source_spectrum[0])**2 / (2. * sigma2))
            
            #plt.plot(freq,gauss,linewidth=1.5)
            #plt.grid()
            #plt.xlabel('Frequency (Hz)',fontsize=16)
            #plt.ylabel('Noise source amplitude (Scaled)',fontsize=16)
            #plt.ylim([0,1.2])
            #plt.xlim([0,0.02])
            #plt.savefig('noise_source_ampspec.eps',format='eps')
            #plt.show()
            S = gauss
        else:
            S = np.ones(n) * inp.source_strength
        
        # Find the relevant frequency indices to save the 'proto'kernel for
        ind_f0 = np.abs(freq-inp.freq_min).argmin()
        ind_f1 = np.abs(freq-inp.freq_max).argmin()
        
        # allocate the kernel array. This could be quite big. 
        # Think about saving the frequency axis for the kernel. Otherwise, it's hard to make sure it is compatible with f(w)
        # Kernel has dimentions (nr. frequencies x nr. source locations+1)
        # The plus one column is for saving the frequency axis
        K = np.zeros(((ind_f1-ind_f0)+2,len(mask_z)+1),dtype='complex')
        
        # Test purposes:
        G1_f0 = np.zeros((3,len(mask_z)+1),dtype='complex')
        G2_f0 = np.zeros((3,len(mask_z)+1),dtype='complex')
        
        K[0,0] = np.nan
        K[1,0] = np.nan
        K[2:,0] = freq[ind_f0:ind_f1]
        K[0,1:] = mask_x 
        K[1,1:] = mask_y
        G1_f0[0:2,:] = K[0:2,:] 
        G2_f0[0:2,:] = K[0:2,:] 
        # ======================================================================
        # Loop over source locations
        # ======================================================================
        count = 0
        for i in mask_i:
            i = int(i)
            z = mask_z[i]
            count += 1
            if count%10000 == 0:
                print 'completed source locations: ',count, ' of ', np.size(mask_z)
            # The format of the filenames for seismograms used here is determined by create_noisemask.py !!
            if inp.input == 'specfem_sac':
                file = 'OUTPUT_FILES/SRC.%08g.' %mask_i[i]
                file1 = os.path.join(inp.database,sta1,file+inp.channel+'.sem.sac')
                file2 = os.path.join(inp.database,sta2,file+inp.channel+'.sem.sac')
                trace1=read(file1)[0]
                trace2=read(file2)[0]
            else:
                srcname = 'SR.%g.S3.' %i
                srcname += inp.channel
                trace1 = tr1[i]#tr1.select(id=srcname)[0]
                trace2 = tr2[i]#tr2.select(id=srcname)[0]
                if int(trace1.stats.station) != i or int(trace2.stats.station) != i :
                    print trace1
                    print trace2
                    print i
                    msg='Inconsistency between Mseed file and source list file.'
                    raise ValueError(msg)
                
            
            
            #### Temporary fix for ugly sac job on daint !!! ####
            for trc in (trace1,trace2):
                trc.stats.starttime=UTCDateTime(2000,01,01)
                trc.detrend('linear')
                trc.taper(type='cosine',max_percentage=0.025)
        
            #    if trc.stats.sampling_rate > 1.:
            #        print 'Downsampling from '+str(trc.stats.sampling_rate)
            #        trc.filter(type='bandpass',freqmin=0.002,freqmax=0.02,\
            #        corners=3,zerophase=True)
            #        trc.decimate(4,no_filter=True)
            #        trc.decimate(4,no_filter=True)
            #        trc.decimate(4,no_filter=True)
            
            #if inp.cutoff == 'R1':
            #    buffer = 20 #samples
            #    distance1=gps2DistAzimuth(lat1,lon1,\
            #    mask_x[i],mask_y[i])[0]/1000.# in km
            #    distance2 = gps2DistAzimuth(lat2,lon2,mask_x[i],mask_y[i])[0]/1000. 
            #    cutoff1 = min(int((distance1 / 3.4 ) * test_seism.stats.sampling_rate + buffer),len(trace1.data))
            #    cutoff2 = min(int((distance2 / 3.4 ) * test_seism.stats.sampling_rate + buffer),len(trace2.data))
            #
            #    win1 = np.zeros(len(trace1.data))
            #    win2 = np.zeros(len(trace2.data))
            #    win1[0:cutoff1] += tukeywin(cutoff1)
            #    win2[0:cutoff2] += tukeywin(cutoff2)
            #    win1[0:cutoff1/2] = 1.
            #    win2[0:cutoff2/2] = 1.
            #    trace1.data = win1*trace1.data
            #    trace2.data = win2*trace2.data
            
            
            
            trace1.filter('bandpass',freqmin=inp.freq_min,freqmax=inp.freq_max,\
            corners=5,zerophase=True)
            trace2.filter('bandpass',freqmin=inp.freq_min,freqmax=inp.freq_max,\
            corners=5,zerophase=True)
            trace1.data = np.ascontiguousarray(trace1.data, np.float32)
            trace2.data = np.ascontiguousarray(trace2.data, np.float32)
            
           
            # Transform to freq. domain. use rfft cause it should be faster than fft
            # Confine to the positive (incl 0) frequency axis, as everything should be Hermitian. Hmm.
           
            pos_spec1 = np.fft.rfft(trace1.data)
            pos_spec2 = np.fft.rfft(trace2.data)
            
            #if i%inp.plot_nsteps == 0:
            #    plt.plot(freq[ind_f0:ind_f1],np.abs(pos_spec1[ind_f0:ind_f1]))
            #    plt.plot(freq[ind_f0:ind_f1],np.abs(pos_spec2[ind_f0:ind_f1]))
            #    plt.show()
            
            # Just for checking once: Store the complex spectral value at f0,
            # at all the locations. 
            
            # calculate G1G2*
            C = pos_spec1 * np.conjugate(pos_spec2)
            C_test = np.correlate(trace1.data,trace2.data,mode='same')
           
            # Store that for the relevant freq. indices as 'proto-kernel'
            K[2:,i+1] += C[ind_f0:ind_f1]
            # Only once for testing: save the complex phase of the two Greens functions at a certain frequency
            G1_f0[2,i+1] = pos_spec1[ind_f0]
            G2_f0[2,i+1] = pos_spec2[ind_f0]
            # Calculate G1G2*Source (where Source is spectrum * spatial weight)
            C *= S * mask_z[i]
            C_test *= mask_z[i]
            if i%inp.plot_nsteps == 0 and i!=0:
                plt.plot(freq,np.abs(C))
                plt.show()
            # Add that.
            correlation += C
            correlation_test += C_test
        
    else: 
        msg = 'This script works only for specfem synthetics in sac or mseed files.\
         Older versions include instaseis. Look up in ARCH folder' 
        raise NotImplementedError(msg)
            
    # Transform back to time domain.
    correlation = np.fft.irfft(correlation)
    print len(correlation)
    # Sort the correlation into positive and negative lag...
    n = len(correlation)
    correlation_temp = np.zeros(n)
    if n%2 == 0:
        print 'correlation length even'
        midsample = n/2
        correlation_temp[0:midsample] = correlation[midsample:]
        correlation_temp[midsample:] = correlation[0:midsample]
    else:
        midsample = int(n/2)
        print 'correlation length odd'
        correlation_temp[0:midsample] = correlation[midsample+1:]
        correlation_temp[midsample+1:] = correlation[0:midsample]
    
    # Test purposes
    np.save(inp.outdir+'test_G1_f0.npy',G1_f0)
    np.save(inp.outdir+'test_G2_f0.npy',G2_f0)
    correlation = correlation_temp
    
    
    # Determine the distance in meters between the two receivers
    dist_meters = gps2DistAzimuth(lat1,lon1,lat2,lon2)[0]
    corr = Trace()
    corr.stats=test_seism.stats
    corr.data = correlation
    corr.stats.sac={}
    corr.stats.sac['dist'] = dist_meters
    corr.stats.sac['stla'] = lat1
    corr.stats.sac['stlo'] = lon1
    corr.stats.sac['evla'] = lat2
    corr.stats.sac['evlo'] = lon2
    
    corr_test = corr.copy()
    corr_test.data = correlation_test
    
    return corr, K, corr_test
