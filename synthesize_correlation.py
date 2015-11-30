import instaseis
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from TOOLS.tukey import tukeywin
from obspy import Trace, read, UTCDateTime
#from obspy.geodetics import gps2dist_azimuth
from UTIL.map_xyz import map_xyz
import INPUT.input_synthetics as inp
#from obspy.signal.util import next_pow_2
from obspy.signal.util import nextpow2
from obspy.core.util import gps2DistAzimuth 

if __name__=='__main__':
    import synthesize_correlation as sc
    print 'Rank '+os.environ[inp.rankvar]+' started at: '+\
    time.strftime('%d.%m.%Y., %H.%M')
    sc.synthesize_correlation()
    print 'Rank '+os.environ[inp.rankvar]+' finished at: '+\
    time.strftime('%d.%m.%Y., %H.%M')
    
def synthesize_correlation():
    
    # =============================================================
    # Preliminaries
    # =============================================================
    rank = int(os.environ[inp.rankvar])
    if rank == 0:
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
    if inp.plot_mask == True:
        map_xyz(mask[1,:,:],mask[2,:,:],mask[3,:,:])
    if rank == 0:
        os.system('cp new_map.png '+inp.outdir)
    # =============================================================
    # Each entry pair: Extract seismograms on the fly, correlate
    # =============================================================
    numpairs = int(len(pairs)/inp.size)
    
    mypairs = pairs[rank*numpairs:(rank+1)*numpairs]
    try:
        mypairs.append(pairs[numpairs*inp.size+rank])
        
    except IndexError:
        pass
    
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
        
        correlation = get_synthetic_correlation(sta1,sta2,lat1,lat2,lon1,lon2,mask)
        
        
        
    # =============================================================
    # Each entry pair: Store resulting correlations
    # =============================================================
        filename = inp.outdir +sta1+'--'+sta2+'.SAC'
        print filename
        correlation.write(filename=filename,format = 'SAC')
           

    
def get_synthetic_correlation(sta1,sta2,lat1,lat2,lon1,lon2,mask):
    

    if inp.method == 'instaseis':
        mask_i = mask[0,:,:]
        mask_x = mask[1,:,:]
        mask_y = mask[2,:,:]
        mask_z = mask[3,:,:]
        db = instaseis.open_db(inp.database)
        source1 = instaseis.ForceSource(latitude=lat1,\
        longitude=lon1,f_r=inp.source_strength)
        source2 = instaseis.ForceSource(latitude=lat2,\
        longitude=lon2,f_r=inp.source_strength)
        
        
        # Initiate correlation
        test_rec = instaseis.Receiver(latitude=0.,longitude=0.,network='X',\
        station='LAE')
        test_seism = db.get_seismograms(source=source1,receiver=test_rec,\
        components='Z')[0]
        
        # Cutoff and zero-pad
        if inp.cutoff == 'R1':
            # longest possible time series: half Earth circumference / group velocity * samplingrate
            t_max_in_samples = 20015. / 3.5 * test_seism.stats.sampling_rate
            # Zeropad for convenience
            t_max_in_samples = next_pow_2(t_max_in_samples)
            correlation = np.zeros(t_max_in_samples*2-1)
        
        else:
            correlation = np.zeros(len(test_seism.data))
        
        # Go through the grid....and calculate seismograms for each grid point
        num_sources=0
        for k in range(np.shape(mask_z)[0]):
            for l in range(np.shape(mask_z)[1]):
                x = mask_x[k,l]
                y = mask_y[k,l]
                z = mask_z[k,l]
                #print x
                #print y
                #print z
                #print '----------'
                receiver = instaseis.Receiver(latitude=x,longitude=y,network='X',\
                station='LAE')
                trace1 = db.get_seismograms(source=source1,receiver=receiver,\
                components='Z')[0]
                trace2 = db.get_seismograms(source=source2,receiver=receiver,\
                components='Z')[0]
                
                # if asked for: determine when Rayleigh wave arrives; taper off afterwards
                if inp.cutoff == 'R1':
                    trace1.data = trace1.data[0:t_max_in_samples]
                    trace2.data = trace2.data[0:t_max_in_samples]
                    distance1 = gps2dist_azimuth(lat1,lon1,x,y)[0]/1000.# in km
                    distance2 = gps2dist_azimuth(lat2,lon2,x,y)[0]/1000. 
                    cutoff1 = min(int((distance1 / 3.5 *2.) * test_seism.stats.sampling_rate),len(trace1.data))
                    cutoff2 = min(int((distance2 / 3.5 *2.) * test_seism.stats.sampling_rate),len(trace2.data))
               
                    win1 = np.zeros(t_max_in_samples)
                    win2 = np.zeros(t_max_in_samples)
                    win1[0:cutoff1] += tukeywin(cutoff1)
                    win2[0:cutoff2] += tukeywin(cutoff2)
                    win1[0:cutoff1/2] = 1.
                    win2[0:cutoff2/2] = 1.
                    trace1.data = win1*trace1.data
                    trace2.data = win2*trace2.data
                
                if inp.noise is not None:
                    trace1.data += np.random.random(len(trace1.data))*inp.noise[1]
                    trace2.data += np.random.random(len(trace2.data))*inp.noise[1]
                    trace1.taper(type='cosine',max_percentage=0.05)
                    trace2.taper(type='cosine',max_percentage=0.05)
                correlation += np.correlate(trace1.data,trace2.data*z,\
                mode='valid')
                
                
    elif inp.method == 'specfem':
        mask_i = mask[0,:]
        mask_x = mask[1,:]
        mask_y = mask[2,:]
        mask_z = mask[3,:]
        test_seism = read(os.path.join(inp.database,sta1,\
        'OUTPUT_FILES/SRC.00000000.MXZ.sem.sac'))[0]
        correlation = np.zeros(len(test_seism.data))
        test_seism.stats.starttime = UTCDateTime(2000,01,01)
        for i in mask_i:
            i = int(i)
            # The format of the filenames here is determined by create_noisemask.py !!
            file = 'OUTPUT_FILES/SRC.%08g.' %mask_i[i]
            z = mask_z[i]
            file1 = os.path.join(inp.database,sta1,file+inp.channel+'.sem.sac')
            file2 = os.path.join(inp.database,sta2,file+inp.channel+'.sem.sac')
            
            trace1=read(file1)[0]
            trace2=read(file2)[0]
            trace1.detrend('linear')
            trace2.detrend('linear')
            #trace1.filter('bandpass',freqmin=0.002,freqmax=0.02,corners=3, zerophase=True)
            #trace2.filter('bandpass',freqmin=0.002,freqmax=0.02,corners=3, zerophase=True)
            
            trace1.data *= inp.source_strength
            trace2.data *= inp.source_strength
            
            correlation += np.correlate(trace1.data,trace2.data*z,\
            mode='same')
           
    
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
    return corr
