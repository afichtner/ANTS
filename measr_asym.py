import os
import pylab as P
import numpy as np
import sys

from glob import glob
from obspy import read
from geographiclib import geodesic,geodesicline
from math import ceil, log
import matplotlib.pyplot as plt

import UTIL.geolib as gl
import measr_inp as minp


if __name__=='__main__':
    import measr_asym as ma
    if len(sys.argv) == 1:
        print 'Usage: python measr_asym.py <argument>'
        print 'Arguments: measure, bin_values, plot_greatcirc'
        print 'Please edit the input file: measr_inp.py\n'
    elif sys.argv[1] == 'measure':
        ma.meas_asym()
    elif sys.argv[1] == 'bin_values':
        ma.seg_measr(plotstyle='points',plot_off=True)
        ma.bin_asym()
    elif sys.argv[1] == 'plot_greatcirc':
        ma.seg_measr(plotstyle='gc')
    else:
        print '\nInvalid input argument!'
        print 'Arguments: measure, bin_values, plot_greatcirc'
        print 'Input file: measr_inp.py\n'

def meas_asym():
    
    """
    input: Provide in input file measr_inp.py
    """
    
    #Initialize output file
    ofid1=open(minp.out_basename+'.msr1.txt','w')
    ofid1.write('Input_pattern:  '+minp.input+'\n')
    ofid1.write('Group_speed/evaluated_phase: '+str(minp.g_speed)+'\n')
    ofid1.write('Window halfwidth: '+str(minp.hw)+' seconds\n')
    ofid1.write('Window_type: '+minp.window+' \n\n')
    
    ofid1.write('Prefilter: '+ str(minp.prefilter)+ '\n\n')
     
        
    ofid2=open(minp.out_basename+'.msr2.txt','w')
    files = glob(minp.input)
    numwins=list()
    if minp.g_speed_msr == None:
        g_speed = minp.g_speed
    else:
        fid = open(minp.g_speed_msr,'r')
        gspeeds = fid.read().split('\n')
        #for entry in gspeeds:
        #    
    
    for file in files:
        
        trace = read(file)[0]
        if minp.ps_nu > 0:
            psfile = file.rstrip('SAC')+'npy'
            # This is a bit messy --- oh well! Let's hope no files have weird names
            psfile = psfile.replace('pcc','pcs',1)
            psfile = psfile.replace('ccc','ccs',1)
            
            try:
                ps = np.load(psfile)
            except IOError:
                print('Could not find phase coherence file. Skipping:')
                print(psfile)
                continue
            ps = np.abs(ps)
            trace.data *= np.power(ps,minp.ps_nu)
       
       
       # -----
        #msr = Nmeasure(trace,minp.g_speed,minp.w1,minp.w2)
        msr, snrc, snra, nw = asym_measr(trace)
        numwins.append(nw)
        id = trace.id + '--' + trace.stats.sac['kuser0'].strip()+\
        '.'+trace.stats.sac['kevnm'].strip()+'.'+\
        trace.stats.sac['kuser1'].strip()+'.'+\
        trace.stats.sac['kuser2'].strip()
        
        ofid1.write( id + '  %g' %nw)
        ofid2.write('%9.4f %9.4f %9.4f %9.4f %12.2f ' \
        %(trace.stats.sac['stla'],trace.stats.sac['stlo'],\
        trace.stats.sac['evla'],trace.stats.sac['evlo'],\
        trace.stats.sac['dist']) )
        
        if minp.doplot == True:
            plot_measr(trace)
        
        if minp.verbose==True:
            print file.split('/')[-1]
        
        ofid2.write('  %3g' %nw)
        ofid2.write('  %10.6f' %snrc)
        ofid2.write('  %10.6f' %snra)
        
        ofid1.write('%10.6f\n' %msr)
        ofid2.write('%10.6f\n' %msr)
        
    #if minp.dohist == True and len(numwins) > 0:
    #    plot_hist(numwins,minp.out_basename)
    
    ofid1.close()
    ofid2.close()
    
def bin_asym():
    dat = open('gmt_scripts/temp/asym_msr.txt','r')
    dat = dat.read().strip().split('\n')
    lons,lats,vals,hits = bin_ind(minp.latmin,minp.latmax,minp.lonmin,\
    minp.lonmax,minp.ddeg_lat,minp.ddeg_lon,dat)
                            
    # write the results
    ofid1 = open('gmt_scripts/temp/vals_xyz.txt','w')
    ofid2 = open('gmt_scripts/temp/hits_xyz.txt','w')
    ofid3 = open('gmt_scripts/temp/info_xyz.txt','w')
    for i in range(0,len(vals)):
        for j in range(0,len(lats)):
            if minp.bin_weight == True:
                areaweight = gl.area_of_sqdeg(lats[j])/gl.area_of_sqdeg(0.)
            else:
                areaweight = 1.
            if hits[i,j]>0:
                ofid1.write("%7.2f %7.2f %7.4f\n" %(lons[i], lats[j],\
                            vals[i,j]/hits[i,j]/areaweight))
            else:
                ofid1.write("%7.2f %7.2f %7.4f\n" %(lons[i], lats[j],\
                            0.))
            ofid2.write("%7.2f %7.2f %7.4f\n" %(lons[i], lats[j],\
            ceil(hits[i,j]/areaweight)))
        
    ofid3.write('input file: '+minp.inp_binning_plotting+'\n')
    ofid3.write('Q = %6.2f \n' %minp.q)
    ofid3.write('freq = %6.2f Hz\n' %minp.f_centr)
    ofid3.write('group v = %6.2f m/s\n' %(minp.g_speed))
    ofid3.write('snr_min: '+str(minp.snr_min)+'\n')
    ofid3.write('sign convention: %6.2f\n' %minp.signconv)
    ofid3.write('great circle segments of %6.2f km length\n' %minp.segper)
    ofid3.write('latmin, lonmin, latmax, lonmax: %6.2f %6.2f %6.2f %6.2f deg\n'\
     %(minp.latmin,minp.lonmin,minp.latmax,minp.lonmax))
    ofid3.write('latitude resolution for binning = %6.2f degree\n'\
     %minp.ddeg_lat)
    ofid3.write('longitude resolution for binning = %6.2f degree\n'\
     %minp.ddeg_lon)
    ofid1.close()
    ofid2.close()
    ofid3.close()



def seg_measr(plotstyle='points',plot_off=False):
    
    print 'Determining great circle segments...'
    infile = open(minp.inp_binning_plotting,'r')
    data = infile.read().split('\n')
    # output files
    ofid1 = open('gmt_scripts/temp/asym_msr.txt','w')
    ofid2 = open('gmt_scripts/temp/asym_stas.txt','w')
    
    # count the valid measurements
    hitcnt = 0
    # count all
    totcnt = 0
    # how many windows, on average...?
    avgwin = 0
    
    for entry in data:
        entry=entry.split()
        
        if len(entry) < 9: continue
        totcnt +=1
        sta_dist = float(entry[4])
        if sta_dist == 0.: continue
        if float(entry[6]) < minp.snr_min and float(entry[7]) < minp.snr_min:
            continue
        if entry[8] == 'nan':
            continue
            
        lat1 = float(entry[0])
        lon1 = float(entry[1])
        lat2 = float(entry[2])
        lon2 = float(entry[3])
        mesr = minp.signconv*float(entry[8])
        hitcnt +=1
        avgwin += int(float(entry[5]))
        
        # write station coordinates to file
        ofid2.write("%8.3f %8.3f \n" %(lon1,lat1))
        ofid2.write("%8.3f %8.3f \n" %(lon2,lat2))
        
        # find midpoint
        mp = gl.get_midpoint(lat1,lon1,lat2,lon2)
        # find antipode of midpoint
        ap = gl.get_antipode(mp[0],mp[1])
        #find distance to antipode of midpoint
        dist = geodesic.Geodesic.WGS84.Inverse(ap[0],ap[1],lat1,lon1)\
        ['s12']/1000.
        dist2 = geodesic.Geodesic.WGS84.Inverse(ap[0],ap[1],lat2,lon2)\
        ['s12']/1000.
        # (half) Nr of segments
        numseg = int(dist/minp.segper)
        if numseg == 0:
            print 'Warning! Zero segments! Setting to 1'
            numseg = 1
        
        # get segments station 1 to antipode
        seg1 = gl.get_gcsegs(lat1,lon1,ap[0],ap[1],numseg,minp.num_max,True,\
        sta_dist,\
        minp.f_centr,minp.q,minp.g_speed)
        # get segments station 2 to antipode
        seg2 = gl.get_gcsegs(lat2,lon2,ap[0],ap[1],numseg,minp.num_max,True,\
        sta_dist,\
        minp.f_centr,minp.q,minp.g_speed)
        
        if plotstyle == 'points':
            for seg in seg1:
                #write
                val = mesr*seg[2]   
                ofid1.write("%7.2f %7.2f  %7.2f\n" %(seg[1],seg[0],-val))
                
            for seg in seg2:
                val = mesr*seg[2]
                ofid1.write("%7.2f %7.2f  %7.2f\n" %(seg[1],seg[0],val))
                
        if plotstyle == 'gc':
            for i in range(len(seg1)-1):
                seg = seg1[i]
                val = mesr*seg[2]
                ofid1.write('> -Z%3.2f\n' %(-val))
                ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],\
                seg1[i+1][1],seg1[i+1][0]))
                
            for i in range(len(seg2)-1):
                seg = seg2[i]
                val = mesr*seg[2]
                ofid1.write('> -Z%3.2f\n' %(val))
                ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],\
                seg2[i+1][1],seg2[i+1][0]))
        
    ofid1.close()
    ofid2.close()
    
    if hitcnt == 0:
        print 'No measurements suit your criteria.'
        return()
    
    if plotstyle == 'points' and plot_off == False:
        os.system('bash gmt_scripts/msr_points.gmt')
    elif plotstyle == 'gc' and plot_off == False:
        os.system('bash gmt_scripts/msr_gcsegs.gmt')
    
    if plot_off == False:
        filename = minp.out_basename+'.jpg'
        os.system('gs -dBATCH -dNOPAUSE -sDEVICE=jpeg -sOutputFile='+filename+\
        ' -r200 gmt_scripts/temp/msr_segments.ps')
   
    print 'Number of successful measurements:'
    print hitcnt
    print 'In percent of total data available:'
    print float(hitcnt)/float(totcnt)*100
    print 'Average number of time windows in each stack:'
    print avgwin/hitcnt

def bin_ind(latmin,latmax,lonmin,lonmax,\
                            ddeg_lat,ddeg_lon,dat):
                            
    lats = np.arange(latmin,latmax+ddeg_lat,ddeg_lat)
    lons = np.arange(lonmin,lonmax+ddeg_lon,ddeg_lon)
    vals = np.zeros((len(lons),len(lats)))
    hits = np.zeros((len(lons),len(lats)))
    
    datalon = np.zeros(len(dat))
    datalat = np.zeros(len(dat))
    dataval = np.zeros(len(dat))
    
    
    for entry in dat:
        if entry.split() == []: continue
        lon = float(entry.split()[0])
        lat = float(entry.split()[1])
        if lon > lonmax: continue
        if lat > latmax: continue
        if lon < lonmin: continue
        if lat < latmin: continue
        val = float(entry.split()[2])
        # Index 1 - longitude index
        i1 = int(round((lon-lonmin)/ddeg_lon))
        if i1 > len(lons)-1: continue
        # Index 2 - latitude index
        i2 = int(round((lat-latmin)/ddeg_lat))
        if i2 > len(lats)-1: continue
        vals[i1,i2] += val
        hits[i1,i2] += 1
        
    return lons,lats,vals,hits
    
def asym_measr(correlation):
    
    if minp.prefilter is not None:
        #correlation.detrend('linear')
        #correlation.detrend('demean')
        correlation.taper(max_percentage=0.02,type='cosine')
        correlation.filter(type='bandpass',freqmin=minp.prefilter[0],\
        freqmax=minp.prefilter[1],corners=minp.prefilter[2],zerophase=True)
        
    win_signl, win_noise, wins = get_wins(correlation)
    
    if wins == True:
        signal = correlation.data*win_signl
        sig_acausal = signal[0:(len(signal)-1)/2]
        sig_causal = signal[(len(signal)-1)/2:len(signal)]
        
        msr = log(np.sum(np.power(sig_causal,2))/np.sum(np.power(sig_acausal,2)))
        
        noise = correlation.data*win_noise
        nse_acausal = noise[0:(len(noise)-1)/2]
        nse_causal = noise[(len(noise)-1)/2:len(noise)]
        
        snrc = np.sum(np.power(sig_causal,2)) / np.sum(np.power(nse_causal,2))
        snra = np.sum(np.power(sig_acausal,2)) / np.sum(np.power(nse_acausal,2))
        
        nw = int(correlation.stats.sac['user0'])
    
        return(msr,snrc,snra,nw)
        
    else:
        if minp.verbose == True:
            print('No measurement windows selected, returning nan.')
            return(np.nan,np.nan,np.nan,0)
            
    
#def plot_hist(numwins,filename):
#        
#        # the histogram of the data with histtype='step'
#        n, bins, patches = P.hist(numwins, 20, histtype='bar')
#        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
#
#        # add a line showing the expected distribution
#        y = P.normpdf( bins, np.mean(numwins), 100)
#        
#        # Labels
#        P.xlabel('Nr. of windows used for stack')
#        P.ylabel('Nr. correlation stacks')
#        P.title(filename.split('/')[:-1])
#        P.savefig(filename + '.hist.png')
#        P.show()
 
 
def plot_measr(correlation):
    
    msr,snrc,snra,nw = asym_measr(correlation)
    win_signl, win_noise, success = get_wins(correlation)
    
    if success == True:
        net = correlation.stats.sac['kuser0'].strip()
        sta = correlation.stats.sac['kevnm'].strip()
        loc = correlation.stats.sac['kuser1'].strip()
        cha = correlation.stats.sac['kuser2'].strip()
        if loc == '-12345':
            loc = ''
        
        id = correlation.id + '--' + net +'.'+sta+'.'+loc+'.'+cha
                
                
        maxlag = (correlation.stats.npts-1)/2/correlation.stats.sampling_rate
        lag = np.linspace(-maxlag,maxlag,len(correlation.data))
        winlen = correlation.stats.sac['user1']
        (x1,y1) = (-maxlag+100,np.min(correlation.data)/2)
        
        plt.plot()
        plt.plot(lag,correlation.data,'k',linewidth=1.7)
        plt.plot(lag,win_signl*np.max(correlation.data),'r--',linewidth=1.5)
        plt.plot(lag,win_noise*np.max(correlation.data)*0.5,'b--',linewidth=1.5)
        
        plt.title(id,fontweight='bold')
        plt.xlabel('Lag (sec)',fontsize=16,fontweight='bold')
        plt.ylabel('Correlation',fontsize=16,fontweight='bold')
        plt.legend(['data','signal window','noise window'])
        plt.annotate('ln(amplitude ratio): %5.4f\ncausal window s/n: %5.4f\
        \nacausal window s/n: %5.4f\nnr. of stacked windows: %g\n\
        window length (s): %g' %(msr,snrc,snra,nw,winlen),\
        xy=(x1,y1),xytext=(x1,y1),bbox=dict(boxstyle="round", fc="0.8"))
                                    
                                    
        plt.xlim([-maxlag,maxlag])
        plt.xticks([-maxlag,-maxlag/2.,0,maxlag/2.,maxlag],fontweight='bold')
        
        plt.show()
     
     
def get_wins(correlation):
    
    # Initialize array for windows
    win_signl = np.zeros(len(correlation.data))
    win_noise = np.zeros(len(correlation.data))
    success = False
    
    # Determine window bounds for signal window
    s_0 = int((len(correlation.data)-1)/2)
    t_lo = int((correlation.stats.sac['dist']/minp.g_speed-minp.hw)*\
    correlation.stats.sampling_rate)
    t_hi = int((correlation.stats.sac['dist']/minp.g_speed+minp.hw)*\
    correlation.stats.sampling_rate)
    w_ind = (s_0-t_hi+1, s_0-t_lo+1, s_0+t_lo, s_0+t_hi)

    if w_ind[2] < w_ind[1] and minp.win_overlap == False:
        if minp.verbose == True:
            print 'No windows found. (Windows overlap) '
        return win_signl, win_noise, success
        
    
    # Construct signal window
    if minp.window == 'boxcar':
         win_signl[w_ind[0]:w_ind[1]] += 1.
         win_signl[w_ind[2]:w_ind[3]] += 1.
    elif minp.window == 'hann':
         win_signl[w_ind[0]:w_ind[1]] += np.hanning(w_ind[1]-w_ind[0])
         win_signl[w_ind[2]:w_ind[3]] += np.hanning(w_ind[3]-w_ind[2])
   
    
    # Determine window bounds for noise window
    noisewinshift = minp.sepsignoise*minp.hw
    t_lo = t_hi + int(noisewinshift*correlation.stats.sampling_rate)
    t_hi = t_lo + int(2*minp.hw*correlation.stats.sampling_rate)
    w_ind = (s_0-t_hi+1, s_0-t_lo+1, s_0+t_lo, s_0+t_hi)
    
    # Out of bounds?
    if w_ind[0] < 0 or w_ind[3] > len(correlation.data):
        if minp.verbose == True:
            print 'No windows found. (Noise window not covered by data)'
        # return two zero arrays - no measurement possible
        return win_noise, win_noise, success

    # Construct noise window
    if minp.window == 'boxcar':
         win_noise[w_ind[0]:w_ind[1]] += 1.
         win_noise[w_ind[2]:w_ind[3]] += 1.
    elif minp.window == 'hann':
         win_noise[w_ind[0]:w_ind[1]] += np.hanning(w_ind[1]-w_ind[0])
         win_noise[w_ind[2]:w_ind[3]] += np.hanning(w_ind[3]-w_ind[2])
    success = True
    
    return win_signl, win_noise, success

    
    