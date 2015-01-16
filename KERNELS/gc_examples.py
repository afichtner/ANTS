import os
import numpy as np
import gc_geolib as gc
from geographiclib import geodesic,geodesicline
from obspy.core.util.geodetics import gps2DistAzimuth

def gc_example():
    p=geodesic.Geodesic.WGS84.Inverse(40.6, -73.8, 1.4, 104)
    l=geodesic.Geodesic.WGS84.Line(p['lat1'],p['lon1'],p['azi1'])
    num=15
    for i in range(num+1):
      b=l.Position(i*p['s12']/num)
  
def midpoint_example():
    rome = (41.99,12.5)
    lat = [21.02, -12.1, -18.13, 47.58, -25.75]
    lon = [105.87, 282.95, 178.42, 237.67, 28.20]
    ofid=open('example_geomidpoints.txt','w')
    
    for (lat,lon) in zip(lat,lon):
        mp = gc.get_midpoint(rome[0],rome[1],lat,lon)
        ofid.write("%7.2f %7.2f %7.2f %7.2f \n" %(lat,lon,mp[0],mp[1]))
    ofid.close()
    os.system('bash KERNELS/midpoint_example.gmt; rm example_geomidpoints.txt')

def segments_example():
    rome = (41.99,12.5)
    lat = [21.02, -12.1, -18.13, 47.58, -25.75]
    lon = [105.87, 282.95, 178.42, 237.67, 28.20]
    ofid1=open('example_seg1.txt','w')
    ofid2=open('example_seg2.txt','w')
    for (lat,lon) in zip(lat,lon):
        ofid1.write("%7.2f %7.2f \n" %(lat,lon))
        segs = gc.get_gcsegs(rome[0],rome[1],lat,lon,20)
        for seg in segs[1:]:
            ofid2.write("%7.2f %7.2f \n" %(seg[0],seg[1]))
    ofid1.close()
    ofid2.close()    
    os.system('bash KERNELS/segments_example.gmt; rm example_seg?.txt')    
    
def seg_example_majarc():
    rome = (41.99,12.5)
    #lat = [21.02, -12.1, -18.13, 47.58, -25.75]
    #lon = [105.87, 282.95, 178.42, 237.67, 28.20]
    lat = [21.02]
    lon = [105.87]
    ofid0=open('example_seg0.txt','w')
    ofid1=open('example_seg1.txt','w')
    ofid2=open('example_seg2.txt','w')
    for (lat,lon) in zip(lat,lon):
        ofid1.write("%7.2f %7.2f \n" %(lat,lon))
        # get midpoint
        mp = gc.get_midpoint(rome[0],rome[1],lat,lon)
        # get antipode
        ap = gc.get_antipode(mp[0],mp[1])
        ofid0.write("%7.2f %7.2f \n" %(ap[0],ap[1]))
        # first station to antipode
        segs1 = gc.get_gcsegs(rome[0],rome[1],ap[0],ap[1],20)
        # second station to antipode
        segs2 = gc.get_gcsegs(lat,lon,ap[0],ap[1],20)
        # write to file
        for seg in segs1[1:]:
            # distance midpoint - segment
            dist = geodesic.Geodesic.WGS84.Inverse(mp[0],mp[1],seg[0],seg[1])['s12']/1000
            ofid2.write("%7.2f %7.2f  %7.2f\n" %(seg[0],seg[1],dist))
        for seg in segs2[1:]:
            dist = geodesic.Geodesic.WGS84.Inverse(mp[0],mp[1],seg[0],seg[1])['s12']/1000
            ofid2.write("%7.2f %7.2f  %7.2f \n" %(seg[0],seg[1],dist)) 
    ofid0.close() 
    ofid1.close()
    ofid2.close()
    os.system('bash KERNELS/segments_example.gmt; rm example_seg?.txt')
    
    
def seg_example_sensi():
    lat1=33.051490
    lon1=-114.827060 
    name1='GLA'
    lat2=38.229000
    lon2=-86.293800
    name2='WCI'
    numseg = 100
    freq = 0.1
    Q = 230
    v = 4600
    mesr = 1
    ofid1 = open('sens_example.txt','w')
    ofid2 = open('sens_stas.txt','w')
    sta_dist, az, baz = gps2DistAzimuth(lat1,lon1,lat2,lon2) 
    
    # find midpoint
    mp = gc.get_midpoint(lat1,lon1,lat2,lon2)
        
    # find antipode of midpoint
    ap = gc.get_antipode(mp[0],mp[1])
    dist = geodesic.Geodesic.WGS84.Inverse(ap[0],ap[1],lat1,lon1)['s12']/1000.
    dist2 = geodesic.Geodesic.WGS84.Inverse(ap[0],ap[1],lat2,lon2)['s12']/1000.

    
    
    # get segments station 1 to antipode
    seg1 = gc.get_gcsegs(lat1,lon1,ap[0],ap[1],numseg,True,sta_dist,freq,Q,v)
    # get segments station 2 to antipode
    seg2 = gc.get_gcsegs(lat2,lon2,ap[0],ap[1],numseg,True,sta_dist,freq,Q,v)
    
    for i in range(len(seg1)-1):
        seg = seg1[i]
        ofid1.write('> -Z %7.2f\n' %(-mesr*seg[2]))
        ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg1[i+1][1],seg1[i+1][0]))
        
    for i in range(len(seg2)-1):
        seg = seg2[i]
        ofid1.write('> -Z %7.2f\n' %(mesr*seg[2]))
        ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg2[i+1][1],seg2[i+1][0]))
        
    ofid1.close()
    ofid2.write("%7.2f %7.2f %s\n" %(lon1,lat1,name1))
    ofid2.write("%7.2f %7.2f %s\n" %(lon2,lat2,name2))
    ofid2.close()

    
    ofid1 = open('sens_example2.txt','w')
    for i in range(len(seg1)-1):
        seg = seg1[i]
        ofid1.write('> -Z %7.2f\n' %(-10))
        ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg1[i+1][1],seg1[i+1][0]))
        
    for i in range(len(seg2)-1):
        seg = seg2[i]
        ofid1.write('> -Z %7.2f\n' %(-10))
        ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg2[i+1][1],seg2[i+1][0]))
    ofid1.close()
    
    os.system('bash KERNELS/example_sens1.gmt')
    
    
    
def seg_example_measr(infid,snr_min=10,nwin_min=1,nwin_max=None,order=1.,minmsr=0.,\
                      Q=750,v=3800.,freq=0.02,thresh=0.,segper=1000.,plotstyle='points'):
    
    infile = open(infid,'r')
    data = infile.read().split('\n')
    # output files
    ofid1 = open('measurement_example.txt','w')
    ofid2 = open('measurement_stas.txt','w')
    
    # count the valid measurements
    hitcnt = 0
    # how many windows, on average...?
    avgwin = 0
    
    
    for entry in data:
        
        entry=entry.split()
        
        if len(entry) == 0: continue
        sta_dist = float(entry[4])
        if sta_dist == 0.: continue
        if int(entry[5]) < nwin_min: continue
        if nwin_max is not None and int(entry[5]) > nwin_max: continue
        if float(entry[6]) < snr_min and float(entry[7]) < snr_min: continue
        
        lat1 = float(entry[0])
        lon1 = float(entry[1])
        lat2 = float(entry[2])
        lon2 = float(entry[3])
        mesr = order*float(entry[8])
        
        if entry[8] == 'nan':
            continue
        if abs(float(entry[8])) < minmsr:
            continue
        
        hitcnt +=1
        avgwin += int(entry[5])
        # write station coordinates to file
        ofid2.write("%8.3f %8.3f \n" %(lon1,lat1))
        ofid2.write("%8.3f %8.3f \n" %(lon2,lat2))
        # find midpoint
        mp = gc.get_midpoint(lat1,lon1,lat2,lon2)
        
        # find antipode of midpoint
        ap = gc.get_antipode(mp[0],mp[1])
        dist = geodesic.Geodesic.WGS84.Inverse(ap[0],ap[1],lat1,lon1)['s12']/1000.
        dist2 = geodesic.Geodesic.WGS84.Inverse(ap[0],ap[1],lat2,lon2)['s12']/1000.
        
        if dist-dist2>10:
            print 'different distances of station-antipodal point!'
            print dist, dist2
        # (half) Nr of segments
        numseg = int(dist/segper)
        
        if numseg == 0:
            print 'Warning! Zero segments! Setting to 1'
            numseg = 1
        
        # get segments station 1 to antipode
        seg1 = gc.get_gcsegs(lat1,lon1,ap[0],ap[1],numseg,True,sta_dist,freq,Q,v)
        # get segments station 2 to antipode
        seg2 = gc.get_gcsegs(lat2,lon2,ap[0],ap[1],numseg,True,sta_dist,freq,Q,v)
        if plotstyle == 'points':
            for seg in seg1:
                #write
                if mesr*seg[2] < -thresh: 
                    val = -thresh
                elif mesr*seg[2] > thresh:
                    val = thresh
                else:
                    val = mesr*seg[2]
                    
                ofid1.write("%7.2f %7.2f  %7.2f\n" %(seg[1],seg[0],-val))
            for seg in seg2:
                #write
                if mesr*seg[2] < -thresh: 
                    val = -thresh
                elif mesr*seg[2] > thresh:
                    val = thresh
                else:
                    val = mesr*seg[2]
                ofid1.write("%7.2f %7.2f  %7.2f\n" %(seg[1],seg[0],val))
        if plotstyle == 'gc':
            for i in range(len(seg1)-1):
                seg = seg1[i]
                if mesr*seg[2] < -thresh: 
                    val = -thresh
                elif mesr*seg[2] > thresh:
                    val = thresh
                else:
                    val = mesr*seg[2]
                ofid1.write('> -Z %7.2f\n' %(-val))
                ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg1[i+1][1],seg1[i+1][0]))
                
            for i in range(len(seg2)-1):
                seg = seg2[i]
                if mesr*seg[2] < -thresh: 
                    val = -thresh
                elif mesr*seg[2] > thresh:
                    val = thresh
                else:
                    val = mesr*seg[2]
                ofid1.write('> -Z %7.2f\n' %(val))
                ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg2[i+1][1],seg2[i+1][0]))
        
    ofid1.close()
    ofid2.close()
    
    
    
    if plotstyle == 'points':
        os.system('bash KERNELS/msr_example_1.gmt')
    elif plotstyle == 'gc':
        os.system('bash KERNELS/msr_example_2.gmt')
    filename = infid.rstrip('.msr2.txt')+'.jpg'
    os.system('gs -dBATCH -dNOPAUSE -sDEVICE=jpeg -sOutputFile='+filename+' -r600 msr_segments.ps')
    #os.system('rm msr_segments.ps')
    #os.system('rm measurement_*.txt')
    print 'Number of successful measurements:'
    print hitcnt
    print 'In percent of total data available:'
    print float(hitcnt)/len(data)*100
    print 'Average number of time windows in each stack:'
    print avgwin/hitcnt
    
def plot_smoothmap(infid,snr_min=10,nwin_min=1,nwin_max=None,order=1.,minmsr=0.,\
                      Q=150,v=3600.,freq=0.005,segper=1000.):
    #generate the input file from the measurement file                           
    #seg_example_measr(infid,snr_min,nwin_min,nwin_max,order,minmsr,\
    #                      Q,v,freq,segper,plotstyle='points')
                          
    dat = open(infid,'r')
    dat = dat.read().split('\n')
    lats=[]
    lons= []
    vals = []
    hits = []
    coords = []
    
    for entry in dat:
        if entry.split()==[]:
            continue
        lon = float(entry.split()[0])
        lat = float(entry.split()[1])
        val = float(entry.split()[2])
        
        if (lat,lon) not in coords:
            lats.append(lat)
            lons.append(lon)
            vals.append(val)
            hits.append(1)
            coords.append((lat,lon))
        else:
            
            i = np.where(np.array(lats)==lat)
            
            for ind in i[0]:
                if lons[ind]==lon:
                    vals[ind]+=val
                    hits[ind]+=1
            
    outf1 = open('addedmsr.txt','w')
    outf2 = open('addedhit.txt','w')
    for i in range(len(lats)):
        outf1.write("%7.2f %7.2f %7.2f \n" %(lons[i],lats[i],vals[i]))
        outf2.write("%7.2f %7.2f %7.2f \n" %(lons[i],lats[i],hits[i]))
        
    outf1.close()
    outf2.close()