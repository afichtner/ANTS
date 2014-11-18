import os
import gc_geolib as gc
from geographiclib import geodesic,geodesicline

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
    
def seg_example_measr(infile,snr_min=10,nwin_min=1,nwin_max=None,order=1.,\
                      Q=750,v=3800.,freq=0.02,plotstyle='points'):
    
    infile = open(infile,'r')
    data = infile.read().split('\n')
    # output files
    ofid1 = open('measurement_example_1.txt','w')
    ofid2 = open('measurement_stas_1.txt','w')
    
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
        
        if mesr == 'nan':
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
        # Now: one dot per thousand km.
        numseg = int(dist/1000)
        
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
                ofid1.write("%7.2f %7.2f  %7.2f\n" %(seg[1],seg[0],-mesr*seg[2]))
            for seg in seg2:
                #write
                ofid1.write("%7.2f %7.2f  %7.2f\n" %(seg[1],seg[0],mesr*seg[2]))
        if plotstyle == 'gc':
            for i in range(len(seg1)-1):
                seg = seg1[i]
                ofid1.write('> -Z %7.2f\n' %(-mesr*seg[2]))
                ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg1[i+1][1],seg1[i+1][0]))
            for i in range(len(seg2)-1):
                seg = seg2[i]
                ofid1.write('> -Z %7.2f\n' %(mesr*seg[2]))
                ofid1.write("%7.2f %7.2f \n %7.2f %7.2f\n" %(seg[1],seg[0],seg2[i+1][1],seg2[i+1][0]))
        
    ofid1.close()
    ofid2.close()
    
    
    
    if plotstyle == 'points':
        os.system('bash KERNELS/msr_example_1.gmt')
    elif plotstyle == 'gc':
        os.system('bash KERNELS/msr_example_2.gmt')
    os.system('gs -dBATCH -dNOPAUSE -sDEVICE=jpeg -sOutputFile=msr_segments.jpg -r300 msr_segments.ps')
    os.system('rm msr_segments.ps')
    os.system('rm measurement_*.txt')
    print 'Number of successful measurements:'
    print hitcnt
    print 'In percent of total data available:'
    print float(hitcnt)/len(data)*100
    print 'Average number of time windows in each stack:'
    print avgwin/hitcnt
        
        