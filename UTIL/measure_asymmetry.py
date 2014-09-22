import os
import pylab as P
import numpy as np

from glob import glob
from obspy import read
from obspy.core.util.geodetics import gps2DistAzimuth

import antconfig as cfg
import TOOLS.read_xml as rxml

def meas_asym(indir, xmldir,g_speed,w1,w2,filename,ps_nu=0,prefilter=None,\
                verbose=False,outdir='RES/asymmetry_measurements/',doplot=False):
    
    """
    A script to plot cross-correlation traces sorted by interstation distance.
    indir: string: path to directory containing results of stack.py
    xmldir: A directory containing the station xml files for all stations (to get the distances). These files must be in xml format and must be called network.station.sta_info.xml
    g_speed: Fastest group speed in the package of frequencies. This will be used to pick the onset of the correlation signal. In km/sec !!
    w1,w2: Cut w1 seconds before and w2 seconds after the predicted 1-D arrival time.
    ps_nu: integer; power to which the weighting by instantaneous phase stack is raised. ps=0 means linear stack.
    prefilter: Tuple of (float,float,integer) which stands for (lower corner, upper corner, order). A butterworth filter to restrict the correlations to a certain frequency band.
    verbose: Boolean, talk or not
    savefig: Boolean, save plot or not
    outdir: string, where to save plot and outfile. default working directory
    """
    
    if os.path.exists(outdir) == False:
        os.system('mkdir '+outdir)
        
    #Initialize output file
    ofid1=open(outdir+filename+'.txt','w')
    ofid1.write('Input_directory:  '+indir+'\n')
    ofid1.write('Group_speed/evaluated_phase: '+str(g_speed)+'\n')
    ofid1.write('Window_start/end,before/after_arrival:  '+str(w1)+str(w2)+'\n')
    
    if prefilter is not None:
        ofid1.write('Prefilter: '+ str(prefilter)+ '\n')
     
        
    ofid2=open(outdir+filename+'.dat','w')
    files = glob(indir+'*.SAC')

    if doplot == True:
        numwins=list()
    
    
    for file in files:
        if verbose: print file
        (dist,lat1,lon1,lat2,lon2) = get_stainf(file)
        if dist == 0.: continue
        
        correlation=read(file)[0]
        
        numwin = correlation.stats.sac['user0']
        if doplot == True:
            numwins.append(numwin)
        
        ofid1.write(file + '  %g\n' %numwin)
       #win = getwin(dist,Fs)
    if doplot == True:
        plot_hist(numwins)
        
        
        
    
    
    
def get_stainf(corrfile):
    
        inf=corrfile.split('/')[-1].split('.')[0:7]
        st = (cfg.datadir+'/stationxml/'+inf[0]+'.'+inf[1]+'*',\
                cfg.datadir+'/stationxml/'+inf[4]+'.'+inf[5]+'*')
        coord=[(0.,0.),(0.,0.)]

        for i in (0,1):
            try:
                stafile=glob(st[i])[0]
            except IndexError:
                if verbose: print 'No station xml found. Trying to retrieve online.'
                try: 
                  rxml.get_staxml(inf[0+4*i],inf[1+4*i])
                  stafile=glob(st[i])[0]
                except:
                    if verbose: print 'No station xml could be retrieved, skip.'
                    return(0.,0.,0.,0.,0.)
                    
            coord[i]=rxml.find_coord(stafile)[1:3]
        
        lat1 = coord[0][0]
        lon1 = coord[0][1]
        lat2 = coord[1][0]
        lon2 = coord[1][1]
        
        (dist,az,baz)=gps2DistAzimuth(lat1, lon2, lat2, lon2)
        dist/=1000.
        
        return(dist,lat1,lon1,lat2,lon2)
        
def plot_hist(numwins):
        
        # the histogram of the data with histtype='step'
        n, bins, patches = P.hist(numwins, 20, histtype='bar')
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

        # add a line showing the expected distribution
        y = P.normpdf( bins, np.mean(numwins), 100)
        l = P.plot(bins, y, 'k--', linewidth=1.5)
        
        #P.show()
        P.savefig('hist.png')