#plotting script
from obspy.iris import Client
from obspy.core import Trace,  read
from obspy.core.util.geodetics import gps2DistAzimuth
import os
import re
from TOOLS.read_xml import read_xml
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

def plot_dist(indir, xmldir, station, channel, corrtype, ps_nu,dist_scaling=10000,maxlag=200, r_speed=None, prefilter=None, verbose=False, savefig=False, outdir=None, indir2=None):
    """
    A script to plot cross-correlation traces sorted by interstation distance.
    indir: string: path to directory containing results of stack.py
    xmldir: A directory containing the station xml files for all stations (to get the distances). These files must be in xml format and must be called network.station.sta_info.xml
    station: The 'reference': All other stations are plotted with respect to the station selected here and distance is distance to this station
    channel: string: A rudimentary selection for channels. BHZ, BHE, BHN or BH* (or LH)
    corrtype: string: ccc or pcc (classical cross correlation, phase cross correlation), or 'both' (comparison will be shown)
    ps_nu: integer; power to which the weighting by instantaneous phase stack is raised. ps=0 means linear stack.
    dist_scaling: The script calculates the interstation distance in meters. This number has to be scaled down for plotting so that the actual cross-correlations become visible. The amount of scaling depends on the farthest distance in the network on is looking at.
    maxlag: x axis is restricted to +- this value
    r_speed: For a quick rough check whether this nice wiggle propagates with a velocity consistent with Rayleigh waves
    prefilter: Tuple of (float,float,integer) which stands for (lower corner, upper corner, order). A butterworth filter to restrict the correlations to a certain frequency band.
    verbose: Boolean, talk or not
    savefig: Boolean, save plot or not
    outdir: string, where to save plot
    
    """
    #Initialize the plot
    fig=plt.figure(figsize=(12, 16))
    fig.hold()
    plt.subplot(111)
    
    if channel=="BH*" or channel=="LH*":
        plt.rc('axes', color_cycle=['b', 'g', 'k'])
    else:
        plt.rc('axes', color_cycle=['k'])
    
    #Is a comparison between different types required?
    if corrtype=='both':
        corrtype1='ccc'
        corrtype2='pcc'
    
    #Find the files containing a correlation of the 'reference' station
    sta=re.compile(station)
    if corrtype=='both':
        ctype=re.compile(corrtype1+'_stack')
    else:
        ctype=re.compile(corrtype+'_stack')
    cha=re.compile(channel)
    
    ls=os.listdir(indir)
    stalist=[]
    
    for filename in ls:
        if sta.search(filename) is not None:
            if ctype.search(filename) is not None:
                if cha.search(filename) is not None:
                    
                    stalist.append(filename)
        else:
            continue
            
    if len(stalist)==0:
        if verbose: print 'No matching file found. Try another reference station.'
        return
    
   

    #For all these files find the interstation distance, and plot accordingly
    dist_old=0
  
    for file in stalist:
        
        if verbose: print file
        inf=file.split('-')[0].split('.')+file.split('-')[1].split('.')
        
        try:
            stafile1=glob(xmldir+'/'+inf[0]+'.'+inf[1]+'*')[0]
            stafile2=glob(xmldir+'/'+inf[4]+'.'+inf[5]+'*')[0]
        except IndexError:
            if verbose: print 'No station xml found. Skipping this correlation.'
            continue
    
    
        
        correlation=read(indir+'/'+file)[0]
        
        if corrtype=='both':
            #print indir2+'/'+file.rstrip('.SAC').rstrip('.MSEED').rstrip('ccc_stack')+'pcc_stack'
            try:
                correlation2=read((glob(indir2+'/'+file.rstrip('.SAC').rstrip('.MSEED').rstrip('ccc_stack')+'pcc_stack.*')[0]))[0]
            except IndexError:
                print 'notfound'
                continue
            
     
        
        
        
        if prefilter is not None:
            correlation.taper(p=0.1)
            correlation.filter('bandpass', freqmin=prefilter[0], freqmax=prefilter[1], corners=prefilter[2], zerophase=True)
            if 'correlation2' in locals():
                correlation2.taper(p=0.1)
                correlation2.filter('bandpass', freqmin=prefilter[0], freqmax=prefilter[1], corners=prefilter[2], zerophase=True)
        
        
        if ps_nu==0:
            stacktype='ls'
        else:
            stacktype='pws'+str(ps_nu)
            fhre=ctype.sub('coherence_stack_real', file)
            fhim=ctype.sub('coherence_stack_imag', file)
            
            phre=read(indir+'/'+fhre)[0]
            phim=read(indir+'/'+fhim)[0]
            if prefilter is not None:
                phre.taper(p=0.1)
                phre.filter('bandpass', freqmin=prefilter[0], freqmax=prefilter[1], corners=prefilter[2], zerophase=True)
                phim.taper(p=0.1)
                phim.filter('bandpass', freqmin=prefilter[0], freqmax=prefilter[1], corners=prefilter[2], zerophase=True)
            phre=phre.data
            phim=phim.data
            
            #Calculate the phase stack
            pstack=np.power(np.sqrt(np.multiply(phre, phre)+np.multiply(phim, phim)), ps_nu)
            
            correlation.data=np.multiply(correlation.data, pstack)
            if 'correlation2' in locals():
                correlation2.data=np.multiply(correlation2.data, pstack)
            
            
        taxis=np.linspace(-(len(correlation.data)-1)/2/correlation.stats.sampling_rate,(len(correlation.data)-1)/2/correlation.stats.sampling_rate, len(correlation.data))
        
        if stafile1==stafile2:
           dist=0
           
        else:
            try:
               inf1=read_xml(stafile1)[1]
               inf2=read_xml(stafile2)[1]
               lat1=float(inf1['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lat'])
               lon1=float(inf1['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lon'])
               lat2=float(inf2['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lat'])
               lon2=float( inf2['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lon'])
               dist=gps2DistAzimuth(lat1, lon1, lat2, lon2)[0]
            except IOError:
                if verbose: print 'No station xml file found, skipping this station.'
                
        
        correlation.data/=np.max(np.abs(correlation.data))
        plt.plot(taxis, correlation.data+dist/dist_scaling)
        if 'correlation2' in locals():
            correlation2.data/=np.max(np.abs(correlation2.data))
            plt.plot(taxis, correlation2.data+dist/dist_scaling, 'r')
        
        if dist_old-dist>=1:
            plt.annotate(inf[1]+'-'+inf[5], xy=(taxis[0], dist/dist_scaling) , xytext=(taxis[0]-100, dist/dist_scaling+0.1), fontsize=12 )
        else:
            plt.annotate(inf[1]+'-'+inf[5], xy=(taxis[0], dist/dist_scaling) , xytext=(taxis[-1:]-100, dist/dist_scaling+0.1), fontsize=12)
        if dist>dist_old: max_dist=dist
        dist_old=dist
            
    plt.xlabel('Lag Time (sec)')
    plt.ylabel("Interstation distance (%g m)" %(dist_scaling))   
  
    plt.title(corrtype+" from "+indir.strip('/').split('/')[-1]+'\nPrefilter: '+str(prefilter))
    if r_speed is not None:
        plt.plot(taxis, np.abs(taxis*r_speed/dist_scaling), linewidth=5.0, color='0.8')
    plt.xlim(-maxlag, maxlag)
    frame1=plt.gca()
    
   
    if savefig==True:
        if outdir==None:
            figname=indir+'/'+indir.strip('/').split('/')[-1]+'.'+station+'.'+channel+'.'+stacktype+'.'+corrtype+'.png'
        else:
            figname=outdir+'/'+indir.strip('/').split('/')[-1]+'.'+station+'.'+channel+'.'+stacktype+'.'+corrtype+'.'+str(prefilter)+'.png'
        
        plt.savefig(figname, format='png', dpi=200)
        
    plt.show()
                     
            
            
            

    
def get_sta_info(indir, filt):
    client=Client()
    if os.path.exists(indir+'station_info')==False: os.mkdir(indir+'station_info')
    
    os.system('ls '+indir+' | grep '+ filt+' > temp.txt')
    fileid=open('temp.txt')
    filelist=fileid.read().split('\n')
    
    for file in filelist:
        try:
            inf=file.split('-')[0].split('.')+file.split('-')[1].split('.')
            outfile=indir+'station_info'+inf[4]+'.'+inf[5]+'.sta_info.xml'
            #Metadata request with obspy
            if os.path.exists(outfile)==False:
                client.station(inf[4], inf[5], filename=outfile)
        except IndexError:
            continue
    os.system('rm temp.txt') 
          



