#plotting script
from obspy.iris import Client
from obspy.core import Trace,  read
from obspy.core.util.geodetics import gps2DistAzimuth
import os
import re
from TOOLS.read_xml import read_xml
import matplotlib.pyplot as plt
import numpy as np

def plot_dist(indir, xmldir, station, channel, corrtype, ps_nu,dist_scaling=10000,maxlag=200, verbose=False, savefig=False, outdir='./'):
    """
    A script to plot cross-correlation traces sorted by interstation distance.
    indir: string: path to directory containing results of stack.py
    xmldir: A directory containing the station xml files for all stations (to get the distances). These files must be in xml format and must be called network.station.sta_info.xml
    station: The 'reference': All other stations are plotted with respect to the station selected here and distance is distance to this station
    channel: string: A rudimentary selection for channels. BHZ, BHE, BHN or BH*
    corrtype: string: ccc or pcc (classical cross correlation, phase cross correlation)
    ps_nu: integer; power to which the weighting by instantaneous phase stack is raised. ps=0 means linear stack.
    dist_scaling: The script calculates the interstation distance in meters. This number has to be scaled down for plotting so that the actual cross-correlations become visible. The amount of scaling depends on the farthest distance in the network on is looking at.
    maxlag: x axis is restricted to +- this value
    """
    #Initialize the plot
    fig=plt.figure(figsize=(10, 18))
    fig.hold()
    plt.subplot(111)
    
    if channel!='BHE' and channel!='BHN' and channel!='BHZ':
        plt.rc('axes', color_cycle=['b', 'g', 'k'])
    else:
        plt.rc('axes', color_cycle=['k'])
    
    
    #Find the files containing a correlation of the 'reference' station
    sta=re.compile(station)
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
    
    #Get time windows of the correlation for information on the plot. 
    mdfile='*'+station+'.*.'+channel+'-'+'*'+station+'.*.'+channel+'.'+corrtype+'.metadata'
    mdfid=open(mdfile, 'r')
    dates=(mdfid.readlines[0], mdfid.readlines[-1])

    #For all these files find the interstation distance, and plot accordingly
    dist_old=0
    for file in stalist:
        if verbose: print file
        inf=file.split('-')[0].split('.')+file.split('-')[1].split('.')
        stafile1=xmldir+'/'+inf[0]+'.'+inf[1]+'.sta_info.xml'
        stafile2=xmldir+'/'+inf[4]+'.'+inf[5]+'.sta_info.xml'
        
        
        correlation=read(indir+'/'+file)[0]
        
        if ps_nu==0:
            stacktype='ls'
        else:
            stacktype='pws'+str(ps_nu)
            fhre=ctype.sub('coherence_stack_real', file)
            fhim=ctype.sub('coherence_stack_imag', file)
            
            phre=read(indir+'/'+fhre)[0].data
            phim=read(indir+'/'+fhim)[0].data
            #Calculate the phase stack
            pstack=np.power(np.sqrt(np.multiply(phre, phre)+np.multiply(phim, phim)), ps_nu)
            
            correlation.data=np.multiply(correlation.data, pstack)
            
        taxis=np.linspace(-(len(correlation.data)-1)/2/correlation.stats.sampling_rate,(len(correlation.data)-1)/2/correlation.stats.sampling_rate, len(correlation.data))
        
        if stafile1==stafile2:
           dist=0
        else:
           inf1=read_xml(stafile1)[1]
           inf2=read_xml(stafile2)[1]
           lat1=float(inf1['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lat'])
           lon1=float(inf1['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lon'])
           lat2=float(inf2['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lat'])
           lon2=float( inf2['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lon'])
           dist=gps2DistAzimuth(lat1, lon1, lat2, lon2)[0]
        correlation.data/=np.max(np.abs(correlation.data))
        
        plt.plot(taxis, correlation.data+dist/dist_scaling)
        
        if dist_old-dist>=1:
            plt.annotate(inf[1]+'-'+inf[5], xy=(taxis[0], dist/dist_scaling) , xytext=(taxis[0]-100, dist/dist_scaling+0.1), fontsize=12 )
        else:
            plt.annotate(inf[1]+'-'+inf[5], xy=(taxis[0], dist/dist_scaling) , xytext=(taxis[-1:]-100, dist/dist_scaling+0.1), fontsize=12 )
        dist_old=dist
            
    plt.xlabel('Lag Time (sec)')
    plt.ylabel("Interstation distance (%g m)" %(dist_scaling))   
    plt.titel(corr_type+" from "+dates[1]+" to "+dates[2]) 
    
    plt.xlim(-maxlag, maxlag)
    frame1=plt.gca()
    #frame1.axes.get_yaxis().set_ticks([])
    
    
    if savefig==True:
        figname=outdir+'/'+station+'.'+channel+'.'+stacktype+'.'+corrtype+'.png'
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
            outfile=indir+'/station_info/'+inf[4]+'.'+inf[5]+'.sta_info.xml'
            #Metadata request with obspy
            if os.path.exists(outfile)==False:
                client.station(inf[4], inf[5], filename=outfile)
        except IndexError:
            continue
    os.system('rm temp.txt') 
          



