#plotting script
from obspy.core import Trace,  read
from obspy.core.util.geodetics import gps2DistAzimuth
import TOOLS.read_xml as rxml
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from glob import glob
import antconfig as cfg
import re


import matplotlib as mpl
#mpl.rcParams['axes.color_cycle'] = ['r', 'k', 'c','b']
mpl.rcParams['axes.color_cycle'] = ['k']

def plot_section(inpattern,ps_nu=0,dist_scaling=10000,\
             maxlag=200,r_speed=None, winwidth=None,excludeauto=True, annotate=False,\
             prefilter=None, scaling=True,verbose=False, savepng=False,\
             savesvg=False,outdir='./'):
    """
    A script to plot cross-correlation traces sorted by interstation distance.
    inpattern: string: path filter containing results of stack.py; containing 
    a pattern like /correlations/dataset/*.pcc.SAC or the like
    ps_nu: integer; power to which the weighting by instantaneous phase stack 
    is raised. ps=0 means linear stack.
    dist_scaling: The script calculates the interstation distance in meters. 
    This number has to be scaled down for plotting so that the actual 
    cross-correlations become visible. The amount of scaling depends on the 
    farthest distance in the network on is looking at.
    maxlag: x axis is restricted to +- this value
    r_speed: For a quick rough check whether this nice wiggle propagates with 
    a velocity consistent with Rayleigh waves
    prefilter: Tuple of (float,float,integer) which stands for (lower corner,
    upper corner, order). A butterworth filter to restrict the correlations 
    to a certain frequency band.
    scaling: Boolean, should all traces be scaled to unity or not
    verbose: Boolean, talk or not
    savefig: Boolean, save plot or not
    outdir: string, where to save plot
    
    """
    # =========================================================================
    # Initialize the plot
    # =========================================================================
    fig=plt.figure(figsize=(12, 7))
    fig.hold()
    ax = fig.add_subplot(111)
    annotpos='r'
    dist_old=0
    
    # =========================================================================
    # what correlations are available?
    # =========================================================================
    corlist = glob(inpattern)
    if len(corlist)==0:
        if verbose: print 'No matching file found. Try another reference station or \'all\'.'
        return
    
    
    # =========================================================================
    # correlation file loop
    # =========================================================================
    for file in corlist:
    
        correlation=read(file)[0]
        if verbose: print file
        
        
    # =========================================================================
    # get station information
    # =========================================================================  
        sta1 = correlation.stats['station']
        sta2 = correlation.stats.sac['kevnm']
        
        lat1 = correlation.stats.sac['stla']
        lon1 = correlation.stats.sac['stlo']
        lat2 = correlation.stats.sac['evla']
        lon2 = correlation.stats.sac['evlo']
        dist = correlation.stats.sac['dist']
       
        if excludeauto == True and dist == 0:
            continue
        
        # =========================================================================
        # phase weighted?
        # =========================================================================
          
        if ps_nu==0:
            stacktype='ls'
        else:
            stacktype='pws'+str(ps_nu)
            psid = re.sub('SAC','npy',file)
            psid = re.sub('\.ccc\.','.ccs.',psid)
            psid = re.sub('\.pcc\.','.pcs.',psid)
            ps = np.load(psid)
            
            ps = np.abs(np.power(ps,ps_nu))
            correlation.data = np.multiply(correlation.data,ps)
            
        if prefilter is not None:
            correlation.taper(max_percentage=0.05, type='cosine')
            correlation.filter('bandpass', freqmin=prefilter[0], \
                        freqmax=prefilter[1], corners=prefilter[2], zerophase=True)
           
        # =========================================================================
        # Scaling and window selection
        # =========================================================================
        
        if scaling==True:
            correlation.data/=np.max(np.abs(correlation.data))
            
        sec1 = (len(correlation.data)/correlation.stats.sampling_rate-1)/2 - maxlag
        correlation.trim(starttime=correlation.stats.starttime+sec1,\
            endtime=correlation.stats.endtime-sec1)
            
        taxis=np.linspace(-(len(correlation.data)-1)/2/correlation.stats.sampling_rate,\
            (len(correlation.data)-1)/2/correlation.stats.sampling_rate, len(correlation.data))
        
    # =========================================================================
    # Plot
    # ========================================================================= 
                
        plt.plot(taxis, correlation.data*dist_scaling+dist)
             
        if annotate==True:
            if dist_old-dist<=5 and annotpos=='l':
                annotpos='r'
            elif dist_old-dist<=5 and annotpos=='r':
                annotpos='l'
                
            if annotpos=='l':
                plt.annotate(sta1+'-'+sta2, xy=(taxis[0], dist+10) , \
                xytext=(taxis[0]-maxlag/6, dist), fontsize=12,fontweight='bold',color='b')
            else:                                                           
                plt.annotate(sta1+'-'+sta2, xy=(taxis[0], dist+10) , \
                xytext=(taxis[-1:], dist), fontsize=12,fontweight='bold',color='b')
        dist_old=dist
    
    
    # =========================================================================
    # Finalize plot with labels
    # =========================================================================
            
    plt.xlabel('Lag Time (s)', fontsize=16, fontweight='bold')
    #plt.ylabel("Interstation distance", fontsize=16, fontweight='bold') 
    xticks = [-maxlag,-3*maxlag/4,-maxlag/2,-maxlag/4,0,maxlag/4,maxlag/2,3*maxlag/4,maxlag] 
    plt.xticks(xticks,fontsize=16,fontweight='bold')
    if scaling==False:
        plt.title("Unscaled Correlations"+'\nFilter: '+str(prefilter),\
         fontsize=16,fontweight='bold')
    else:
        plt.title("Scaled Correlations"+'\nFilter: '+str(prefilter),\
         fontsize=16,fontweight='bold')   
    plt.yticks([],fontsize=16)
    plt.xlim([-maxlag, maxlag])
    plt.axis('tight')
    plt.grid()
    
    # =========================================================================
    # Add a line for estimated Rayleigh wave speed
    # =========================================================================
    
    if r_speed is not None:
        plt.plot(taxis,taxis*r_speed,'r--',linewidth=1.5)
    if winwidth is not None:
        plt.plot(taxis,taxis*r_speed+winwidth,'b--',linewidth=1.2)
        plt.plot(taxis,taxis*r_speed-winwidth,'b--',linewidth=1.2)
    
    # =========================================================================
    # save plot
    # =========================================================================
    if scaling==False:
        sc='.unsc.'
    else:
        sc='.sc.'
        
    
    figname=outdir+'/'+inpattern.strip('/').split('/')[-1]+'.'+\
    stacktype+sc+str(prefilter)
   
    if savesvg==True:
        plt.savefig(figname+'.svg', format='svg', dpi=200)
    if savepng==True:
        plt.savefig(figname+'.png', format='png', dpi=200)
    
    plt.show()