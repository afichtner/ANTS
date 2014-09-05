#plotting script
from obspy.core import Trace,  read, UTCDateTime
from obspy.core.util.geodetics import gps2DistAzimuth
import os
import TOOLS.read_xml as rxml
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import UTIL.shoot as sht
from glob import glob
import antconfig as cfg

def plot_paths(indir,station,prepname,savefig=True, annotate=True, plotpaths=True, verbose=False):
    """
    A script to plot cross-correlation traces sorted by interstation distance.
    indir: string: path to directory containing results of stack.py
    xmldir: A directory containing the station xml files for all stations (to get the distances). These files must be in xml format and must be called network.station.sta_info.xml
    station: The 'reference': All other stations are plotted with respect to the station selected here and distance is distance to this station
    prepname: string: the 'stamp' of the dataset
    stacktype: string: ls or pws; show the linear or the phase weighted stack (I need to implement first that actually the type is specified on the output file.
    """
    
    #Initialize the map
    # create new figure, axes instances.
    fig=plt.figure(figsize=(16, 8))
    fig.hold()
    ax=fig.add_axes([0.1,0.1,0.8,0.8])
    # setup mercator map projection.
    #m = Basemap(llcrnrlon=-180.,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=80.,\
    #            rsphere=(6378137.00,6356752.3142),\
    #            resolution='c',projection='merc', \
    #            lat_0=0., lon_0=0., lat_ts=20.)
    m = Basemap(rsphere=(6378137.00,6356752.3142),\
                resolution='c',projection='moll', \
                lat_0=0., lon_0=0.)
    m.drawcoastlines()
    m.fillcontinents(lake_color='aqua')
    #m.etopo()
    # draw parallels
    #m.drawparallels(np.arange(-90,90,90))
    # draw meridians
    #m.drawmeridians(np.arange(-180,180,180))
    
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
    m.drawmapboundary(fill_color='aqua')
    
    #==========================================================================================
    #= Find files 
    #==========================================================================================
    stalist=list()
    prepnames=prepname.split(' ')
    for pre in prepnames:
        if station=='all':
            lst=glob(indir+'/*.'+pre+'.*')
        else:
            lst=glob(indir+'/*.'+station+'.*.'+pre+'.*')
        if len(lst)>0:
            stalist.extend(lst)
        
    if len(stalist)==0:
        if verbose: print 'No matching file found.'
        return
    
    
    #==========================================================================================
    #= Find coordinates
    #==========================================================================================
    coordlist=list()
    
    for file in stalist:
        inf=file.split('/')[-1].split('.')
        print inf
        
        #- Find the distance:
        print inf[0],inf[1],inf[4],inf[5]
        coords=rxml.get_coord_dist(inf[0],inf[1],inf[4],inf[5])
        coordlist.append((inf[1],inf[5])+coords)
        
        print cfg.datadir+'/stationxml/'+inf[0]+'.'+inf[1]+".xml"
        print cfg.datadir+'/stationxml/'+inf[4]+'.'+inf[5]+".xml"
        
        stafile1=glob(cfg.datadir+'/stationxml/'+inf[0]+'.'+inf[1]+".xml")[0]
        stafile2=glob(cfg.datadir+'/stationxml/'+inf[4]+'.'+inf[5]+".xml")[0]
       
        if stafile1==stafile2:
           continue
        
        if verbose: 
            print 'Plot path between '+inf[1]+' and '+inf[5]
        
        
        
    #==========================================================================================
    #= Plot great circles
    #==========================================================================================   
    
        (lat1, lon1, lat2, lon2, dist, az, baz)=coords
        #- draw great circle route between sta1, sta2
        if plotpaths==True:
            line,=m.drawgreatcircle(lon1,lat1,lon2,lat2,del_s=50,color='b',linewidth=0.8)
            p = line.get_path()
        # find the index which crosses the dateline (the delta is large)
            cut_point = np.where(np.abs(np.diff(p.vertices[:, 0])) > 10000000)[0]
            print cut_point
            if cut_point:
                cut_point = cut_point[0]
                # create new vertices with a nan inbetween and set those as the path's vertices
                new_verts = np.concatenate(
                                           [p.vertices[:cut_point, :], 
                                            [[np.nan, np.nan]], 
                                            p.vertices[cut_point+1:, :]]
                                               )
                p.codes = None
                p.vertices = new_verts
    
    #==========================================================================================
    #= Plot stations
    #==========================================================================================
    sta_mark=list()
    for stapair in coordlist:
        (sta1,sta2,lat1, lon1, lat2, lon2, dist, az, baz)=stapair   
        x1, y1=m(lon1, lat1)
        x2, y2=m(lon2, lat2)
        ax.plot(x1, y1, 'r.',linewidth=0.8)
        ax.plot(x2, y2, 'r.',linewidth=0.8)
        if annotate==True and sta1 not in sta_mark:
            ax.annotate(sta1, xy=(x1, y1), xytext=(x1, y1))
            sta_mark.extend(sta1)
        if annotate==True and sta2 not in sta_mark:
            ax.annotate(sta2, xy=(x2, y2), xytext=(x2, y2))
            sta_mark.extend(sta2)
        
        
    if station=='all':    
        ax.set_title('Correlations in dataset '+prepname)
    else:
        ax.set_title('Correlations in dataset '+prepname+' w.r.t. station '+station)   
    
    if savefig:
        figname=cfg.testdir+'/Maps/'+station+'.map.'+UTCDateTime().strftime("%Y%m%d")+'.png'
        plt.savefig(figname, format='png', dpi=300)
       
    plt.show()
                     