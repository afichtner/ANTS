import xml.etree.ElementTree as et
import os
import obspy as obs
from warnings import warn
from glob import glob
from obspy import read_inventory
try:
    from obspy.core.util import gps2DistAzimuth
except ImportError:
    from obspy.geodetics import gps2dist_azimuth
try:
    from obspy.fdsn import Client
except ImportError:
    from obspy.clients.fdsn import Client


import antconfig as cfg

#==============================================================================================

def read_xml(filename):
    
    def recursive_dict(element):
        return element.tag, \
            dict(map(recursive_dict, element)) or element.text
    
    
    doc = et.parse(filename)
    
    return recursive_dict(doc.getroot())

#==============================================================================================

def find_coord(path_to_xml):
    sta = path_to_xml.split('/')[-1].split('.')[1]
    
    try:
        #inf = read_xml(path_to_xml)
        inv = read_inventory(path_to_xml)
    except IOError, TypeError: 
        get_staxml(path_to_xml.split('/')[-1].split('.')[0], sta)
        msg = 'stationxml file not found, trying to download...'
        warn(msg)
        inv = read_inventory(path_to_xml)
        sta = inv[0][0]
        staname = sta.code
        lat=sta.latitude
        lon=sta.longitude
        return str(staname), float(lat),float(lon)
    
    except KeyError: 
        msg='Faulty stationxml file: Could not retrieve coordinates.'
        warn(msg)       
        return '000',0,0

    
#==============================================================================================
    
def get_staxml(net,sta):
    client=Client()
    outfile=cfg.datadir+'/stationxml/'+net+'.'+sta+'.xml'
    print outfile
    # Metadata request with obspy
    if os.path.exists(outfile)==False:
        client.get_stations(network=net,station=sta,filename=outfile)

#==============================================================================================

def get_coord_staxml(net1, sta1):

    try:
        stafile1=glob(cfg.datadir+'/stationxml/'+net1+'.'+sta1+'*.xml')[0]
        (staname1,lat1,lon1)=find_coord(cfg.datadir+'/stationxml/'+net1+'.'+sta1+'.xml')
    except IndexError:
        print 'Trying to download stationxml nr. 1...'
        try:
            get_staxml(net1,sta1)
            stafile1=glob(cfg.datadir+'/stationxml/'+net1+'.'+sta1+'*.xml')[0]
            (staname1,lat1,lon1)=find_coord(cfg.datadir+'/stationxml/'+net1+'.'+sta1+'.xml')
        except:
            return(0,0,0,0)
            
    
    return (lat1,lon1)

#==============================================================================
def get_geoinf(x1,y1,x2,y2,inp='coord'):
    
    if inp == 'coord':
        dist=gps2DistAzimuth(x1, y1, x2, y2)[0]
        az=gps2DistAzimuth(x1, y1, x2, y2)[1]
        baz=gps2DistAzimuth(x1, y1, x2, y2)[2]
   
    return (x1, y1, x2, y2, dist, az, baz)

#==============================================================================================
    
def get_sta_info(stationlist,verbose=False):
    
    fid = open(stationlist,'r')
    idlist = fid.read().split('\n')
    
    for id in idlist:
        if verbose: print id
        if id == '': continue
        
        network = id.split('.')[0]
        sta = id.split('.')[1]
        try:
            get_staxml(network,sta)
        except:
            if verbose: print '\n ================== \
                               No xml downloaded for station \
                               ====================\n'+id
            continue

def get_antip_pt(lon,lat):
    """
    Return the coordinates (lon, lat) of the antipode.
    """
    
    if lon <= 0.:
        lon_new = lon + 180.
    else:
        lon_new = lon - 180.
        
    lat_new = -1.*lat
    
    return lon_new,lat_new
    