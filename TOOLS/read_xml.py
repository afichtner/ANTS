import xml.etree.ElementTree as et
import os
import obspy as obs
from obspy.fdsn import Client
from obspy.core.util.geodetics import gps2DistAzimuth
from glob import glob
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
    tree=et.parse(path_to_xml)
    root=tree.getroot()
    
    sta=path_to_xml.split('/')[-1].split('.')[1]

    lat=root.find('*//Station').find('Latitude').text
    lon=root.find('*//Station').find('Longitude').text
    return sta, float(lat),float(lon)
    
#==============================================================================================
    
def get_staxml(net,sta):
    client=Client()
    outfile=cfg.datadir+'/stationxml/'+net+'.'+sta+'.xml'

    # Metadata request with obspy
    if os.path.exists(outfile)==False:
        client.get_stations(network=net,station=sta,filename=outfile)
        os.system('UTIL/ch_rootelem.sh '+outfile)

#==============================================================================================

def get_coord_staxml(net1, sta1, net2, sta2):


    try:
        stafile1=glob(cfg.datadir+'/stationxml/'+net1+'.'+sta1+'*.xml')[0]
    
    except IndexError:
        print 'Trying to download stationxml nr. 1...'
        try:
            get_staxml(net1,sta1)
        except:
            return(0,0,0,0)
        
    
    try:
        stafile2=glob(cfg.datadir+'/stationxml/'+net2+'.'+sta2+'*.xml')[0]    
        
    except IndexError:
        print 'Trying to download stationxml nr. 2...'
        try:
            get_staxml(net2,sta2)
        except:
            return(0,0,0,0)
            
    (staname1,lat1,lon1)=find_coord(cfg.datadir+'/stationxml/'+net1+'.'+sta1+'.xml')
    (staname2,lat2,lon2)=find_coord(cfg.datadir+'/stationxml/'+net2+'.'+sta2+'.xml')
    return (lat1,lon1,lat2,lon2)

#==============================================================================
def get_geoinf(lat1,lon1,lat2,lon2):
    
    dist=gps2DistAzimuth(lat1, lon1, lat2, lon2)[0]
    az=gps2DistAzimuth(lat1, lon1, lat2, lon2)[1]
    baz=gps2DistAzimuth(lat1, lon1, lat2, lon2)[2]
   
    return (lat1, lon1, lat2, lon2, dist, az, baz)

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
    