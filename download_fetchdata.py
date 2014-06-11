import TOOLS.read_xml as rxml
from obspy import UTCDateTime
import os
import sys
import shutil
import download_fetchdata as fd
import antconfig as cfg


if __name__=='__main__':
    xmlin=str(sys.argv[1])
    fd.download_fetchdata(xmlin)
    

def download_fetchdata(xmlinput):
    
    """
    
    Tool for the download continuous seismic data from a collection of stations and/or networks.

    The download is based on the iris DMC FetchData script and takes as input an xml file that specifies the download parameters.

    """
    
    datadir=cfg.datadir
    
     #- read input file ============================================================================

    dat=rxml.read_xml(xmlinput)[1]
   
    
    # Verbose?
    if dat['verbose']=='1':
        v=True
        vfetchdata='-v '
    else:
        vfetchdata=''
        
    # Directory where executable is located
    exdir=dat['exdir']
    
    # network, channel, location and station list
    stafile=dat['ids']

    # time interval of request
    t1=dat['time']['starttime']
    t1str=UTCDateTime(t1).strftime('%Y.%j.%H.%M.%S')
    t2=dat['time']['endtime']
    t2str=UTCDateTime(t2).strftime('%Y.%j.%H.%M.%S')
    
    # minimum length
    minlen=dat['time']['minlen']

    # geographical region
    lat_min=dat['region']['lat_min']
    lat_max=dat['region']['lat_max']
    lon_min=dat['region']['lon_min']
    lon_max=dat['region']['lon_max']
     
    
    # storage of the data
    targetloc=datadir+'raw/latest/'
    respfileloc=datadir+'resp/'
    
    if os.path.isdir(targetloc)==False:
        cmd='mkdir '+targetloc
        os.system(cmd)   
    
    if os.path.isdir(respfileloc)==False:
        cmd='mkdir '+respfileloc
        os.system(cmd)
        
    fh=open(stafile, 'r')
    ids=fh.read().split('\n')
    for id in ids:
        if id=='': continue
        
        #-Formulate a polite request
        filename=targetloc+id+'.'+t1str+'.'+t2str+'.mseed'
        if os.path.exists(filename)==False:
            network=id.split('.')[0]
            station=id.split('.')[1]
            location=id.split('.')[2]
            channel=id.split('.')[3]
            #print network, station, location, channel
            
            reqstring=exdir+'/FetchData '+vfetchdata+' -N '+network+ ' -S '+station+' -L '+location+' -C '+channel+' -s '+t1+' -e '+t2+' -msl '+minlen+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -o '+filename
            os.system(reqstring)
            
            
    # Clean up (some files come back with 0 data)
    cmd=('./UTIL/cleandir.sh '+targetloc)       
    os.system(cmd)
    return
