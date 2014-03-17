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
    
    Tool for the download of response information, to be stored in RESP files.

    The download is based on the iris DMC FetchData script and takes as input an xml file that specifies the download parameters.

    """
    
    datadir=cfg.datadir
    
     #- read input file ============================================================================

    datainput=rxml.read_xml(xmlinput)
    dat1=datainput[1]
    
    
    # Verbose?
    if dat1['verbose']=='1':
        v=True
        vfetchdata='-v '
    else:
        vfetchdata=''
        
    # Directory where executable is located
    exdir=dat1['exdir']
    
    # network, channel, location and station list
    networks=dat1['data']['networks'].strip().split(' ')
    channels=dat1['data']['channels'].strip().split(' ')
    stafile=dat1['data']['stations']

    # time interval of request
    t1=dat1['time']['starttime']
    t1str=UTCDateTime(t1).strftime('%Y.%j.%H.%M.%S')
    t2=dat1['time']['endtime']
    t2str=UTCDateTime(t2).strftime('%Y.%j.%H.%M.%S')
    
    # geographical region
    lat_min=dat1['region']['lat_min']
    lat_max=dat1['region']['lat_max']
    lon_min=dat1['region']['lon_min']
    lon_max=dat1['region']['lon_max']
     
    
    for network in networks:
        # storage of the data
        targetloc=datadir+'raw/'+network
        respfileloc=datadir+'resp/'
        
        if os.path.isdir(targetloc)==False:
            cmd='mkdir '+targetloc
            os.system(cmd)   
    
        if os.path.isdir(respfileloc)==False:
            cmd='mkdir '+respfileloc
            os.system(cmd)
            
            
        for channel in channels:
            
            
            if stafile=="*":
                if os.path.exists(filename)==False:
                    reqstring='./FetchData '+vfetchdata+' -N '+network+ '-C '+channel+' -s '+t1+' -e '+t2+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -rd '+respfileloc
                    print(reqstring)
                    os.system(reqstring)
            
            else:
                fh=open(stafile, 'r')
                stations=fh.read().split('\n')
                for station in stations:
                    if station=='': continue
                    reqstring=exdir+'/FetchData '+vfetchdata+' -N '+network+ ' -S '+station + ' -C '+channel+' -s '+t1+' -e '+t2+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -rd '+respfileloc
                    print(reqstring)
                    os.system(reqstring)
    
        return
