import TOOLS.read_xml as rxml
from obspy import UTCDateTime
import os
import sys
import shutil
import download_resp as dr
import antconfig as cfg


if __name__=='__main__':
    xmlin=str(sys.argv[1])
    dr.download_resp(xmlin)
    

def download_resp(xmlinput):
    
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
    stafile=dat1['ids']
    
    # geographical region
    lat_min=dat1['region']['lat_min']
    lat_max=dat1['region']['lat_max']
    lon_min=dat1['region']['lon_min']
    lon_max=dat1['region']['lon_max']
     
    

    respfileloc=datadir+'resp/'
        
    if os.path.isdir(respfileloc)==False:
        cmd='mkdir '+respfileloc
        os.system(cmd)
            
    fh=open(stafile, 'r')
    ids=fh.read().split('\n')
    
    
    for id in ids:
        if id=='': continue
        
        network=id.split('.')[0]
        station=id.split('.')[1]
        channel=id.split('.')[3]
        
        print '\n Downloading response information from: '+id+'\n'
        reqstring=exdir+'/FetchData '+vfetchdata+' -N '+network+ ' -S '+station+' -C '+channel+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -rd '+respfileloc
        os.system(reqstring)
    
    return
