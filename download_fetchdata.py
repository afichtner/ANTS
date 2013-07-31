import TOOLS.read_xml as rxml
from obspy import UTCDateTime
import os

def download_fetchdata(xmlinput):
    
    """
    
    Tool for the download continuous seismic data from a collection of stations and/or networks.

    The download is based on the iris DMC FetchData script and takes as input an xml file that specifies the download parameters.

    """
    
     #- read input file ============================================================================

    datainput=rxml.read_xml(xmlinput)
    dat1=datainput[1]
    
    # Verbose?
    if dat1['verbose']=='1':
        v=True
        vfetchdata='-v '
    else:
        vfetchdata=''

    # network, channel, location and station list
    networks=dat1['data']['networks'].strip().split(' ')
    channels=dat1['data']['channels'].strip().split(' ')
    locations=dat1['data']['locations'].strip().split(' ')
    stafile=dat1['data']['stations']

    # time interval of request
    t1=dat1['time']['starttime']
    t1str=UTCDateTime(t1).strftime('%Y%m%d%H%M%S')
    t2=dat1['time']['endtime']
    t2str=UTCDateTime(t2).strftime('%Y%m%d%H%M%S')

    # geographical region
    lat_min=dat1['region']['lat_min']
    lat_max=dat1['region']['lat_max']
    lon_min=dat1['region']['lon_min']
    lon_max=dat1['region']['lon_max']
     
    # storage of the data
    targetloc=dat1['storage']['downloadloc']
    respfileloc=dat1['storage']['respfileloc']
    
    if os.path.isdir(targetloc)==False:
        cmd='mkdir '+targetloc
        os.system(cmd)   
    
    if os.path.isdir(respfileloc)==False:
        cmd='mkdir '+respfileloc
        os.system(cmd)
    
    
    for network in networks:
        for channel in channels:
            for location in locations:

                if stafile=="*":
                    filename=targetloc+'/'+network+'.'+location+'.'+channel+'.'+t1str+'.'+t2str+'.mseed'
                    reqstring='./FetchData '+vfetchdata+' -N '+network+' -C '+channel+' -s '+t1+' -e '+t2+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -o '+filename+' -rd '+respfileloc
                    print reqstring
                    os.system(reqstring)

                else:
                    fh=open(stafile, 'r')
                    stations=fh.read().split('\n')
                    for station in stations:
                        if station=='': continue
                        print network + station + location + channel
                        #-Formulate a polite request
                        filename=targetloc+'/'+network+'.'+station+'.'+location+'.'+channel+'.'+t1str+'.'+t2str+'.mseed'
                        reqstring='./FetchData '+vfetchdata+' -N '+network+' -S '+station + ' -C '+channel+' -s '+t1+' -e '+t2+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -o '+filename+' -rd '+respfileloc
                        print reqstring
                        os.system(reqstring)