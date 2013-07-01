import TOOLS.read_xml as rxml
from obspy import UTCDateTime
import sys
import os


def download_data(inputfile):

    """
    
    Tool for the download continuous seismic data from a collection of stations and/or networks.

    The download is based on the obspyDMT and takes as input an xml file that specifies the download parameters.

    """
        
    #==============================================================================================
    #- initialisations
    #==============================================================================================

    #- read input file ============================================================================

    datainput=rxml.read_xml(inputfile)
    dat1=datainput[1]

    # network, channel, location and station list
    networks=dat1['data']['networks'].strip().split(' ')
    channels=dat1['data']['channels'].strip().split(' ')
    location=dat1['data']['location']
    stafile=dat1['data']['stations']

    print networks
    print channels

    # time interval of request
    t1=' --min_date '+dat1['time']['starttime']
    t2=' --max_date '+dat1['time']['endtime']

    # geographical region
    lat_min=dat1['region']['lat_min']
    lat_max=dat1['region']['lat_max']
    lon_min=dat1['region']['lon_min']
    lon_max=dat1['region']['lon_max']
    region=' --station_rect '+lon_min+'/'+lon_max+'/'+lat_min+'/'+lat_max

    # storage of the data
    downloadloc=' --datapath '+dat1['storage']['downloadloc']
    format=' --'+dat1['storage']['format']+' '

    #- specify data centre (exclude ArcLink when centre='iris') ====================================

    if dat1['centre']=='iris':
        centre=" --arc 'N' "
        update=' --iris_update '
    else:
        centre=''
        update=' --update_all '
        
    #==============================================================================================
    #- downloads
    #==============================================================================================

    # loop over networks and channels
    for network in networks:
        for channel in channels:

            #- download data when no station file is provided. Download all available stations ====

            if stafile=='*':

                identity=' --identity '+network+'.*.'+location+'.'+channel 

                # download when directory already exists
                if os.path.isdir(dat1['storage']['downloadloc']):
                    cmd='obspyDMT --continuous --ic_no '+update+dat1['storage']['downloadloc']+identity+t1+t2+region+format
                # create new directory
                else:
                    cmd='obspyDMT --continuous --ic_no '+identity+t1+t2+region+centre+format+downloadloc
            
                os.system(cmd)
       
            #- download data when station file is provided ========================================

            else:
        
                # open the station file
                f=open(stafile, 'r')
                stations=f.read().split('\n')
        
                # loop through all the stations
                for station in stations:           
                    identity = ' --identity ' + network + '.'+station+'.'+location + '.' + channel 

                    # download when directory already exists
                    if os.path.isdir(dat1['storage']['downloadloc']):
                        cmd='obspyDMT --continuous --ic_no '+update+dat1['storage']['downloadloc']+identity+t1+t2+region+format
                    # create new directory
                    else:
                        cmd='obspyDMT --continuous --ic_no '+identity+t1+t2+region+centre+format+downloadloc
                
                    os.system(cmd)        
