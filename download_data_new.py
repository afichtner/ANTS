import TOOLS.read_xml as rxml
import TOOLS.renamer as renamer
from obspy import UTCDateTime
import obspy as obs
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
    
    #-Remove any previous temporary download directory=============================================
    if os.path.isdir('download_temp'):
        os.system('rm -rf download_temp')
    
    os.system('mkdir download_temp')
    

    #- read input file ============================================================================

    datainput=rxml.read_xml(inputfile)
    dat1=datainput[1]
    
    # Verbose?
    if dat1['verbose']=='1':
        v=True

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
    
    # format
    format=' --'+dat1['storage']['format']+' '
     
    # storage of the data
    targetloc=dat1['storage']['downloadloc']
    respfileloc=dat1['storage']['respfileloc']
    
    if os.path.isdir(targetloc)==False:
        cmd='mkdir '+targetloc
        os.system(cmd)   

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
    upd=0
    
    for network in networks:
        for channel in channels:
            #- download data when no station file is provided. Download all available stations ====

            if stafile=='*':

                identity=' --identity '+network+'.*.'+location+'.'+channel
                if upd==0:
                    cmd='obspyDMT --continuous --datapath download_temp --ic_no '+identity+t1+t2+region+centre+format
                    upd=1
                else:
                    cmd='obspyDMT ' +update+' download_temp '+identity
                os.system(cmd)
                
            #- download data when station file is provided ========================================

            else:
        
                # open the station file
                f=open(stafile, 'r')
                stations=f.read().split('\n')
                print stations        
                # loop through all the stations
                for station in stations: 
                    if station=='':
                        continue
                    identity = ' --identity ' + network + '.'+station+'.'+location + '.' + channel
                    
                    if upd==0:
                        cmd='obspyDMT --continuous --datapath download_temp  --ic_no '+identity+t1+t2+region+centre+format
                        upd=1
                    else:
                        cmd='obspyDMT ' +update+' download_temp '+identity
                    os.system(cmd)
                        
        
    #- Get a list of the nice data that was just downloaded====================================
    os.system('bash general/get_datalist.sh download_temp/*')
    
    #- Rename data ====================================
    fh=open('channellist.txt', 'r')
    filelist=fh.read().split('\n')
    
    for file in filelist:
        if file=='':
            continue
        stream=obs.read(file)
        renamer.rename_seismic_data(stream[0], targetloc, False, v)
    
    #-Move all the resp files to the respfileloc folder
    cmd='bash general/save_resp.sh download_temp/* '+respfileloc
    os.system(cmd)
    if v:
        print '===================================================================='
        print 'Moved Response files.'
    
    #- Remove temporary directory ====================================
    os.system('rm -rf download_temp')
    if v:
        print 'Removed temporary download folder.'
        print 'Done.'
