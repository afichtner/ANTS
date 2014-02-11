import TOOLS.read_xml as rxml
from obspy import UTCDateTime
import os
import sys
from shutil import copy

from mpi4py import MPI
from math import ceil

if __name__=='__main__':
    import par_download as pd
    xmlin=str(sys.argv[1])
    pd.download_fetchdata(xmlin)
    

def download_fetchdata(xmlinput):
    
    """
    
    Tool for the download continuous seismic data from a collection of stations and/or networks.

    The download is based on the iris DMC FetchData script and takes as input an xml file that specifies the download parameters.

    """
    #==============================================================================================
    # preliminaries
    #==============================================================================================
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size=comm.Get_size()

    #==============================================================================================
    #- MASTER process:
    #- reads in xmlinput
    #- creates output directory
    #- creates a list of input files
    #==============================================================================================
    if rank==0:
        dat1=rxml.read_xml(xmlinput)[1]
        network=dat1['data']['networks'].strip().split(' ')[0]
        # storage of the data
        targetloc='./DATA/raw/'+network
        respfileloc='./DATA/resp/'
        
        if os.path.isdir(targetloc)==False:
            cmd='mkdir '+targetloc
            os.system(cmd)   
        
        if os.path.isdir(respfileloc)==False:
            cmd='mkdir '+respfileloc
            os.system(cmd)
   
        
    #- Read the input station list==================================================================
        stafile=dat1['data']['stations'].strip().split(' ')[0]
        if stafile=="*":
            print 'Wildcarding stations is not allowed for parallel download. Provide a station list file.'
            return
        else:
            fh=open(stafile, 'r')
            stations=fh.read().split('\n')
            
        
    elif rank!=0:
        dat1=list()
        stations=list()  


    #==============================================================================================
    #- broadcast the input; and the station list
    stations=comm.bcast(stations, root=0)
    dat1=comm.bcast(dat1, root=0)
    #==============================================================================================
    
    
    # Verbose?
    if dat1['verbose']=='1':
        v=True
        vfetchdata='-v '
    else:
        vfetchdata=''
    # Directory where executable is located
    exdir=dat1['exdir']
    
    
    # network, channel, location and station list
    network=dat1['data']['networks'].strip().split(' ')[0]
    channels=dat1['data']['channels'].strip().split(' ')
    locations=dat1['data']['locations'].strip().split(' ')
    
    # Folders
    targetloc='./DATA/raw/'+network
    respfileloc='./DATA/resp/'
    
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
     
    #==============================================================================================
    #- Assign each rank its own chunk of stations
    #==============================================================================================
    clen=int(ceil(float(len(stations))/float(size)))
    chunk=(rank*clen, (rank+1)*clen)
    mystations=stations[chunk[0]:chunk[1]]
    if v: 
        print 'Hi I am process nr ', rank, 'and I request data from these stations:'
        print mystations
    
    for channel in channels:
        for location in locations:
            for station in mystations:
                if station=='': continue
                #-Formulate a polite request
                filename=targetloc+'/'+network+'.'+station+'.'+location+'.'+channel+'.'+t1str+'.'+t2str+'.mseed'
                reqstring=exdir+'/FetchData '+vfetchdata+' -N '+network+' -S '+station+' -L '+location + ' -C '+channel+' -s '+t1+' -e '+t2+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -o '+filename+' -rd '+respfileloc
                if v: print reqstring
                os.system(reqstring)
