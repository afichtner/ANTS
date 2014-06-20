# A script to download ambient vibration records
import os
import sys
import shutil
from obspy import UTCDateTime
from math import ceil

from mpi4py import MPI

import TOOLS.read_xml as rxml 
import antconfig as cfg


if __name__=='__main__':
    import par_download as pd
    xmlin=str(sys.argv[1])
    print 'XML input file: '+ xmlin
    pd.par_download(xmlin)


def par_download(xmlinput):
    
    """
    
    Parallel download from IRIS DMC
    
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
    
       datadir=cfg.datadir
       dat=rxml.read_xml(xmlinput)[1]
      
       # network, channel, location and station list
       stafile=dat['ids']
       fh=open(stafile, 'r')
       ids=fh.read().split('\n')
       
    #==============================================================================================
    #- All processes:
    #- receive the input; and the list of files
    #- read variables from broadcasted input
    #==============================================================================================
    
    else:
        ids=list()
        dat=list()    
       
    ids=comm.bcast(ids, root=0)
    dat=comm.bcast(dat, root=0)
    
    datadir=cfg.datadir
    targetloc=datadir+'raw/latest/rank'+str(rank)+'/'
    
    if os.path.isdir(targetloc)==False:
        cmd='mkdir '+targetloc
        os.system(cmd)
    
    
    # Directory where executable is located
    exdir=dat['exdir']
    
    # Verbose?
    if dat['verbose']=='1':
        v=True
        vfetchdata='-v '
    else:
        vfetchdata=''
        
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
    
    #==============================================================================================
    #- Assign each rank its own chunk of input
    #==============================================================================================

    clen=int(float(len(ids))/float(size))
    chunk=(rank*clen, (rank+1)*clen)
    myids=ids[chunk[0]:chunk[1]]
    if rank==size-1:
        myids=ids[chunk[0]:]
    
    #==================================================================================
    # Input files loop
    #==================================================================================
      
    for id in myids:
        
        if id=='': continue
        
        #-Formulate a polite request
        filename=targetloc+id+'.'+t1str+'.'+t2str+'.mseed'
        if os.path.exists(filename)==False:
            network=id.split('.')[0]
            station=id.split('.')[1]
            channel=id.split('.')[3]
            #print network, station, location, channel
            print '\n Rank '+str(rank)+'\n'
            print '\n Attempting to download data from: '+id+'\n'
            reqstring=exdir+'/FetchData '+vfetchdata+' -N '+network+ ' -S '+station+' -C '+channel+' -s '+t1+' -e '+t2+' -msl '+minlen+' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -o '+filename
            os.system(reqstring)
      
     
    # Clean up (some files come back with 0 data)
    stafile=dat['ids']
    t1s=t1str.split('.')[0]+'.'+t1str.split('.')[1]
    t2s=t2str.split('.')[0]+'.'+t2str.split('.')[1]
    cleanupinfo=targetloc+stafile.split('/')[-1].split('.')[0]+'.'+t1s+'.'+t2s
    
    cmd=('./UTIL/cleandir.sh '+targetloc+' '+cleanupinfo)      
    os.system(cmd)
      #return