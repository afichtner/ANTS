# A script to download ambient vibration records
from __future__ import print_function
import os
import sys
import shutil
from obspy import UTCDateTime
from math import ceil

from mpi4py import MPI
from glob import glob
import TOOLS.read_xml as rxml 
import antconfig as cfg


if __name__=='__main__':
    import par_download as pd
    
    pd.par_download()
    

def par_download():
    
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
       dat=rxml.read_xml(os.path.join(cfg.inpdir,'input_download.xml'))[1]
      
       # network, channel, location and station list
       stalist=os.path.join(cfg.inpdir,'downloadlist.txt')
       fh=open(stalist,'r')
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
        
    # Quality?
    quality = dat['quality']
        
    # time interval of request
    t1=dat['time']['starttime']
    t1str=UTCDateTime(t1).strftime('%Y.%j.%H.%M.%S')
    t2=dat['time']['endtime']
    t2str=UTCDateTime(t2).strftime('%Y.%j.%H.%M.%S')
 
    # data segment length
    if dat['time']['len']==None:
        winlen=UTCDateTime(t2)-UTCDateTime(t1)
    else:
        winlen = int(dat['time']['len'])
    
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
        
        t = UTCDateTime(t1)
        while t < UTCDateTime(t2):
            
            tstart = UTCDateTime(t).strftime('%Y-%m-%d,%H:%M:%S')
            tstartstr = UTCDateTime(t).strftime('%Y.%j.%H.%M.%S')
            
            tstep = min((UTCDateTime(t)+winlen),UTCDateTime(t2)).\
            strftime('%Y-%m-%d,%H:%M:%S')
            tstepstr = min((UTCDateTime(t)+winlen),UTCDateTime(t2)).\
            strftime('%Y.%j.%H.%M.%S')
            
            
            #-Formulate a polite request
            filename=targetloc+id+'.'+tstartstr+'.'+tstepstr+'.mseed'
            if os.path.exists(filename)==False:
                network=id.split('.')[0]
                station=id.split('.')[1]
                channel=id.split('.')[3]
                #print network, station, location, channel
                print('\n Rank '+str(rank)+'\n',file=None)
                print('\n Attempting to download data from: '+id+'\n',file=None)
                print(filename,None)
                
                reqstring=exdir+'/FetchData '+vfetchdata+' -N '+network+ \
                 ' -S '+station+' -C '+channel+' -s '+tstart+' -e '+tstep+ \
                 ' -msl '+minlen+' --lat '+lat_min+':'+lat_max+ \
                ' --lon '+lon_min+':'+lon_max+' -o '+filename+' -Q '+quality
                  
                os.system(reqstring)
            t += winlen
          
     
    # Clean up (some files come back with 0 data)
    stafile=dat['ids']
    t1s=t1str.split('.')[0]+'.'+t1str.split('.')[1]
    t2s=t2str.split('.')[0]+'.'+t2str.split('.')[1]
    
    cmd=('./UTIL/cleandir.sh '+targetloc)     
    os.system(cmd)
    os.system('mv '+targetloc+'* '+targetloc+'/..')
    os.system('rmdir '+targetloc)
    
    
    # Download resp files for all epochs!
    respfileloc=datadir+'resp/'
        
    if os.path.isdir(respfileloc)==False:
        cmd='mkdir '+respfileloc
        os.system(cmd)
    
    
    for id in myids:
        if id=='': continue
        
        network=id.split('.')[0]
        station=id.split('.')[1]
        channel=id.split('.')[3]
        
        print('\n Downloading response information from: '+id+'\n')
        reqstring=exdir+'/FetchData '+vfetchdata+' -N '+network+ ' -S '+station+' -C '+channel+\
        ' --lat '+lat_min+':'+lat_max+' --lon '+lon_min+':'+lon_max+' -rd '+respfileloc
        os.system(reqstring)
        
    comm.Barrier()

    
    if rank==0:
        outfile=os.path.join(cfg.datadir,'raw/latest/download_report.txt')
        outf=open(outfile,'w')
        
        print('Attempting to download data from stations: \n',file=outf)
        print('****************************************** \n',file=outf)
        for id in ids:
            print(id,file=outf)
        print('****************************************** \n',file=outf)
        stalist=os.path.join(cfg.inpdir,'downloadlist.txt')
        fh=open(stalist,'r')
        ids=fh.read().split('\n')
        
        noreturn=[]
        
        for id in ids:
            if id=='': continue
            fls=glob(os.path.join(cfg.datadir,'raw/latest',id+'*'))
            if fls != []:
                print('Files downloaded for id: '+id,file=outf)
                print('First file: '+fls[0],file=outf)
                print('Last file: '+fls[-1],file=outf)
                print('****************************************** \n',file=outf)    
            else: 
                noreturn.append(id)
            
        if noreturn != []:
            print('NO files downloaded for: \n',file=outf)
            
            print(noreturn,file=outf)
     
        print('****************************************** \n',file=outf)
        print('Download parameters were: \n',file=outf)
        print('****************************************** \n',file=outf)
        outf.close()
        
        os.system('cat '+cfg.inpdir+'/input_download.xml >> '+outfile)
        