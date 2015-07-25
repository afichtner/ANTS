# A script to process ambient vibration records
from __future__ import print_function
# Use the print function to be able to switch easily between stdout and a file
import os
import sys
import shutil
import time

from math import ceil
from obspy import read, Stream,  Trace, UTCDateTime
from obspy.signal import filter
from mpi4py import MPI
from glob import glob

import matplotlib.pyplot as plt
import numpy as np

import TOOLS.processing as proc
import TOOLS.read_xml as rxml 
import TOOLS.mergetraces as mt

import antconfig as cfg
import INPUT.input_correction as inp

if __name__=='__main__':
    import par_ic as pic
    pic.ic()


def ic(content=None):
    
    """
    
    This script preprocesses the MSEED files at the path specified as command line argument 2.
    Command line argument 1 must be xml input file.
    
    """

    #==============================================================================================
    # preliminaries
    #==============================================================================================
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size=comm.Get_size()
    t0=time.time()
   
    #==============================================================================================
    #- MASTER process:
    #- reads in xmlinput
    #- creates output directory
    #- creates a list of input files
    #==============================================================================================
    
    if rank==0:
    
       datadir=cfg.datadir
       verbose=inp.verbose
       update=inp.update
       check=inp.check
       prepname=inp.prepname
       #startyr=int(inp1['input']['startyr'][0:4])
       #endyr=int(inp1['input']['endyr'][0:4])
      
       
       #- copy the input xml to the output directory for documentation ===============================
       inname=datadir+'/processed/input/ic.'+prepname+'.txt'
       
       if os.path.exists(inname)==True and update == False:
           print('Name tag already in use! New generic name tag chosen. \
           Please review tag later to avoid overwriting.',file=None)
           prepname = UTCDateTime().strftime('proc%Y-%j')
           inname=datadir+'/processed/input/ic.'+prepname+'.txt'
           print('New tag is '+prepname,file=None)    
       
       if update == False:
           shutil.copy(os.path.join(cfg.inpdir,'input_correction.py'),inname)
       
       
       #for i in range(startyr-1,endyr+1):
    #       if os.path.exists(datadir+'/processed/'+str(i)+'/')==False:
    #       os.mkdir(datadir+'/processed/'+str(i))
       if os.path.exists(datadir+'processed/'+prepname)==False:
           os.mkdir(datadir+'processed/'+prepname)
       
       #- check what input is, list input from different directories =================================
       if content==None:
           indirs=inp.indirs
           content=list()
           for indir in indirs:
               print(indir)
               content.extend(glob(indir+'/*'))
           
       elif type(content)==str: 
           filename=content
           content=list()
           content.append(filename)
       
       content.sort()
           
       #- If only a check run is performed, then only a couple of files are preprocessed
       if check==True and len(content)>4:
           content=[content[0],content[1],content[len(content)-2],\
           content[len(content)-1]]
           
   
    #==============================================================================================
    #- All processes:
    #- receive the input; and the list of files
    #- read variables from broadcasted input
    #==============================================================================================
    
    else:
        content=list()
        prepname=''
       
    t1=time.time()-t0
    content=comm.bcast(content, root=0)
    prepname=comm.bcast(prepname, root=0)
    t2=time.time()-t0-t1
    
    verbose=inp.verbose
    update=inp.update
    datadir=cfg.datadir
    
    if update ==True:
        ofid=open(datadir+'/processed/out/update.'+prepname+'.rank_'+\
        str(rank)+'.txt','w')
    else:
        ofid=open(datadir+'/processed/out/proc.'+prepname+'.rank_'+str(rank)+\
        '.txt','w')
    
    check=inp.check
    respdir=inp.respdir
    unit=inp.unit
    freqs=inp.freqs
    wl=inp.waterlevel
    seglen=inp.length_in_sec
    minlen=inp.min_length_in_sec
    mergegap=inp.maxgaplen
    Fs_original=inp.Fs_old
    Fs_new=inp.Fs_new
    Fs_new.sort() # Now in ascending order
    Fs_new=Fs_new[::-1] # Now in descending order
    
    #==============================================================================================
    #- Assign each rank its own chunk of input
    #==============================================================================================
    
    nfiles = int(len(content) / size)
    restfiles = len(content) % size
    
    mycontent=content[rank * nfiles : (rank + 1) * nfiles]
    if rank < restfiles:
        mycontent.append(content[size * nfiles + rank])
    del content
    t3=time.time()-t0-t2
    
    
    #- Print some nice comments to output file ---------------------------------------- 
    if rank == 0:      
        if verbose:
            print('Time at start was '+UTCDateTime().strftime('%Y.%m.%d, %H:%M')+\
            ' GMT\n',file=None)
            print('Rank 0 took '+str(t1)+' seconds to read in input\n',file=None)
            print('Broadcasting took '+str(t2)+' seconds \n',file=None)
            
    print('\nHi I am rank number %d and I am processing the following files for you:\
           \n' %rank,file=ofid)
    for fname in mycontent:
        ofid.write(fname+'\n')
    
    #==============================================================================================
    #- Input file loop
    #==============================================================================================
    mydir=datadir+'processed/'+prepname+'/rank'+str(rank)
    if os.path.exists(mydir)==False:
        os.mkdir(mydir)
        
    
    for filepath in mycontent:
        if verbose==True:
            print('===========================================================',\
            file=ofid)
            print('* opening file: '+filepath+'\n',file=ofid)
            
        #- read data
        try:
            data=read(filepath)
            
        except (TypeError, IOError):
            if verbose==True: print('** file wrong type or not found, skip.',file=ofid)
            continue
        except:
            if verbose: print('** unexpected read error, skip.',file=ofid)
            continue    
        #- check if this file contains data
        if len(data) == 0:
            print('File contains no data!',file=None)
            print('File contains no data!',file=ofid)
            continue
        
        #- clean the data merging segments with gap shorter than a specified number of seconds:
        print(data[0].stats.starttime)
        data=mt.mergetraces(data,Fs_original,mergegap)
        print(data[0].stats.starttime)
        data.split()
        
        #- initialize stream to 'recollect' the split traces
        colloc_data=Stream()
        
      
        #- split traces into shorter segments======================================================
        if inp.split_do == True:
            print('slice')
            print(data[0].stats.starttime)
            data=proc.slice_traces(data,seglen,minlen,verbose,ofid)
            print(data[0].stats.starttime)
        n_traces=len(data)
            
        if verbose==True:
            print('* contains '+str(n_traces)+' trace(s)',file=ofid)
            
        #- trim ===============================================================================
        
        if inp.trim == True:
            data=proc.trim_next_sec(data,verbose,ofid)
        
        
        #==================================================================================
        # trace loop
        #==================================================================================
        for k in np.arange(n_traces):
            trace=data[k]
            
            if check==True:
                ctr=trace.copy()
                ctr.stats.network='Original Data'
                ctr.stats.station=''
                ctr.stats.location=''
                ctr.stats.channel=''
                cstr=Stream(ctr)
            
            if update == True:
                if len(glob(getfilepath(mydir,trace.stats,prepname,True))) > 0:
                    print('File already processed, proceeding...',file=ofid)
                    break
                else:
                    print('Updating...',file=ofid)
            
            if verbose==True: print('-----------------------------------------\
            ------------------',file=ofid)
    
            #==================================================================================
            # basic quality checks
            #==================================================================================
    
            #- check NaN
            if True in np.isnan(trace.data):
                if verbose==True: print('** trace contains NaN, discarded',\
                file=ofid)
                continue
    
            #- check infinity
            if True in np.isinf(trace.data):
                if verbose==True: print('** trace contains infinity, discarded',\
                file=ofid)
                continue
    
            if verbose==True: print('* number of points: '+str(trace.stats.npts)+\
            '\n',file=ofid)
    
            #==================================================================================
            # processing (detrending, filtering, response removal, decimation)
            #==================================================================================
                              
            #- demean============================================================================
            if inp.detrend == True:
    
                trace=proc.detrend(trace,verbose,ofid)
                
            if inp.demean == True:
    
                trace=proc.demean(trace,verbose,ofid)
            
    
            #- taper edges ========================================================================
    
            if inp.taper_do == True:
    
                trace=proc.taper(trace,inp.taper_width,verbose,ofid)
                
          
            
            #- downsampling =======================================================================
            k=0
            while k<len(Fs_new):
                if trace.stats.sampling_rate>Fs_new[k]:
                    print('Decimate')
                    print(data[0].stats.starttime)
                    trace=proc.downsample(trace,Fs_new[k],verbose,ofid)
                    print(data[0].stats.starttime)
                k+=1
            newtrace = trace.copy()
            del trace
               
               
            #- remove instrument response =========================================================
    
            if inp.remove_response == True:
                print('IC')
                print(data[0].stats.starttime)
                removed,newtrace=proc.remove_response(newtrace,respdir,unit,\
                freqs,wl,verbose,ofid)
                print(data[0].stats.starttime)
                if removed==False:
                    print('** Instrument response could not be removed! \
                    Trace discarded.',file=ofid)
                    continue
                    
                if True in np.isnan(newtrace):
                    print('** Deconvolution seems unstable! Trace discarded.',\
                    file=ofid)
                    continue
        
            if check==True:
                ctr = newtrace.copy()
                ctr.stats.network='After IC to '+unit
                ctr.stats.station=''
                ctr.stats.location=''
                ctr.stats.channel=''
                cstr.append(ctr)
                cstr.plot(outfile=datadir+'/processed/out/'+\
                filepath.split('/')[-1]+'.'+prepname+'.png',equal_scale=False)
                cstr.trim(endtime=cstr[0].stats.starttime+3600)
                cstr.plot(outfile=datadir+'/processed/out/'+\
                filepath.split('/')[-1]+'.'+prepname+'.1hr.png',equal_scale=False)
                
            #- merge all into final trace =========================================================
            colloc_data+=newtrace
             
            #- flush buffer of output file ========================================================
            ofid.flush()
            
            del newtrace
        
        if len(colloc_data) == 0: 
            print('*** NO data returned from this file: '+filepath.split('/')[-1])
            continue
        colloc_data=mt.mergetraces(colloc_data,Fs_new,mergegap,ofid)
        colloc_data._cleanup()

        for k in range(len(colloc_data)):
            if ((inp.remove_response==True) and \
            (removed==1)) or \
                inp.remove_response==False:
                
                filepathnew = getfilepath(mydir,colloc_data[k].stats,prepname)
                
                #- write to file
                colloc_data[k].write(filepathnew,\
                format=colloc_data[k].stats._format)
                       
                if verbose==True:
                    print('* renamed file: '+filepathnew,file=ofid)
        
        del colloc_data
        del data
        
    if ofid:
        print("Rank %g has completed processing." %rank,file=None)
        ofid.close()
    os.system('mv '+mydir+'/* '+mydir+'/../')
    os.system('rmdir '+mydir)    
        
def getfilepath(mydir,stats,prepname,startonly=False):
    
    network=stats.network
    station=stats.station
    location=stats.location
    channel=stats.channel
    format=stats._format
        
    t1=stats.starttime.strftime('%Y.%j.%H.%M.%S')
    t2=stats.endtime.strftime('%Y.%j.%H.%M.%S')
    yr=str(t1[0:4])
    
    if startonly == False:
        filepathnew=mydir+'/'+network+'.'+station+'.'+location+'.'+\
        channel+'.' + t1 + '.' +t2+'.'+prepname+'.'+format 
    else:
        filepathnew=mydir+'/'+network+'.'+station+'.'+location+'.'+\
        channel+'.' + t1 + '..'+prepname+'.'+format
    return filepathnew
        
    
