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
import TOOLS.normalisation as nrm
import TOOLS.read_xml as rxml 
import TOOLS.mergetraces as mt
import antconfig as cfg

if __name__=='__main__':
    import par_ic as pic
    xmlin=str(sys.argv[1])
    print('XML input file: '+ xmlin,file=None)
    pic.ic(xmlin)


def ic(xmlinput,content=None):
    
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
       inp1=rxml.read_xml(xmlinput)[1]
       
       verbose=bool(int(inp1['verbose']))
       update=bool(int(inp1['update']))
       check=bool(int(inp1['check']))
       prepname=inp1['prepname']
       #startyr=int(inp1['input']['startyr'][0:4])
       #endyr=int(inp1['input']['endyr'][0:4])
      
       
       #- copy the input xml to the output directory for documentation ===============================
       xmlinname=datadir+'/processed/xmlinput/ic.'+prepname+'.xml'
       
       if os.path.exists(xmlinname)==True and update == False:
           print('Name tag already in use! New generic name tag chosen. Please review tag later to avoid overwriting.',file=None)
           prepname = UTCDateTime().strftime('proc%Y-%j')
           xmlinname=datadir+'/processed/xmlinput/ic.'+prepname+'.xml'
           print('New tag is '+prepname,file=None)
       
       if update == False:
           shutil.copy(xmlinput,xmlinname)
       
       
       #for i in range(startyr-1,endyr+1):
    #       if os.path.exists(datadir+'/processed/'+str(i)+'/')==False:
    #       os.mkdir(datadir+'/processed/'+str(i))
       if os.path.exists(datadir+'processed/'+prepname)==False:
           os.mkdir(datadir+'processed/'+prepname)
       
       #- check what input is, list input from different directories =================================
       if content==None:
           indirs=inp1['input']['indirs'].strip().split(' ')
           content=list()
           for indir in indirs:
               content.extend(glob(indir+'/*'))
           
       elif type(content)==str: 
           filename=content
           content=list()
           content.append(filename)
       
       content.sort()
           
       #- If only a check run is performed, then only a couple of files are preprocessed
       if check and len(content)>4:
           content=[content[0],content[1],content[len(content)-2],content[len(content)-1]]
           
   
    #==============================================================================================
    #- All processes:
    #- receive the input; and the list of files
    #- read variables from broadcasted input
    #==============================================================================================
    
    else:
        content=list()
        inp1=list()
        prepname=''
       
    t1=time.time()-t0
    content=comm.bcast(content, root=0)
    inp1=comm.bcast(inp1, root=0)
    prepname=comm.bcast(prepname, root=0)
    t2=time.time()-t0-t1
    
    verbose=bool(int(inp1['verbose']))
    update=bool(int(inp1['update']))
    datadir=cfg.datadir
    if update ==True:
        ofid=open(datadir+'/processed/out/update.'+prepname+'.rank_'+str(rank)+'.txt','w')
    else:
        ofid=open(datadir+'/processed/out/proc.'+prepname+'.rank_'+str(rank)+'.txt','w')
        
    respdir=inp1['processing']['instrument_response']['respdir']
    unit=inp1['processing']['instrument_response']['unit']
    freqs=inp1['processing']['instrument_response']['freqs']
    wl=inp1['processing']['instrument_response']['waterlevel']
    seglen=float(inp1['processing']['split']['length_in_sec'])
    minlen=float(inp1['quality']['min_length_in_sec'])
    mergegap=float(inp1['quality']['maxgaplen'])
    Fs_original=inp1['processing']['decimation']['Fs_old'].split(' ')
    Fs_old=list()
    for fs in Fs_original:
        Fs_old.append(float(fs))
    Fs_down=inp1['processing']['decimation']['Fs_new'].split(' ')
    Fs_new=list()
    for fs in Fs_down:
        Fs_new.append(float(fs))
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
            print('Time at start was '+UTCDateTime().strftime('%Y.%m.%d, %H:%M')+' GMT\n',file=None)
            print('Rank 0 took '+str(t1)+' seconds to read in input\n',file=None)
            print('Broadcasting took '+str(t2)+' seconds \n',file=None)
            
    print('\nHi I am rank number %d and I am processing the following files for you: \n' %rank,file=ofid)
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
            print('===========================================================',file=ofid)
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
    
        #- clean the data merging segments with less than a specified number of seconds:
        data=mt.mergetraces(data,Fs_old,mergegap)
        data.split()
        
        #- initialize stream to 'recollect' the split traces
        colloc_data=Stream()
        
      
        #- split traces into shorter segments======================================================
        if inp1['processing']['split']['doit']=='1':
            data=proc.slice_traces(data,seglen,minlen,verbose,ofid)
        n_traces=len(data)
        if verbose==True:
            print('* contains '+str(n_traces)+' trace(s)',file=ofid)
            
        #- trim ===============================================================================
        
        if inp1['processing']['trim']=='1':
            data=proc.trim_next_sec(data,verbose,ofid)
        
        
        #==================================================================================
        # trace loop
        #==================================================================================
        for k in np.arange(n_traces):
            trace=data[k]
            
            if update == True:
                if len(glob(getfilepath(trace.stats,prepname,True))) > 0:
                    print('File already processed, proceeding...',file=ofid)
                    break
                else:
                    print('Updating...',file=ofid)
            
            if verbose==True: print('-----------------------------------------------------------',file=ofid)
    
            #==================================================================================
            # basic quality checks
            #==================================================================================
    
            #- check NaN
            if True in np.isnan(trace.data):
                if verbose==True: print('** trace contains NaN, discarded',file=ofid)
                continue
    
            #- check infinity
            if True in np.isinf(trace.data):
                if verbose==True: print('** trace contains infinity, discarded',file=ofid)
                continue
    
            if verbose==True: print('* number of points: '+str(trace.stats.npts)+'\n',file=ofid)
    
            #==================================================================================
            # processing (detrending, filtering, response removal, decimation)
            #==================================================================================
                              
            #- demean============================================================================
            if inp1['processing']['detrend']=='1':
    
                trace=proc.detrend(trace,verbose,ofid)
                
            if inp1['processing']['demean']=='1':
    
                trace=proc.demean(trace,verbose,ofid)
            
    
            #- taper edges ========================================================================
    
            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose,ofid)
                
          
            
            #- downsampling =======================================================================
            k=0
            while k<len(Fs_new):
                if trace.stats.sampling_rate>Fs_new[k]:
                    trace=proc.downsample(trace,Fs_new[k],verbose,ofid)
                k+=1
            newtrace = trace.copy()
            del trace
               
               
            #- remove instrument response =========================================================
    
            if inp1['processing']['instrument_response']['doit']=='1':
    
                removed,newtrace=proc.remove_response(newtrace,respdir,unit,freqs,wl,verbose,ofid)
                if removed==False:
                    print('** Instrument response could not be removed! Trace discarded.',file=ofid)
                    continue
                    
                if True in np.isnan(newtrace):
                    print('** Deconvolution seems unstable! Trace discarded.',file=ofid)
                    continue
          
            #- merge all into final trace =========================================================
            colloc_data+=newtrace
             
            #- flush buffer of output file ========================================================
            ofid.flush()
            
            del newtrace

        colloc_data=mt.mergetraces(colloc_data,Fs_new,mergegap,ofid)
        colloc_data._cleanup()

        for k in range(len(colloc_data)):
            if ((inp1['processing']['instrument_response']['doit']=='1') and (removed==1)) or \
                inp1['processing']['instrument_response']['doit']!='1':
                
                filepathnew = getfilepath(mydir,colloc_data[k].stats,prepname)
                
                #- write to file
                colloc_data[k].write(filepathnew,format=colloc_data[k].stats._format)
                       
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
        filepathnew=mydir+'/'+network+'.'+station+'.'+location+'.'+channel+'.' + t1 + '.' +t2+'.'+prepname+'.'+format 
    else:
        filepathnew=mydir+'/'+network+'.'+station+'.'+location+'.'+channel+'.' + t1 + '.*.'+prepname+'.'+format
    return filepathnew
        
    
