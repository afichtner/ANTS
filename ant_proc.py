# A script to process ambient vibration records
from __future__ import print_function
# Use the print function to be able to switch easily between stdout and a file
from mpi4py import MPI
import os
import sys
import shutil
import time

from math import ceil
from obspy import read, Stream,  Trace, UTCDateTime
from glob import glob

import matplotlib.pyplot as plt
import numpy as np

import TOOLS.processing as proc
import TOOLS.read_xml as rxml 
import TOOLS.mergetraces as mt
import TOOLS.event_excluder as ee

import antconfig as cfg
import INPUT.input_correction as inp

if __name__=='__main__':
    import ant_proc as pp

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    inname=cfg.datadir+'/processed/input/ic.'+inp.prepname+'.txt'
    if rank==0 and os.path.exists(inname)==True and inp.update == False:
        sys.exit("Choose a new name tag or set update to True. Aborting")
        MPI.COMM_WORLD.Abort(1)
    pp.ic(rank,size)


def ic(rank,size):
    
    """
    
    This script preprocesses the MSEED files at the path specified as command line argument 2.
    Command line argument 1 must be xml input file.
    
    """
    datadir=cfg.datadir
    verbose=inp.verbose
    update=inp.update
    check=inp.check
    prepname=inp.prepname
    datadir=cfg.datadir
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
     
    try:
        os.mkdir(datadir+'processed/'+prepname)
    except OSError:
        pass
    
     
    #- copy the input xml to the output directory for documentation =========== 
    if rank==0:
        inname=datadir+'/processed/input/ic.'+prepname+'.txt'
        if update == False:
            shutil.copy(os.path.join(cfg.inpdir,'input_correction.py'),inname)
       
        #- check what input is, list input from different directories =================================
    
    indirs=inp.indirs
    content=list()
    for indir in indirs:
        print(indir)
        content.extend(glob(indir+'/*'))
        
    content.sort()
           
       #- If only a check run is performed, then only a couple of files are preprocessed
    if check==True and len(content)>4:
        content=[content[0],content[1],content[len(content)-2],\
        content[len(content)-1]]
    
    if update ==True:
        ofid=open(datadir+'/processed/out/update.'+prepname+'.rank_'+\
        str(rank)+'.txt','w')
    else:
        ofid=open(datadir+'/processed/out/proc.'+prepname+'.rank_'+str(rank)+\
        '.txt','w')
     #==============================================================================================
    #- Assign each rank its own chunk of input
    #==============================================================================================
    
    nfiles = int(len(content) / size)
    restfiles = len(content) % size
    
    mycontent=content[rank * nfiles : (rank + 1) * nfiles]
    if rank < restfiles:
        mycontent.append(content[size * nfiles + rank])
    del content
    
    
    
    #- Print some nice comments to output file ---------------------------------------- ------ 
    print('\nHi I am rank number %d and I am processing the following files for you:\
           \n' %rank,file=ofid)
    for fname in mycontent:
        ofid.write(fname+'\n')
    
    if check==True and inp.debugfile is not None:
        dfile=open(inp.debugfile,'w') #==============================================================================================
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
        data=mt.mergetraces(data,Fs_original,mergegap)
        data.split()
        
        #- initialize stream to 'recollect' the split traces
        colloc_data=Stream()
        
      
        #- split traces into shorter segments======================================================
        if inp.split_do == True:
            data=proc.slice_traces(data,seglen,minlen,verbose,ofid)
        n_traces=len(data)
        if verbose==True:
            print('* contains '+str(n_traces)+' trace(s)',file=ofid)
            
        #- trim ===============================================================================
        
        if inp.trim == True:
            data=proc.trim_next_sec(data,verbose,ofid)
        
        
        #==================================================================================
        # trace loop
        #==================================================================================
        for trace_index in np.arange(n_traces):
            
            trace=data[trace_index]
            if trace.stats.npts / inp.Fs_new[-1] < minlen:
                continue
            if trace.stats.npts / inp.Fs_new[-1] < 39:
                continue    
                
            if check==True:
                ctr=trace.copy()
                ctr.stats.network='Original Data'
                ctr.stats.station=''
                ctr.stats.location=''
                ctr.stats.channel=''
                cstr=Stream(ctr)
                print(trace,file=dfile)
                dfile.write('-----------------------------------------------\n')
                dfile.write('Original\n')
                print(trace.data[0:20],file=dfile)
                dfile.write('\n')
            
            
            if update == True:
                if len(glob(getfilepath(mydir,trace.stats,prepname,True))) > 0:
                    print('File already processed, proceeding...',file=ofid)
                    print(trace)
                    print('File already processed, proceeding...',file=None)
                    
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
                if verbose: print('** trace contains infinity, discarded',\
                file=ofid)
                continue
    
            if verbose: print('* number of points: '+str(trace.stats.npts)+\
            '\n',file=ofid)
    
            #==================================================================================
            # processing (detrending, filtering, response removal, decimation)
            #==================================================================================
                              
            #- demean============================================================================
            if inp.detrend:
    
                trace=proc.detrend(trace,verbose,ofid)
                
                if check:
                    dfile.write('Detrended\n')
                    print(trace.data[0:20],file=dfile)
                    dfile.write('\n')
                
            if inp.demean:
    
                trace=proc.demean(trace,verbose,ofid)
                
                if check:
                    dfile.write('Mean removed\n')
                    print(trace.data[0:20],file=dfile)
                    dfile.write('\n')
         
         
#- event exclusion based on energy levels.. ========================================================================                    
            # This should operate directly on the trace.
            if inp.exclude_events:
                ee.event_exclude(trace,inp.exclude_windows,inp.exclude_n,\
                inp.exclude_freq,inp.exclude_level)
                
    
            #- taper edges ========================================================================
    
            if inp.taper_do == True:
    
                trace=proc.taper(trace,inp.taper_width,verbose,ofid)
                
                if check == True:
                    dfile.write('Tapered\n')
                    print(trace.data[0:20],file=dfile)
                    dfile.write('\n')
            
            #- downsampling =======================================================================
            sampling_rate_index=0
            while sampling_rate_index<len(Fs_new):
                if trace.stats.sampling_rate>Fs_new[sampling_rate_index]:
                    trace=proc.downsample(trace,Fs_new[sampling_rate_index],\
                    verbose,ofid)
                sampling_rate_index+=1
            newtrace = trace.copy()
            del trace
               
            if check == True:
                dfile.write('(Downsampled), copied\n')
                print(newtrace.data[0:20],file=dfile)
                dfile.write('\n')   
            #- remove instrument response =========================================================
    
            if inp.remove_response == True:
    
                removed,newtrace=proc.remove_response(newtrace,respdir,unit,\
                freqs,wl,verbose,ofid)
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
                if inp.rankvariable != 'ALPS_APP_PE':
                    cstr.append(ctr)
                    cstr.plot(outfile=datadir+'/processed/out/'+\
                    filepath.split('/')[-1]+'.'+prepname+'.png',equal_scale=False)
                    cstr.trim(endtime=cstr[0].stats.starttime+3600)
                    cstr.plot(outfile=datadir+'/processed/out/'+\
                    filepath.split('/')[-1]+'.'+prepname+'.1hr.png',equal_scale=False)
                dfile.write('Instrument response removed\n')
                print(newtrace.data[0:20],file=dfile)
                dfile.write('\n')
                
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

        for trace_index_2 in range(len(colloc_data)):
            if ((inp.remove_response==True) and \
            (removed==1)) or \
                inp.remove_response==False:
                
                filepathnew = getfilepath(mydir,colloc_data[trace_index_2].stats,prepname)
                
                #- write to file
                colloc_data[trace_index_2].write(filepathnew,\
                format=colloc_data[trace_index_2].stats._format)
                       
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
        channel+'.' + t1 + '.*.'+prepname+'.'+format
    return filepathnew
        
    
