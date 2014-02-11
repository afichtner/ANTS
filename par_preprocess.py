# A script to process ambient vibration records

import os
import sys
import shutil
import time
from glob import glob

from obspy.core import read
from obspy import Stream,  Trace
from obspy.signal import filter

import numpy as np
from mpi4py import MPI
from math import ceil

import TOOLS.processing as proc
import TOOLS.normalisation as nrm
import TOOLS.read_xml as rxml 
import TOOLS.renamer as rn
import TOOLS.testinput as ti
import TOOLS.psd_estimate as psd


import par_preprocess as p

if __name__=='__main__':
    xmlin=str(sys.argv[1])
    try:
        rawdata=str(sys.argv[2])
    except IndexError:
        rawdata=None
    p.prep(xmlin, rawdata)


def prep(xmlinput,content=None):
    
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
        
        inp1=rxml.read_xml(xmlinput)[1]
        if ti.testinput(inp1)==False:
            print 'Problems in xmlinput, forced to interrupt.' 
            return
        startyr=int(inp1['input']['startyr'][0:4])
        endyr=int(inp1['input']['endyr'][0:4])
        prepname=inp1['prepname']
        #- copy the input xml to the output directory for documentation ===============================
        if os.path.exists('DATA/processed/xmlinput/proc.'+prepname+'.xml')==True:
            print '\n\nChoose a new name or delete inputfile of the name: proc.'+prepname+'.xml in ./xmlinput. Be aware this may cause chaos. Aborting.\n\n'
            return
        else:       
            shutil.copy(xmlinput,'DATA/processed/xmlinput/proc.'+prepname+'.xml')
        
        #- create yearly directories if not present===================================================
        for i in range(startyr-1,endyr+1):
            if os.path.exists('DATA/processed/'+str(i)+'/')==False:
                os.mkdir('DATA/processed/'+str(i))
                
        

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
            
  #==============================================================================================
  #- All processes:
  #- broadcast the input; and the list of files
  #- read variables from broadcasted input
  #==============================================================================================       
    else:
        content=list()
        inp1=list()    
       
    t1=time.time()-t0
    content=comm.bcast(content, root=0)
    inp1=comm.bcast(inp1, root=0)
    t2=time.time()-t0-t1
    
    verbose=bool(int(inp1['verbose']))
    check=bool(int(inp1['check']))
    prepname=inp1['prepname']
    saveplot=bool(int(inp1['saveplot']))
    startyr=int(inp1['input']['startyr'][0:4])
    endyr=int(inp1['input']['endyr'][0:4])
    ofid=open('DATA/processed/out/proc.'+prepname+'.rank_'+str(rank)+'.txt','w')
    
    #==============================================================================================
    #- Assign each rank its own chunk of input
    #==============================================================================================

    clen=int(ceil(float(len(content))/float(size)))
    chunk=(rank*clen, (rank+1)*clen)
    mycontent=content[chunk[0]:chunk[1]]
    if check==True:
        mycontent=[mycontent[0]]
    t3=time.time()-t0-t2   

    #==================================================================================
    # Input files loop
    #==================================================================================
    
    
    #- Print some nice comments to output file ----------------------------------------       
    if verbose:
        ofid.write('Time at start was '+str(t0)+'\n')
        ofid.write('Rank 0 took '+str(t1)+' seconds to read in input\n')
        ofid.write('Broadcasting took '+str(t2)+' seconds \n')
        ofid.write('I got my task assigned in '+str(t3)+' seconds \n')
        
        ofid.write('\nHi I am rank number %d and I am processing the following files for you: \n' %rank)
        for fname in mycontent:
            ofid.write(fname+'\n')
        
    for filepath in mycontent:
        
        filename=filepath.split('/')[-1]
  
        if verbose==True:
            ofid.write('\n========================================================================================\n')
            ofid.write('opening file: '+filepath+'\n')
            t=time.time()-t0
            ofid.write('Time elapsed since start: '+str(t)+'\n')
        #- read data
        try:
            data=read(filepath)
        except TypeError:
            if verbose==True: 
                ofid.write('file could not be opened, skip.')
            continue
        except IOError:
            if verbose: 
                ofid.write('file could not be opened, skip.')
            continue
        
        #- initialize a stream that is plotted in case check option is set
        if check==True:
            cstr=Stream()
        #- initialize the stream that recollects the trace segments
        colloc_data=Stream()
        
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- decimate-first routine +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if inp1['first_step']=='decimate':
                
            #- trim ===============================================================================
            
            if inp1['processing']['trim']=='1':
                data=proc.trim_next_sec(data,verbose, ofid)
                
                
            #- downsampling =======================================================================
            
            if inp1['processing']['decimation']['doit']=='1':
                new_fs=inp1['processing']['decimation']['new_sampling_rate'].split(' ')
                
                for fs in new_fs:
                    data=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose,ofid)
                    data=proc.downsample(data,float(fs),verbose,ofid)
            
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
        #- split traces into shorter segments======================================================
        if inp1['processing']['split']['doit']=='1':
            data=proc.split_traces(data,float(inp1['processing']['split']['length_in_sec']),float(inp1['quality']['min_length_in_sec']),verbose,ofid)
        
        if check==True:
            n_traces=min(3,len(data))
        else:
            n_traces=len(data)
        
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- split-first routine ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if inp1['first_step']=='split':
            #- trim ===============================================================================
            
            if inp1['processing']['trim']=='1':
                data=proc.trim_next_sec(data,verbose,ofid)
            
        #- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
        
        if verbose==True:
            ofid.write('contains '+str(n_traces)+' trace(s)\n')
    
        
        #==================================================================================
        # trace loop
        #==================================================================================
        for k in np.arange(n_traces):
    
            t_start=time.time()

            trace=data[k]
            
            if check:
                ctr=Trace(data=trace.data)
                ctr.stats.network='Original Data'
                cstr.append(ctr)
                
            if verbose==True: 
                ofid.write('- trace '+str(k+1)+' ----------------------------------------------------\n')
    
            #==================================================================================
            # basic quality checks
            #==================================================================================
    
            #- check NaN
            if True in np.isnan(trace.data):
                if verbose==True:
                    ofid.write('** trace contains NaN, discarded\n')
                continue
    
            #- check infinity
            if True in np.isinf(trace.data):
                if verbose==True:
                    ofid.write('** trace contains infinity, discarded\n')
                continue
    
            if verbose==True:
                ofid.write('* number of points: '+str(trace.stats.npts)+'\n')

            #==================================================================================
            # processing (detrending, filtering, response removal, decimation)
            #==================================================================================
                
            #- demean============================================================================
    
            if inp1['processing']['demean']=='1':
    
                trace=proc.demean(trace,verbose,ofid)
             
            #- detrend ============================================================================
          
            if inp1['processing']['detrend']=='1':
    
                trace=proc.detrend(trace,verbose,ofid)

            #- taper edges ========================================================================
    
            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose,ofid)
                
            #- bandpass, first stage ==============================================================
    
            if inp1['processing']['bandpass_1']['doit']=='1':
    
                trace=proc.bandpass(trace,int(inp1['processing']['bandpass_1']['corners']),float(inp1['processing']['bandpass_1']['f_min']),float(inp1['processing']['bandpass_1']['f_max']),verbose,ofid)
                
                if check:
                    ctr=Trace(data=trace.data)
                    ctr.stats.network='After Bandpass 1'
                    cstr.append(ctr)
                    
            
                
            #- remove instrument response =========================================================
    
            #ta=time.time()
            if inp1['processing']['instrument_response']['doit']=='1':
    
                removed,trace=proc.remove_response(trace,inp1['processing']['instrument_response']['respdir'],inp1['processing']['instrument_response']['unit'],inp1['processing']['instrument_response']['waterlevel'],verbose,ofid)
                if ((True in np.isnan(trace)) or (removed==0)):
                    ofid.write('Deconvolution seems unstable or instrument response was without succes. Trace discarded.')
                    continue
                if check:
                    ctr=Trace(data=trace.data)
                    ctr.stats.network='After IC to '+inp1['processing']['instrument_response']['unit']
                    cstr.append(ctr)
                
            #- split-first routine ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if inp1['first_step']=='split':
            #- downsampling =======================================================================
            
                if inp1['processing']['decimation']['doit']=='1':
                    new_fs=inp1['processing']['decimation']['new_sampling_rate'].split(' ')
                    for fs in new_fs:
                        data=proc.downsample(data,float(fs),verbose,ofid)
                
            #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            
            #ofid.write('+++++++++ '+str(time.time()-ta)+'\n')
            #- bandpass, second stage =============================================================
    
            if inp1['processing']['bandpass_2']['doit']=='1':
    
                trace=proc.bandpass(trace,int(inp1['processing']['bandpass_2']['corners']),float(inp1['processing']['bandpass_2']['f_min']),float(inp1['processing']['bandpass_2']['f_max']),verbose,ofid)
                
                if check:
                    ctr=Trace(data=trace.data)
                    ctr.stats.network='After Bandpass 2'
                    cstr.append(ctr)

            #- taper edges ========================================================================

            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose,ofid)
               
            #======================================================================================
            # normalisations (whitening, averages, ...)
            #======================================================================================
    
            #- one-bit normalisation ==============================================================
    
            if inp1['normalisation']['onebit']=='1':
            
                trace=nrm.onebit(trace,verbose)
    
            #- rms clipping =======================================================================
    
            if inp1['normalisation']['rms_clipping']=='1':
            
                trace=nrm.clip(trace,verbose)
    
            #- running average normalisation ======================================================
    
            if inp1['normalisation']['ram_normal']['doit']=='1':
    
                trace=nrm.ram_normal(trace,float(inp1['normalisation']['ram_normal']['window_length']),verbose)
    
            #- iterative clipping above a multiple of the rma (waterlevel normalisation) ==========
    
            if inp1['normalisation']['waterlevel']['doit']=='1':
    
                trace=nrm.waterlevel(trace,float(inp1['normalisation']['waterlevel']['level']),verbose)
    
            #- spectral whitening =================================================================
    
            if inp1['normalisation']['whitening']['doit']=='1':
    
                trace=nrm.whiten(trace,float(inp1['normalisation']['whitening']['smoothing']),verbose)

            #======================================================================================
            # timing and storage of results
            #======================================================================================
            if check:
                ctr=Trace(data=trace.data)
                ctr.stats.network='After preprocessing'
                cstr.append(ctr)
                cstr.plot(outfile='DATA/processed/out/'+filepath.split('/')[-1]+'.'+str(k)+'.'+prepname+'.png',equal_scale=False)
                cstr.trim(endtime=cstr[0].stats.starttime+3600)
                cstr.plot(outfile='DATA/processed/out/'+filepath.split('/')[-1]+'.'+str(k)+'.'+prepname+'.1hr.png',equal_scale=False)
                
                
            if verbose==True:
                t_end=time.time()
                ofid.write('time per trace: '+str(t_end-t_start)+' s\n')

            #- merge all into final trace
            colloc_data+=trace
        
        colloc_data._cleanup()
        
        
        if verbose==True:
            t=time.time()-t0
            ofid.write('- rename and store ------------------------------------------------------\n')
            ofid.write('time elapsed since start: '+str(t)+'\n')
            
        if len(colloc_data)>0:
            if (inp1['saveprep']=='1'):
                for k in range(len(colloc_data)):
                    rn.rename_seismic_data(colloc_data[k], prepname, verbose, ofid)
                    
    ofid.close()
    
