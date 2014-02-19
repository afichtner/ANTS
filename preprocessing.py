# A script to process ambient vibration records
from __future__ import print_function
# Use the print function to be able to switch easily between stdout and a file
import os
import sys
import shutil


from obspy.core import read
from obspy import Stream,  Trace
from obspy.signal import filter
import matplotlib.pyplot as plt
import numpy as np

import TOOLS.processing as proc
import TOOLS.normalisation as nrm
import TOOLS.read_xml as rxml 
import TOOLS.renamer as rn
import TOOLS.testinput as ti
from  TOOLS.plot_spectra import plotpsdwpete
import TOOLS.psd_estimate as psd



if __name__=='__main__':
    import preprocessing as p
    xmlin=str(sys.argv[1])
    print('XML input file: '+ xmlin,file=None)
    p.prep(xmlin)


def prep(xmlinput,content=None):
    
    """
    
    This script preprocesses the MSEED files at the path specified as command line argument 2.
    Command line argument 1 must be xml input file.
    
    """

    #==============================================================================================
    # preliminaries
    #==============================================================================================
   
    #- read input files ===========================================================================
    
    
    inp1=rxml.read_xml(xmlinput)
    inp1=inp1[1]
        
    if ti.testinput(inp1)==False:
        print('Problems in xmlinput, forced to interrupt.', file=None)
        return
    
    
    verbose=bool(int(inp1['verbose']))
    check=bool(int(inp1['check']))
    outf=bool(int(inp1['outfile']))
    prepname=inp1['prepname']
    saveplot=bool(int(inp1['saveplot']))
    startyr=int(inp1['input']['startyr'][0:4])
    endyr=int(inp1['input']['endyr'][0:4])
    
    
    
    #- copy the input xml to the output directory for documentation ===============================
    if os.path.exists('DATA/processed/xmlinput/proc.'+prepname+'.xml')==True:
        print('\n\nChoose a new name or delete inputfile of the name: proc.'+prepname+'.xml in ./xmlinput. Be aware this may cause chaos. Aborting.\n\n',file=None)
        return
        
    shutil.copy(xmlinput,'DATA/processed/xmlinput/proc.'+prepname+'.xml')
    
    if check==True or outf==True:
        outfile=open('DATA/processed/out/proc.'+prepname+'.txt','w')
    else:
        outfile=None
    
    for i in range(startyr-1,endyr+1):
        if os.path.exists('DATA/processed/'+str(i)+'/')==False:
            os.mkdir('DATA/processed/'+str(i))
    
    #- check what input is, list input from different directories =================================
    if content==None:
        indirs=inp1['input']['indirs'].strip().split(' ')
        content=list()
        for indir in indirs:
            content.extend(os.listdir(indir))
            for k in range(len(content)):
                content[k]=indir+'/'+content[k]

    elif type(content)==str: 
        filename=content
        content=list()
        content.append(filename)
        
    #- If only a check run is performed, then only a couple of files are preprocessed
    if check:
        content=[content[0],content[1],content[len(content)-2],content[len(content)-1]]
        
   
    #==================================================================================
    # Input files loop
    #==================================================================================
    
    for filepath in content:
        
        
        if verbose==True:
            print('===========================================================',file=outfile)
            print('* opening file: '+filepath+'\n',file=outfile)
            
        
        #- read data
        try:
            data=read(filepath)
        except TypeError:
            if verbose==True: print('** file could not be opened, skip.',file=outfile)
            continue
        except IOError:
            if verbose: print('** file could not be opened, skip.',file=outfile)
            continue
        
        #- initialize stream to 'recollect' the split traces

        colloc_data=Stream()
        
        #- decimate-first routine +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if inp1['first_step']=='decimate':
            #- trim ===============================================================================
            
            if inp1['processing']['trim']=='1':
                data=proc.trim_next_sec(data,verbose,outfile)
                
            #- downsampling =======================================================================
            
            if inp1['processing']['decimation']['doit']=='1':
                new_fs=inp1['processing']['decimation']['new_sampling_rate'].split(' ')
                for fs in new_fs:
                    data=proc.taper(data,float(inp1['processing']['taper']['taper_width']),verbose,outfile)
                    data=proc.downsample(data,float(fs),verbose,outfile)
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
        
    
        #- split traces into shorter segments======================================================
        if inp1['processing']['split']['doit']=='1':
            data=proc.split_traces(data,float(inp1['processing']['split']['length_in_sec']),float(inp1['quality']['min_length_in_sec']),verbose,outfile)
        
        if check:
            n_traces=min(3,len(data))
        else:
            n_traces=len(data)
        
        
        #- split-first routine ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if inp1['first_step']=='split':
            #- trim ===============================================================================
            
            if inp1['processing']['trim']=='1':
                data=proc.trim_next_sec(data,verbose,outfile)
        #- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            
            
        
        if verbose==True:
            print('* contains '+str(n_traces)+' trace(s)',file=outfile)
    
        
        
        #==================================================================================
        # trace loop
        #==================================================================================
        for k in np.arange(n_traces):
    
            trace=data[k]
            
            if check:
                ctr=trace.copy()
                ctr.stats.network='Original Data'
                ctr.stats.station=''
                ctr.stats.location=''
                ctr.stats.channel=''
                cstr=Stream(ctr)
                
            
            if verbose==True: print('-----------------------------------------------------------',file=outfile)
    
            #==================================================================================
            # basic quality checks
            #==================================================================================
    
            #- check NaN
            if True in np.isnan(trace.data):
                if verbose==True: print('** trace contains NaN, discarded',file=outfile)
                continue
    
            #- check infinity
            if True in np.isinf(trace.data):
                if verbose==True: print('** trace contains infinity, discarded',file=outfile)
                continue
    
            if verbose==True: print('* number of points: '+str(trace.stats.npts)+'\n',file=outfile)
    
            #==================================================================================
            # processing (detrending, filtering, response removal, decimation)
            #==================================================================================
                              
            #- demean============================================================================
    
            if inp1['processing']['demean']=='1':
    
                trace=proc.demean(trace,verbose,outfile)
            if inp1['processing']['detrend']=='1':
    
                trace=proc.detrend(trace,verbose,outfile)
                    
    
            #- taper edges ========================================================================
    
            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose,outfile)
                
          
            #- bandpass, first stage ==============================================================
    
            if inp1['processing']['bandpass_1']['doit']=='1':
    
                trace=proc.bandpass(trace,int(inp1['processing']['bandpass_1']['corners']),float(inp1['processing']['bandpass_1']['f_min']),float(inp1['processing']['bandpass_1']['f_max']),verbose,outfile)
                
                if check:
                    ctr=trace.copy()
                    ctr.stats.network='After Bandpass 1'
                    ctr.stats.station=''
                    ctr.stats.location=''
                    ctr.stats.channel=''
                    cstr.append(ctr)
                
            #- split-first routine ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if inp1['first_step']=='split':
            #- downsampling =======================================================================
            
                if inp1['processing']['decimation']['doit']=='1':
                    trace=proc.downsample(trace,inp1['processing']['decimation']['new_sampling_rate'],verbose,outfile)
                
            #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
                
            #- remove instrument response =========================================================
    
            if inp1['processing']['instrument_response']['doit']=='1':
    
                removed,trace=proc.remove_response(trace,inp1['processing']['instrument_response']['respdir'],inp1['processing']['instrument_response']['unit'],inp1['processing']['instrument_response']['waterlevel'],verbose,outfile)
                if True in np.isnan(trace):
                    print('** Deconvolution seems unstable! Trace discarded.',file=outfile)
                    continue
                if check:
                    ctr=trace.copy()
                    ctr.stats.network='After IC to '+inp1['processing']['instrument_response']['unit']
                    ctr.stats.station=''
                    ctr.stats.location=''
                    ctr.stats.channel=''
                    cstr.append(ctr)
              
            #- bandpass, second stage =============================================================
    
            if inp1['processing']['bandpass_2']['doit']=='1':
    
                trace=proc.bandpass(trace,int(inp1['processing']['bandpass_2']['corners']),float(inp1['processing']['bandpass_2']['f_min']),float(inp1['processing']['bandpass_2']['f_max']),verbose,outfile)
                if check:
                    ctr=trace.copy()
                    ctr.stats.network='After Bandpass 2'
                    ctr.stats.station=''
                    ctr.stats.location=''
                    ctr.stats.channel=''
                    cstr.append(ctr)
                    
            #- taper edges ========================================================================
    
            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose,outfile)
               
            #======================================================================================
            # normalisations (whitening, averages, ...)
            #======================================================================================
    
            #- one-bit normalisation ==============================================================
    
            if inp1['normalisation']['onebit']=='1':
            
                trace=nrm.onebit(trace,verbose,outfile)
    
            #- rms clipping =======================================================================
    
            if inp1['normalisation']['rms_clipping']=='1':
            
                trace=nrm.clip(trace,verbose,outfile)
    
            #- running average normalisation ======================================================
    
            if inp1['normalisation']['ram_normal']['doit']=='1':
    
                trace=nrm.ram_normal(trace,float(inp1['normalisation']['ram_normal']['window_length']),verbose,outfile)
    
            #- iterative clipping above a multiple of the rma (waterlevel normalisation) ==========
    
            if inp1['normalisation']['waterlevel']['doit']=='1':
    
                trace=nrm.waterlevel(trace,float(inp1['normalisation']['waterlevel']['level']),verbose,outfile)
    
            #- spectral whitening =================================================================
    
            if inp1['normalisation']['whitening']['doit']=='1':
    
                trace=nrm.whiten(trace,float(inp1['normalisation']['whitening']['smoothing']),verbose,outfile)
            
            #======================================================================================
            # plot and store results
            #======================================================================================
            if check:
                ctr=trace.copy()
                ctr.stats.network='After preprocessing'
                ctr.stats.station=''
                ctr.stats.location=''
                ctr.stats.channel=''
                cstr.append(ctr)
                cstr.plot(outfile='DATA/processed/out/'+filepath.split('/')[-1]+'.'+prepname+'.png',equal_scale=False)
                cstr.trim(endtime=cstr[0].stats.starttime+3600)
                cstr.plot(outfile='DATA/processed/out/'+filepath.split('/')[-1]+'.'+prepname+'.1hr.png',equal_scale=False)
            #merge all into final trace
            colloc_data+=trace
        
        colloc_data._cleanup()
        print('* Done Preprocessing. Output are:\n',len(colloc_data),' traces.',file=outfile)
        
                    
                   
        if (inp1['saveprep']=='1'):
            for k in range(len(colloc_data)):
                if ((inp1['processing']['instrument_response']['doit']=='1') and (removed==1)) or (inp1['processing']['instrument_response']['doit']!='1',outfile):
                    rn.rename_seismic_data(colloc_data[k],prepname,verbose,outfile)
    if outfile:
        outfile.close()                  