# A script to process ambient vibration records

import os
import sys
import shutil

from obspy.core import read
from obspy import Stream,  Trace
import matplotlib.pyplot as plt
import numpy as np

import TOOLS.processing as proc
import TOOLS.normalisation as nrm
import TOOLS.read_xml as rxml 
import TOOLS.renamer as rn

def preprocessing(processing_input):

    #==============================================================================================
    # preliminaries
    #==============================================================================================

    #- read input files ===========================================================================

    inp1=rxml.read_xml(processing_input)
    inp1=inp1[1]

    indirs=inp1['directories']['indirs'].strip().split(' ')
    outdir=inp1['directories']['outdir']

    verbose=bool(int(inp1['verbose']))
    plot=bool(int(inp1['plot']))

    #- make target directory if it does not exist =================================================
    
    if os.path.exists(outdir)==False:
        os.mkdir(outdir)

    #- copy the input xml to the output directory for documentation ===============================
    shutil.copy(processing_input,outdir)
    
    #- loop through all files =====================================================================

    #- loop through input directories -------------------------------------------------------------
    for indir in indirs:

        content=os.listdir(indir)
        
        #==================================================================================
        # Input files loop
        #==================================================================================
        for filename in content:
          
            filepath=indir+'/'+filename
        
            if verbose==True:
                print '\nopening file: '+filepath
                
            #- read data
            try:
                data=read(filepath)
            except TypeError:
                if verbose==True: print 'file could not be opened, skip.'
                continue
            except IOError:
                if verbose: print 'file could not be opened, skip.'
                continue
            
            #- Lowpass filter==========================
            data=proc.taper(data,float(inp1['processing']['taper']['taper_width']),verbose)
            data=proc.lowpass(data,4,float(inp1['processing']['decimation']['new_sampling_rate'])*0.25,verbose )
            
            #- split traces into shorter segments==============================================
            if inp1['processing']['split']['doit']=='1':
                data=proc.split_traces(data,float(inp1['processing']['split']['length_in_sec']),float(inp1['quality']['min_length_in_sec']),verbose)
            n_traces=len(data)
            
            
            if verbose==True:
                print 'contains '+str(n_traces)+' trace(s)'

          
            
            #==================================================================================
            # trace loop
            #==================================================================================
            for k in np.arange(n_traces):

                trace=data[k]
                
                if verbose==True: print '-----------------------------------------------------------'

                if plot==True:
                    print '* trace before processing'
                    trace.plot()

                #==================================================================================
                # basic quality checks
                #==================================================================================

                #- check NaN
                if True in np.isnan(trace.data):
                    if verbose==True: print '** trace contains NaN, discarded'
                    continue

                #- check infinity
                if True in np.isinf(trace.data):
                    if verbose==True: print '** trace contains infinity, discarded'
                    continue

                if verbose==True: print '* number of points: '+str(trace.stats.npts)

                #==================================================================================
                # processing (detrending, filtering, response removal, decimation)
                #==================================================================================
                #- taper edges ========================================================================

                if inp1['processing']['taper']['doit']=='1':
                    trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose)
                    
                #- trim ===========================================================================
                
                if inp1['processing']['trim']=='1':
                    trace=proc.trim_next_sec(trace,verbose)
                    
                    
                #- detrend ============================================================================

                if inp1['processing']['detrend']=='1':
                    trace=proc.detrend(trace,verbose)
                    
                    
                #- demean============================================================================

                if inp1['processing']['demean']=='1':
                    trace=proc.demean(trace,verbose)
                    
                   
                #- taper edges ========================================================================

                if inp1['processing']['taper']['doit']=='1':
                   trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose)
                   
            
                #- bandpass, first stage ==============================================================

                if inp1['processing']['bandpass_1']['doit']=='1':
                    trace=proc.bandpass(trace,int(inp1['processing']['bandpass_1']['corners']),float(inp1['processing']['bandpass_1']['f_min']),float(inp1['processing']['bandpass_1']['f_max']),verbose)
              
                    
                    
                #- downsampling =======================================================================
    
                if inp1['processing']['bandpass_1']['doit']=='1' and inp1['processing']['decimation']['doit']=='1':
                    trace=proc.downsample(trace,inp1['processing']['decimation']['new_sampling_rate'],verbose)
                elif inp1['processing']['bandpass_1']['doit']!='1' and inp1['processing']['decimation']['doit']=='1':
                    print "It is necessary to apply an anti-alias filter first (see bandpass first stage in input file).\n No filter has been defined, therefore a fourth-order lowpass filter with fc=0.25*new sampling rate is now imposed."
                    trace=proc.lowpass(trace, 4, 0.25*float(inp1['processing']['decimation']['new_sampling_rate']), verbose)
                    trace=proc.downsample(trace,inp1['processing']['decimation']['new_sampling_rate'],verbose)
       
                   
                #- remove instrument response =========================================================

                if inp1['processing']['instrument_response']['doit']=='1':

                    removed,trace=proc.remove_response(trace,inp1['processing']['instrument_response']['respdir'],inp1['processing']['instrument_response']['unit'],inp1['processing']['instrument_response']['waterlevel'],verbose)

                #- bandpass, second stage =============================================================

                if inp1['processing']['bandpass_2']['doit']=='1':

                    trace=proc.bandpass(trace,int(inp1['processing']['bandpass_2']['corners']),float(inp1['processing']['bandpass_2']['f_min']),float(inp1['processing']['bandpass_2']['f_max']),verbose)
                    trace.plot()

                #- taper edges ========================================================================

                if inp1['processing']['taper']['doit']=='1':

                    trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose)

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
                # plot and store results
                #======================================================================================
            
                
            
                if len(trace)>0:
                    if plot:
                        print '* trace after processing'
                        trace.plot()

                        
                    
                    if (inp1['saveprep']=='1'):
                        if ((inp1['processing']['instrument_response']['doit']=='1') and (removed==1)) or (inp1['processing']['instrument_response']['doit']!='1'):
                            rn.rename_seismic_data(trace, outdir, True, verbose)
                    
