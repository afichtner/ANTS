# A script to process ambient vibration records

import os
import sys
import shutil

from obspy.core import read

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

        #- loop through input files ---------------------------------------------------------------
        for filename in content:
            if filename=='.DS_Store': continue
            
            filepath=indir+'/'+filename
        
            if verbose==True:
                print '\nopening file: '+filepath
                
            #- read data
            try:
                data=read(filepath)
            except TypeError:
                if verbose==True: print 'file could not be opened, abort'
                continue

            #- split traces into shorter segments
            if inp1['processing']['split']['doit']=='1':
                data=proc.split_traces(data,float(inp1['processing']['split']['length_in_sec']),float(inp1['quality']['min_length_in_sec']),verbose)

            n_traces=len(data)

            if verbose==True:
                print 'contains '+str(n_traces)+' trace(s)'

            data_original=data.copy()

            #- loop through traces ----------------------------------------------------------------
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
                #- trim ===========================================================================
            
                if inp1['processing']['trim']=='1':
                    trace=proc.trim_next_sec(trace,verbose)

                #- detrend ============================================================================

                if inp1['processing']['detrend']=='1':

                    trace=proc.detrend(trace,verbose)

                #- taper edges ========================================================================

                if inp1['processing']['taper']['doit']=='1':

                    trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose)
                    
                #- demean============================================================================

                if inp1['processing']['demean']=='1':

                    trace=proc.demean(trace,verbose)

                #- bandpass, first stage ==============================================================

                if inp1['processing']['bandpass_1']['doit']=='1':

                    trace=proc.bandpass(trace,int(inp1['processing']['bandpass_1']['corners']),float(inp1['processing']['bandpass_1']['f_min']),float(inp1['processing']['bandpass_1']['f_max']),verbose)

                #- downsampling =======================================================================

                if inp1['processing']['decimation']['doit']=='1':

                    trace=proc.downsample(trace,inp1['processing']['decimation']['new_sampling_rate'],verbose)

                #- remove instrument response =========================================================

                if inp1['processing']['instrument_response']['doit']=='1':

                    removed,trace=proc.remove_response(trace,inp1['processing']['instrument_response']['respdir'],inp1['processing']['instrument_response']['unit'],verbose)

                #- bandpass, second stage =============================================================

                if inp1['processing']['bandpass_2']['doit']=='1':

                    trace=proc.bandpass(trace,int(inp1['processing']['bandpass_2']['corners']),float(inp1['processing']['bandpass_2']['f_min']),float(inp1['processing']['bandpass_2']['f_max']),verbose)

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
            
                if plot==True:
                    print '* trace after processing'
                    trace.plot()
                    
                if (inp1['saveprep']=='1') & (removed==1):
                    rn.rename_seismic_data(trace, outdir, True, verbose)
                
