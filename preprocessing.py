# A script to process ambient vibration records

import os
import sys

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

    #- make target directory if it does not exist =================================================
    
    if os.path.exists(outdir)==False:
        os.mkdir(outdir)
    
    #- loop through all files =====================================================================

    # number of successfully processed seismograms
    count_success=0
    # number of seismograms where processing failed
    count_fail=0

    #- loop through input directories
    for indir in indirs:

        content=os.listdir(indir)

        #- loop through input files
        for filename in content:
        
            filepath=indir+'/'+filename
        
            if verbose==True:
                print '\nopening file: '+filepath
                
            #obtain trace containing this file's data
            data=read(filepath)[0]

            data_original=data.copy()

            #======================================================================================
            # processing (detrending, filtering, response removal, decimation)
            #======================================================================================

            #- detrend ============================================================================

            if inp1['processing']['detrend']=='1':

                data=proc.detrend(data,verbose)

            #- taper edges ========================================================================

            if inp1['processing']['taper']['doit']=='1':

                data=proc.taper(data,float(inp1['processing']['taper']['taper_width']),verbose)

            #- bandpass, first stage ==============================================================

            if inp1['processing']['bandpass_1']['doit']=='1':

                data=proc.bandpass(data,int(inp1['processing']['bandpass_1']['corners']),float(inp1['processing']['bandpass_1']['f_min']),float(inp1['processing']['bandpass_1']['f_max']),verbose)

            #- downsampling =======================================================================

            if inp1['processing']['decimation']['doit']=='1':

               data=proc.downsample(data,inp1['processing']['decimation']['new_sampling_rate'],verbose)

            #- remove instrument response =========================================================

            if inp1['processing']['instrument_response']['doit']=='1':

                success,data=proc.remove_response(data,inp1['processing']['instrument_response']['respdir'],inp1['processing']['instrument_response']['unit'],verbose)
                count_success=count_success+success
                count_fail=count_fail+1-success

            #- bandpass, second stage =============================================================

            if inp1['processing']['bandpass_2']['doit']=='1':

                data=proc.bandpass(data,int(inp1['processing']['bandpass_2']['corners']),float(inp1['processing']['bandpass_2']['f_min']),float(inp1['processing']['bandpass_2']['f_max']),verbose)


            #======================================================================================
            # normalisations (whitening, averages, ...)
            #======================================================================================

            #- one-bit normalisation ==============================================================

            if inp1['normalisation']['onebit']=='1':
                
                data=nrm.onebit(data,verbose)

            #- rms clipping =======================================================================

            if inp1['normalisation']['rms_clipping']=='1':
                
                data=nrm.clip(data,verbose)

            #- running average normalisation ======================================================

            if inp1['normalisation']['ram_normal']['doit']=='1':

                data=nrm.ram_normal(data,float(inp1['normalisation']['ram_normal']['window_length']),verbose)

            #- iterative clipping above a multiple of the rma (waterlevel normalisation) ==========

            if inp1['normalisation']['waterlevel']['doit']=='1':

                data=nrm.waterlevel(data,float(inp1['normalisation']['waterlevel']['level']),verbose)

            #- spectral whitening =================================================================

            if inp1['normalisation']['whitening']['doit']=='1':

                data=nrm.whiten(data,float(inp1['normalisation']['whitening']['smoothing']),verbose)

            #======================================================================================
            # store results
            #======================================================================================

            #- rename files =======================================================================

            if inp1['rename']=='1':
                if inp1['saveprep']=='1':
                    #- write processed seismograms
                    rn.rename_seismic_data(data, outdir, True, verbose)
                    
                #- write original seismograms
                rn.rename_seismic_data(data_original, outdir, False, verbose)
            else:
                #- write processed seismograms
                data.write(filepath+'_prep',data.stats._format)


            #- test spectral whitening
        
            #(freq, spec)=fourier_spectrum(data[0].data, Fs)
            #plt.plot(freq[1000:5000], abs(spec[1000:5000]))
            #plt.xlabel('Amplitude [?]')    
            #plt.xlabel('Frequency [Hz]')
            #plt.show()
            #psdplot(data[0].data, Fs)
