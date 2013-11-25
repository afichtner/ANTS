# A script to process ambient vibration records

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
import TOOLS.psd_estimate as psd




if __name__=='__main__':
    
    import preprocessing as p
    xmlin=str(sys.argv[1])
    print 'XML input file: '+ xmlin
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
        print 'Problems in xmlinput, forced to interrupt.'
        return
    
    outdir=inp1['directories']['outdir']

    verbose=bool(int(inp1['verbose']))
    check=bool(int(inp1['check']))
    saveplot=bool(int(inp1['saveplot']))
    

    #- make target directory if it does not exist =================================================
    
    if os.path.exists(outdir)==False:
        os.mkdir(outdir)

    #- copy the input xml to the output directory for documentation ===============================
    
    shutil.copy(xmlinput,outdir)
    
    #- check what input is, list input from different directories =================================
    if content==None:
        indirs=inp1['directories']['indirs'].strip().split(' ')
        content=list()
        for indir in indirs:
            content.extend(os.listdir(indir))
            for k in range(len(content)):
                content[k]=indir+'/'+content[k]

    elif type(content)==str: 
        filename=content
        content=list()
        content.append(filename)
   
    #==================================================================================
    # Input files loop
    #==================================================================================
    
    for filepath in content:
        print filepath
        filename=filepath.split('/')[-1]
        
    
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
        
        
        #- decimate-first routine +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if inp1['first_step']=='decimate':
            #- taper edges ========================================================================
    
#            if inp1['processing']['taper']['doit']=='1':
#                data=proc.taper(data,float(inp1['processing']['taper']['taper_width']),verbose)
                
            #- trim ===============================================================================
            
            if inp1['processing']['trim']=='1':
                data=proc.trim_next_sec(data,verbose)
                
            #- downsampling =======================================================================
            
            if inp1['processing']['decimation']['doit']=='1':
                new_fs=inp1['processing']['decimation']['new_sampling_rate'].split(' ')
                #data.filter('lowpassCheby2', freq=float(new_fs[-1])*0.25, maxorder=8)
                data=proc.lowpass(data,4,float(new_fs[-1])*0.2,verbose )
                
                for fs in new_fs:
                    data=proc.downsample(data,float(fs),verbose)
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
        
    
        #- split traces into shorter segments======================================================
        if inp1['processing']['split']['doit']=='1':
            data=proc.split_traces(data,float(inp1['processing']['split']['length_in_sec']),float(inp1['quality']['min_length_in_sec']),verbose)
        n_traces=len(data)
        
        
        #- split-first routine ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if inp1['first_step']=='split':
            #- trim ===============================================================================
            
            if inp1['processing']['trim']=='1':
                data=proc.trim_next_sec(data,verbose)
        #- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            
            
        
        if verbose==True:
            print 'contains '+str(n_traces)+' trace(s)'
    
        #data_original=data.copy()
        colloc_data=Stream()
        
        #==================================================================================
        # trace loop
        #==================================================================================
        for k in np.arange(n_traces):
    
            trace=data[k]
            if check:
                trace_original=trace.copy()
            
            if verbose==True: print '-----------------------------------------------------------'
    
            if check==True:
                print '* trace before processing'
                trace.plot()
                print '* spectrum before processing'
                psd.plot_psd(trace,100, trace_original)
    
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
                
            #- detrend ============================================================================
    
#            if inp1['processing']['detrend']=='1':
#    
#                trace=proc.detrend(trace,verbose)
#                
            #- demean============================================================================
    
            if inp1['processing']['demean']=='1':
    
                trace=proc.demean(trace,verbose)
            if inp1['processing']['detrend']=='1':
    
                trace=proc.detrend(trace,verbose)
                    
    
            #- taper edges ========================================================================
    
            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose)
                
          
            #- bandpass, first stage ==============================================================
    
            if inp1['processing']['bandpass_1']['doit']=='1':
    
                trace=proc.bandpass(trace,int(inp1['processing']['bandpass_1']['corners']),float(inp1['processing']['bandpass_1']['f_min']),float(inp1['processing']['bandpass_1']['f_max']),verbose)
                if check==True:
                    print '* spectrum after bandpass nr 1'
                    psd.plot_psd(trace,100, trace_original)
                
                
            #- split-first routine ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if inp1['first_step']=='split':
            #- downsampling =======================================================================
            
                if inp1['processing']['decimation']['doit']=='1':
                    trace=proc.downsample(trace,inp1['processing']['decimation']['new_sampling_rate'],verbose)
                
            #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
                
            #- remove instrument response =========================================================
    
            if inp1['processing']['instrument_response']['doit']=='1':
    
                removed,trace=proc.remove_response(trace,inp1['processing']['instrument_response']['respdir'],inp1['processing']['instrument_response']['unit'],inp1['processing']['instrument_response']['waterlevel'],verbose)
                if True in np.isnan(trace):
                    print 'Deconvolution seems unstable! Trace discarded.'
                    continue
                if check==True:
                    print '* spectrum after Instrument correction'
                    psd.plot_psd(trace,100, trace_original)
    
            #- bandpass, second stage =============================================================
    
            if inp1['processing']['bandpass_2']['doit']=='1':
    
                trace=proc.bandpass(trace,int(inp1['processing']['bandpass_2']['corners']),float(inp1['processing']['bandpass_2']['f_min']),float(inp1['processing']['bandpass_2']['f_max']),verbose)
    
            #- taper edges ========================================================================
    
            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose)
               
                if check==True:
                    print '* spectrum after bandpass nr 2 and taper'
                    psd.plot_psd(trace,100, trace_original)
                    
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
            if check==True:
                print '* spectrum after processing'
                psd.plot_psd(trace,100, trace_original)
                print '* single trace after processing:'
                trace.plot()
                
            #merge all into final trace
            colloc_data+=trace
        
        colloc_data._cleanup()
        print '* traces after processing are:\n', colloc_data
        
        if len(colloc_data)>0:
            if check:
                colloc_data.plot()
                if saveplot:
                    figname=outdir+'/'+colloc_data[0].stats.station+'.'+colloc_data[0].stats.channel+'.png'
                    colloc_data.plot(outfile=figname)
                    print '* plot saved to ', figname
                psd.plot_psd(colloc_data[0],100)
                    
                    
            if (inp1['saveprep']=='1'):
                for k in range(len(colloc_data)):
                    if ((inp1['processing']['instrument_response']['doit']=='1') and (removed==1)) or (inp1['processing']['instrument_response']['doit']!='1'):
                        rn.rename_seismic_data(colloc_data[k], outdir, True, verbose)
                    
    
    
