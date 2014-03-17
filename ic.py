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

import antconfig as cfg


if __name__=='__main__':
    import ic
    xmlin=str(sys.argv[1])
    print('XML input file: '+ xmlin,file=None)
    ic.ic(xmlin)


def ic(xmlinput,content=None):
    
    """
    
    This script preprocesses the MSEED files at the path specified as command line argument 2.
    Command line argument 1 must be xml input file.
    
    """

    #==============================================================================================
    # preliminaries
    #==============================================================================================
   
    #- read input files ===========================================================================
    
    datadir=cfg.datadir
    inp1=rxml.read_xml(xmlinput)[1]
    if testinput(inp1,ofid=None)==False:
        print('Problems in xmlinput, forced to interrupt.', file=None)
        return
    
    
    verbose=bool(int(inp1['verbose']))
    check=bool(int(inp1['check']))
    outf=bool(int(inp1['outfile']))
    prepname=inp1['prepname']
    saveplot=bool(int(inp1['saveplot']))
    startyr=int(inp1['input']['startyr'][0:4])
    endyr=int(inp1['input']['endyr'][0:4])
    respdir=inp1['processing']['instrument_response']['respdir']
    unit=inp1['processing']['instrument_response']['unit']
    freqs=inp1['processing']['instrument_response']['freqs']
    wl=inp1['processing']['instrument_response']['waterlevel']
    seglen=float(inp1['processing']['split']['length_in_sec'])
    minlen=float(inp1['quality']['min_length_in_sec'])
    
    #- copy the input xml to the output directory for documentation ===============================
    if os.path.exists(datadir+'/processed/xmlinput/ic.'+prepname+'.xml')==True:
        print('\n\nChoose a new name or delete inputfile of the name: ic.'+prepname+'.xml in ./xmlinput. Be aware this may cause chaos. Aborting.\n\n',file=None)
        return
        
    shutil.copy(xmlinput,datadir+'/processed/xmlinput/ic.'+prepname+'.xml')
    
    if check==True or outf==True:
        outfile=open(datadir+'/processed/out/ic.'+prepname+'.txt','w')
    else:
        outfile=None
    
    for i in range(startyr-1,endyr+1):
        if os.path.exists(datadir+'/processed/'+str(i)+'/')==False:
            os.mkdir(datadir+'/processed/'+str(i))
    
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
    if check and len(content)>4:
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
            if verbose==True: print('** file could not be opened, TypeError, skip.',file=outfile)
            continue
        except IOError:
            if verbose: print('** file could not be opened, IOError, skip.',file=outfile)
            continue
        except KeyError:
            if verbose: print('** file could not be opened, KeyError, skip.',file=outfile)
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
                    data=proc.downsample(data,float(fs),verbose,outfile)
        #- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
        
    
        #- split traces into shorter segments======================================================
        if inp1['processing']['split']['doit']=='1':
            data=proc.slice_traces(data,seglen,minlen,verbose,outfile)
        
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
            if inp1['processing']['detrend']=='1':
    
                trace=proc.detrend(trace,verbose,outfile)
            
            if inp1['processing']['demean']=='1':
    
                trace=proc.demean(trace,verbose,outfile)
            
    
            #- taper edges ========================================================================
    
            if inp1['processing']['taper']['doit']=='1':
    
                trace=proc.taper(trace,float(inp1['processing']['taper']['taper_width']),verbose,outfile)
                if check:
                    ctr=trace.copy()
                    ctr.stats.network='After Detrending'
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
    
                removed,trace=proc.remove_response(trace,respdir,unit,freqs,wl,verbose,outfile)
                if True in np.isnan(trace):
                    print('** Deconvolution seems unstable! Trace discarded.',file=outfile)
                    continue
                if check:
                    ctr=trace.copy()
                    ctr.stats.network='After IC to '+unit
                    ctr.stats.station=''
                    ctr.stats.location=''
                    ctr.stats.channel=''
                    cstr.append(ctr)
              
         
            if check:
                cstr.plot(outfile=datadir+'/processed/out/'+filepath.split('/')[-1]+'.'+prepname+'.png',equal_scale=False)
                cstr.trim(endtime=cstr[0].stats.starttime+3600)
                cstr.plot(outfile=datadir+'/processed/out/'+filepath.split('/')[-1]+'.'+prepname+'.1hr.png',equal_scale=False)
            #merge all into final trace
            colloc_data+=trace
        
        colloc_data._cleanup()
                  
        if (inp1['saveprep']=='1'):
            for k in range(len(colloc_data)):
                if ((inp1['processing']['instrument_response']['doit']=='1') and (removed==1)) or (inp1['processing']['instrument_response']['doit']!='1',outfile):
                    rn.rename_seismic_data(colloc_data[k],prepname,verbose,outfile)
    if outfile:
        outfile.close()
        
        
def testinput(inp,ofid):
    niceinput=True
    
    if not inp['prepname']:
        print('A valid name must be provided for this processing run.',file=ofid)
        niceinput=False
        
    if not inp['comment']:
        print('It is necessary that you provide a comment on this preprocessing run.',file=ofid)
        niceinput=False
        
    
    if inp['verbose']!='1' and inp['verbose']!='0':
        print('Verbose must be 1 or 0.',file=ofid)
        niceinput=False
        
    if inp['check']!='1' and inp['check']!='0':
        print('check must be 1 or 0.',file=ofid)
        niceinput=False
        
    try:
        sm=int(inp['input']['startyr'])
        if sm<1970:
            print('Start year must be format YYYY, no earlier than 1970.',file=ofid)
            niceinput=False
         
    except ValueError:
        print('Start year must be an integer, format YYYY.',file=ofid)
        niceinput=False
       
        
    try:
        em=int(inp['input']['endyr'])
        if em<1970:
            print('End year year must be format YYYY, no earlier than 1970.',file=ofid)
            niceinput=False
            
    
    except ValueError:
        print('End year must be an integer, format YYYY.',file=ofid)
        niceinput=False
        
    try:
        ml=int(inp['quality']['min_length_in_sec'])
    except ValueError:
        print('min_length_in_sec must be a number.',file=ofid)
        niceinput=False
       
    if inp['first_step']!='split' and inp['first_step']!='decimate':
        print('Invalid choice for first_step: Must be split or decimate.',file=ofid)
        niceinput=False
        
    
    if inp['processing']['split']['doit']!='0' and inp['processing']['split']['doit']!='1':
        print('Choice for split must be 0 or 1.',file=ofid)
        niceinput=False
    
    try:
        int(inp['processing']['split']['length_in_sec'])
    except ValueError:
        print('length_in_sec must be a number.',file=ofid)
        niceinput=False
        
        
    if inp['processing']['detrend']!='0' and inp['processing']['detrend']!='1':
        print('Choice for detrend must be 0 or 1.',file=ofid)
        niceinput=False    
        
        
    if inp['processing']['demean']!='0' and inp['processing']['demean']!='1':
        print('Choice for demean must be 0 or 1.',file=ofid)
        niceinput=False

     
    if inp['processing']['trim']!='0' and inp['processing']['trim']!='1':
        print('Choice for trim must be 0 or 1.',file=ofid)
        niceinput=False  
    
    if inp['processing']['taper']['doit']!='0' and inp['processing']['taper']['doit']!='1':
        print('Choice for taper must be 0 or 1.',file=ofid)
        niceinput=False
    
    try:
        tw=float(inp['processing']['taper']['taper_width'])
    except ValueError:
        print('taper width must be a number.',file=ofid)
        niceinput=False

    
    if inp['processing']['decimation']['doit']!='0' and inp['processing']['decimation']['doit']!='1':
        print('Choice for decimation must be 0 or 1.',file=ofid)
        niceinput=False
        
        
    try:
        fs_new=inp['processing']['decimation']['new_sampling_rate'].split(' ')
        for fs in fs_new:
            fs=float(fs)
    except ValueError:
        print('Sampling frequency must be number.',file=ofid)
    
    return niceinput
