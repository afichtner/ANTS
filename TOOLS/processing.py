from __future__ import print_function
import os as os
import numpy as np

from obspy.core import read

from scipy.interpolate import interp1d
from obspy.core import Trace, Stream, UTCDateTime
from scipy.signal import cheby2,  cheb2ord,  filtfilt
from glob import glob
from gc import collect


#==================================================================================================
# SPLIT TRACES INTO SHORTER SEGMENTS
#==================================================================================================

def split_traces(s, length_in_sec, min_len, verbose, ofid):
    """
    Split an ObsPy stream object with multiple traces into a stream with traces of a predefined
    maximum length.
    """
            
    s_new=Stream()
    
    #- loop through traces ------------------------------------------------------------------------
    
    for k in np.arange(len(s)):
           
        #- set initial start time
        start=s[k].stats.starttime
        
        #- march through the trace until the endtime is reached
        while start<s[k].stats.endtime-min_len:
            s_copy=s[k].copy()
            s_copy.trim(start,start+length_in_sec-1/(s[k].stats.sampling_rate))
            s_new.append(s_copy)
            del s_copy
            collect()
            start+=length_in_sec

    return s_new

#==================================================================================================
# SLICE TRACES (looking for a less memory intensive way....)
#==================================================================================================

def slice_traces(s,length_in_sec, min_len, verbose, ofid):
    """
    Slice an ObsPy stream object with multiple traces; The stream of new (sliced) traces merely contains
    references to the original trace.
    """
    s_new=Stream()
    
    #- loop through traces ------------------------------------------------------------------------
    
    for k in np.arange(len(s)):
           
        #- set initial start time
        start=s[k].stats.starttime
        
        #- march through the trace until the endtime is reached
        while start<s[k].stats.endtime-min_len:
            s_part=s[k].slice(start,start+length_in_sec-1/(s[k].stats.sampling_rate))
        
            s_new.append(s_part)
            
            start+=length_in_sec

    return s_new
    

#==================================================================================================
# TAPER
#==================================================================================================

def taper(data,width,verbose,ofid,ttype='cosine'):
    try:
        data.taper(max_percentage=width/2.0,type=ttype)
        if verbose==True:
            print('* taper '+str(100*width)+' percent of trace\n',file=ofid)
    except TypeError:
        # Means you are on Rothorn
        data.taper()
        if verbose==True: 
            print('* taper with default taper\n',file=ofid)
                
        
    
    
    return data


#==================================================================================================
# DETREND
#==================================================================================================

def detrend(data,verbose, ofid):

    if verbose==True:
        print('* detrend\n',file=ofid)
        
    data.detrend('linear')
    
    return data
    
#==================================================================================================
# DEMEAN
#==================================================================================================

def demean(data,verbose, ofid):

    if verbose==True: 
        print('* demean\n',file=ofid)
        
    data.detrend('demean')

    return data

#==================================================================================================
# BANDPASS FILTER
#==================================================================================================

def bandpass(data,corners,f_min,f_max,verbose, ofid):
    
    if verbose==True: 
        print('* bandpass between '+str(f_min)+' and '+str(f_max)+' Hz\n',file=ofid)
        print ('* filter order '+str(corners)+' \n',file=ofid)
        
    data.filter('bandpass',freqmin=f_min, freqmax=f_max,corners=corners,zerophase=True)
    
    return data

#==================================================================================================
# CHEBYCHEFF LOWPASS FILTER
#==================================================================================================

def antialias(data, freqmax, verbose, ofid=None):
    if isinstance(data,Trace):
        zerophase_chebychev_lowpass_filter(data, freqmax, verbose,ofid)
        if verbose: 
            print('* Applied low-pass Chebychev filter with corner freq. '+str(freqmax)+'\n',file=ofid)
    else:
        msg = 'I only take traces for antialiasing. Regards, your chebycheff filter.'
        raise TypeError(msg)
           
    return data
            

def zerophase_chebychev_lowpass_filter(trace, freqmax, verbose, ofid=None):
    """
    Custom Chebychev type two zerophase lowpass filter useful for decimation
    filtering.

    This filter is stable up to a reduction in frequency with a factor of 10.
    If more reduction is desired, simply decimate in steps.

    Partly based on a filter in ObsPy.

    :param data: Input trace
    :param freqmax: The desired lowpass frequency.

    Will be replaced once ObsPy has a proper decimation filter.
    """ 
        
    # rp - maximum ripple of passband, rs - attenuation of stopband
    rp, rs, order = 1, 96, 1e99
    ws = freqmax / (trace.stats.sampling_rate * 0.5)  # stop band frequency
    wp = ws  # pass band frequency

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = cheb2ord(wp, ws, rp, rs, analog=0)

    b, a = cheby2(order, rs, wn, btype="low", analog=0, output="ba")

    # Apply twice to get rid of the phase distortion.
    trace.data = filtfilt(b, a, trace.data)
    
#==================================================================================================
# BUTTERWORTH LOWPASS FILTER
#==================================================================================================

def lowpass(data,corners,f_max,verbose, ofid):

    if verbose==True:
        print('* lowpass below '+str(f_max)+' Hz\n',file=ofid)
        
    data.filter('lowpass', freq=f_max,corners=corners,zerophase=False)
    

    return data

#==================================================================================================
# REMOVE INSTRUMENT RESPONSE
#==================================================================================================

def remove_response(data,respdir,unit,freqs,waterlevel,verbose, ofid):

    """
    Remove instrument response located in respdir from data. Unit is displacement (DIS), velocity (VEL) or acceleration (ACC).

    Return 1 if successful. Return 0 otherwise.
    """

    #- RESP file ==================================================================================

    resp_file=respdir+'/RESP.'+data.stats.network+'.'+data.stats.station+'.'+data.stats.location+'.'+data.stats.channel

    if verbose==True:
        print('* RESP file: '+resp_file+'\n',file=ofid)

    #- try to remove response if the RESP file exists =============================================

    if os.path.exists(resp_file):

        success=1

        if verbose==True:
            print('* remove instrument response, unit='+unit+'\n',file=ofid)
                    
        resp_dict = {"filename": resp_file, "units": unit, "date": data.stats.starttime}
        
        try:
            data.simulate(seedresp=resp_dict, water_level=float(waterlevel),\
            nfft_pow2=True, simulate_sensitivity=False,pre_filt=tuple(freqs),\
            pitsasim=False,sacsim=True)
        except ValueError:
            if verbose==True: 
                print('** could not remove instrument response\n',file=ofid)
            success=0

    #- response cannot be removed because RESP file does not exist

    else:
        if verbose==True: 
            print('** could not find correct RESP file\n',file=ofid)
            
        success=0

    return success, data



#==================================================================================================
# TRIM TO NEXT FULL SECOND
#==================================================================================================

def trim_next_sec(data,verbose,ofid):
    
    """ 
    Trim data to the full second. Ensures that recordings start and end with a full second.

    data=trim_next_sec(data,verbose)

    data:       Is an obspy stream or trace. The returned stream/trace is a bit shorter.
    verbose:    Talk or not.
    ofid: Output file in 'check' mode; otherwise, 'None' is the default inherited from main.
    
    """

    

    if isinstance(data,Trace):
        sec_to_remove=data.stats.starttime.microsecond/1e6
        sec_to_add=1.-sec_to_remove
        if sec_to_add > data.stats.delta:
            if verbose: 
                    print('* Trimming to full second.\n',file=ofid)
            data.trim(starttime=data.stats.starttime+sec_to_add,nearest_sample=True)
        
    elif isinstance(data,Stream):
        for tr in data:
            starttime=tr.stats.starttime
            sec_to_remove=tr.stats.starttime.microsecond/1e6
            sec_to_add=1.-sec_to_remove
            if sec_to_add > tr.stats.delta:
                if verbose: 
                        print('* Trimming to full second.\n',file=ofid)
                tr.trim(starttime=tr.stats.starttime+sec_to_add,nearest_sample=True)
            
    return data
            
 
        

    


#==================================================================================================
# DOWNSAMPLING
#==================================================================================================

def downsample(data, Fsnew, verbose, ofid):

    """
    Downsample data to a new sampling rate.

    data_new=downsample(data, Fsnew, verbose)

    data:       ObsPy stream with original data.
    Fsnew:      New sampling rate in samples per second.
    verbose:    Talk or not.

    data_new:   Downsampled data.
    
    A lowpassfilter is applied to avoid aliasing.

    """
    try:
        Fs=float(data.stats.sampling_rate)
    except AttributeError:
        Fs=float(data[0].stats.sampling_rate)
    
    Fsnew=np.double(Fsnew)
    

    #- Check if data already have the desired sampling rate =======================================

    if Fs==Fsnew:
        if verbose==True:
            print('* Current and new sampling rate are equal. No downsampling performed.\n',file=ofid)

    #- Downsampling ===============================================================================

    else:
        
        dec=int(Fs/Fsnew)
        data=antialias(data,Fsnew*0.4,verbose,ofid)
        data.decimate(dec, no_filter=True, strict_length = False)
        
        try:
            tl=len(data.data)
        except AttributeError:
            tl=len(data[0].data)
        
        if verbose==True:
            print('* downsampling by factor '+str(dec)+', length of new trace: '+str(tl)+'\n',file=ofid)

    return data
    
    
