import os as os
import numpy as np

from obspy.core import read
from obspy.signal import bandpass
from scipy.interpolate import interp1d
from obspy.core import Trace, Stream, UTCDateTime
from scipy.signal import cheby2,  cheb2ord,  filtfilt


#==================================================================================================
# SPLIT TRACES INTO SHORTER SEGMENTS
#==================================================================================================

def split_traces(s, length_in_sec, min_len, verbose, ofid=None):
    """
    Split an ObsPy stream object with multiple traces into a stream with traces of a predefined
    maximum length.
    """

    #- compute length of segments that produces a number of samples close to 2**n -----------------

    power=np.floor(np.log2(length_in_sec*s[0].stats.sampling_rate))
    n_samples=2**power

    new_length_in_sec=n_samples*s[0].stats.delta

    #- print information --------------------------------------------------------------------------

    if verbose==True:
        if ofid==None:
            print '* split into traces of '+str(new_length_in_sec)+' s length'
        else:
            ofid.write('* split into traces of '+str(new_length_in_sec)+' s length\n')
            
    s_new=Stream()
    
    #- loop through traces ------------------------------------------------------------------------
    
    for k in np.arange(len(s)):
           
        #- set initial start time
        start=s[k].stats.starttime
        
        #- march through the trace until the endtime is reached
        while start<s[k].stats.endtime-min_len:
            s_copy=s[k].copy()
            s_copy.trim(start,start+new_length_in_sec-1/(s[k].stats.sampling_rate))
            
            if (s_copy.stats.endtime-s_copy.stats.starttime)<min_len:
                if verbose==True:
                    if ofid==None:
                        print '** trace too short, discarded'
                    else:
                        ofid.write('** trace too short, discarded\n')
            else:
                s_new.append(s_copy)
           
            
            start+=length_in_sec

    return s_new


#==================================================================================================
# TAPER
#==================================================================================================

def taper(data,width,verbose, ofid=None):

    if verbose==True:
        if ofid==None:
            print '* taper '+str(100*width)+' percent of trace'
        else:
            ofid.write('* taper '+str(100*width)+' percent of trace\n')
    data.taper('cosine',p=width)
    
    return data


#==================================================================================================
# DETREND
#==================================================================================================

def detrend(data,verbose, ofid=None):

    if verbose==True:
        if ofid==None:
            print '* detrend'
        else:
            ofid.write('* detrend\n')
        
    data.detrend('linear')
    
    return data
    
#==================================================================================================
# DEMEAN
#==================================================================================================

def demean(data,verbose, ofid=None):

    if verbose==True: 
        if ofid==None:
            print '* demean'
        else:
            ofid.write('* demean\n')
    data.detrend('demean')

    return data

#==================================================================================================
# BANDPASS FILTER
#==================================================================================================

def bandpass(data,corners,f_min,f_max,verbose, ofid=None):
    
    if verbose==True: 
        if ofid==None:
            print '* bandpass between '+str(f_min)+' and '+str(f_max)+' Hz'
        else:
            ofid.write('* bandpass between '+str(f_min)+' and '+str(f_max)+' Hz\n')
    data.filter('bandpass',freqmin=f_min, freqmax=f_max,corners=corners,zerophase=False)
    
    return data

#==================================================================================================
# CHEBYCHEFF LOWPASS FILTER
#==================================================================================================

def antialias(data, freqmax, verbose, ofid=None):
    if isinstance(data,Trace):
        zerophase_chebychev_lowpass_filter(data, freqmax, verbose)
    elif isinstance(data,Stream):
        for trace in data:
            zerophase_chebychev_lowpass_filter(trace, freqmax, verbose)
    
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
    if verbose: 
        if ofid==None:
            print '* Applied low-pass Chebychev filter with corner freq. ', freqmax
        else:
            ofid.write('* Applied low-pass Chebychev filter with corner freq. '+str(freqmax)+'\n')

#==================================================================================================
# BUTTERWORTH LOWPASS FILTER
#==================================================================================================

def lowpass(data,corners,f_max,verbose, ofid=None):

    if verbose==True:
        if ofid==None:
            print '* lowpass below '+str(f_max)+' Hz'
        else:
            ofid.write('* lowpass below '+str(f_max)+' Hz\n')
    data.filter('lowpass', freq=f_max,corners=corners,zerophase=False)
    

    return data

#==================================================================================================
# REMOVE INSTRUMENT RESPONSE
#==================================================================================================

def remove_response(data,respdir,unit,waterlevel,verbose, ofid=None):

    """
    Remove instrument response located in respdir from data. Unit is displacement (DIS), velocity (VEL) or acceleration (ACC).

    Return 1 if successful. Return 0 otherwise.
    """

    #- RESP file ==================================================================================

    resp_file=respdir+'/RESP.'+data.stats.network+'.'+data.stats.station+'.'+data.stats.location+'.'+data.stats.channel

    if verbose==True:
        if ofid==None:
            print '* RESP file: '+resp_file
        else:
            ofid.write('* RESP file: '+resp_file+'\n')

    #- try to remove response if the RESP file exists =============================================

    if os.path.exists(resp_file):

        success=1

        if verbose==True:
            if ofid==None:
                print '* remove instrument response, unit='+unit
            else:
                ofid.write('* remove instrument response, unit='+unit+'\n')
                    
        resp_dict = {"filename": resp_file, "units": unit, "date": data.stats.starttime}

        try:
            data.simulate(seedresp=resp_dict, water_level=float(waterlevel),nfft_pow2=True)
        except ValueError:
            if verbose==True: 
                if ofid==None:
                    print '** could not remove instrument response'
                else:
                    ofid.write('** could not remove instrument response\n')
            success=0

    #- response cannot be removed because RESP file does not exist

    else:
        if verbose==True: 
            if ofid==None:
                print '** could not find correct RESP file'
            else:
                ofid.write('** could not find correct RESP file\n')
            
        success=0

    return success, data



#==================================================================================================
# REMOVE INSTRUMENT RESPONSE
#==================================================================================================

#def ic_sac(data,respdir,unit,verbose):
#    import subprocess
#
#   """
#   Remove instrument response located in respdir from data. Unit is displacement (DIS), velocity (VEL) or acceleration (ACC).
#
#   Return 1 if successful. Return 0 otherwise.
#   """
#
#   #- RESP file ==================================================================================
#   resp_file=respdir+'/RESP.'+data.stats.network+'.'+data.stats.station+'.'+data.stats.location+'.'+data.stats.channel
#
#   if verbose==True: print '* RESP file: '+resp_file
#
#   #- try to remove response if the RESP file exists =============================================
#
#   if os.path.exists(resp_file):
#
#       success=1
#
#       if verbose==True: print '* remove instrument response, unit='+unit
#
#        p = subprocess.Popen(['sac'],
#                             stdout = subprocess.PIPE,
#                             stdin  = subprocess.PIPE,
#                             stderr = subprocess.STDOUT )
#                             
#        s = \
#        'setbb resp ../Resp/' + resp_file.split('/')[-1] + '\n' + \
#        'read ../BH_RAW/' + trace.split('/')[-1] + '\n' + \
#        'rtrend' + '\n' + \
#        'taper' + '\n' + \
#        'rmean' + '\n' + \
#        'trans from evalresp fname %resp to ' + unit_sac + ' freqlim ' + freqlim + '\n' + \
#        'write ' + unit.lower() + '.' + trace_info[1] + '.' + trace_info[2] + \
#                                            '.' + trace_info[3] + '\n' + \
#        'quit\n'
#
#
#   #- response cannot be removed because RESP file does not exist
#
#   else:
#       if verbose==True: print '** could not find correct RESP file'
#       success=0
#
#   return success, data
#



#==================================================================================================
# TRIM TO NEXT FULL SECOND
#==================================================================================================

def trim_next_sec(data,verbose, ofid=None):
    
    """ 
    Trim data to the full second and add a little buffer. Ensures that recordings start and end with a full second.
    A little 10 s buffer of zeroes is added at both ends. Data should be tapered before this operation.

    data=trim_next_sec(data,verbose)

    data:       Is an obspy stream or trace. The returned stream/trace is a bit shorter.
    verbose:    Talk or not.
    
    """

    if verbose: 
        if ofid==None:
            print '* Trimming to full second.'#' Add small buffer around the edges.'
        else:
            ofid.write('* Trimming to full second.\n')

    if isinstance(data,Trace):
        starttime=data.stats.starttime
        endtime=data.stats.endtime
    
        fullsecondtime_start=starttime.strftime('%Y%m%d%H%M%S')
        fullsecondtime_start=UTCDateTime(fullsecondtime_start)#-10
    
        fullsecondtime_end=endtime.strftime('%Y%m%d%H%M%S')
        fullsecondtime_end=UTCDateTime(fullsecondtime_end)#+10
    
        data.trim(starttime=fullsecondtime_start, pad=True, fill_value=0.0)
        data.trim(endtime=fullsecondtime_end, pad=True, fill_value=0.0)
        
    elif isinstance(data,Stream):
        for k in range(len(data)):
            starttime=data[k].stats.starttime
            endtime=data[k].stats.endtime
    
            fullsecondtime_start=starttime.strftime('%Y%m%d%H%M%S')
            fullsecondtime_start=UTCDateTime(fullsecondtime_start)#-10
        
            fullsecondtime_end=endtime.strftime('%Y%m%d%H%M%S')
            fullsecondtime_end=UTCDateTime(fullsecondtime_end)#+10
        
            data[k].trim(starttime=fullsecondtime_start, pad=True, fill_value=0.0)
            data[k].trim(endtime=fullsecondtime_end, pad=True, fill_value=0.0)

    return data
    


#==================================================================================================
# DOWNSAMPLING
#==================================================================================================

def downsample(data, Fsnew, verbose, ofid=None):

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
    
    Fsnew=float(Fsnew)
    f_max=Fsnew/4
    
    data_new=data.copy()

    #- Check if data already have the desired sampling rate =======================================

    if Fs==Fsnew:
        if verbose==True:
            if ofid==None:
                print '* Current and new sampling rate are equal. No downsampling performed.'
            else:
                ofid.write('* Current and new sampling rate are equal. No downsampling performed.\n')

    #- Downsampling ===============================================================================

    else:
        
        dec=int(Fs/Fsnew)
        data_new.decimate(dec, no_filter=True)
        
        try:
            tl=len(data_new.data)
        except AttributeError:
            tl=len(data_new[0].data)
        
        if verbose==True:
            if ofid==None:
                print '* downsampling by factor '+str(dec)+', length of new trace: '+str(tl)
            else:
                ofid.write('* downsampling by factor '+str(dec)+', length of new trace: '+str(tl)+'\n')

    return data_new
    
    
