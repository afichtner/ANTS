import os as os
import numpy as np
from obspy.core import read
from obspy.core import stream
from obspy.signal import bandpass
from scipy.interpolate import interp1d
from obspy.core import trace, stream, UTCDateTime

#==================================================================================================
# SPLIT TRACES INTO SHORTER SEGMENTS
#==================================================================================================

def split_traces(s,length_in_sec,min_length_in_sec,verbose):
    """
    Split an ObsPy stream object with multiple traces into a stream with traces of a predefined
    maximum length.
    """

    if verbose==True: print '* split into traces of '+str(length_in_sec)+' s length'

    s_new=stream.Stream()

    #- loop through traces ------------------------------------------------------------------------
    for k in np.arange(len(s)):
        #- set initial start time
        start=s[k].stats.starttime
        #- march through the trace until the endtime is reached
        while start<s[k].stats.endtime:
            s_copy=s.copy()
            s_copy.trim(start,start+length_in_sec)

            start+=length_in_sec
            if (s_copy[0].stats.endtime-s_copy[0].stats.starttime)<min_length_in_sec:
                if verbose==True: print '** trace too short, discarded'
                continue
            else:
                s_new.append(s_copy[0])

    return s_new

#==================================================================================================
# TAPER
#==================================================================================================

def taper(data,width,verbose):

    if verbose==True: print '* taper '+str(100*width)+' percent of trace'
    data.taper('cosine',p=width)

    return data


#==================================================================================================
# DETREND
#==================================================================================================

def detrend(data,verbose):

    if verbose==True: print '* detrend'
    data.detrend('linear')

    return data
    
#==================================================================================================
# DEMEAN
#==================================================================================================

def demean(data,verbose):

    if verbose==True: print '* demean'
    data.detrend('demean')

    return data

#==================================================================================================
# BANDPASS FILTER
#==================================================================================================

def bandpass(data,corners,f_min,f_max,verbose):

    if verbose==True: print '* bandpass between '+str(f_min)+' and '+str(f_max)+' Hz'
    data.filter('bandpass',freqmin=f_min, freqmax=f_max,corners=corners,zerophase=False)
    

    return data

#==================================================================================================
# LOWPASS FILTER
#==================================================================================================

def lowpass(data,corners,f_max,verbose):

    if verbose==True: print '* lowpass below '+str(f_max)+' Hz'
    data.filter('lowpass', freqmax=f_max,corners=corners,zerophase=False)
    

    return data


#==================================================================================================
# REMOVE INSTRUMENT RESPONSE
#==================================================================================================

def remove_response(data,respdir,unit,verbose):

    """
    Remove instrument response located in respdir from data. Unit is displacement (DIS), velocity (VEL) or acceleration (ACC).

    Return 1 if successful. Return 0 otherwise.
    """

    #- RESP file ==================================================================================
    resp_file=respdir+'/RESP.'+data.stats.network+'.'+data.stats.station+'.'+data.stats.location+'.'+data.stats.channel

    if verbose==True: print '* RESP file: '+resp_file

    #- try to remove response if the RESP file exists =============================================

    if os.path.exists(resp_file):

        success=1

        if verbose==True: print '* remove instrument response, unit='+unit
        resp_dict = {"filename": resp_file, "units": unit, "date": data.stats.starttime}

        try:
            data.simulate(seedresp=resp_dict)
        except ValueError:
            if verbose==True: print '** could not remove instrument response'
            success=0

    #- response cannot be removed because RESP file does not exist

    else:
        if verbose==True: print '** could not find correct RESP file'
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

def trim_next_sec(data,verbose):
    
    """ 
    Trim data to the full second and add a little buffer. Ensures that recordings start and end with a full second.
    A little 10 s buffer of zeroes is added at both ends. Data should be tapered before this operation.

    data=trim_next_sec(data,verbose)

    data:       Is an obspy stream or trace. The returned stream/trace is a bit shorter.
    verbose:    Talk or not.
    
    """

    if verbose: print '* Trimming to full second. Add small buffer around the edges.'

    try:
        starttime=data.stats.starttime
        endtime=data.stats.endtime
    except AttributeError:        
        starttime=data[0].stats.starttime
        endtime=data[0].stats.endtime
    
    fullsecondtime_start=starttime.strftime('%Y%m%d%H%M%S')
    fullsecondtime_start=UTCDateTime(fullsecondtime_start)-10

    fullsecondtime_end=endtime.strftime('%Y%m%d%H%M%S')
    fullsecondtime_end=UTCDateTime(fullsecondtime_end)+10

    data.trim(starttime=fullsecondtime_start, pad=True, fill_value=0.0)
    data.trim(endtime=fullsecondtime_end, pad=True, fill_value=0.0)

    return data
    


#==================================================================================================
# DOWNSAMPLING AND INTERPOLATION
#==================================================================================================

def downsample(data, Fsnew, verbose):

    """
    Downsample data to a new sampling rate.

    data_new=downsample(data, Fsnew, verbose)

    data:       ObsPy stream with original data.
    Fsnew:      New sampling rate in samples per second.
    verbose:    Talk or not.

    data_new:   Downsampled data.
    
    A lowpassfilter is applied to avoid aliasing.

    """

    Fs=float(data.stats.sampling_rate) 
    Fsnew=float(Fsnew)
    f_max=Fsnew/4

    data_new=data.copy()

    #- Check if data already have the desired sampling rate =======================================

    if Fs==Fsnew:
        if verbose==True: 
            print '* Current and new sampling rate are equal. No downsampling performed.'

    #- Downsampling ===============================================================================

    else:
        data_new.filter('lowpass', freq=f_max, corners=4,zerophase=False)
        dec=int(Fs/Fsnew)
        data_new.decimate(dec, no_filter=True)
        if verbose: print '* Lowpass filter with f_c = ', f_max, ' before downsampling.'
        if verbose==True:  print '* downsampling by factor '+str(dec)

    return data_new
    
    
