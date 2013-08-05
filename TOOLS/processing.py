import os as os
import numpy as np
from obspy.core import read
from scipy.interpolate import interp1d
from obspy.core import trace, stream, UTCDateTime


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

	if verbose==True: print '* demean and detrend'
	data.detrend('linear')
	data.detrend('demean')

	return data


#==================================================================================================
# BANDPASS FILTER
#==================================================================================================

def bandpass(data,corners,f_min,f_max,verbose):

	if verbose==True: print '* bandpass between '+str(f_min)+' and '+str(f_max)+' Hz'
	data.filter('lowpass',freq=f_max,corners=corners,zerophase=False)
	data.filter('highpass',freq=f_min,corners=corners,zerophase=False)

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
# TRIM TO NEXT FULL SECOND
#==================================================================================================

def trim_next_sec(data,verbose):
    
    """ 
    Trim data to the full second and add a little buffer. Ensures that recordings start and end with a full second.
    A little 10 s buffer of zeroes is added at both ends. Data should be tapered before this operation.

  	data=trim_next_sec(data,verbose)

    data: 		Is an obspy stream or trace. The returned stream/trace is a bit shorter.
    verbose:	Talk or not.
    
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

	data:		ObsPy stream with original data.
	Fsnew:		New sampling rate in samples per second.
	verbose:	Talk or not.

	data_new:	Downsampled data.

	"""

	Fs=float(data.stats.sampling_rate) 
	Fsnew=float(Fsnew)

	data_new=data.copy()

	#- Check if data already have the desired sampling rate =======================================

	if Fs==Fsnew:
		if verbose==True: 
			print '* Current and new sampling rate are equal. No downsampling performed.'

	#- Downsampling ===============================================================================

	else:
		dec=int(Fs/Fsnew)
		data_new.decimate(dec, no_filter=True)
		if verbose==True:  print '* downsampling by factor '+str(dec)

	return data_new
    
    
