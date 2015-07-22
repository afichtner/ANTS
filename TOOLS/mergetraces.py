from __future__ import print_function
import numpy as np
from obspy import Stream

def mergetraces(data, Fs, maxgap=10.0, ofid=None):
    
    """
    Small script to merge traces over short gaps. Gaps are filled with zeroes.
    Intended to close short gaps of a couple of seconds.
    
    data: obspy stream, where all the traces you would like to merge are collected already
    Fs: List of original sampling rates; it is checked whether the traces have some weird deviating sampling rate
    maxgap: float, maximum length of gap to be filled with zeros. Default is 10 seconds.
    Every gap longer than this value will remain a gap.
    
    
    """

    # check sampling rates and dtypes
    # Did not run merge Checks: Do not return if just two traces would cause problems, because other traces might still be merged
    # order matters!
    # It is good to have the order functionality here because cleanup would not order if mergeChecks throws exception
    data.sort(keys=['network', 'station', 'location', 'channel', 'starttime', 'endtime'])

    # build up dictionary with lists of traces with same ids
    traces_dict = {}
    # using pop() and try-except saves memory
    try:
        while True:
            # I think this could also directly be done as data.pop(0)
            trace = data.pop(0)
            # Skip empty traces (if any)
            if trace.stats.npts <= 0:
                continue

            trace.stats.sampling_rate = round(trace.stats.sampling_rate, 4)
            # Throw data with the wrong sampling rate out.
            if trace.stats.sampling_rate not in Fs:
                print('Bad sampling rate: '+str(trace.stats.sampling_rate), file=ofid)
                continue

            # add trace to respective list or create that list
            traces_dict.setdefault(trace.id, []).append(trace)
            
    except IndexError:
        pass

    # 'data' contains no traces now; fill it with merged content
    
    # loop through ids
    for id in list(traces_dict.keys()):
        #print('===================')
        
        trace_list = traces_dict[id]
        cur_trace = Stream(trace_list.pop(0))

        # work through all traces of same id
        while trace_list:
            trace = trace_list.pop(0)
            #print('===================')
            print(cur_trace)
            #print('Next trace:')
            print(trace)
            # Case 1: Overlap =================================================
            # Overlaps should be removed in any case, as we don't want to correlate the data twice; especially if they differ between traces (which is most likely the case due to different tapering etc.)
            # Case 2: Perfectly adjacent ======================================
            # Case 3: Short gap ===============================================
            if trace.stats.starttime - cur_trace[-1].stats.endtime <= maxgap:
                #print('Overlap or adjacent or short gap.')
                
                cur_trace += trace
                cur_trace.merge(method=1,interpolation_samples=0,\
                fill_value=0)
            # Case 4: Long gap ================================================
            else:
                #print('Long gap.')
                # Add to output stream
                data += cur_trace
                # Start new stream
                cur_trace = Stream(trace)
            
        # Add the last trace of this id
        # This is actually a stream, but that is ok, adding streams just adds their traces
        data += cur_trace

    return data

