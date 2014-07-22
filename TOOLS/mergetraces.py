from __future__ import print_function
from obspy import Stream

def mergetraces(data,Fs,maxgap=10.0,ofid=None):
    
    """
    Small script to merge traces over short gaps. Gaps are filled with zeroes.
    Intended to close short gaps of a couple of seconds.
    
    data: obspy stream, where all the traces you would like to merge are collected already
    Fs: List of original sampling rates; it is checked whether the traces have some weird deviating sampling rate
    maxgap: float, maximum length of gap to be filled with zeros. Default is 10 seconds.
    Every gap longer than this value will remain a gap.
    
    
    """
    
    #Clean up
    data._cleanup()
    if len(data)==1:
        data[0].stats.sampling_rate = round(data[0].stats.sampling_rate,6)
        return data

    
    
        
    else:
        i=0
        newstream=Stream()
        while i<(len(data)-1):
            data[i].stats.sampling_rate = round(data[i].stats.sampling_rate,6)
            # Check the sampling rate
            if data[i].stats.sampling_rate not in Fs:
                print('Bad sampling rate: '+str(Fs),file=ofid)
                i+=1
                continue
           
            # Check if IDs match
            if data[i+1].id!=data[i].id:
                # Append trace i, go to the next trace
                newstream+=data[i]
                i+=1
            else:
                # Check how long the gap is
                if data[i+1].stats.starttime-data[i].stats.endtime<=maxgap:
                # Merge traces
                    s=Stream()
                    try:
                        s=data[i:i+2].merge(method=1,fill_value=0,interpolation_samples=0)
                        newstream+=s
                    except:
                        newstream+=data[i]
                    i+=1
                
                else:
                     # Append trace i and go to the next trace
                    newstream+=data[i]
                    i+=1
                    
        # remove the overlapping segments
        newstream._cleanup()
    
    return newstream

            
    
