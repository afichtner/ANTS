from obspy import Stream

def mergetraces(data,maxgap=10.0):
    
    """
    Small script to merge traces over short gaps. Gaps are filled with zeroes.
    Intended to close short gaps of a couple of seconds.
    
    data: obspy stream, where all the traces you would like to merge are collected already
    maxgap: float, maximum length of gap to be filled with zeros. Default is 10 seconds.
    Every gap longer than this value will remain a gap.
    
    
    """
    
    #Clean up
    
    data._cleanup()

    i=0
    
    newstream=Stream()
    
    while i<(len(data)-1):
        
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
                s=data[i:i+2].merge(method=1,fill_value=0,interpolation_samples=0)
                newstream+=s
                i+=1
            
            else:
            # Do nothing and go to the next trace
                newstream+=data[i]
                i+=1
                
    # remove the overlapping segments
    newstream._cleanup()
    
    return newstream
                
            
            
    
