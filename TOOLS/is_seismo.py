def nice_seismogram(data, threshold=100):
    
    
    """
    A small function to determine wheter a time series looks like a seismogram in the sense that it has many local extrema (ups and downs)
    data: Time series to be checked, a numpy array
    threshold (integer): This many local extrema are needed for the time series to qualify as nice seismogram
    returns True if number of extrema 
    """
    cnt=0
    threshold=int(threshold)
    
    for i in range(1, len(data)-2):
        
        if data[i]<data[i-1] and data[i]<data[i+1]:
            cnt+=1
        elif data[i]>data[i-1] and data[i]>data[i+1]:
            cnt+=1
            
        if cnt>=threshold:
            return True
    
    return False
    
