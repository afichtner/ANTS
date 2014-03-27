
from math import ceil,  log

def wl_adjust(win_len, Fs):
    """Get the nearest power of two in terms of samples; then determine the corresponding window length in seconds. 
    win_len: Integer, User-defined window length in seconds
    Fs: Integer, Sampling rate """
    #current window length
    cwl=Fs*win_len;
    nwl=2**ceil(log(cwl)/log(2))
    
    #convert the new window length to seconds
    win_len=nwl/Fs
    
    
    return win_len

    
