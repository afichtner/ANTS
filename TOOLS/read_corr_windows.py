import numpy as np
import os

def read_corr_windows(inputfile,size_of_float=4,nbytes_stringhead=256,nbytes_windowname=24):
    """
    Intermediate correlation windows from ant_corr are saved in a specific binary format.
    The header contains: sampling rate (4 byte), number of samples per trace (4 byte), 
    number of subwindows per substack (4 byte), endianness (256 byte), preprocessing string(256 byte)
    Then the data follow in 4-byte floats.
    """
    
    f_in = open(inputfile,'rb')
    
    Fs = np.fromfile(f_in,dtype='f'+str(size_of_float),count=1)[0]
    npts = np.fromfile(f_in,dtype='f'+str(size_of_float),count=1)[0]
    nsub = np.fromfile(f_in,dtype='f'+str(size_of_float),count=1)[0]
    
    endianness = f_in.read(nbytes_stringhead)
    endianness = str(endianness.decode('utf-8')).strip()
    preproc = f_in.read(nbytes_stringhead)
    preproc = str(preproc.decode('utf-8')).strip()
    
    nbytes_header = nbytes_stringhead * 2 + size_of_float * 3
    ntraces = (os.path.getsize(inputfile)-nbytes_header)/\
    (size_of_float*npts+nbytes_windowname)
    
    print "This file contains %g traces of sampling rate %g Hz, each of which contains a stack of %g subtrace(s). \
The byte order is %s-endian and the following preprocessing operations have been applied prior to correlation: %s."\
%(ntraces,Fs,nsub,endianness,preproc)
    
    traces = np.zeros((ntraces,npts))
    
    i = 0
    winnames = dict()
    while i < ntraces:
        winname = f_in.read(nbytes_windowname)
        winname = str(winname.decode('utf-8')).strip()
        winnames[winname] = i
        traces[i,:] += np.fromfile(f_in,dtype='f'+str(size_of_float),count=npts)
        i += 1
        
    return winnames, traces