import numpy as np
import os

import INPUT.input_kernels as inp

from glob import glob
from scipy.interpolate import griddata
from obspy import read
from mpi4py import MPI

if __name__ == '__main__':
    import noise_source_kernel as nsk
    comm  = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    nsk.obtain_gradient(rank, size)
    
def obtain_kernel(rank, size):
    kernelfiles = glob(inp.kernelfiles)
    mykernelfiles = kernelfiles[rank:len(pairs):size]
    
    # For all 'proto kernel' files:
    for kfile in mykernelfiles:
    # Figure out where the adjoint source is located
        k = kfile.split('/')[-1]
        k = k.replace('.pKern.npy','')
        id1 = k.split('--')[0]
        id2 = k.split('--')[1]
        id = id1 + '.' +inp.channel +'--'+\
        id2 +'.'+inp.channel
        #id = id1 +'--'+id2 # added for synthetic tests
        print id
        adjfile = inp.adj_dir + id +'.'+inp.msr_type+'.SAC'
        if not os.path.exists(adjfile):
            print '+++ Could not find adjoint source for ', id
            continue
        print kfile
        print adjfile
    # Calculate the kernel
        x,y,K = integrate_kernel_sources(kfile,adjfile,inp.mask,inp.f0,inp.f1)

    # Save the plot and the kernel
        kfilenew = inp.outdir + id + '.' + str(inp.f0) + '.' + \
        str(inp.f1) +'.'+inp.msr_type +'.kernel.npy'
        np.save(kfilenew,np.vstack((x,y,K)))
        
        if inp.cap is not None:
            cap = inp.cap * np.max(np.max(np.abs(np.real(K))))
            K = np.clip(K,-cap,cap)
    # Add the kernel to the gradient
        if 'K_tot' in locals():
            K_tot += K 
        else:
            K_tot = K 
    # Collect the misfit (somehow)
    # Save the kernel
    kfilenew_t = inp.outdir + str(inp.f0) + '.' + \
    str(inp.f1) +'.'+inp.msr_type+'.grad'+str(rank)+'.npy'
    np.save(kfilenew_t,np.vstack((x,y,K_tot)))

def integrate_kernel_sources(kernelfile,adjstf,mask,f0,f1,Fs=None):
    #===============================================================================
    # Green's functions
#===============================================================================

    # load kernel
    K = np.load(kernelfile)
    
    # load source location coordinates
    mask_y = np.real(K[0,1:])
    mask_x = np.real(K[1,1:])
    
    # load kernel frequencies
    freq_kern = np.real(K[2:,0])
    
    
    
    if f0 < freq_kern[0]:
        print 'Low frequency outside available range, reset lowest frequency to '\
        +str(freq_kern[0])
        f0 = freq_kern[0]
    
    if f1 > freq_kern[-1]:
        print 'High frequency outside available range, reset highest frequency to '\
        +str(freq_kern[-1])
        f1 = freq_kern[-1]
    
    #===============================================================================
    # Adjoint source
#===============================================================================

    try:
        f = np.load(adjstf)
        if Fs == None:
            msg = 'Fs has to be set when using npy file for adjoint source.'
            raise ValueError(msg)
    except IOError:
        f = read(adjstf)[0]
        Fs = f.stats.sampling_rate
        f = f.data
    
    
    # Transform adjoint source time function to freq. domain
    F = np.fft.rfft(f)
    # Testing: Try smoothing F
    
    
    # Get also the frequency axis.
    if len(f)%2 == 0:
        n = len(f)/2 + 1
    else:
        n = (len(f)+1)/2
        
    freq = np.fft.fftfreq(len(f),d = 1./Fs)
    freq = freq[0:n]
    dfreq = abs(freq[1]-freq[0])
    
    
    # find the frequencies relevant to K
    ind_f0 = np.abs(freq-f0).argmin()
    ind_f1 = np.abs(freq-f1).argmin()
    ind_f0_kern = np.abs(freq_kern-f0).argmin()
    ind_f1_kern = np.abs(freq_kern-f1).argmin()
    
    if ind_f1_kern-ind_f0_kern != ind_f1-ind_f0:
        msg = 'Something went wrong. Are you sure these\
         Green\'s functions and this adjoint source belong together?'
        raise IndexError(msg)
        
    # Allocate K_int
    K_int = np.zeros(len(K[0,:])-1,dtype='complex')
    # multiply K*f
    i = 0

    while freq[ind_f0+i] < freq[ind_f1]:
        print np.real(K[ind_f0_kern+i+2,0])
        print freq[ind_f0+i]
        print '-------------------------'
        #K[ind_f0_kern+i+2,1:] *= F[ind_f0+i]
        

    # Sum up
        K_int += K[ind_f0_kern+i+2,1:] * F[ind_f0+i] * dfreq
        i+=1

    #return kernel
    return mask_x,mask_y,K_int

