import numpy as np
from scipy.signal import hann
from obspy.signal.invsim import cosTaper
from math import ceil
import matplotlib.pyplot as plt

def event_exclude(trace,windows,n_compare,min_freq,\
factor_enrg=1.,taper_perc=0.05,thresh_stdv=1.,undo_taper=True):

    # " Undo " the taper temporarily. Needed so that taper doesn't appear like an energy increase.
    #if undo_taper:
    #    twin = cosTaper(trace.stats.npts,taper_perc*2.)
    #    twin[0] = twin[1]
    #    twin[-1]= twin[-2]
    #    trace.data/=twin
        
        
    weight = np.ones(trace.stats.npts)
    windows.sort() # make sure ascending order
    # length of taper dep on minimum frequency that should be available
    # check this again, not sure about it
    n_hann = int(trace.stats.sampling_rate / min_freq)
    tpr = hann(2*n_hann)
    
    for win in windows: 
# initialize an array of subtrace values (energy, standard deviation) for each window length; maybe just use an expandable list (max. a couple of hundred values)
        enrg = []
        stdv = []
        t = []
        marker = []
        weighting_trace = np.ones(trace.stats.npts)
# fill those arrays
        t0 = trace.stats.starttime
        while t0 < trace.stats.endtime-win:
            
            enrg.append(np.sum(np.power(trace.slice(starttime=t0,endtime=t0+win-1).data,2))/win)
            subwin = int(win/3)
            [a,b,c] = [np.std(trace.slice(starttime=t0,endtime=t0+win-1).data[0:subwin]),
                        np.std(trace.slice(starttime=t0,endtime=t0+win-1).data[subwin:2*subwin]),
                    np.std(trace.slice(starttime=t0,endtime=t0+win-1).data[2*subwin:3*subwin])]  
                              #np.std(trace.slice(starttime=t0,endtime=t0+win-1).data[3*subwin:4*subwin])]
                        #np.std(trace.slice(starttime=t0,endtime=t0+win-1).data[4*subwin:5*subwin])]
            #    
            stdv.append(np.max([a,b,c])/np.min([a,b,c]))
            t.append((t0+win/2).strftime('%s'))
            t0 += win
        
        # count how many windows are excluded on counter
        # step through the array enrg, stdv; this should be relatively fast as array are not particularly long
        winsmp = int(ceil(win*trace.stats.sampling_rate))
         
        for i in range(2,len(enrg)):
            i0 = i - n_compare if i>=n_compare else 0
            i1 = i0 + n_compare 
            if i1 >= len(enrg):
                i0 = i0 - (len(enrg)-i1)
                i1 = len(enrg)
                
            mean_enrg = np.mean(enrg[i0:i1])
            
            if enrg[i] > factor_enrg * mean_enrg: #and stdv[i] > thresh_stdv:
                
                marker.append(1)
                weighting_trace[i*winsmp:(i+1)*winsmp] *= 0.
                if i*winsmp > n_hann:
                    weighting_trace[i*winsmp-n_hann:i*winsmp] *= 1-tpr[0:n_hann]
                else:
                    weighting_trace[0:i*winsmp] *= 1-tpr[n_hann-i*winsmp:n_hann]
                
                if (i+1)*winsmp+n_hann <=len(weighting_trace):
                    weighting_trace[(i+1)*winsmp:(i+1)*winsmp+n_hann] *= 1-tpr[n_hann:]
                else:
                    weighting_trace[(i+1)*winsmp:] *= 1-tpr[n_hann:n_hann+len(weighting_trace)-(i+1)*winsmp]
                # build in that if taper is longer than trace itself, it gets shortened.
            else:
                marker.append(0)
        #plt.plot(trace.times()+int(trace.stats.starttime.strftime('%s')),weighting_trace)
        weight *= weighting_trace        
        marker=np.array(marker)
        enrg=np.array(enrg)
        

# plot the trace; original and processed; plot only part if trace v. long

   
    #plt.plot(trace.times()+int(trace.stats.starttime.strftime('%s')),trace.data,color='0.5')

    trace.data *= weight
    # percentage of data that was cut:
    pctg_cut=float(np.sum(weight==0.))/float(trace.stats.npts)*100
    # display a summary of how much was kept 
    print('cut %g percent of data from trace: ' %pctg_cut)
    print trace.id
    
    #plt.plot(trace.times()+int(trace.stats.starttime.strftime('%s')),trace.data,'k')
 #   plt.show()
#    
    #spec_excl = np.fft.fft(trace.data)
    #f = np.fft.fftfreq(trace.stats.npts,d=trace.stats.delta)
    
    #plt.plot(f,spec)
    #plt.plot(f,spec_excl) #plt.plot(trace.times()+int(trace.stats.starttime.strftime('%s')),weight,linewidth=2.)
  #  plt.show()


# don't have to return trace -- operates directly on trace.
    return()
    
    
# discarded:

        #stdv=np.array(stdv) 
        #plt.plot(t,enrg/np.max(np.abs(enrg)))
        #plt.plot(t,stdv/np.max(np.abs(stdv)),'--') 
        #plt.plot(t,marker,'d')
            
            
            #stdv.append(np.max(np.abs(trace.slice(starttime=t0,endtime=t0+win-1).data)))        
# walk through making decisions about keeping or throwing out
# Throw out by setting respective window in weighting_trace to
# 1 - Hann[0:half window],0,1-Hann[half window to end] (or similar; look at window)
        


