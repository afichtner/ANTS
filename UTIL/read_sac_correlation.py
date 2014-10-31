from __future__ import print_function

from obspy import read, Trace
from obspy.core.trace import Stats
from obspy.core import AttribDict
import numpy as np
from obspy.noise.correlation import Correlation, CorrelationStream

from glob import glob
import fnmatch
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
font={'size':14,  'family':'normal'}
mpl.rc('font', **font)



def read_corr_sac(filepath,corrtype,verbose=False):
    
    fls = glob(filepath)
    
    cstream = CorrelationStream()
    
    for file in fls:
        corrtype=file.split('/')[-1].split('.')[8]
        # The tag should become obsolete here once the correlation 
        # metainformation is all contained in the correlation file.
        tag = file.split('/')[-1].split('.')[9]
        
        stats_a=Stats()
        stats_b=Stats()
        
        tr = read(file)[0]
        st = tr.stats
        
        stats_a = st
        stats_a.network = str(stats_a.network)
        stats_a.station = str(stats_a.station)
        stats_a.location = str(stats_a.location)
        stats_a.channel = str(stats_a.channel)
        
        # hide the tag somewhere
        stats_a.sac['kt2'] = tag
        
        stats_b.network=str(st.sac['kuser0']).strip()
        stats_b.station=str(st.sac['kevnm']).strip()
        stats_b.location=str(st.sac['kuser1']).strip()
        stats_b.channel=str(st.sac['kuser2']).strip()
        stats_b.sampling_rate=float(st.sampling_rate)
        
        if str(stats_b.location)=='-12345':
            stats_b.location=''

        
        corr = Correlation(stats_a,stats_b,tr.data,max_lag=st.sac['e'],\
                                correlation_type=corrtype,min_lag=0,n_stack=stats_a.sac['user0'])
        if verbose==True:
            print(corr)
        cstream+=corr
    return cstream


def write_corr_sac(corrstack,corrtype,stackname,outdir,verbose):
#==============================================================================
    #- Write metadata info to sac header
    #- Store results
#==============================================================================
    print(corrstack.id)
    id = corrstack.id.split('--')[0]+'.'+corrstack.id.split('--')[1]
    #- open file and write correlation function
    filename=outdir+id+'.'+corrtype+'.'+stackname+'.SAC'
    
    if os.path.exists(filename):
        msg = 'Choose a new directory or name to save, otherwise old data ' +\
                'would be overwritten.'
        raise ValueError(msg)
        
    tr = Trace(data=corrstack.correlation)
    tr.stats = corrstack.stats_a
    tr.write(filename,format='SAC')
    



def plot_sect(cstream,prefilter=None,annotate=True,\
              station1=None,station2=None,exagg=1.,save_as=None):
#==============================================================================
    #- Plot the traces contained in the correlation stream
    #- offset = interstation distance
#==============================================================================
              
    fig = plt.figure(figsize=(10, 10))
    fig.hold()
    plt.subplot(111)
    
    ids = []
    
    
    for i in range(len(cstream.select(station1,station2))):
        
        corr = cstream[i]
        
        net1 = corr.stats_a.network
        sta1 = corr.stats_a.station
        net2 = corr.stats_b.network
        sta2 = corr.stats_b.station
        
        id = net1+'.'+sta1+'-'+\
             net2+'.'+sta2
        if id not in ids:
            ids.append(id)
            
            data = corr.correlation.copy()
            
            dist=float(corr.stats_a.sac['dist'])/10000
            lag = corr.stats_a.sac['e']
            
            data /= corr.stats_a.sac['user0']
            
            if prefilter is not None:
                tr = Trace(data=data)
                tr.stats.sampling_rate = corr.stats_a.sampling_rate
                tr.filter('bandpass',freqmin = prefilter[0],\
                            freqmax = prefilter[1], corners = prefilter[2],\
                                zerophase = True)
                data = tr.data
            
            
            wcount=' ('+str(int(corr.stats_a.sac['user0']))+' windows)'
            
            plt.plot(np.linspace(-lag,lag,len(corr.correlation)),\
                exagg*data+dist,'k')
            
            plt.annotate(id, xy=(lag, dist+1),xytext=(-lag-3000,dist+1),\
                         fontsize=12,color='k')
            plt.annotate(wcount, xy=(lag, dist+0.5),xytext=(lag,dist),\
                         fontsize=12,color='k')
    plt.xlabel('Lag (s)')
    plt.ylabel('Station-station dist (1000 km)')
    plt.ylim([0,2000])
    
    if save_as is not None:
        plt.savefig(save_as,format='eps')
    plt.show()
    

#rrelationStream(object):
#
#__init__(self,correlations=None,cor_typ='ccc',max_lag=None,win_len=None,\
#         overlap=None,Fs=None,bandpass=None):
#
#self.__correlations = []
#
#if isinstance(correlations, Correlation):
#    self.__correlations.append(correlations)
## This thing with the list of correlation objects....not sure.
#if isinstance(correlations,list):
#    self.__correlations.extend(correlations)
# 
##self.cstats=AttribDict()
##self.cstats.update({'cor_typ':cor_typ})
##self.cstats.update({'max_lag':max_lag})
##self.cstats.update({'win_len':win_len})
##self.cstats.update({'overlap':overlap})
##self.cstats.update({'Fs':Fs})
##self.cstats.update({'bandpass':bandpass})
#
#__add__(self,other):
#if isinstance(other, Correlation):
#    other=CorrelationStack(correlations=other)
#if not isintance(other, CorrelationStack):
#    raise TypeError
#correlations = self.__correlations + other.__correlations
#return self.__class__(correlations=correlations)
#
#__iadd__(self,other):
#if isinstance(other, Correlation):
#    other=CorrelationStack(correlations=other)
#if not isintance(other, CorrelationStack):
#    raise TypeError
#self.__correlations += other.__correlations
#
# __iter__(self):
# return list(self.__correlations).__iter__()
#__getitem__(self, index):
#    """
#    __getitem__ method
#    :return: Correlation objects
#    """
#    if isinstance(index, slice):
#        return self.__class__(correlations=\
#                self.__correlations.__getitem__(index))
#    else:
#        return self.__correlations.__getitem__(index)
# 
#__str__(self,extended=False):
#out = 'Contains %g correlation(s):\n' %len(self.__correlations)
##Don't print for all, it s too long!
#if len(self.__correlations) <= 10 or extended is True:
#    out = out + "\n".join([_i.__str__() for _i in self])
#else: 
#    out = out + "\n" + self.__correlations[0].__str__() + "\n" + \
#        '...\n(%i other correlations)\n...\n' % (len(self.__correlations) - 2) + \
#        self.__correlations[-1].__str__() + '\n\n[Use "print(' + \
#        'CorrelationStack.__str__(extended=True))" to print all correlations]'
#return out
#
# bandpassfilter(self,filter=None):
# 
# for corr in self.__correlations:
#append(self,other):
#if isinstance(other, Correlation):
#    self.__correlations.append(other)
#
#stack(self,tr1,tr2):
#
##Check if necessary input is available
#if self.cstats['win_len']==None:
#    msg='No window length specified.'
#    raise ValueError(msg)
#
## Check if these are traces
#if not isinstance(tr1,Trace) or not isinstance(tr2,Trace):
#    msg='Input must be obspy Trace object.'
#    raise TypeError(msg)
#
## Check sampling rate
#if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
#    msg='Sampling rates must be the same for input traces.'
#    raise ValueError(msg)
#    
#if self.cstats['Fs'] is None:
#    Fs=tr1.stats.sampling_rate
#else:
#    Fs=self.cstats['Fs']
#    
#if tr1.stats.sampling_rate != Fs or\
#    tr2.stats.sampling_rate != Fs:
#    msg='Incompatible sampling rate.'
#    raise ValueError(msg) 
#
## Check trace length
#if len(tr1.data) <= 2*Fs*self.cstats['max_lag'] or\
#    len(tr2.data) <= 2*Fs*self.cstats['max_lag']:
#    msg='Trace(s) is/are too short.'
#    raise ValueError(msg)
#
#
#
## Cut the trace in windows
## Window lenght in seconds
#wlen = int(self.cstats['win_len'])
#t1=max(tr1.stats.starttime,tr2.stats.starttime)
#t2=t1+wlen
#stack=list()
#
#while t2<min(tr1.stats.endtime,tr2.stats.endtime):
#    
#    tr_a=tr1.slice(starttime=t1,endtime=t2)
#    tr_b=tr2.slice(starttime=t1,endtime=t2)
#    
#    # Correlate them
#    corr = correlate_trace(tr_a,tr_b,self.cstats['max_lag'],\
#                           self.cstats['cor_typ'])
#    corr.stats_a.sac={}
#    corr.stats_a.sac['user0']=1
#    corr.stats_a.sac['b']=-self.cstats['max_lag']
#    corr.stats_a.sac['e']=self.cstats['max_lag']
#    
#    stack.append(corr)
#    
#    # Update window
#    t1=t2
#    t2+=wlen
#    
#return CorrelationStack(correlations=stack)
#
#
#save_linstack_sac(self,filename=None):
#
##if num_corrs and num_wins:
##    msg('Window and correlation number are mutually '+\
##        'exclusive stack length designations.')            
#
##if num_corrs is not None:
##    N = num_corrs
##elif num_wins is not None:
##    N = num_wins
##else:
##    N = len(self.select(station1,station2,tag).__correlations)
#    
## implement the above...
#
#
#stack,sdv,cnt,ncnt=self.lin_stack()
#
#if filename == None:
#    id = "%(network)s.%(station)s.%(location)s.%(channel)s"
#    filename = outdir+'/'+ id % self[0].stats_a + '.' + id % self[0].stats_b +\
#               '.' + self[0].correlation_type + '.' +\
#               tag+'.SAC'
#
#if stack is not None:
#    t = Trace(data=stack)
#    t.stats=self[0].stats_a.copy()
#    t.stats.sac['user0'] = ncnt
#    t.write(filename,format = 'SAC')
#
#plot(self,station1=None,station2=None):
#
#fig = plt.figure(figsize=(14, 6))
#fig.hold()
#plt.subplot(211)
#
#for tr in self.select(station1,station2).__correlations:
#    maxlag = tr.max_lag
#    npts = len(tr.correlation)
#    lag = np.linspace(-maxlag,maxlag,npts)
#    plt.plot(lag,tr.correlation/tr.stats_a.sac['user0'],color='0.75')
#stack,sdv,cnt,wcnt = self.lin_stack(station1,station2)
#
#plt.plot(lag,stack/wcnt,'r',linewidth=1.5)
#plt.grid()
#plt.ylabel('Correlations')
#plt.xlim([-maxlag,maxlag])
#
#plt.subplot(212)
#plt.plot(lag,sdv,'k--')
#plt.grid()
#plt.xlim([-maxlag,maxlag])
#plt.xlabel('Lag (s)')
#plt.ylabel('Standard deviation of stack')
#plt.show()
#
#plot_sect(self,annotate=True,sdv=False,\
#          station1=None,station2=None,exagg=1.,save_as=None):
#          
#fig = plt.figure(figsize=(10, 10))
#fig.hold()
#plt.subplot(111)
#
#ids = []
#
#for corr in self.select(station1,station2).__correlations:
#    
#    net1 = corr.stats_a.network
#    sta1 = corr.stats_a.station
#    net2 = corr.stats_b.network
#    sta2 = corr.stats_b.station
#    
#    id = net1+'.'+sta1+'-'+\
#         net2+'.'+sta2
#    if id not in ids:
#        ids.append(id)
#        dist=float(corr.stats_a.sac['dist'])/10000
#        lag = corr.stats_a.sac['e']
#        
#        stack,sdv,cnt,wcnt = self.lin_stack(sta1,sta2)
#        stack/=wcnt
#        
#        wcount=' ('+str(int(corr.stats_a.sac['user0']))+' windows)'
#        plt.plot(np.linspace(-lag,lag,len(corr.correlation)),\
#            exagg*stack+dist,'k')
#        #plt.plot(np.linspace(-lag,lag,len(corr.correlation)),\
#        #    stack+sdv+dist,'--',color='0.5')
#        #plt.plot(np.linspace(-lag,lag,len(corr.correlation)),\
#        #    stack-sdv+dist,'--',color='0.5')
#        plt.annotate(id, xy=(lag, dist+1),xytext=(-lag-3000,dist+1),\
#                     fontsize=12,color='k')
#        plt.annotate(wcount, xy=(lag, dist+0.5),xytext=(lag,dist),\
#                     fontsize=12,color='k')
#plt.xlabel('Lag (s)')
#plt.ylabel('Station-station dist (1000 km)')
#plt.ylim([0,2000])
#
#if save_as is not None:
#    plt.savefig(save_as,format='eps')
#plt.show()
#
#        
#
#lin_stack(self,station1=None,station2=None,location1=None,\
#          location2=None,tag=None,n_corr=None,\
#          n_win=None,npts=None):
## Implement a 'selection step' -- 
##- stations
##- start-, endtime
##- exclude stacking things that are based on different 
##- processing, i. e., the 'correlation stats' or at least the ids
##- have to be the same
#stack = None
#sdv = None
#cnt=0
#wcnt=0
#
#for tr in self.select(station1,station2,location1,location2,tag):
#    if stack==None:
#        stack=tr.correlation.copy()
#        stats_a=tr.stats_a
#        stats_b=tr.stats_b
#        max_lag=tr.stats_a.sac['e']
#        corrtype=tr.correlation_type
#        cnt=1
#        wcnt=tr.stats_a.sac['user0']
#    else:
#        # Find a better way to handle sampling rates better
#        if npts is not None and len(tr.correlation)!=npts:
#            continue
#        
#        if n_corr is not None and cnt>n_corr: break
#        if n_win is not None and cnt>n_win: break
#        
#        stack+=tr.correlation       
#        cnt+=1
#        wcnt+=tr.stats_a.sac['user0']
#        
## Maybe make this a bit more concise.
#if stack != None and cnt>1:
#    n=0.
#    wn=0
#    sdv = np.zeros(len(stack))
#    for tr in self.select(station1,station2):
#        if len(tr.correlation)!=len(stack): continue
#        if n_corr is not None and n>n_corr: break
#        if n_win is not None and wn>n_win: break
#        
#        sdv+=np.power(tr.correlation-stack/cnt,2)
#        n+=1.
#        wn+=tr.stats_a.sac['user0']
#        
#    sdv = np.sqrt(1./(n-1.)*sdv)
#    
#else:
#    sdv = None
#    
#return stack, sdv, cnt, wcnt
## 'cleaner': return a new correlationTrace, but where to put? 
## Should it overwite the old stack?
#
# 
#select(self,station1=None,station2=None,location1=None,\
#       location2=None,tag=None,\
#       endbefore=None, startafter=None):
#            
#if endbefore != None or startafter != None:
#    msg = "Not implemented, sorry."
#    raise NotImplementedError(msg)
#       
#correlations = []
#for corr in self.__correlations:
#    # skip trace if any given criterion is not matched
#    
#    if station1 is not None and station2 is None:
#        if not fnmatch.fnmatch(corr.stats_a.station.upper(),
#                               station1.upper())\
#           and not fnmatch.fnmatch(corr.stats_b.station.upper(),
#                               station1.upper()):                    
#            continue
#    if station1 is not None and station2 is not None:
#         if not fnmatch.fnmatch(corr.stats_a.station.upper(),
#                                station1.upper()):
#             continue
#    if station2 is not None:
#        if not fnmatch.fnmatch(corr.stats_b.station.upper(),
#                               station2.upper()):
#            continue
#    
#    if location1 is not None:
#        if not corr.stats_a.location in location1:
#            continue
#    if location2 is not None:
#        if not corr.stats_b.location in location2:
#            continue
#    if tag is not None:
#        if not fnmatch.fnmatch(corr.stats_a.sac['kt2'].upper(),
#                               tag.upper()):
#            continue
#        print(corr.stats_a.sac['kt2'].upper())
#            
#    correlations.append(corr)
#return self.__class__(correlations=correlations)
#        
#            
#
