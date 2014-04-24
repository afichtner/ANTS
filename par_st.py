from __future__ import print_function
import time
import sys
import os
 
import obspy as obs
import TOOLS.read_xml as rxml
import antconfig as cfg
import numpy as np
import TOOLS.correlations as corr

from mpi4py import MPI
from glob import glob
from obspy.core import Stats
from obspy.signal.tf_misfit import cwt

if __name__=='__main__':
    import par_st as pst
    xmlin=str(sys.argv[1])
    pst.par_st(xmlin)

import matplotlib.pyplot as plt

def par_st(xmlinput):
    """
    Script to parallely stack cross-correlation functions
    
    input: 
    
    xmlinput, string: Path to xml input file
    
    output: 
    
    None
    
    """
    
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size=comm.Get_size()
    #t0=time.time()
    
    #==============================================================================================
    #- MASTER process:
    #- reads in xmlinput
    #- gets the list of correlation pairs
    #- broadcasts both to workers
    #==============================================================================================
    if rank==0:
    
        print('The size is '+str(size),file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
        #- Read the input from xml file------------------------------------------------------------
        inp=rxml.read_xml(xmlinput)
        inp=inp[1]
        
        print('Read Input\n')
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
        #- copy the input xml to the output directory for documentation ---------------------------
        if os.path.exists(cfg.datadir+'/correlations/xmlinput/'+inp['corrname']+'.xml')==True:
            print('\n\nA new name was chosen to avoid overwriting an older run.\n\n',file=None)
            corrname=inp['corrname']+'_1'
            print('New name: '+corrname,file=None)
        else:
            corrname=inp['corrname']
            
        os.system('cp '+xmlinput+' '+cfg.datadir+'/correlations/xmlinput/'+corrname+'.xml')
        print('Copied input file',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
        #- Get list of correlation pairs-----------------------------------------------------------
        idpairs=parlistpairs(inp['data']['idfile'],int(inp['data']['npairs']),inp['data']['mix_channels'])
        print('Obtained list with correlations',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
    else:
        inp=dict()
        idpairs=list()
        corrname=''
    #==============================================================================================
    #- ALL processes:
    #- receive input and list of possible station pairs
    #- loop through 'block list'
    #==============================================================================================
    
    #- Broadcast the input and the list of pairs---------------------------------------------------
    idpairs=comm.bcast(idpairs, root=0)
    inp=comm.bcast(inp, root=0)
    corrname=comm.bcast(corrname,root=0)
    
    if rank==0:
        print('Broadcasted input and correlation pair list',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
    
    #- Create your own directory
    dir=cfg.datadir+'correlations/'+inp['timethings']['startdate'][0:4]+'/rank'+str(rank)+'/'
    if os.path.exists(dir)==False:
        os.mkdir(dir)
   
    if rank==0:
        print('Created own directory for output',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
       
    #- Each rank determines the part of data it has to work on ------------------------------------
    n1=int(len(idpairs)/size)
    n2=len(idpairs)%size
    ids=list()
    
    for i in range(0,n1):
        ids.append(idpairs[i*size+rank])
    if rank<n2:
        ids.append(idpairs[n1*size+rank])
    
    #- Print info to outfile of this rank --------------------------------------------------------
    if bool(int(inp['verbose']))==True:
        ofid=open(cfg.datadir+'/correlations/out/'+corrname+'.rank'+str(rank)+'.txt','w')
        print('\nRank number %d is correlating: \n' %rank,file=ofid)
        for block in ids:
            for tup in block:
                print(str(tup)+'\n',file=ofid)
    else:
        ofid=None
    
    if rank==0:
        print('Got task assigned, ready to start correlating',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
    #- Run correlation for blocks ----------------------------------------------------------------
    for block in ids:
        if rank==0:
            print('Starting a block of correlations',file=None)
            print(time.strftime('%H.%M.%S')+'\n',file=None)
        corrblock(inp,block,dir,corrname,ofid,bool(int(inp['verbose'])))
        if rank==0:
            print('Finished a block of correlations',file=None)
            print(time.strftime('%H.%M.%S')+'\n',file=None)
    
        
def corrblock(inp,block,dir,corrname,ofid=None,verbose=False):
    """
    Receives a block with station pairs
    Loops through those station pairs
    Checks if the required station data are already in memory by checking the list of ids in memory;
    if not, read the data and enter them into the list of ids in memory
    then assigns the two traces in question to dat1, dat2 which then are passed on to stacking routine
    it gets back the correlation stack and writes this stack to sac file. Metadata is written to sac header.
    
    input:
    
    inp, python dict object: a dictionary containing the input read in from the xmlinput file
    block, python list object: a list of tuples where every tuple is two station ids of 
    stations that should be correlated
    dir: Directory ro write to; needed so that every rank can write to its own directory
    ofid: output file id
    verbose: talk or shut up
    
    ouput:
    
    None
    
    """
   
    idlist=dict()
    datstr=obs.Stream()
    ntr=0

    
    for pair in block:
        id1=pair[0]
        id2=pair[1]
        endday=obs.UTCDateTime(inp['timethings']['enddate'])
        #==============================================================================================
        #- check if data for first station is in memory
        #- if it isn't, it needs to be read in
        #- typically if should be filtered
        #==============================================================================================
        
        if id1 in idlist:
            # Split the trace at its gaps
            str1=datstr[idlist[id1]].split()
        else:
            traces=glob(inp['selection']['indir']+'/'+id1+'.*.'+inp['selection']['prepname']+'.*')
            traces.sort()
            
            if len(traces)==0: continue
            
            #- collect a trace with masked samples where gaps are.
            #- This is convenient because we then get only one trace that can be handled with an index 
            #- inside the datstr objects more easily (rather than having a stream with variable number of traces)
            for filename in traces:
     
                (ey,em)=filename.split('/')[-1].split('.')[4:6]
                ef=obs.UTCDateTime(ey+','+em)
                if ef>endday:
                    continue

                newtr=obs.read(filename)[0]
                
                #- Check if at least one window contained
                if len(newtr.data)<int(inp['timethings']['winlen']):
                    continue
                
                if bool(int(inp['bandpass']['doit']))==True:
                    freqmin=float(inp['bandpass']['f_min'])
                    freqmax=float(inp['bandpass']['f_max'])
                    corners=int(inp['bandpass']['corners'])
                    newtr.filter('bandpass',freqmin=freqmin,freqmax=freqmax,corners=corners,zerophase=True)
                    
                if 'colltr' in locals():
                    colltr+=newtr
                else:
                    colltr=newtr.copy()    
            #- add this entire trace (which contains all data of this station that are available in this directory) to datstr and update the idlist
            if 'colltr' in locals():
                datstr+=colltr
                idlist.update({id1:ntr})
                ntr+=1
                if verbose:
                    print('Read in Traces for station '+id1,file=ofid)
                del newtr
                del colltr
                #- Split the trace at its gaps
                str1=datstr[idlist[id1]].split()
        
        #==============================================================================================
        #- Same thing for the second station, unless it's identical to the first
        #- check if data is in memory
        #- if it isn't, it needs to be read in
        #- typically if should be filtered
        #==============================================================================================
        if id2==id1:
            str2=str1
        else:
            if id2 in idlist:
                str2=datstr[idlist[id2]].split()
            else:
                traces=glob(inp['selection']['indir']+'/'+id2+'.*.'+inp['selection']['prepname']+'.*')
                traces.sort()
                
                if len(traces)==0: continue
                
                for filename in traces:
                    
                    (ey,em)=filename.split('/')[-1].split('.')[4:6]
                    ef=obs.UTCDateTime(ey+','+em)
                    if ef>endday:
                        continue
                        
                    newtr=obs.read(filename)[0]
                
                    #- Check if at least one window contained
                    if len(newtr.data)<int(inp['timethings']['winlen']):
                        continue
                
                    if bool(int(inp['bandpass']['doit']))==True:
                        freqmin=float(inp['bandpass']['f_min'])
                        freqmax=float(inp['bandpass']['f_max'])
                        corners=int(inp['bandpass']['corners'])
                        newtr.filter('bandpass',freqmin=freqmin,freqmax=freqmax,corners=corners,zerophase=True)
                    
                    if 'colltr' in locals():
                        colltr+=newtr
                    else:
                        colltr=newtr.copy()    
                #- add this entire trace (which contains all data of this station that are available in this directory) to datstr and update the idlist
                if 'colltr' in locals():
                    datstr+=colltr
                    idlist.update({id2:ntr})
                    ntr+=1
                    if verbose:
                        print('Read in Traces for station '+id2,file=ofid)
                    del newtr
                    del colltr
                    #- Split the trace at its gaps
                    str2=datstr[idlist[id2]].split()
        if len(datstr)==0:
		if verbose:
			print('No matching files found',file=ofid)
		return() 
        #==============================================================================================
        #- Get some information needed for the cross correlation
        #==============================================================================================
        startday=obs.UTCDateTime(inp['timethings']['startdate'])
        winlen=int(inp['timethings']['winlen'])
        maxlag=int(inp['correlations']['max_lag'])
        pccnu=int(inp['correlations']['pcc_nu'])
        verbose=bool(int(inp['verbose']))
        check=bool(int(inp['check']))
        Fs=float(inp['timethings']['Fs'])
        
        if ('freqmax' in locals())==False:
            freqmax=0.5*Fs
        if ('freqmin' in locals())==False:
            freqmax=0.001
        
        
        #==============================================================================================
        #- Run cross correlation
        #==============================================================================================
        
        (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1,str2,winlen,maxlag,pccnu,startday,endday,Fs,freqmin,freqmax,corrname,check,verbose)
        
        if npcc==0 and nccc==0:
            if verbose: print('No windows stacked for stations '+id1+' and '+id2,file=ofid)
            continue
        
        if verbose: print('Correlated traces from stations '+id1+' and '+id2,file=ofid)
        
        #==============================================================================================
        #- Write metadata info to sac header
        #- Store results
        #==============================================================================================
        (lat1, lon1, lat2, lon2, dist, az, baz)=rxml.get_coord_dist(id1.split('.')[0],id1.split('.')[1],id2.split('.')[0],id2.split('.')[1])
        tr_ccc=obs.Trace()
        tr_pcc=obs.Trace()
        
        for tr in (tr_ccc,tr_pcc):
            
            
            tr.stats.sampling_rate=Fs
            tr.stats.starttime=obs.UTCDateTime(2000, 01, 01)-maxlag
            tr.stats.network=id1.split('.')[0]
            tr.stats.station=id1.split('.')[1]
            tr.stats.location=id1.split('.')[2]
            tr.stats.channel=id1.split('.')[3]
            
            tr.stats.sac={}
            
            tr.stats.sac['user1']=winlen 
            tr.stats.sac['b']=-maxlag
            tr.stats.sac['e']=maxlag
            tr.stats.sac['kt0']=startday.strftime('%Y%j')
            tr.stats.sac['kt1']=endday.strftime('%Y%j')
            tr.stats.sac['iftype']=5
            tr.stats.sac['stla']=lat1
            tr.stats.sac['stlo']=lon1
            tr.stats.sac['kevnm']=id2.split('.')[1]
            tr.stats.sac['evla']=lat2
            tr.stats.sac['evlo']=lon2
            tr.stats.sac['dist']=dist
            tr.stats.sac['az']=az
            tr.stats.sac['baz']=baz
            tr.stats.sac['kuser0']=id2.split('.')[0]
            tr.stats.sac['kuser1']=id2.split('.')[2]
            tr.stats.sac['kuser2']=id2.split('.')[3]
        
        tr_ccc.data=ccc
        tr_pcc.data=pcc  
        tr_ccc.stats.sac['kcmpnm']='CCC'
        tr_pcc.stats.sac['kcmpnm']='PCC'
        tr_ccc.stats.sac['user0']=nccc
        tr_pcc.stats.sac['user0']=npcc
        
        
        #- open file and write correlation function
        fileid_ccc=dir+id1+'.'+id2+'.ccc.'+corrname+'.SAC'
        tr_ccc.write(fileid_ccc,format='SAC')      
        
        #- open file and write correlation function
        fileid_pcc=dir+id1+'.'+id2+'.pcc.'+corrname+'.SAC'
        tr_pcc.write(fileid_pcc,format='SAC')
        
        #- write coherence: As numpy array datafile (can be conveniently loaded again)
        fileid_ccc_cwt=dir+id1+'.'+id2+'.ccs.'+corrname+'.npy'
        np.save(fileid_ccc_cwt,cstack_ccc)
        
        fileid_pcc_cwt=dir+id1+'.'+id2+'.pcs.'+corrname+'.npy'
        np.save(fileid_pcc_cwt,cstack_pcc)
        
    

def parlistpairs(infile,nf,mix_cha,auto=False):
    """
    Find the 'blocks' to be processed by a single node.
    
    input:
    
    infile: Ascii file, contains all station ids to be correlated
    nf: number of pairs that should be in one block (to be held in memory and processed by one node)
    
    output:
    
    idpairs, python list object: list of tuples where each tuple contains two station ids
    
    """
    fid=open(infile,'r')
    idlist=fid.read().split('\n')
    
    idpairs=list()
    idcore=list()
    pcount=0
    
    for i in range(0,len(idlist)):
        if idlist[i]=='':
            continue
        
        for j in range(0,i+1):
            
            if idlist[j]=='': continue
            if idlist[i]==idlist[j] and auto==False:
                continue
            
            if pcount<nf:
                if idlist[i]<=idlist[j]:
                    idcore.append((idlist[i].split()[0],idlist[j].split()[0]))
                else:
                    idcore.append((idlist[j].split()[0],idlist[i].split()[0]))
                pcount+=1
            else:
                idpairs.append(idcore)
                idcore=list()
                if idlist[i]<=idlist[j]:
                    idcore.append((idlist[i].split()[0],idlist[j].split()[0])) 
                else:
                    idcore.append((idlist[j].split()[0],idlist[i].split()[0]))
                pcount=1
    idpairs.append(idcore) 
     
    return idpairs
    
    
def corr_pairs(str1,str2,winlen,maxlag,nu,startday,endday,Fs,fmin,fmax,corrname,check,verbose):
    """
    Step through the traces in the relevant streams and correlate whatever overlaps enough.
    
    input:
    
    str1,str2, obspy stream objects: Streams holding traces for station 1 and station 2
    winlen, int: window length in seconds
    maxlag, int: maximum lag for correlation in seconds
    nu, int: pcc nu, exponent for phase cross correlation
    startday, UTCDateTime object: Time where stack should start (if data available)
    endday, UTCDateTime object: Maximum time until where stacking should be carried out
    Fs, float: Sampling rate of data in Hz
    fmin: Minimum frequency the signal was filtered to (or 0.001 Hz)
    fmax: Maximum frequency the signal was filtered to (or 0.5*sampling frequency)
    verbose, boolean: loud or quiet
    
    output:
    
    cccstack, numpy array: classical cross correlation in time domain
    pccstack, numpy array: phase cross correlation
    ccccnt, int: Number of windows stacked for cccstack
    pcccnt, int: Number of windows stacked for pccstack
    
    
    """
    
    pcccnt=0
    ccccnt=0
    n1=0
    n2=0
    t1=startday
    
    while n1<len(str1) and n2<len(str2):
        
        
        # Check if the end of one of the traces has been reached
        if str1[n1].stats.endtime-t1<winlen-1:
            n1+=1
            continue
        elif str2[n2].stats.endtime-t1<winlen-1:
            n2+=1
            continue
        
        # Check the starttime of the potentially new trace
        t1=max(t1,str1[n1].stats.starttime,str2[n2].stats.starttime)
        t2=t1+winlen
        
        # Check if the end of the desired stacking window is reached
        if t2>endday: break
        
        tr1=str1[n1].slice(starttime=t1,endtime=t2-1/Fs)
        tr2=str2[n2].slice(starttime=t1,endtime=t2-1/Fs)
        
        # Correlate the traces
        if len(tr1.data)==winlen and len(tr2.data)==winlen:
            
            pcc=corr.phase_xcorrelation(tr1, tr2, maxlag, nu, varwl=True)
            ccc=corr.xcorrelation_td(tr1, tr2, maxlag)
            
        # Obtain the phase coherence
            npts = (2*maxlag+1)
            dt = 1.0/Fs
            coh_pcc=cwt(pcc, dt, 8.0, fmin, fmax)
            tol=np.mean(np.abs(coh_pcc))/1000.0
            coh_pcc=coh_pcc/(np.abs(coh_pcc)+tol)
            coh_ccc=cwt(ccc, dt, 8.0, fmin, fmax)
            tol=np.mean(np.abs(coh_ccc))/1000.0
            coh_ccc=coh_ccc/(np.abs(coh_ccc)+tol)
            
            
    
        # Write intermediate results if check flag is set
        # Include option to write phase coherence?
            if check==True:
                id1=str1[n1].id.split('.')[0]+'.'+str1[n1].id.split('.')[1]
                id2=str2[n2].id.split('.')[0]+'.'+str2[n2].id.split('.')[1]
                if os.path.exists(cfg.datadir+'/correlations/interm/'+id1+'_'+id2+'/')==False:
                    os.mkdir(cfg.datadir+'/correlations/interm/'+id1+'_'+id2+'/')
                id_corr=cfg.datadir+'/correlations/interm/'+id1+'_'+id2+'/'+str1[n1].id.split('.')[1]+'.'+str2[n2].id.split('.')[1]+t1.strftime('.%Y.%j.%H.%M.%S.')+corrname+'.SAC'
                trace_pcc=obs.Trace(data=pcc)
                trace_pcc.stats=Stats({'network':'pcc','station':str1[n1].id.split('.')[1],'location':str2[n2].id.split('.')[1],'sampling_rate':Fs})
                trace_pcc.write(id_corr, format='SAC')
                
                id_cwt=cfg.datadir+'/correlations/interm/'+id1+'_'+id2+'/'+str1[n1].id.split('.')[1]+'.'+str2[n2].id.split('.')[1]+t1.strftime('.%Y.%j.%H.%M.%S.')+corrname+'.npy'
                np.save(id_cwt,coh_pcc)
                
            #Add to the stack
            if 'pccstack' in locals():
                pccstack+=pcc
            else:
                pccstack=pcc
            pcccnt+=1
            if 'cccstack' in locals():
                cccstack+=ccc
            else:
                cccstack=ccc
            ccccnt+=1
            if 'cstack_ccc' in locals():
                cstack_ccc+=coh_ccc
            else:
                cstack_ccc=coh_ccc
            if 'cstack_pcc' in locals():
                cstack_pcc+=coh_pcc
            else:
                cstack_pcc=coh_pcc
            
        #Update starttime
        t1=t2
    
    return(cccstack,pccstack,cstack_ccc,cstack_pcc,ccccnt,pcccnt)
            
            
                
        
       
    
    
