from __future__ import print_function
import time
import sys
import os
 
import obspy as obs
import TOOLS.read_xml as rxml
import antconfig as cfg
import numpy as np
import TOOLS.processing as proc
import TOOLS.rotationtool as rt

from mpi4py import MPI
from glob import glob
from obspy.core import Stats, Trace
from obspy.noise.correlation import Correlation
from obspy.noise.correlation_functions import phase_xcorr
from obspy.signal import cross_correlation
from obspy.signal.tf_misfit import cwt
from scipy.signal import hilbert

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
    
    
#==============================================================================
    #- MASTER process:
    #- reads in xmlinput
    #- gets the list of correlation pairs
    #- broadcasts both to workers    
#==============================================================================

    if rank==0:
    
        print('The size is '+str(size),file=None)
        
        #- Read the input from xml file----------------------------------------
        inp=rxml.read_xml(xmlinput)
        inp=inp[1]
        
        
        #- copy the input xml to the output directory for documentation -------
        if os.path.exists(cfg.datadir+'/correlations/xmlinput/'+\
            inp['corrname']+'.xml') == False or\
                 bool(int(inp['update'])) == True:
            corrname=inp['corrname']
        else:
            msg = 'Avoiding to overwrite data: choose a new name tag or change\
              update mode to True. Name tag %s already taken.' %inp['corrname']
            raise ValueError(msg)
        
        if os.path.exists(cfg.datadir+'correlations/'+corrname) == False:
            os.mkdir(cfg.datadir+'correlations/'+corrname)
        
           
        os.system('cp '+xmlinput+' '+cfg.datadir+'/correlations/xmlinput/'+\
            corrname+'.xml')
        print(corrname+'\n',file=None)
        print('Copied input file',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
        #- Get list of correlation pairs----------------------------------------
        idpairs=parlistpairs(inp['data']['idfile'],int(inp['data']['npairs']),\
            corrname,bool(int(inp['correlations']['autocorr'])))
        print('Obtained list with correlations',file=None)
        print('Number of possible correlations: '+str(len(idpairs*int(inp['data']['npairs']))))
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
    else:
        inp=dict()
        idpairs=list()
        corrname=''
    
#==============================================================================
    #- ALL processes:
    #- receive input and list of possible station pairs
    #- loop through 'block list'    
#==============================================================================
    
    #- Broadcast the input and the list of pairs-------------------------------
    idpairs=comm.bcast(idpairs, root=0)
    inp=comm.bcast(inp, root=0)
    corrname=comm.bcast(corrname,root=0)
    
    if rank==0:
        print('Broadcasted input and correlation pair list',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
    
    #- Create your own directory
    dir=cfg.datadir+'correlations/'+corrname+'/rank'+str(rank)+'/'
    if os.path.exists(dir)==False:
        os.mkdir(dir)
   
    if rank==0:
        print('Created output directory',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
       
    #- Each rank determines the part of data it has to work on ----------------
    n1=int(len(idpairs)/size)
    n2=len(idpairs)%size
    ids=list()
    
    for i in range(0,n1):
        ids.append(idpairs[i*size+rank])
    if rank<n2:
        ids.append(idpairs[n1*size+rank])
    
    #- Print info to outfile of this rank --------------------------------------
    if bool(int(inp['verbose']))==True:
        ofid=open(cfg.datadir+'/correlations/out/'+corrname+'.rank'+str(rank)+\
            '.txt','w')
        print('\nRank number %d is correlating: \n' %rank,file=ofid)
        for block in ids:
            for tup in block:
                print(str(tup),file=ofid)
    else:
        ofid=None
    
    if rank==0:
        print('Station pairs assigned, start correlating',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
    #- Run correlation for blocks ----------------------------------------------
    for block in ids:
        
        corrblock(inp,block,dir,corrname,ofid,bool(int(inp['verbose'])))
        if rank==0:
            print('Finished a block of correlations',file=None)
            print(time.strftime('%H.%M.%S'),file=None)
        
        # Flush the outfile buffer every now and then...
        if bool(int(inp['verbose']))==True:
	    ofid.flush()
    
    
    os.system('mv '+dir+'/* '+cfg.datadir+'correlations/'+corrname+'/')
    os.system('rmdir '+dir)
    print('Rank %g finished correlations.' %rank,file=None)
        
def corrblock(inp,block,dir,corrname,ofid=None,verbose=False):
    """
    Receives a block with station pairs
    Loops through those station pairs
    Checks if the required station data are already
    in memory by checking the list of ids in memory;
    if not, read the data and enter them into the list of ids in memory
    then assigns the two traces in question to dat1, dat2 
    These are passed on to stacking routine, which passes
    back the correlation stack and writes this stack to sac file. 
    Metadata is written to sac header.
    
    input:
    inp, python dict object: a dictionary containing the input read in from 
    the xmlinput file
    block, python list object: a list of tuples where every tuple is two 
    station ids of stations that should be correlated
    dir: Directory ro write to; needed so that every rank can write to its own 
    directory
    ofid: output file id
    verbose: talk or shut up
    
    ouput:
    None
    
    """
   
    datstr=obs.Stream()
    idlist=list()
    
    
#==============================================================================
    #- Get some information needed for the cross correlation   
#============================================================================== 
    
    startday=obs.UTCDateTime(inp['timethings']['startdate'])
    endday=obs.UTCDateTime(inp['timethings']['enddate'])
    Fs=inp['timethings']['Fs'].split(' ')
    Fs = [float(fs) for fs in Fs] 
    winlen=int(inp['timethings']['winlen'])
    maxlag=int(inp['correlations']['max_lag'])
    overlap=int(inp['timethings']['olap'])
    pccnu=int(inp['correlations']['pcc_nu'])
    verbose=bool(int(inp['verbose']))
    check=bool(int(inp['check']))
    cha=inp['data']['channel']
    comp=inp['data']['components']
    mix_cha=bool(int(inp['data']['mix_cha']))
    indir=inp['selection']['indir']
    corrtype=inp['correlations']['corrtype']
    prepname=inp['selection']['prepname']
    tfpws = bool(int(inp['correlations']['tfpws']))
    
    if inp['bandpass']['doit']=='1':
        freqmin=float(inp['bandpass']['f_min'])
        freqmax=float(inp['bandpass']['f_max'])
        prefilt=(freqmin,freqmax,float(inp['bandpass']['corners']))
    else:
        prefilt=None
        freqmax=0.5*min(Fs)
        freqmin=0.001

    for pair in block:
        str1=obs.Stream()
        str2=obs.Stream()
        id1 = pair[0]
        id2 = pair[1]
        
        if comp=='Z':
            
            id1 = [id1+cha+'Z']
            id2 = [id2+cha+'Z']
            
        elif comp=='EN':
            id1=[id1+cha+'E', id1+cha+'N', id1+cha+'1', id1+cha+'2']
            id2=[id2+cha+'E', id2+cha+'N', id2+cha+'1', id2+cha+'2']
            
        
#==============================================================================
        #- check if data for first station is in memory
        #- if it isn't, it needs to be read in
        #- typically if should be filtered 
#==============================================================================
        for id in id1:
            
            station = id.split('.')[1]
            channel = id.split('.')[-1]
            
            if id in idlist:
                str1 += datstr.select(station=station, channel=channel).split()
            else:
                (colltr,readsuccess) = \
                                addtr(id,indir,prepname,winlen,endday,prefilt)
        
                #- add this entire trace (which contains all data of this 
                #- station that are available in this directory) to datstr and 
                #- update the idlist
                if readsuccess == True:
                    datstr += colltr
                    str1 += colltr.split()
                    idlist += id
                    
                   # if verbose:
                   #     print('Read in traces for channel '+id,file=ofid)
                   # del colltr
                else:
                   # if verbose:
                   #     print('No traces found for channel '+id,file=ofid)
                    continue
            
        
#==============================================================================
        #- Same thing for the second station, unless it's identical to the 1st
        #- check if data is in memory
        #- if it isn't, it needs to be read in
        #- typically if should be filtered        
#==============================================================================
        if id2 == id1:
            str2 = str1
        else:
            for id in id2:
                station = id.split('.')[1]
                channel = id.split('.')[-1]
                
                if id in idlist:
                    str2 += datstr.select(station=station, \
                        channel=channel).split()
                else:
                    (colltr,readsuccess)=addtr(id,indir,prepname,\
                        winlen,endday,prefilt)
                    
                    if readsuccess == True:
                        datstr += colltr
                        str2 += colltr.split()
                        idlist += id
                        
                        
                        if verbose:
                            print('Read in traces for channel '+id,file=ofid)
                        del colltr
                    else:
                        if verbose:
                            print('No traces found for channel '+id,file=ofid)
                        continue
                    
        
#==============================================================================
        #- No files found?
        
#==============================================================================
           
        if len(str1) == 0 or len(str2) == 0:
            
            if verbose==True:
                print('No data found for one or both of:\n',file=ofid)
                print(str(id1)+str(id2),file=ofid)
                continue
            else:
                continue    
            
		       
            
        
        
#==============================================================================
        #- Rotate horizontal traces        
#==============================================================================

        #- Geoinf: (lat1, lon1, lat2, lon2, dist, az, baz)
        geoinf=rxml.get_coord_dist(id1[0].split('.')[0],id1[0].split('.')[1],\
            id2[0].split('.')[0],id2[0].split('.')[1])
        
        if comp=='EN':
            try:
                str1 = rt.rotate_streams(str1.select(component='N'),\
                    str1.select(component='E'),geoinf[6],'NE->RT')
                str2 = rt.rotate_streams(str2.select(component='N'),\
                    str2.select(component='E'),geoinf[6],'NE->RT')
            
            except ValueError:
                print('East and North traces do not cover same time span,not\
                    rotated',file=None)
                continue
            
            str1_T=str1.select(channel=cha+'T').split()
            str1_R=str1.select(channel=cha+'R').split()
            str2_T=str2.select(channel=cha+'T').split()
            str2_R=str2.select(channel=cha+'R').split()
            
            if len(str1_T) == 0 or len(str1_R) == 0 \
                or len(str2_T) == 0 or len(str2_R) == 0:
                        
                if verbose==True: 
                    print('No rotated traces: Probably original traces \
                    are components 1, 2. Rotation not implemented yet.',file=ofid)
                continue
        
        
        
#==============================================================================
        #- Run cross correlation
        
#==============================================================================
        
        #- Case: Mix channels True or false and channel==z: Nothing special 
        if comp=='Z':
            id_1=str1[0].id
            id_2=str2[0].id
            if verbose == True:
                print(id_1,file=ofid)
                print(id_2,file=ofid)
            
            (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1,str2,\
                winlen,overlap,maxlag,pccnu,tfpws,startday,endday,Fs,freqmin,\
                freqmax,corrname,corrtype,check,verbose)
            
            if npcc!=0 or nccc!=0:
                savecorrsac(ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc,id_1,\
                    id_2,geoinf,winlen,overlap,maxlag,pccnu,tfpws,startday,\
                    endday,Fs,freqmin,freqmax,corrname,corrtype,dir,check,\
                    verbose)
                if verbose: print('Correlated traces from stations '+id1[0]+\
                ' and '+id2[0],file=ofid)
            
        elif comp=='EN':
            id1_T=str1_T[0].id
            id1_R=str1_R[0].id
            id2_T=str2_T[0].id
            id2_R=str2_R[0].id
            
            
            (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_T,\
                str2_T,winlen,overlap,maxlag,pccnu,tfpws,startday,endday,Fs,\
                freqmin,freqmax,corrname,corrtype,check,verbose)
                
            if npcc != 0 or nccc != 0:
                savecorrsac(ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc,id1_T,\
                    id2_T,geoinf,winlen,overlap,maxlag,pccnu,tfpws,startday,\
                    endday,Fs,freqmin,freqmax,corrname,corrtype,dir,check,\
                    verbose)
            del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc
            
            # Component RR
            (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_R,\
                str2_R,winlen,overlap,maxlag,pccnu,tfpws,startday,endday,Fs,\
                freqmin,freqmax,corrname,corrtype,check,verbose)
            if npcc != 0 or nccc != 0:
                savecorrsac(ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc,id1_R,\
                    id2_R,geoinf,winlen,overlap,maxlag,pccnu,tfpws,startday,\
                    endday,Fs,freqmin,freqmax,corrname,corrtype,dir,check,\
                    verbose)
            del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc
            
            if mix_cha == True:
            # Get the remaining component combinations
            # Component T1R2
                (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_T,\
                    str2_R,winlen,overlap,maxlag,pccnu,tfpws,startday,endday,Fs,\
                    freqmin,freqmax,corrname,corrtype,check,verbose)
                if npcc != 0 or nccc != 0:
                    savecorrsac(ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc,id1_T,\
                    id2_R,geoinf,winlen,overlap,maxlag,pccnu,tfpws,startday,\
                    endday,Fs,freqmin,freqmax,corrname,corrtype,dir,check,\
                    verbose)
                del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc
            # Component R1T2
                (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_R,\
                    str2_T,winlen,overlap,maxlag,pccnu,tfpws,startday,endday,Fs,\
                    freqmin,freqmax,corrname,corrtype,check,verbose)
                if npcc != 0 or nccc != 0:
                    savecorrsac(ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc,id1_R,\
                    id2_T,geoinf,winlen,overlap,maxlag,pccnu,tfpws,startday,\
                    endday,Fs,freqmin,freqmax,corrname,corrtype,dir,check,\
                    verbose)
                del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc


def parlistpairs(infile,nf,corrname,auto=False):
    """
    Find the 'blocks' to be processed by a single node.
    
    input:
    infile: Ascii file, contains all station ids to be correlated
    nf: number of pairs that should be in one block (to be held in memory and 
    processed by one node)
    auto: whether or not to calculate autocorrelation
    
    output:
    idpairs, python list object: list of tuples where each tuple contains two 
    station ids
    
    """
    fid=open(infile,'r')
    ids=fid.read().split('\n')
    idlist=list()
    
    for item in ids:
        #- Sort out empty lines
        if item=='': continue
        #- Sort out doubles
        if item not in idlist:
            idlist.append(item)
    
    idpairs=list()
    idcore=list()
    pcount=0
    
    for i in range(0,len(idlist)):
       
        for j in range(0,i+1):
            
            #- In update mode: Check if the correlation is there already
            if idlist[i]<=idlist[j]:
                fileid = cfg.datadir + 'correlations/' + corrname + '/' +\
                idlist[i] + '???.' + idlist[j] + '???.?cc.' + corrname + '.SAC'
            else:
                fileid = cfg.datadir + 'correlations/' + corrname + '/' + \
                idlist[j] + '???.' + idlist[i] + '???.?cc.' + corrname + '.SAC'
                
            if glob(fileid) != []:
                continue
            
            #- Autocorrelation?
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
    
    
def corr_pairs(str1,str2,winlen,overlap,maxlag,nu,tfpws,startday,endday,Fs_new,\
                    fmin,fmax,corrname,corrtype,check,verbose):
    """
    Step through the traces in the relevant streams and correlate whatever 
    overlaps enough.
    
    input:
    
    str1,str2, obspy stream objects: Streams holding traces for station 1 and 
    station 2
    winlen, int: window length in seconds
    overlap, int: overlap in seconds
    maxlag, int: maximum lag for correlation in seconds
    nu, int: pcc nu, exponent for phase cross correlation
    tfpws, boolean: type of phase weighted stack (if true time-frequency pws 
    is 
    calculated, otherwise time domain)
    startday, UTCDateTime object: Time where stack should start (if data 
    available)
    endday, UTCDateTime object: Maximum time until where stacking should be 
    carried out
    Fs, float: Sampling rate of data in Hz
    fmin: Minimum frequency the signal was filtered to (or 0.001 Hz)
    fmax: Maximum frequency the signal was filtered to (or 0.5*sampling 
    frequency)
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
    tlen=int(maxlag*Fs_new[-1])*2+1
    cccstack=np.zeros(tlen)
    pccstack=np.zeros(tlen)
    cstack_ccc=np.zeros(tlen,dtype=np.complex)
    cstack_pcc=np.zeros(tlen,dtype=np.complex)
    
    
    while n1<len(str1) and n2<len(str2):
    
        # Check if the end of one of the traces has been reached
        if str1[n1].stats.endtime-t1<winlen-1:
            n1+=1
            #print('No more windows in trace 1..',file=None)
            continue
        elif str2[n2].stats.endtime-t1<winlen-1:
            n2+=1
            #print('No more windows in trace 2..',file=None)
            continue
        
        # Check the starttime of the potentially new trace
        t1=max(t1,str1[n1].stats.starttime,str2[n2].stats.starttime)
        #print(t1,file=None)
        t2=t1+winlen
        
        # Check if the end of the desired stacking window is reached
        if t2>endday: 
            #print('At end of correlation time',file=None)
            break
        
        
        tr1=str1[n1].slice(starttime=t1,endtime=t2-1/Fs_new[-1])
        tr2=str2[n2].slice(starttime=t1,endtime=t2-1/Fs_new[-1])
        
        #- Downsample
        if len(tr1.data)>40 and len(tr2.data)>40:
            k=0
            while k<len(Fs_new):
                if Fs_new[k]<tr1.stats.sampling_rate:
                    tr1=proc.trim_next_sec(tr1,False,None)
                    tr1=proc.downsample(tr1,Fs_new[k],False,None)
                if Fs_new[k]<tr2.stats.sampling_rate:
                    tr2=proc.trim_next_sec(tr2,False,None)
                    tr2=proc.downsample(tr2,Fs_new[k],False,None)
                k+=1
        else:
            t1 = t2 - overlap
            continue   
    
        
        if len(tr1.data) == len(tr2.data):
            mlag = maxlag / tr1.stats.delta
            mlag=int(mlag)
            
            # Check if the traces are both long enough
            if len(tr1.data)<=2*mlag or len(tr2.data)<=2*mlag:
                t1 = t2 - overlap
                continue
            if tr1.data.any()==np.nan or tr2.data.any()==np.nan:
                t1 = t2 - overlap
                continue
            if tr1.data.any()==np.inf or tr2.data.any()==np.inf:
                t1 = t2 - overlap
                continue
                
            #- Classical correlation part =====================================
            if corrtype == 'ccc' or corrtype == 'both':
                ccc=classic_xcorr(tr1, tr2, mlag)
                coh_ccc = hilbert(ccc)
                tol = np.max(coh_ccc)/10000
                coh_ccc = coh_ccc/(np.absolute(coh_ccc)+tol)
                cccstack+=ccc
                ccccnt+=1
                cstack_ccc+=coh_ccc
                    
                if check==True:
                    id1=str1[n1].id.split('.')[0]+'.'+str1[n1].id.split('.')[1]
                    id2=str2[n2].id.split('.')[0]+'.'+str2[n2].id.split('.')[1]
                    
                    if os.path.exists(cfg.datadir+'/correlations/interm/'+id1+\
                        '_'+id2+'/')==False:
                        os.mkdir(cfg.datadir+'/correlations/interm/'+\
                            id1+'_'+id2+'/')
                            
                    id_corr=cfg.datadir+'/correlations/interm/'+id1+'_'+id2+\
                        '/'+str1[n1].id.split('.')[1]+'.'+\
                        str2[n2].id.split('.')[1]+\
                        t1.strftime('.%Y.%j.%H.%M.%S.')+\
                        corrname+'.ccc.SAC'
                        
                    trace_ccc=obs.Trace(data=pcc)
                    trace_ccc.stats=Stats({'network':'pcc',\
                        'station':str1[n1].id.split('.')[1],\
                        'location':str2[n2].id.split('.')[1],\
                        'sampling_rate':Fs})
                        
                    trace_ccc.write(id_corr, format='SAC')
                    
                    id_coh=cfg.datadir+'/correlations/interm/'+id1+'_'+id2+\
                        '/'+str1[n1].id.split('.')[1]+'.'+\
                        str2[n2].id.split('.')[1]+\
                        t1.strftime('.%Y.%j.%H.%M.%S.')+corrname+'.ccc.npy'
                        
                    np.save(id_coh,coh_pcc)
                    
            #- Phase correlation part =========================================
            if corrtype == 'pcc' or corrtype == 'both':
                pcc=phase_xcorr(tr1, tr2, mlag, nu)
                coh_pcc = hilbert(pcc)
                tol = np.max(coh_pcc)/10000
                coh_pcc = coh_pcc/(np.absolute(coh_pcc)+tol)
                pccstack+=pcc
                pcccnt+=1
                cstack_pcc+=coh_pcc
                    
                if check==True:
                    id1=str1[n1].id.split('.')[0]+'.'+str1[n1].id.split('.')[1]
                    id2=str2[n2].id.split('.')[0]+'.'+str2[n2].id.split('.')[1]
                    
                    if os.path.exists(cfg.datadir+'/correlations/interm/'+id1+\
                        '_'+id2+'/')==False:
                        os.mkdir(cfg.datadir+'/correlations/interm/'+\
                            id1+'_'+id2+'/')
                            
                    id_corr=cfg.datadir+'/correlations/interm/'+id1+'_'+id2+\
                        '/'+str1[n1].id.split('.')[1]+'.'+\
                        str2[n2].id.split('.')[1]+\
                        t1.strftime('.%Y.%j.%H.%M.%S.')+\
                        corrname+'.pcc.SAC'
                        
                    trace_pcc=obs.Trace(data=pcc)
                    trace_pcc.stats=Stats({'network':'pcc',\
                        'station':str1[n1].id.split('.')[1],\
                        'location':str2[n2].id.split('.')[1],\
                        'sampling_rate':Fs})
                        
                    trace_pcc.write(id_corr, format='SAC')
                    
                    id_coh=cfg.datadir+'/correlations/interm/'+id1+'_'+id2+\
                        '/'+str1[n1].id.split('.')[1]+'.'+\
                        str2[n2].id.split('.')[1]+\
                        t1.strftime('.%Y.%j.%H.%M.%S.')+corrname+'.pcc.npy'
                        
                    np.save(id_coh,coh_pcc)
                    
                    
        #Update starttime
            t1 = t2 - overlap
        else:
            #print('Traces have unequal length!',file=None)
            t1 = t2 - overlap

    return(cccstack,pccstack,cstack_ccc,cstack_pcc,ccccnt,pcccnt)
    
    
    
    
    
    

def addtr(id,indir,prepname,winlen,endday,prefilt):
    
    """
    Little reader.
    
    Needs to read all available data of one channel with a specified tag in a 
    directory.
    
    """
 
    traces=glob(indir+'/'+id+'.*.'+prepname+'.*')
    traces.sort()
    readit=False
    
    
    if len(traces) == 0:
        return (Trace(),False)
             
    #- collect a trace with masked samples where gaps are.
    #- This is convenient because we then get only one trace that can be 
    #- handled with an index 
    #- inside the datstr objects more easily (rather than having a stream with 
    #- variable number of traces)
    for filename in traces: 
        
        (ey,em)=filename.split('/')[-1].split('.')[4:6]
        ef=obs.UTCDateTime(ey+','+em)
        
        if ef>endday:
            continue
        
        try:
            newtr=obs.read(filename)
        except:
            print('Problems opening data file:\n',file=ofid)
            print(tr,file=ofid)
            continue
        
        for tr in newtr:
            #- Check if at least one window contained
            if len(tr.data)*tr.stats.delta<winlen-tr.stats.delta:
                continue
            
            #- Bandpass filter
            if prefilt is not None:
                tr.filter('bandpass',freqmin=prefilt[0],freqmax=prefilt[1],\
                corners=prefilt[2],zerophase=True)
            
            if readit==False:
                colltr=tr.copy()
                #print(colltr,file=None)
                readit=True
                #print('Started coll. trace',file=None)
                #print(colltr,file=None)
            else:
                try:
                    colltr+=tr
                except TypeError:
                    continue   
     
    if readit == True:           
        return (colltr,readit)
    else:
        return(Trace(),readit)
    
   
def savecorrsac(ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc,id1,id2,geoinf,winlen,\
    overlap,maxlag,pccnu,tfpws,startday,endday,Fs,freqmin,freqmax,corrname,\
    corrtype,outdir,check,verbose):
    
    
    
#==============================================================================
    #- Write metadata info to sac header
    #- Store results
#==============================================================================

    traces = list()
    
    if corrtype =='pcc' or corrtype == 'both':
        tr_pcc=obs.Trace(data=pcc)
        tr_pcc.stats.sac={}
        tr_pcc.stats.sac['user0']=npcc
        traces.append(tr_pcc)
        
        
    elif corrtype =='ccc' or corrtype == 'both':
        tr_ccc=obs.Trace(data=ccc)
        tr_ccc.stats.sac={}
        tr_ccc.stats.sac['user0']=nccc
        traces.append(tr_ccc)
        
    (lat1, lon1, lat2, lon2, dist, az, baz)=geoinf
    if type(Fs) == list: Fs=Fs[-1]
    
    for tr in traces:

        tr.stats.sampling_rate=Fs
        tr.stats.starttime=obs.UTCDateTime(2000, 01, 01)-maxlag*Fs
        tr.stats.network=id1.split('.')[0]
        tr.stats.station=id1.split('.')[1]
        tr.stats.location=id1.split('.')[2]
        tr.stats.channel=id1.split('.')[3]
        
        tr.stats.sac['user1']=winlen
        tr.stats.sac['user2']=overlap
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

    
    
    if corrtype == 'pcc' or corrtype =='both':
        #- open file and write correlation function
        fileid_pcc=outdir+id1+'.'+id2+'.pcc.'+corrname+'.SAC'
        tr_pcc.write(fileid_pcc,format='SAC')
        
        fileid_pcc_cwt=outdir+id1+'.'+id2+'.pcs.'+corrname+'.npy'
        np.save(fileid_pcc_cwt,cstack_pcc)
        
    if corrtype == 'ccc' or corrtype =='both':
        #- open file and write correlation function
        fileid_ccc=outdir+id1+'.'+id2+'.ccc.'+corrname+'.SAC'
        tr_ccc.write(fileid_ccc,format='SAC')      
        
        #- write coherence: As numpy array datafile (can be conveniently 
        #- loaded again)
        fileid_ccc_cwt=outdir+id1+'.'+id2+'.ccs.'+corrname+'.npy'
        np.save(fileid_ccc_cwt,cstack_ccc)
        
    
    
def classic_xcorr(data1, data2, max_lag):
   
    xcorr = cross_correlation.xcorr(data1.data * 1000, data2.data * 1000,\
        max_lag, True)[2]
    
    return xcorr
