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
import INPUT.input_correlation as inp

from math import sqrt
from glob import glob
from obspy.core import Stats, Trace
from obspy.noise.correlation import Correlation
from obspy.noise.correlation_functions import phase_xcorr
from obspy.signal import xcorr
from obspy.signal.util import nextpow2
from obspy.signal.tf_misfit import cwt
from scipy.signal import hilbert

if __name__=='__main__':
    import par_st_daint as pst
    rank = int(os.environ['ALPS_APP_PE'])
    size = int(sys.argv[1])
    
    if inp.update == False and os.path.exists(cfg.datadir+\
    '/correlations/input/'+inp.corrname+'.txt') == True:
        sys.exit('Choose a new correlation name tag or set update=True.\
         Aborting.')
    pst.par_st(size, rank)
    

def par_st(size,rank):
    """
    Script to parallely stack cross-correlation functions
    
    input: 
    
    number of tasks, rank of this task
    
    output: 
    
    None
    
    """

    
    
#==============================================================================
    #- MASTER process:
   
    #- gets the list of correlation pairs
    #- broadcasts both to workers    
#==============================================================================
    corrname = inp.corrname
    try:
        os.mkdir(cfg.datadir+'correlations/'+corrname)
    except OSError:
        pass
    
    if rank==0:
        print('The size is '+str(size),file=None)
        os.system('cp INPUT/input_correlation.py '+cfg.datadir+\
        '/correlations/input/'+corrname+'.txt')
        print(corrname+'\n',file=None)
        print('Copied input file',file=None)
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
    #- Get list of correlation pairs----------------------------------------
    idpairs=parlistpairs(corrname)
    
    if rank == 0:
        print('Obtained list with correlations',file=None)
        print('Approx. number of possible correlations: '+str(len(idpairs)*inp.npairs))
        print(time.strftime('%H.%M.%S')+'\n',file=None)
        
#==============================================================================
    #- ALL processes:
    #- receive input and list of possible station pairs
    #- loop through 'block list'    
#==============================================================================

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
    if inp.verbose==True:
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
        
        corrblock(block,dir,corrname,ofid)
        if rank==0:
            print('Finished a block of correlations',file=None)
            print(time.strftime('%H.%M.%S'),file=None)
        
        # Flush the outfile buffer every now and then...
        if inp.verbose==True:
	    ofid.flush()
    
    
    os.system('mv '+dir+'/* '+cfg.datadir+'correlations/'+corrname+'/')
    os.system('rmdir '+dir)
    print('Rank %g finished correlations.' %rank,file=None)
        
def corrblock(block,dir,corrname,ofid=None):
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
    verbose=inp.verbose
    
#==============================================================================
    #- Get some information needed for the cross correlation   
#============================================================================== 
    
    
    cha=inp.channel
    comp=inp.components
    mix_cha=inp.mix_cha
    

    for pair in block:
        str1=obs.Stream()
        str2=obs.Stream()
        id1 = pair[0]
        id2 = pair[1]
        
        if comp=='Z':
            
            id1 = [id1+cha+'Z']
            id2 = [id2+cha+'Z']
            
        elif comp=='RT' or comp=='R' or comp=='T':
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
                (colltr,readsuccess) = addtr(id)
        
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
                    (colltr,readsuccess)=addtr(id)
                    
                    if readsuccess == True:
                        datstr += colltr
                        str2 += colltr.split()
                        idlist += id
                        
                        
                        if inp.verbose:
                            print('Read in traces for channel '+id,file=ofid)
                        del colltr
                    else:
                        if inp.verbose:
                            print('No traces found for channel '+id,file=ofid)
                        continue
                    
        
#==============================================================================
        #- No files found?
        
#==============================================================================
           
        if len(str1) == 0 or len(str2) == 0:
            
            if inp.verbose==True:
                print('No data found for one or both of:\n',file=ofid)
                print(str(id1)+str(id2),file=ofid)
                continue
            else:
                continue    
            
	    
#==============================================================================
        #- Rotate horizontal traces        
#==============================================================================
        #- Get information on the geography of the two traces
        #- This is all not very beautiful, could be done up sometime
        try: 
            lat1 = str1[0].stats['lat']
            lon1 = str1[0].stats['lon']
            lat2 = str2[0].stats['lat']
            lon2 = str2[0].stats['lon']
        except KeyError:
            (lat1,lon1,lat2,lon2) = rxml.get_coord_staxml(id1[0].\
            split('.')[0],id1[0].split('.')[1],\
            id2[0].split('.')[0],id2[0].split('.')[1])
            if (lat1,lon1,lat2,lon2) ==(0,0,0,0):
                print('Problems with metadata: No station coordinates \
found, setting distance to zero for station pair:')
                print(id1[0].\
            split('.')[0],id1[0].split('.')[1],\
            id2[0].split('.')[0],id2[0].split('.')[1]+'\n')
            
        #- Geoinf: (lat1, lon1, lat2, lon2, dist, az, baz)
        geoinf=rxml.get_geoinf(lat1,lon1,lat2,lon2)
        
        if comp=='RT' or comp=='R' or comp=='T':
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
                        
                if inp.verbose==True: 
                    print('No rotated traces: Original traces \
                    are components 1, 2. Rotation not implemented yet.',file=ofid)
                continue
        
        
        
#==============================================================================
        #- Run cross correlation
        
#==============================================================================
        
        # Z components ========================================================
        #- Case: Mix channels True or false and channel==z: Nothing special
        #======================================================================
        if comp=='Z':
            id_1=str1[0].id
            id_2=str2[0].id
            if inp.verbose == True:
                print(id_1,file=ofid)
                print(id_2,file=ofid)
            
            (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1,str2,\
                corrname,geoinf)
            
            if npcc != 0:
                savecorrs(pcc,cstack_pcc,npcc,id_1,\
                    id_2,geoinf,corrname,'pcc',dir)
            if nccc != 0:
                savecorrs(ccc,cstack_ccc,nccc,id_1,\
                    id_2,geoinf,corrname,'ccc',dir)
                if inp.verbose: print('Correlated traces from stations '+id1[0]+\
                ' and '+id2[0],file=ofid)
        
        # Hor componets =======================================================
        #- Through all channels
        #======================================================================   
        elif comp=='RT':
            id1_T=str1_T[0].id
            id1_R=str1_R[0].id
            id2_T=str2_T[0].id
            id2_R=str2_R[0].id
            
            
            (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_T,\
                str2_T,corrname,geoinf)
            
            # Component TT    
            if npcc != 0:
                savecorrs(pcc,cstack_pcc,npcc,id1_T,\
                    id2_T,geoinf,corrname,'pcc',dir)
            if nccc != 0:
                savecorrs(ccc,cstack_ccc,nccc,id1_T,\
                    id2_T,geoinf,corrname,'ccc',dir)
                
            del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc
            
            # Component RR
            (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_R,\
                str2_R,corrname,geoinf)
            if npcc != 0:
                savecorrs(pcc,cstack_pcc,npcc,id1_R,\
                    id2_R,geoinf,corrname,'pcc',dir)
            if nccc != 0:
                savecorrs(ccc,cstack_ccc,nccc,id1_R,\
                    id2_R,geoinf,corrname,'ccc',dir)
            del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc
            
            if mix_cha == True:
            # Get the remaining component combinations
            # Component T1R2
                (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_T,\
                    str2_R,corrname,geoinf)
                if npcc != 0 or nccc != 0:
                    savecorrs(ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc,id1_T,\
                    id2_R,geoinf,corrname,dir)
                del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc
            # Component R1T2
                (ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc)=corr_pairs(str1_R,\
                    str2_T,corrname,geoinf)
                if npcc != 0:
                    savecorrs(pcc,cstack_pcc,npcc,id1_R,\
                    id2_T,geoinf,corrname,'pcc',dir)
                if nccc != 0:
                    savecorrs(ccc,cstack_ccc,nccc,id1_R,\
                    id2_T,geoinf,corrname,'ccc',dir)
                del ccc,pcc,cstack_ccc,cstack_pcc,nccc,npcc


def parlistpairs(corrname):
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
    # input...
    infile=inp.idfile
    nf=inp.npairs
    corrtype=inp.corrtype
    auto=inp.autocorr
    
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
                idlist[i] + '???.' + idlist[j] + '???.'+corrtype+'.' + corrname + '.SAC'
            else:
                fileid = cfg.datadir + 'correlations/' + corrname + '/' + \
                idlist[j] + '???.' + idlist[i] + '???.'+ corrtype + '.' + corrname + '.SAC'
            
            if glob(fileid) != [] and inp.update == True:
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
    
    
def corr_pairs(str1,str2,corrname,geoinf):
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

    startday, UTCDateTime object: Time where stack should start (if data 
    available)
    endday, UTCDateTime object: Maximum time until where stacking should be 
    carried out
    Fs, float: Sampling rate of data in Hz
    fmin: Minimum frequency the signal was filtered to (or 0.001 Hz)
    fmax: Maximum frequency the signal was filtered to (or 0.5*sampling 
    frequency)
    onebit: Boolean, do one-bitting or not
    verbose, boolean: loud or quiet
    
    output:
    
    cccstack, numpy array: classical cross correlation in time domain
    pccstack, numpy array: phase cross correlation
    ccccnt, int: Number of windows stacked for cccstack
    pcccnt, int: Number of windows stacked for pccstack
    
    
    """
    
    
   
    startday=obs.UTCDateTime(inp.startdate)
    endday=obs.UTCDateTime(inp.enddate)
    t1=startday
    Fs_new=inp.Fs
    tlen=int(inp.max_lag*Fs_new[-1])*2+1
    
    # Initialize arrays and variables
    pcccnt=0
    ccccnt=0
    n1=0
    n2=0
    cccstack=np.zeros(tlen)
    pccstack=np.zeros(tlen)
    cstack_ccc=np.zeros(tlen)
    cstack_pcc=np.zeros(tlen)
    
    while n1<len(str1) and n2<len(str2):
    
        # Check if the end of one of the traces has been reached
        if str1[n1].stats.endtime-t1<inp.winlen-1:
            n1+=1
            #print('No more windows in trace 1..',file=None)
            continue
        elif str2[n2].stats.endtime-t1<inp.winlen-1:
            n2+=1
            #print('No more windows in trace 2..',file=None)
            continue
        
        # Check the starttime of the potentially new trace
        t1=max(t1,str1[n1].stats.starttime,str2[n2].stats.starttime)
        #print(t1,file=None)
        t2=t1+inp.winlen
        
        # Check if the end of the desired stacking window is reached
        if t2>endday: 
            #print('At end of correlation time',file=None)
            break
        
        
        tr1=str1[n1].slice(starttime=t1,endtime=t2-1/Fs_new[-1])
        tr2=str2[n2].slice(starttime=t1,endtime=t2-1/Fs_new[-1])
        
        if tr1.stats.npts != tr2.stats.npts:
            t1 = t2 - inp.olap
            continue
        #==============================================================================
        #- Data treatment        
        #==============================================================================
        
        #- Downsampling ===============================================================
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
            t1 = t2 - inp.olap
            continue   

        #- Glitch correction ==========================================================
        if inp.cap_glitches == True:
            std1 = np.std(tr1.data*1.e6)
            gllow = inp.glitch_thresh * -std1
            glupp = inp.glitch_thresh * std1
            tr1.data = np.clip(tr1.data*1.e6,gllow,glupp)/1.e6
            
            std2 = np.std(tr2.data*1.e6)
            gllow = inp.glitch_thresh * -std2
            glupp = inp.glitch_thresh * std2
            tr2.data = np.clip(tr2.data*1.e6,gllow,glupp)/1.e6
            
        
            
#        #- Whitening            ==================================================================
        
        if inp.apply_white == True:
            df = 1/(tr1.stats.npts*tr1.stats.delta)
            print(df)
            freqaxis=np.fft.fftfreq(tr1.stats.npts,tr1.stats.delta)
            ind_fw1 = int(round(inp.white_freqs[0]/df))
            ind_fw2 = int(round(inp.white_freqs[1]/df))
            print(ind_fw1,ind_fw2)
            print(freqaxis[ind_fw1])
            print(freqaxis[ind_fw2])
            
            length_taper = int(round((inp.white_freqs[1]-inp.white_freqs[0])*\
            inp.white_tape/df))
            
            taper_left = np.linspace(0.,np.pi/2,length_taper)
            taper_left = np.square(np.sin(taper_left))
            
            taper_right = np.linspace(np.pi/2,np.pi,length_taper)
            taper_right = np.square(np.sin(taper_right))
            
            taper = np.zeros(tr1.stats.npts)
            taper[ind_fw1:ind_fw2] += 1.
            taper[ind_fw1:ind_fw1+length_taper] = taper_left
            taper[ind_fw2-length_taper:ind_fw2] = taper_right
            
            tr1.taper(max_percentage=0.05, type='cosine')
            tr2.taper(max_percentage=0.05, type='cosine')
            
            spec1 = np.fft.fft(tr1.data)
            spec2 = np.fft.fft(tr2.data)
            print(len(taper))
            print(len(spec1))
            spec1 /= np.abs(spec1)
            spec1 *= taper
            spec2 /= np.abs(spec2)
            spec2 *= taper
            
            
            tr1.data = np.real(np.fft.ifft(spec1,n=len(tr1.data)))
            tr2.data = np.real(np.fft.ifft(spec2,n=len(tr2.data)))
            
        #- One-bitting ================================================================
        
        if inp.apply_onebit == True:
            tr1.data = np.sign(tr1.data)
            tr2.data = np.sign(tr2.data)
        
            
        #==============================================================================
        #- Checks     
        #============================================================================== 
        
        if len(tr1.data) == len(tr2.data):
            mlag = inp.max_lag / tr1.stats.delta
            mlag=int(mlag)
            
            # Check if the traces are both long enough
            if len(tr1.data)<=2*mlag or len(tr2.data)<=2*mlag:
                t1 = t2 - inp.olap
                continue
            if tr1.data.any()==np.nan or tr2.data.any()==np.nan:
                t1 = t2 - inp.olap
                continue
            if tr1.data.any()==np.inf or tr2.data.any()==np.inf:
                t1 = t2 - inp.olap
                continue
        #==============================================================================
        #- Correlations proper 
        #==============================================================================         #- Taper
            if inp.taper_traces == True:
                tr1.taper(type='cosine',max_percentage=inp.perc_taper)
                tr2.taper(type='cosine',max_percentage=inp.perc_taper)
        
        #- Classical correlation part =====================================
            if inp.corrtype == 'ccc' or inp.corrtype == 'both':
                
                #ccc=classic_xcorr(tr1, tr2, mlag)
                (ccc, params) = cross_covar(tr1.data, \
                tr2.data, mlag)
                
                # normalization by trace energy
                en1 = params[2]
                en2 = params[3]
                if inp.normalize_correlation == True:
                    ccc/=(sqrt(en1)*sqrt(en2))
                
                cccstack+=ccc
                ccccnt+=1
                
                # Make this faster by zero padding
                
                
                if inp.get_pws == True:
                    coh_ccc = np.zeros(nextpow2(len(ccc)),dtype=np.complex)
                    coh_ccc[0:len(ccc)] += ccc*np.hanning(len(ccc))
                    coh_ccc = hilbert(ccc)
                    tol = np.max(coh_ccc)/1000.
                    if tol < 1e-9:
                        tol = 1e-9
                    coh_ccc = coh_ccc/(np.absolute(coh_ccc)+tol)
                    coh_ccc = coh_ccc[:len(ccc)]
                    cstack_ccc+=coh_ccc
                else: 
                    coh_ccc = None
                    cstack_ccc = None
                    
                if inp.write_all==True:
                    id1=str1[n1].id.split('.')[0]+'.'+str1[n1].id.split('.')[1]
                    id2=str2[n2].id.split('.')[0]+'.'+str2[n2].id.split('.')[1]
                    win_dir = cfg.datadir+'/correlations/interm/'+id1+\
                        '_'+id2+'/'
                    
                    if os.path.exists(win_dir)==False:
                        os.mkdir(win_dir)
                    
                    timestring = tr1.stats.starttime.strftime('.%Y.%j.%H.%M.%S')
                    savecorrs(ccc,coh_ccc,1,tr1.id,tr2.id,geoinf,\
                    corrname,'ccc',win_dir,params,timestring,startday=t1,endday=t2)
                            
                    
            #- Phase correlation part =========================================
            # To be implemented: Getting trace energy
            
            if inp.corrtype == 'pcc' or inp.corrtype == 'both':
                pcc=phase_xcorr(tr1, tr2, mlag, inp.pcc_nu)
                pccstack+=pcc
                pcccnt+=1
                
                
                
                if inp.get_pws == True:
                    coh_pcc = np.zeros(nextpow2(len(pcc)),dtype=np.complex)
                    coh_pcc[0:len(pcc)] += pcc*np.hanning(len(pcc))
                    coh_pcc = hilbert(pcc)
                    tol = np.max(coh_pcc)/1000.
                    if tol < 1e-9:
                        tol = 1e-9
                    coh_pcc = coh_pcc/(np.absolute(coh_pcc)+tol)
                    coh_pcc = coh_pcc[:len(pcc)]
                    cstack_pcc+=coh_pcc
                else: 
                    coh_pcc = None
                    cstack_pcc = None
                if inp.write_all==True:
                    id1=str1[n1].id.split('.')[0]+'.'+str1[n1].id.split('.')[1]
                    id2=str2[n2].id.split('.')[0]+'.'+str2[n2].id.split('.')[1]
                    win_dir = cfg.datadir+'/correlations/interm/'+id1+\
                        '_'+id2+'/'
                    
                    if os.path.exists(win_dir)==False:
                        os.mkdir(win_dir)
                    timestring = tr1.stats.starttime.strftime('.%Y.%j.%H.%M.%S')  
                    savecorrs(pcc,coh_pcc,1,tr1.id,tr2.id,geoinf,\
                    corrname,'pcc',win_dir,None,timestring,startday=t1,endday=t2)
                       
               
        #Update starttime
            t1 = t2 - inp.olap
        else:
            #print('Traces have unequal length!',file=None)
            t1 = t2 - inp.olap


    return(cccstack,pccstack,cstack_ccc,cstack_pcc,ccccnt,pcccnt)
    
    
    
    
    
    

def addtr(id):
    
    """
    Little reader.
    
    Needs to read all available data of one channel with a specified tag in a 
    directory.
    
    """
 
    traces=glob(inp.indir+'/'+id+'.*.'+inp.prepname+'.*')
    traces.sort()
    readone=False
    endday=obs.UTCDateTime(inp.enddate)
    
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
            if len(tr.data)*tr.stats.delta<inp.winlen-tr.stats.delta:
                continue
            
            
            #- Bandpass filter
            if inp.apply_bandpass == True:
                tr.filter('bandpass',freqmin=inp.filter[0],freqmax=inp.filter[1],\
                corners=inp.filter[2],zerophase=True)
                
            
            if readone==False:
                colltr=tr.copy()
                readone=True
                
            else:
                try:
                    colltr+=tr
                except TypeError:
                    continue   
     
    if readone == True:           
        return (colltr,readone)
    else:
        return(Trace(),readone)
    
   
def savecorrs(correlation,phaseweight,n_stack,id1,id2,geoinf,\
    corrname,corrtype,outdir,params=None,timestring='',startday=None,\
    endday=None):
    
    
    
#==============================================================================
    #- Write metadata info to sac header
    #- Store results
#==============================================================================

    
    tr=obs.Trace(data=correlation)
    tr.stats.sac={}
    
    if startday == None:
        startday=obs.UTCDateTime(inp.startdate)
    if endday == None:
        endday=obs.UTCDateTime(inp.enddate)
        
    (lat1, lon1, lat2, lon2, dist, az, baz)=geoinf
    
    
    # Add a preprocessing string:
    prepstring = ' '
    if inp.cap_glitches == True: prepstring += 'g'
    if inp.apply_white == True: prepstring += 'w'
    if inp.apply_onebit == True: prepstring += 'o'

    tr.stats.sampling_rate=inp.Fs[-1]
    tr.stats.starttime=obs.UTCDateTime(2000, 01, 01)-inp.max_lag*inp.Fs[-1]
    tr.stats.network=id1.split('.')[0]
    tr.stats.station=id1.split('.')[1]
    tr.stats.location=id1.split('.')[2]
    tr.stats.channel=id1.split('.')[3]
    
    tr.stats.sac['kt2']=prepstring
    tr.stats.sac['kt8']=corrtype
    tr.stats.sac['user0']=n_stack
    tr.stats.sac['user1']=inp.winlen
    tr.stats.sac['user2']=inp.olap
    tr.stats.sac['b']=-inp.max_lag
    tr.stats.sac['e']=inp.max_lag
    tr.stats.sac['kt0']=startday.strftime('%Y%j')
    tr.stats.sac['kt1']=endday.strftime('%Y%j')
    tr.stats.sac['iftype']=1
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
    
    if params is not None:
        tr.stats.sac['user3']=params[0]
        tr.stats.sac['user4']=params[1]
        tr.stats.sac['user5']=params[2]
        tr.stats.sac['user6']=params[3]
        tr.stats.sac['user7']=params[4]
        tr.stats.sac['user8']=params[5]
    
    
    #- open file and write correlation function
    fileid=outdir+id1+'.'+id2+'.'+corrtype+'.'+\
    corrname+timestring+'.SAC'
    tr.write(fileid,format='SAC')
    
    if phaseweight is not None:
        if corrtype == 'pcc': corr_type = 'pcs'
        if corrtype == 'ccc': corr_type = 'ccs'
        
        fileid_cwt=outdir+id1+'.'+id2+'.'+corrtype+\
        '.'+corrname+timestring+'.npy'
        np.save(fileid_cwt,phaseweight)
            
    
    
    
def classic_xcorr(trace1, trace2, max_lag_samples):
   
    x_corr = xcorr(trace1.data, trace2.data,\
        max_lag_samples, True)[2]
    
    return x_corr
    
def cross_covar(data1, data2, max_lag_samples, normalize_traces=True):
    
    # Remove mean and normalize; this should have no effect on the energy-normalized correlation result, but may avoid precision issues if trace values are very small
    if normalize_traces == True:
        data1-=np.mean(data1)
        data2-=np.mean(data2)
        data1/=np.max(np.abs(data1))
        data2/=np.max(np.abs(data2))
    
    # Make the data more convenient for C function np.correlate
    data1 = np.ascontiguousarray(data1, np.float32)
    data2 = np.ascontiguousarray(data2, np.float32)
    
    # Get the signal energy; most people normalize by the square root of that
    ren1 = np.correlate(data1,data1,mode='valid')[0]
    ren2 = np.correlate(data2,data2,mode='valid')[0]
    
    
    # Get the window rms
    rms1 = sqrt(ren1 / len(data1))
    rms2 = sqrt(ren2 / len(data2)) 
    
    # A further parameter to 'see' impulsive events: range of standard deviations
    nsmp = len(data1)/4
    std1 = np.zeros(4)
    std2 = np.zeros(4)
    for i in range(4):
        std1[i] = np.std(data1[i*nsmp:(i+1)*nsmp])
        std2[i] = np.std(data2[i*nsmp:(i+1)*nsmp])

    rng1 = np.max(std1)/np.mean(std1)
    rng2 = np.max(std2)/np.mean(std2)
    
    
    # Obtain correlation via FFT and IFFT
    ccv = np.correlate(data1,data2,mode='full')
    
    # Cut out the desired samples from the middle...
    i1 = (len(ccv) - (2*max_lag_samples+1))/2
    i2 = i1 + 2 * max_lag_samples + 1
    
    params = (rms1,rms2,ren1,ren2,rng1,rng2)
    
    return ccv[i1:i2], params
    
    

