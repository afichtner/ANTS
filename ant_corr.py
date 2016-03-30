from __future__ import print_function
from mpi4py import MPI

import time
import sys
import os
import numpy as np 
 
#import obspy as obs
from ANTS.TOOLS import read_xml as rxml
from ANTS import antconfig as cfg
from ANTS.TOOLS import processing as proc
from ANTS.TOOLS import rotationtool as rt
from ANTS.INPUT import input_correlation as inp

from math import sqrt
from glob import glob
from obspy.core import Stats, Trace, Stream, UTCDateTime, read
#from obspy.noise.correlation import Correlation
#from obspy.noise.correlation_functions import phase_xcorr
from obspy.signal.cross_correlation import xcorr
from obspy.signal.filter import envelope
from obspy.signal.util import nextpow2
from obspy.signal.tf_misfit import cwt
from scipy.signal import hilbert
from scipy import fftpack

if __name__=='__main__':
    from ANTS import ant_corr as pc
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    #rank = int(os.environ[inp.rankvariable])
    #size = int(sys.argv[1])
    
    if rank==0 and inp.update == False and os.path.exists(cfg.datadir+\
    '/correlations/input/'+inp.corrname+'.txt') == True:
        print('Choose a new correlation name tag or set update=True.\
         Aborting all processes.',file=None)
        # sys.exit did NOT work properly here. 
        MPI.COMM_WORLD.Abort(1)
    pc.par_st(size, rank)
    

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
    #n1=int(len(idpairs)/size)
    #n2=len(idpairs)%size
    #ids=list()
    
    #for i in range(0,n1):
    #    ids.append(idpairs[i*size+rank])
    #if rank<n2:
    #    ids.append(idpairs[n1*size+rank])
    ids = idpairs[rank:len(idpairs):size]
    
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
        
        corrblock(block,dir,corrname,rank,ofid)
        if rank==0:
            print('Finished a block of correlations',file=None)
            print(time.strftime('%H.%M.%S'),file=None)
        
        # Flush the outfile buffer every now and then...
        if inp.verbose==True:
	    ofid.flush()
    
    print('\nTrying to move computed calculations from: ',file=None)
    print(dir+'* ',file=None)
    print('to:',file=None)
    print(cfg.datadir+'correlations/'+corrname+'/\n\n',file=None)
    
    dir1 = os.path.join(dir,'*')
    dir2 = os.path.join(cfg.datadir,'correlations',corrname)
    os.system('mv '+dir1+' '+dir2)
    os.system('rmdir '+dir)
    print('Rank %g finished correlations.' %rank,file=None)
        
def corrblock(block,dir,corrname,rank,ofid=None):
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
    print('Rank %g: Working on a block of station pairs...\n' %rank,file=None)
    
    
    datstr=Stream()
    idlist=list()
    verbose=inp.verbose
    
#==============================================================================
    #- Get some information needed for the cross correlation   
#============================================================================== 
    
    
    cha=inp.channel
    comp=inp.components
    mix_cha=inp.mix_cha
    

    for pair in block:
        str1=Stream()
        str2=Stream()
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
                (colltr,readsuccess) = addtr(id,rank)
        
                #- add this entire trace (which contains all data of this 
                #- station that are available in this directory) to datstr and 
                #- update the idlist
                if readsuccess == True:
                    datstr += colltr
                    str1 += colltr.split()
                    idlist += id
                    
                    if verbose:
                        print('Read in traces for channel '+id,file=ofid)
                    del colltr
                else:
                    if verbose:
                        print('No traces found for channel '+id,file=ofid)
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
                    (colltr,readsuccess)=addtr(id,rank)
                    
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
            (lat1,lon1) = rxml.get_coord_staxml(id1[0].\
            split('.')[0],id1[0].split('.')[1])
            (lat2,lon2) = rxml.get_coord_staxml(id2[0].\
            split('.')[0],id2[0].split('.')[1])
                        
            if (lat1,lon1) == (0,0) or (lat2,lon2) == (0,0):
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
            if nccc != 0 or npcc != 0 and inp.verbose:
                print('Correlated traces from stations '+id1[0]+\
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
                fileid1 = cfg.datadir + 'correlations/' + corrname + '/rank*/' +\
                idlist[i] + '???.' + idlist[j] + '???.'+corrtype+'.' + corrname + '.SAC'
            else:
                fileid = cfg.datadir + 'correlations/' + corrname + '/' + \
                idlist[j] + '???.' + idlist[i] + '???.'+ corrtype + '.' + corrname + '.SAC'
                fileid1 = cfg.datadir + 'correlations/' + corrname + '/rank*/' + \
                idlist[j] + '???.' + idlist[i] + '???.'+ corrtype + '.' + corrname + '.SAC'
            
            if glob(fileid) != [] and inp.update == True:
                print('Correlation already available, continuing...')
                continue
            if glob(fileid1) != [] and inp.update == True:
                print('Correlation already available, continuing...')
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
    
    print('Computing correlation stack for:',file=None)
    print('-------------',file=None)
    print(str1[0].id)
    print(str2[0].id)
    print('-------------',file=None)
    
   
    startday=UTCDateTime(inp.startdate)
    endday=UTCDateTime(inp.enddate)
    t1=startday
    Fs_new=inp.Fs
    tlen=int(inp.max_lag*Fs_new[-1])*2+1
    
    # Initialize arrays and variables
    pcccnt=0
    ccccnt=0
    n1=0
    n2=0
    cccstack=np.zeros(tlen)
    pccstack=np.zeros(tlen,dtype=np.float64)
    cstack_ccc=np.zeros(tlen,dtype=np.complex128)
    cstack_pcc=np.zeros(tlen,dtype=np.complex128)
    
    # Collect intermediate traces in a binary file.
    if inp.write_all:
        if inp.get_pws:
            msg = 'Saving intermediate windows of phase weighted stack\
is not implemented yet. Intermediate windows will be saved as linear stack.'
        
        # set the size of float
        
        # set size of character array for processing string
        
        # format of file:
        # header consisting of float, float, float, string of 256, string of 256
        # traces, traces, traces...
        # header values: Sampling rate Fs, number of samples in each trace,  Nr. of windows in intermediate stacks, Endianness, preprocessing string
        interm_fs = Fs_new[-1]
        interm_nsam = tlen
        interm_nwin = inp.interm_nstack
        
        if cccstack.dtype.byteorder == '=':
            interm_endian = sys.byteorder
        elif cccstack.dtype.byteorder == '<':
            interm_endian = 'little'
        elif cccstack.dtype.byteorder == '>':
            interm_endian = 'big'
            
        interm_preproc = get_prepstring()
        
        # open the file(s)
        if inp.corrtype in ['both','pcc','ccc']:
            
            outdir = os.path.join(cfg.datadir,'correlations',inp.corrname)
            interm_file=os.path.join(outdir,str1[0].id+'.'+str2[0].id+'.'+inp.corrtype+'.'+\
            inp.corrname+'.windows.bin')
            interm_file = open(interm_file,'wb')
            header_1 = np.array([interm_fs,interm_nsam,interm_nwin],dtype='f4')
            header_2 = np.array([interm_endian,interm_preproc],dtype='S256')
            header_1.tofile(interm_file)
            header_2.tofile(interm_file)
            
        else:
            print('Correlation type not recognized. Correlation types are:\
ccc, pcc or both.')
            MPI.COMM_WORLD.Abort(1)
            
            
         
    while n1<len(str1) and n2<len(str2):
    
        # Check if the end of one of the traces has been reached
        if str1[n1].stats.endtime-t1<inp.winlen-1:
            n1+=1
            print('No more windows in trace 1..',file=None)
            continue
        elif str2[n2].stats.endtime-t1<inp.winlen-1:
            n2+=1
            print('No more windows in trace 2..',file=None)
            continue
        
        # Check the starttime of the potentially new trace
        t1=max(t1,str1[n1].stats.starttime,str2[n2].stats.starttime)
        #print(t1,file=None)
        t2=t1+inp.winlen
        print(t1)
        # Check if the end of the desired stacking window is reached
        if t2>endday: 
            #print('At end of correlation time',file=None)
            break
        
        
        tr1=str1[n1].slice(starttime=t1,endtime=t2-1/Fs_new[-1])
        tr2=str2[n2].slice(starttime=t1,endtime=t2-1/Fs_new[-1])
        
        if tr1.stats.npts != tr2.stats.npts:
            t1 = t2 - inp.olap
            
            continue
       # tr1.plot()
        #tr2.plot()
        
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
        #==============================================================================
        #- Checks     
        #============================================================================== 
        if tr1.data.any()==np.nan or tr2.data.any()==np.nan:
            t1 = t2 - inp.olap
            print('Encountered nan, skipping this trace pair...',file=None)
            continue
        if tr1.data.any()==np.inf or tr2.data.any()==np.inf:
            t1 = t2 - inp.olap
            print('Encountered inf, skipping this trace pair...',file=None)
            continue
            
        if len(tr1.data) == len(tr2.data):
            mlag = inp.max_lag / tr1.stats.delta
            mlag=int(mlag)
            
        # Check if the traces are both long enough
        if len(tr1.data)<=2*mlag or len(tr2.data)<=2*mlag:
            t1 = t2 - inp.olap
            print('One or both traces too short',file=None)
            continue
        # Check if too many zeros
        # I use epsilon for this check. That is convenient but not strictly right. It seems to do the job though. min doesn't work.
        
        if np.sum(np.abs(tr1.data)<sys.float_info.epsilon) > 0.1*tr1.stats.npts or \
        np.sum(np.abs(tr2.data)<sys.float_info.epsilon) > 0.1*tr2.stats.npts:
            t1 = t2 - inp.olap
            
            if inp.verbose: print('More than 10% of trace equals 0, skipping.',file=None)
            continue
        
         #==============================================================================
        #- Data treatment        
        #==============================================================================


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
        
        if inp.apply_white:
            tr1 = whiten(tr1)
            tr2 = whiten(tr2)
            
        #- One-bitting ================================================================
        
        if inp.apply_onebit:
            tr1.data = np.sign(tr1.data)
            tr2.data = np.sign(tr2.data)
        
#- RAM normalization...who wants to do all this stuff!! ================================================================
        if inp.apply_ram:
            tr1 = ram_norm(tr1,inp.ram_window,prefilt=inp.ram_filter)
            tr2 = ram_norm(tr2,inp.ram_window,prefilt=inp.ram_filter)
        
       
        #==============================================================================
        #- Correlations proper 
        #==============================================================================        #- Taper
        if inp.taper_traces == True:
            tr1.taper(type='cosine',max_percentage=inp.perc_taper)
            tr2.taper(type='cosine',max_percentage=inp.perc_taper)
        
    #-   Classical correlation part =====================================
        if inp.corrtype == 'ccc' or inp.corrtype == 'both':
            #ccc=classic_xcorr(tr1, tr2, mlag)
            (ccc, params) = cross_covar(tr1.data, \
            tr2.data, mlag,inp.normalize_correlation)
            
            
            
            if ccc.any() == np.nan:
                msg='NaN encountered, omitting correlation from stack.'
                warn(msg)
                print(tr1)
                print(tr2)
                t1 = t2 - inp.olap
                continue
                
            # normalization by trace energy
            en1 = params[2]
            en2 = params[3]
            if inp.normalize_correlation:
                ccc/=(sqrt(en1)*sqrt(en2))
            
            cccstack+=ccc
            ccccnt+=1
           
            print('Finished a correlation window',file=None)
            # Make this faster by zero padding
            
            
            if inp.get_pws == True:
                coh_ccc = np.zeros(nextpow2(len(ccc)))
                startindex = int(0.5*(len(coh_ccc) - len(ccc)))
                coh_ccc[startindex:startindex+len(ccc)] += ccc*np.hanning(len(ccc))
                coh_ccc = hilbert(coh_ccc)
                tol = np.max(coh_ccc)/10000.
                #if tol < 1e-9:
                #    tol = 1e-9
                coh_ccc = coh_ccc/(np.absolute(coh_ccc)+tol)
                coh_ccc = coh_ccc[startindex:startindex+len(ccc)]
                cstack_ccc+=coh_ccc
                
            else: 
                coh_ccc = None
                cstack_ccc = None
                
            if inp.write_all==True and ccccnt % inp.interm_nstack == 0:
                trcname = t2.strftime("end%Y.%j.%H.%M.%S")
                trcname = np.array([trcname],dtype='S24')
                trcname.tofile(interm_file)
                ccc = np.array(ccc,dtype='f4')
                ccc.tofile(interm_file)
                
                #id1=str1[n1].id.split('.')[0]+'.'+str1[n1].id.split('.')[1]
                #id2=str2[n2].id.split('.')[0]+'.'+str2[n2].id.split('.')[1]
                #win_dir = cfg.datadir+'/correlations/interm/'+id1+\
                #    '_'+id2+'/'
                
                #if os.path.exists(win_dir)==False:
                #    os.mkdir(win_dir)
                
                #timestring = tr1.stats.starttime.strftime('.%Y.%j.%H.%M.%S')
                #savecorrs(ccc,coh_ccc,1,tr1.id,tr2.id,geoinf,\
                #corrname,'ccc',win_dir,params,timestring,startday=t1,endday=t2)
                        
                
                #- Phase correlation part =========================================
                # To be implemented: Getting trace energy
        
        elif inp.corrtype == 'pcc' or inp.corrtype == 'both':
            pcc=phase_xcorr(tr1.data, tr2.data, mlag, inp.pcc_nu)
            pccstack+=pcc
            pcccnt+=1
            
            
            
            if inp.get_pws == True:
                coh_pcc = np.zeros(nextpow2(len(pcc)))
                startindex = int(0.5*(len(coh_pcc) - len(pcc)))
                coh_pcc[startindex:startindex+len(pcc)] += pcc*np.hanning(len(pcc))  # Tapering and zero padding to make hilbert trafo faster
                coh_pcc = hilbert(coh_pcc)
                tol = np.max(coh_pcc)/10000.
                #if tol < 1e-9:
                #    tol = 1e-9
                coh_pcc = coh_pcc/(np.absolute(coh_pcc)+tol)
                coh_pcc = coh_pcc[startindex:startindex+len(pcc)]
                cstack_pcc+=coh_pcc
            else: 
                coh_pcc = None
                cstack_pcc = None
            if inp.write_all==True:
                trcname = t2.strftime("end%Y.%j.%H.%M.%S")
                trcname = np.array([trcname],dtype='S24')
                trcname.tofile(interm_file)
                pcc = np.array(pcc,dtype='f4')
                pcc.tofile(interm_file)
                #id1=str1[n1].id.split('.')[0]+'.'+str1[n1].id.split('.')[1]
                #id2=str2[n2].id.split('.')[0]+'.'+str2[n2].id.split('.')[1]
                #win_dir = cfg.datadir+'/correlations/interm/'+id1+\
                #    '_'+id2+'/'
                #
                #if os.path.exists(win_dir)==False:
                #    os.mkdir(win_dir)
                #timestring = tr1.stats.starttime.strftime('.%Y.%j.%H.%M.%S')  
                #savecorrs(pcc,coh_pcc,1,tr1.id,tr2.id,geoinf,\
                #corrname,'pcc',win_dir,None,timestring,startday=t1,endday=t2)
                   
           
        #Update starttime
        t1 = t2 - inp.olap
    if 'interm_file' in locals():  
        interm_file.close()
    return(cccstack,pccstack,cstack_ccc,cstack_pcc,ccccnt,pcccnt)
    
    
    
    
    
    

def addtr(id,rank):
    
    """
    Little reader.
    
    Needs to read all available data of one channel with a specified tag in a 
    directory.
    
    """
    print('Rank %g: Reading noise data...\n' %rank,file=None)
    traces=glob(inp.indir+'/'+id+'.*.'+inp.prepname+'.*')
    traces.sort()
    readone=False
    endday=UTCDateTime(inp.enddate)
    startday=UTCDateTime(inp.startdate)
    if len(traces) == 0:
        return (Trace(),False)
             
    #- collect a trace with masked samples where gaps are.
    #- This is convenient because we then get only one trace that can be 
    #- handled with an index 
    #- inside the datstr objects more easily (rather than having a stream with 
    #- variable number of traces)
    for filename in traces: 
        
        (sy,sm)=filename.split('/')[-1].split('.')[4:6]
        (ey,em)=filename.split('/')[-1].split('.')[9:11]
        sf=UTCDateTime(sy+','+sm)
        ef=UTCDateTime(ey+','+em)
        if sf>endday or ef<startday:
            continue
        
        try:
            newtr=read(filename)
            print(newtr[0].stats.starttime)
        except:
            print('Problems opening data file:\n',file=None)
            print(tr,file=None)
            continue
        
        for tr in newtr:
            #- Check if at least one window contained
            if len(tr.data)*tr.stats.delta<inp.winlen-tr.stats.delta:
                continue
            
            
            #- Bandpass filter
            if inp.apply_bandpass:
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

    
    tr=Trace(data=correlation)
    tr.stats.sac={}
    
    if startday == None:
        startday=UTCDateTime(inp.startdate)
    if endday == None:
        endday=UTCDateTime(inp.enddate)
        
    (lat1, lon1, lat2, lon2, dist, az, baz)=geoinf
    
    
    # Add a preprocessing string:
    prepstring = get_prepstring()
    

    tr.stats.sampling_rate=inp.Fs[-1]
    tr.stats.starttime=UTCDateTime(2000, 01, 01)-inp.max_lag*inp.Fs[-1]
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
        
        fileid_cwt=outdir+id1+'.'+id2+'.'+corrtype+\
        '.'+corrname+timestring+'.npy'
        np.save(fileid_cwt,phaseweight)
            
    
    
    
def classic_xcorr(trace1, trace2, max_lag_samples):
   
    x_corr = xcorr(trace1.data, trace2.data,\
        max_lag_samples, True)[2]
    
    return x_corr
    
def cross_covar(data1, data2, max_lag_samples, normalize_traces):
    
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
    nsmp = int(len(data1)/4)
    std1 = np.zeros(4)
    std2 = np.zeros(4)
    for i in range(4):
        std1[i] = np.std(data1[i*nsmp:(i+1)*nsmp])
        std2[i] = np.std(data2[i*nsmp:(i+1)*nsmp])

    rng1 = np.max(std1)/np.min(std1)
    rng2 = np.max(std2)/np.min(std2)
    
    
    # Obtain correlation via FFT and IFFT
    ccv = np.correlate(data1,data2,mode='same')
    
    # Cut out the desired samples from the middle...
    i1 = (len(ccv) - (2*max_lag_samples+1))/2
    i2 = i1 + 2 * max_lag_samples + 1
    
    params = (rms1,rms2,ren1,ren2,rng1,rng2)
    
    
    return ccv[i1:i2], params
    
    
    
def whiten(tr):
    df = 1/(tr.stats.npts*tr.stats.delta)
    freqaxis=np.fft.fftfreq(tr.stats.npts,tr.stats.delta)
    
    ind_fw1 = int(round(inp.white_freqs[0]/df))
    ind_fw2 = int(round(inp.white_freqs[1]/df))
    
    
    length_taper = int(round((inp.white_freqs[1]-inp.white_freqs[0])*\
    inp.white_tape/df))
    
    taper_left = np.linspace(0.,np.pi/2,length_taper)
    taper_left = np.square(np.sin(taper_left))
    
    taper_right = np.linspace(np.pi/2,np.pi,length_taper)
    taper_right = np.square(np.sin(taper_right))
    
    taper = np.zeros(tr.stats.npts)
    taper[ind_fw1:ind_fw2] += 1.
    taper[ind_fw1:ind_fw1+length_taper] = taper_left
    taper[ind_fw2-length_taper:ind_fw2] = taper_right
    
    tr.taper(max_percentage=0.05, type='cosine')
    
    
    
    spec = fftpack.fft(tr.data)
    
    # Don't divide by 0
    tol = np.max(np.abs(spec)) / 1e5
    spec /= np.abs(spec+tol)
    spec *= taper
   
    
    tr.data = np.real(fftpack.ifft(spec,n=len(tr.data)))
    return tr
    
    
def ram_norm(trace,winlen,prefilt=None):
    
    trace_orig = trace.copy()
    hlen = int(winlen*trace.stats.sampling_rate/2.)
    weighttrace = np.zeros(trace.stats.npts)
    
    if prefilt is not None:
        trace.filter('bandpass',freqmin=prefilt[0],freqmax=prefilt[1],\
        corners=prefilt[2],zerophase=True)
        
    envlp = envelope(trace.data)

    for n in xrange(hlen,trace.stats.npts-hlen):
        weighttrace[n] = np.sum(envlp[n-hlen:n+hlen+1]/(2.*hlen+1))
        
    weighttrace[0:hlen] = weighttrace[hlen]
    weighttrace[-hlen:] = weighttrace[-hlen-1]
    
    trace_orig.data /= weighttrace
    return(trace_orig)

def get_prepstring():
    
    prepstring = ''
    if inp.apply_bandpass: 
        prepstring +='b'
    else:
        prepstring += '-'
    if inp.cap_glitches: 
        prepstring += 'g'
    else:
        prepstring += '-'    
    if inp.apply_white: 
        prepstring += 'w'
    else:
        prepstring += '-'
    if inp.apply_onebit: 
        prepstring += 'o'
    else:
        prepstring += '-'
    if inp.apply_ram: 
        prepstring += 'r'
    else:
        prepstring += '-'
    
    return prepstring