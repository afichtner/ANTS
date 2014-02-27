#Stacking correlation windows

import os
import TOOLS.read_xml as rxml
import antconfig as cfg

from obspy import read,  Trace
from glob import glob

def st_res(indirs, corrname, outdir,num_win, format='SAC',stackall=False,verbose=False):
    
    """
    Cross-correlation and phase stack results from specified directories are added up to a 
    specified correlation time length and saved in a new directory
    indirs: python list object, containing the directories where to look for correlation files
    corrname: string, the correlation run which to filter for
    outdir: string, name of directory to store the result in 
    num_win: integer, number of single-window correlations that should be stacked.
    stackall: Boolean, if False, the script exits after obtaining only one stack of the desired window length.
    """
    
    #- output directory ======================================================
    if os.path.exists(outdir)==False:
        os.mkdir(outdir)
    #- Make sure stacks don't get added twice! If indir and outdir are not the same and files aren't reopened for stacking, this shouldnt be a problem though
    
    
    #- loop through the list of directories ==================================
    #- For now I am ignoring the phase weights. Too complicated. Get there later.
    
    for dir in indirs:
        
        pccs=glob(dir+'/*.pcc.'+corrname+'.'+format)
        pccs.sort()
        print pccs
        
        print 'Number of pcc files: '
        print len(pccs)
        
        cccs=glob(dir+'/*.ccc.'+corrname+'.'+format)
        cccs.sort()
        
        print 'Number of ccc files: '
        print len(cccs)
        
    if len(pccs)>=num_win:
        stackup(pccs,num_win,outdir,stackall)
        
    if len(cccs)>=num_win:
        stackup(cccs,num_win,outdir,stackall)
        
        
        
def stackup(filelist,num_win,outdir,stackall):
    
    tcnt=0
    
    while tcnt<=len(filelist)-num_win:
        
        cnt=0
    
        print '========Stacking:=============================================================='
        
        
        while cnt<num_win:
            
            try:
                tr=read(filelist[tcnt])
                
                print filelist[tcnt]
                
                if cnt==0:
                    filename=outdir+'/'+filelist[tcnt].split('/')[-1].rstrip('SAC')+str(num_win)+'.SAC'
                    stack=tr[0]
                else:
                    stack.data+=tr[0].data
                cnt+=1
                tcnt+=1
                
            except (IOError,TypeError):
                print 'Some Error occured.'
                tcnt+=1  
        
        print num_win
        stack.data/num_win
        stack.write(filename=filename,format='SAC')
        if stackall==False: return()
            
# a very very simple get-the-convergence script, just to save some lines on plotting
# Each stack length only contained once in the directory, not several times     
            
def convergetest(indir,corrname,start,dt):
    import numpy as np
    import matplotlib.pyplot as plt
    from obspy.core import UTCDateTime
    
    
    
    lens=[]
    mf=[]
    ccoef=[]
    
    
    filelist=glob(indir+'*'+corrname+'*.SAC')
    for file in filelist:
        lens.append(int(file.split('.')[-2]))
      
    ref=filelist[np.argmax(lens)]
    print 'reference: ',ref
    
    reftr=read(ref)[0]
    print reftr.stats.starttime
    print reftr.stats.endtime
    print start+dt
    reftr.trim(starttime=UTCDateTime(1970,01,01)+start,endtime=UTCDateTime(1970,01,01)+start+dt)
    reftr=reftr.data
    
    filelist=[]
    filelist=glob(indir+'*'+corrname+'*.SAC')
    for file in filelist:
        tr=read(file)[0]
        tr.trim(starttime=UTCDateTime(1970,01,01)+start,endtime=UTCDateTime(1970,01,01)+start+dt)
        tr=tr.data
        #plt.plot(reftr)
        #plt.plot(tr)
        #plt.show()
        mf.append(np.sum(np.power(np.abs(tr-reftr),2)))
        ccoef.append(np.corrcoef(tr,reftr)[0,1])
       
    
    plt.plot(lens,mf,'*',linewidth=2)
    plt.xlabel('Number of time windows')
    plt.ylabel('L2 misfit')
    plt.show()
    
    plt.plot(lens,ccoef,'d',linewidth=2)
    plt.xlabel('Number of time windows')
    plt.ylabel('Correlation coefficient')
    plt.show()
    
    
    
    
        
    
    
    
    
    
    
            
            
            
            
            
            
            
        
       
        
        
        
