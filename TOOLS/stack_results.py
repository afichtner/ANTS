#Stacking monthly stacks....
from obspy import read,  Trace
import os
import sys
import re
from shutil import copy
from glob import glob

def st_res(list_dirs, fileform, outdir, verbose=False):
    
    """
    Cross-correlation and phase stack results from specified directories are added and saved in a new directory
    list_dirs: list, full path all the directories to be searched (strings)
    outdir: string, name of directory to store the result in 
    fileform: format of data files to look for (e. g. MSEED)
    """
    #dirfile=open(list_dirs, 'r')
    #dirs=dirfile.read().split('\n')
    
    #- output directory ======================================================
    if os.path.exists(outdir)==False:
        os.mkdir(outdir)
    else:
        print 'Target directory exists already! Stacks may not be added twice. Aborting.'
        return
        
    #- Logfile ===============================================================
    #logf=open(outdir+'/'+list_dirs.split('/')[-1]+'.txt', 'w')
    #- loop through the list of directories ==================================
    
    for dir in list_dirs:
        
        content=os.listdir(dir)
        
    
    #- for each file there check if it is in other directories, too ==========
    
        for file in content:
            #if verbose: print 'Opening file '+file
            #Only start with correlation files, not metadata or coherence stacks
            
            inf=file.split('.')
            
            #rudimentary check to exclude some files that are also there, like plot pngs
            if 'metadata' in inf or 'MSEED' in inf:
                
                
                if verbose: 
                    print '---------------------------------------------------------'
                    print 'Stacking ', file
                
                infile=os.path.join(dir, file)
                outfile=os.path.join(outdir, file)
                
                if os.path.exists(outfile)==False:
                    copy(infile, outdir)
                else:
                    if inf[7]=='metadata':
                        continue
                    
                    elif inf[8]=='MSEED':  
                        try:
                            trace1=read(infile)[0]
                            trace2=read(outfile)[0]
                        except (IOError, TypeError):    
                            if verbose: print 'Wrong file type. Skipping this file.'
                            continue
              
                        #Add traces       
                        trace1.data+=trace2.data
                        try:
                            
                            trace1.write(outfile, format='MSEED')
                            
                        except NotImplementedError:
                            print 'Warning: NonImplementedError occured.\n One or both files may contain masked array.'
                            #logf.write('A problem occured with the stacking of '+dir+'/'+file+'\n')
                        
                    elif inf[8]=='metadata':
                        #Add windows
                        cmd=('cat '+infile+ ' | grep \'^[0-9, ]*$\'' + ' >> '+outfile)
                        os.system(cmd)
                        
                        
                
def cl_corr(indir, form, merge_loc=False, verbose=False):
    files=glob(indir+'/*.?cc_stack.'+form)
   
    for file in files:
        if verbose:
            print '---------------------------------------------------'
        
        corrtype=file.split('/')[-1].split('.')[-2].rstrip('stack').rstrip('_')
        corr=file.split('/')[-1].rstrip(form).rstrip('.'+corrtype+'_stack.')
        net1=file.split('/')[-1].split('.')[0]
        net2=file.split('/')[-1].split('.')[3].split('-')[1]
        sta1=file.split('/')[-1].split('.')[1]
        sta2=file.split('/')[-1].split('.')[4]
        loc1=file.split('/')[-1].split('.')[2]
        loc2=file.split('/')[-1].split('.')[5]
        cha1=file.split('/')[-1].split('.')[3].split('-')[0]
        cha2=file.split('/')[-1].split('.')[6]
       
        if verbose:
            print sta1,  sta2
        if sta1==sta2: continue
        
        #- If different locations should be merged, wildcard the location
        comp=[]
        if merge_loc==True:
            corr_loc=net1+'.'+sta1+'.*.'+cha1+'-'+net2+'.'+sta2+'.*.'+cha2
            corr_rev=net2+'.'+sta2+'.*.'+cha2+'-'+net1+'.'+sta1+'.*.'+cha1
            comp+=(glob(indir+'/'+corr_loc+'.'+corrtype+'_stack.'+form))
            comp+=(glob(indir+'/'+corr_rev+'.'+corrtype+'_stack.'+form))
            
        else:
            corr_rev=corr.split('-')[1]+'-'+corr.split('-')[0]
            comp+=(glob(indir+'/'+corr_rev+'.'+corrtype+'_stack.'+form))
      
        
        
        
        
        
        for file2 in comp:
            
            if file2==file:
                continue
            
            if os.path.exists(file)==False: 
                continue
            
            compname=file2.split('/')[-1].rstrip(form).rstrip('.'+corrtype+'_stack.')
            
            if verbose:
                print 'Merging '+corr+' and '+compname
            
            #- Merge the data themselves
            tr1=read(file)[0]
            tr2=read(file2)[0]
            
            if file2.split('/')[-1].split('.')[1]==sta2:
                tr1.data+=tr2.data[::-1]
            elif file2.split('/')[-1].split('.')[1]==sta1:
                tr1.data+=tr2.data
            else:
                print 'An Error occured. Go to lunch.'
                continue
                
            tr1.write(file, form)
            os.system('rm '+file2)
            
            #- Merge the coherence stacks
            coh_re1=read(indir+'/'+corr+'.coherence_stack_real.'+form)[0]
            coh_im1=read(indir+'/'+corr+'.coherence_stack_imag.'+form)[0]
            
            coh_re2=read(indir+'/'+compname+'.coherence_stack_real.'+form)[0]
            coh_im2=read(indir+'/'+compname+'.coherence_stack_imag.'+form)[0]
            
            coh_re1.data+=coh_re2.data[::-1]
            coh_im1.data+=coh_im2.data[::-1]
            
            coh_re1.write(indir+'/'+corr+'.coherence_stack_real.'+form, form)
            coh_im1.write(indir+'/'+corr+'.coherence_stack_imag.'+form, form)
            
            os.system('rm '+indir+'/'+compname+'.coherence_stack_real.'+form)
            os.system('rm '+indir+'/'+compname+'.coherence_stack_imag.'+form)
            
            #- Merge the metadata
            md1=indir+'/'+corr+'.'+corrtype+'.metadata'
            md2=indir+'/'+compname+'.'+corrtype+'.metadata'
            cmd=('cat '+md2+ ' | grep \'^[0-9, ]*$\'' + ' >> '+md1)
            os.system(cmd)
            os.system('rm '+md2)
         
            
            
            
            
            
            
            
            
            
            
        
       
        
        
        
