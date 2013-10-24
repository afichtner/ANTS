#Stacking monthly stacks....
import obspy as obs
import os
import sys
import re
from shutil import copy


def st_res(list_dirs, outdir, fileform, verbose=False):
    
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
                            trace1=obs.read(infile)[0]
                            trace2=obs.read(outfile)[0]
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
                        cmd=('cat '+infile+' >> '+outfile)
                        os.system(cmd)
                        
                        
                
