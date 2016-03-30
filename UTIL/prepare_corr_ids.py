from glob import glob
from ANTS import antconfig as cfg
import os
import sys

if __name__=='__main__':
    from UTIL import prepare_corr_ids as pp
   
    indir = sys.argv[1]
    prepname = sys.argv[2]
    
    pp.prepare_corr_ids(indir,prepname)
    
def prepare_corr_ids(indir,prepname):
 
    
    outfile=os.path.join(cfg.inpdir,'correlationlist_new.txt')
    fid = open(outfile,'w')
    
    idlist = glob(indir+'/*.'+prepname+'.mseed')
    idlist += glob(indir+'/*.'+prepname+'.MSEED')
    idlist += glob(indir+'/*.'+prepname+'.sac')
    idlist += glob(indir+'/*.'+prepname+'.SAC')
    
    ids_un = list()
    
    print 'I found the following suitable channel IDs for correlation:'
    
    for id in idlist:
        
        id = id.split('/')[-1]
        
        id = id.split('.')[0]+'.'+id.split('.')[1]+'.'+id.split('.')[2]+'.'
        
        if id not in ids_un:
            print id
            fid.write(id+'\n')
            ids_un.append(id)
           
    print 'They are saved under:'
    print outfile
    fid.close()
        
        
    
    
    
    
    


