from glob import glob
import antconfig as cfg

def prepare_corr_ids(indir,prepname):
    
    outfile='INPUT/ID_LISTS/ids_'+prepname+'_corr.txt'
    fid = open(outfile,'w')
    
    idlist = glob(indir+'/*.'+prepname+'.mseed')
    idlist += glob(indir+'/*.'+prepname+'.MSEED')
    idlist += glob(indir+'/*.'+prepname+'.sac')
    idlist += glob(indir+'/*.'+prepname+'.SAC')
    
    ids_un = list()
    
    for id in idlist:
        
        id = id.split('/')[-1]
        
        id = id.split('.')[0]+'.'+id.split('.')[1]+'.'+id.split('.')[2]+'.'
        
        if id not in ids_un:
            print id
            fid.write(id+'\n')
            ids_un.append(id)
            
    fid.close()
    return()
        
        
    
    
    
    
    


