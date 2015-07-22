#- Script to quickly get a table of successfully processed data
from glob import glob
from obspy import UTCDateTime
import collections
import time
import sys

if __name__=='__main__':
    import listproc as lp
    if len(sys.argv)==4:
        lp.listproc(sys.argv[1],sys.argv[2],sys.argv[3])
    elif len(sys.argv)==5:
        lp.listproc(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    elif len(sys.argv)==6:
        lp.listproc(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    elif len(sys.argv)==7:        
        lp.listproc(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
        
    
def listproc(indir,oname,procname,station='*',location='*',channel='*'):
    """
    indir: Directory where processed data is located (string)
    oname: The information will be provided as ascii file in this location (string)
    procname: processing run name to check (string)
    """
    ofid=open(oname,'w')
    ofid.write('Nr ID               t1                   t2                     gap (hours) \n' )
    ofid.write('===========================================================================\n')
    files=glob(indir+'/*.'+procname+'.*')
    files.sort()
    toc=dict()
    foundgap = False
    foundover = False
    
    for file in files:
        file=file.split("/")[-1]
        file=file.split('.')
        
        if station!='*' and file[1]!=station:
            continue
        if location!='*' and file[2]!=location:
            continue
        if channel!='*' and file[3]!=channel:
            continue
        
        id=file[0]+'.'+file[1]+'.'+file[2]+'.'+file[3]
        t1_new=UTCDateTime(file[4]+','+file[5]+','+file[6]+','+file[7]+','+file[8])
        t2_new=UTCDateTime(file[9]+','+file[10]+','+file[11]+','+file[12]+','+file[13])
        
        if id in toc:
            gap=toc[id][2]
            
            if t1_new>t2_old:
                if t1_new-t2_old>1:
                    gap+=t1_new-t2_old-1
                    foundgap = True
                toc[id][1] = t2_new
                toc[id][2] = gap 
                
            elif t2_new<t1_old:
                if t1_old-t2_new>1:
                    gap+=t1_old-t2_new-1

                toc[id][0] = t1_new
                toc[id][2] = gap
            
            elif t1_new < t2_old and t1_new > t1_old:
                foundover = True
                toc[id][1] = t2_new
                
        else:
            toc.update({id:[t1_new,t2_new,0]})
          
        t1_old = t1_new
        t2_old = t2_new
    
          
    toc=collections.OrderedDict(sorted(toc.items()))   
    k=1
    
    for item in toc:
        ofid.write("\n%g%15s  %s    %s    %g" %(k,item,toc[item][0].strftime('%Y.%m.%d.%H.%M.%S'),\
        toc[item][1].strftime('%Y.%m.%d.%H.%M.%S'),round(toc[item][2]/3600,4)))
        k+=1
    ofid.close()
    
    if foundgap == True and foundover == True:
        print '\nChecked %g files, found gaps and overlaps, see output\
file %s for details.\n' %(len(files),oname)
    elif foundover == True:
        print '\nChecked %g files, found overlaps, see output file\
%s for details.\n' %(len(files),oname)   
    elif foundgap == True:
        print '\nChecked %g files, found gaps, see output file %s for\
details.\n' %(len(files),oname)   
    else:
        print '\nChecked %g files.\n' %len(files)      
            
