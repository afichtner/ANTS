#- Script to quickly get a table of successfully processed data
from glob import glob
from obspy import UTCDateTime
import collections
import time
import sys

if __name__=='__main__':
    import listproc as lp
    lp.listproc(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    
def listproc(indir,oname,procname,network='*',channel='*',location='*'):
    """
    indir: Directory where processed data is located (string)
    oname: The information will be provided as ascii file in this location (string)
    procname: processing run name to check (string)
    """
    ofid=open(oname,'w')
    ofid.write('Nr  ID          t1                   t2                  gap (hours) \n' )
    ofid.write('====================================================================\n')
    files=glob(indir+'/*.'+procname+'.*')
    files.sort()
    toc=dict()
    
    for file in files:
        print file
        file=file.split("/")[-1]
        file=file.split('.')
        
        if network!='*' and file[0]!=network:
            continue
        if channel!='*' and file[2]!=channel:
            continue
        if location!='*' and file[0]!=location:
            continue
        
        id=file[0]+'.'+file[1]+'.'+file[2]+'.'+file[3]
        t1_new=UTCDateTime(file[4]+','+file[5]+','+file[6]+','+file[7]+','+file[8])
        t2_new=UTCDateTime(file[9]+','+file[10]+','+file[11]+','+file[12]+','+file[13])
        
        if id in toc:
            t1_old=toc[id][0]
            t2_old=toc[id][1]
            gap=toc[id][2]
            print gap
            
            if t1_new>t2_old:
                if t1_new-t2_old>1:
                    gap+=t1_new-t2_old-1
                print gap
                toc[id]=(t1_old,t2_new,gap) 
                
            elif t2_new<t1_old:
                if t1_old-t2_new>1:
                    gap+=t1_old-t2_new-1
                print gap
                toc[id]=(t1_new,t2_old,gap)
                
        else:
            toc.update({id:(t1_new,t2_new,0)})
          
          
    toc=collections.OrderedDict(sorted(toc.items()))   
    k=1
    
    for item in toc:
        ofid.write(str(k)+'  '+item+'  '+toc[item][0].strftime('%Y.%m.%d.%H.%M.%S')+'  '+toc[item][1].strftime('%Y.%m.%d.%H.%M.%S')+'  '+str(toc[item][2]/3600)+'  ')
        ofid.write("\n")
        k+=1
    ofid.close()
            
                
            
            
