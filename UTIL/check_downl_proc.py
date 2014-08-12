import sys
import collections
import os
from obspy import UTCDateTime

if __name__=='__main__':
    import check_downl_proc as cdp
    
    cdp.check_downl_proc(indir=str(sys.argv[1]),outfile=str(sys.argv[2]))
    
def check_downl_proc(indir=None,infile=None,outfile=None):
    
    if indir != None:
        
        filelist=os.listdir(indir)
        
    elif infile != None:
        
        filelist = open(infile,'r')
        filelist = filelist.read().split('\n')
    else:
        msg = 'Must provide either file or directory.'
        raise Exception(msg)
        
    data = dict()
    
    if filelist is not None:
        
        for file in filelist:
            
            if file == '': continue
            if file.split('.')[-1] not in ['MSEED','mseed','SAC','sac']:
                continue
            
            inf = file.split('.')
            
            id = inf[0] + '.' + inf[1] + '.' + \
                inf[2]+ '.' + inf[3]
            
            
            if id not in data:
                data.update({id: [inf[4]+inf[5],inf[9]+inf[10]]})
                
            else:
                if inf[4]+inf[5] < data[id][0]:
                    data[id][0] = inf[4]+inf[5]
                    
                if inf[9]+inf[10] > data[id][1]:
                    data[id][1] = inf[9]+inf[10]
                    
        data = collections.OrderedDict(sorted(data.items()))
        
    if outfile ==  None:    
        return data
    else:
        ofid=open(outfile,'w')
        for entry in data:
            ofid.write("%16s  |  %8s  |  %8s\n" %(entry, data[entry][0], data[entry][1]))
                
        