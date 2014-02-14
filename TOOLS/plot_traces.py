#Plot several preprocessed traces in a directory for convenient short check
import obspy as obs
from glob import glob


def plot_traces(indir,keyw='*',num_traces=10):
    """
    Little convenience function to obspy-plot all the traces in a directory the name of which contains the key expression keyw.
    indir: Directory where to look for traces
    keyw: string, an expression to filter by
    num_traces: How many traces should be plotted in one panel
    """
    
    files=glob(indir+'/*'+keyw+'*')
    stream=obs.Stream()
    legstr=['...']
    cnt=0
    
    for file in files:
        
        if cnt<num_traces:
           try:
               stre=obs.read(file)
               print stre
               
               if len(stream)==0:
                   stre[0].stats.channel+='.'+str(cnt)+str(stre[0].stats.starttime.strftime('.%Y.%j'))
                   stream=stre
                   cnt+=1
                   
               else: 
                   for tr in stre:
                 
                       tr.stats.channel+='.'+str(cnt)+str(tr.stats.starttime.strftime('.%Y.%j'))
                       tr.stats.starttime=stream[0].stats.starttime
                       stream.append(tr)
                       
                       cnt+=1
                   

           except IOError, TypeError:
               print 'Cannot read file.'
               continue
               
        elif cnt==num_traces:            
           stream.plot(equal_scale=False)
           stream=obs.Stream()
           cnt=0
           
        continue
        
   
