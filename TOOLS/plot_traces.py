#Plot several preprocessed traces in a directory for convenient short check
import obspy as obs
import os

def plot_traces(indir):
    
    files=os.listdir(indir)
    stream=obs.Stream()
  
    
    
    for file in files:
        
        
        try:
            str=obs.read(os.path.join(indir+file))
            
            if len(stream)==0: 
                stream=str
                
            elif (str[0].stats.station==stream[0].stats.station):
                stream.append(str[0])
                
            else:
                stream.plot()
                stream=obs.Stream()
            
        except (IOError, TypeError):
            print 'Cannot read file.'
            continue
            
   
