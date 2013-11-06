#Plot several preprocessed traces in a directory for convenient short check
import obspy as obs
import os


def plot_traces(indir, num_traces=None):
    
    files=os.listdir(indir)
    stream=obs.Stream()
    
    if type(num_traces)==int:
            cnt=1
    
    for file in files:
        
        
        if cnt>num_traces:
            try:
                stream.plot()
            except ValueError:
                print 'WTF?'
            
            stream.clear()
            cnt=1
            #raw_input('Wanna see more? Hit Enter ')
            
        
        try:
            str=obs.read(os.path.join(indir+file))
            
            if len(stream)==0: 
                stream=str
                
            elif (str[0].stats.station==stream[0].stats.station):
                stream.append(str[0])
                
            else:
                try:
                    stream.plot()
                except ValueError:
                    print 'WTF?'
                stream.clear()
            
            cnt+=1
            
        except (IOError, TypeError):
            print 'Cannot read file.'
            continue
            
   
