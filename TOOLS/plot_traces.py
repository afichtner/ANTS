#Plot several preprocessed traces in a directory for convenient short check
import obspy as obs
import os


def plot_traces(indir, num_traces=None):
    
    files=os.listdir(indir)
    stream=obs.Stream()
    
    cnt=1
    
    for file in files:
        
        if type(num_traces)==int:
            if cnt==num_traces:
                try:
                    stream.plot()
                except ValueError:
                    print 'WTF?'
                stream.clear()
                cnt=1
                continue
   
        
        try:
            str=obs.read(os.path.join(indir+file))
            print str[0]
            
            if len(stream)==0: 
                stream=str
                
            elif (str[0].stats.station==stream[0].stats.station):
                stream.append(str[0])
                
            else:
                try:
                    stream.plot()
                except ValueError:
                    print 'WTF?'
                stream=str
            
            cnt+=1
            
        except (IOError, TypeError):
            print 'Cannot read file.'
            try:
                stream.plot()
            except ValueError, IndexError:
                print 'WTF?'
            stream.clear()
            
            
            continue
            
   
