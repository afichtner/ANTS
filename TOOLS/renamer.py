import obspy

def rename_seismic_data(data,targetdir,processed,verbose, ofid=None):

    """
    rename file to network.station.location.channel_starttime_endtime_samplingrate_format

    input:
        data: ObsPy data stream
        targetdir: directory where files are written
        processed: True for processed data, False for original data
        verbose: print screen output (True)
        ofid: file id of an open file that output can be written to instead of stdout
    """
    
    #- read statistics
    if isinstance(data, obspy.Stream) and len(data)>1:
        starttime=data[0].stats.starttime
        endtime=data[len(data)-1].stats.endtime
        stats=data[0].stats
    elif isinstance(data, obspy.Stream) and len(data)==1:
        starttime=data[0].stats.starttime
        endtime=data[0].stats.endtime
        stats=data[0].stats
    elif isinstance(data, obspy.Trace):
        starttime=data.stats.starttime
        endtime=data.stats.endtime
        stats=data.stats
    else:
        if ofid==None:
            print 'Nothing saved: Object is not a trace or stream.'
        else:
            ofid.write('Nothing saved: Object is not a trace or stream.')
        return
        
        
    samplerate=stats.sampling_rate
    network=stats.network
    station=stats.station
    location=stats.location
    channel=stats.channel
    format=stats._format
        
    #- convert time objects to string
    t1=starttime.strftime('%Y%m%d%H%M%S')
    t2=endtime.strftime('%Y%m%d%H%M%S')

    if processed==True:
        filepathnew=targetdir+'/'+network+'.'+station+'.'+location+'.'+channel+'.' + t1 + '.' +t2+'.proc.'+ format
    else:
        filepathnew=targetdir+'/'+network+'.'+station+'.'+location+'.'+channel+'.' + t1 + '.' +t2 +'.'+ format
    
    #- write to file
    data.write(filepathnew,format)
    if verbose==True:
        if ofid==None:
            print '* renamed file: '+filepathnew 
        else:
            ofid.write('* renamed file: '+filepathnew+'\n')
