from __future__ import print_function
import obspy

def rename_seismic_data(data,prepname,verbose, ofid):

    """
    rename file to network.station.location.channel_starttime_endtime_samplingrate_format

    input:
        data: ObsPy data stream
        targetdir: directory where files are written
        prepname: name of preprocessing run
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
        print('** Nothing saved: Object is not a trace or stream.',file=ofid)
        
        
    samplerate=stats.sampling_rate
    network=stats.network
    station=stats.station
    location=stats.location
    channel=stats.channel
    format=stats._format
        
    #- convert time objects to string
    t1=starttime.strftime('%Y.%j.%H.%M.%S')
    yr=str(t1[0:4])
    t2=endtime.strftime('%Y.%j.%H.%M.%S')

    filepathnew='DATA/'+yr+'/'+network+'.'+station+'.'+location+'.'+channel+'.' + t1 + '.' +t2+'.'+prepname+'.'+format
    
    #- write to file
    data.write(filepathnew,format)
    if verbose==True:
        print('* renamed file: '+filepathnew,file=ofid)
        