import obspy

def rename_seismic_data(data,targetdir,processed,verbose):

    """
    rename file to network.station.location.channel_starttime_endtime_samplingrate_format

    input:
        data: ObsPy data stream
        targetdir: directory where files are written
        processed: True for processed data, False for original data
        verbose: print screen output (True)
    """

    #- read statistics
    starttime=data.stats.starttime
    endtime=data.stats.endtime
    samplerate=data.stats.sampling_rate
    network=data.stats.network
    station=data.stats.station
    location=data.stats.location
    channel=data.stats.channel
    format=data.stats._format
        
    #- convert time objects to string
    t1=starttime.strftime('%Y%m%d%H%M%S')
    t2=endtime.strftime('%Y%m%d%H%M%S')

    #- new filename
    filepathnew=targetdir+'/'+network+'.'+station+'.'+location+'.'+channel+'.' + t1 + '.' +t2 + '.' + str(int(samplerate))+'.'+ format
        
    if processed==True:
        filepathnew=filepathnew+'.prep'

    
    #- write to file
    data.write(filepathnew,format)
    if verbose==True: print '* renamed file: '+filepathnew   
