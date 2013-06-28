def download_data(inputfile):
    from obspy import UTCDateTime
    import sys
    sys.path.append("./general")
    sys.path.append("./download")
    import os
    from read_xml import read_xml
        
    datainput = read_xml(inputfile)
    dat1 = datainput[1]
    
    #parameters of data request
    t1 = (dat1['starttime'])
    t2 =(dat1['endtime'])
    
    #How to request several networks?
    network = dat1['network']
    channel = dat1['channel']
    location = dat1['location']
    
    downloadloc = dat1['downloadloc']
    
    if dat1['centre']=='iris':
        centre=' --arc N '
    else:
        centre=''
        
    format= ' --' +dat1['format']+' '
    
    
    stafile=dat1['stations']
    if stafile=='*':
        if os.path.isdir(downloadloc):
            if dat1['centre']=='iris':
                cmd='obspyDMT --iris_update '+downloadloc+' --identity '+ network + ".*."+location + '.' + channel + ' --min_date ' + t1 + ' --max_date ' + t2 + format
                print cmd
                os.system(cmd)
            else:
                cmd='obspyDMT --arc_update '+downloadloc+' --identity '+ network + ".*."+location + '.' + channel + ' --min_date ' + t1 + ' --max_date ' + t2 + format
                print cmd
                os.system(cmd)
        else:
            cmd='obspyDMT --continuous --identity '+ network + ".*."+location + '.' + channel + ' --min_date ' + t1 + ' --max_date ' + t2 +centre+ format + ' --datapath '+downloadloc
            print cmd
            os.system(cmd)
        cmd='obspyDMT --ic_all '+downloadloc+'/'+t1+"*"+'/'+"*/"+'BH_RAW'
        print cmd
        os.system(cmd)
 
        
    else:
        fhst=open(stafile, 'r')
        for line in fhst:
            endl=os.linesep
            station=line.strip(endl)
            identity =  network + '.'+station+'.'+location + '.' + channel 
            if os.path.isdir(downloadloc):
                if dat1['centre']=='iris':
                    cmd='obspyDMT --iris_update '+downloadloc+' --identity '+ identity + ' --min_date ' + t1 + ' --max_date ' + t2 + format
                    os.system(cmd)
                else:
                    cmd='obspyDMT --arc_update '+downloadloc+' --identity '+ identity + ' --min_date ' + t1 + ' --max_date ' + t2 + format 
                    os.system(cmd)
            else:
                cmd='obspyDMT --continuous --identity '+ identity + ' --min_date ' + t1 + ' --max_date ' + t2 +centre+ format + ' --datapath '+downloadloc
                print cmd
                os.system(cmd)
            cmd='obspyDMT --ic_all '+downloadloc+'/'+t1+"*"+'/'+"*/"+'BH_RAW'
            print cmd
            os.system(cmd)    
