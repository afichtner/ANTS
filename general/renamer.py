#sequentially open and rename data files
def rename_seismic_data(indir, targetdir):
    import os
    import sys
    import obspy
    import shutil
    
    sys.path.append("./general")
    
    if os.path.exists(targetdir)==False:
        os.mkdir(targetdir)
    
    cmd='bash get_datalist.sh '+indir
    os.system(cmd)
    
    filelist = open('channellist.txt', 'r')
    filelist = filelist.readlines()
    
    channellist=open('av_channellist.txt', 'w')


    if len(filelist)>0:
        print 'Files copied to: \n-----------------------'
        for filepath in filelist:
            endl=os.linesep
            filepath=filepath.strip(endl)
            
            
            
            #read start and endtime from file
            data=obspy.read(filepath)[0]
            starttime=data.stats.starttime
            endtime=data.stats.endtime
            samplerate=data.stats.sampling_rate
            network=data.stats.network
            format=data.stats._format
            
            #convert time objects to string
            t1=starttime.strftime('%Y%m%d%H%M%S')
            t2=endtime.strftime('%Y%m%d%H%M%S')
           
            #rename files and put in target folder
            filename=filepath.split('/')[-1:]
            #check if the files have been renamed already. dont do it twice.
            if filename[0].split('.')[0]!='dis':
                filepathnew=targetdir +'/'+ filename[0]
                print filepathnew
                channellist.write(filepathnew)
                channellist.write('\n')
            else:
                filename=filename[0].strip('dis.')
                filepathnew=targetdir +'/'+ network + '.'+ filename+ '_' + t1 + '_' +t2 + '_' + str(int(samplerate))+'_'+ format
                print filepathnew
                channellist.write(filepathnew)
                channellist.write('\n')
            shutil.copyfile(filepath, filepathnew)
            
    else:
        print 'No files found.'
            
            
            
            
            
