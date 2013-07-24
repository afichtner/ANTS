from obspy import UTCDateTime
import subprocess
import os

def sortbyday(indir, firstday, lastday, secondstolerance=10.0, verbose=False):

    firstday=UTCDateTime(firstday)
    lastday=UTCDateTime(lastday)
    
    datetime=firstday
    while datetime<=lastday:
        
        if verbose:
            print 'Searching data for day:'
            print datetime
       
        cmd='find ' + indir+' -maxdepth 1 -type f -not -name .DS_Store > sortfilelist.txt'
        os.system(cmd)
    
        
        filelist=open('sortfilelist.txt', 'r')
        filelist=filelist.read().split('\n')
        if len(filelist)==0: return()
        
        for file in filelist:
            file=file.split('/')[-1:][0]
            if file=='': continue
    
            
            startdate=UTCDateTime(file.split('.')[4])
           
            
            if abs(startdate-datetime)<secondstolerance:
                dirname=indir+'/'+datetime.strftime('%Y%m%d')
                if os.path.isdir(dirname)==False:
                    cmd=('mkdir '+dirname)
                    os.system(cmd)
                cmd=('mv '+indir+'/'+file+' '+dirname)
                os.system(cmd)
        
        datetime+=3600*24
            
    os.system('rm sortfilelist.txt')
    
