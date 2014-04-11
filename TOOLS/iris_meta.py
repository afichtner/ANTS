# A script for producing iris-like metadata
import os
from glob import glob
from TOOLS.read_xml import read_xml, find_coord


def write_stationlst(indir, xmldir, outdir,corrname):
    
    files=os.listdir(indir)
    lst=list()
    outfid=open(outdir+corrname+'.station.lst', 'w')
    
    
    
    for file in files:
        try:
            (net, sta)=(file.split('.')[0], file.split('.')[1])
            if ((net, sta)) not in lst:
                lst.append((net, sta))
        except IndexError: 
            continue
                
    
     
    for entry in lst:

        xmlfile=glob(xmldir+'/'+entry[0]+'.'+entry[1]+'*.xml')
        
        if len(xmlfile)==1:
           
           (staname,lat,lon)=find_coord(xmlfile[0])
           
           outfid.write(entry[0].ljust(5)+entry[1].ljust(8)+str(lat).ljust(11)+str(lon).ljust(10)+'\n')
        
