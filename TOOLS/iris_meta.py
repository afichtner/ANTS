# A script for producing iris-like metadata
import os
from glob import glob
from TOOLS.read_xml import read_xml


def write_stationlst(indir, xmldir, outdir):
    
    files=os.listdir(indir)
    lst=list()
    outfid=open(outdir+'/station.lst', 'w')
    
    
    
    for file in files:
        
        (net, sta)=(file.split('.')[0], file.split('.')[1])
        if ((net, sta)) not in lst:
            lst.append((net, sta))
            
    
     
    for entry in lst:

        xmlfile=glob(xmldir+'/'+entry[0]+'.'+entry[1]+'*.xml')
        
        if len(xmlfile)==1:
           
           inf=read_xml(xmlfile[0])[1]
           lat=float(inf['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lat'])
           lon=float(inf['{http://www.data.scec.org/xml/station/}Network']['{http://www.data.scec.org/xml/station/}Station']['{http://www.data.scec.org/xml/station/}StationEpoch']['{http://www.data.scec.org/xml/station/}Lon'])
           
           outfid.write(entry[0].ljust(5)+entry[1].ljust(8)+str(lat).ljust(11)+str(lon).ljust(10)+'\n')
        
