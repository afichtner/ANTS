# A script to process ambient vibration records
import obspy as obs
import numpy as np
import os
import sys
sys.path.append("./general")
sys.path.append("./download")
sys.path.append("./processing")
import matplotlib.pyplot as plt
import avp_pre as pre
from read_xml import read_xml
from download_data import download_data
from downsample import downsample
from fourier_tools import fourier_spectrum


#Input file handling--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read input
inp1 = read_xml('avp_input.xml')
inp1 = inp1[1]
#Also read the specifications for data:
dat1=read_xml('avp_data.xml')
dat1=dat1[1]


#Load data if requested. If no data is downloaded, the latest channel list will be used as it is.----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if inp1['datadownload']=='1':
    from renamer import rename_seismic_data
    download_data('avp_data.xml')
    downloadloc=dat1['downloadloc']
    #a list of the downloaded data has to be obtained
    cmd='./general/get_datalist.sh ' +downloadloc+'/'+dat1['starttime']+"*"
    os.system(cmd)
    
    targetdir=dat1['targetdir']
    datadir=downloadloc + '/' +dat1['starttime']+'_'+dat1['endtime']+'*'
    rename_seismic_data(datadir, targetdir)

# the read command of obspy can be used to read common file types like MiniSEED, sac and so on. It contains a number of traces. 
#Each trace consists of stats object containing meta information and the data itself which is stored in a numpy ndarray.

#Processing routine------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Open file after file for processing
filelist = open('av_channellist.txt', 'r')
filelist = filelist.readlines()


for filename in filelist:
    endl=os.linesep
    filename=filename.strip(endl)
    
    datastream=obs.read(filename)
    data=datastream.copy()
     
    #Have to get sampling rate for different purposes 
    Fs=data[0].stats['sampling_rate']
    
    #Step 0: Necessary preliminaries
    #0 a: Downsampling and lowpassfilter--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if dat1['rundown']=='1':
        Fsnew=int(dat1['downs'])
        data=downsample(data, filename, Fs, Fsnew, data[0].stats._format)
        
    #0 b: Demean and detrend
    data.detrend('demean')
    data.detrend('linear')
    #print data[0].data.mean()
    data.plot()
    
    #Step 1: Preprocessing---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #1 a: Time domain normalization
    if inp1['onebit']=='1':
        data[0].data = pre.onebit(data[0].data)
        
    if inp1['clip']=='1':
        data[0].data = pre.clip(data[0].data)
        
    if inp1['ram_normal']=='1':
        windowlength=int(inp1['winlength'])
        data[0].data = pre.ram_normal(data[0].data, windowlength)
        
    if inp1['waterlevel']=='1':
        level=int(inp1['level'])
        data[0].data = pre.waterlevel(data[0].data,level)
        
    #1 b: Whitening
    if inp1['whitening']=='1':
        #Testversion: plot the PSD before and after whitening. So as to see whether it behaves the way it should.
        from psd_show import psdplot
        psdplot(data[0].data, Fs)
        
        whitesmooth=int(inp1['whitesmooth'])
        data[0].data = np.real(pre.whiten(data[0].data,whitesmooth, Fs))
        
        #(freq, spec)=fourier_spectrum(data[0].data, Fs)
        #plt.plot(freq[1000:5000], abs(spec[1000:5000]))
        #plt.xlabel('Amplitude [?]')    
        #plt.xlabel('Frequency [Hz]')
        #plt.show()
        psdplot(data[0].data, Fs)
    
    #Write the data to a file. In the same directory as unprocessed file.
    if inp1['saveprep']=='1':
        filename=filename+'_prep'
        data.write(filename, data[0].stats._format)
