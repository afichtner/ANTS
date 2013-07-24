import numpy as np
from obspy.core import read
from obspy.core import stream
from obspy.core import  UTCDateTime
import os

import matplotlib.pyplot as plt
import TOOLS.read_xml as rxml
import TOOLS.correlations as corr

#==================================================================================================
# Go through day after day and stack
#==================================================================================================
def stack(stack_input):
    #Read the input
    inp1=rxml.read_xml(stack_input)
    inp1=inp1[1]

    indir=inp1['directories']['indir']
    outdir=inp1['directories']['outdir']
    verbose=bool(int(inp1['verbose']))
    plotting=bool(int(inp1['plotting']))
    
    
    startday=inp1['timethings']['startdate']
    endday=UTCDateTime(inp1['timethings']['enddate'])
    print 'endday = ',  endday
    
    #Force maximum lag to be integer. In seconds:
    maxlag=int(inp1['correlations']['max_lag'])
    
    #Use directory of the first day in order to find the channel combinations
    if os.path.exists(indir+'/'+startday+'/.DS_Store'):
        os.remove(indir+'/'+startday+'/.DS_Store')
    
    listch=os.listdir(indir+'/'+startday)
    #Find out what are the relevant channel combinations. This returns a list of tuples with identifiers that are to be correlated (e. g. ('G.ECH.00.BHE','G.CAN.00.BHE'))
    corr_ch=matchchannels(listch)
    
    #-Loop over station pairs==========================================================================
    for chpair in corr_ch:
        if verbose:
            print 'Stacking correlations for:'
            print chpair[0]
            print chpair[1]
            
        #open those for the first day, to get the sampling rate
        dat1=read(indir+'/'+startday+'/'+chpair[0]+'*')[0]
        dat2=read(indir+'/'+startday+'/'+chpair[1]+'*')[0]
        #Prepare stacking loop
        #xcorrstack=np.zeros((2*maxlag*dat1.stats.sampling_rate+1, 1))
        
        xcorrstack=corr.xcorrelation_td(dat1, dat2, maxlag)
        n=1
        day=UTCDateTime(startday)
        
        #Loop over days, update stack
        while day<=endday-86400:
            #Obtain correlation
            #xcorrstack+=corr.xcorrelation_td(dat1, dat2, maxlag)
            #n+=1
            
            #Step day
            day+=86400
            crntday=day.strftime('%Y%m%d')
            #get new data
            dat1=read(indir+'/'+crntday+'/'+chpair[0]+'*')[0]
            dat2=read(indir+'/'+crntday+'/'+chpair[1]+'*')[0]
            
            if len(dat1) < 86000 or len(dat2)<86000: continue
            
            xcorrstack+=corr.xcorrelation_td(dat1, dat2, maxlag)
            n+=1
            
            
            if plotting:
                x=np.linspace(-maxlag*dat1.stats.sampling_rate,maxlag*dat1.stats.sampling_rate, len(xcorrstack))
                plt.plot(x, xcorrstack, linewidth=1.8)
                plt.show()
      
        
        #Linear stack is easy, what about phase weighted stack? Two-step procedure, the phase stack needs to be determined first. 
        xcorrstack/=n
        #write the stacked functions---> to an ascii file? to an mseed? What would be the headers then?
        
    
    




def matchchannels(channels):
    
    ccpairs=[]
    
    for i in range(len(channels)):
        for j in range(len(channels)):
            if i<j: continue
            inf1=channels[i].split('.')
            inf2=channels[j].split('.')
            #station
            sta1=inf1[1]
            sta2=inf2[1]
            #channel
            cha1=inf1[3]
            cha2=inf2[3]
            #sampling rate
            Fs1=inf1[6]
            Fs2=inf2[6]
            
            chcomb=()
            if  sta1!=sta2 and cha1==cha2 and Fs1==Fs2:
                name1=inf1[0]+'.'+inf1[1]+'.'+inf1[2]+'.'+inf1[3]
                name2=inf2[0]+'.'+inf2[1]+'.'+inf2[2]+'.'+inf2[3]
                chcomb=(name1,name2)
                ccpairs.append(chcomb)
    
    return ccpairs
            
    #what would we like here? I guess a python list of tuples 
    # I d like to have an extra function that takes care of this
    
    
    #Allocate variables: each channel combination should get a structure where to put the correlation array and some meta-information (sample rate; correlated channels; number of correlated days; anything else?)
    #Where to put this metadata in the final file? In the header? In the 
    



#==================================================================================================
# Run the selected type of correlation on a number of files in a specific directory
#==================================================================================================
def xcorr(indir, lag   , verbose=False):
    
    #Create filelist
    cmd='ls '+indir+' > filelist.txt'
    os.system(cmd)
    #Open filelist
    filelist = open('filelist.txt', 'r')
    filelist=filelist.read().split('\n')[:-1]
    
    #Cross-correlation loop over files
    for i in range(len(filelist)):
       for j in range(len(filelist)):
            if i<j: continue
            
            filename1=indir+'/'+filelist[i]
            filename2=indir+'/'+filelist[j]
            dat1=read(filename1)[0]
            dat2=read(filename2)[0]        
           
            if  dat1.stats.station!=dat2.stats.station and dat1.stats.channel==dat2.stats.channel:
                
                if dat1.stats.sampling_rate!=dat2.stats.sampling_rate:
                    print 'Unequal sampling rates. Cannot crosscorrelate.'
                    continue
        
                if verbose:
                    print '---------------------------------------------------------------------'
                    print 'Crosscorrelation of...'
                    print filename1 
                    print filename2
                
                #- Classical cross-correlation time domain============================================================================
                if inp1['correlations']['cc_td']=='1':
                    if verbose:
                        print '* Time-domain classical cross-correlation'
                    cc_td=corr.xcorrelation_td(dat1, dat2, max_lag)
                    taxis=np.linspace(-max_lag, max_lag, len(cc_td), True)
                    plt.plot(taxis, cc_td, linewidth=2.0)
                
                
                #- Classical cross-correlation frequency domain=======================================================================
                if inp1['correlations']['cc_fd']=='1':
                    if verbose:
                        print '* Frequency-domain classical cross-correlation'
                    cc_fd=corr.xcorrelation_fd(dat1, dat2)
                    taxis=np.linspace(0,len(cc_fd)/dat1.stats.sampling_rate, len(cc_fd), True)
                    plt.plot(taxis, cc_fd, linewidth=2.0)
                    
                #- Classical cross-correlation frequency domain=======================================================================
                if inp1['correlations']['pcc']['doit']=='1':
                    if verbose:
                        print '* Phase cross-correlation'
                    nu=int(inp1['correlations']['pcc']['nu'])
                    pcc=corr.phase_xcorrelation(dat1, dat2, max_lag, nu)
                    taxis=np.linspace(-max_lag, max_lag, len(pcc), True)
                    plt.plot(taxis, pcc, linewidth=2.0)
                plt.xlabel('Lage time seconds')
                plt.ylabel('Cross-Correlation')
                plt.legend(('Time domain CC', 'Freq. domain CC', 'Phase CC'))
                plt.show()

    #Clean up
    os.system('rm filelist.txt')
