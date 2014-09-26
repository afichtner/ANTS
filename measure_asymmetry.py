import os
import pylab as P
import numpy as np
import sys

from glob import glob
from obspy import read
from KERNELS.noisemeasurement import Nmeasure

import antconfig as cfg



if __name__=='__main__':
    import measure_asymmetry as ma
    ma.meas_asym(indir=str(sys.argv[1]))

def meas_asym(indir, xmldir='DATA/stationxml/',g_speed=3000.,w1=200.,w2=200.,\
               ps_nu=0,prefilter=None,\
                verbose=False,outdir='RES/asymmetry_measurements/',doplot=True):
    
    """
    A script to plot cross-correlation traces sorted by interstation distance.
    indir: string: path to directory containing results of stack.py
    xmldir: A directory containing the station xml files for all stations (to get the distances). These files must be in xml format and must be called network.station.sta_info.xml
    g_speed: Fastest group speed in the package of frequencies. This will be used to pick the onset of the correlation signal. In km/sec !!
    w1,w2: Cut w1 seconds before and w2 seconds after the predicted 1-D arrival time.
    ps_nu: integer; power to which the weighting by instantaneous phase stack is raised. ps=0 means linear stack.
    prefilter: Tuple of (float,float,integer) which stands for (lower corner, upper corner, order). A butterworth filter to restrict the correlations to a certain frequency band.
    verbose: Boolean, talk or not
    savefig: Boolean, save plot or not
    outdir: string, where to save plot and outfile. default working directory
    """
    
    if os.path.exists(outdir) == False:
        os.system('mkdir '+outdir)
    if indir[-1]!='/':
        indir=indir+'/'
    
    filename = indir.split('/')[-2]
      
    #Initialize output file
    ofid1=open(outdir+filename+'.msr1.txt','w')
    ofid1.write('Input_directory:  '+indir+'\n')
    ofid1.write('Group_speed/evaluated_phase: '+str(g_speed)+'\n')
    ofid1.write('Window_start/end,before/after_arrival:  '+str(w1)+'/'+str(w2)+'\n')
    
    if prefilter is not None:
        ofid1.write('Prefilter: '+ str(prefilter)+ '\n')
     
        
    ofid2=open(outdir+filename+'.msr2.txt','w')
    files = glob(indir+'*.pcc.*.SAC')


    numwins=list()
    
    for file in files:
        
        msr = Nmeasure(file,2,1,3000,200,prefilter=(0.005,0.04,3))
        #msr.plot()
        numwins.append(msr.nw)
        ofid1.write(msr.id + '  %g' %msr.nw)
        ofid2.write('%9.4f %9.4f %9.4f %9.4f %12.2f ' \
                    %msr.geoinf())
        asym = msr.take_measurement()
        print msr.id
        print asym
        ofid1.write('%10.6f\n' %asym)
        ofid2.write('%10.6f\n' %asym)
        
    if doplot == True:
        plot_hist(numwins,filename,outdir)
    
    ofid1.close()
    ofid2.close()
    



def plot_hist(numwins,filename,outdir):
        
        # the histogram of the data with histtype='step'
        n, bins, patches = P.hist(numwins, 20, histtype='bar')
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

        # add a line showing the expected distribution
        y = P.normpdf( bins, np.mean(numwins), 100)
        
        # Labels
        P.xlabel('Nr. of windows used for stack')
        P.ylabel('Nr. correlation stacks')
        P.title(filename.split('/')[:-1])
        
        P.show()
        P.savefig(outdir+ 'hist.png')