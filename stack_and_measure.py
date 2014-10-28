import measure_asymmetry as ma
import UTIL.stack_results as sr
import os

sr.create_idlist('INPUT/ID_LISTS/ids_medi_corr.txt','pcc',
                    prefix='DATA/correlations/medi_z',cha1='LHZ',cha2='LHZ',\
                    autocorr=False,outfile='dummyids.txt')

sr.stack_results('dummyids.txt','test1', 'pcc', \
                    outdir='DATA/stacks/',verbose=True)
                    
#os.system('rm dummyids.txt')

#ma.meas_asym(input='DATA/correlations/hum_selection/stacks/*.jul.*',\
#            filename='hum_select.msr.hann.txt',g_speed=3590.,w1=325,\
 #           w2=275,window='hann',verbose=True)