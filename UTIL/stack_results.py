# A script that

# - obtain a list of the correlations to be stacked
# - reading is HEAVY here
# - Try to do it serial, see how long it takes
# - go through each item of the list, load data, stack, write.
# - phase weighted stack?


import UTIL.read_sac_correlation as readsac
import os
from par_st import parlistpairs



def stack_results(inp_ids, stackname, corrtype, \
                   outdir='DATA/stacks/',verbose=False):
    
    fid = open(inp_ids,'r')
    pats = fid.read().split('\n')
    
    
    for entry in pats:
        if entry == ' ': continue
        cstream = readsac.read_corr_sac(entry, corrtype, verbose)
        stk = cstream.stack(location1='00',location2='00')
            
        if stk.n_stack > 0:
            readsac.write_corr_sac(stk,corrtype,stackname,outdir,verbose)
        
        
def create_idlist(corr_ids,corrtype,prefix,cha1='LHZ',cha2='LHZ',\
                    autocorr=False,outfile='dummyids.txt'):
    
    fid_in = open(corr_ids, 'r')
    fid_out = open(outfile,'w')
   
    stations = fid_in.read().split('\n')
    fid_in.close()
    
    ids = list()
    stalist = list()
    
    for sta in stations:
        #- Sort out empty lines
        if sta=='': continue
        #- Sort out doubles
        if sta not in stalist:
            stalist.append(sta)
        
            
    for i in range(len(stalist)):
        for j in range(i+1):
            sta1 = stalist[i].split('.')[1]
            sta2 = stalist[j].split('.')[1]
            
            if sta1 == sta2 and autocorr == False:
                continue
            
            if stalist[i] <= stalist[j]:
                fid_out.write(prefix+'/*.'+\
                            sta1+'.*.'+cha1+'.*.'+sta2+'.*.'+\
                                cha2+'.'+corrtype+'*.SAC\n')
            else:
                fid_out.write(prefix+'/*.'+\
                            sta2+'.*.'+cha2+'.*.'+sta1+'.*.'+\
                                cha1+'.'+corrtype+'*.SAC\n')
                 
    
    fid_out.close()