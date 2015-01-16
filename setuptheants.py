# setup script to create all the directories that PANOCO needs
import os

import antconfig as cfg

datadir=cfg.datadir
if os.path.exists(datadir)==False:
    os.mkdir(datadir)
    os.mkdir(datadir+'resp/')
    os.mkdir(datadir+'raw/')
    os.mkdir(datadir+'raw/latest/')
    os.mkdir(datadir+'processed/')
    os.mkdir(datadir+'correlations/')
    os.mkdir(datadir+'stationxml/')
    
    os.mkdir(datadir+'processed/out/')
    os.mkdir(datadir+'processed/xmlinput/')
    os.mkdir(datadir+'correlations/interm/')
    os.mkdir(datadir+'correlations/out/')
    os.mkdir(datadir+'correlations/xmlinput/')

testdir=cfg.testdir
if os.path.exists(testdir)==False:
    os.mkdir(testdir)
    
resdir=cfg.resdir
if os.path.exists(resdir)==False:
    os.mkdir(resdir)
    os.mkdir(resdir+'maps/')
    os.mkdir(resdir+'correlations/')
 