# setup script to create all the directories that PANOCO needs
import os
import antconfig as cfg

datadir=cfg.datadir
if os.path.exists(datadir)==False:
    os.mkdir(datadir)
    os.mkdir(datadir+'resp/')
    os.mkdir(datadir+'raw/')
    os.mkdir(datadir+'processed/')
    os.mkdir(datadir+'correlations/')
    os.mkdir(datadir+'stationxml/')
    
    os.mkdir(datadir+'processed/out/')
    os.mkdir(datadir+'processed/xmlinput/')
    os.mkdir(datadir+'correlations/interm/')
    os.mkdir(datadir+'correlations/xmlinput/')
    
