# setup script to create all the directories that ANTS needs
import os
import antconfig as cfg

datadir=cfg.datadir
if os.path.exists(datadir)==False:
    os.mkdir(datadir)
if os.path.exists(datadir+'resp/')==False:    
    os.mkdir(datadir+'resp/')
if os.path.exists(datadir+'raw/')==False:
    os.mkdir(datadir+'raw/')
    os.mkdir(datadir+'raw/latest/')
if os.path.exists(datadir+'processed/')==False:
    os.mkdir(datadir+'processed/')
    os.mkdir(datadir+'processed/out/')
    os.mkdir(datadir+'processed/xmlinput/')
    os.mkdir(datadir+'processed/input/')
if os.path.exists(datadir+'correlations/')==False:
    os.mkdir(datadir+'correlations/')
    os.mkdir(datadir+'correlations/out/')
    os.mkdir(datadir+'correlations/input/')
if os.path.exists(datadir+'stationxml/')==False:
    os.mkdir(datadir+'stationxml/')
