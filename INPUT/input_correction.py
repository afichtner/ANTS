
#*******************************************************************************
# Name tag. Must be unique to this dataset unless update is set to True
# ******************************************************************************
prepname = 'noisy'
#*******************************************************************************
# Output file format. Set to 'mseed','sac' or 'same'. 'same' will save the ouput
# files in the same format as the input files were.
# ******************************************************************************
# UPDATE Needed: Implement format change if asked for!!
# outform = 'mseed'
#*******************************************************************************
# Add an optional description in this comment line
# Apply a trigger on the swiss data -- for comparison if something changes.
# ...
# ******************************************************************************
verbose = True
outfile = True
# If check is set to true, only the first up to four files are preprocessed and images of several stages saved to processed/out in you data directory
check = False
debugfile = 'test_noisy.txt'
# If update is set to true, script first controls which files have already been processed, and skips these.
update = False

#*******************************************************************************
# Input paths
# ******************************************************************************
indirs = ['<Your Directory>/DATA/raw/latest/rank0/']


#*******************************************************************************
# Quality/gaps
# ******************************************************************************

min_length_in_sec = 0.
maxgaplen = 3600.

#*******************************************************************************
# Event excluder
# ******************************************************************************
exclude_events = False
exclude_windows = [900,3600,7200]
exclude_std = 2.
exclude_n = 6
exclude_freq = 0.001
exclude_level = 1.5
#*******************************************************************************
# Preprocessing
# ******************************************************************************
split_do = False
length_in_sec = 131072
trim = True
detrend = True
demean = True
taper_do = True
taper_width = 0.05
# downsampling:
#  The original sampling rates are needed, because traces with a wrong sampling rate will to be sorted out. If the new sampling rate is different from the old one, data will be downsampled. A lowpass filter is hardcoded into downsampling! It has a corner frequency of TWO FIFTHS (0.4) of the new sampling frequency. When a large amont of downsampling is performed, do it in two steps! e. g. instead of by factor 20 decimate by 5 then 4. IMPORTANT: Since correlation windows are cut by start and end time in seconds, don't decimate below 1 Hz in the preprocessing stage!!

Fs_old = [1.0]
Fs_new = [1.0,0.5]

#*******************************************************************************
# Instrument correction
# ******************************************************************************
remove_response = True
unit  = 'VEL' #'DIS','ACC'
respdir = '<Your Directory>/DATA/resp/'
# This is a prefilter (bandpass) that works like in SAC and restricts the data to the reliable passband of the instrument.
freqs = [0.003, 0.007, 0.1, 0.2]
waterlevel = 0.0
