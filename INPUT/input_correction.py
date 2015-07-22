prepname = 'noisy'
#*******************************************************************************
# Add an optional description in this comment line
# This is a preprocessing example
# ...
# ******************************************************************************
verbose = True
outfile = True
# If check is set to true, only the first up to four files are preprocessed and images of several stages saved to processed/out in you data directory
check = False
# If update is set to true, script first controls which files have already been processed, and skips these.
update = False

#*******************************************************************************
# Input paths
# ******************************************************************************
indirs = ['DATA/raw/noisy']


#*******************************************************************************
# Quality/gaps
# ******************************************************************************

min_length_in_sec = 3600
maxgaplen = 300

#*******************************************************************************
# Preprocessing
# ******************************************************************************
split_do = True
length_in_sec = 131072
trim = True
detrend = True
demean = True
taper_do = True
taper_width = 0.05
# downsampling:
#  The original sampling rates are needed, because traces with a wrong sampling rate will to be sorted out. If the new sampling rate is different from the old one, data will be downsampled. A lowpass filter is hardcoded into downsampling! It has a corner frequency of TWO FIFTHS (0.4) of the new sampling frequency. When a large amont of downsampling is performed, do it in two steps! e. g. instead of by factor 20 decimate by 5 then 4. IMPORTANT: Since correlation windows are cut by start and end time in seconds, don't decimate below 1 Hz in the preprocessing stage!!

Fs_old = [1.0]
Fs_new = [1.0]

#*******************************************************************************
# Instrument correction
# ******************************************************************************
remove_response = True
unit  = 'VEL' #'DIS','ACC'
respdir = 'DATA/resp/'
# This is a prefilter (bandpass) that works like in SAC and restricts the data to the reliable passband of the instrument.
freqs = [0.003, 0.007 ,0.25, 0.4 ]
waterlevel = 0.0
