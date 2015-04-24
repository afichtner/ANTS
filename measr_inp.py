# Input parameters for rapid source distribution estimate
# ============================================================================
# Input paths
# ============================================================================
# Measurement input
input='DATA/examples/*.SAC'
# Binning input file
inp_binning='example_asymmetry.msr2.txt'
# Plotting input file
inp_plotting = 'example_asymmetry.msr2.txt'

out_basename='example_asymmetry'
verbose=True
doplot=True
dohist=False

g_speed=2900.
# If you have group speed measurements in a file, provide its path here. Otherwise set to None
g_speed_msr = None
f_centr=0.15
q=120
hw=40

# Prefilter (format: (freq_low,freq_high,corners) or None)
prefilter=(0.1,0.2,3)

# window type
window='hann'

# Are windows allowed to overlap? If yes, set to True
win_overlap = False

# Separation of noise window from signal window: 1 --> 1*halfwidth, 2 --> 2*halfwidth etc.
sepsignoise = 2.  

# Phase weighted stack? Set to value > 0 to include
ps_nu=0

# Minimum signal to noise energy ratio
snr_min=10.

# Sign convention: Obspy ccc 1, obspy pcc -1
signconv=1.

# ============================================================================
# Inputs relevant for binning
# ============================================================================ 

# Segment per ? km
segper=100.
# Binning: Min lon
lonmin=-10
# Binning: Max lon
lonmax=30
# Binning: Min lat
latmin=30
# Binning: Max lat
latmax=60
# Binning: d_deg lon
ddeg_lon=0.5
# Binning: d_deg lat
ddeg_lat=0.5


