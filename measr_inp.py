# Input parameters for rapid source distribution estimate
# ============================================================================
# Input paths
# ============================================================================
# Measurement input
input='DATA/stacks/G.CAN*.SAC'
# Input file for both binning and great circle path plotting
inp_binning_plotting='example_asymmetry.msr2.txt'


# ============================================================================
# Input relating to the measurement
# ============================================================================
out_basename='test'
verbose=True
doplot=True
dohist=False

g_speed=3700.
# If you have group speed measurements in a file, provide its path here. Otherwise set to None
g_speed_msr = None
f_centr=0.007
q=138
hw=200

# Prefilter (format: (freq_low,freq_high,corners) or None)
prefilter=(0.004,0.01,3)

# window type
window='hann'

# Are windows allowed to overlap? If yes, set to True
win_overlap = True

# Separation of noise window from signal window: 1 --> 1*halfwidth, 2 --> 2*halfwidth etc.
sepsignoise = 4.  

# Phase weighted stack? Set to value > 0 to include
ps_nu=0

# ============================================================================
# Inputs relating to plotting and binning
# ============================================================================

# length of great circle segments in km
segper=100

# Maximum number of such segments to consider. If set to 'None' they go all
# the way to the antipodal point 
num_max = 50

# Minimum signal to noise energy ratio
snr_min=5.

# Sign convention: Obspy ccc 1, obspy pcc -1
signconv=1.

# Segment per ? km
segper=100.

# ============================================================================
# Inputs relating to binning
# ============================================================================ 


# Binning: Min lon
lonmin=-180
# Binning: Max lon
lonmax=180
# Binning: Min lat
latmin=-80
# Binning: Max lat
latmax=80
# Binning: d_deg lon
ddeg_lon=2.
# Binning: d_deg lat
ddeg_lat=2.
# Bin weighting by approximate size of square
bin_weight=False


