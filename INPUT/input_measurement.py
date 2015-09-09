# Input parameters for rapid source distribution estimate
# ============================================================================
# Input paths
# ============================================================================
# Measurement input or input file for great circle plotting and binning.
input='DATA/correlations/ch_jan/CH.FUORN*ccc.*.SAC'
#'../time_dependent_sources/120_months_of_hum/data/hum_05_ccc_feb.006.008.ccc.hann.msr2.txt'
out_basename='dummy'
# ============================================================================
# Input relating to the measurement
# ============================================================================
verbose=True
doplot=True

g_speed=2800.
# If you have group speed measurements in a file, provide its path here. Otherwise set to None
g_speed_msr = None
f_centr=0.15
q=120.
hw=40

# Prefilter (format: (freq_low,freq_high,corners) or None)
prefilter=(0.1,0.2,3)

# window type
window='hann'

# Are windows allowed to overlap? If yes, set to True
win_overlap = True

# Separation of noise window from signal window: 1 --> 1*halfwidth, 2 --> 2*halfwidth etc.
sepsignoise = 2.  

# Phase weighted stack? Set to value > 0 to include
ps_nu=0

# ============================================================================
# Inputs relating to plotting and binning
# ============================================================================

# length of great circle segments in km
segper=100.

# Maximum number of such segments to consider. If set to 'None' they go all
# the way to the antipodal point 
num_max = None

# Minimum signal to noise energy ratio
snr_min=1.

# Sign convention: Obspy ccc 1, obspy pcc -1
signconv=1.



# ============================================================================
# Inputs relating to binning
# ============================================================================ 


# Binning: Min lon
lonmin=-180
# Binning: Max lon
lonmax=180
# Binning: Min lat
latmin=-90
# Binning: Max lat
latmax=90
# Binning: d_deg lon
ddeg_lon=2.
# Binning: d_deg lat
ddeg_lat=2900.
# Bin weighting by approximate size of square
bin_weight=False


