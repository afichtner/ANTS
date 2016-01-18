# Input parameters for rapid source distribution estimate
# ============================================================================
# Input paths
# ============================================================================
# Measurement input or input file for great circle plotting and binning.
#
input='hum_rays.msr2.txt'
#'RES/hum/hum_months/hum_jul.006.008.ccc.hann.msr2.txt'
#'
#input='hum_2014.msr2.txt'
#../Results/plots_sourcemapping/hum/measurement_files/hum_jan.006.008.ccc.hann.msr2.txt'
#input='/Volumes/cowpox/insta_iso/*.SAC'
#input='/Volumes/cowpox/insta_homogeneous_8degrees/*.SAC'
out_basename='hum_rays'
# ============================================================================
# Input relating to the measurement and bin-plotting
# ============================================================================
verbose=True
doplot=False

g_speed=3700.
# If you have group speed measurements in a file, provide its path here. Otherwise set to None
g_speed_msr = None
f_centr=0.0055
q=120.
hw=400

# Prefilter (format: (freq_low,freq_high,corners) or None)
prefilter=(0.003,0.08,3)

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
segper=20.

# Maximum number of such segments to consider. If set to 'None' they go all
# the way to the antipodal point 
num_max = None

# Minimum signal to noise energy ratio
snr_min=7.5

# Sign convention: Obspy ccc 1, obspy pcc -1
signconv=1.

# Plot type: Specify 'py' if you'd like to have a matplotlib-friendly output
plot_type='py'#'py'

# ============================================================================
# Inputs relating to binning
# ============================================================================ 

# Binning: Min lon
lonmin=-180.
# Binning: Max lon
lonmax=180.
# Binning: Min lat
latmin=-80.
# Binning: Max lat
latmax=80.
# Binning: d_deg lon
ddeg_lon=1.0
# Binning: d_deg lat
ddeg_lat=1.0
# Bin weighting by approximate size of square
bin_weight=True


