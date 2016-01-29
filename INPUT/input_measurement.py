# Input parameters for rapid source distribution estimate
# ============================================================================
# Input paths
# ============================================================================
# Measurement input or input file for great circle plotting and binning.
#
input='bhutan_pcc_jan.msr2.txt'
#'/Volumes/cowpox/DATA/correlations/bhutan_pcc_jan/rank*/*.SAC'
#
#'/Volumes/cowpox/DATA/correlations/bhutan_noquakes_jan_1/*.SAC'
#'bhutan_noquakes_jul_1_mindist.msr2.txt'
#
#'bhutan_noquakes_jan_1.msr2.txt'
#
#'bhutan_noquakes_jan_1.msr2.txt'
#
#'bhutan_noquakes_jan_1.msr2.txt'
#'/Volumes/cowpox/DATA/correlations/bhutan_noquakes_jan_1/*.SAC'
#

#'RES/hum/hum_months/hum_jul.006.008.ccc.hann.msr2.txt'
#'
#input='hum_2014.msr2.txt'
#../Results/plots_sourcemapping/hum/measurement_files/hum_jan.006.008.ccc.hann.msr2.txt'
#input='/Volumes/cowpox/insta_iso/*.SAC'
#input='/Volumes/cowpox/insta_homogeneous_8degrees/*.SAC'
out_basename='bhutan_pcc_jan'
# ============================================================================
# Input relating to the measurement and bin-plotting
# ============================================================================
verbose=True
doplot=False

g_speed=2900.
# If you have group speed measurements in a file, provide its path here. Otherwise set to None
g_speed_msr = None
f_centr=0.067
q=120.
hw=30

# Prefilter (format: (freq_low,freq_high,corners) or None)
prefilter=(0.05,0.1,3)

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
segper=10.

# Maximum number of such segments to consider. If set to 'None' they go all
# the way to the antipodal point 
num_max = 1000

# Minimum signal to noise energy ratio
snr_min=5.0

# Sign convention: Obspy ccc 1, obspy pcc -1
signconv=-1.

# Plot type: Specify 'py' if you'd like to have a matplotlib-friendly output
plot_type='py'#'py'

# ============================================================================
# Inputs relating to binning
# ============================================================================ 

# Binning: Min lon
lonmin=60.
# Binning: Max lon
lonmax=115.
# Binning: Min lat
latmin=5.
# Binning: Max lat
latmax=45.
# Binning: d_deg lon
ddeg_lon=0.5
# Binning: d_deg lat
ddeg_lat=0.5
# Bin weighting by approximate size of square
bin_weight=True


