#!/bin/env python
# encoding: utf-8
#
# This is the main script to be launched for reading data, performing burst
# search and plotting results.
#
# The script has been designed to work inside the ipython shell in order to
# give access to interactive data manipulation and plotting.
# 
# To run this file, uncomment only the data filename that you want to process,
# launch "ipython -pylab" (in the same dir as the script) and type:
#
# >>> run -i burst
#
# now in the 'd' variable you have all the data (ph timetag, burst search, etc)
# relative to the data file.
#
# For all the 4-channel plots we call the function plot_mburst_share() 
# passing the fuction name of the desired plot type. For ex. to plot the fret
# histogram type:
#
# >>> plot_mburstm_share(d, fun=hist_fret)
#
# The list of plot functions that can be passed can be found in burst_plot.py
#
# For further details, read the code! It's easy!
# 

from glob import glob
import cPickle as pickle
from path_def import *

try:
    ip = get_ipython()
except:
    ip = _ip

## Loading the fuctions to read the data files
#ip.magic("run -i dataload.py")
from dataload.multi_ch_reader import *

## Loading the lib functions related to the burst search
ip.magic("run -i burstlib.py")

## Loading the plotting functions
ip.magic("run -i burst_plot.py")

swap_D_A = False

## - - - - - - - - - - - - - - - - - - - - - - - - - - - FILE NAME START -.
## Data file path                                                          |
#fname = 'Dualview_8x1_BS_200nm_100X_36uW_0.dat'; BT = 0.09
#fname = 'Dualview_8x1_BS_SSDNA_100pM_30mW_2.dat'; BT = 0.09
#fname = 'Dualview_8x1_Cy3Cy5Dicr_SSDNA_500pM_80mW_3.dat'; BT = 0.09
#fname = 'Dualview_8x1_Cy3Cy5Dicr_SSDNA_100pM_80mW_1.dat'; L= 20; BT = 0.09
#fname = 'DV_8x1_Dicr_DSDNA_tmr_atto_100pM_80mW_2.dat'; BT = 0.09
#fname = '2011-05-04/DV_8x1_Dicr_dsdna_100pM_160mW_12.dat'; BT = 0.09
#fname = '2011-05-04/DV_8x1_Dicr_dsdna_100pM_240mW_14.dat'; BT = 0.09
#fname = '2011-05-06/DV_8x1_Dicr_dsdna_10pM_75mW_2.dat'; BT = 0.09
#fname = '2011-05-16/DV_8x1_Dicr_dsdna_10pM_1.dat'; BT = 0.09
#fname = '2011-05-17/DV_8x1_Dicr_ssdna_q100pM_3.dat'; BT = 0.09
#fname = '2011-05-18/DV_8x1_Dicr_ssdna_q100pM_sample2_4.dat'; BT = 0.09
#fname = '2011-05-18/DV_8x1_Dicr_ssdna_tmr_alexa_100pM_2.dat'; BT = 0.09
#fname = '2011-05-18/DV_8x1_Dicr_high_fret_no1_100pM_1.dat'; BT = 0.11
#fname = '2011-05-18/DV_8x1_Dicr_high_fret_no1_10pM_2.dat'; BT = 0.11

#fname = '2011-05-27/DV_8x1_Dicr_xFRET3b_100pM_0.dat'; BT = 0.11
#fname = '2011-05-27/DV_8x1_Dicr_xFRET7b_100pM_1.dat'; BT = 0.11
#fname = '2011-06-06/DV_8x1_Dicr_xFRET5b_100pM_0.dat'; BT = 0.11
#fname = '2011-06-06/DV_8x1_Dicr_xFRET10b_100pM_1.dat'; BT = 0.11 # BIOS 10bp

#fname = '2011-06-07/DV_8x1_Dicr_xFRET10b_33pM_0.dat'; BT = 0.11
fname = '2011-06-07/DV_8x1_Dicr_xFRET3b_33pM_3.dat'; BT = 0.11 # BIOS 3bp
#fname = '2011-06-07/DV_8x1_Dicr_xFRET5b_33pM_2.dat'; BT = 0.11 # BIOS 5bp
#fname = '2011-06-07/DV_8x1_Dicr_xFRET7b_33pM_1.dat'; BT = 0.11

#fname = '2011-06-15/DV_8x1_Dicr_xFRET5b_100pM_ND0_0.dat'; BT = 0.11
#fname = '2011-06-15/DV_8x1_Dicr_xFRET5b_100pM_ND0_3_1.dat'; BT = 0.11
#fname = '2011-06-15/DV_8x1_Dicr_xFRET5b_100pM_ND0_6_2.dat'; BT = 0.11
#fname = '2011-06-16/DV_8x1_Dicr_xFRET7b_25pM_ND0_0.dat'; BT = 0.11
#fname = '2011-06-16/DV_8x1_Dicr_xFRET7b_25pM_ND0_4_1.dat'; BT = 0.11
#fname = '2011-06-16/DV_8x1_Dicr_xFRET10b_25pM_ND0_6_2.dat'; BT = 0.11

## From 2011-06-17 on all measurements w/ TE 250mM NaCl IF not specified.

#fname = '2011-06-17/DV_8x1_Dicr_xFRET7b_100pM_ND0_6_0.dat'; BT = 0.11
#fname = '2011-06-17/DV_8x1_Dicr_xFRET7b_100pM_ND0_1.dat'; BT = 0.11
#fname = '2011-06-17/DV_8x1_Dicr_xFRET7b_25pM_ND0_3.dat'; BT = 0.11 # BIOS 7bp
#fname = '2011-06-17/DV_8x1_Dicr_xFRET7b_25pM_ND0_6_4.dat'; BT = 0.11

## Ron samples
#fname = '2011-09-27/rg1+rgc1_tubeSM_ND0_0.dat'; BT = 0.11; swap_D_A=True
#fname = '2011-09-27/rg1+rgc2_1o5tube_ND0_2.dat'; BT = 0.11; swap_D_A=True
#fname = '2011-09-28/rg1+rgc1+Don1-Acc1-100x_0.dat'; BT = 0.11; swap_D_A=True
#fname = '2011-09-29/rg1+rgc2_1-1_200_0.dat'; BT = 0.11; swap_D_A=True

## DCR and TEST
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_DCR_light_off_real_2.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_DCR_light_on_real_3.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_1x6mW_pix2-6-off_10.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_1x6mW_4pix_9.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_1x6mW_pix1-5_ok_8.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_1x6mW_pix1-5_7.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_1x6mW_pix2-6_6.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_1x6mW_pix4-8_5.dat'
#fname = '2011-05-04/DV_8x1_Dicr_Alexa594_1x6mW_pix3-7_4.dat'
#fname = '2011-05-09/Detector_Noise_3.dat'
# data_dir = './data/' -> Now defined in path_def.py
fname = data_dir+fname

## Example using globbing for file name
#fname = glob(data_dir+'2011-05-27/*_0.dat')[0]; BT = 0.11

## Example using GUI to select a file
#from utils.gui import gui_fname
#fname = gui_fname(data_dir)

if len(sys.argv) > 1: fname = sys.argv[1]

short_fname = shorten_fname(fname) # path with only the last subfolder
#                                                                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - FILE NAMES END -'

## From 2011-05-04/DV_8x1_Dicr_Alexa594_DCR_light_off_real_2.dat
dcr = array([ 2019.40481586, 712.32578196, 91.51662338, 45.27361157,
      51.2865131, 157.27864018, 1108.45678316, 38.08170973])

## Parameters
clock_period_ns = 12.5 # -> 80MHz
clk_p = clock_period_ns*1e-9

## Read data
pprint('\n >>>>>>>>>> PROCESSING: %s\n\n' % fname)

# OLD DATA LOAD FUNCTIONS
#ph_times, detector = load_data_ordered(fname, 300*2**20)
#ph_times_m, red = merge_green_red(ph_times, detector)
#masks = [(detector==x) for x in range(1,9)]

# NEW DATA LOAD FUNCTIONS (2x faster)
ph_times_m, red, ph_times_m8 = load_data_ordered2(fname, -1, swap_D_A)

# Background calculation
#rates8, rates4 = background_rates_calc_m(ph_times_m8, short_fname,
#        nocache=False,
#        # fit_tau keywords:
#        clk_p=clk_p, use_linregr=True, min_delay_ms=1., auto_min_delay=True,
#        debug=False) 
#rate_dd = rates8[:self.nch]
#rate_da = rates8[self.nch:]
#rate_m = rates4
#
rate_m, rate_dd, rate_ad, rate_aa = background_rates_calc_multi_c(
        fname=short_fname, ph_times_m=ph_times_m, A_em=red, clk_p=clk_p, 
        use_linregr=True, min_delay_ms=1., auto_min_delay=True,
        no_cache=False, debug=False)

# Save everything in a variable of type Data (see burstlib.py)
d = Data(fname=short_fname, nch=4, clk_p=clk_p, BT=BT, ALEX=False, 
        ph_times_m=ph_times_m, A_em=red, 
        rate_m=rate_m, rate_dd=rate_dd, rate_ad=rate_ad, rate_aa=rate_aa,
        )

# Burst search
L = 20; m = 3
d.burst_search(L=L, m=m, P_user=0.05)

# Compute corrections
d.background_correction()
d.bleed_through_correction()

# Compute FRET efficiency and stochiometry
d.calc_fret()

# Run a (optional) script for printing stats, logging and auto-plotting
ip.magic("run -i only_burst.py")

