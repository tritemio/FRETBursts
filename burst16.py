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
# For all the 8-channel plots we call the function dplot() 
# passing the fuction name of the desired plot type. For ex. to plot the fret
# histogram type:
#
# >>> dplot(d, fun=hist_fret)
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
#from dataload.multi_ch_reader import *

## Loading the lib functions related to the burst search
ip.magic("run -i burstlib.py")
#ip.magic("run -i burst_selection.py")

## Loading the plotting functions
ip.magic("run -i burst_plot.py")

## Loading default plotting style
ip.magic("run -i style.py")

## - - - - - - - - - - - - - - - - - - - - - - - - - - - FILE NAME START -.

## NOTE: 2012-10-24 MEMO on detector numbering
# A SPAD array is connected to 1-8ch of the FPGA breakout board
# D SPAD array is connected to the 9-16ch
# '-> THUS we need to swap D and A to maintain the convention D first
#
# There is correspondence from the FPGA breakout board and the LV program
# 
# Until 2012-10-03 (included) SPAD A was mislabeled (8->1, 7->2, etc...).
# Initially, up to 2012-10-03, was connected in 8..1 order so the REAL order
# was correct and NO remap is needed for those measurement.
#
# On 2012-10-04 the A SPAD ch were inverted to match the (wrong) labeling.
# THUS we need to remap_A=True for those measurement.
#
# On 2012-10-25 the A SPAD was reconnected as initially (and relabeled).
# From this date on we don't need to remap anymore (remap_A=False)
#
##

## Data file path                         
# Conventional order: 1-8 donor ch, 9-16 acceptor ch

swap_D_A=True; BT=0.045

## NOTE FROM 2012-10-02 to 10-04 use swap_D_A=True and remap_A=False
#remap_A=False # USE until 2012-10-03
#fname = '2012-10-02/12d_DSDNA_250pM_0'
#fname = '2012-10-03/12d_DSDNA_125pM_0'; BT=0.065
#fname = '2012-10-03/12d_DSDNA_125pM_slight_green_realign_1'
#fname = 'det_num_test_0'


## NOTE FROM 2012-10-04 to 10-25 swap_D_A=True and remap_A=True
#remap_A=True # USE starting from 2012-10-04
#fname = '2012-10-11/12d_dsDNA_2-5nM_HC_0'
#fname = '2012-10-11/12d_dsDNA_2-5nM_SM_1'
#fname = '2012-10-11/12d_dsDNA_2-5nM_SM_2'
#fname = '2012-10-11/12d_dsDNA_2-5nM_nochamber_SM_ND1_3'
#fname = '2012-10-11/12d_dsDNA_2-5nM_nochamber_SM_5'
#fname = '2012-10-11/12d_dsDNA_2-5nM_nochamber_SM_ND0.4_6'
#fname = '2012-10-11/12d_dsDNA_2-5nM_nochamber2_SM_ND0.4_7'
#fname = '2012-10-11/12d_dsDNA_2-5nM_nochamber3_SM_ND0_10'

#fname = '2012-10-18/R6G_40p_Sucrose_1'
#fname = '2012-10-18/R6G_40p_Sucrose_ND0-4_2'
#fname = '2012-10-18/R6G_40p_Sucrose_maxPower_4'
#fname = '2012-10-18/Cy3B_40p_Sucrose_1000X_maxPower_5'
#fname = '2012-10-18/Cy3B_40p_Sucrose_2000X_ND1_6'
#fname = '2012-10-18/Cy3B_40p_Sucrose_2000X_ND0.4_7'
#fname = '2012-10-24/_12d_2n_TX_12p_Sucr_ND0.4_8'
#fname = '2012-10-24/_12d_2n_TX_12p_Sucr_ND1'
#fname = '2012-10-24/_12d_2n_TX_12p_Sucr_ND0_10'; BT=0.001

#fname = '2012-10-26/12d_TX_TE50_ND0.4_3'; BT=0.05
#fname = '2012-10-26/12d_TX_TE50_ND0.4_4'; BT=0.035
#fname = '2012-10-26/12d_TX_TE50_ND0_5'; BT=0.05

#fname = '2012-11-27/22d_TE50_ND0-4___0'
#fname = '2012-11-27/22d_TE50_ND0___1'
#fname = '2012-11-27/22d_TE50_ND1___2'


#fname = '2012-11-28/17d_TE50_ND1___2.dat'
#fname = '2012-11-28/17d_TE50_ND0-4___3.dat'
#fname = '2012-11-28/17d_TE50_ND1_z2___5.dat'
#fname = '2012-11-28/17d_TE50_ND0-4_z2___6.dat' # BEST 17d
#fname = '2012-11-28/22d_TE50_ND1_z2___9.dat'
#fname = '2012-11-28/22d_TE50_ND0-4_z3___10.dat'
#fname = '2012-11-28/22d_TE50_ND0_z3___11.dat'

fname = '2012-12-13/12d ND04_6.dat'
fname = '2012-12-13/27d ND04_0.dat'

## NOTE FROM 2012-10-25 onward use swap_D_A=True and no remapping
#remap_A=False

#fname = data_dir+fname
#fname = glob(fname+'*')[0]

## Example using GUI to select a file
from utils.gui import gui_fname
fname = gui_fname(data_dir)

if len(sys.argv) > 1: fname = sys.argv[1]

short_fname = shorten_fname(fname) # path with only the last subfolder
#                                                                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - FILE NAMES END -'

## From 2011-05-04/DV_8x1_Dicr_Alexa594_DCR_light_off_real_2.dat
dcr = array([ 2019.40481586, 712.32578196, 91.51662338, 45.27361157,
      51.2865131, 157.27864018, 1108.45678316, 38.08170973])

## Parameters
clk_p = 12.5e-9 # -> 80MHz
BT300_EM = r_[0.03653,0.03668,0.03702,0.03696,0.03391,0.03011,0.03319,0.03316]

## Read data
pprint('\n >>>>>>>>>> PROCESSING: %s\n\n' % fname)

d = Data(fname=fname, clk_p=clk_p, nch=8, BT=BT, gamma=1.)
#d26 = Data(fname=fname, clk_p=clk_p, nch=8, BT=BT300_EM, gamma=0.5)
d.load_multispot_cache(bytes_to_read=-1, swap_D_A=swap_D_A)
d.calc_bg_cache(bg_calc_exp, time_s=5, tail_min_p=0.1)
#d.calc_bg_cache(bg_calc_gauss, time_s=10)

d.burst_search_t(L=10,m=10,P=None,F=6, ph_sel='DA')
#df = d.fuse_bursts(ms=0.0)
#dfs = Select_bursts(df, select_bursts_nda, th1=50, gamma=1.)
#dfs.update_gamma(0.5)
#dfs.fit_E()

#dplot(dfs, hist_fret, bins=r_[-0.2:1.21:0.04])
#xlim(-0.2,1.1)
#ylim(0,300)
 
# Filters out donor-only and acceptor-only population (Prob{BG>=na} < P)
#d2 = Select_bursts(d, select_bursts_na_bg_p, P=0.05)
#d3 = Select_bursts(d2, select_bursts_nd_bg_p, P=0.05)

def split(d, t=60, n=4):
    Sel, sel_time = Select_bursts, select_bursts_time
    D = [Sel(d, sel_time, time_s1=t*i, time_s2=t*(i+1)) for i in range(n)]
    return D
    
def hack_plot(D):
    for i,d in enumerate(D):
        dplot(d, hist_fret, bins=r_[-0.2:1.21:0.04])
        #xlim(-0.2,1.1); ylim(0,300)
        #savefig(d.name()+"_TOT_slice%d.png"%(i+1))
        xlim(-0.2,1.1); ylim(0,80)
        savefig(d.name()+"_TEST_slice%d.png"%(i+1))
        close()
