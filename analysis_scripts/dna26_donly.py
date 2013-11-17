# -*- coding: utf-8 -*-
import scipy.stats as S
from path_def import *

try: ip = get_ipython()
except: ip = _ip
ip.magic("run -i burst_selection.py")
ip.magic("run -i burst_plot.py")
ip.magic("run -i style.py")

fo = '2013-03-26/DO_300p_TE50_TX_320mW_24.dat'
#fo = '2013-03-26/DO_300p_TE50_TX_800mW_25.dat'

fo = data_dir+fo

clk_p = 12.5e-9
experim = True

do = Data(fname=fo, clk_p=clk_p, nch=8, BT=0.0, gamma=1.)
do.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
do.calc_bg_cache(bg_calc_exp, time_s=10, tail_min_p=0.1)
do.burst_search_t(L=10,m=10,P=None,F=6,nofret=True)
do.cal_ph_num()

if experim:
    bt_prebg = do.calc_bt_from_donly(debug=True)
    do_cut = Sel(do, select_bursts_for_bt_fit,BT=bt_prebg[:,0],nofret=True)
    bt_prebg2 = do_cut.calc_bt_from_donly(debug=True)
    do_cut2 = Sel(do, select_bursts_for_bt_fit,BT=bt_prebg2[:,0],nofret=True)
    bt_prebg3 = do_cut2.calc_bt_from_donly(debug=True)
    do_cut3 = Sel(do, select_bursts_for_bt_fit,BT=bt_prebg3[:,0],nofret=True)

## Apply BG correction to nd, na, nt...
do.background_correction_t()
if experim: do_cut3.background_correction_t()

## Fit BT after BG subtraction
bt = do.calc_bt_from_donly(debug=True)
if experim:
    bt_cut3 = do_cut3.calc_bt_from_donly(debug=True)
    bt_mean = mean([bt, bt_cut3], axis=0)

b = bt[:,0]
if experim:
    b = bt_cut3[:,0]
    b = bt_mean[:,0]

## Apply the BT correction
do.update_bt(BT=b) 

## Finally compute FRET
do.dither(lsb=1) # apply dithering because the threshold is very low
do.calc_fret(count_ph=False, corrections=True) # corrections will not be
                                               # reapplied

dos = Sel(do, select_bursts_nt, th1=30)
#dos.calc_fret(count_ph=False, corrections=False)

#dos.update_bt(BT=bt[:,0])
dplot(do, hist_fret)
dplot(dos, hist_fret)

def acprint(a, decimals=5):
    s = ",".join(["%.5f" % x for x in a])
    return "r_["+s+"]"

BT800 = r_[0.03560,0.03514,0.03839,0.03520,0.03343,0.03189,0.03401,0.03056]
BT800_E = r_[0.04049,0.04318,0.04144,0.03854,0.03657,0.03614,0.03676,0.03359]
BT800_EM = r_[0.03804,0.03916,0.03991,0.03687,0.03500,0.03402,0.03538,0.03207]

BT300 = r_[0.03370,0.03122,0.03536,0.03516,0.03173,0.02717,0.03149,0.03163]
BT300_E = r_[0.03935,0.04213,0.03868,0.03875,0.03609,0.03305,0.03488,0.03470]
BT300_EM = r_[0.03653,0.03668,0.03702,0.03696,0.03391,0.03011,0.03319,0.03316]

