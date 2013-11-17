# -*- coding: utf-8 -*-
import scipy.stats as S
from path_def import *

try: ip = get_ipython()
except: ip = _ip
ip.magic("run -i burst_selection.py")
ip.magic("run -i burst_plot.py")
ip.magic("run -i style.py")

#fo = '2013-03-26/DO_300p_TE50_TX_320mW_24.dat'
fo = '2013-03-26/DO_300p_TE50_TX_800mW_25.dat'

fo = data_dir+fo

clk_p = 12.5e-9
gamma = 0.45
dither = False

#do = Data(fname=fo, clk_p=clk_p, nch=8, BT=0.0, gamma=1.)
#do.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
#do.calc_bg_cache(bg_calc_exp, time_s=10, tail_min_p=0.1)
do.burst_search_t(L=10,m=10,P=None,F=6,nofret=True)
do.cal_ph_num()

## - Fit the slop pre-bg correction (bt_prebg)
## - Select a symmetric population around the slope (do_cut)
## - Fit again the (more) symmetric population to reduce skew (bt_prebg2)
## - Repeat again (bt_prebg3)
bt_prebg = do.calc_bt_from_donly(debug=True)
do_cut = Sel(do, select_bursts_for_bt_fit,BT=bt_prebg[:,0],nofret=True)
bt_prebg2 = do_cut.calc_bt_from_donly(debug=True)
do_cut2 = Sel(do, select_bursts_for_bt_fit,BT=bt_prebg2[:,0],nofret=True)
bt_prebg3 = do_cut2.calc_bt_from_donly(debug=True)
do_cut3 = Sel(do, select_bursts_for_bt_fit,BT=bt_prebg3[:,0],nofret=True)


#figure(); title('slope')
#plot(bt_prebg[:,0])
#plot(bt_prebg2[:,0])
#figure(); title('intercept')
#plot(bt_prebg[:,1])
#plot(bt_prebg2[:,1])

do.background_correction_t()
do_cut.background_correction_t()
do_cut2.background_correction_t()
do_cut3.background_correction_t()

bt = do.calc_bt_from_donly(debug=True)
bt_cut = do_cut.calc_bt_from_donly(debug=True)
bt_cut2 = do_cut2.calc_bt_from_donly(debug=True)
bt_cut3 = do_cut3.calc_bt_from_donly(debug=True)
bt_mean = mean([bt, bt_cut3], axis=0)
figure(); title('slope')
plot(bt[:,0])
plot(bt_cut[:,0])
plot(bt_cut2[:,0])
plot(bt_cut3[:,0])
plot(bt_mean[:,0], '--')

#bt = do.calc_bt_from_donly(debug=True)
#b = bt[:,0]
#AX,_ = dplot(do, scatter_da)
#dplot(do, plot_bt_overlay, bt=b, AX=AX)
#xlim(0,600); ylim(0,42)
#do_copy = do.copy()
#do_copy.update_bt(BT=b) # correction not applied yet
#print "Correction on do_copy"
#do_copy.bleed_through_correction()
#dplot(do_copy, scatter_da)
#xlim(0,600); ylim(0,42)

do.dither(lsb=1)
do.calc_fret(count_ph=False, corrections=False)
dos = Sel(do, select_bursts_nt, th1=30, nofret=True)
dos.calc_fret(count_ph=False, corrections=False)

#figure(10)
#title("Slope")
#plot(bt[:,0], label='all bursts')
#plot(bt_cut[:,0], label='selection')
#legend()
#figure(11)
#title("intercept")
#plot(bt[:,1], label='all bursts')
#plot(bt_cut[:,1], label='selection')
#legend()

#AX,_ = dplot(dos, scatter_da)
#dplot(dos, plot_bt_overlay, bt=btsb[:,0], AX=AX)

#bt=btsb[:,0]

dos.update_bt(BT=bt[:,0])
dplot(dos, hist_fret)

dos.update_bt(BT=bt_cut3[:,0])
dplot(dos, hist_fret)

dos.update_bt(BT=bt_mean[:,0])
dplot(dos, hist_fret)
