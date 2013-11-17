# -*- coding: utf-8 -*-
from path_def import * # loads data_dir variable

try: ip = get_ipython()
except: ip = _ip
ip.magic("run -i burst_selection.py")
ip.magic("run -i burst_plot.py")
ip.magic("run -i style.py")

fo = '2013-05-09/12d_DO_150p_320mW_28.dat'
fos = '2013-05-09/12d_DO_150p_320mW_steering_29.dat'
fo = data_dir+fo
fos = data_dir+fos

clk_p = 12.5e-9

do = Data(fname=fo, clk_p=clk_p, nch=8, BT=0.0, gamma=1.)
do.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
do.calc_bg_cache(bg_calc_exp_cdf, time_s=10, tail_min_p=0.1)
do.burst_search_t(L=10,m=10,P=None,F=6,nofret=True, ph_sel='D')
do.cal_ph_num()

do60 = Sel(do, select_bursts_nda, th1=60, nofret=1)
BT_ML = do60.calc_bt_from_donly_ML()
BT_LS = do60.calc_bt_from_donly_LS()
BT_LR = do60.calc_bt_from_donly_LR()

dos = Data(fname=fos, clk_p=clk_p, nch=8, BT=0.0, gamma=1.)
dos.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
dos.calc_bg_cache(bg_calc_exp_cdf, time_s=10, tail_min_p=0.1)
dos.burst_search_t(L=10,m=10,P=None,F=6,nofret=True, ph_sel='D')
dos.cal_ph_num()

dos60 = Sel(dos, select_bursts_nda, th1=60, nofret=1)
BTs_ML = dos60.calc_bt_from_donly_ML()
BTs_LS = dos60.calc_bt_from_donly_LS()
BTs_LR = dos60.calc_bt_from_donly_LR()