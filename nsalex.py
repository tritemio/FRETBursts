from path_def import *
from dataload.spcreader import load_spc
from utils import gui_fname

try: ip = get_ipython()
except: ip = _ip
ip.magic("run -i burst_selection.py")
ip.magic("run -i burst_plot.py")


donor_ch, accept_ch = 0, 2
clk_p = 50e-9

gain = 2.; range_ = 60e-9; bin_width_s = (range_/gain)/4096

fname = '2012-1-3/po_dsdna7_green100_red0.spc'
fname = nsalex_data_dir+fname
fname = gui_fname(nsalex_data_dir)

ddd = load_spc(fname)
ph_times = ddd[0]
det = ddd[1]
ntime = ddd[2]

d_em = (det == donor_ch)
a_em = (det == accept_ch)
assert (d_em+a_em).sum() == det.size
ph_times_m = [ph_times]
A_em = [a_em]
D_em = [d_em]

d = Data(fname=fname, clk_p=clk_p, nch=1, BT=0.1, gamma=1.)
d.add(ph_times_m=ph_times_m, A_em=A_em, nanotime=ntime, ALEX=False)

pprint("BG Calculation ... ")
d.calc_bg(bg_calc_exp, time_s=1)
pprint("DONE\n")

d.burst_search_t(L=10,m=10,P=None,F=4)
#df = d.fuse_bursts(ms=0.1)
#dfs = Select_bursts(df, select_bursts_nt, th1=15)

if 1:
    b = 8
    H = histogram(d.nanotime[d_em], bins=r_[0:4096:b])
    nt = H[0]
    bins_ns = H[1][:-1]*bin_width_s*1e9*b
    plot(bins_ns, nt, lw=1.5, alpha=0.7, label=d.name())

