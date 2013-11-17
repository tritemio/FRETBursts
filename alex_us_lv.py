#
# Load SM data file written by LV in the us-ALEX setup.
#
# Run this file from the burst/ folder as follows:
#
# >>> run -i sm/alex_us_lv.py
#

from glob import glob
from path_def import *
from dataload.smreader import load_sm

#try:
#    from burstlib import *
#except:
#    print "WARNING: Importing burstlib failed."

debug = True
compute_bg_rates = True

MONITOR_HEADER = 177
STD_HEADER = 166

## Encoding and setup parameters
donor_ch, accept_ch = 1, 0
monitor_ch = 2 # None
clock_period_ns = 12.5   # (ns) -> 80MHz
clk_p = clock_period_ns*1e-9

# Build masks for the alternating periods
def select_outer_range(times, period, edges):
    return ((times%period) > edges[0]) + ((times%period) < edges[1])

def select_inner_range(times, period, edges):
    return ((times%period) > edges[0]) * ((times%period) < edges[1])


def process_data(fname, header=STD_HEADER, **kwargs):
    ## Load data file
    print " - Loading '%s' ... " % fname
    ph_times_t, det_t = load_sm(fname,header=header) 
    #ph_times_t, det_t = load_sm(fname, header=0xAB) # 2012-09-10 measurements
    print " [DONE]\n"
    
    d = process_data_noload(ph_times_t, det_t, fname=fname, **kwargs)
    return d


def process_data_noload(ph_times_t, det_t, fname, BT=0, m=10,L=10,P=0.05, ms=2, 
        F=1., corr=True, swap_DA=False, plot_bg=False):
    period = 2*switch_window
    
    if monitor_ch is not None:
        m_ = (det_t != monitor_ch)
        det_t = det_t[m_]
        ph_times_t = ph_times_t[m_]
        
    ## Create a "short" fname string
    short_fname = shorten_fname(fname) # path with only the last subfolder

    if swap_DA: print '- Swaping D and A'
    if swap_DA: D_ch, A_ch = accept_ch, donor_ch
    else: D_ch, A_ch = donor_ch, accept_ch

    # Create masks for donor and acceptor channels
    donor_t_mask = (det_t==D_ch)
    accept_t_mask = (det_t==A_ch)
    print "#donor: %d  #acceptor: %d \n" % (donor_t_mask.sum(), 
            accept_t_mask.sum())

    # Build masks for excitation windows
    donor_ex_t = select_outer_range(ph_times_t, period, DONOR_ON)
    accept_ex_t = select_inner_range(ph_times_t, period, ACCEPT_ON)

    # Safety check: each ph is either D or A ex (not both)
    assert (donor_ex_t*accept_ex_t == False).any() 
    
    ## - - - - - - - - - - - - - - - - - 
    ## Ph. stream for burst search: uncomment only one of (1), (2) and (3)
    ##
    # (1) All photon burst search
    #mask = ones(ph_times_t.size, dtype=bool)
    #ph_burst_search = "all_ph"  # TODO: use broader EXCITATION selection in
    #                            #       this case

    # (2) Burst search on donor_exitation
    #mask = donor_excitation
    #ph_burst_search = "donor_ex"

    # (3) Burst search on donor_excitation+acceptor_excitation
    mask = donor_ex_t + accept_ex_t
    
    #ph_burst_search = "donor_accept_ex"

    # Assign the new ph selection mask
    ph_times = ph_times_t[mask]
    d_em = donor_t_mask[mask]
    a_em = accept_t_mask[mask]
    d_ex = donor_ex_t[mask]
    a_ex = accept_ex_t[mask]
    
    assert d_em.sum() + a_em.sum() == ph_times.size
    assert (d_em * a_em).any() == False
    assert a_ex.size == a_em.size == d_ex.size == d_em.size == ph_times.size
    
    ph_times_m, A_ex, D_ex, A_em, D_em = [ph_times],[a_ex],[d_ex],[a_em],[d_em]
    d = Data(fname=fname, nch=1, clk_p=clk_p, BT=BT, gamma=1., ALEX=True,
            D_ON=DONOR_ON, A_ON=ACCEPT_ON, switch_window=switch_window,
            ph_times_m=ph_times_m,
            D_em=D_em, A_em=A_em, D_ex=D_ex, A_ex=A_ex,
            # Needed for alternation hist
            ph_times_t=[ph_times_t], D_em_t = [donor_t_mask]
            )
    assert d.ph_times_m[0].size == d.A_em[0].size
    assert d.switch_window*2 == period

    # Safety check: each ph is either D or A ex (not both)
    assert (d.A_ex[0]*d.D_ex[0]).any() == False
    assert (d.A_ex[0]+d.D_ex[0]).all() # D + A ex are ALL ph

    # Safety check: D_ex ph only in D_ON and A_ex ph only in A_ON
    assert (d.ph_times_m[0][d.A_ex]%period > d.A_ON[0]).all()
    assert (d.ph_times_m[0][d.A_ex]%period < d.A_ON[1]).all()
    
    ph_wrap = d.ph_times_m[0][d.D_ex]%period
    assert ((ph_wrap > d.D_ON[0]) + (ph_wrap < d.D_ON[1])).all()
    assert (ph_wrap < switch_window*2).all()
    
    # each ph is either D or A em (not both)
    assert (d.A_em[0]*d.D_em[0]).any() == False
    assert (d.A_em[0]+d.D_em[0]).all() # D + A em are ALL ph
    
    # Compute BG rates
    #d.compute_bg_rates_mch()
    d.calc_bg(bg_calc_exp, time_s=10, tail_min_p=0.25)

    # Burst search
    #d.burst_search(L=L,m=m,P_user=P, nofret=True)
    d.burst_search_t(L=L, m=m, P=P, F=F, nofret=True)

    # (Optionally) fuse bursts less distant than X ms
    if ms >= 0: d = d.fuse_bursts(ms=ms, process=False)
    
    # Compute FRET efficiency, stochiometry and E-S histogram
    d.calc_fret(count_ph=True, corrections=corr)

    return d

if __name__ == '__main__':
    try:
        ip = get_ipython()
    except:
        ip = _ip
    #ip.magic("run -i burstlib.py")
    ip.magic("run -i burst_selection.py")
    ip.magic("run -i burst_plot.py")
    ip.magic("run -i alex_plot.py")

    ## SM Data files
    
    # Configuration 1: for files until 2012-07-31
    #DONOR_ON = (1720, 450); ACCEPT_ON = (700, 1400); 
    #switch_window = 1000; swap_DA = False
    
    #fname = '2012-05-09/001-xfret-10bp-25pm.sm'; 
    #fname = '2012-05-11/002-xfret-10bp-6pm.sm' # NO SYNC
    #fname = '2012-05-11/003-xfret-only Buffer.sm'
    #fname = '2012-05-14/001-2nd try 6m.sm'
    #fname = '2012-05-14/002-buffer-aom-IN.sm'
    #fname = '2012-05-14/003-buffer-aom-OUT.sm'
    #fname = '2012-05-16/001-6pm2nd try.sm'
    #fname = '2012-06-06/001166green-69red-XFRET10bp-6.5pm.sm'
    #fname = '2012-06-06/001sanyoonsample.sm'
    #fname = '2012-06-06/002150green-92red-25pm-10bp.sm'
    #fname = '2012-06-06/003150green-92red-25pm-10bp.sm'
    #fname = '2012-06-06/004150green-92red-25pm-10bp.sm'
    #fname = '2012-06-06/005150green-92red-25pm-10bp.sm'
    #fname = '2012-06-06/006150green-red-25pm-10bp.sm'
    #fname = '2012-06-06/00738-85green-38-96red-25pm-10bp.sm'
    #fname = '2012-06-06/00885green-90red-25pm-10bp.sm'
    
    # Configuration 2
    #DONOR_ON = (1530, 240); ACCEPT_ON = (540, 1250); 
    #switch_window = 1000; swap_DA = False 
    
    #fname = '2012-07-03/001test1.sm'
    #fname = '2012-07-03/002test2.sm'
    #fname = '2012-07-03/009Sangyoon sample_new2.sm'
    #fname = '2012-07-03/010Sangyoon sample_new3.sm'
    #fname = '2012-07-05/001test3_SY.sm'
    #fname = '2012-07-05/002test3_SY.sm'
    #fname = '2012-07-05/021Sangyoon sample_new3.sm'
    #fname = '2012-07-05/022test2.sm'
    #fname = '2012-07-25/001_Cy3b_DNA-test.sm'
    #fname = '2012-07-25/002_Cy3b_1mw-DNA-test.sm'
    #fname = '2012-07-25/009cy3b-1test.sm'
    #fname = '2012-07-25/010cy3b-2test.sm'
    #fname = '2012-07-25/011cy3b-3test.sm'
    #fname = '2012-07-25/cy3b-1test'
    #fname = '2012-07-30/005cy3bDNA-run1.sm'
    #fname = '2012-07-30/006cy3bDNA-run2.sm'
    #fname = '2012-07-30/008cy3bDNA-20Micro-PD.sm'
    #fname = '2012-07-30/010cy3bDNA-40Micro-PD.sm'
    #fname = '2012-07-30/011cy3bDNA-60Micro-PD.sm'
    #fname = '2012-07-31/001cy3bDNA-PD-CENTER.sm'
    #fname = '2012-07-31/003cy3bDNA-PD-50micro.sm'
    #fname = '2012-07-31/004cy3bDNA-PD-80micro.sm'
    #fname = '2012-07-31/005cy3bDNA-PD-90micro.sm'
    #fname = '2012-07-31/006cy3bDNA-PD-110micro.sm'
    #fname = '2012-07-31/007cy3bDNA-PD-120micro.sm'
    #fname = '2012-07-31/008cy3bDNA-PD-140micro.sm'
    #fname = '2012-07-31/009cy3bDNA-PD-150micro.sm'
    #fname = '2012-07-31/010cy3bDNA-PD-150micro.sm'
    
    #fname = '2012-08-01/001new sample- test1.sm' # Conf 2
    #fname = '2012-08-01/002new sample- test2.sm'
    #fname = '2012-08-01/003bleeching test .sm'
    #fname = '2012-08-01/004_Cy3b_bleeching test2.sm.sm' # Conf 2
    
    # Configuration XM-SPCM
    #DONOR_ON = (2580, 330); ACCEPT_ON = (550, 2300); 
    #switch_window = 1000; swap_DA = True
    
    # Selected meaningful files for 2012-08-24
    #fname='2012-08-24/002_12bp-Atto550-Atto647N-50nM-100uW-532nm-50uW-638nm.sm';BT=0.10#OK
    #fname = '2012-08-24/004_12bp-Atto550-Atto647N-50nM-200uW-532nm-100uW-638nm.sm';BT=0.10
    #fname = '2012-08-24/005_12bp-Atto550-Atto647N-50nM-100uW-532nm-50uW-638nm.sm';BT=0.10
    
    # Selected meaningful files for 2012-08-27
    #fname='2012-08-27/014_12bp-A550-A647N-50nM-200uW-532nm-100uW-638nm.sm.sm';BT=0.10#OK
    #fname='2012-08-27/023_12bp-A550-A647N-50nM-120uW-532nm-80uW-638nm-PinH3.sm.sm';BT=0.07
    
    # Sample degradation in (005, 006), 001 last measurement of the day
    #fname='2012-08-28/005_12bp-A550-A647N-165nM-100uW-532nm-50uW-638nm-PinH2.sm';BT=0.1
    #fname ='2012-08-28/006_12bp-A550-A647N-165nM-100uW-532nm-25uW-638nm-PinH2.sm';BT=0.1
    #fname ='2012-08-28/007_12bp-A550-A647N-165nM-TE250-100uW-532nm-25uW-638nm-PinH2.sm';BT=0.1
    #fname ='2012-08-28/008_12bp-A550-A647N-50nMno1-TE500-100uW-532nm-25uW-638nm-PinH2.sm';BT=0.1#OK
    #fname = '2012-08-28/009_12bp-A550-A647N-155nM-TE250-100uW-532nm-50uW-638nm-PinH2.sm';BT=0.1
    #fname ='2012-08-28/001_12bp-A550-A647N-50pM-TE50-100uW-532nm-50uW-638nm-PinH2_PDM.sm.sm';BT=0.1#OK
    
    # DATA FILES NOT tested
    #fname = '2012-09-10/001_dsDNA-ATTO647N_100uW_532nm_60uW_635nm.sm'
    #fname = '2012-09-10/002_OnlyDye-ATTO647N-TE50_100uW_532nm_60uW_635nm.sm'
    #fname = '2012-09-10/003_OnlyDye-ATTO647N-TE50_100uW_532nm_60uW_635nm.sm'
    #fname = '2012-09-10/004_OnlyDye-ATTO647N-TE50_100uW_532nm_60uW_635nm.sm'
    #fname = '2012-09-11/001DSDNA-ATTO647N-100uw-green-60uw-Red.sm'
    #fname = '2012-09-11/002onlydye-10^-4-ATTO647N-TE50-100uw-green-60uw-Red.sm'
    #fname = '2012-09-11/003onlydye-10^-4-ATTO647N-TE50-100uw-green-60uw-Red.sm'
    #fname = '2012-09-11/004onlydye-10^-3-ATTO647N-TE50-100uw-green-60uw-Red.sm'
    #fname = '2012-09-11/004onlydye-10^-4-ATTO647N-TE50-100uw-green-60uw-Red.sm'
    #fname = '2012-09-11/005onlydye-10^-4-ATTO647N-TE50-100uw-green-60uw-Red.sm'
    #fname = '2012-09-11/006onlydye-10^-4-ATTO647N-TE50-100uw-green-60uw-Red.sm'
    #fname = '2012-09-11/007onlydye-10^-4-ATTO647N-TE50-100uw-green-60uw-Red.sm'
    
    #fname='2012-10-05/001_12d_dsDNA_A550_A647N_125pM_100uW_532nm_60uW_635nm.sm';BT=0.065
    #fname ='2012-10-05/002_TE50_z10um_100uW_532nm_60uW_635nm.sm';BT=0.065
    #fname ='2012-10-05/003_12d_dsDNA_A550_A647N_80pM_z10um_100uW_532nm_60uW_635nm.sm';BT=0.065#OK
    #fname ='2012-10-05/004_TE50_z10um_100uW_532nm_60uW_635nm.sm';BT=0.065
    #fname='2012-10-11/003_12d_sdDNA_2-5nM_from08-28-tube-prep10-11-_100uW_532nm_40uW_635nm.sm';BT = 0.1
    #fname='2012-10-11/004_12d_sdDNA_2-5nM_chamber_from08-28-tube-prep10-11-_100uW_532nm_40uW_635nm.sm';BT = 0.1
    #fname='2012-10-11/005_12d_sdDNA_2-5nM_chamber2_from08-28-tube-prep10-11-_100uW_532nm_40uW_635nm.sm';BT = 0.1

    #fname='2012-10-15/001_12d_dsDNA_6nM_TE250_100uW_532nm_40uW_635nm.sm';BT=0.1
    #fname='2012-10-15/002_12d_dsDNA_6nM_TE250_100uW_532nm_40uW_635nm.sm';BT=0.1
    #fname='2012-10-15/003_12d_dsDNA_6nM_TE250_100uW_532_blocked_red.sm';BT=0.1
    #fname='2012-10-15/004_12d_dsDNA_6nM_TE250_200uW_532nm_40uW_635nm.sm';BT=0.1
    
    #fname='N/2012-10-16/001_12d_dsDNA_TE50_100uW_532nm_40uW_635nm.sm';BT=0.07
    #fname='N/2012-10-16/002_12d_dsDNA_C2_100uW_532nm_40uW_635nm.sm';BT=0.07
    #fname='N/2012-10-16/003_12d_dsDNA_C2_TE50_500uW_532nm_NO_RED.sm';BT=0.07
    #fname='N/2012-10-16/004_12d_dsDNA_C2_TE50_1000uW_532nm_NO_RED.sm';BT=0.07
    #fname='N/2012-10-16/006_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED.sm';BT=0.07
    #fname='N/2012-10-16/007_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-2.sm';BT=0.07
    #fname='N/2012-10-16/008_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-3.sm';BT=0.07
    #fname='N/2012-10-16/009_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-4.sm';BT=0.07
    #fname='N/2012-10-16/010_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-5.sm';BT=0.07
    #fname='N/2012-10-16/011_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-6.sm';BT=0.07
    #fname='N/2012-10-16/012_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-7.sm';BT=0.07
    #fname='N/2012-10-16/013_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-8.sm';BT=0.07

    #fname='Yazan/2012-10-19/001_'; BT=0.07
    
    #fname='2012-10-24/006_12d_dsDNA_2nM_TE50_TX_12p_Suc_'; BT = 0.07
    #fname='2012-10-24/008_12d_dsDNA_2nM_AFTER_MS'; BT = 0.07;
    
    ## Uncomment these 2 lines if you selected a filename above
    #fname = alex_data_dir+fname
    #fname = glob(fname+'*')[0].replace(os.path.sep, '/')

    ## Uncomment these lines to select a file using a GUI
    from utils import gui_fname
    fname = gui_fname(alex_data_dir)
    # Configuration XM-SPCM
    DONOR_ON = (2850, 580); ACCEPT_ON = (930, 2580); 
    switch_window = 2000; swap_DA = True
    BT=0.1
     
    ## Eitan CONF
    #DONOR_ON = (1530, 280); ACCEPT_ON = (580, 1275); switch_window = 1000;
    #swap_DA = True; BT = 0.06
    
    ## ATTO647N-Only conf
    #DONOR_ON = (3000, 330); ACCEPT_ON = (800, 2500); 
    #switch_window = 2000; swap_DA = True
    
    ## Meas 2012-11-9
    #DONOR_ON = (2800, 580); ACCEPT_ON = (840, 2600); 
    #switch_window = 2000; swap_DA = True; BT = 0.06

    #A_frac = (ACCEPT_ON[1]-ACCEPT_ON[0])/(2.*switch_window)
    #D_frac = (2*switch_window-DONOR_ON[0]+DONOR_ON[1])/(2.*switch_window)

    ### Short name measurements
    #short_fname = shorten_fname(fname)
    #name = short_fname.replace('/','_')[:-3]
    ##if not os.path.exists(fig_dir+name): os.makedirs(fig_dir+name)

    ### Load data file
    #print " - Loading '%s' ... " % fname
    #ph_times_t, det_t = load_sm(fname, header=0xA6) # usual header length
    #print " [DONE]\n"
    
    ## Data manipulation
    if 0:
        ip.magic("run -i sim/dcr.py")
        ch = accept_ch
        DCR = 500
        dcr = gen_dcr(DCR, 600, clk_p).astype(ph_times_t.dtype)
        n = dcr.size
        ph_times_t = hstack([ph_times_t, dcr])
        det_t = hstack([det_t, ((ch+swap_DA)%2)*ones(n, dtype=det_t.dtype)])
        argsort = ph_times_t.argsort()
        ph_times_t = ph_times_t[argsort]
        det_t = det_t[argsort]
    #ph_times_t, det_t = ph_times_t[::2], det_t[::2]
    
    #d = process_data_noload(ph_times_t, det_t, BT=BT,m=3,L=10,swap_DA=swap_DA)
    
    d = process_data(fname, header=STD_HEADER, BT=BT, swap_DA=swap_DA,
                     m=3, L=40, F=1., P=0.95, ms=0, corr=True, plot_bg=False) 
    
    #d = process_data(fname, BT=BT, m=10, L=10, P=0.05, ms=0, corr=False,
    #        swap_DA=swap_DA) 
    #do = d
    #do.calc_fret(corrections=True)
    #d = Select_bursts(do, select_bursts_nt_bg_p, P=0.01, F=3)
    #do.calc_fret(corrections=True)
    
    #print print_burst_stats(d)
    
    #d2 = Select_bursts(d, select_bursts_fret_na_bg_p, P=0.01)
    #d3 = Select_bursts(d2, select_bursts_fret_nd_bg_p, P=0.01)
    #dplot(d3, hist_bleaching, normalize=True)
    
    #dplot(d, hist_fret)
    #savefig(fig_dir+name+'/'+name+'_hist_fret'+'.png')
    
    #dplot(d, scatter_da)
    #savefig(fig_dir+name+'/'+name+'_scatter_da'+'.png')
    
    #dplot(d, hist2d_alex)
    #savefig(fig_dir+name+'/'+name+'_hist2d_alex'+'.png')
    
    #dplot(d, hist_bleaching, bins=arange(-0.6,0.61,0.05))
    #xlim(-0.6,0.6)
    #savefig(fig_dir+name+'/'+name+'_hist_bleaching'+'.png')
    
    
    #df = Select_bursts(d, select_bursts_fret_na_bg_p, P=0.05)
    #dff = Select_bursts(df, select_bursts_fret_nd_bg_p, P=0.05)
    
    #dplot(dff, hist_bleaching, bins=arange(-0.6,0.61,0.05))
    #xlim(-0.6,0.6)
    #savefig(fig_dir+name+'/'+name+'_hist_bleaching_fret'+'.png')

    #dplot(dff, scatter_da)
    #savefig(fig_dir+name+'/'+name+'_scatter_da_fret'+'.png')

    ## Plot the alternation histogram
    #plot_alternation_hist(d)

    ## Plot the 1D FRET histogram
    #dplot(d, hist_fret)

    ## Plot the scatter plot acceptor vs donor
    #dplot(d, scatter_da)

    ## Plot the burst size distribution
    #dplot(d, hist_size)

    ## Plot the skewness index (bleaching)
    #dplot(d, hist_bleaching)

    ## Plot the ALEX histogram
    #AH = AlexHist(d, bin_step=0.05)
    #alexhist_imshow(AH, vmax=80)
