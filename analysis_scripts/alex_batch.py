# Configuration XM-SPCM
DONOR_ON = (2580, 330); ACCEPT_ON = (550, 2300); switch_window = 2000; swap_DA = True
    
fname_list = [
    # Selected meaningful files for 2012-08-24
    dict(fname='2012-08-24/002_12bp-Atto550-Atto647N-50nM-100uW-532nm-50uW-638nm.sm',BT=0.10),#OK
    dict(fname='2012-08-24/004_12bp-Atto550-Atto647N-50nM-200uW-532nm-100uW-638nm.sm',BT=0.10),
    dict(fname='2012-08-24/005_12bp-Atto550-Atto647N-50nM-100uW-532nm-50uW-638nm.sm',BT=0.10),
    
    # Selected meaningful files for 2012-08-27
    dict(fname='2012-08-27/014_12bp-A550-A647N-50nM-200uW-532nm-100uW-638nm.sm.sm',BT=0.10),#OK]
    dict(fname='2012-08-27/023_12bp-A550-A647N-50nM-120uW-532nm-80uW-638nm-PinH3.sm.sm',BT=0.07),
    
    # Sample degradation in (005, 006), 001 last measurement of the day
    dict(fname='2012-08-28/005_12bp-A550-A647N-165nM-100uW-532nm-50uW-638nm-PinH2.sm',BT=0.1),
    dict(fname='2012-08-28/006_12bp-A550-A647N-165nM-100uW-532nm-25uW-638nm-PinH2.sm',BT=0.1),
    dict(fname='2012-08-28/007_12bp-A550-A647N-165nM-TE250-100uW-532nm-25uW-638nm-PinH2.sm',BT=0.1),
    dict(fname='2012-08-28/008_12bp-A550-A647N-50nMno1-TE500-100uW-532nm-25uW-638nm-PinH2.sm',BT=0.1),#OK]
    dict(fname='2012-08-28/009_12bp-A550-A647N-155nM-TE250-100uW-532nm-50uW-638nm-PinH2.sm',BT=0.1),
    dict(fname='2012-08-28/001_12bp-A550-A647N-50pM-TE50-100uW-532nm-50uW-638nm-PinH2_PDM.sm.sm',BT=0.1),#OK]
    
    dict(fname='2012-10-05/001_12d_dsDNA_A550_A647N_125pM_100uW_532nm_60uW_635nm.sm',BT=0.065),
    dict(fname='2012-10-05/003_12d_dsDNA_A550_A647N_80pM_z10um_100uW_532nm_60uW_635nm.sm',BT=0.065),#OK]
]

fname_list = [
    dict(fname='N/2012-10-16/001_12d_dsDNA_TE50_100uW_532nm_40uW_635nm.sm',BT=0.07),
    dict(fname='N/2012-10-16/002_12d_dsDNA_C2_100uW_532nm_40uW_635nm.sm',BT=0.07),
    dict(fname='N/2012-10-16/003_12d_dsDNA_C2_TE50_500uW_532nm_NO_RED.sm',BT=0.07),
    dict(fname='N/2012-10-16/004_12d_dsDNA_C2_TE50_1000uW_532nm_NO_RED.sm',BT=0.07),
    dict(fname='N/2012-10-16/006_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED.sm',BT=0.07),
    dict(fname='N/2012-10-16/007_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-2.sm',BT=0.07),
    dict(fname='N/2012-10-16/008_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-3.sm',BT=0.07),
    dict(fname='N/2012-10-16/009_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-4.sm',BT=0.07),
    dict(fname='N/2012-10-16/010_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-5.sm',BT=0.07),
    dict(fname='N/2012-10-16/011_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-6.sm',BT=0.07),
    dict(fname='N/2012-10-16/012_12d_dsDNA_C2_TE50_200uW_532nm_NO_RED-7.sm',BT=0.07),
]

for line in fname_list:
    fname = line['fname']
    BT = line['BT']

    fname = alex_data_dir+fname
    short_fname = shorten_fname(fname)
    name = short_fname.replace('/','_')[:-3]

    d = process_data(fname, BT=BT, m=3, L=20, swap_DA=swap_DA)
    #print print_burst_stats(d)
    
    #d2 = Select_bursts(d, select_bursts_fret_na_bg_p, P=0.02)
    #d3 = Select_bursts(d2, select_bursts_fret_nd_bg_p, P=0.02)
    
    if not os.path.exists(fig_dir+name):
        pprint('\n - Creating the figure folder "%s" ...' % (fig_dir+name))
        os.mkdir(fig_dir+name)
        pprint('[DONE]\n')
    else:
        pprint('\n - Figure folder already present.\n')
    
    #dplot(d3, hist_bleaching, normalize=True)
    #savefig(fig_dir+name+'/'+name+'_hist_bleaching'+'.png')

    #dplot(d, hist_fret)
    #savefig(fig_dir+name+'/'+name+'_hist_fret'+'.png')

    #dplot(d, scatter_da)
    #savefig(fig_dir+name+'/'+name+'_scatter_da'+'.png')
    
    #a = 0.2
    #if d.nt[0].size > 20000: a = 0.1
    #dplot(d, hist2d_alex, scatter_alpha=a, figsize=(7,5))
    #savefig(fig_dir+name+'/'+name+'_hist2d_alex'+'.png')
    
    #dplot(d, scatter_alex)
    #savefig(fig_dir+name+'/'+name+'_scatter_alex'+'.png')

    dplot(d, hist_size, pbar=False, vmax=600)
    savefig(fig_dir+name+'/'+name+'_hist_size_nt'+'.png')
    
    df = Select_bursts(d, select_bursts_fret_na_bg_p, P=0.05)
    dff = Select_bursts(df, select_bursts_fret_nd_bg_p, P=0.05)
    
    dplot(dff, hist_bleaching, bins=arange(-0.6,0.61,0.05))
    xlim(-0.6,0.6)
    savefig(fig_dir+name+'/'+name+'_hist_bleaching_fret'+'.png')

    del d
    #del d2
    #del d3

    close('all')

    


    


