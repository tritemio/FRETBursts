# -*- coding: utf-8 -*-
import scipy.stats as S

ip = get_ipython()
ip.magic("run -i burstlib.py")
ip.magic("run -i burst_plot.py")
ip.magic("run -i style.py")

f7 = '2013-05-15/7d_New_150p_320mW_steer_3.dat'
f12 = '2013-05-15/12d_New_30p_320mW_steer_3.dat'
f17 = '2013-05-15/17d_100p_320mW_steer_1.dat'
f22 = '2013-05-15/22d_30p_320mW_steer_1.dat'
f27 = '2013-05-15/27d_50p_320mW_steer_1.dat'
fo = '2013-05-15/DO12_No2_50p_320mW_steer_1.dat'

f7, f12, f17, f22, f27, fo = [data_dir+f for f in [f7,f12,f17,f22,f27,fo]]

clk_p = 12.5e-9
gamma = 1
dither = False
fuse_ms = -1
BT = 0.036
F=6

reload_data = 0
burst_search = 0
delete_ph_times = False
gamma_calib = False
bg_plot = 0
timetrace_plot = 0
fret_plot = 0
fret_dist_plot = 1
model_plot = 0
deviance_plot = 1
save_figure = 0
nosuptitle = False
hist_weights = None

rc('savefig', dpi=250)
fret_bins = r_[-0.2:1.201:0.02]+0.01
bg_fun = bg_calc_exp_cdf

def correct_E(Er, gamma):
    return Er/(gamma-gamma*Er+Er)

## 7d
if reload_data:
    d7 = Data(fname=f7, clk_p=clk_p, nch=8, BT=BT, gamma=gamma)
    d7.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    d7.calc_bg_cache(bg_fun, time_s=10, tail_min_p=0.1)
if burst_search:
    d7.burst_search_t(L=10,m=10,P=None,F=F,dither=dither)
    df7 = d7.fuse_bursts(ms=fuse_ms)
    if delete_ph_times: (d7.delete('ph_times_m'), df7.delete('ph_times_m'))
df7.update_gamma(gamma); df7.update_bt(BT)

ds = Sel(df7, select_bursts_nda, th1=20, gamma1=0.5)
#dsc = Select_bursts(ds, select_bursts_E, E1=0.6)
#dsc.fit_E(fit_fun=gaussian_fit)
#ds.fit_from(dsc)
ds.fit_E_generic(E1=0.6, fit_fun=gaussian_fit_hist, weights=hist_weights)

dfs7 = ds.copy()
if gamma_calib:
    dfs7.calc_per_ch_gamma()
    dfs7.update_gamma(gamma*dfs7.gamma)
    dfs7c = Select_bursts(dfs7, select_bursts_E, E1=0.80)
    dfs7c.fit_E()
    dfs7.fit_from(dfs7c)
E_fit7 = dfs7.E_fit
if timetrace_plot:
    dplot(d7, timetrace_da, nosuptitle=nosuptitle); xlim(0,1); ylim(-100,100)
    if save_figure: savefig("7d_timetrace_"+d7.name()+".png"); close()
if bg_plot:
    dplot(d7, timetrace_bg, nosuptitle=nosuptitle); ylim(ymax=20e3)
    if save_figure: savefig("7d_BG_timetrace_"+d7.name()+".png"); close()
    AX, s = dplot(d7, hist_bg_fit, bp=0, scale=False, nosuptitle=nosuptitle)
    if save_figure: savefig("7d_BG_"+d7.name()+".png"); close()
if fret_plot: 
    dplot(dfs7, hist_fret, bins=fret_bins, show_fit=True, nosuptitle=nosuptitle)
    if save_figure: savefig("7d_FRET_"+d7.name()+".png"); close()

## 12d
if reload_data:
    d12 = Data(fname=f12, clk_p=clk_p, nch=8, BT=BT, gamma=gamma)
    d12.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    d12.calc_bg_cache(bg_fun, time_s=10, tail_min_p=0.1)
if burst_search:
    d12.burst_search_t(L=10,m=10,P=None,F=F,dither=dither)
    df12 = d12.fuse_bursts(ms=fuse_ms)
    if delete_ph_times: (d12.delete('ph_times_m'), df12.delete('ph_times_m'))
df12.update_gamma(gamma); df12.update_bt(BT)

ds = Sel(df12, select_bursts_nda, th1=20, gamma1=0.5)
#dsc = Select_bursts(ds, select_bursts_E, E1=0.2)
#dsc.fit_E(fit_fun=gaussian_fit)
#ds.fit_from(dsc)
ds.fit_E_generic(E1=0.2, fit_fun=gaussian_fit_hist, weights=hist_weights)

dfs12 = ds.copy()
if gamma_calib:
    dfs12.calc_per_ch_gamma()
    dfs12.update_gamma(gamma*dfs12.gamma)
    dfs12c = Select_bursts(dfs12, select_bursts_E, E1=0.45)
    dfs12c.fit_E()
    dfs12.fit_from(dfs12c)
E_fit12 = dfs12.E_fit
if timetrace_plot:
    dplot(d12, timetrace_da); xlim(0,1); ylim(-100,100)
    if save_figure: savefig("12d_timetrace_"+d12.name()+".png"); close()
if bg_plot:
    dplot(d12, timetrace_bg)
    if save_figure: savefig("12d_BG_timetrace_"+d12.name()+".png"); close()
    dplot(d12, hist_bg_fit, bp=0, scale=False)
    if save_figure: savefig("12d_BG_"+d12.name()+".png"); close()
if fret_plot: 
    dplot(dfs12, hist_fret, bins=fret_bins, show_fit=True, nosuptitle=nosuptitle)
    if save_figure: savefig("12d_FRET_"+d12.name()+".png"); close()

## 17d
if reload_data:
    d17 = Data(fname=f17, clk_p=clk_p, nch=8, BT=BT, gamma=gamma)
    d17.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    d17.calc_bg_cache(bg_fun, time_s=10, tail_min_p=0.1)
if burst_search:
    d17.burst_search_t(L=10,m=10,P=None,F=F,dither=dither)
    df17 = d17.fuse_bursts(ms=fuse_ms)
    if delete_ph_times: (d17.delete('ph_times_m'), df17.delete('ph_times_m'))
df17.update_gamma(gamma); df17.update_bt(BT)

ds = Sel(df17, select_bursts_nda, th1=30, gamma1=0.5)
#dsc = Select_bursts(ds, select_bursts_E, E1=0.05)
#dsc.fit_E(fit_fun=gaussian_fit)
#ds.fit_from(dsc)
ds.fit_E_generic(E1=0.05, fit_fun=gaussian_fit_hist, weights=hist_weights)

dfs17 = ds.copy()
if gamma_calib:
    dfs17.calc_per_ch_gamma()
    dfs17.update_gamma(gamma*dfs17.gamma)
    dfs17c = Select_bursts(dfs17, select_bursts_E, E1=0.15)
    dfs17c.fit_E()
    dfs17.fit_from(dfs17c)
E_fit17 = dfs17.E_fit
if timetrace_plot:
    dplot(d17, timetrace_da); xlim(0,1); ylim(-100,100)
    if save_figure: savefig("17d_timetrace_"+d17.name()+".png"); close()
if bg_plot:
    dplot(d17, timetrace_bg)
    if save_figure: savefig("17d_BG_timetrace_"+d17.name()+".png"); close()
    dplot(d17, hist_bg_fit, bp=0, scale=False)
    if save_figure: savefig("17d_BG_"+d17.name()+".png"); close()
if fret_plot: 
    dplot(dfs17, hist_fret, bins=fret_bins, show_fit=True, nosuptitle=nosuptitle)
    if save_figure: savefig("17d_FRET_"+d17.name()+".png"); close()

## 22d
if reload_data:
    d22 = Data(fname=f22, clk_p=clk_p, nch=8, BT=BT, gamma=gamma)
    d22.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    d22.calc_bg_cache(bg_fun, time_s=10, tail_min_p=0.1)
if burst_search:
    d22.burst_search_t(L=10,m=10,P=None,F=F,dither=dither)
    df22 = d22.fuse_bursts(ms=fuse_ms)
    if delete_ph_times: (d22.delete('ph_times_m'), df22.delete('ph_times_m'))
df22.update_gamma(gamma); df22.update_bt(BT)
dfs22 = Sel(df22, select_bursts_nda, th1=59, gamma=0.5)
#dfs22.fit_E(fit_fun=gaussian_fit)
dfs22.fit_E_generic(fit_fun=gaussian_fit_hist, weights=hist_weights)
if gamma_calib:
    dfs22.calc_per_ch_gamma()
    dfs22.update_gamma(gamma*dfs22.gamma)
    dfs22.fit_E()
E_fit22 = dfs22.E_fit
if timetrace_plot:
    dplot(d22, timetrace_da); xlim(0,1); ylim(-100,100)
    if save_figure: savefig("22d_timetrace_"+d22.name()+".png"); close()
if bg_plot:
    dplot(d22, timetrace_bg)
    if save_figure: savefig("22d_BG_timetrace_"+d22.name()+".png"); close()
    dplot(d22, hist_bg_fit, bp=0, scale=False)
    if save_figure: savefig("22d_BG_"+d22.name()+".png"); close()
if fret_plot: 
    dplot(dfs22, hist_fret, bins=fret_bins, show_fit=True, nosuptitle=nosuptitle)
    if save_figure: savefig("22d_FRET_"+d22.name()+".png"); close()

## 27d
if reload_data:
    d27 = Data(fname=f27, clk_p=clk_p, nch=8, BT=BT, gamma=gamma)
    d27.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    d27.calc_bg_cache(bg_fun, time_s=10, tail_min_p=0.1)
if burst_search:
    d27.burst_search_t(L=10,m=10,P=None,F=F,dither=dither)
    df27 = d27.fuse_bursts(ms=fuse_ms)
    if delete_ph_times: (d27.delete('ph_times_m'), df27.delete('ph_times_m'))
df27.update_gamma(gamma); df27.update_bt(BT)
dfs27 = Sel(df27, select_bursts_nda, th1=80, gamma=0.5)
#dfs27.fit_E(fit_fun=gaussian_fit)
dfs27.fit_E_generic(fit_fun=gaussian_fit_hist, weights=hist_weights)
if gamma_calib:
    dfs27.calc_per_ch_gamma()
    dfs27.update_gamma(gamma*dfs27.gamma)
    dfs27.fit_E()
E_fit27 = dfs27.E_fit
if timetrace_plot:
    dplot(d27, timetrace_da); xlim(0,1); ylim(-100,100)
    if save_figure: savefig("27d_timetrace_"+d27.name()+".png"); close()
if bg_plot:
    dplot(d27, timetrace_bg)
    if save_figure: savefig("27d_BG_timetrace_"+d27.name()+".png"); close()
    dplot(d27, hist_bg_fit, bp=0, scale=False)
    if save_figure: savefig("27d_BG_"+d27.name()+".png"); close()
if fret_plot: 
    dplot(dfs27, hist_fret, bins=fret_bins, show_fit=True, nosuptitle=nosuptitle)
    if save_figure: savefig("27d_FRET_"+d27.name()+".png"); close()

## Donly
if reload_data:
    do = Data(fname=fo, clk_p=clk_p, nch=8, BT=BT, gamma=gamma)
    do.load_multispot_cache(bytes_to_read=-1, swap_D_A=True)
    do.calc_bg_cache(bg_fun, time_s=10, tail_min_p=0.1)
if burst_search:
    do.burst_search_t(L=10,m=10,P=None,F=F,dither=dither)
    dfo = do.fuse_bursts(ms=fuse_ms)
    if delete_ph_times: (do.delete('ph_times_m'), dfo.delete('ph_times_m'))
dfo.update_gamma(gamma); dfo.update_bt(BT)
dfso = Sel(dfo, select_bursts_nda, th1=80, gamma=0.5)
#dfso.fit_E(fit_fun=gaussian_fit)
dfso.fit_E_generic(fit_fun=gaussian_fit_hist, weights=hist_weights)
if gamma_calib:
    dfso.calc_per_ch_gamma()
    dfso.update_gamma(gamma*dfso.gamma)
    dfso.fit_E()
E_fito = dfso.E_fit
if timetrace_plot:
    dplot(do, timetrace_da); xlim(0,1); ylim(-100,100)
    if save_figure: savefig("Donly_timetrace_"+do.name()+".png"); close()
if bg_plot:
    dplot(do, timetrace_bg)
    if save_figure: savefig("Donly_BG_timetrace_"+do.name()+".png"); close()
    dplot(do, hist_bg_fit, bp=0, scale=False)
    if save_figure: savefig("Donly_BG_"+do.name()+".png"); close()
if fret_plot: 
    dplot(dfso, hist_fret, bins=fret_bins, show_fit=True, nosuptitle=nosuptitle)
    if save_figure: savefig("Donly_FRET_"+do.name()+".png"); close()

## Plot 27d FRET + D-only
dfsoc = Select_bursts(dfso, select_bursts_time, time_s1=100, time_s2=160)
#dfsoc.fit_E()
dfsoc.fit_E_generic(fit_fun=gaussian_fit_hist, weights=hist_weights)
if fret_plot:
    AX, s = dplot(dfsoc, hist_fret, bins=fret_bins, lw=2, edgecolor='grey',color='grey', alpha=.3,
            histtype='stepfilled', nosuptitle=nosuptitle, show_fit=True, no_text=True)

    dplot(dfs27, hist_fret, bins=fret_bins, show_fit=True, AX=AX, color='b', alpha=0.3, 
          edgecolor='white', nosuptitle=nosuptitle)
    ylim(0,200)
    if save_figure: savefig("27d_FRETvsDonly_"+d27.name()+".png")

if fret_dist_plot:
    ## GLOBAL Plots
    font1 = {'fontname':'Liberation Sans','fontsize':16}
    font2 = {'fontname':'Arial','fontsize':16}
    rcParams["font.family"] = font2['fontname']
    fontsize = 'large'
    rcParams['xtick.labelsize']=fontsize
    rcParams['ytick.labelsize']=fontsize
    #figure(figsize=(8,5.975))
    figure(figsize=(9.1,6.6))

    dist_s_bp = [7,12,17,22,27]
    EE = correct_E(vstack([E_fit7,E_fit12,E_fit17,E_fit22,E_fit27]), 0.43)
    #EE = vstack([E_fit7,E_fit12,E_fit17,E_fit22,E_fit27])
    plot(dist_s_bp, EE, '+', lw=2, mew=1.2, ms=8, zorder=4)
    ylim(0,1); xlim(0,30); grid(True)
    Eus = r_[0.92,0.749,0.455,0.23,0.144] # See plot/ALEX sample comp.py
    Eus = r_[0.93,0.749,0.46,0.225,0.13]

    plot(dist_s_bp, Eus, '-', lw=3, mew=0, alpha=0.5, color='r', zorder=3)
    xlabel('Distance in base-pairs', fontsize=fontsize); 
    ylabel('E', fontsize=fontsize)
    legend(['CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8', u'μsALEX'], 
    #legend(['CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8'], 
           loc='best', fancybox=True, prop={'size':fontsize})
           #loc='upper right', ncol=1, frameon=False, prop={'size':fontsize})
    #title('Multi-spot smFRET dsDNA, Gamma = %.2f' % gamma)
    if deviance_plot:
        dfun = lambda E_fit: E_fit.std()*100
        dfun = lambda E_fit: 100*(E_fit.max()-E_fit.min())    
        
        text(1.2,0.05,r"$\Delta\,=\,(E_{MAX} - E_{MIN})\cdot100$", 
                ha='left',va='bottom', 
                bbox=dict(facecolor='#DCE6F2', pad=15), 
                zorder=1, fontsize=fontsize)

        D7, D12, D17, D22, D27 = dfun(EE[0]), dfun(EE[1]), dfun(EE[2]), \
                dfun(EE[3]),dfun(EE[4]),
        text(6.4,0.92,r"$\Delta$"+" = %.1f%%" %D7, ha='right',va='center', 
                bbox=dict(facecolor='white', edgecolor='white'), 
                zorder=-2, fontsize=fontsize)
        text(11.4,0.74,r"$\Delta$"+" = %.1f%%" %D12, ha='right',va='center', 
                bbox=dict(facecolor='white', edgecolor='white'), 
                zorder=-2, fontsize=fontsize)
        text(16.4,0.46,r"$\Delta$"+" = %.1f%%" %D17, ha='right',va='center', 
                bbox=dict(facecolor='white', edgecolor='white'), 
                zorder=-2, fontsize=fontsize)
        text(21.4,0.22,r"$\Delta$"+" = %.1f%%" %D22, ha='right',va='center', 
                bbox=dict(facecolor='white', edgecolor='white'), 
                zorder=-2, fontsize=fontsize)
        text(26.4,0.11,r"$\Delta$"+" = %.1f%%" %D27, ha='right',va='center', 
                bbox=dict(facecolor='white', edgecolor='white'), 
                zorder=-2, fontsize=fontsize)

    if model_plot:
        bp_rot = 2*pi/10. # rotation per bp (36 degree)
        dna_radius = 2.4e-9*0.5
        bp_pitch = 0.33e-9
        dist_lin_bp = r_[0:30:0.2]
        dist_lin = dist_lin_bp*bp_pitch

        dist_c = sqrt((dist_lin)**2 + (2*dna_radius*sin(0.5*bp_rot*dist_lin_bp))**2)
        dist_c_bp = dist_c/bp_pitch

        ## R^6 Model
        #R0 = 6.5e-9
        #plot(dist_lin_bp, 1./(1+(dist_lin/R0)**6), lw=2, alpha=0.4, color='grey', 
        #        zorder=-1)
        #text(15,0.82,r"$R_0$ = 6.5nm", ha='center',va='center', rotation=-54, 
        #        color='grey', bbox=dict(facecolor='white', edgecolor='white'), 
        #        zorder=1, fontsize=fontsize)
        
        ## Clegg R0 = 6.5nm
        R0 = 6.5e-9
        plot(dist_lin_bp, 1./(1+(dist_c/R0)**6), lw=2, alpha=0.4, color='grey', 
             zorder=-1)
        text(13.7,0.82,r"$R_0$ = 6.5nm", ha='center',va='center', rotation=-56, 
                color='grey', bbox=dict(facecolor='white', edgecolor='white'), 
                zorder=1, fontsize=fontsize)

        ## Clegg R0 = 5nm
        R0 = 5e-9
        plot(dist_lin_bp, 1./(1+(dist_c/R0)**6), lw=2, alpha=0.4,
                color='grey', zorder=-1)
        text(15.6,0.33,r"$R_0$ = 5nm", ha='center',va='center', rotation=-48.5,
                color='grey', bbox=dict(facecolor='white', edgecolor='white'), 
                zorder=1, fontsize=fontsize)
        #text(18,0.22,r"$R_0$ = 5nm", ha='center',va='center', rotation=-33, 
        #         color='grey', bbox=dict(facecolor='white', edgecolor='white'),
        #         zorder=1, fontsize=fontsize)
        
        ## Text for the 6nm curve
        #text(23.5,0.1,r"$R_0$ = 6nm", ha='center',va='center', rotation=0, 
        #         color='grey', bbox=dict(facecolor='white', edgecolor='white'),
        #         fontsize=fontsize)

    gca().set_axisbelow(True)
    if save_figure: savefig("FRET vs distance.png")

    # defined in grant.py
    #AX = plot_mburstm_8ch_global([dfs7, dfs12, dfs17, dfs22,dfs27],bins=r_[-0.2:1.2:0.05]+0.025, normed=True, show_fit=True)
    #ylim(0,4.8)

##
# Save to file
# 

#EE = correct_E(vstack([E_fit7,E_fit12,E_fit17,E_fit22,E_fit27]), 1.)
#savetxt('8spot 5samples fret fit - Gamma=1.txt', EE)

#EE = correct_E(vstack([E_fit7,E_fit12,E_fit17,E_fit22,E_fit27]), .43)
#savetxt('8spot 5samples fret fit - Gamma=0_43.txt', EE)

#savetxt('usALEX 5samples fret fit.txt', Eus)
