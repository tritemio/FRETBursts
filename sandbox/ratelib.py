
def rate_ph_da(fname, kcps=True):
    d = Data(fname=fname, clk_p=12.5e-9, nch=8, BT=0., gamma=1.)
    d.load_multispot(bytes_to_read=-1, swap_D_A=True)
    rates_a, rates_d = zeros(8), zeros(8)
    for i,(ph,a) in enumerate(zip(d.ph_times_m,d.A_em)):
        if a.any():
            rates_a[i] = 1.*ph[a].size/(ph[a][-1]-ph[a][0])
    for i,(ph,a) in enumerate(zip(d.ph_times_m,d.A_em)):
        if (-a).any(): 
            rates_d[i] = 1.*ph[-a].size/(ph[-a][-1]-ph[-a][0])
    rates_d /= d.clk_p
    rates_a /= d.clk_p
    if kcps:
        rates_d /= 1e3
        rates_a /= 1e3
    return rates_d, rates_a

if __name__ == "__main__":
    try: ip = get_ipython()
    except: ip = _ip
    ip.magic("cd C:\Data\Antonio\software\burst")
    ip.magic("run -i burst_selection.py")
    ip.magic("run -i burst_plot.py")
    fname = data_dir+'2013-05-21/Orig-Pos_LCOS_43_5-2_0.dat'
    rd,ra = rate_ph_da(fname)
    ch = r_[1:9]   
    plot(ch,rd, 'g', lw=2)
    plot(ch,ra, 'r', lw=2)
    ylim(0); ylabel("kcp")
     