def hist_green_red(i,b,d, max_size=100, bins=None, **kwargs):
    if bins is None: bins=arange(d.L,max_size,2)
    mask_green, mask_red = mask_burst_green_red(d, i)
    fret_eff = fret_efficiency(d, i)
    print " Fret eff CH %d: %1.2f" % (i+1, fret_eff)
    ng1 = sort(d.nd[i][mask_green]); ng1 = ng1[ng1<max_size]
    ng2 = sort(d.nd[i][mask_red]); ng2 = ng2[ng2<max_size]
    histog1, bins1 = histogram(ng1, bins=bins, **kwargs)
    histog2, bins2 = histogram(ng2, bins=bins, **kwargs)
    bar(bins1[:-1], histog1, bins1[1]-bins1[0], color='green', alpha=0.7)
    bar(bins2[:-1], histog2, bins2[1]-bins2[0], color='red', alpha=0.7)
    text(0.6,0.8,"Fret eff. = %1.2f" % fret_eff, transform = gca().transAxes)
    #hist(d.nd[i][mask_green], bins=nbins, color='green', alpha=0.5)
    #hist(d.nd[i][mask_red], bins=nbins, color='red', alpha=0.5)
    xlabel('Burst intensity'); ylabel('# bursts')
def scatter_da_c(i,b,d):
    mask_green, mask_red = mask_burst_green_red(d, i)
    plot(d.nd[i][mask_green],d.na[i][mask_green],'o', mew=0,ms=3, 
            alpha=0.7, color='green')
    plot(d.nd[i][mask_red],d.na[i][mask_red],'o', mew=0,ms=3, 
            alpha=0.7, color='red')
    ngreen = mask_green.sum()
    nred = mask_red.sum()
    #text(0.4,0.9,"#donor = %d, #red= %d" %(ngreen, nred), 
    #        transform = gca().transAxes)
    xlabel('Number of donor ph.'); ylabel('Number of acceptor ph.')
    xlim(-15,250); ylim(-15,130)

def timetrace_merged(i,b,d, bin_width=1e-3, **plot_kw):
    trace, time = binning(d.ph_times_m[i],bin_width_ms=bin_width*1e3)
    plot(clk_to_s(time[1:]), trace, 'b', lw=1.2, alpha=0.8, **plot_kw)
    axhline( d.m/d.T[i]*bin_width, color='b')
    xlabel('Time (s)'); ylabel('# ph')

def timetrace_bars(i,b,d,ax,binw):
    x1,x2 = ax.get_xlim()
    x1ck,x2ck = [x/d.clk_p for x in (x1,x2)]
    mask = (b[:,0] > x1ck)*(b[:,0]<x2ck)
    sel_b = b[mask,:]
    bstart = clk_to_s(sel_b[:,0])
    bwidth = clk_to_s(sel_b[:,1])
    green_val = d.nd[i][mask]/(bwidth/binw)
    red_val = d.na[i][mask]/(bwidth/binw)
    ax.bar(bstart, green_val,bwidth,alpha=0.3,linewidth=0,color='green')
    ax.bar(bstart, -red_val,bwidth,alpha=0.3,linewidth=0,color='red')
    for start,width in zip(bstart,bwidth):
        ax.axvline(start+width*.5,0.75,1, lw=2, color='gray')

def twin_burst_stats(i,b,d,ax1,ax2,stats,binw):
    assert stats.size == b.shape[0]
    ax2.set_ylim(-1,1)
    #ax2.set_ylim(-100,100)
    x1,x2 = ax.get_xlim()
    x1ck,x2ck = [x/d.clk_p for x in (x1,x2)]
    mask = (b[:,0] > x1ck)*(b[:,0]<x2ck)
    sel_b = b[mask,:]
    bstart = clk_to_s(sel_b[:,0])
    ax2.plot(bstart,stats[mask],'d',color='black', mew=0, alpha=0.7)
    ax2.set_xlim(x1,x2)
    draw()


def hist_fret_inset(i,b,d, ifs=9, tpos=0.11, twidth=0.9, plegend=True):
    hist_fret(i,b,d,smoothing=False, clip=False, bins=20, label='CH %d' % (i+1))
    if plegend: 
        leg = legend(fancybox=True, loc='upper right')
        leg.get_frame().set_alpha(0.5)
        leg.get_frame().set_visible(False)
    iax = mai.inset_axes(gca(), 1,1, loc=9)
    timetrace(i,b,d, bin_width=1e-3)
    iax.set_xlabel('s', fontsize=ifs); iax.set_ylabel('')
    iax.set_ylim(-29,29); iax.set_xlim(tpos,tpos+twidth)
    setp(iax.get_yticklabels(), fontsize=ifs)
    setp(iax.get_xticklabels(), fontsize=ifs)

def hist_fret_split(i,b,d,clip=False, bins=None):
    if bins is None: bins = arange(-0.2,1.21,0.05)
    if clip: mask = (d.E[i]>0)*(d.E[i]<=1)
    else: mask = ones(d.E[i].size, dtype=bool)
    mask_green, mask_red = mask_burst_green_red(d, i) 
    hist(d.E[i][mask*mask_green], bins=bins, color='green', alpha=0.7)
    hist(d.E[i][mask*mask_red], bins=bins, color='red', alpha=0.7)
    xlabel('FRET Efficiency'); ylabel('# Bursts')
    ylim(ymin=0);xlim(-0.2,1.2) 

def hist_intensity_split(i,b,d, bins=None, max_size=100):
    if bins is None: bins=arange(d.L,max_size,2)    
    mask_green, mask_red = mask_burst_green_red(d, i) 
    nfr = mask_green.sum()
    fr = mask_red.sum()
    #assert fr+nfr == d.E[i].size
    print "NF: %d  F: %d" % (nfr, fr)
    I1 = sort(b[(mask_green,2)]); I1 = I1[I1<=max_size]
    I2 = sort(b[(mask_red,2)]); I2 = I2[I2<=max_size]
    hist(I1, bins=bins, color='green', alpha=0.5)
    hist(I2, bins=bins, color='red', alpha=0.5)
    xlabel('# Ph.'); ylabel('# Bursts')
    xlim(d.L/2,max_size); ylim(ymin=0)

def plot_mburstm_8ch_super(d, fun=scatter_width_size, sharey=True,
        scroll=False,pgrid=True, figsize=(4,3), name=None, **kwargs):
    if name is None: name = d.fname
    fig, ax = subplots(1,1,figsize=figsize)
    subplots_adjust(left=0.18, right=0.92, top=0.93, bottom=0.18)
    ax.grid(pgrid)
    for i in xrange(d.nch):
        b = d.mburst[i] if hasattr(d, 'mburst') else None   
        if (not fun in [timetrace, ratetrace]) and b.size is 0: continue
        fun(i,b,d,**kwargs)
    ax.autoscale(enable=True,axis='y')
    s=None
    if scroll: s=ScrollingToolQT(fig)
    return AX,s

def plot_mburstm_8ch_vstack(d, fun=scatter_width_size, sharey=True,
        scroll=False,pgrid=True, figsize=(3,9), name=None, **kwargs):
    if name is None: name = d.fname
    fig, AX = subplots(8,1,figsize=figsize, sharex=True, sharey=sharey)
    subplots_adjust(left=0.18, right=0.96, top=0.93, bottom=0.07)
    for i in xrange(d.nch):
        b = d.mburst[i] if hasattr(d, 'mburst') else None   
        if (not fun in [timetrace, ratetrace]) and b.size is 0: continue
        ax = AX.ravel()[i]
        ax.grid(pgrid)
        sca(ax)
        fun(i,b,d,**kwargs)
    setp([a.get_xticklabels() for a in AX[:-1]], visible=False)
    #setp([a.get_yticklabels() for a in AX], visible=False)
    [a.set_xlabel('') for a in AX[:-1]]
    [a.set_ylabel('') for a in AX]
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    if sharey:
        ax.autoscale(enable=True,axis='y')
    s=None
    if scroll: s=ScrollingToolQT(fig)
    return AX,s

def plot_mburstm(d, fun=scatter_width_size, figsize=(12,9), **kwargs):
    p = int(sqrt(len(d.mburst))) # p is 1 or 2    
    fig, AX = subplots(p,p,figsize=figsize, squeeze=False); AX=AX.ravel()
    for i,b in enumerate(d.mburst): 
        if b.size is 0 and fun is not timetrace: continue
        #plot_num = i+1
        #if not len(d.mburst) == 1: AX.append(subplot(2,2,plot_num))
        sca(AX[i])
        if i == 0: suptitle(d.fname+", L = %d, m = %d"%(d.L,d.m))
        title('CH %d, BG = %d cps, T = %1.2f ms, #bursts = %d' %\
                (i+1, d.rate_m[i], d.T[i]*1e3, b.shape[0]), fontsize=12)
        grid(True)
        fun(i,b,d,**kwargs)
    return AX

def plot_mburstm_share(d, fun=scatter_width_size, sharex=True,sharey=True,
        scroll=False,pgrid=True,
        figsize=(12,9), name=None, **kwargs):
    if name is None: name = d.fname
    fig, AX = subplots(2,2,figsize=figsize, sharex=sharex, sharey=sharey)
    for i,b in enumerate(d.mburst): 
        if b.size is 0 and fun is not timetrace: continue
        ax = AX.ravel()[i]

        #ax.axvline(0.8, lw=2, ls='--', color='gray') BIOS-2012
        #ax.axvline(0.85, lw=2, ls='--', color='gray')
        #ax.axvline(0.9, lw=2, ls='--', color='gray')
        #ax.axvline(0.95, lw=2, ls='--', color='gray')
        if i == 0: suptitle(name+", L = %d, m = %d"%(d.L,d.m))
        
        ax.set_title('CH %d, BG = %d cps, T = %1.2f ms, #bursts = %d' %\
                (i+1, d.rate_m[i], d.T[i]*1e3, b.shape[0]), fontsize=12)
        ax.grid(pgrid)
        sca(ax)
        fun(i,b,d,**kwargs)
    if sharex:
        setp([a.get_xticklabels() for a in AX[0,:]], visible=False)
        [a.set_xlabel('') for a in AX[0,:]]
        fig.subplots_adjust(hspace=0.1)
    if sharey:
        setp([a.get_yticklabels() for a in AX[:,1]], visible=False)
        [a.set_ylabel('') for a in AX[:,1]]
        fig.subplots_adjust(wspace=0.1)
        ax.autoscale(enable=True,axis='y')
    s=None
    if scroll: s=ScrollingToolQT(fig)
    return AX,s

def plot_mburstm_8ch_twin(d, fun=scatter_width_size, sharex=True,sharey=True,
        scroll=False, pgrid=True, figsize=(12,9), name=None, **kwargs):
    if name is None: name = d.fname
    fig1, AX1 = subplots(2,2,figsize=figsize, sharex=sharex, sharey=sharey)
    fig2, AX2 = subplots(2,2,figsize=figsize, sharex=sharex, sharey=sharey)
    AX = concatenate([AX1,AX2])
    for i,b in enumerate(d.mburst): 
        if b.size is 0 and fun is not timetrace: continue
        ax = AX.ravel()[i]
        if i==0 or i==4: suptitle(d.status())
        ax.set_title(u'CH%d, BG=%dcps, T=%dμs, #bu=%d' %\
                (i+1, d.rate_m[i], d.T[i]*1e6, b.shape[0]), fontsize=12)
        ax.grid(pgrid)
        sca(ax)
        fun(i,b,d,**kwargs)
    if sharex:
        setp([a.get_xticklabels() for a in AX[0,:]], visible=False)
        [a.set_xlabel('') for a in concatenate([AX1[0,:],AX2[0,:]]).ravel()]
        fig1.subplots_adjust(hspace=0.1); fig2.subplots_adjust(hspace=0.1)
    if sharey:
        setp([a.get_yticklabels() for a in AX[:,1]], visible=False)
        [a.set_ylabel('') for a in concatenate([AX1[:,1],AX[:,1]]).ravel()]
        fig1.subplots_adjust(wspace=0.1); fig2.subplots_adjust(wspace=0.1)
        ax.autoscale(enable=True,axis='y')
    s1,s2=None,None
    if scroll: 
        s1=ScrollingToolQT(fig1)
        s2=ScrollingToolQT(fig2)
    return AX,s1,s2

def plot_tau_fit(hh, AA, Tau):
    import pylab as P
    for i,h in enumerate(hh):
        us = 1e-6
        Tau_us = array(Tau)*clock_period_ns*us
        xd = h[1][:-1]*clock_period_ns*us
        yd = h[0]
        P.semilogy(xd,yd,'--', lw=1.5, color=co[i], alpha=0.3)
        P.semilogy(xd,AA[i]*P.exp(-(xd/Tau_us[i])), 
                lw=2.5, color=co[i], alpha=0.8, 
                label=u"Ch. %d, τ=%1.2f ms, A=%3.1e" % (i+1, Tau_us[i], AA[i]))
    P.xlabel(u'Delay between ph. (μs)')
    P.ylabel('Number of ph.')
    P.legend(); P.grid(True); P.title(fname)

def plot_hist_diff(d, bins_us=1000):
    if d.ALEX:
        hh = [histogram(diff(ph[m]), bins=bins_us*1e-6/d.clk_p)
                for ph,m in zip(d.ph_times_m, d.D_ex)]
    else:
        hh = [histogram(diff(ph), bins=bins_us*1e-6/d.clk_p) 
                for ph in d.ph_times_m]
  
    import pylab as P
    for det,h in enumerate(hh):
        P.semilogy(h[1][:-1]*d.clk_p*1e6,h[0], label='CH %d' %(det+1))
    #P.xlabel(u'Delay between ph. (ms)')
    P.xlabel(ur'Delay between ph. (μs)')
    P.legend(); P.grid(True); P.title(d.fname)


