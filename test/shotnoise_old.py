from scipy.stats.distributions import binom
from scipy.interpolate import interp1d

#from background import cal_bg_per_burst

## Data() method, temporary position
#    def calc_per_burst_BG(self, relax_bg_t=True):
#        pprint("   - Computing background in eahc burst.\n")
#        bbg_dd, bbg_ad, bbg_t = [0]*self.nch, [0]*self.nch, [0]*self.nch
#        for ich, mb in enumerate(self.mburst):
#            if mb.size == 0: continue
#            width = mb[:,iwidth]*self.clk_p
#            period = self.bp[ich]
#            bbg_dd[ich] = self.bg_dd[ich][period] * width
#            bbg_ad[ich] = self.bg_ad[ich][period] * width
#            if relax_bg_t:
#                bbg_t[ich] = self.bg[ich][period] * width
#            else:
#                bbg_t[ich] = bbg_dd[ich] + bbg_ad[ich]
#            if self.ALEX:
#                bgg_aa[ich] = self.bg_aa[ich][period] * width
#                bgg_t[ich] += bbg_aa[ich]
#        return bbg_dd, bbg_ad, bbg_t

def calc_burst_size_distr(d, ich=0):
    BS = d.nd[ich] + d.na[ich] # total size with BG subtraction
    BSD, bbins = histogram(BS,bins=arange(BS.min(),BS.max()+2))
    BurstSizes = bbins[:-1]

    return BurstSizes, BSD

def shotnoise_hist(d, ich=0, fret_eff=0.8, interp='nearest'):
    """Analytical shotnoise simulation"""
    BurstSizes, BSD = calc_burst_size_distr(d, ich=0)
    step = 0.02
    xx = arange(0,1+step, step) # FRET values axis
    total_hist = zeros(xx.size)
    for num_bursts,burst_size in zip(BSD,BurstSizes):
        Bi = binom(burst_size, fret_eff)
        x = arange(burst_size+2) # axis of possible acceptor counts
        y = Bi.pmf(x)            # prob. of each acceptor count value
        #yy = interp1d(x/float(x.max()-1),y, kind=interp)(xx)
        yy = interp1d(1.*x/x.max(),y, kind=interp)(xx)
        total_hist += num_bursts*(yy)
    return total_hist, xx

def shotnoise_hist_bg_experim(d, ich=0, fret_eff=0.8, interp='nearest'):
    """Analytical shotnoise simulation with background"""
    period = d.bp[ich]
    nd, na, nt, mb = d.nd[ich], d.na[ich], d.nt[ich], d.mburst[ich] 
    BurstSizes, BSD = calc_burst_size_distr(d, ich=0)
    width = b_width(mb)*d.clk_p
    bg_d, bg_a = d.bg_dd[ich][period]*width, d.bg_ad[ich][period]*width
    step = 0.001
    E_ax = arange(0,1+2*step, step)
    total_hist = zeros(E_ax.size)
    da_hist = zeros(E_ax.size)
    for i,bs in enumerate(nt):
        #burst_size = round(bs)
        #Bi = binom(burst_size, fret_eff)
        #a_bi_ax = arange(burst_size+1) # axis of possible acceptor counts
        #P_bi = Bi.pmf(a_bi_ax)         # prob. of each acceptor count value
        
        ## Poisson distributions for current burst background
        Po_d, Po_a = poisson(bg_d[i]), poisson(bg_a[i])
        ## 2D space of possible D and A values due to bg_d and bg_a
        D_bg_ax, A_bg_ax = mgrid[:5*bg_d[i],:5*bg_a[i]]
        ## 2D Probability of each (D,A) coordinate
        P_bg = Po_d.pmf(D_ax) * Po_a.pmf(A_ax)

        da_hist = 0
        for d,a,p in zip(burst_size-a_bi_ax, a_bi_ax, P_bi):
            FRET = (d+D_bg_ax)/(d+D_bg_ax + a+A_bg_ax)
            da_hist += p*histogram(FRET, weights=P_bg, bins=E_ax)
        total_hist += da_hist
    return total_hist, E_ax

def shotnoise_hist_bg_experim2(d, ich=0, fret_eff=0.8, step=0.02):
    """Analytical shotnoise simulation with background.
    """
    if hasattr(d, "bg_corrected") and d.bg_corrected:
        print "ERROR: Data is already BG corrected."
        return
    ## These values (nd,na,nt) must be uncorrected for BG
    nd, na, nt, mb = d.nd[ich], d.na[ich], d.nt[ich], d.mburst[ich] 
    period = d.bp[ich]
    BurstSizes, BSD = calc_burst_size_distr(d, ich=0)
    ## Compute the BG for each burst
    width = b_width(mb)*d.clk_p
    bg_d, bg_a = d.bg_dd[ich][period]*width, d.bg_ad[ich][period]*width
    bg_t = d.bg[ich][period]*width
    ## Compute reasonable max BG count in each burst (must be smaller that nt)
    max_bg_count = bg_t + 5*sqrt(bg_t) + 1
    max_bg_count[max_bg_count > nt] = nt[max_bg_count > nt]
    print 'width', width
    print "per-burst Bg rates (d,a,t):", bg_d, bg_a, bg_t
    ## Setup the histogram
    E_ax = arange(-0.4,1.41, step)+0.5*step
    total_hist = zeros(E_ax.size-1)
    T_bg_ax = arange(max_bg_count.max()+1)
    ## Loop over each burst
    for i,bs in enumerate(nt):
        #if i == 1: continue
        #print "\nBurst %d, nt=%d, nd=%d, na=%d" %(i,bs,nd[i],na[i])

        ## Compute the range of possible BG values inside the burst
        #t_bg_ax = arange(max_bg_count[i]+1)
        t_bg_ax = T_bg_ax[:max_bg_count[i]+1]

        ## Compute the Poisson distribution pmf for the 3 BG (D, A, Total)
        p_bg_t = poisson(bg_t[i]).pmf(t_bg_ax)
        # Restrict axis to values with non-negligible probability
        mask = p_bg_t > 1e-4
        t_bg_ax, p_bg_t = t_bg_ax[mask], p_bg_t[mask]
        #assert p_bg_t.sum() > 0.995
        #print "p_bg_t.sum()",p_bg_t.sum()
        Po_d, Po_a = poisson(bg_d[i]), poisson(bg_a[i])
        bg_ax = arange(t_bg_ax.max()+1) # preallocation
        bg_a_pmf = Po_a.pmf(bg_ax)      # pre-compute on large axis 
        bg_d_pmf = Po_d.pmf(bg_ax)
        
        ## For each total background value (t_bg) in t_bg_ax the sum of donor
        ## and accept. BG is t_bg, so possible values of A BG are
        ## bg_ax[:t_bg+1] and of donor BG are bg_ax[t_bg::-1].
        ## For each donor/accept. BG partition the probability has been
        ## computed in bg_a_pmf and bg_d_pmf.
        ## We compute a vector of possible FRET values for each t_bg (due to
        ## different D/A partion) and a vector of probabilities for each FRET
        ## value.
        FRET = [1.*(na[i]-bg_ax[:t_bg+1])/(nd[i]+na[i]-t_bg)
                for t_bg in t_bg_ax]
        P_FRET = [bg_a_pmf[:t_bg+1]*bg_a_pmf[t_bg::-1]
                for t_bg in t_bg_ax]
        
        ## Alternative implementation with a loop
        #-for t_bg in t_bg_ax:
            #print " - Tot bg", t_bg, p
            #-p_ad_bg = bg_a_pmf[:t_bg+1]*bg_a_pmf[t_bg::-1]
            #print "p_ad_bg.sum()",p_ad_bg.sum()
            #Nd = nd[i]-bg_ax[t_bg::-1]
            #Na = na[i]-bg_ax[:t_bg+1]
            #assert (Nd+Na == bs-t_bg).all()
            #FRET.append(1.*Na/(Nd+Na))
            #-FRET.append(1.*(na[i]-bg_ax[:t_bg+1])/(nd[i]+na[i]-t_bg))
            #-P_FRET.append(p_ad_bg) 
        FRET, P_FRET = concatenate(FRET),concatenate(P_FRET)
        total_hist += histogram(FRET, weights=P_FRET, bins=E_ax)[0]
    #print "Shapes", FRET.shape, P_FRET.shape
    return total_hist, E_ax

def shotnoise_hist_mc(d, ich=0, fret_eff=0.93, N=1e5, bins=None):
    """Monte carlo shotnoise simulation"""
    BurstSizes, BSD = calc_burst_size_distr(d, ich)
    if bins is None: bins = arange(-0.2,1.21,0.05)
    total_hist = zeros(bins.size-1)
    for num_bursts,burst_size in zip(BSD,BurstSizes):
        Bi = binom(burst_size, fret_eff)
        fret_samples = Bi.rvs(N)/float(burst_size)
        h,_ = histogram(fret_samples, bins=bins)
        total_hist += num_bursts*(h/float(N))
    return total_hist, bins

def ShotNoise(d, fret_eff=0.93, N=1e6, bins=None):
    His, Bins = [], []
    for i in range(len(d.ng)):
        #his, bins = shotnoise_hist_mc(d,i,fret_eff=fret_eff,N=N,bins=bins)
        his, bins = shotnoise_hist(d,i,fret_eff=fret_eff)
        His.append(his); Bins.append(bins)
    return His, Bins

def plot_overlay(AX, His, Bins, **kwargs):
    for i,ax in enumerate(AX.ravel()):
        sca(ax)
        plot(Bins[i],His[i],lw=1.6, **kwargs)

def plot_overlay_mc(AX, His, Bins):
    for i,ax in enumerate(AX.ravel()):
        sca(ax); b, h = Bins[i], His[i]
        print b.shape, h.shape
        bar(b[:-1],h,width=b[1]-b[0], alpha=0.6,color='r')

#d2 = Select_bursts(d, select_bursts_fret_value, F=0.95, th=0.0001)
#BS = d.mburst[0][:,2]
#BSD, Bins = histogram(BS,bins=arange(BS.min(),BS.max()+2)); Sizes = Bins[:-1]

#'2011-06-07/DV_8x1_Dicr_xFRET3b_33pM_3.dat'
#d2 = Select_bursts(d, select_bursts_fret, alpha=4)
#H,B=ShotNoise(d2)
#AX,_ = plot_mburstm_share(d2, hist_fret)
#plot_overlay(AX, H,B)

if __name__ == '__main__':
    clk_p = 12.5e-9
    w = 1e-3/clk_p
    ## Make fake data
    mb = [r_[[[0,w,30,0,0,0],[1,w,60,0,0,0],[0,w,30,0,0,0]]]]
    d = Data(mburst=mb,nt=[r_[50,90,50]], na=[r_[17,27,17]], nd=[r_[33,63,33]], 
            bp=[r_[0,0,0]], bg=[r_[1e4]], bg_dd=[r_[3e3]], bg_ad=[r_[7e3]],
            clk_p=clk_p, nch=1, BT=0., gamma=0.)
    
    step = 0.02
    E,E_ax = shotnoise_hist_bg_experim2(d, 0, 0.3, step=step)
    bar(E_ax[:-1],E,step, alpha=0.5); xlim(-0.2,1.2)
