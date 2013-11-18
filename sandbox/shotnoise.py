"""
Compute shot-noise FRET peak broadening.
"""


from scipy.stats.distributions import binom
from scipy.interpolate import interp1d
import numexpr as NE

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
        yy = interp1d(1.*x/x.max(),y, kind=interp)(xx)
        #yy = interp1d(1.*x/burst_size,y, kind=interp)(xx) # UNTESTED BUT
        #                                                  # SHOULD BE CORRECT!
        total_hist += num_bursts*(yy)
    return total_hist, xx

def shotnoise_hist_bg(d, ich=0, step=0.02):
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
    Bg_d, Bg_a = d.bg_dd[ich][period]*width, d.bg_ad[ich][period]*width
    Bg_t = d.bg[ich][period]*width
    ## Compute reasonable max BG count in each burst (must be smaller that nt)
    Max_bg_count = ceil(Bg_t + 5*sqrt(Bg_t))
    Max_bg_count[Max_bg_count > nt] = nt[Max_bg_count > nt]
    print 'width', width
    print "per-burst Bg rates (d,a,t):", Bg_d, Bg_a, Bg_t
    ## Setup the histogram
    E_ax = arange(-0.4,1.41, step)+0.5*step
    total_hist = zeros(E_ax.size-1)
    T_bg_ax = arange(Max_bg_count.max()+1)
    ## Loop over each burst
    for i,bs in enumerate(nt):
        #if i == 1: continue
        #print "\nBurst %d, nt=%d, nd=%d, na=%d" %(i,bs,nd[i],na[i])

        ## Compute the range of possible BG values inside the burst
        #t_bg_ax = arange(max_bg_count[i]+1)
        bg_max = Max_bg_count[i]
        bg_t_ax = T_bg_ax[:bg_max+1]
        
        ## Compute the Poisson distribution pmf for the 3 BG (D, A, Total)
        bg_t_p = poisson(Bg_t[i]).pmf(bg_t_ax)
        
        # Restrict axis to values with non-negligible probability
        # (disable for array implementation)
        #mask = bg_t_p > 1e-4
        #bg_t_ax, bg_t_p = bg_t_ax[mask], bg_t_p[mask]
        
        #assert bg_t_p.sum() > 0.995
        #print "bg_t_p.sum()",bg_t_p.sum()
        poiss_d, poiss_a = poisson(Bg_d[i]), poisson(Bg_a[i])
        
        bg_ax = arange(bg_t_ax.max()+1)    # preallocation
        # pre-compute on large axis 
        bg_a_pmf = poiss_a.pmf(bg_ax)
        bg_d_pmf = poiss_d.pmf(bg_ax)      
        bg_d_pmfi = bg_d_pmf[::-1]      
        
        ## For each total background value (bg_t) in bg_t_ax the sum of donor
        ## and accept. BG is bg_t, so possible values of A BG are
        ## bg_ax[:bg_t+1] and of donor BG are bg_ax[bg_t::-1].
        ## For each donor/accept. BG partition the probability has been
        ## computed in bg_a_pmf and bg_d_pmf.
        ## We compute a vector of possible FRET values for each t_bg (due to
        ## different D/A partion) and a vector of probabilities for each FRET
        ## value.

        # FRET with no gamma        
        #FRET = [1.*(na[i]-bg_ax[:bg_t+1])/(nd[i]+na[i]-bg_t)
        #       for bg_t in bg_t_ax]
        #P_FRET = [bg_a_pmf[:bg_t+1]*bg_d_pmf[bg_t::-1]
        #        for bg_t in bg_t_ax]
 
        # FRET with gamma (loop implementation, list)
        #FRET1 = [[]]*bg_t_ax.size
        #P_FRET1 = [[]]*bg_t_ax.size
        #for i1, bg_t in enumerate(bg_t_ax):
        #    FRET1[i1] = 1.*(na[i]-bg_ax[:bg_t+1])/(na[i]-bg_ax[:bg_t+1] + \
        #            1.*(nd[i]-bg_ax[bg_t::-1]))
        #    P_FRET1[i1] = bg_a_pmf[:bg_t+1]*bg_d_pmf[bg_t::-1]
        #for ix in range(bg_t_ax.size):
        #    assert (FRET[ix] == FRET1[ix]).all()
        #    assert (P_FRET[ix] == P_FRET1[ix]).all()

        # FRET with gamma (more explicit loop implementation, list)
        #FRET1 = [[]]*bg_t_ax.size
        #P_FRET1 = [[]]*bg_t_ax.size
        #for i1, bg_t in enumerate(bg_t_ax):
        #    bg_a_ax = bg_ax[:bg_t+1]
        #    bg_a_p = bg_a_pmf[:bg_t+1]
        #    bg_d_ax = bg_ax[bg_t::-1]
        #    bg_d_p = bg_d_pmf[bg_t::-1]
        #    assert ((bg_a_ax + bg_d_ax) == bg_t).all()
        #    FRET1[i1] = 1.*(na[i]-bg_a_ax)/(na[i]-bg_a_ax+ 1.*(nd[i]-bg_d_ax)) 
        #    P_FRET1[i1] = bg_a_p*bg_d_p
        #for ix in range(bg_t_ax.size):
        #    assert (FRET[ix] == FRET1[ix]).all()
        #    assert (P_FRET[ix] == P_FRET1[ix]).all()
        
        # FRET with gamma (loop implementation, array)
        FRET1 = zeros((bg_ax.size,bg_ax.size))
        P_FRET1= zeros((bg_ax.size,bg_ax.size))
        bg_d_p = zeros(bg_ax.size)
        nd_i, na_i = nd[i], na[i]
        for i1, bg_t in enumerate(bg_t_ax):
            bg_a_ax = bg_ax
            bg_a_p = bg_a_pmf.copy()
            bg_a_p[bg_t+1:] = 0
            bg_d_ax = bg_ax[::-1] - (bg_max-bg_t)
            #bg_d_p = poiss_d.pmf(bg_d_ax)
            bg_d_p[:bg_t+1] = bg_d_pmfi[bg_max-bg_t:]
            #assert (bg_d_p1 == bg_d_p).all()
            #assert ((bg_a_ax + bg_d_ax) == bg_t).all()
            FRET1[i1] = 1.*(na_i-bg_a_ax)/(na_i-bg_a_ax + 1.*(nd_i-bg_d_ax)) 
            P_FRET1[i1] = bg_a_p*bg_d_p
        #for ix in range(bg_t_ax.size):
        #    m = (P_FRET1[ix] > 0)
        #    # When bg_t == nt (even with negligible probability) then
        #    # denominator of FRET is 0 and if na == bg_a_ax then we habe a 0/0
        #    # that result in nan. Two nan are never equal, so I set them to 99.
        #    FRET[ix][isnan(FRET[ix])] = 99 
        #    FRET1[ix][isnan(FRET1[ix])] = 99 
        #    assert (FRET[ix] == FRET1[ix][m]).all()
        #    assert (P_FRET[ix] == P_FRET1[ix][m]).all()

        # FRET with gamma
        #FRET = [1.*(na[i]-bg_ax[:bg_t+1])/(na[i]-bg_ax[:bg_t+1] + \
        #        1.*(nd[i]-bg_ax[bg_t::-1])) for bg_t in bg_t_ax]
        #for ix in range(bg_t_ax.size):
        #    assert (FRET[ix] == FRET1[ix]).all()
        
        FRET, P_FRET = FRET1, P_FRET1
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
    d_ = Data(mburst=mb,nt=[r_[50,90,50]],na=[r_[17,27,17]],nd=[r_[33,63,33]], 
            bp=[r_[0,0,0]], bg=[r_[1e4]], bg_dd=[r_[3e3]], bg_ad=[r_[7e3]],
            clk_p=clk_p, nch=1, BT=0., gamma=0.)
    
    step = 0.02
    E,E_ax = shotnoise_hist_bg(d_, 0, step=step)
    #bar(E_ax[:-1],E,step, alpha=0.5); xlim(-0.2,1.2)
