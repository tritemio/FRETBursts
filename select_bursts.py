"""
Functions to select bursts according to different criteria
"""


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  BURSTS SELECTION FUNCTIONS
#

def str_G(gamma, gamma1):
    """Return a string to indicate if and how gamma or gamma1 were used."""
    if gamma1 is None: s = "G%.1f" % gamma
    else: s = "G1_%.1f" % gamma1
    return s

## Selection on E or S values
def E(d, ich=0, E1=-0.2, E2=1.2):
    """Select the burst with E between E1 and E2."""
    burst_mask = (d.E[ich] >= E1)*(d.E[ich] <= E2)
    return burst_mask, ''
def S(d, ich=0, S1=-0.2, S2=1.2):
    """Select the burst with S between S1 and S2."""
    burst_mask = (d.S[ich] >= S1)*(d.S[ich] <= S2)
    return burst_mask, ''
def ES(d, ich=0, E1=-0.2, E2=1.2, S1=-0.2, S2=1.2):
    """Select the burst with E between E1 and E2 and S between S1 and S2."""
    burst_mask = (d.S[ich] >= S1)*(d.S[ich] <= S2) * \
            (d.E[ich] >= E1)*(d.E[ich] <= E2)
    return burst_mask, ''
def ESe(d, ich=0, E1=-0.2, E2=1.2, S1=-0.2, S2=1.2):
    """Select the burst with E-S inside an ellipsis inscribed in E1,E2,S1,S2"""
    def ellips(x,y,x1,x2,y1,y2):
        rx, ry = 0.5*abs(x2-x1), 0.5*abs(y2-y1)
        return ((x-mean([x1,x2]))/rx)**2 + ((y-mean([y1,y2]))/ry)**2
    burst_mask = (ellips(d.E[ich],d.S[ich],E1,E2,S1,S2) <= 1)
    return burst_mask, ''

## Selection on static burst size, width or period
def period(d, ich=0, bp1=0, bp2=None):
    """Select the burst from period bp1 to period bp2 (included)."""
    if bp2 is None: bp2 = d.bp[ich].max()
    burst_mask = (d.bp[ich] >= bp1)*(d.bp[ich] <= bp2)
    return burst_mask, ''
def time(d, ich=0, time_s1=0, time_s2=None):
    """Select the burst starting from time_s1 to time_s2 (in seconds)."""
    burst_start = b_start(d.mburst[ich])*d.clk_p
    if time_s2 is None: time_s2 = burst_start.max()
    burst_mask = (burst_start >= time_s1)*(burst_start <= time_s2)
    return burst_mask, ''

def nd(d, ich=0, th1=20, th2=1000):
    """Select bursts with (nd >= th1) and (nd <= th2)."""
    bursts_mask = (d.nd[ich] >= th1)*(d.nd[ich] <= th2)
    return bursts_mask, ''
def na(d, ich=0, th1=20, th2=1000):
    """Select bursts with (na >= th1) and (na <= th2)."""
    bursts_mask = (d.na[ich] >= th1)*(d.na[ich] <= th2)
    return bursts_mask, ''
def naa(d, ich=0, th1=20, th2=1000):
    """Select bursts with (naa >= th1) and (naa <= th2)."""
    bursts_mask = (d.naa[ich] >= th1)*(d.naa[ich] <= th2)
    return bursts_mask, ''
def nda(d, ich=0, th1=20, th2=1000, gamma=1., gamma1=None,
                      add_naa=False):
    """Select bursts with (nd+na >= th1) and (nd+na <= th2).
    If `gamma` or `gamma1` is specified burst size is computed as:
        nd+na/gamma  (so th1 is the min. burst size for donly bursts)
        nd*gamma1+na (so th1 is the min. burst size for high FRET bursts)
    If data is ALEX and `add_naa` is True, `naa` is added to burst-size.
    """
    if gamma1 is not None:
        burst_size = (1.*d.nd[ich]*gamma1 + d.na[ich])
    else:
        burst_size = (d.nd[ich] + 1.*d.na[ich]/gamma)
    if d.ALEX and add_naa:
        burst_size += d.naa[ich]
    bursts_mask = (burst_size >= th1)*(burst_size <= th2)
    s = "nda_th%d" % th1
    if th2 < 1000: s +="_th2_%d" % th2
    return bursts_mask, s+str_G(gamma, gamma1)
def nda_percentile(d, ich=0, q=50, low=False,
        gamma=1., gamma1=None):
    """Select bursts with SIZE >= q-percentile (or <= if `low` is True)
    `gamma` and `gamma1` are used to compute SIZE like in nda()
    """
    if gamma1 is not None:
        burst_size = (1.*d.nd[ich]*gamma1 + d.na[ich])
    else:
        burst_size = (d.nd[ich] + 1.*d.na[ich]/gamma)
    q_percentile = percentile(burst_size, q=q)
    if low: bursts_mask = (burst_size <= q_percentile)
    else: bursts_mask = (burst_size >= q_percentile)
    return bursts_mask, 'perc%d' % q
def topN_nda(d, ich=0, N=500, gamma=1., gamma1=None):
    """Select the N biggest bursts in the channel.
    `gamma` and `gamma1` are used to compute SIZE like in nda()
    """
    if gamma1 is not None:
        burst_size = (1.*d.nd[ich]*gamma1 + d.na[ich])
    else:
        burst_size = (d.nd[ich] + 1.*d.na[ich]/gamma)
    index_sorted = burst_size.argsort()
    burst_mask = zeros(burst_size.size, dtype=bool)
    burst_mask[index_sorted[-N:]] = True
    return burst_mask, 'topN%d%s' % (N, str_G(gamma,gamma1))

def width(d, ich=0, th1=0.5, th2=1):
    """Select bursts with width between th1 and th2 (ms)."""
    th1, th2 = th1*1e-3/d.clk_p, th2*1e-3/d.clk_p
    burst_width = d.mburst[ich][:,iwidth]
    bursts_mask = (burst_width >= th1)*(burst_width <= th2)
    return bursts_mask, ''


# Selection on burst rate
def max_rate(d, ich=0, min_rate_p=0.1):
    min_rate = d.max_rate[ich].max()*min_rate_p
    mask = (d.max_rate[ich] >= min_rate)
    return mask


## Selection on burst time (nearby, overlapping or isolated bursts)
def single(d, ich=0, th=1):
    """Select bursts that are at least th millisec apart from the others."""
    th = th*1e-3/d.clk_p
    burst_start = d.mburst[ich][:,itstart]
    burst_end = burst_start + d.mburst[ich][:,iwidth]
    gap_mask = (burst_start[1:] - burst_end[:-1]) >= th
    bursts_mask = hstack([gap_mask,False])*hstack([False,gap_mask])
    return bursts_mask
def attached(d, ich=0):
    """Select the first burst of consecutive bursts."""
    burst_mask = (burst_separation(d, ich=ich) <= 0)
    return hstack([burst_mask, (False,)])
def attached2(d, ich=0):
    """Select the second burst of consecutive bursts."""
    burst_mask = (burst_separation(d, ich=ich) <= 0)
    return hstack([(False,), burst_mask])
def nearby(d, ich=0, ms=0.2, clk_p=12.5e-9):
    """Select the first burst of bursts disting less than "ms" millisec."""
    burst_mask = (burst_separation(d, ich=ich) <= (ms*1e-3)/clk_p)
    return hstack([burst_mask, (False,)])
def nearby2(d, ich=0, ms=0.2, clk_p=12.5e-9):
    """Select the second burst of bursts disting less than "ms" millisec."""
    burst_mask = (burst_separation(d, ich=ich) <= (ms*1e-3)/clk_p)
    return hstack([(False,), burst_mask])


## Selection on burst size vs BG
def nd_bg(d, ich=0, F=5):
    """Select bursts with (nd >= bg_dd*F)."""
    bg_burst = d.bg_dd[ich][d.bp[ich]]*b_width(d.mburst[ich])*d.clk_p
    bursts_mask = (d.nd[ich] >= F*bg_burst)
    return bursts_mask
def na_bg(d, ich=0, F=5):
    """Select bursts with (na >= bg_ad*F)."""
    bg_burst = d.bg_ad[ich][d.bp[ich]]*b_width(d.mburst[ich])*d.clk_p
    bursts_mask = (d.na[ich] >= F*bg_burst)
    return bursts_mask
def naa_bg(d, ich=0, F=5):
    """Select bursts with (naa >= bg_aa*F)."""
    bg_burst = d.bg_aa[ich][d.bp[ich]]*b_width(d.mburst[ich])*d.clk_p
    bursts_mask = (d.naa[ich] >= F*bg_burst)
    return bursts_mask
def nt_bg(d, ich=0, F=5):
    """Select bursts with (nt >= bg*F)."""
    bg_burst = d.bg[ich][d.bp[ich]]*b_width(d.mburst[ich])*d.clk_p
    bursts_mask = (d.nt[ich] > F*bg_burst)
    return bursts_mask

## Selection on burst size vs BG (probabilistic)
def na_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ AD signal using P{F*BG>=na} < P."""
    accept_ch_bg_rate = d.rate_ad[ich]
    bursts_width = clk_to_s(d.mburst[ich][:,iwidth])
    max_num_bg_ph = poisson(F*accept_ch_bg_rate*bursts_width).isf(P)
    #print "Min num. ph = ", max_num_bg_ph
    bursts_mask = (d.na[ich] >= max_num_bg_ph)
    return bursts_mask
def nd_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ DD signal using P{F*BG>=nd} < P."""
    donor_ch_bg_rate = d.rate_dd[ich]
    bursts_width = clk_to_s(d.mburst[ich][:,iwidth])
    max_num_bg_ph = poisson(F*donor_ch_bg_rate*bursts_width).isf(P)
    #print "Min num. ph = ", max_num_bg_ph
    bursts_mask = (d.nd[ich] >= max_num_bg_ph)
    return bursts_mask
def naa_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ AA signal using P{F*BG>=naa} < P."""
    A_em_ex_bg_rate = d.rate_aa[ich]
    bursts_width = clk_to_s(d.mburst[ich][:,iwidth])
    max_num_bg_ph = poisson(F*A_em_ex_bg_rate*bursts_width).isf(P)
    #print "Min num. ph = ", max_num_bg_ph
    bursts_mask = (d.naa[ich] >= max_num_bg_ph)
    return bursts_mask
def nt_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ signal using P{F*BG>=nt} < P."""
    bg_rate = d.rate_m[ich]
    bursts_width = clk_to_s(d.mburst[ich][:,iwidth])
    max_num_bg_ph = poisson(F*bg_rate*bursts_width).isf(P)
    #print "Min num. ph = ", max_num_bg_ph
    #print "burst width (ms) = ", bursts_width*1e3
    #print "Poisson rate = ", bg_rate*bursts_width
    #print "rate = ", bg_rate
    bursts_mask = (d.nt[ich] >= max_num_bg_ph)
    return bursts_mask


## Selection on burst skeness
def centered(d, ich=0, th=0.1):
    """Select bursts with absolute value of skewness index less than th."""
    skew_index,_,_ = bleaching(d, ich=ich, exclude_nan=False)
    bursts_mask = (skew_index <= th)*(skew_index >= -th)
    return bursts_mask
def skewness(d, ich=0, th1=0.2, th2=1, **kwargs):
    """Select bursts with skeness between th1 and th2."""
    skew_index,_,_ = bleaching(d, ich=ich, **kwargs)
    bursts_mask = (skew_index <= th2)*(skew_index >= th1)
    return bursts_mask


## Old selection functions
def for_bt_fit(d, ich=0, BT=None):
    """Select bursts for more accurate BT fitting (select before BG corr.)"""
    assert size(BT) == d.nch
    # Selection to be applied as a second-step fit after a first BT fit.
    bursts_mask = (d.na[ich] <= 2*BT[ich]*d.nd[ich])
    return bursts_mask, ''
def size_noise(d, ich=0, th=2):
    """Select bursts w/ size th times above the noise on both D and A ch."""
    burst_width = d.mburst[ich][:,iwidth]
    noise_d, noise_a = burst_width*d.rate_dd[ich], burst_width*d.rate_ad[ich]
    bursts_mask = (d.nd[ich] >= th*noise_d)*(d.na[ich] >= th*noise_a)
    return bursts_mask
def size_noise_or(d, ich=0, th=2):
    """Select bursts w/ size th times above the noise on D or A ch."""
    burst_width = d.mburst[ich][:,iwidth]
    noise_d, noise_a = burst_width*d.rate_dd[ich], burst_width*d.rate_ad[ich]
    bursts_mask = (d.nd[ich] >= th*noise_d)+(d.na[ich] >= th*noise_a)
    return bursts_mask
def no_bg(d, ich=0, P=0.005, NF=1.):
    """Select bursts with prob. to be from BG < P."""
    burst_prob = prob_to_be_bg(d, ich=ich, NF=NF)
    bursts_mask = (burst_prob < P)
    return bursts_mask

def fret_value(d, ich=0, F=0.5, P_th=0.01):
    """Select bursts with prob. > P_th to have fret of F."""
    bsizes = around(d.nd[ich]+d.na[ich]).astype(uint16)
    bursts_mask = zeros(bsizes.size, dtype=bool)
    for burst_size in range(bsizes.min(), bsizes.max()+1):
        indexes = find(bsizes == burst_size)
        RV = binom(burst_size, F)
        #accept_num = arange(burst_size+1)
        #y = RV.cdf(accept_num)
        #min_accept_num = interp(th, y,accept_num)
        min_accept_num = RV.ppf(P_th) # ppf: percent point function (cdf^-1)
        bursts_mask[indexes] = (d.na[ich][indexes] > min_accept_num)
    return bursts_mask



