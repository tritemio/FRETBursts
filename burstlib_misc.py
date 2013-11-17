"""
Misc functions moved from burstlib.py. These functions may be useful in future
or may be deleted.
"""

def fret_stats(d, th=None):
    FR, NFR, EFF = [], [], []
    for ich in range(d.nch):
        mask_green, mask_red = mask_burst_green_red(d,ich,th=th)
        NFR.append(mask_green.sum())
        FR.append(mask_red.sum())
        EFF.append(fret_efficiency(d, ich))
    return FR, NFR, EFF

def print_fret_stats(d, th=None):
    FR, NFR, EFF = fret_stats(d, th=th)
    nch = len(d.mburst)
    s = "\n\nFRET STATISTICS"
    s += "\nPixel:           "+"%7d "*nch % tuple(range(1,nch+1))
    s += "\n#burst w/  FRET: "+"%7d "*nch % tuple(FR)
    s += "\n#burst w/o FRET: "+"%7d "*nch % tuple(NFR)
    s += "\n#FRET eff.(don): "+"%7.2f "*nch % tuple(EFF)
    return s

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  BURST ANALISYS  
#
    
def bleaching1(d, ich, ib, use_median=True, normalize=True):
    """Compute bleaching of a single burst."""
    if use_median: fun = median
    else: fun = mean
    b = d.mburst[ich][ib]
    burst_slice = slice(b_istart(atleast_2d(b)), b_iend(atleast_2d(b))+1)
    ph = d.ph_times_m[ich][burst_slice] 
    a_mask = d.A_em[ich][burst_slice]
    d_mask = -a_mask
    md = fun(ph[d_mask].astype(float))
    ma = fun(ph[a_mask].astype(float))
    mt = fun(ph.astype(float))
    if normalize: skew = (ma-md)/float(b_width(atleast_2d(b)))
    else: skew = 1e3*clk_to_s(ma-md) # milli-seconds
    return skew

def bleaching(d, ich=0, exclude_nan=True, use_median=True, normalize=True):
    if use_median: fun = median
    else: fun = mean
    mg, mr, mt = stat_burst(d, ich, fun)
    #sg, sr, st = stat_burst(d, ich, std)
    if normalize:
        bleach = (mr-mg)/d.mburst[ich][:,1]
        s1='(%s(red)-%s(green))/width' % (fun.__name__, fun.__name__)
    else:
        bleach = 1e3*clk_to_s(mr-mg)
        s1='(%s(red)-%s(green)) [ms]' % (fun.__name__, fun.__name__)
    #bleach = (mr-mg)/st; s1='(mean(red)-mean(green))/std(width)'
    if isnan(bleach).any(): print "WARN: NaN in bleaching function."
    if exclude_nan: bleach = bleach[-isnan(bleach)]
    bleach_val = (sum(bleach<0)-sum(bleach>0))/float(bleach.size)
    s2 = "CH %d, #burst = %5d/%5d, Bleaching = %6.3f %%"  %\
            (ich+1, bleach.size, d.mburst[ich].shape[0], bleach_val*100)
    print s2
    return bleach, s1, s2

def Bleaching(d, **kwargs):
    BL = []
    for ich in range(d.nch):
        BL.append(bleaching(d, ich, **kwargs)[0])
    return BL

def asymmetry(d, ich=0):
    mg, mr, mt = stat_burst(d, ich, median)
    bsize = b_size(d.mburst[ich])
    bwidth = b_width(d.mburst[ich])
    asym = (mr-mg)*bsize/bwidth; s='(m(red)-m(green))*bsize/bwidth'
    return asym, s

def prob_to_be_bg(d, ich=0, NF=1., clk_p=12.5e-9):
    """Returns the prob. that a burst is due to bg, for each burst in ich."""
    burst_width = d.mburst[ich][:,iwidth]
    burst_size = d.mburst[ich][:,inum_ph]
    bg_rate = d.ratem[ich]
    #prob = SS.poisson(bg_rate*(burst_width*clk_p)).sf(d.L-1)
    prob = SS.poisson(bg_rate*NF*(burst_width*clk_p)).sf(burst_size-1)
    return array(prob)
def Prob_to_be_bg(d, **kwargs):
    Prob = []
    for ich in range(d.nch):
        Prob.append(prob_to_be_bg(d, ich, **kwargs))
    return Prob

## 
# MISC FUCNTIONS
#

def binning(times, bin_width_ms=1, max_num_bins=1e5, clk_p=12.5e-9):
    """Return the binned histogram of array times."""
    bin_width_clk = (bin_width_ms*1e-3)/clk_p
    num_bins = min(times.max()/bin_width_clk,max_num_bins)
    h = histogram(times[times<(num_bins*bin_width_clk)], bins=num_bins) 
    return h

def mbinning(times,det, bin_width_ms=1, max_num_bins=1e5, clk_p=12.5*1e9):
    """Return the binned arrival times for all channels."""
    num_det = det.max()
    H = [binning(times[det==d], bin_width_ms, max_num_bins, clk_p)
            for d in range(1,num_det+1)]
    return H

# Maybe will delete these since they are of little use:
def burst_start(d, ich=0):
    return d.mburst[ich][:,itstart]
def burst_end(d, ich=0):
    return d.mburst[ich][:,itstart] + d.mburst[ich][:,iwidth]
def burst_width(d, ich=0):
    return d.mburst[ich][:,iwidth]
def burst_size(d, ich=0):
    return d.mburst[ich][:,inum_ph]
def burst_separation(d, ich=0):
    return burst_start(d,ich)[1:]-burst_end(d,ich)[:-1]


