#
# Load Becker & Hikl SPC data file from Ron's ALEX setup.
#
# Run this file from the burst/ folder as follows:
#
# >>> run -i spc/alex_pulsed_bh.py
#

import os
from spcreader import load_spc

ip = get_ipython()
ip.magic("run -i burstlib.py")
ip.magic("run -i burst_plot.py")

donor_ch, accept_ch = 4, 6
clk_period = 50e-9  # (seconds), macrotime discretization period

## ALEX Parameters
D_EX_START = 160    # first donor excitation start in the utime hist (0-4096)
A_EX_START = 660    # first accept excitation start in the utime hist (0-4096)
alt_freq = 1000     # donor-acceptor laser alternating freq (0-4096)
W_size = 400        # window size for donor and acceptor excitation (0-4096)
Nap = 3             # Number of donor-accept alternation to collect

def get_utime_hist(utime, det, mask=None):
    """Returns two 4096-bins histog of utime for donor and accept channels."""
    if mask is None: mask = ones(utime.shape, dtype=bool)
    hd,_ = histogram(utime[(det==donor_ch)*mask],bins=arange(4096))
    ha,_ = histogram(utime[(det==accept_ch)*mask],bins=arange(4096))
    return hd, ha

def get_excit_mask(nanotime):
    ut = nanotime
    donor_ex_mask = zeros(ut.shape, dtype=bool)
    accept_ex_mask = zeros(ut.shape, dtype=bool)
    for i in range(Nap):
        donor_ex_mask += \
                (ut>i*alt_freq+D_EX_START)*(ut<i*alt_freq+D_EX_START+W_size)
        accept_ex_mask += \
                (ut>i*alt_freq+A_EX_START)*(ut<i*alt_freq+A_EX_START+W_size)
    return donor_ex_mask, accept_ex_mask

def plot_range():
    for i in range(Nap):
        donor_ex1 = D_EX_START + i*alt_freq
        donor_ex2 = donor_ex1 + W_size
        accept_ex1 = A_EX_START + i*alt_freq
        accept_ex2 = accept_ex1 + W_size
        axvspan(accept_ex1,accept_ex2,color='red',alpha=0.4)
        axvspan(donor_ex1,donor_ex2,color='green',alpha=0.4)

def check_window_selection(nanotime,det):
    dem,aem = get_excit_mask(nanotime)
    #assert ((det[dem]==donor_ch)+(det[dem]==accept_ch)).sum() == 0
    #assert ((det[aem]!=donor_ch)+(det[aem]!=accept_ch)).sum() == 0

    hd,ha = get_utime_hist(nanotime, det)
    plot(hd,'-g',ha,'-r', alpha=0.3)
    
    hdd,had = get_utime_hist(nanotime,det,dem)
    plot(hdd,'-g',had,'-r')
    
    hda,haa = get_utime_hist(nanotime,det,aem)
    plot(hda,'-g',haa,'-r')
    
    plot_range()

def split_donor_acceptor_ex(d):
    d.dem,d.aem = get_excit_mask(d.tc) # donor and accept excitation masks
    assert sum(-((d.dem!=donor_ch)+(d.dem!=accept_ch))) == 0
    assert sum(-((d.aem!=donor_ch)+(d.aem!=accept_ch))) == 0
    d.ph_times_d = d.ph_times_t[d.dem]
    d.ph_times_a = d.ph_times_t[d.aem]
    d.det_d = d.det_t[d.dem]
    d.det_a = d.det_t[d.aem]

def get_donor_ex(d):
    """Returns a new Data() obj containing ph in the donor excitation."""
    if not hasattr(d,'det_d'): split_donor_acceptor_ex(d)
    D = Data(ph_times=d.ph_times_d,det=d.det_d,ph_times_m=[d.ph_times_d],
            red=[d.det_d==accept_ch],ratem=[0],rates=[0,0],clk_p=d.clk_p,nch=1)
    return D
def get_accept_ex(d):
    """Returns a new Data() obj containing ph in the acceptor excitation."""
    if not hasattr(d,'det_d'): split_donor_acceptor_ex(d)
    D = Data(ph_times=d.ph_times_a,det=d.det_a,ph_times_m=[d.ph_times_a],
            red=[d.det_a==accept_ch],ratem=[0],rates=[0,0],clk_p=d.clk_p,nch=1)
    return D

def bg_raw_rates(d):
    ph_times, red = d.ph_times_m[0], d.red[0]
    rate_d = (-red).sum()/(d.ph_times_m[0][-red].max()*d.clk_p)
    rate_a = (red.sum())/(d.ph_times_m[0][red].max()*d.clk_p)
    return rate_d, rate_a

def bg_rates(d):
    donor_mask, accept_mask = (d.det == donor_ch), (d.det == accept_ch)
    tau_a = fit_tau(d.ph_times[accept_mask], use_linregr=True, clk_p=d.clk_p,
        debug=True)
    tau_d = fit_tau(d.ph_times[donor_mask], use_linregr=True, clk_p=d.clk_p,
        debug=True)
    d.rate_a = array([1/tau_a])
    d.rate_d = array([1/tau_d])
    raw_rate_d, raw_rate_a = bg_raw_rates(d)
    if d.rate_a > raw_rate_a:
        print "WARN: Using raw rate for acceptor."
        d.rate_a = raw_rate_a
    if d.rate_d > raw_rate_d:
        print "WARN: Using raw rate for donor."
        d.rate_d = raw_rate_d
    d.ratem = [d.rate_d+d.rate_a]
    d.rates = [d.rate_d, d.rate_a]

def burst_search_bg(d,L=20,m=3,T_us=1000):
    d.L, d.m, d.T = L, m, [T_us*1e-6]
    bg_rates(d)
    #T_ch_index, T_ch = get_adaptative_T(d.ratem,m,clk_p=clk_p)
    T_ch_index = int(round((T_us*1e-6)/clk_period))
    burst = ba(d.ph_times, L, m, T_ch_index)
    mburst = [burst]
    ng,nr = ph_num_green_red(mburst,d.ph_times_m,d.red)
    d.add(mburst=mburst,ng=ng,nr=nr)
    
def donor_only_analisys(fname, clk_period=clk_period):
    d = load_spc_obj(fname, clk_period=clk_period)
    #Container()
    #d.fname = fname
    #d.clk_p = clk_period 
    #d.ph_times_t, d.det_t, d.tc, dat = load_spc(fname)
    
    check_window_selection(d.tc,d.det)

    # Assign the attributes to use only donor excitation data
    dd = get_donor_ex(d)

    burst_search_bg(dd,L=20,m=3,T_us=1000)
    bleed_through_correct(dd, 0.17)
    background_subtract(dd)

    dd.eff = calculate_fret_eff(dd)
    print print_burst_stats(dd)
    return dd

def load_spc_obj(fname, clk_period=50e-9):
    ph_times_tot, det_tot, nanotime = load_spc(fname)
    d = Data(fname=os.path.basename(fname), clk_p=clk_period, 
            ph_times_tot=ph_times_tot, det_tot=det_tot, nanotime=nanotime)
    return d

def bg_rates_alex(d):
    
    def get_raw_rate(ph_times, ex_mask, ch_mask, clk_p=d.clk_p):
        return ch_mask[ex_mask].sum()/(ph_times[ch_mask*ex_mask].max()*clk_p)
    raw_rate_dd = get_raw_rate(d.ph_times, d.donor_ex, d.donor)
    raw_rate_da = get_raw_rate(d.ph_times, d.donor_ex, -d.donor)
    raw_rate_ad = get_raw_rate(d.ph_times, -d.donor_ex, d.donor)
    raw_rate_aa = get_raw_rate(d.ph_times, -d.donor_ex, -d.donor)

    def get_rate(ph_times):
        return 1./array([fit_tau(ph_times, use_linregr=True, clk_p=d.clk_p, 
            debug=True)])
    d.rate_dd = get_rate(d.ph_times[d.donor*d.donor_ex])
    d.rate_da = get_rate(d.ph_times[(-d.donor)*d.donor_ex])
    d.rate_ad = get_rate(d.ph_times[d.donor*(-d.donor_ex)])
    d.rate_aa = get_rate(d.ph_times[(-d.donor)*(-d.donor_ex)])
    
    if d.rate_dd > raw_rate_dd:
        print "WARN: Using raw rate for donor_ex, donor_ch."
        d.rate_dd = raw_rate_dd
    if d.rate_da > raw_rate_da:
        print "WARN: Using raw rate for donor_ex, accept_ch."
        d.rate_da = raw_rate_da
    if d.rate_ad > raw_rate_ad:
        print "WARN: Using raw rate for accept_ex, donor_ch."
        d.rate_ad = raw_rate_ad
    if d.rate_aa > raw_rate_aa:
        print "WARN: Using raw rate for accept_ex, accept_ch."
        d.rate_aa = raw_rate_aa
    d.rateD = (d.rate_dd+d.rate_da)*0.5
    d.rateA = (d.rate_ad+d.rate_aa)*0.5

def alex_analisys(d):
    #check_window_selection(d.nanotime,d.det_t)

    # Calculate masks to select ph in valid excitation periods
    DEM, AEM = get_excit_mask(d.nanotime)
    EM = DEM+AEM
    ph_times = d.ph_times_tot[EM]
    DEM, AEM = DEM[EM], AEM[EM]
    det = d.det_tot[EM]
    
    # Flags to identify ph: donor or accept. ch; donor or accept. excitation
    DonorCh = (det == donor_ch)
    DonorEx = DEM

    # Safety check and variables cleanups to save RAM
    assert (DEM+AEM).sum() == EM.sum()
    assert DonorCh.shape == DonorEx.shape
    d.delete('det_tot', 'ph_times_tot', 'nanotime'); del det
    d.add(ph_times=ph_times, donor=DonorCh, donor_ex=DonorEx)

    # Calculate the bg rates and store it in d.rate_XX 
    bg_rates_alex(d)
    
def DCBS(d, m=3, L=20):
    """Dual-channel burst search (DCBS) with adaptative threshold.
    
    Ref. for DCBS: Nir et al. 2006. J.Phys.Chem.B.
    """
    get_T = lambda r: get_adaptative_T(r, m=m, clk_p=d.clk_p)[0][0]

    burstD = ba(d.ph_times[d.donor_ex],L,m, get_T(d.rateD)) 
    burstA = ba(d.ph_times[-d.donor_ex],L,m, get_T(d.rateA)) 
    print " DONOR EX.  rate: %d \t T: %d" % (d.rateD, get_T(d.rateD))  
    print " ACCEPT EX. rate: %d \t T: %d" % (d.rateA, get_T(d.rateA))  
    d.add(burstA=burstA, burstD=burstD)

    d.ndd, d.nda = count_burst_donor_accept(burstD, d.ph_times[d.donor_ex])
    d.nad, d.naa = count_burst_donor_accept(burstA, d.ph_times[-d.donor_ex])
    
    # Masks for bursts
    D_burst_num_mask = get_ph_burst_num_mask(d.ph_times, d.donor_ex, burstD)
    A_burst_num_mask = get_ph_burst_num_mask(d.ph_times, -d.donor_ex, burstA)
    
    assert max(D_burst_num_mask) == burstD.shape[0]-1
    assert max(A_burst_num_mask) == burstA.shape[0]-1
    overlap = (D_burst_num_mask>0)*(A_burst_num_mask>0)
    
    Bursts_iD = unique(D_burst_num_mask[overlap])
    Bursts_iA = unique(A_burst_num_mask[overlap])
    # NOTE: Bursts_iD.shape != Bursts_iA
    return  Bursts_iD, Bursts_iA, overlap, D_burst_num_mask,A_burst_num_mask
    

    # TODO: check bg rate calculation (seems quite a bit off)

    #bleed_through_correct(d, 0.17)
    #bg_subtract(d)

    #d.eff = calculate_fret_eff(d)
    #print print_burst_stats(d)
    #return d

def get_ph_burst_mask(ph_times, bursts):
    assert (itstart, iwidth, inum_ph, iistart) == (0,1,2,3)
    mask = zeros(ph_times.size, dtype=bool)
    for b in bursts:
        mask[b[iistart]:b[iistart]+b[inum_ph]] = True
    return mask

def get_ph_burst_num_mask(ph_times, mask, bursts):
    """Return a numeric burst mask for ph in ph_times.
    """ 
    assert (itstart, iwidth, inum_ph, iistart) == (0,1,2,3)
    ph_index = arange(ph_times.size)[mask]
    num_mask = zeros(ph_times.size, dtype=uint32)
    for i,b in enumerate(bursts):
        i_start = ph_index[b[iistart]]
        i_end = ph_index[b[iistart]+b[inum_ph]] 
        num_mask[i_start:i_end] = i
        #if i%100 == 0: pprint("%d %%" % (100*float(i)/bursts.shape[0])) 
    return num_mask

def coincident_bursts(d):
    bursts = []
    ph_index_D = arange(d.ph_times.size)[d.donor_ex]
    ph_index_A = arange(d.ph_times.size)[-d.donor_ex]
    for iD, bD in enumerate(d.burstD):
        btendD = bD[itstart]+bD[iwidth]
        for iA, bA in enumerate(d.burstA):
            btendA = bA[itstart]+bA[iwidth]
                
            btstart, btend = min(bD[itstart], bA[itstart]), max(btendD, btendA)
            if btend > btstart:
                if btendD > bA[itstart] or bD[itstart] < btendA:
                    bwidth = btend - btstart
                    b_num_D, b_num_A = iD, iA
                    b_index_start_D = ph_index_D[bD[iistart]]
                    b_index_start_A = ph_index_A[bA[iistart]]
                    b_index_end_D = ph_index_D[bD[iistart]+bD[inum_ph]]
                    b_index_end_A = ph_index_A[bA[iistart]+bA[inum_ph]]
                    
                    BurstDA = array([btstart, bwidth, b_num_D, b_num_A,
                            b_index_start_D, b_index_start_A,
                            b_index_end_D, b_index_end_A])
                    bursts.append(BurstDA)
        pprint("%d " %iD)
    return array(bursts)


if __name__ == '__main__':
    #fname = '10bp(cy3b+atto647n)(110-57uw)(green+red).spc' #
    #fname = '10bp(cy3b+atto647n)(65-24uw)(green).spc' #
    
    #fname = '20dt+20da(tmr+alexa)(250mm-nacl)(100-40uw)(+bsa)(green).spc' #
    #fname = '20dt+20da(tmr+alexa)(250mm-nacl)(100-40uw)(green+red).spc' #
    
    #fname = '20dt(250mm-nacl)(110-57uw)(green+red).spc'
    
    #fname = '20dt(tmr+alexa)(250mm-nacl)(100-40uw)(+bsa)(green+red).spc' #
    #fname = '20dt(tmr+alexa)(250mm-nacl)(100-40uw)(+bsa)(green).spc' #
    
    #fname = '20dt(tmr+alexa)(65-24uw)(green+red).spc' #
    #fname = '20dt(tmr+alexa)(65-24uw)(green).spc' #
    
    fname = '3bp(cy3b+atto647n)(100-40uw)(green+red).spc'
    #fname = '5bp(cy3b+atto647n)(100-40uw)(green+red).spc'
    #fname = '5bp(cy3b+atto647n)(110-57uw)(green+red).spc'
    
    #fname = '7bp(cy3b+atto647n)(110-57uw)(green+red).spc' #
    #fname = '7bp(cy3b+atto647n)(65-24uw)(green).spc' #

    data_dir = '/home/anto/ucla/burst/data/Ron/spc files/'
    fname = data_dir+fname

    ## DONOR EXCITATION - - - - - - - - - - -
    d = donor_only_analisys(data_dir+fname)

    ## ALEX - - - - - - - - - - - - - - - - - 
    # Load and pre-filter the data
    #d = load_spc_obj(fname)
    #alex_analisys(d)

    # Burst search
    #bd,ba,o,md,ma = DCBS(d,L=30)
    
    #start = 100000
    #end = start + 50000
    #t = arange(start,end)*d.clk_p*1e3
    #plot(t,(1)*(md[start:end]>0),color='green')
    #plot(t,(-1)*(ma[start:end]>0),color='red')
    #plot(t,(0.8)*(o[start:end]>0),color='k', lw=2)
    #xlabel('ms')


    #pprint("\n  >>>>>>>>>>>>>>>>>>>>>>> %s\n" % os.path.basename(fname))
    #d = donor_only_analisys(fname, clk_period=clk_period)
   
    #check_window_selection(d.tc,d.det_t)

