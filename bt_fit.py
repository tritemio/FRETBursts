from scipy.stats import poisson, binom, linregress
from scipy.misc import factorial
from scipy.optimize import minimize_scalar, leastsq
from scipy.special import gammaln
from numpy import r_, log, round

from burstsearch.bs import b_width

#from burst_selection import Sel, select_bursts_nda
# TODO: Solve the problem of imports without using run -i

def likelihood_bruteforce(Kbt, nd, na, bg_a):
    """-Likelihood function for BT. Computed brute-force from poisson PMF."""
    P = 1.
    for ndi, nai, bg_di, bg_ai in zip(nd,na,bg_d,bg_a):
        lam_ai = bg_ai + Kbt*(ndi-bg_di)
        P *= poisson(lam_ai).pmf(nai)
    return -P

def likelihood_raw_log(Kbt, nd, na, bg_d, bg_a):
    """-Likelihood function for BT. Computed from raw log-likelihood."""
    logP = 0
    for ndi, nai, bg_di, bg_ai in zip(nd,na,bg_d,bg_a):
        lam_ai = bg_ai + Kbt*(ndi-bg_di)
        #assert log(factorial(nai)) == gammaln(nai+1)
        logP += log(lam_ai**nai) - gammaln(nai+1) - lam_ai
    P = exp(logP)
    return -P

def loglikelihood_func_raw(Kbt, nd, na, bg_d, bg_a):
    """True -Log-Likelihood function for BT."""
    logP = 0
    for ndi, nai, bg_di, bg_ai in zip(nd,na,bg_d,bg_a):
        lam_ai = bg_ai + Kbt*(ndi-bg_di)
        logP += log(lam_ai**nai) - gammaln(nai+1) - lam_ai
    return -logP

def loglikelihood_func_opt(Kbt, nd, na, bg_d, bg_a):
    """Opt. -Log-Likelihood func. for BT. Faster to compute but w/ same min.
    """
    logP = 0
    for ndi, nai, bg_di, bg_ai in zip(nd,na,bg_d,bg_a):
        lam_ai = bg_ai + Kbt*(ndi-bg_di)
        #logP += log(lam_ai**nai) - lam_ai
        logP += nai*log(lam_ai) - lam_ai
    #logP -= (bg_a+Kbt*(nd-bg_d)).sum()
    return -logP

def LogLikelihood_binom(E, nd, na):
    """Likelihood function for (nd,na) to be from a binom with p=E (no BG)."""
    return -(log(binom.pmf(na, na+nd, E))).sum()

def fit_ML_s(nd, na, bg_d, bg_a):
    """Fit the BT using MLE method."""
    res = minimize_scalar(loglikelihood_func_opt, bounds=(1e-4,0.2),
            method='bounded', args=(nd,na,bg_d,bg_a), 
            options={'disp':0, 'xtol': 1e-8})
    print res
    return res.x
def fit_LS_s(nd, na, bg_d, bg_a):
    """Fit the BT using LS regression to a line (w/ intercept=0)."""
    err_fun = lambda kb, nd, na, bg_d, bg_a: (nd-bg_d)*kb - (na-bg_a)
    res = leastsq(err_fun, x0=0.05, args=(nd,na,bg_d,bg_a), full_output=0,
            xtol=1e-8)
    return res[0][0]
def fit_LR(nd, na, bg_d, bg_a, ret_intercept=False):
    """Fit the BT using linear regression (NOTE: intercept is not set to 0)."""
    res = linregress(x=(nd-bg_d),y=(na-bg_a))
    slope, intercept = res[0], res[1]
    print "  - Intercept: ", intercept
    if ret_intercept: return slope, intercept
    else: return slope

def fit_ML_binom_alt(d, ich):
    """Fit the BT using MLE method with a binomial distribution.
    This version fits the BG-corrected (nd,na), trying to avoid negative values
    """
    assert not (d.bg_corrected*d.bt_corrected)
    
    # NOTE: Values are with no correction here
    nd, na, bg_d, bg_a = expand(d,ich)

    nd, na = round(nd-bg_d).astype(int), round(na-bg_a).astype(int)
    
    # The binomial distribution can not handle negative values
    # so we must find a way to "remove" the negativa values
    # a. remove bursts with neg. values, but we can skew the fit
    # b. remove bursts with nd < nd[na<0].max(), but few bursts left
    # c. remove bursts with na+nd < nd[na<0].max()
    pos_bursts = (nd>=0)*(na>=0)
    if (-pos_bursts).any():
        # a. remove bursts with neg. values
        #nd, na = nd[pos_bursts], na[pos_bursts]
        
        # b. Cut all the part with na<0 to have a less skewed estimation
        nd_min = nd[na<0].max()
        #nd, na = nd[nd>nd_min], na[nd>nd_min]
        
        # c. remove bursts with na+nd < nd[na<0].max()
        nd, na = nd[nd+na>nd_min], na[nd+na>nd_min]
        
    # Note a max p=0.5 means a max k of p/(1-p)=1 (100%)
    res = minimize_scalar(LogLikelihood_binom, bounds=(1e-6,0.5),
            method='bounded', args=(nd,na), options={'disp':0})
    p = res.x
    k = p/(1-p)
    return k

def _fit_ML_binom(d, ich):
    """Fit the BT using MLE method with a binomial distribution.
    It fits the raw nd, na values, the value must be corrected after.
    """
    assert not (d.bg_corrected*d.bt_corrected)
    # Note a max p=0.5 means a max k of p/(1-p)=1 (100%)
    res = minimize_scalar(LogLikelihood_binom, bounds=(1e-6,0.5),
            method='bounded', args=expand(d,ich)[:2], options={'disp':0})
    p = res.x
    k = p/(1-p)
    return k

def fit_ML_binom(d, ich):
    """Fit the BT using MLE method with a binomial distribution.
    It fits the raw nd, na values, then corrects the value.
    """
    nd,na,bg_d,bg_a,w = expandw(d,ich)
    kbt = _fit_ML_binom(d,ich)    
    kd = fit_slope(nd,bg_d)
    ka = fit_slope(na,bg_d)

    return  kbt*(1-ka)/(1-kd)

def fit_slope(x, y):
    """Fit the slope using Least-Square regression to a line (w/intercept=0)."""
    err_fun = lambda k, x, y: y - k*x
    res = leastsq(err_fun, x0=1, args=(x,y), full_output=0)
    return res[0][0]

def fit_ML(d, ich):
    """Fit the BT using MLE method."""
    assert not (d.bg_corrected*d.bt_corrected)
    #ds = Sel(d, select_bursts_nda, th1=th, nofret=1)
    res = minimize_scalar(loglikelihood_func_opt, bounds=(1e-4,0.2),
            method='bounded', args=expand(d,ich), options={'disp':0})
    return res.x
def fit_LS(d, ich):
    """Fit the BT using Least-Square regression to a line (w/ intercept=0)."""
    assert not (d.bg_corrected*d.bt_corrected)
    err_fun = lambda kb, nd, na, bg_d, bg_a: (nd-bg_d)*kb - (na-bg_a)
    #ds = Sel(d, select_bursts_nda, th1=th, nofret=1)
    res = leastsq(err_fun, x0=0.05, args=expand(d,ich), full_output=0)
    return res[0][0]
def fit_LR(d, ich, ret_intercept=False):
    """Fit the BT using linear regression (NOTE: intercept is not set to 0)."""
    assert not (d.bg_corrected*d.bt_corrected)
    nd, na, bg_d, bg_a = expand(d, ich)
    res = linregress(x=(nd-bg_d),y=(na-bg_a))
    slope, intercept = res[0], res[1]
    print "  - Intercept: ", intercept
    if ret_intercept: return slope, intercept
    else: return slope

def expand(d, ich):
    period = d.bp[ich]
    w = b_width(d.mburst[ich])*d.clk_p
    bg_a = d.bg_ad[ich][period]*w
    bg_d = d.bg_dd[ich][period]*w
    return d.nd[ich], d.na[ich], bg_d, bg_a
def expandw(d, ich):
    period = d.bp[ich]
    w = b_width(d.mburst[ich])*d.clk_p
    bg_a = d.bg_ad[ich][period]*w
    bg_d = d.bg_dd[ich][period]*w
    return d.nd[ich], d.na[ich], bg_d, bg_a, w


## NOTE: Plot function moved to burst_plot

if __name__ == '__main__':
    if 1:
        dir = '/home/anto/Documents/ucla/src/burst/data/'
        fname = dir+'2013-04-05/do_320mW_xd_pix2-6_0y_3.dat'

        clk_p = 12.5e-9
        d = Data(fname=fname, clk_p=clk_p, nch=8, BT=1., gamma=1.)
        d.load_multispot_cache(bytes_to_read=-1, swap_D_A=1)
        d.calc_bg_cache(bg_calc_exp, time_s=5, tail_min_p=0.1)
    d.burst_search_t(L=10,m=10,P=None,F=6, ph_sel='D', nofret=True)
    d.cal_ph_num()

    #compare_bt_fit(d, th=60, ich=0)
    
    for ich in range(d.nch):
        figure()
        test_LS_bt_fit(d, ich=ich, Th=r_[10:250:5])
        test_ML_bt_fit(d, ich=ich, Th=r_[10:250:5])
        legend(loc='best')
        title('CH%d: %s'%(ich+1,d.name())); grid(1)
    

    
