"""
Try to subtract the BG from the rate histogram using an Eralang distrib.
"""

from scipy.stats import erlang, invgamma
from scipy.special import gamma, gammainc
from scipy.optimize import leastsq
from scipy.signal import gaussian

# Quick functions to calculate rate-trace from ph_times (dense)
ph_rate_d = lambda m,ph: 1.*m/(ph[m-1:]-ph[:ph.size-m+1])     # rate
ph_rate_d_t = lambda m, ph: 0.5*(ph[m-1:]+ph[:ph.size-m+1])   # time for rate

# Quick functions to calculate rate-trace from ph_times (sparse)
ph_rate_s = lambda m, ph: 1.*m/diff(ph[::m])
ph_rate_t_s = lambda m, ph: 0.5*(ph[::m][:-1]+ph[::m][1:])   # time for rate

gauss = lambda M,std: gaussian(M,std)/gaussian(M,std).sum()

def reinterp(x,y, num=1e4):
    """Reinterpolate x,y on a regular grid with `num` points."""
    xi = linspace(x.min(),x.max(),num=num)
    yi = interp(xi, x,y)
    return xi, yi

def deriv(x, y):
    """Derivative of x,y normalized to have area=1."""
    dx = diff(x)
    dy = diff(y)
    return x[:-1], (1./y.size)*dy/dx

def test_pdf(d, ich, m):
    ph, rate0, r, tr, D = get_ph_rate(d, ich, m)
    D.sort()
    ecdf_x = D*1e3
    ecdf_y = arange(ecdf_x.size)/float(ecdf_x.size-1)
    
    x,y = reinterp(ecdf_x,ecdf_y,1e5)
    
    px, py = deriv(x,y)
    #plot(px,py)
    yc = convolve(y, gauss(400,60), mode='same')
    plot(*deriv(x,yc),lw=2)
    #yc = convolve(y, gauss(800,120), mode='same')
    #plot(*deriv(x,yc),lw=2)
    yc = convolve(y, gauss(1600,240), mode='same')
    plot(*deriv(x,yc),lw=2)
    yc = convolve(y, gauss(3200,480), mode='same')
    plot(*deriv(x,yc),lw=2)
    #yc = convolve(y, gauss(6400,960), mode='same')
    #plot(*deriv(x,yc),lw=2)
    ylim(0,py.max()*1.3)

def test_pdf2(d, ich, m, bins=r_[:15:0.01]):
    ph, rate0, r, tr, D = get_ph_rate(d, ich, m)
    H = histogram(D*1e3, bins=bins, normed=True)
    epdf_x = H[1][:-1]
    epdf_y = H[0]
    plot(epdf_x, epdf_y)
    plot(epdf_x, convolve(epdf_y, gauss(15,2),'same'), lw=2)
    plot(epdf_x, 0.8*erlang.pdf(epdf_x, m, scale=1e3/rate0), lw=2, color='k')


if 0:
    ## TEST Erlang
    _rate = 0.5
    _k = 3.
    x = arange(0,20,0.01)
    # Parameters: (delay or waiting time, k or number of events,scale or 1/rate)
    plot(x, erlang.cdf(x, _k, scale=1./_rate))
    plot(x, gammainc(_k, _rate*x), lw=3, alpha=0.3)
    p95 = erlang.ppf(0.95, _k, scale=1./_rate)
    axhline(0.95)
    axvline(p95)

## Assumes a measurement loaded in d
#name = "d.ph_time_m[%d] m=%d" % (ich,m)

def get_ph_rate(d, ich, m, bp=0):
    ph = d.ph_times_m[ich]
    #ph = ph[ph<d.bg_time_s/d.clk_p]
    rate0 = d.bg[ich][bp]
    r = ph_rate_s(m,ph)/clk_p
    tr = ph_rate_t_s(m,ph)*clk_p
    D = 1.*m/r
    return ph, rate0, r, tr, D

def fit_delays_dist(d, ich=0, m=3):
    ph, rate0, r, tr, D = get_ph_rate(d, ich, m)
    
    ## ECDF for the experimental delays distrib. (BG+Bursts)
    D.sort()
    if 0:
        # Do some averaging to make numerical derivative smoother
        n_avg = 20; nsize = (D.size/n_avg)*n_avg
        ecdf_x = D[:nsize].reshape(nsize/n_avg,n_avg).mean(-1)
    else:
        ecdf_x = D
    ecdf_y = arange(ecdf_x.size)/float(ecdf_x.size-1)

    ## CDF for only BG
    _x = linspace(0,D.max(), num=1e5)
    bg_pdf = erlang.pdf(_x, m, scale=1./rate0)
    #bg_cdf = erlang.cdf(_x, m, scale=1./rate0)

    dx = diff(ecdf_x)
    epdf_x = ecdf_x[:-1]
    epdf_y = (1./ecdf_y.size)/dx

    ## Fitting the empirical CDF
    fit_fun = lambda x, a, rate_b: a*erlang.pdf(x, m, scale=1./rate0) +\
            (1-a)*erlang.pdf(x, m, scale=1./rate_b)
    errfunc = lambda p, x, y: fit_fun(x, p[0], p[1]) - y
    p,v = leastsq(errfunc, x0=[0.8, 1e5], args=(epdf_x,epdf_y))

def test_delays_dist(d, ich=0, m=3):
    """Build the distrib. of sum of delays of m ph (Erlang)."""
    ph, rate0, r, tr, D = get_ph_rate(d, ich, m)
    
    ## ECDF for the experimental delays distrib. (BG+Bursts)
    D.sort()
    if 0:
        # Do some averaging to make numerical derivative smoother
        n_avg = 20; nsize = (D.size/n_avg)*n_avg
        ecdf_x = D[:nsize].reshape(nsize/n_avg,n_avg).mean(-1)
    else:
        ecdf_x = D
    ecdf_y = arange(ecdf_x.size)/float(ecdf_x.size-1)

    ## CDF for only BG
    _x = linspace(0,D.max(), num=1e5)
    bg_pdf = erlang.pdf(_x, m, scale=1./rate0)
    bg_cdf = erlang.cdf(_x, m, scale=1./rate0)

    # Fix a threshold
    # Fit amplitude of CDF

    dx = diff(ecdf_x)
    epdf_x = ecdf_x[:-1]
    epdf_y = (1./ecdf_y.size)/dx

    title("Delays dist.")
    plot(epdf_x, epdf_y)
    plot(_x, bg_pdf, lw=3, color='k')

    #figure()
    #plot(ecdf_x, ecdf_y)
    #plot(_x, bg_cdf*0.56, lw=3, alpha=0.3)

def test_rates_dist(d, ich=0, m=3):
    """Build the distrib. of rates of m ph (Inverse Gamma)."""
    ph, rate0, r, tr, D = get_ph_rate(d, ich, m)
    
    ## ECDF for the experimental delays distrib. (BG+Bursts)
    r.sort()
    # Do some averaging to make numerical derivative smoother
    n_avg = 20; nsize = (r.size/n_avg)*n_avg
    ecdf_x = r[:nsize].reshape(nsize/n_avg,n_avg).mean(-1)
    ecdf_y = arange(ecdf_x.size)/float(ecdf_x.size-1)

    ## CDF for only BG
    _x = linspace(0,r.max(), num=1e5)
    # NOTE: invgamma gives the distribution of inverse of delays
    # to get the rate we still have to multiply it by m
    bg_pdf = invgamma.pdf(_x, m, scale=m*rate0)
    bg_cdf = invgamma.cdf(_x, m, scale=m*rate0)

    # Fix a threshold
    # Fit amplitude of CDF

    dx = diff(ecdf_x)
    epdf_x = ecdf_x[:-1]
    epdf_y = (1./ecdf_y.size)/dx

    title("Rates dist. (m=%d)" % m)
    plot(epdf_x, epdf_y)
    plot(_x, bg_pdf, lw=3, color='k')

    #figure()
    #plot(ecdf_x, ecdf_y)
    #plot(_x, bg_cdf*0.56, lw=3, alpha=0.3)
    
