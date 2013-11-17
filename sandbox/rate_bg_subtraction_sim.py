"""
Try to subtract the BG from the rate histogram using an Eralang distrib.
"""

from scipy.stats import erlang, expon
from scipy.special import gamma, gammainc

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

## Simulate the delays for pure Poisson BG
m = 3
rate0 = 5e3
sim_delays = expon(scale=1./rate0).rvs(1e6)
D = erlang(m, scale=1./rate0).rvs(1e6)
d_exp = expon(scale=1./rate0).rvs(1e6*m)
D2 = d_exp.reshape(1e6, m).sum(axis=-1)

## ECDF for the experimental delays distrib. (BG+Bursts)
D.sort()
if 1:
    # Do some averaging to make numerical derivative smoother
    n_avg = 100; nsize = (D.size/n_avg)*n_avg
    ecdf_x = D[:nsize].reshape(nsize/n_avg,n_avg).mean(-1)
else:
    ecdf_x = D
ecdf_y = arange(ecdf_x.size)/float(ecdf_x.size-1)

## CDF for only BG
_x = linspace(0,D.max(), num=1e3)
bg_pdf = erlang.pdf(_x, m, scale=1./rate0)
bg_cdf = erlang.cdf(_x, m, scale=1./rate0)

# Fix a threshold
# Fit amplitude of CDF

dx = diff(ecdf_x)
epdf_x = ecdf_x[:-1]
epdf_y = (1./ecdf_y.size)/dx

plot(epdf_x, epdf_y)
plot(_x, bg_pdf, lw=3, color='k')

figure()
plot(ecdf_x, ecdf_y)
plot(_x, bg_cdf, lw=3, alpha=0.3)

