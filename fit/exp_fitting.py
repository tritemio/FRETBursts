import numpy as N
import scipy.optimize as O
import scipy.stats as S
import numpy.random as R
from scipy.special import erf
from scipy.optimize import leastsq


def expon_fit(s, s_min=0):
    """Eponential fit of samples s using MLE.
    All samples < s_min are discarded (s_min must be >= 0).
    Returns the lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0    
    Tau = s.mean() - s_min
    return 1./Tau

def expon_fit_cdf(s, s_min=0):
    """Eponential fit of samples s using a curve fit to the empirical CDF.
    All samples < s_min are discarded (s_min must be >= 0).
    Returns the lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0    
    ## Empirical CDF
    ecdf = [N.sort(s), N.arange(0.5,s.size+0.5)*1./s.size]
    decr_line = N.log(1-ecdf[1])
    L = S.linregress(ecdf[0], decr_line)
    Lambda = -L[0]
    return Lambda


def exp_tail_fit(sample, auto_min_val=0.2, min_val=None, max_val=None,
        end_discard=1., debug=False, plot=False):
    """Returns the eximated tau for the tail of the distribution of sample.
    """
    def plot_tau_fit_single(xd, yd, A, tau, m):
        """Plots xd vs yd and the exponential fit A*exp(-x/tau)"""
        import pylab as P
        l, = P.semilogy(xd, yd, '.', alpha=0.5)
        c = l.get_color()
        xx = P.arange(0, xd.max(), xd.max()/200.)
        P.plot(xx, A*P.exp(-xx/tau), lw=2, alpha=0.4, color=c, 
                label="Tau = %.3f, Rate = %.3f" % (tau, 1./tau))
        P.axvspan(xd[m][0], xd[m][-1], color=c, alpha=0.2)
        P.legend()
        P.xlabel(u'Samples'); P.grid(True)

    ssample = sort(sample)[:-end_discard]

    xd = ssample
    yd = arange(xd.size,0,-1)
    
    if max_val is None: max_val = xd[-1]
    
    if min_val is None:
         min_val = max_val*auto_min_val
         print "AUTO MIN VALUE: min_val %.2f*max_val = %.2f"%\
                 (auto_min_val, min_val)

    # Selection mask for the range to fit
    m = (xd > min_val)*(xd < max_val)
    
    ## Fitting using linear regression 
    x1, y1 = xd[m], yd[m]
    if debug: print " - Fitting on %d points." % x1.size
    assert (x1.size > 2) and (y1.size > 2)
    ly1 = log(y1)
    lr = linregress(x1, ly1)
    tau = -(1.0/lr[0])
    if debug or plot:
        X = [xd[m][0], xd[m][-1]]
        Y = [yd[m][0], yd[m][-1]]
        tau2 = (1.0*X[1]-X[0])/(log(Y[0]/Y[1]))
    if debug:
        print "m: %f y0: %f, r-value: %f, p-value: %f, stderr: %f" % lr
        print " Tau linregress: %.4e  -  Tau 2-points fit: %.4e" %(tau,tau2)
    
    A = exp(lr[1])
    if debug or plot:
        print "A = ", A
        plot_tau_fit_single(xd, yd, A, tau, m)
    
    return tau, A

def exp_tail_fit2(sample, min_val_p=0.1, min_val=None, max_val=None, 
        debug=False, plot=False):
    """Returns the eximated tau for the tail of the distribution of sample.
    """
    def plot_tau_fit_single(xd, yd, A, tau, m):
        """Plots xd vs yd and the exponential fit A*exp(-x/tau)"""
        import pylab as P
        l, = P.semilogy(xd, yd, '.', alpha=0.5)
        c = l.get_color()
        xx = P.arange(0, xd.max(), xd.max()/200.)
        P.plot(xx, A*P.exp(-xx/tau), lw=2, alpha=0.4, color=c, 
                label="Tau = %.3f, Rate = %.3f" % (tau, 1./tau))
        P.axvspan(xd[m][0], xd[m][-1], color=c, alpha=0.2)
        P.legend()
        P.xlabel(u'Samples'); P.grid(True)

    M = 1000
    hh = histogram(sample, bins=arange(M))

    xd = hh[1][:-1]+0.5
    yd = hh[0]
    
    # Exclude the last 10 points
    i = find(cumsum(yd[::-1]) == 10)
    # Exclude the range for which there is 1 or 0 events
    i = find(yd[::-1] == 2)[0]
    if max_val is None: max_val = xd[-i]

    # Min delay adaptively set to a fraction of max_delay
    if min_val is None: min_val = max_val*min_val_p
    print "AUTO MIN VAL: min_val = %.2f*MAX = %.2f"%(min_val_p,max_val)

    # Selection mask for the range to fit
    m = (xd > min_val)*(xd <= max_val)*(yd>0)
    
    ## Fitting using linear regression 
    x1, y1 = xd[m], yd[m]
    
    if debug: print " - Fitting on %d points." % x1.size
    assert (x1.size > 2) and (y1.size > 2)
    ly1 = log(y1)
    lr = linregress(x1, ly1)
    tau = -(1.0/lr[0])
    if debug or plot:
        X = [xd[m][0], xd[m][-1]]
        Y = [yd[m][0], yd[m][-1]]
        tau2 = (1.0*X[1]-X[0])/(log(Y[0]/Y[1]))
    if debug:
        print "m: %f y0: %f, r-value: %f, p-value: %f, stderr: %f" % lr
        print " Tau linregress: %.4e  -  Tau 2-points fit: %.4e" %(tau,tau2)
    
    A = exp(lr[1])
    if debug or plot:
        print "A = ", A
        plot_tau_fit_single(xd, yd, A, tau, m)
    
    return tau, A

def exp_tail_fit3(sample, min_val_p=0.1, min_val=None, max_val=None, 
        debug=False, plot=False, end_discard=1):
    """Returns the eximated tau for the tail of the distribution of sample.
    """
    def plot_tau_fit_single(xd, yd, A, tau, min_val):
        """Plots xd vs yd and the exponential fit A*exp(-x/tau)"""
        import pylab as P
        l, = P.semilogy(xd, yd, '.', alpha=0.5)
        c = l.get_color()
        xx = P.arange(0, xd.max(), xd.max()/200.)
        P.plot(xx, A*P.exp(-xx/tau), lw=2, alpha=0.4, color=c, 
                label="Tau = %.3f, Rate = %.3f" % (tau, 1./tau))
        P.axvspan(min_val, 1000, color=c, alpha=0.2)
        P.legend()
        P.xlabel(u'Samples'); P.grid(True)

    for i in range(end_discard):
        sample = sample[sample<sample.max()]
        
    if max_val is None:
        ssample = sort(sample)
        max_val = ssample[-1]
    
    if min_val is None: min_val = max_val*min_val_p

    m = (sample>min_val)
    
    tau = sample[m].mean() - min_val
    step = 1
    A = step*m.sum()/tau * exp(min_val/tau)

    if debug or plot:
        print "A = ", A
        hh = histogram(sample, bins=arange(0,500,step))
        xd = hh[1][:-1]
        yd = hh[0]
        plot_tau_fit_single(xd, yd, A, tau, min_val)
    
    return tau, A

