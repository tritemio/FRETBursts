#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module provides functions to fit gaussian distributions and gaussian
distribution mixtures (2 components). These functions can be used directly,
or more often, in a typical FRETBursts workflow they are passed to higher
level methods like :meth:`fretbursts.burstlib.Data.fit_E_generic`.

Single Gaussian distribution fit:

    * :func:`gaussian_fit_hist`
    * :func:`gaussian_fit_cdf`
    * :func:`gaussian_fit_pdf`

For 2-Gaussians fit we have the following models:

    * :func:`two_gauss_mix_pdf`: *PDF of 2-components Gaussians mixture*
    * :func:`two_gauss_mix_ab`: *linear combination of 2 Gaussians*

Main functions for mixture of 2 Gaussian distribution fit:

    * :func:`two_gaussian_fit_hist` *histogram fit using `leastsq`*
    * :func:`two_gaussian_fit_hist_min` *histogram fit using `minimize`*
    * :func:`two_gaussian_fit_hist_min_ab` *the same but using _ab model*
    * :func:`two_gaussian_fit_cdf` *curve fit of the CDF*
    * :func:`two_gaussian_fit_EM` *Expectation-Maximization fit*
    * :func:`two_gaussian_fit_EM_b` *the same with boundaries*

Also, some functions to fit 2-D gaussian distributions and mixtures are
implemented but not thoroughly tested.

The reference documentation for **all** the functions follows.
"""

from __future__ import division, print_function, absolute_import

import numpy as np
import numpy.random as R
import scipy.optimize as O
import scipy.stats as S
from scipy.special import erf
from scipy.optimize import leastsq, minimize
import scipy.ndimage as ndi

#from scipy.stats import gaussian_kde
from .weighted_kde import gaussian_kde_w  # this version supports weights


def normpdf(x, mu=0, sigma=1.):
    """Return the normal pdf evaluated at `x`."""
    assert sigma > 0
    u = (x-mu)/sigma
    y = 1/(np.sqrt(2*np.pi)*sigma)*np.exp(-u*u/2)
    return y

##
# Single gaussian distribution
#

def gaussian_fit_curve(x, y, mu0=0, sigma0=1, a0=None, return_all=False,
        **kwargs):
    """Gaussian fit of curve (x,y).
    If a0 is None then only (mu,sigma) are fitted (to a gaussian density).
    `kwargs` are passed to the leastsq() function.

    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq is returned.
    """
    if a0 is None:
        gauss_pdf = lambda x, m, s: np.exp(-((x-m)**2)/(2*s**2))/\
                                    (np.sqrt(2*np.pi)*s)
        err_fun = lambda p, x, y: gauss_pdf(x, *p) - y
        res = leastsq(err_fun, x0=[mu0, sigma0], args=(x, y), **kwargs)
    else:
        gauss_fun = lambda x, m, s, a: a*np.sign(s)*np.exp(-((x-m)**2)/(2*s**2))
        err_fun = lambda p, x, y: gauss_fun(x, *p) - y
        res = leastsq(err_fun, x0=[mu0, sigma0, a0], args=(x, y), **kwargs)

    if 'full_output' in kwargs:
        return_all = kwargs['full_output']
    mu, sigma = res[0][0], res[0][1]
    if return_all: return res
    return mu, sigma

def get_epdf(s, smooth=0, N=1000, smooth_pdf=False, smooth_cdf=True):
    """Compute the empirical PDF of the sample `s`.

    If smooth > 0 then apply a gaussian filter with sigma=smooth.
    N is the number of points for interpolation of the CDF on a uniform range.
    """
    ecdf = [np.sort(s), np.arange(0.5, s.size+0.5)*1./s.size]
    #ecdf = [np.sort(s), np.arange(s.size)*1./s.size]
    _x = np.linspace(s.min(), s.max(), N)
    ecdfi = [_x, np.interp(_x, ecdf[0], ecdf[1])]
    if smooth_cdf and smooth > 0:
        ecdfi[1] = ndi.filters.gaussian_filter1d(ecdfi[1], sigma=smooth)
    epdf = [ecdfi[0][:-1], np.diff(ecdfi[1])/np.diff(ecdfi[0])]
    if smooth_pdf and smooth > 0:
        epdf[1] = ndi.filters.gaussian_filter1d(epdf[1], sigma=smooth)
    return epdf

def gaussian_fit_pdf(s, mu0=0, sigma0=1, a0=1, return_all=False,
        leastsq_kwargs={}, **kwargs):
    """Gaussian fit of samples s using a fit to the empirical PDF.
    If a0 is None then only (mu,sigma) are fitted (to a gaussian density).
    `kwargs` are passed to get_epdf().
    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq and the PDF curve is returned.
    """
    ## Empirical PDF
    epdf = get_epdf(s, **kwargs)

    res = gaussian_fit_curve(epdf[0], epdf[1], mu0, sigma0, a0, return_all,
            **leastsq_kwargs)
    if return_all: return res, epdf
    return res

def gaussian_fit_hist(s, mu0=0, sigma0=1, a0=None, bins=np.r_[-0.5:1.5:0.001],
        return_all=False, leastsq_kwargs={}, weights=None, **kwargs):
    """Gaussian fit of samples s fitting the hist to a Gaussian function.
    If a0 is None then only (mu,sigma) are fitted (to a gaussian density).
    kwargs are passed to the histogram function.
    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq and the histogram is returned.
    `weights` optional weights for the histogram.
    """
    histogram_kwargs = dict(bins=bins, density=True, weights=weights)
    histogram_kwargs.update(**kwargs)
    H = np.histogram(s, **histogram_kwargs)
    x, y = 0.5*(H[1][:-1] + H[1][1:]), H[0]
    #bar(H[1][:-1], H[0], H[1][1]-H[1][0], alpha=0.3)

    res = gaussian_fit_curve(x, y, mu0, sigma0, a0, return_all,
            **leastsq_kwargs)
    if return_all: return res, H, x, y
    return res

def gaussian_fit_cdf(s, mu0=0, sigma0=1, return_all=False, **leastsq_kwargs):
    """Gaussian fit of samples s fitting the empirical CDF.
    Additional kwargs are passed to the leastsq() function.
    If return_all=False then return only the fitted (mu,sigma) values
    If return_all=True (or full_output=True is passed to leastsq)
    then the full output of leastsq and the histogram is returned.
    """
    ## Empirical CDF
    ecdf = [np.sort(s), np.arange(0.5, s.size+0.5)*1./s.size]

    ## Analytical Gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(np.sqrt(2)*sigma)))

    ## Fitting the empirical CDF
    err_func = lambda p, x, y: y - gauss_cdf(x, p[0], p[1])
    res = leastsq(err_func, x0=[mu0, sigma0], args=(ecdf[0], ecdf[1]),
            **leastsq_kwargs)
    if return_all: return res, ecdf
    return res[0]

def gaussian_fit_ml(s, mu_sigma_guess=[0.5, 1]):
    """Gaussian fit of samples s using the Maximum Likelihood (ML method).
    Didactical, since scipy.stats.norm.fit() implements the same method.
    """
    n = s.size
    ## Log-likelihood (to be maximized)
    log_l = lambda mu, sig: -n/2.*np.log(sig**2) - \
                             1./(2*sig**2)*np.sum((s-mu)**2)

    ## Function to be minimized
    min_fun = lambda p: -log_l(p[0], p[1])

    res = O.minimize(min_fun, [0, 0.5], method='powell',
                     options={'xtol': 1e-6, 'disp': True, 'maxiter': 1e9})

    print(res)
    mu, sigma = res['x']
    return mu, sigma

##
# Two-component gaussian mixtures
#

def two_gauss_mix_pdf(x, p):
    """PDF for the distribution of a mixture of two Gaussians."""
    mu1, sig1, mu2, sig2, a = p
    return a*normpdf(x, mu1, sig1) + (1-a)*normpdf(x, mu2, sig2)

def two_gauss_mix_ab(x, p):
    """Mixture of two Gaussians with no area constrain."""
    mu1, sig1, a1, mu2, sig2, a2 = p
    return a1*normpdf(x, mu1, sig1) + a2*normpdf(x, mu2, sig2)

def reorder_parameters(p):
    """Reorder 2-gauss mix params to have the 1st component with smaller mean.
    """
    if p[0] > p[2]:
        p = p[np.array([2, 3, 0, 1, 4])]  # swap (mu1, sig1) with (mu2, sig2)
        p[4] = 1 - p[4]                   # "swap" the alpha of the mixture
    return p

def reorder_parameters_ab(p):
    """Reorder 2-gauss mix params to have the 1st component with smaller mean.
    """
    if p[0] > p[3]:
        p = p[np.array([3, 4, 5, 0, 1, 2])]
    return p


def bound_check(val, bounds):
    """Returns `val` clipped inside the interval `bounds`."""
    if bounds[0] is not None and val < bounds[0]:
        val = bounds[0]
    if bounds[1] is not None and val > bounds[1]:
        val = bounds[1]
    return val

def two_gaussian_fit_curve(x, y, p0, return_all=False, verbose=False, **kwargs):
    """Fit a 2-gaussian mixture to the (x,y) curve.
    `kwargs` are passed to the leastsq() function.

    If return_all=False then return only the fitted paramaters
    If return_all=True then the full output of leastsq is returned.
    """
    if kwargs['method'] == 'leastsq':
        kwargs.pop('method')
        kwargs.pop('bounds')
        def err_func(p, x, y):
            return (y - two_gauss_mix_pdf(x, p))
        res = leastsq(err_func, x0=p0, args=(x, y), **kwargs)
        p = res[0]
    else:
        def err_func(p, x, y):
            return ((y - two_gauss_mix_pdf(x, p))**2).sum()
        res = minimize(err_func, x0=p0, args=(x, y), **kwargs)
        p = res.x

    if verbose:
        print(res, '\n')
    if return_all: return res
    return reorder_parameters(p)

def two_gaussian_fit_KDE_curve(s, p0=[0, 0.1, 0.6, 0.1, 0.5], weights=None,
                               bandwidth=0.05, x_pdf=None, debug=False,
                               method='SLSQP', bounds=None,
                               verbose=False, **kde_kwargs):
    """Fit sample `s` with two gaussians using a KDE pdf approximation.

    The 2-gaussian pdf is then curve-fitted to the KDE pdf.

    Arguments:
        s (array): population of samples to be fitted
        p0 (sequence-like): initial parameters [mu0, sig0, mu1, sig1, a]
        bandwidth (float): bandwidth for the KDE algorithm
        method (string): fit method, can be 'leastsq' or one of the methods
            accepted by scipy `minimize()`
        bounds (None or 5-element list): if not None, each element is a
            (min,max) pair of bounds for the corresponding parameter. This
            argument can be used only with L-BFGS-B, TNC or SLSQP methods.
            If bounds are used, parameters cannot be fixed
        x_pdf (array): array on which the KDE PDF is evaluated and curve-fitted
        weights (array): optional weigths, same size as `s` (for ex.
            1/sigma^2 ~ nt).
        debug (bool): if True perfoms more tests and print more info.

    Additional kwargs are passed to scipy.stats.gaussian_kde().

    Returns:
        Array of parameters for the 2-gaussians (5 elements)
    """
    if x_pdf is None: x_pdf = np.linspace(s.min(), s.max(), 1000)

    ## Scikit-learn KDE estimation
    #kde_skl = KernelDensity(bandwidth=bandwidth, **kde_kwargs)
    #kde_skl.fit(x)[:, np.newaxis])
    ## score_samples() returns the log-likelihood of the samples
    #log_pdf = kde_skl.score_samples(x_pdf)[:, np.newaxis])
    #kde_pdf = np.exp(log_pdf)

    ## Weighted KDE estimation
    kde = gaussian_kde_w(s, bw_method=bandwidth, weights=weights)
    kde_pdf = kde.evaluate(x_pdf)

    p = two_gaussian_fit_curve(x_pdf, kde_pdf, p0=p0, method=method,
                               bounds=bounds, verbose=verbose)
    return p


def two_gaussian_fit_EM_b(s, p0=[0, 0.1, 0.6, 0.1, 0.5], weights=None,
                          bounds=[(None, None,)]*5,
                          max_iter=300, ptol=1e-4, debug=False):
    """
    Fit the sample s with two gaussians using Expectation Maximization.

    This version allows setting boundaries for each parameter.

    Arguments:
        s (array): population of samples to be fitted
        p0 (sequence-like): initial parameters [mu0, sig0, mu1, sig1, a]
        bound (tuple of pairs): sequence of (min, max) values that constrain
            the parameters. If min or max are None, no boundary is set.
        ptol (float): convergence condition. Relative max variation of any
            parameter.
        max_iter (int): max number of iteration in case of non convergence.
        weights (array): optional weigths, same size as `s` (for ex.
            1/sigma^2 ~ nt).

    Returns:
        Array of parameters for the 2-gaussians (5 elements)
    """
    assert np.size(p0) == 5
    if weights is None: weights = np.ones(s.size)
    assert weights.size == s.size
    weights *= (1.*weights.size)/weights.sum() # Normalize to (#samples)
    #weights /= weights.sum()  # Normalize to 1
    if debug: assert np.abs(weights.sum() - s.size) < 1e-6
    bounds_mu = [bounds[0], bounds[2]]
    bounds_sig = [bounds[1], bounds[3]]
    bounds_pi0 = bounds[4]

    # Initial guess of parameters and initializations
    mu = np.array([p0[0], p0[2]])
    sig = np.array([p0[1], p0[3]])
    pi_ = np.array([p0[4], 1-p0[4]])

    gamma = np.zeros((2, s.size))
    N_ = np.zeros(2)
    p_new = np.array(p0)

    # EM loop
    counter = 0
    stop_iter, converged = False, False
    while not stop_iter:
        # Compute the responsibility func. (gamma) and the new parameters
        for k in [0, 1]:
            gamma[k, :] = weights*pi_[k]*normpdf(s, mu[k], sig[k]) / \
                            two_gauss_mix_pdf(s, p_new)
            N_[k] = gamma[k, :].sum()
            mu[k] = np.sum(gamma[k]*s)/N_[k]
            mu[k] = bound_check(mu[k], bounds_mu[k])
            sig[k] = np.sqrt( np.sum(gamma[k]*(s-mu[k])**2)/N_[k] )
            sig[k] = bound_check(sig[k], bounds_sig[k])
            if k < 1:
                pi_[k] = N_[k]/s.size
                pi_[k] = bound_check(pi_[k], bounds_pi0)
            else:
                pi_[k] = 1 - pi_[0]
        p_old = p_new
        p_new = np.array([mu[0], sig[0], mu[1], sig[1], pi_[0]])
        if debug:
            assert np.abs(N_.sum() - s.size)/float(s.size) < 1e-6
            assert np.abs(pi_.sum() - 1) < 1e-6

        # Convergence check
        counter += 1
        relative_delta = np.abs(p_new - p_old)/p_new
        converged = relative_delta.max() < ptol
        stop_iter = converged or (counter >= max_iter)

    if debug:
        print("Iterations: ", counter)
    if not converged:
        print("WARNING: Not converged, max iteration (%d) reached." % max_iter)
    return reorder_parameters(p_new)

def two_gaussian_fit_EM(s, p0=[0, 0.1, 0.6, 0.1, 0.5], max_iter=300, ptol=1e-4,
                        fix_mu=[0, 0], fix_sig=[0, 0], debug=False,
                        weights=None):
    """
    Fit the sample s with two gaussians using Expectation Maximization.

    This vesion allows to optionally fix mean or std. dev. of each component.

    Arguments:
        s (array): population of samples to be fitted
        p0 (sequence-like): initial parameters [mu0, sig0, mu1, sig1, a]
        bound (tuple of pairs): sequence of (min, max) values that constrain
            the parameters. If min or max are None, no boundary is set.
        ptol (float): convergence condition. Relative max variation of any
            parameter.
        max_iter (int): max number of iteration in case of non convergence.
        weights (array): optional weigths, same size as `s` (for ex.
            1/sigma^2 ~ nt).

    Returns:
        Array of parameters for the 2-gaussians (5 elements)
    """
    assert np.size(p0) == 5
    if weights is None: weights = np.ones(s.size)
    assert weights.size == s.size
    weights *= (1.*weights.size)/weights.sum() # Normalize to (#samples)
    #weights /= weights.sum()  # Normalize to 1
    if debug: assert np.abs(weights.sum() - s.size) < 1e-6

    # Initial guess of parameters and initializations
    mu = np.array([p0[0], p0[2]])
    sig = np.array([p0[1], p0[3]])
    pi_ = np.array([p0[4], 1-p0[4]])

    gamma = np.zeros((2, s.size))
    N_ = np.zeros(2)
    p_new = np.array(p0)

    # EM loop
    counter = 0
    stop_iter, converged = False, False
    while not stop_iter:
        # Compute the responsibility func. (gamma) and the new parameters
        for k in [0, 1]:
            gamma[k, :] = weights*pi_[k]*normpdf(s, mu[k], sig[k]) / \
                    two_gauss_mix_pdf(s, p_new)
            ## Uncomment for SCHEME2
            #gamma[k, :] = pi_[k]*normpdf(s, mu[k], sig[k]) / \
            #        two_gauss_mix_pdf(s, p_new)
            N_[k] = gamma[k, :].sum()
            if not fix_mu[k]:
                mu[k] = np.sum(gamma[k]*s)/N_[k]
                ## Uncomment for SCHEME2
                #mu[k] = np.sum(weights*gamma[k]*s)/N_[k]
            if not fix_sig[k]:
                sig[k] = np.sqrt( np.sum(gamma[k]*(s-mu[k])**2)/N_[k] )
            pi_[k] = N_[k]/s.size
        p_old = p_new
        p_new = np.array([mu[0], sig[0], mu[1], sig[1], pi_[0]])
        if debug:
            assert np.abs(N_.sum() - s.size)/float(s.size) < 1e-6
            assert np.abs(pi_.sum() - 1) < 1e-6

        # Convergence check
        counter += 1
        fixed = np.concatenate([fix_mu, fix_sig, [0]]).astype(bool)
        relative_delta = np.abs(p_new[-fixed] - p_old[-fixed])/p_new[-fixed]
        converged = relative_delta.max() < ptol
        stop_iter = converged or (counter >= max_iter)

    if debug:
        print("Iterations: ", counter)
    if not converged:
        print("WARNING: Not converged, max iteration (%d) reached." % max_iter)
    return reorder_parameters(p_new)

def two_gaussian_fit_hist(s, bins=np.r_[-0.5:1.5:0.001], weights=None,
        p0=[0.2,1,0.8,1,0.3], fix_mu=[0,0], fix_sig=[0,0], fix_a=False):
    """Fit the sample s with 2-gaussian mixture (histogram fit).

    Uses scipy.optimize.leastsq function. Parameters can be fixed but
    cannot be constrained in an interval.

    Arguments:
        s (array): population of samples to be fitted
        p0 (5-element list or array): initial guess or parameters
        bins (int or array): bins passed to `np.histogram()`
        weights (array): optional weights passed to `np.histogram()`
        fix_a (tuple of bools): Whether to fix the amplitude of the gaussians
        fix_mu (tuple of bools): Whether to fix the mean of the gaussians
        fix_sig (tuple of bools): Whether to fix the sigma of the gaussians

    Returns:
        Array of parameters for the 2-gaussians (5 elements)
    """
    assert np.size(p0) == 5
    fix = np.array([fix_mu[0], fix_sig[0], fix_mu[1], fix_sig[1], fix_a],
                dtype=bool)
    p0 = np.array(p0)
    p0_free = p0[-fix]
    p0_fix = p0[fix]

    H = np.histogram(s, bins=bins, weights=weights, density=True)
    x, y = 0.5*(H[1][:-1] + H[1][1:]), H[0]
    assert x.size == y.size

    ## Fitting
    def err_func(p, x, y, fix, p_fix, p_complete):
        p_complete[-fix] = p
        p_complete[fix] = p_fix
        return y - two_gauss_mix_pdf(x, p_complete)

    p_complete = np.zeros(5)
    p, v = leastsq(err_func, x0=p0_free, args=(x, y, fix, p0_fix, p_complete))
    p_new = np.zeros(5)
    p_new[-fix] = p
    p_new[fix] = p0_fix
    return reorder_parameters(p_new)


def two_gaussian_fit_hist_min(s, bounds=None, method='L-BFGS-B',
        bins=np.r_[-0.5:1.5:0.001], weights=None,  p0=[0.2,1,0.8,1,0.3],
        fix_mu=[0,0], fix_sig=[0,0], fix_a=False, verbose=False):
    """Fit the sample `s` with 2-gaussian mixture (histogram fit). [Bounded]

    Uses scipy.optimize.minimize allowing constrained minimization.

    Arguments:
        s (array): population of samples to be fitted
        method (string): one of the methods accepted by scipy `minimize()`
        bounds (None or 5-element list): if not None, each element is a
            (min,max) pair of bounds for the corresponding parameter. This
            argument can be used only with L-BFGS-B, TNC or SLSQP methods.
            If bounds are used, parameters cannot be fixed
        p0 (5-element list or array): initial guess or parameters
        bins (int or array): bins passed to `np.histogram()`
        weights (array): optional weights passed to `np.histogram()`
        fix_a (tuple of bools): Whether to fix the amplitude of the gaussians
        fix_mu (tuple of bools): Whether to fix the mean of the gaussians
        fix_sig (tuple of bools): Whether to fix the sigma of the gaussians
        verbose (boolean): allows printing fit information

    Returns:
        Array of parameters for the 2-gaussians (5 elements)
    """
    assert np.size(p0) == 5
    fix = np.array([fix_mu[0], fix_sig[0], fix_mu[1], fix_sig[1], fix_a],
                dtype=bool)
    p0 = np.array(p0)
    p0_free = p0[-fix]
    p0_fix = p0[fix]

    H = np.histogram(s, bins=bins, weights=weights, density=True)
    x, y = 0.5*(H[1][:-1] + H[1][1:]), H[0]
    assert x.size == y.size

    ## Fitting
    def err_func(p, x, y, fix, p_fix, p_complete):
        p_complete[-fix] = p
        p_complete[fix] = p_fix
        return ((y - two_gauss_mix_pdf(x, p_complete))**2).sum()

    p_complete = np.zeros(5)
    res = minimize(err_func, x0=p0_free, args=(x, y, fix, p0_fix, p_complete),
                    method=method, bounds=bounds)
    if verbose: print(res)
    p_new = np.zeros(5)
    p_new[-fix] = res.x
    p_new[fix] = p0_fix
    return reorder_parameters(p_new)

def two_gaussian_fit_hist_min_ab(s, bounds=None, method='L-BFGS-B',
        bins=np.r_[-0.5:1.5:0.001], weights=None,  p0=[0.2,1,0.8,1,0.3],
        fix_mu=[0,0], fix_sig=[0,0], fix_a=[0,0], verbose=False):
    """Histogram fit of sample `s` with 2-gaussian functions.

    Uses scipy.optimize.minimize allowing constrained minimization. Also
    each parameter can be fixed.

    The order of the parameters is: mu1, sigma1, a1, mu2, sigma2, a2.

    Arguments:
        s (array): population of samples to be fitted
        method (string): one of the methods accepted by scipy `minimize()`
        bounds (None or 6-element list): if not None, each element is a
            (min,max) pair of bounds for the corresponding parameter. This
            argument can be used only with L-BFGS-B, TNC or SLSQP methods.
            If bounds are used, parameters cannot be fixed
        p0 (6-element list or array): initial guess or parameters
        bins (int or array): bins passed to `np.histogram()`
        weights (array): optional weights passed to `np.histogram()`
        fix_a (tuple of bools): Whether to fix the amplitude of the gaussians
        fix_mu (tuple of bools): Whether to fix the mean of the gaussians
        fix_sig (tuple of bools): Whether to fix the sigma of the gaussians
        verbose (boolean): allows printing fit information

    Returns:
        Array of parameters for the 2-gaussians (6 elements)
    """
    nparams = 6
    assert np.size(p0) == nparams
    fix = np.array([fix_mu[0], fix_sig[0], fix_a[0],
                    fix_mu[1], fix_sig[1], fix_a[1]], dtype=bool)
    p0 = np.asarray(p0)
    p0_free = p0[-fix]
    p0_fix = p0[fix]

    counts, bins = np.histogram(s, bins=bins, weights=weights, density=True)
    x = bins[:-1] + 0.5*(bins[1] - bins[0])
    y = counts
    assert x.size == y.size

    ## Fitting
    def err_func(p, x, y, fix, p_fix, p_complete):
        p_complete[-fix] = p
        p_complete[fix] = p_fix
        return ((y - two_gauss_mix_ab(x, p_complete))**2).sum()

    p_complete = np.zeros(nparams)
    res = minimize(err_func, x0=p0_free, args=(x, y, fix, p0_fix, p_complete),
                    method=method, bounds=bounds)
    if verbose: print(res)
    p_new = np.zeros(nparams)
    p_new[-fix] = res.x
    p_new[fix] = p0_fix
    return reorder_parameters_ab(p_new)

def two_gaussian_fit_cdf(s, p0=[0., .05, .6, .1, .5],
                         fix_mu=[0, 0], fix_sig=[0, 0]):
    """Fit the sample s with two gaussians using a CDF fit.

    Curve fit 2-gauss mixture Cumulative Distribution Function (CDF) to the
    empirical CDF for sample `s`.

    Note that with a CDF fit no weighting is possible.

    Arguments:
        s (array): population of samples to be fitted
        p0 (5-element list or array): initial guess or parameters
        fix_mu (tuple of bools): Whether to fix the mean of the gaussians
        fix_sig (tuple of bools): Whether to fix the sigma of the gaussians

    Returns:
        Array of parameters for the 2-gaussians (5 elements)
    """
    assert np.size(p0) == 5
    fix = np.array([fix_mu[0], fix_sig[0], fix_mu[1], fix_sig[1], 0],
                    dtype=bool)
    p0 = np.array(p0)
    p0_free = p0[-fix]
    p0_fix = p0[fix]

    ## Empirical CDF
    ecdf = [np.sort(s), np.arange(0.5, s.size+0.5)*1./s.size]
    x, y = ecdf

    ## Analytical gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(np.sqrt(2)*sigma)))
    def two_gauss_mix_cdf(x, p):
        return p[4]*gauss_cdf(x, p[0], p[1]) + (1-p[4])*gauss_cdf(x, p[2], p[3])

    ## Fitting the empirical CDF
    def err_func(p, x, y, fix, p_fix, p_complete):
        p_complete[-fix] = p
        p_complete[fix] = p_fix
        return y - two_gauss_mix_cdf(x, p_complete)

    p_complete = np.zeros(5)
    p, v = leastsq(err_func, x0=p0_free, args=(x, y, fix, p0_fix, p_complete))
    p_new = np.zeros(5)
    p_new[-fix] = p
    p_new[fix] = p0_fix
    return reorder_parameters(p_new)

def test_two_gauss():
    m01 = 0.
    m02 = 0.6
    s01 = 0.05
    s02 = 0.1
    alpha = 0.
    p_real = [m01, s01, m02, s02, alpha]

    N = 500
    si1 = round(alpha*N)
    si2 = round((1-alpha)*N)
    s1 = R.normal(size=si1, loc=m01, scale=s01)
    s2 = R.normal(size=si2, loc=m02, scale=s02)
    s = np.r_[s1,s2]

    pc = two_gaussian_fit_cdf(s, fix_mu=[1,0], p0=[-0.01,0.05,0.5,0.2,0.4])
    ph = two_gaussian_fit_hist(s, fix_mu=[1,0], p0=[-0.01,0.05,0.5,0.2,0.4])
    pe = two_gaussian_fit_EM(s, fix_mu=[1,0], p0=[-0.01,0.05,0.5,0.2,0.4])

    hist(s, bins=40, normed=True)

    x = np.r_[s.min()-1:s.max()+1:200j]
    plot(x, a*normpdf(x,mu1,sig1), lw=2)
    plot(x, (1-a)*normpdf(x,mu2,sig2), lw=2)
    plot(x, two_gauss_mix_pdf(x, p0), lw=2)

    axvline(m01, lw=2, color='k', alpha=0.3)
    axvline(m02, lw=2, color='gray', alpha=0.3)
    axvline(mu1, lw=2, ls='--', color='k', alpha=0.3)
    axvline(mu2, lw=2, ls='--', color='gray', alpha=0.3)
    axvline(mu1h, lw=2, ls='--', color='r', alpha=0.3)
    axvline(mu2h, lw=2, ls='--', color='r', alpha=0.3)

def compare_two_gauss():
    m01 = 0.
    m02 = 0.5
    s01 = 0.08
    s02 = 0.15
    alpha = 0.7
    p_real = [m01, s01, m02, s02, alpha]

    N = 1000
    si1 = round(N*alpha)
    si2 = round((1-alpha)*N)

    p0 = [-0.01,0.05,0.6,0.2,0.4]
    fix_mu = [0,0]

    n = 500
    PC, PH, PE = np.zeros((n,5)), np.zeros((n,5)), np.zeros((n,5))
    for i in xrange(n):
        s1 = R.normal(size=si1, loc=m01, scale=s01)
        s2 = R.normal(size=si2, loc=m02, scale=s02)
        s = np.r_[s1,s2]

        pc = two_gaussian_fit_cdf(s, fix_mu=fix_mu, p0=p0)
        ph = two_gaussian_fit_hist(s, fix_mu=fix_mu, p0=p0)
        pe = two_gaussian_fit_EM(s, fix_mu=fix_mu, p0=p0)

        PC[i], PH[i], PE[i] = pc, ph, pe

    Label = ['Mu1', 'Sig1', 'Mu2', 'Sig2', 'Alpha']
    ftype = 'png'

    for i in range(5):
        figure()
        title(Label[i])
        vmin = min([PC[:,i].min(), PH[:,i].min(), PE[:,i].min()])
        vmax = max([PC[:,i].max(), PH[:,i].max(), PE[:,i].max()])
        b = np.r_[vmin:vmax:80j]
        if vmax == vmin: b = np.r_[vmin-.1:vmax+.1:200j]
        hist(PC[:,i], bins=b, alpha=0.3, label='CDF')
        hist(PH[:,i], bins=b, alpha=0.3, label='Hist')
        hist(PE[:,i], bins=b, alpha=0.3, label='EM')
        legend(loc='best')
        axvline(p_real[i], color='k', lw=2)
        #savefig('Two-gaussian Fit Comp - %s.png' % Label[i])

def gaussian2d_fit(sx, sy, guess=[0.5,1]):
    """2D-Gaussian fit of samples S using a fit to the empirical CDF."""
    assert sx.size == sy.size

    ## Empirical CDF
    ecdfx = [np.sort(sx), np.arange(0.5,sx.size+0.5)*1./sx.size]
    ecdfy = [np.sort(sy), np.arange(0.5,sy.size+0.5)*1./sy.size]

    ## Analytical gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(np.sqrt(2)*sigma)))

    ## Fitting the empirical CDF
    fitfunc = lambda p, x: gauss_cdf(x, p[0], p[1])
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    px,v = leastsq(errfunc, x0=guess, args=(ecdfx[0],ecdfx[1]))
    py,v = leastsq(errfunc, x0=guess, args=(ecdfy[0],ecdfy[1]))
    print("2D Gaussian CDF fit", px, py)

    mux, sigmax = px[0], px[1]
    muy, sigmay = py[0], py[1]
    return mux, sigmax, muy, sigmay

def test_gaussian2d_fit():
    mx0 = 0.1
    my0 = 0.9
    sigx0 = 0.4
    sigy0 = 0.25

    Size = 500
    sx = R.normal(size=Size, loc=mx0, scale=sigx0)
    sy = R.normal(size=Size, loc=my0, scale=sigy0)

    mux, sigmax, muy, sigmay = gaussian2d_fit(sx, sy)

    plot(sx, sy, 'o', alpha=0.2, mew=0)

    X,Y = np.mgrid[sx.min()-1:sx.max()+1:200j, sy.min()-1:sy.max()+1:200j]

    def gauss2d(X,Y, mx, my, sigx, sigy):
        return np.exp(-((X-mx)**2)/(2*sigx**2))*np.exp(-((Y-my)**2)/(2*sigy**2))

    contour(X,Y,gauss2d(X,Y,mux,muy,sigmax,sigmay))

    plot(mx0,my0, 'ok', mew=0, ms=10)
    plot(mux,muy, 'x', mew=2, ms=10, color='green')


def two_gaussian2d_fit(sx, sy, guess=[0.5,1]):
    """2D-Gaussian fit of samples S using a fit to the empirical CDF."""
    ## UNFINISHED (I have 2 alphas unp.sign the xy projections)
    assert sx.size == sy.size

    ## Empirical CDF
    ecdfx = [np.sort(sx), np.arange(0.5,sx.size+0.5)*1./sx.size]
    ecdfy = [np.sort(sy), np.arange(0.5,sy.size+0.5)*1./sy.size]

    ## Analytical gaussian CDF
    gauss_cdf = lambda x, mu, sigma: 0.5*(1+erf((x-mu)/(np.sqrt(2)*sigma)))

    gauss2d_cdf = lambda X,Y,mx,sx,my,sy: gauss_cdf(X,mx,sx)*gauss_cdf(Y,my,sy)

    two_cdf = lambda x, m1, s1, m2, s2, a:\
        a*gauss_cdf(x,m1,s1)+(1-a)*gauss_cdf(x,m2,s2)

    two2d_cdf = lambda X,Y, mx1, sx1, mx2, sx2, my1, sy1, my2, sy2, a:\
        a*gauss2d_cdf(X,Y,mx1,sx1,my1,sy1)+(1-a)*gauss_cdf(X,Y,mx2,sx2,my2,sy2)

    ## Fitting the empirical CDF
    fitfunc = lambda p, x: two_cdf(x, *p)
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    fitfunc2d = lambda p, X,Y: two2d_cdf(X,Y, *p)
    errfunc2d = lambda p, X,Y,Z: fitfunc2d(p, X,Y) - Z

    px,v = leastsq(errfunc, x0=guess, args=(ecdfx[0],ecdfx[1]))
    py,v = leastsq(errfunc, x0=guess, args=(ecdfy[0],ecdfy[1]))
    print("2D Two-Gaussians CDF fit", px, py)

    mux1, sigmax1, mux2, sigmax2, alphax = px
    muy1, sigmay1, muy2, sigmay2, alphay = py
    return mu1, sigma1, mu2, sigma2, alpha

def test_gaussian_fit():
    m0 = 0.1
    s0 = 0.4
    size = 500

    s = R.normal(size=size, loc=m0, scale=s0)
    #s = s[s<0.4]
    mu, sig = gaussian_fit(s)
    mu1, sig1 = S.norm.fit(s)
    mu2, sig2 = gaussian_fit_ml(s)

    print("ECDF ", mu, sig)
    print("ML         ", mu1, sig1)
    print("ML (manual)", mu2, sig2)

    H = np.histogram(s, bins=20, density=True)
    h = H[0]
    bw = H[1][1] - H[1][0]
    #bins_c = H[1][:-1]+0.5*bw
    bar(H[1][:-1], H[0], bw, alpha=0.3)

    x = np.r_[s.min()-1:s.max()+1:200j]
    plot(x, normpdf(x,m0,s0), lw=2, color='grey')
    plot(x, normpdf(x,mu,sig), lw=2, color='r', alpha=0.5)
    plot(x, normpdf(x,mu1,sig1), lw=2, color='b', alpha=0.5)

if __name__ == '__main__':
    #compare_two_gauss()
    #test_gaussian2d_fit()
    #test_gaussian_fit()
    #show()
    pass
