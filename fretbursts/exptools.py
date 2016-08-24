# -*- coding: utf-8 -*-
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
Tools for handling and testing exponential-like samples.

References:

[1] http://itl.nist.gov/div898/handbook/eda/section3/eda35e.htm

[2] Stephens, M. A. (1974). EDF Statistics for Goodness of Fit and Some
    Comparisons, Journal of the American Statistical Association,
    Vol. 69, pp. 730-737.

[3] Stephens, M. A. (1976). Asymptotic Results for Goodness-of-Fit Statistics
    with Unknown Parameters, Annals of Statistics, 4, pp. 357-369.

[4] Parr, W. C. and Schucany, W. R. (1980). Minimum Distance and Robust
    Estimation, Journal of the American Statistical Association, Vol. 75,
    No. 371 (Sep., 1980), pp. 616-624

[5] Stephens, M. A. (1978). Goodness of fit test with special reference
    to test for exponentiality. Technical Report No 262.
    Department of Statistics, Stanford University.

"""


from __future__ import division

import numpy as np
import math

# Critical points for different statistics
critical_points = dict(
    KS=[0.926, 0.990, 1.094, 1.190, 1.308],
    CM=[0.149, 0.177, 0.224, 0.273, 0.337],
    Wa=[0.112, 0.130, 0.161, 0.191, 0.230],
    AD=[0.922, 1.078, 1.341, 1.606, 1.957],
)

# Significance level in percentage for `critical_points`
significance = [15, 10, 5, 2.5, 1]


def weighted_median(data, weights=None):
    data = np.asfarray(data)
    if weights is None:
        return np.median(data)
    weights = np.asfarray(weights)
    data_argsort = data.argsort()
    sorted_data = data[data_argsort]
    sorted_weights = weights[data_argsort]

    # See http://en.wikipedia.org/wiki/Percentile#Weighted_percentile
    Sn = np.cumsum(sorted_weights)
    prob_n = (Sn - 0.5*sorted_weights)/Sn[-1]
    return np.interp(0.5, prob_n, sorted_data)


def estimate_tau(sample, median=False, weights=None):
    """Estimate the `tau` parameter from an exponentially-distributed `sample`.

    Arguments:
        sample (array): the exponetially-distributed samples
        median (bool): if False thes mean estimator is mean(sample). If True
            uses median(samples)/ln(2) for mean esitmator (more robust).
        weights (array or None): optional array of sample weights.

    Returns:
        An estimation of the `tau` parameter (the mean, inverse of rate).
    """
    if median == False:
        fitted_tau = np.average(sample, weights=weights)
    else:
        fitted_tau = weighted_median(sample, weights)/np.log(2)
        if weights is None:
            assert np.allclose(fitted_tau*np.log(2), np.median(sample))
    return fitted_tau

def tail_mean(sample, threshold=0, weights=None, median=False,
              return_ci=False):
    """Estimate `tau` and num. samples of the exponetial tail of `samples`.
    """
    assert hasattr(sample, '__array__'), \
        'The first argument must be array-like, not type %s.' % type(sample)
    mask = sample >= threshold

    weights_th = None
    if weights is not None:
        weights_th = weights[mask]

    sample_tail = sample[mask] - threshold
    mean_ = np.nan
    if mask.sum() > 0:
        mean_ = estimate_tau(sample_tail, median=median, weights=weights_th)

    num_samples = mask.sum()
    if weights_th is not None and weights_th.size > 0:
        num_samples = weights_th.sum()

    if return_ci:
        # See http://en.wikipedia.org/wiki/Exponential_distribution#Confidence_intervals
        delta = 1.96/np.sqrt(num_samples)
        ci = (mean_/(1 + delta), mean_/(1 - delta))
        return mean_, num_samples, ci
    else:
        return mean_, num_samples

def select_tail(sample_tot, threshold):
    return sample_tot[sample_tot >= threshold] - threshold

def zeta_values(sorted_sample, median=False):
    assert (sorted_sample >= 0).all()
    fitted_tau = estimate_tau(sorted_sample, median=median)
    zeta = 1 - np.exp(-sorted_sample/fitted_tau)
    assert (zeta >= 0).all()
    return zeta


def kolgomorv_stat(zeta):
    n = zeta.size
    i = np.arange(1, n + 1)
    Dplus = np.max(i/n - zeta)
    Dminus = np.max(zeta - (i - 1)/n)
    D = max(Dplus, Dminus)
    return D

def kolgomorv_stat_n(zeta):
    D = kolgomorv_stat(zeta)
    n = zeta.size
    sqrt_n = math.sqrt(n)
    Dstar = (D - 0.2/n)*(sqrt_n + 0.26 + 0.5/sqrt_n)
    return Dstar


def cramervonmises_stat(zeta):
    n = zeta.size
    i = np.arange(1, n + 1)
    term2 = (2*i - 1)/(2*n)
    W2 = np.sum((zeta -  term2)**2) + 1/(12*n)
    return W2

def cramervonmises_stat_n(zeta):
    W2 = cramervonmises_stat(zeta)
    n = zeta.size
    W2star = W2*(1 + 0.16/n)
    return W2star


def watson_stat(zeta):
    n = zeta.size
    W2 = cramervonmises_stat(zeta)
    U2 = W2 - n*(np.mean(zeta) - 0.5)**2
    return U2

def watson_stat_n(zeta):
    U2 = watson_stat(zeta)
    n = zeta.size
    U2star = U2*(1 + 0.16/n)
    return U2star


def andersondarling_stat(zeta):
    n = zeta.size
    i = np.arange(1, n + 1)
    f1 = (2*i - 1)
    log1 = np.log(zeta)
    log2 = np.log(1 - zeta[::-1])
    A2 = -np.sum(f1*(log1 + log2))/n - n
    return A2

def andersondarling_stat_n(zeta):
    A2 = andersondarling_stat(zeta)
    n = zeta.size
    A2star = A2*(1 + 0.6/n)
    return A2star


def exp_test_stat(sample, threshold, median=False, metric='KS',
                  asymptotic=True):
    """Return the specified statistic for exponetiality test.

    Arguments:
        sample (array): supposedly exponential distributed samples.
        threshold (float): theshold used to leselt the sample "tail",
            i.e. `sample[sample > threshold]`.
        median (bool): if False, estimate the sample mean using the emirical
            mean. If True, the sample mean is estimated using the median.
        metric (string): the particular metric to compute. Valid values are:
            'KS' Kolmogorov-Smirnov D statistic,
            'CM' Crames von Mises W^2 statistic,
            'Wa' Watson U^2 statistic,
            'AD' Anderson-Darling A^2 statistic.
        asymptotic (bool): if True use the asymptotic form for the statistic,
            if False use the finite-sample correction.

    Return:
        The requested statistic (a float) for `sample[sample > threshold]`.
    """
    metric_names = {'KS': kolgomorv_stat, 'CM': cramervonmises_stat,
                    'AD': andersondarling_stat, 'Wa': watson_stat}
    metric_names_n = {'KS': kolgomorv_stat_n, 'CM': cramervonmises_stat_n,
                      'AD': andersondarling_stat_n, 'Wa': watson_stat_n}
    assert metric in metric_names
    sample_tail = select_tail(sample, threshold)
    zeta = zeta_values(np.sort(sample_tail), median=median)
    if asymptotic:
        stat_func = metric_names[metric]
    else:
        stat_func = metric_names_n[metric]
    return stat_func(zeta)


def exp_tail_stats(sample, thresholds, metric, asymptotic, median,
                   weights=None):
    """Compute statistics on the sample tail for an array of thresholds.

    Return arrays of mean_fit, conf. intervals, num of samples in the tail
    and the specified exponetiality test statistic.
    """
    mean_fit = np.zeros(thresholds.size)*np.nan
    mean_ci = np.zeros((thresholds.size, 2))*np.nan
    num_samples = np.zeros(thresholds.size)*np.nan
    stats = np.zeros(thresholds.size)*np.nan
    for idx, th in enumerate(thresholds):
        mean_fit[idx], num_samples[idx], mean_ci[idx] = \
            tail_mean(sample, threshold=th,
                      median=median, return_ci=True)
        if num_samples[idx] == 0: break
        stats[idx] = exp_test_stat(sample,
                                   threshold=th, metric=metric,
                                   median=median, asymptotic=asymptotic)
    return mean_fit, mean_ci, num_samples, stats

def exp_dist_amplitude(meantau_th, meantau_th_ci, num_samples_th,
                       thresholds, mean_fitrange):
    """Compute total exponential distribution parameters from the tail.

    Most input arguments are quantities computed on a sample sliced
    with and array of thresholds.

    Arguments:
        meantau_th (array): estimated means as a function of the `thresholds`.
        meantau_th_ci (array): estimated 95% CI for `meantau_th`.
        num_samples (array): number of samples in the tail for each threshold.
        thresholds (array): the different thresholds used to slice the sample.
        meanfit_range (array): min-max values used to fit a single mean
            from `meantau_th`.

    Returns:
        meantau (float): the mean fitted tau for the total distribution.
        meantau_ci (array): 95% confidence interval for meantau.
        alpha (array): `thresholds / meantau`.
        alpha_ci (array): 95% confidence intervals for `alpha`.
        N_fit (int): number of sample in the total distribution.
        N_ci (array): 95% confidence intervals for `N_fit`.
    """
    mean_fitmask = (thresholds >= mean_fitrange[0]) * \
                   (thresholds <= mean_fitrange[1])
    meantau = meantau_th[mean_fitmask].mean()
    meantau_ci = meantau_th_ci[mean_fitmask].mean(0)
    alpha = thresholds/meantau
    alpha_ci = thresholds[:, np.newaxis]/meantau_ci[np.newaxis, :]
    N_fit = num_samples_th*np.exp(alpha)
    N_ci = num_samples_th[:, np.newaxis]*np.exp(alpha_ci)
    return meantau, meantau_ci, alpha, alpha_ci, N_fit, N_ci
