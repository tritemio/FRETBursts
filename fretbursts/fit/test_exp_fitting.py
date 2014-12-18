"""
Unittest for exp_fitting.py
"""

from __future__ import print_function
import pytest

import numpy as np
import scipy.stats as SS

from fretbursts.fit.exp_fitting import expon_fit, expon_fit_cdf, expon_fit_hist


sample_size = 5000
sample_tau = 3.
sample_min = 2.
max_relative_error = 0.05


@pytest.fixture(scope="module")
def sample():
    np.random.seed(1)
    return SS.expon.rvs(size=sample_size, scale=sample_tau)


def test_expon_fit(sample):
    lambda_fit, resid, x_resid, size = expon_fit(sample, s_min=sample_min)
    tau_fit = 1./lambda_fit
    relative_error = np.abs(tau_fit-sample_tau)/sample_tau
    print('\n [expon_fit] Fit (tau): %.2f  - Relative error: %.2f %%' % \
            (tau_fit, relative_error*100))
    assert relative_error < max_relative_error

def test_expon_fit_cdf(sample):
    lambda_fit, resid, x_resid, size = expon_fit_cdf(sample, s_min=sample_min)
    tau_fit = 1./lambda_fit
    relative_error = np.abs(tau_fit-sample_tau)/sample_tau
    print('\n [expon_fit_cdf] Fit (tau): %.2f  - Relative error: %.2f %%' % \
            (tau_fit, relative_error*100))
    assert relative_error < max_relative_error

def test_expon_fit_hist(sample):
    binw = sample_tau/20.
    bins = np.arange(0, sample_tau*6, binw)
    lambda_fit, resid, x_resid, size = expon_fit_hist(sample, s_min=sample_min, bins=bins)
    tau_fit = 1./lambda_fit
    relative_error = np.abs(tau_fit-sample_tau)/sample_tau
    print('\n [expon_fit_hist] Fit (tau): %.2f  - Relative error: %.2f %%' % \
            (tau_fit, relative_error*100))
    assert relative_error < max_relative_error

def test_expon_fit_histw(sample):
    binw = sample_tau/20.
    bins = np.arange(0, sample_tau*6, binw)
    lambda_fit, resid, x_resid, size = expon_fit_hist(sample, s_min=sample_min, bins=bins,
                                weights='hist_counts')
    tau_fit = 1./lambda_fit
    relative_error = np.abs(tau_fit-sample_tau)/sample_tau
    print('\n [expon_fit_hist] Fit (tau): %.2f  - Relative error: %.2f %%' % \
            (tau_fit, relative_error*100))
    assert relative_error < max_relative_error

if __name__ == '__main__':
    pytest.main("-x -v -s fretbursts/fit/test_exp_fitting.py")
