"""
Unittest for exp_fitting.py
"""

import pytest

import numpy as np
import scipy.stats as SS

from fit.exp_fitting import expon_fit, expon_fit_cdf, expon_fit_hist


sample_size = 1000
sample_tau = 3.
sample_min = 2.
max_relative_error = 0.05


@pytest.fixture(scope="module")
def sample():
    np.random.seed(1)
    return SS.expon.rvs(size=sample_size, scale=sample_tau)
    

def test_expon_fit(sample):
    lambda_fit = expon_fit(sample, s_min=sample_min)
    tau_fit = 1./lambda_fit
    relative_error = np.abs(tau_fit-sample_tau)/sample_tau
    print '\n [expon_fit] Fit (tau): %.2f  - Relative error: %.2f %%' % \
            (tau_fit, relative_error)
    assert relative_error < max_relative_error

def test_expon_fit_cdf(sample):
    lambda_fit = expon_fit_cdf(sample, s_min=sample_min)
    tau_fit = 1./lambda_fit
    relative_error = np.abs(tau_fit-sample_tau)/sample_tau
    print '\n [expon_fit_cdf] Fit (tau): %.2f  - Relative error: %.2f %%' % \
            (tau_fit, relative_error)
    assert relative_error < max_relative_error

def test_expon_fit_hist(sample):
    lambda_fit = expon_fit_hist(sample, s_min=sample_min, bins=0.1e-3)
    tau_fit = 1./lambda_fit
    relative_error = np.abs(tau_fit-sample_tau)/sample_tau
    print '\n [expon_fit_hist] Fit (tau): %.2f  - Relative error: %.2f %%' % \
            (tau_fit, relative_error)
    assert relative_error < max_relative_error

if __name__ == '__main__':
    pytest.main("-x -v -s fit/test_exp_fitting.py")
