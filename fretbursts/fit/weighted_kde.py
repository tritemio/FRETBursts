"""
This module implements a weighted Kernel Density Estimation (KDE).

The code is a modified version of Scipy gaussian_kde that adds support for
weights.

Only direct set of bandwidth is supported.

This code borrows from a mods done by Zachary Pincus (zachary.pincus@yale.edu):
http://mail.scipy.org/pipermail/scipy-user/2013-May/034580.html
"""

from builtins import range, zip

# Scipy imports.
#from scipy import linalg, special
from numpy import (atleast_2d, reshape, zeros, newaxis, dot, exp, pi, sqrt,)
                   #ravel, power, atleast_1d, squeeze, sum, transpose)
import numpy as np



class gaussian_kde_w(object):
    """
    Representation of a kernel-density estimate using Gaussian kernels.

    Modified version including weights.
    """
    def __init__(self, dataset, bw_method, weights=None):
        self.dataset = atleast_2d(dataset)
        if not self.dataset.size > 1:
            raise ValueError("`dataset` input should have multiple elements.")

        self.d, self.n = self.dataset.shape
        self.set_bandwidth(bw_method=bw_method)
        if weights is None:
            weights = np.ones(dataset.shape, dtype=float)
        else:
            weights = np.asarray(weights, dtype=float)
        self.weights = weights / weights.sum()
        #self.inv_cov = inv_cov
        #self._norm_factor = norm_factor

    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError : if the dimensionality of the input points is different than
                     the dimensionality of the KDE.

        """
        points = atleast_2d(points)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % \
                        (d, self.d)
                raise ValueError(msg)

        result = zeros((m,), dtype=np.float)

        if m >= self.n:
            # there are more points than data, so loop over data
            for i in range(self.n):
                diff = self.dataset[:, i, newaxis] - points
                tdiff = dot(self.inv_cov, diff)
                energy = np.sum(diff * tdiff,axis=0) / 2.0
                #result = result + exp(-energy)                 # weights mod
                result = result + self.weights[i]*exp(-energy)  # weights mod
        else:
            # loop over points
            for i in range(m):
                diff = self.dataset - points[:, i, newaxis]
                tdiff = dot(self.inv_cov, diff)
                energy = np.sum(diff * tdiff, axis=0) / 2.0
                #result[i] = np.sum(exp(-energy), axis=0)              # weights mod
                result[i] = np.sum(self.weights*exp(-energy), axis=0)  # weights mod

        result = result / self._norm_factor
        return result

    __call__ = evaluate

    def set_bandwidth(self, bw_method=None):
        """Compute the estimator bandwidth with given method.

        The new bandwidth calculated after a call to `set_bandwidth` is used
        for subsequent evaluations of the estimated density.

        Parameters
        ----------
        bw_method : scalar
            The method sets the estimator bandwidth.

        """
        if np.isscalar(bw_method):
            self._bw_method = 'use constant'
            self.covariance_factor = lambda: bw_method
        else:
            msg = "`bw_method` should be a scalar "
            raise ValueError(msg)

        self._compute_covariance()

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor().
        """
        self.factor = self.covariance_factor()
        self.covariance = self.factor**2
        self.inv_cov = 1. / self.covariance
        self._norm_factor = sqrt(2*pi*self.covariance) #* self.n

