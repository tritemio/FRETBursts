#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
This model provides a class for fitting multi-channel data
(:class:`MultiFitter`) and a series of predefined functions for common
models used to fit E or S histograms.
"""

from __future__ import division, print_function, absolute_import
from builtins import range, zip

import numpy as np
import pandas as pd
import lmfit
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
    has_matplotlib = False
else:
    has_matplotlib = True

from .fit import gaussian_fitting as gf


## Utility functions
def find_max(x, y, xmin=None, xmax=None):
    """Find peak position of a curve (x, y) between `xmin` and `xmax`.
    """
    if xmin is None:
        xmin = x.min()
    if xmax is None:
        xmax = x.max()

    mask = np.where((x >= xmin)*(x <= xmax))
    return x[mask][y[mask].argmax()]


## Base function shapes
def gaussian(x, center, sigma, amplitude=1):
    """A gaussian function. Parameters: `center`, `sigma`, `amplitude`.
    """
    return amplitude * np.exp(-(x - center)**2/(2*sigma**2))

def asym_gaussian(x, center, sigma1, sigma2, amplitude):
    """A asymmetric gaussian function composed by two gaussian halves.

    This function is composed from two gaussians joined at their peak, so that
    the left and right side decay with different sigmas.

    Arguments:
        x (array): 1-D array for the independent variable
        center (float): function peak position
        sigma1 (float): sigma of the left-side gaussian (for x < center)
        sigma2 (float): sigma of the right-side gaussian (for x > center)
        amplitude (float): maximum value reach for x = center.

    Return:
        An array (same shape as `x`) with the function values.
    """

    left_half = x < center
    sigma = np.zeros_like(x)
    sigma[left_half] = sigma1
    sigma[-left_half] = sigma2
    u = (x - center)/sigma
    return amplitude*np.exp(-u*u)

def bridge_function(x, center1, center2, sigma1, sigma2, amplitude):
    """A "bridge" function, complementary of two gaussian peaks.

    Let `g` be a Gaussian function (with amplitude = 1), the bridge function
    is defined as::

        amplitude * (1 - g(x, center1, sigma1) - g(x, center2, sigma2))

    for `center1 < x < center2`. The function is 0 otherwise.

    Arguments:
        x (array): 1-D array for the independent variable
        center1 (float): center of the first gaussian (left side)
        center2 (float): center of the second gaussian (right side)
        sigma1 (float): sigma of the left-side gaussian
        sigma2 (float): sigma of the right-side gaussian
        amplitude (float): maximum (asymptotic) value of the bridge (plateau)

    Return:
        An array (same shape as `x`) with the function values.
    """
    assert center1 <= center2
    assert amplitude >= 0
    mask = (x > center1)*(x < center2)
    y = np.zeros_like(x)
    y[mask] = 1.
    y[mask] -= gaussian(x[mask], center1, sigma1)
    y[mask] -= gaussian(x[mask], center2, sigma2)
    y *= amplitude
    return y

def bridge_function2(x, center1, center2, sigma1, sigma2, amplitude):
    """A "bridge" function, geometric complementary of two gaussian peaks.

    Let `g` be a Gaussian function (with amplitude = 1), the bridge function
    is defined as:

        amplitude * (1 - g(x, center1, sigma1)) * (1 - g(x, center2, sigma2))

    for center1 < x < center2. The function is 0 otherwise.

    Arguments:
        x (array): 1-D array for the independent variable
        center1 (float): center of the first gaussian (left side)
        center2 (float): center of the second gaussian (right side)
        sigma1 (float): sigma of the left-side gaussian
        sigma2 (float): sigma of the right-side gaussian
        amplitude (float): maximum (asymptotic) value of the bridge (plateau)

    Return:
        An array (same shape as `x`) with the function values.
    """
    assert center1 <= center2
    assert amplitude >= 0
    mask = (x > center1)*(x < center2)
    y = np.zeros_like(x)
    y[mask] = amplitude
    y[mask] *= 1 - gaussian(x[mask], center1, sigma1)
    y[mask] *= 1 - gaussian(x[mask], center2, sigma2)
    y *= amplitude
    return y


##
# Factory functions that return initialized `lmfit.Model` objects
#
def factory_gaussian(center=0.1, sigma=0.1, amplitude=1):
    """Return an lmfit Gaussian model that can be used to fit data.

    Arguments are initial values for the model parameters.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    model = lmfit.models.GaussianModel()
    model.set_param_hint('center', value=center, min=-1, max=2)
    model.set_param_hint('sigma', value=sigma)
    model.set_param_hint('amplitude', value=amplitude)
    return model

def factory_asym_gaussian(center=0.1, sigma1=0.1, sigma2=0.1, amplitude=1):
    """Return a lmfit Asymmetric Gaussian model that can be used to fit data.

    For the definition of asymmetric Gaussian see :func:`asym_gaussian`.
    Arguments are initial values for the model parameters.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    model = lmfit.model.Model(asym_gaussian)
    model.set_param_hint('center', value=center, min=-1, max=2)
    model.set_param_hint('sigma1', value=sigma1)
    model.set_param_hint('sigma2', value=sigma2)
    model.set_param_hint('amplitude', value=amplitude)
    return model

def factory_two_gaussians(add_bridge=False, p1_center=0.1, p2_center=0.9,
                          p1_sigma=0.03, p2_sigma=0.03):
    """Return a 2-Gaussian + (optional) bridge model that can fit data.

    The optional "bridge" component (i.e. a plateau between the two peaks)
    is a function that is non-zero only between `p1_center` and `p2_center`
    and is defined as::

        br_amplitude * (1 - g(x, p1_center, p1_sigma) - g(x, p1_center, p2_sigma))

    where `g` is a gaussian function with amplitude = 1 and `br_amplitude`
    is the height of the plateau and the only additional parameter introduced
    by the bridge. Note that both centers and sigmas parameters in the bridge
    are the same ones of the adjacent Gaussian peaks. Therefore a
    2-Gaussian + bridge model has 7 free parameters: 3 for each Gaussian and
    an additional one for the bridge.
    The bridge function is implemented in :func:`bridge_function`.

    Arguments:
        p1_center, p2_center (float): initial values for the centers of the
            two Gaussian components.
        p1_sigma, p2_sigma (float): initial values for the sigmas of the
            two Gaussian components.
        add_bridge (bool): if True adds a bridge function between the two
            gaussian peaks. If False the model has only two Gaussians.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    peak1 = lmfit.models.GaussianModel(prefix='p1_')
    peak2 = lmfit.models.GaussianModel(prefix='p2_')
    model = peak1 + peak2

    model.set_param_hint('p1_center', value=p1_center, min=-1, max=2)
    model.set_param_hint('p2_center', value=p2_center, min=-1, max=2)
    model.set_param_hint('p1_sigma', value=p1_sigma, min=0.01, max=0.2)
    model.set_param_hint('p2_sigma', value=p2_sigma, min=0.01, max=0.2)
    model.set_param_hint('p1_amplitude', value=1, min=0.01)
    model.set_param_hint('p2_amplitude', value=1, min=0.01)
    name = '2-gaussians'

    if add_bridge:
        bridge = lmfit.model.Model(bridge_function, prefix='br_')
        model += bridge
        model.set_param_hint('br_amplitude', value=0.0001, min=0)
        model.set_param_hint('br_center1', min=0, expr='p1_center')
        model.set_param_hint('br_center2', min=0, expr='p2_center')
        model.set_param_hint('br_sigma1', min=0, expr='p1_sigma')
        model.set_param_hint('br_sigma2', min=0, expr='p2_sigma')
        name += '-bridge'
    model.name = name
    return model

def factory_three_gaussians(p1_center=0., p2_center=0.5, p3_center=1,
                            sigma=0.05):
    """Return a 3-Gaussian model that can fit data.

    The other arguments are initial values for the `center` for each
    Gaussian component plus an single `sigma` argument that is used
    as initial sigma for all the Gaussians. Note that during the fitting
    the sigma of each Gaussian is varied independently.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    peak1 = lmfit.models.GaussianModel(prefix='p1_')
    peak2 = lmfit.models.GaussianModel(prefix='p2_')
    peak3 = lmfit.models.GaussianModel(prefix='p3_')
    model = peak1 + peak2 + peak3

    model.set_param_hint('p1_center', value=p1_center, min=0, max=1)
    model.set_param_hint('p2_center', value=p2_center, min=0, max=1)
    model.set_param_hint('p3_center', value=p3_center, min=0, max=1)
    model.set_param_hint('p1_sigma', value=sigma, min=0.02, max=0.2)
    model.set_param_hint('p2_sigma', value=sigma, min=0.02, max=0.2)
    model.set_param_hint('p3_sigma', value=sigma, min=0.02, max=0.2)
    model.set_param_hint('p1_amplitude', value=1, min=0.01)
    model.set_param_hint('p2_amplitude', value=1, min=0.01)
    model.set_param_hint('p3_amplitude', value=1, min=0.01)
    model.name = '3-gaussians'
    return model

## lmfit composite model used for fitting
def factory_two_asym_gaussians(add_bridge=False, p1_center=0.1, p2_center=0.9,
                               p1_sigma=0.03, p2_sigma=0.03):
    """Return a 2-Asym-Gaussians + (optional) bridge model that can fit data.

    The Asym-Gaussian function is :func:`asym_gaussian`.

    Arguments:
        add_bridge (bool): if True adds a bridge function between the two
            gaussian peaks. If False the model has only two Asym-Gaussians.

    The other arguments are initial values for the model parameters.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    peak1 = lmfit.model.Model(asym_gaussian, prefix='p1_')
    peak2 = lmfit.model.Model(asym_gaussian, prefix='p2_')
    model = peak1 + peak2

    model.set_param_hint('p1_center', value=p1_center, min=-1, max=2)
    #model.set_param_hint('p2_center', value=p2_center, min=-1, max=2)
    model.set_param_hint('p1_sigma1', value=p1_sigma, min=0.01, max=0.2)
    model.set_param_hint('p1_sigma2', value=p1_sigma, min=0.01, max=0.2)
    model.set_param_hint('p2_sigma1', value=p2_sigma, min=0.01, max=0.2)
    model.set_param_hint('p2_sigma2', value=p2_sigma, min=0.01, max=0.2)
    model.set_param_hint('p1_amplitude', value=1, min=0.01)
    model.set_param_hint('p2_amplitude', value=1, min=0.01)

    model.set_param_hint('pos_delta', value=0.3, min=0)
    model.set_param_hint('p2_center', min=-1, expr='p1_center + pos_delta')

    name = '2-asym-gaussians'

    if add_bridge:
        bridge = lmfit.model.Model(bridge_function, prefix='br_')
        model += bridge
        model.set_param_hint('br_amplitude', value=0.0001, min=0)
        model.set_params_hint('br_center1', min=-1, expr='p1_center')
        model.set_params_hint('br_center2', min=-1, expr='p2_center')
        model.set_params_hint('br_sigma1', min=0.01, expr='p1_sigma2')
        model.set_params_hint('br_sigma2', min=0.01, expr='p2_sigma1')
        name += '-bridge'
    model.name = name
    return model


##
# Classes to perform fit-related operations on multi-channel data
#
class FitterBase(object):
    """Base class for histogramming a dataset.

    To set weights assign attribute an array with same size as `data` to
    the attribute `.weights`.
    """
    def __init__(self, data):
        self.data = data
        self.weights = None
        self._hist_computed = False

    def histogram(self, **kwargs):
        """Compute the histogram of the data.

        All the kwargs are passed to `numpy.histogram`.
        """
        kwargs.update(density=False)
        if self.weights is not None:
            kwargs.update(weights=self.weights)
        hist_counts, bins = np.histogram(self.data, **kwargs)
        self._set_hist_data(hist_counts, bins)

    def _set_hist_data(self, hist_counts, bins):
        self.hist_bins = bins
        self.hist_binwidth = (bins[1] - bins[0])
        self.hist_axis = bins[:-1] + 0.5*self.hist_binwidth
        self.hist_counts = np.array(hist_counts)
        self.hist_pdf = np.array(hist_counts, dtype=np.float)
        self.hist_pdf /= self.hist_counts.sum(1)[:, np.newaxis]
        self.hist_pdf /= self.hist_binwidth
        self._hist_computed = True

    @property
    def x_axis(self):
        if not hasattr(self, '_x_axis'):
            self._x_axis = np.linspace(self.hist_axis[0],
                                       self.hist_axis[-1], 1000)
        return self._x_axis


class MultiFitter(FitterBase):
    """A class handling a list of 1-D datasets for histogramming, KDE, fitting.

    This class takes a list of 1-D arrays of samples (such as E values
    per burst). The list contains one 1-D array for each channel in
    a multispot experiment. In single-spot experiments the list contains only
    one array of samples.
    For each dataset in the list, this class compute histograms, KDEs and
    fits (both histogram fit and KDE maximum). The list of datasets is
    stored in the attribute `data_list`.
    The histograms can be fitted with an arbitrary model (lmfit.Model).
    From KDEs the peak position in a range can be estimated.

    Optionally weights can be assigned to each element in a dataset.
    To assign weights a user can assign the `.weights` attribute with a list
    of arrays; corresponding arrays in `.weights` and `.data_list` must have
    the same size.

    Alternatively a function returning the weights can be used. In this case,
    the method `.set_weights_func` allows to set the function to be called
    to generate weights.
    """
    def __init__(self, data_list):
        self.data_list = data_list
        self.ndata = len(data_list)
        self._weights = [None]*self.ndata
        self._hist_computed = False
        self._kde_computed = False

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, values):
        try:
            self._weights = [v for v in values]
        except TypeError:
            # values is not iterable
            self._weights = [values]*self.ndata

    def histogram(self, binwidth=0.03, bins=None, verbose=False, **kwargs):
        """Compute the histogram of the data for each channel.

        If `bins` is None, `binwidth` determines the bins array (saved in
        `self.hist_bins`). If `bins` is not None, `binwidth` is ignored and
        `self.hist_binwidth` is computed from `self.hist_bins`.

        The kwargs and `bins` are passed to `numpy.histogram`.
        """
        # Check if the same histogram is already computed
        if self._hist_computed:
            if bins is None and self.hist_binwidth == binwidth:
                return
            elif bins is not None and np.array_equal(self.hist_bins, bins):
                return

        if verbose:
            print(" - Computing histogram.")
        if bins is None:
            bins = np.r_[-0.2 : 1.2 : binwidth]
        kwargs.update(bins=bins, density=False)
        hist_counts = []
        for data, weights in zip(self.data_list, self.weights):
            if weights is not None:
                kwargs.update(weights=weights)
            counts, bins = np.histogram(data, **kwargs)
            hist_counts.append(counts)
        self._set_hist_data(hist_counts, bins)

    def set_weights_func(self, weight_func, weight_kwargs=None):
        """Setup of the function returning the weights for each data-set.

        To compute the weights for each dataset the `weight_func` is called
        multiple times. Keys in `weight_kwargs` are arguments of
        `weight_func`. Values in `weight_kwargs` are either scalars, in which
        case they are passed to `weight_func`, or lists. When an argument
        is a list, only one element of the list is passed each time.

        Arguments:
            weight_func (function): function that returns the weights
            weight_kwargs (dict): keyword arguments to be passed to
                `weight_func`.
        """
        self.weights = []
        for i in range(self.ndata):
            weight_kw_i = {k: v[i] if isinstance(v, list) else v
                           for k, v in weight_kwargs.items()}
            self.weights.append(weight_func(**weight_kw_i))

    def fit_histogram(self, model=None, pdf=True, **fit_kwargs):
        """Fit the histogram of each channel using the same lmfit model.

        A list of `lmfit.Minimizer` is stored in `.fit_res`.
        The fitted values for all the parameters and all the channels are
        save in a Pandas DataFrame `.params`.

        Arguments:
            model (lmfit.Model object): lmfit Model with all the parameters
                already initialized used for fitting.
            pdf (bool): if True fit the normalized histogram (.hist_pdf)
                otherwise fit the raw counts (.hist_counts).
            fit_kwargs (dict): keyword arguments passed to `model().fit`.
        """
        if model is not None:
            self.model = model
        if not self._hist_computed:
            self.histogram()

        data_list = self.hist_pdf if pdf else self.hist_counts
        self.params = pd.DataFrame(index=range(self.ndata),
                                   columns=sorted(self.model.param_names))
        self.fit_res = []
        #init_params = copy.deepcopy(self.model.params)
        for i, data in enumerate(data_list):
            self.fit_res.append(
                self.model.fit(data, x=self.hist_axis, **fit_kwargs)
                )
            self.params.iloc[i] = pd.Series(self.fit_res[-1].values)

    def calc_kde(self, bandwidth=0.03):
        """Compute the list of kde functions and save it in `.kde`.
        """
        if self._kde_computed and self.kde_bandwidth == bandwidth:
            return

        self.kde_bandwidth = bandwidth
        self.kde = []
        for data, weights_i in zip(self.data_list, self.weights):
            self.kde.append(
                gf.gaussian_kde_w(data, bw_method=bandwidth,
                                  weights=weights_i))
        self._kde_computed = True

    def find_kde_max(self, x_kde, xmin=None, xmax=None):
        """Finds the peak position of kde functions between `xmin` and `xmax`.

        Results are saved in the list `.kde_max_pos`.
        """
        if not hasattr(self, 'kde'):
            self.calc_kde()
        self.kde_max_pos = np.zeros(self.ndata)
        for ich, kde in enumerate(self.kde):
            self.kde_max_pos[ich] = \
                    find_max(x_kde, kde(x_kde), xmin=xmin, xmax=xmax)


def plot_mfit(fitter, ich=0, residuals=False, ax=None, plot_kde=False,
              plot_model=True, return_fig=False, bins=None):
    """Plot data histogram and fitted model from a :class:`MultiFitter` object.

    Assumes data between 0 and 1 and a two peaks model with parameters
    `p1_center` and `p2_center`.

    Return
        A matplotlib figure object.
    """
    if not has_matplotlib:
        raise RuntimeError('Matplotlib not installed.')
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4.5))
    else:
        fig = ax.figure

    if not fitter._hist_computed or bins is not None:
        fitter.histogram(bins=bins)
    ax.plot(fitter.hist_axis, fitter.hist_pdf[ich], '-o', alpha=0.5)
    num_bursts = fitter.data_list[ich].size
    ax.set_title('CH = %d #Bursts %d' % (ich, num_bursts))
    ax.set_xlim(-0.19, 1.19)
    ax.grid(True)
    red = '#E41A1C'

    if hasattr(fitter, 'fit_res') and plot_model:
        fit_res = fitter.fit_res[ich]
        x = fitter.x_axis
        ax.plot(x, fit_res.model.eval(x=x, **fit_res.values), 'k', alpha=0.8)
        if  fit_res.model.components is not None:
            for component in fit_res.model.components:
                ax.plot(x, component.eval(x=x, **fit_res.values), '--k',
                        alpha=0.8)
        for param in fitter.params:
            if param.endswith('center'):
                ax.axvline(fitter.params[param][ich], ls='--', color=red)

        if residuals:
            ax.xaxis.set_ticklabels([])
            divider = make_axes_locatable(ax)
            ax_resid = divider.append_axes("bottom", size=1.2, pad=0.1,
                                           sharex=ax)
            ax_resid.plot(fitter.hist_axis,
                          fitter.hist_pdf[ich] - fit_res.best_fit)
            ax_resid.set_yticks(np.r_[-2:2:0.5])
            ax_resid.set_ylim(-0.49, 0.49)
            ax_resid.set_xlim(-0.2, 1.2)
    elif hasattr(fitter, 'kde') and plot_kde:
        kde = fitter.kde[ich]
        x = fitter.x_axis
        ax.plot(x, kde(x), color='gray', alpha=0.8)
        if hasattr(fitter, 'kde_max_pos'):
            ax.axvline(fitter.kde_max_pos[ich], ls='--', color=red)

    if return_fig:
        return fig
