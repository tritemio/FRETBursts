"""
Fit multiple data populations (multi-channel data). Each population/channel
is treated independently.
"""

from __future__ import division
from functools import partial
import numpy as np
import pandas as pd
import lmfit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import fretbursts.fit.gaussian_fitting as gf
import fretbursts.burstlib_ext as bext


## Base function shapes
def gaussian(x, center, sigma, amplitude=1):
    """A gaussian function. Parameters: `center`, `sigma`, `amplitude`.
    """
    return amplitude*np.exp(-(x - center)**2/(2*sigma**2))

def asym_gaussian(x, center, sigma1, sigma2, amplitude):
    """A asymmetric gaussian function composed by two gaussian halves.

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

    Let `g` be a Gaussian function (with apmplitude = 1), the bridge function
    is defined as:

        amplitude * (1 - g(x, center1, sigma1) - g(x, center2, sigma2))

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
    y[mask] = 1.
    y[mask] -= gaussian(x[mask], center1, sigma1)
    y[mask] -= gaussian(x[mask], center2, sigma2)
    y *= amplitude
    return y

def bridge_function2(x, center1, center2, sigma1, sigma2, amplitude):
    """A "bridge" function, geometric complementary of two gaussian peaks.

    Let `g` be a Gaussian function (with apmplitude = 1), the bridge function
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
def factory_gaussian():
    """Return a Gaussian model that can fit data.

    Arguments:
        add_bridge (bool): if True adds a bridge function between the two
            gaussian peaks. If False the model has only two Gaussians.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    model = lmfit.models.GaussianModel()
    model.set_param('center', 0.1)
    model.set_param('sigma', 0.1)
    model.set_param('amplitude', 1)
    return model

def factory_two_gaussians(add_bridge=False, p1_center=0., p2_center=0.5,
                          p1_sigma=0.03, p2_sigma=0.08):
    """Return a 2-Gaussian + (optional) bridge model that can fit data.

    Arguments:
        add_bridge (bool): if True adds a bridge function between the two
            gaussian peaks. If False the model has only two Gaussians.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    peak1 = lmfit.models.GaussianModel(prefix='p1_')
    peak2 = lmfit.models.GaussianModel(prefix='p2_')
    model = peak1 + peak2

    model.set_param('p1_center', p1_center, min=-0.15, max=0.15)
    model.set_param('p2_center', p2_center, min=0.2, max=1)
    model.set_param('p1_sigma', p1_sigma, min=0.02, max=0.2)
    model.set_param('p2_sigma', p2_sigma, min=0.02, max=0.2)
    model.set_param('p1_amplitude', 1, min=0.01)
    model.set_param('p2_amplitude', 1, min=0.01)

    if add_bridge:
        bridge = lmfit.model.Model(bridge_function, prefix='br_')
        model += bridge
        bridge.set_param('amplitude', 0.0001, min=0)
        model.params['br_center1'].expr = 'p1_center'
        model.params['br_center2'].expr = 'p2_center'
        model.params['br_sigma1'].expr = 'p1_sigma'
        model.params['br_sigma2'].expr = 'p2_sigma'
    return model

## lmfit composite model used for fitting
def factory_two_asym_gaussians(add_bridge=False):
    """Return a 2-Asym-Gaussians + (optional) bridge model that can fit data.

    The Asym-Gaussian function is :func:`asym_gaussian`.

    Arguments:
        add_bridge (bool): if True adds a bridge function between the two
            gaussian peaks. If False the model has only two Asym-Gaussians.

    Returns
        An `lmfit.Model` object with all the parameters already initialized.
    """
    peak1 = lmfit.model.Model(asym_gaussian, prefix='p1_')
    peak2 = lmfit.model.Model(asym_gaussian, prefix='p2_')
    bridge = lmfit.model.Model(bridge_function, prefix='br_')

    peak1.set_param('center', 0.0, min=-0.15, max=0.15)
    peak1.set_param('sigma1', 0.03, min=0.02, max=0.5)
    peak1.set_param('sigma2', 0.03, min=0.02, max=0.5)
    peak1.set_param('amplitude', 1, min=0.01)
    peak2.set_param('center', 0.8, min=0.2, max=1)
    peak2.set_param('sigma1', 0.08, min=0.02, max=0.5)
    peak2.set_param('sigma2', 0.08, min=0.02, max=0.5)
    peak2.set_param('amplitude', 1, min=0.01)

    model = peak1 + peak2
    if add_bridge:
        model += bridge
        bridge.set_param('amplitude', 0.0001, min=0)
        model.params['br_center1'].expr = 'p1_center'
        model.params['br_center2'].expr = 'p2_center'
        model.params['br_sigma1'].expr = 'p1_sigma2'
        model.params['br_sigma2'].expr = 'p2_sigma1'
    return model


##
# Classes to perform fit-related operations on multi-channel data
#
class MultiFitterBase(object):
    """Base class for fitting multiple dataset contained in a list.
    """
    def __init__(self, data_list):
        self.data_list = data_list
        self.ndata = len(data_list)
        self.set_weights()

    def set_weights(self, weight_func=None, weight_kwargs=None):
        """Set weights for each data-set.

        To compute the weights fo each dataset the `weight_func` is called
        multiple times. Keys in `weight_kwargs` are arguments of
        `weight_func`. Values in `weight_kwargs` are either scalars, in which
        case they are passed to `weight_func`, or lists. When an argument
        is a list, only one element of the list is passed each time.

        Arguments:
            weight_func (function): function that returns the weights
            weight_kwargs (dict): keywork arguments to be passed to
                `weight_func`.
        """
        if weight_func is None:
            self.weights = [None]*self.ndata
        else:
            self.weights = []
            for i in range(self.ndata):
                weight_kw_i = {k: v[i] if np.size(v) > 1 else v
                                    for k, v in weight_kwargs.items()}
                self.weights.append( weight_func(**weight_kw_i) )

class MultiFitterKDE(MultiFitterBase):
    """
    Performs KDE and fitting of peak position on a list of data populations
    (one population per channel).

    Arguments:
        data_list (list): a list of data populations. For example
            a list of per-channel E or S values (one value per burst).
    """
    def __init__(self, data_list):
        super(MultiFitterKDE, self).__init__(data_list)

    def calc_kde(self, bandwidth=0.03):
        """Compute the list of kde functions and save it in `.kde`.
        """
        self.kde_bandwidth = bandwidth
        self.kde = []
        for data, weights_i in zip(self.data_list, self.weights):
            self.kde.append(
                gf.gaussian_kde_w(data, bw_method=bandwidth,
                                  weights=weights_i))

    def find_kde_max(self, x_kde, xmin=None, xmax=None):
        """Finds the peak position of kde functions between `xmin` and `xmax`.

        Results are saved in the list `.kde_max_pos`.
        """
        if not hasattr(self, 'kde'):
            self.calc_kde()
        self.kde_max_pos = np.zeros(self.ndata)
        for ich, kde in enumerate(self.kde):
            self.kde_max_pos[ich] = \
                    bext.find_max(x_kde, kde(x_kde), xmin=xmin, xmax=xmax)


class MultiFitterModel(MultiFitterBase):
    """
    Performs histogramming and fitting to a model on a list of data
    populations (one population per channel).

    Arguments:
        data_list (list): a list of data populations. For example
            a list of per-channel E or S values (one value per burst).
    """
    def __init__(self, data_list):
        super(MultiFitterModel, self).__init__(data_list)

    def histogram(self, **kwargs):
        """Compute the histogram of the data for each channel.

        All the kwargs are passed to `numpy.histogram`.
        """
        kwargs.update(density=False)
        hist_counts = []
        for data, weights in zip(self.data_list, self.weights):
            if weights is not None:
                kwargs.update(weights=weights)
            counts, bins = np.histogram(data, **kwargs)
            hist_counts.append( counts )
        self.hist_bins = bins
        self.hist_bin_width = (bins[1] - bins[0])
        self.hist_axis = bins[:-1] + 0.5*self.hist_bin_width
        self.hist_counts = np.array(hist_counts)
        self.hist_pdf = np.array(hist_counts, dtype=np.float)
        self.hist_pdf /= self.hist_counts.sum(1)[:, np.newaxis]
        self.hist_pdf /= self.hist_bin_width

    @property
    def model_class(self):
        return self._model_class

    @model_class.setter
    def model_class(self, model_class):
        self._model_class = model_class
        self.model_class_kwargs = {}

    @property
    def model_class_kwargs(self):
        return self._model_class_kwargs

    @model_class_kwargs.setter
    def model_class_kwargs(self, kwargs):
        self._model_class_kwargs = kwargs

    def fit_histogram(self, pdf=True, **fit_kwargs):
        """Fit the histogram of each channel using the same lmfit model.

        A list of `lmfit.Minimizer` is stored in `.fit_obj`.
        The fitted values for all the parameters and all the channels are
        save in a Pandas DataFrame `.fit_params`.

        Arguments:
            model_func (function): when called it returns an lmfit Model with
                all the parameters already initialized.
            pdf (bool): if True fit the normalized histogram (.hist_pdf)
                otherwise fit the raw counts (.hist_counts).

            fit_kwargs (dict): keyword arguments passed to `model().fit`.
        """
        model_class = partial(self.model_class, **self.model_class_kwargs)

        data_list = self.hist_pdf if pdf else self.hist_counts
        self.fit_params = pd.DataFrame(index=range(self.ndata),
                                       columns=model_class().param_names)
        self.fit_obj = []
        for i, data in enumerate(data_list):
            self.fit_obj.append(model_class().fit(data, x=self.hist_axis,
                                **fit_kwargs) )
            self.fit_params.iloc[i] = pd.Series(self.fit_obj[-1].values)

    @property
    def x_axis(self):
        if not hasattr(self, '_x_axis'):
            self._x_axis = np.linspace(self.hist_axis[0],
                                       self.hist_axis[-1], 1000)
        return self._x_axis


def plot_mfit(fitter, ich=0, residuals=False, ax=None):
    """Plot data histogram and fitted model from a `MultiFiter` object.

    Assumes data between 0 and 1 and a two peaks model with parameters
    `p1_center` and `p2_center`.

    Return
        A matplotlib figure object.
    """
    if ax is None:
        fig = plt.figure(figsize=(7, 5.5))
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    ax.plot(fitter.hist_axis, fitter.hist_pdf[ich], '-o', alpha=0.5)
    num_bursts = fitter.data_list[ich].size
    ax.set_title('CH = %d #Bursts %d' % (ich, num_bursts))
    ax.set_xlim(-0.19, 1.19)

    if not hasattr(fitter, 'fit_obj'):
        return fig

    fit_obj = fitter.fit_obj[ich]
    x = fitter.x_axis
    red = '#E41A1C'

    ax.plot(x, fit_obj.model.eval(x=x, **fit_obj.values), 'k', alpha=0.8)
    if  fit_obj.model.components is not None:
        for component in fit_obj.model.components:
            ax.plot(x, component.eval(x=x, **fit_obj.values), '--k',
                    alpha=0.8)
    for param in ['center', 'p1_center', 'p2_center']:
        if param in fitter.fit_params:
            ax.axvline(fitter.fit_params[param][ich], ls='--', color=red)

    if residuals:
        ax.xaxis.set_ticklabels([])
        divider = make_axes_locatable(ax)
        ax_resid = divider.append_axes("bottom", size=1.2, pad=0.1,
                                       sharex=ax)
        ax_resid.plot(fitter.hist_axis,
                      fitter.hist_pdf[ich] - fit_obj.best_fit)
        ax_resid.set_yticks(np.r_[-2:2:0.5])
        ax_resid.set_ylim(-0.49, 0.49)
        ax_resid.set_xlim(-0.2, 1.2)
    return fig


