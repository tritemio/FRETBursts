"""
Fit multiple data populations (multi-channel data). Each population/channel
is treated independently.
"""

import numpy as np
import pandas as pd
import lmfit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import fretbursts.fit.gaussian_fitting as gf
import fretbursts.burstlib_ext as bext


## Base function shapes
def _gaussian(x, center, sigma, amplitude=1):
    return amplitude*np.exp(-(x - center)**2/(2*sigma**2))

def _asym_gaussian(x, center, sigma1, sigma2, amplitude):
    left_half = x < center
    sigma = np.zeros_like(x)
    sigma[left_half] = sigma1
    sigma[-left_half] = sigma2
    u = (x - center)/sigma
    return amplitude*np.exp(-u*u)

def _bridge_function(x, center1, center2, sigma1, sigma2, amplitude):
    assert center1 <= center2
    assert amplitude >= 0
    mask = (x > center1)*(x < center2)
    y = np.zeros_like(x)
    y[mask] = 1.
    y[mask] -= _gaussian(x[mask], center1, sigma1)
    y[mask] -= _gaussian(x[mask], center2, sigma2)
    y *= amplitude
    return y


## lmfit composite model used for fitting
def two_gaussians_model(add_bridge=False):
    peak1 = lmfit.models.GaussianModel(prefix='p1_')
    peak2 = lmfit.models.GaussianModel(prefix='p2_')
    bridge = lmfit.model.Model(_bridge_function, prefix='br_')

    peak1.set_param('center', 0.0, min=-0.15, max=0.15)
    peak1.set_param('sigma', 0.03, min=0.02, max=0.5)
    peak1.set_param('amplitude', 1, min=0.01)
    peak2.set_param('center', 0.8, min=0.2, max=1)
    peak2.set_param('sigma', 0.08, min=0.02, max=0.5)
    peak2.set_param('amplitude', 1, min=0.01)

    model = peak1 + peak2
    if add_bridge:
        model += bridge
        bridge.set_param('amplitude', 0.0001, min=0)
        model.params['br_center1'].expr = 'p1_center'
        model.params['br_center2'].expr = 'p2_center'
        model.params['br_sigma1'].expr = 'p1_sigma'
        model.params['br_sigma2'].expr = 'p2_sigma'
    return model


## lmfit composite model used for fitting
def two_asym_gaussians_model(add_bridge=False):
    peak1 = lmfit.model.Model(_asym_gaussian, prefix='p1_')
    peak2 = lmfit.model.Model(_asym_gaussian, prefix='p2_')
    bridge = lmfit.model.Model(_bridge_function, prefix='br_')

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

## Class that performs fit-related operations on multi-channel data
class MultiFitter(object):
    """
    Performs distribution (KDE, histogramming) and fitting operations
    on a list of data populations (one population per channel).

    Arguments:
        data_list (list): a list of data populations. For example
            a list of per-channel E or S values (one value per burst).
    """
    def __init__(self, data_list):
        self.data_list = data_list
        self.ndata = len(data_list)

    def calc_kde(self, bandwidth=0.03, weight_func=None, weight_kwargs=None):
        """Compute the list of kde functions and save it in `.kde`.
        """
        if weight_func is None:
            weights = [None]*self.ndata
        else:
            weights = []
            for i in range(self.ndata):
                weight_kw_i = {k: v[i] if np.size(v) > 1 else v
                                    for k, v in weight_kwargs.items()}
                weights.append( weight_func(**weight_kw_i) )

        self.kde_bandwidth = bandwidth
        self.kde_weights = weights
        self.kde = []
        for data, weights_i in zip(self.data_list, weights):
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

    def histogram(self, **kwargs):
        """Compute the histogram of the data for each channel.

        All the kwargs are passed to `numpy.histogram`.
        """
        kwargs.update(density=False)
        hist_counts = []
        for data in self.data_list:
            counts, bins = np.histogram(data, **kwargs)
            hist_counts.append( counts )
        self.hist_bins = bins
        self.hist_bin_width = (bins[1] - bins[0])
        self.hist_axis = bins[:-1] + 0.5*self.hist_bin_width
        self.hist_counts = np.array(hist_counts)
        self.hist_pdf = np.array(hist_counts, dtype=np.float)
        self.hist_pdf /= self.hist_counts.sum(1)[:, np.newaxis]
        self.hist_pdf /= self.hist_bin_width

    def fit_hist(self, model_func, pdf=True, **kwargs):
        """Fit the histogram of each channel using the same lmfit model.

        A list of `lmfit.Minimizer` is stored in `.fit_res_obj`.
        The fitted values for all the parameters and all the channels are
        save in a Pandas DataFrame `.fit_res`.

        Arguments:
            model (function): when called it returns an lmfit Model with all
                the paramenters already initialized.
            pdf (bool): if True fit the normalized histogram (.hist_pdf)
                otherwise fit the raw counts (.hist_counts).

            kwargs (dict): keyword arguments passed to `model().fit`.


        """
        data_list = self.hist_pdf if pdf else self.hist_counts
        self.fit_res = pd.DataFrame(index=range(self.ndata),
                                    columns=model_func().param_names)
        self.fit_res_obj = []
        for i, data in enumerate(data_list):
            self.fit_res_obj.append(model_func().fit(data, x=self.hist_axis,
                                    **kwargs) )
            self.fit_res.iloc[i] = pd.Series(self.fit_res_obj[-1].values)

def plot_mfit(fitter, ich=0, residuals=False, ax=None):
    if not hasattr(plot_mfit, 'x'):
        plot_mfit.x = np.linspace(fitter.hist_axis[0],
                                  fitter.hist_axis[-1], 1000)
    x = plot_mfit.x
    fit_obj = fitter.fit_res_obj[ich]
    red = '#E41A1C'

    if ax is None:
        fig = plt.figure(figsize=(7, 5.5))
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    ax.plot(fitter.hist_axis, fitter.hist_pdf[ich], '-o', alpha=0.5)
    ax.plot(x, fit_obj.model.eval(x=x, **fit_obj.values), 'k', alpha=0.8)
    for component in fit_obj.model.components:
        ax.plot(x, component.eval(x=x, **fit_obj.values), '--k', alpha=0.8)
    ax.axvline(fitter.fit_res['p1_center'][ich], ls='--', color=red)
    ax.axvline(fitter.fit_res['p2_center'][ich], ls='--', color=red)
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(0, 5)
    num_bursts = fitter.data_list[ich].size
    ax.set_title('CH = %d #Bursts %d' % (ich, num_bursts))

    if residuals:
        ax.xaxis.set_ticklabels([])
        divider = make_axes_locatable(ax)
        ax_resid = divider.append_axes("bottom", size=1.2, pad=0.1, sharex=ax)
        ax_resid.plot(fitter.hist_axis, fitter.hist_pdf[ich] - fit_obj.best_fit)
        ax_resid.set_yticks(np.r_[-2:2:0.5])
        ax_resid.set_ylim(-0.49, 0.49)
        ax_resid.set_xlim(-0.2, 1.2)
    return fig


