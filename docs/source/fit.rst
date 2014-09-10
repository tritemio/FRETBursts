.. _fit-section:

Fit framework
=============

For histogram fittings, FRETBursts make use of the `lmfit <>`_ library.
Since lmfit is a pure python package the installation requires just::

    pip install lmfit

FRETBursts requires lmfit version 0.8rc4 or higher.

Fitting E or S histograms
-------------------------

The module :module:`mfit` provides a class :class:`fretbursts.mfit.MultiFitter`
that allow to build histograms and KDE on a multi-channel sample population
(typically E or S values for each burst). The MultiFitter class can find
the max peak position of a KDE or fit the histogram with an abitrary model.
A set of predefined models is provided to handle common cases.
While sensible defaults are applied the user can control
every detail of the fit by setting initial values, parameter bounds
(min, max), algebric constrains and so on. New models can be created by
composing simpler models (by using `+` operator).

A convenience function :func:fretbursts.burstlib_ext.burst_fitter` can be
used to create a `MultiFitter` object to fit either E or S. As an example
let suppose to have a measurent loaded in the variable `d`. To create a
fitter object and compute the FRET histogram we do::

    bext.burst_fitter(d)    # Creates d.E_fitter
    d.E_fitter.histogram()  # Compute the histogram for all the channels

Now to fit with a 2-Gaussians models::

    d.E_fitter.fit_histogram(mfit.factory_two_gaussians)

And to plot the histogram and the model with fitted parameters::

    dplot(d, hist_fret, show_model=True)

More detailed example can be found in the `tutorials <>`_
in notebooks on
`us-ALEX analysis <>`_.

Lmfit introduction
------------------
Lmfit provides a simple and flexible interface for non-linear least squares
and other minimization methods. All the model parameters can be fixed/varied,
have bounds (min, max) or constrained to an algebric expression.

Moreover lmfit provides a Model class and a set of built-in models
that allows to express curve-fitting problems in an compact and expressive
form. Basic models (such as a Gaussian peak) and be composed allowing
an easy definitions of a variety of models (2 or 3 Gaussians).

For more information refer to the official
`lmfit documentation <>`_.


Legacy Fit functions
====================

These are legacy functions used in versions of FRETBursts < 0.4.
This function are retained for backward compatibility.

A set of generic fit function is provided in `fretbursts/fit`.
These are low-level (i.e. generic) fit functions to fit gaussian or
exponential models.

.. toctree::
    :maxdepth: 3

    gaussian_fitting
    exp_fitting
