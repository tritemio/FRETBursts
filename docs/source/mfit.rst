.. currentmodule:: fretbursts.mfit

MultiFitter reference documentation
===================================

.. automodule:: fretbursts.mfit

.. contents::


The MultiFitter class
---------------------

.. autoclass:: MultiFitter
    :members:

Model factory functions
-----------------------

In this section you find the documentation for the factory-functions
that return pre-initialized models for fitting E and S data.

.. autofunction:: factory_gaussian

.. autofunction:: factory_asym_gaussian

.. autofunction:: factory_two_gaussians

.. autofunction:: factory_two_asym_gaussians

.. autofunction:: factory_three_gaussians


Utility functions
-----------------

The following functions are utility functions used to build the
the model functions (i.e. the "factory functions") for the fitting.

.. autofunction:: bridge_function

.. autofunction:: asym_gaussian
