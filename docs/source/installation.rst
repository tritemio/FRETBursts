.. _installation:

FRETBursts Installation
=======================

FRETBursts can be installed as a standard python package either via `conda`
or PIP (see below). Being written in python, FRETBursts runs on OS X,
Windows and Linux.

For updates on the latest FRETBursts version please refer to the
`Release Notes (What's new?) <https://github.com/tritemio/FRETBursts/releases>`__.

.. _package_install:

Installing latest stable version
--------------------------------

The preferred way to to install and keep FRETBursts updated is through
`conda`, a package manager used by Anaconda scientific python distribution.
If you haven't done it already, please install the python3 version of
`Continuum Anaconda distribution <https://www.continuum.io/downloads>`__
(legacy python2 works too but is less updated).
Then, you can install or upgrade FRETBursts with::

    conda install fretbursts -c conda-forge

After the installation, it is recommended that you download and run the
`FRETBursts notebooks <https://github.com/tritemio/FRETBursts_notebooks/archive/master.zip>`__
to get familiar with the workflow. If you don't know what a Jupyter Notebooks is
and how to launch it please see:

* `Jupyter/IPython Notebook Quick Start Guide <http://jupyter-notebook-beginner-guide.readthedocs.org/en/latest/>`__

See also the FRETBursts documentation section: :ref:`running_fretbursts`.

Alternative methods: using PIP
------------------------------

Users that prefer using `PIP <https://pypi.python.org/pypi/pip>`__, have to
make sure that all the non-pure python dependencies are properly installed
(i.e. numpy, scipy, pandas, matplotlib, pyqt, pytables), then use the
usual::

    pip install fretbursts --upgrade

The previous command installs or upgrades FRETBursts to the latest stable release.


.. _source_install:

Install latest development version
----------------------------------

As a rule, all new development takes place on separate "feature branches".
The master branch should always be stable and releasable.
The advantage of installing from the master branch is that you can
get updates without waiting for a formal release.
If there are some errors you can always roll back to the latest
released version to get your job done. Since you have the full version
down to the commit level printed in the notebook you will know which version
works and which does not.

You can install the latest development version directly from GitHub with::

    pip install git+git://github.com/tritemio/FRETBursts.git

.. note ::
    Note that the previous command fails if `git <http://git-scm.com/>`__
    is not installed.

Alternatively you can do an "editable" installation, i.e. executing
FRETBursts from the source folder. In this case, modifications in the source
files are immediately available on the next FRETBursts import.
To do so, clone FRETBursts and install it as follows::

    git clone https://github.com/tritemio/FRETBursts.git
    cd FRETBursts
    pip install -e .

It is recommended that you install `cython <http://cython.org/>`__ before
FRETBursts so that the optimized C routines are installed as well.
Also, make sure you have `lmfit` and `seaborn` installed before running
FRETBursts.


Install FRETBursts in a separate environment
--------------------------------------------

If you want to install multiple versions of FRETBursts, you can create separate
`environments with conda <https://conda.io/docs/using/envs.html>`__.
Each conda environments can contain
a totally different set of packages, so you can have an environment with the
latest released FRETBursts and one with the latest master version, for example.

FRETBursts is not in the `default conda channel <https://docs.continuum.io/anaconda/pkg-docs>`__,
but in the `conda-forge channel <https://conda-forge.github.io/ >`__.
You can add conda-forge to the channel list with::

    conda config --append channels conda-forge

This **appends** `conda-forge` to the channel list, with a lower
priority than the default channel. It means that a package available,
with the same version, in both conda-forge and the default channel,
will be installed from default.

To make a new environment called `fbmaster` containing python 3.6 and
fretbursts::

    conda create -n fbmaster python=3.6 frebursts

The environment needs to be activated::

    . activate fbmaster

(on windows remove the leading "dot").

Once the environment is activated you can install/remove more packages in it.
For example you can replace the stable FRETBursts with the version from github master using
``pip install -e .`` in the same terminal where the environment has been activated.
Installing the stable FRETBursts first allows installing all the dependencies through conda.
Conda adds the environment to the notebook menu. So when you open a notebook, you can go to the
menu *Kernel* -> *Change kernel* and select *fbmaster* instead of default (or vice versa).
The latest used kernel is saved in the notebook so you don't have to switch every time.

Environments help to be more reproducible in computations. They can be "saved"
or exported to a text file for recreation on a different machine. For example,
you can have the analysis for an old paper that fails to run or gives different
results on an updated python installation. If you saved the environment file,
you can restore the old environment with the exact version of all packages.
It saves you time, trouble and makes the analysis more reproducible.

Refer to the conda documentation
`Managing environments <https://conda.io/docs/using/envs.html>`__ for details.
