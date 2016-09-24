.. _installation:

FRETBursts Installation
=======================

FRETBursts can be installed as a standard python package either via `conda`
or PIP (see below). Being written in python, FRETBursts runs on OS X,
Windows and Linux.

For updates on the latest FRETBursts version please refer to the
`Release Notes (What's new?) <https://github.com/tritemio/FRETBursts/releases>`_.

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

Users that prefer to use PIP (the standard python package manager), have to
make sure that all the non-pure python dependencies are properly installed
(i.e. numpy, scipy, pandas, matplotlib, pyqt, pytables), then use the
usual::

    pip install fretbursts --upgrade

The previous installs or upgrades FRETBursts to the latest stable release.


.. _source_install:

Install latest development version
----------------------------------

You can install the latest development version directly from GitHub with::

    pip install git+git://github.com/tritemio/FRETBursts.git

.. note ::
    Note that the previous command fails if `git <http://git-scm.com/>`_
    is not installed.

Alternatively you can clone FRETBursts git repository and run from the
source folder the following commands::

    python setup.py build
    pip install .

The optimized C extensions are installed in both cases. Make sure that
the dependencies `lmfit` and `seaborn` have been installed.

Note that to do an "editable" or "source" installation, i.e. executing
FRETBursts from the source folder you need to add ``-e`` to the lasted command::

    pip install . -e

In this case, modifications in the source files would be immediately available
on the next FRETBursts import.
