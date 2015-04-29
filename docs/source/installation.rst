.. _installation:

FRETBursts Installation
=======================

FRETBursts can be installed as a standard python package either via `conda`
or PIP (see below). Being written in python, FRETBursts runs on Mac OSX,
Windows and Linux.

In order to easily install the different dependencies, you need to
install a **scientific python distribution**. If you don't know where
to start just install
`Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`_
(free version).

After installing FRETBursts, it is recommended that you download and run the
`FRETBursts notebooks <https://github.com/tritemio/FRETBursts_notebooks/archive/master.zip>`__
to get started. See :ref:`running_fretbursts`.

For updates the latest version please refer to the
`Release Notes (What's new?) <https://github.com/tritemio/FRETBursts/releases>`_.

.. _package_install:

Package installation (conda)
----------------------------

`conda` is the package manager used by Anaconda that simplifies installation
and upgrades. If you installed the Anaconda distribution, this is
the preferred method to install FRETBursts and keep it updated.

As a first step, add the "channel" containing the FRETBursts
packages::

    conda config --add channels tritemio

Then install FRETBursts with::

    conda install fretbursts

The previous command will upgrade to the newest FRETBursts version
whenever available. See also the
`Release Notes (What's new?) <https://github.com/tritemio/FRETBursts/releases>`_
to stay updated on new features and changes.

Package installation (PIP)
--------------------------

FRETBursts can be also installed through the standard python package
manager (PIP).

In this case, make sure that all the non-pure python dependencies are already
installed (i.e. numpy, scipy, pandas, matplotlib, pyqt, pytables).

Then install or upgrade FRETBursts with::

    pip install fretbursts --upgrade

Use the previous command every time you want to upgrade to the latest stable
release.

Install from GitHub (development)
---------------------------------

You can install the latest development version directly from GitHub with::

    pip install git+git://github.com/tritemio/FRETBursts.git

.. note ::
    Note that the previous command fails when `git <http://git-scm.com/>`_
    is not installed.

Alternatively you can clone FRETBursts git repository and from the
source folder run::

    python setup.py build
    python setup.py install

The optimized C extensions are installed in both cases. Make sure that
the dependencies `lmfit` and `seaborn` have been installed.

.. _source_install:

Source installation (development)
---------------------------------

To perform a "source installation", i.e. executing FRETBursts from the source
folder, download the code from GitHub, and execute the ipython notebook
`FRETBursts Installation`.

In the case of source installation, instead of the normal import, we need to
use a little helper script (:ref:`load_fretbursts.py <load_fretbursts>`) that
finds the sources and import FRETBursts.
Once the script `load_fretbursts.py` is copied in your notebook
folder you can load/import FRETBursts with::

    %run load_fretbursts --nogui -- sources

A copy of `load_fretbursts.py` can be found under `notebooks` in the
source tree.
