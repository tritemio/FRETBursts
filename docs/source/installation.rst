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

Install FRETBursts in new environment
----------------------------------

As a rule, all new development takes place on separate branches. In priniciple, the master branch is always stable and releasable. The advantage of installing from the master branch is that you can easily get immidiate updates without waiting for a formal release. If there are some errors you can always roll back to the latest released version to get your job done. Since you have the full version down to the commit level printed in the notebook you will know which version works and which does not.

Another even safer approach is creating a separate environment with conda. You can create a new environment in which you install a
totally different set of packages. You can have an environment with the latest released FRETBursts and one with the latest master version.
    
Use conda-forge to install a new environment for FRETBursts::

    conda config --append channels conda-forge 
    
This appends the conda-forge channel to the conda channels. 

The default environment is called root. To make a new environment called fbmaster:

    conda create -n fbmaster python=3.6 frebursts
    
This installs python 3.6 and names the environment 'fbmaster.'

Then activate the environment:

    activate fbmaster

See `Managing environments <https://conda.io/docs/using/envs.html>`__ for documentation. 

Once the environment is activated you can install/remove more packages in it. For example you can replace the stable FRETBursts with the version from github master using pip install -e . in the same terminal where the environment has been activated. Installing the stable FRETBursts first, allows installing all the dependencies through conda. Conda adds all the environment to the notebook menu. So when you open a notebook, you can go to the menu Kernel -> Change kernel and select fbmaster instead of default (or vice versa). The latest used kernel is saved in the notebook so you don't have to switch every time.

Environments help to be more reproducible in computations. They can be "saved" or exported to a text file so that they can be recreated on a different machine. For example, you can have an environment for a paper. 1 or 2 years later, an updated python installation may have some incompatibilities causing the old analysis to fail or to give different results. If you saved the environment file, you can restore the old environment with the exact version of all packages with one command. Having a saved environment saves you the trouble of finding and fixing these incompatibilities (they are almost always trivial but it takes time).
