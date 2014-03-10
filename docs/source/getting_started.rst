
Getting started
===============

Installation
------------

FRETBursts is distributed as source code. To run FRETBursts you need to

-  install the dependencies
-  download FRETBursts sources from GitHub
-  execute the notebooks in the notebook dir

These 3 steps are described below.

Install python dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~

On all the main platforms, the preferred way to install python and all
the scientific packages is using a python distribution like
`Anaconda <https://store.continuum.io/cshop/anaconda/>`__ or
`Canopy <https://www.enthought.com/products/canopy/>`__.

By installing a python distribution all the python dependencies are
fulfilled.

FRETBursts is tested on Anaconda 1.9 or newer.

    **List of python dependencies:**

    For those interested the list of used python packages is: numpy,
    scipy, matplotlib, IPython, pytables (optional), cython (optional).

    Unit tests are run with `py.test <http://pytest.org/latest/>`__.

    The documentation is built using
    `Sphinx <http://sphinx-doc.org/>`__.

Installing a compilers (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some core burst-search core functions can be optionally compiled to gain
significant execution speed. This process requires a compiler to be
installed.

On **Linux** the preferred compiler is GCC, that is easily available for
any distribution.

On **Windows**, the MS Visual Studio compiler is preferred. To install
it search on internet for the files VS2008ExpressWithSP1ENUX1504728.iso
and GRMSDKX\_EN\_DVD.iso.

On **Mac OSX** you should install the LLVM compiler included in Xcode
(untested).

Installing a revision control system (suggested)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FRETBursts uses `Git <http://git-scm.com/>`__ as revision control
system. Installing Git is suggested for all users, because FRETBursts
notebooks can keep track of the software revision. Furthermore, Git will
make easier to download any future update.

Unless you are familiar with Git we suggest to install a graphical
interface like `SourceTree <http://www.sourcetreeapp.com/>`__.

On **Windows**, install SourceTree and, when asked, select the
single-user installation and to download the embedded Git.

On **Mac OSX**, install SourceTree and configure it to use a system-wide
Git installation. Git can be installed system-wide using the
`homebrew <http://brew.sh/>`__ package manager.

On **Linux** Git is easily installed and usually comes with the **gitk**
graphical interface.

Downloading FRETBursts
~~~~~~~~~~~~~~~~~~~~~~

Once Git is installed you can download the latest FRETBursts version by
copying the **clone URL** and using either SourceTree or the command
line.

From SourceTree, click on *Clone/New* -> *Clone Repository* and paste
the **clone URL** in *Source Path/URL*. You can choose where to put the
sources

From the command line, type (FRETBursts will be downloaded in the
current folder):

::

    git clone clone_URL

Download the data-samples
~~~~~~~~~~~~~~~~~~~~~~~~~

Data samples to run FRETBursts are provided separately from the source.
You can download the datasets `here <>`__.

Notebook-based workflow
-----------------------

The user guide for **FRETBursts** is provided in the form of IPython
Notebooks. These notebooks are a sort of enhanced scripts that mix
(live) code, execution results and a rich HTML description in one single
document. Using an IPython Notebook a single document can contain the
analysis code, interleaved with descriptive content and analysis
results.

The preferred way to use the software is copy one of the provided
notebooks and execute and modify it to perform the desired analysis.

The FRETBursts is "revision-control aware", meaning that the exact
FRETBursts revision used during each execution is stored (and displayed)
at load time. Saving the software revision together with analysis
commands and results allows long term reproducibility and provides a
lightweight approach for regression testing.

Configure the notebook workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A typical FRETBursts notebook starts with the line:

::

    %run load_fretbursts

This command switched from notebook folder to FRETBursts source folder
and loads the software.

We must tell the script where the FRETBursts folder is. You can either
paste the folder name in ``load_fretbursts.py`` or set an environment
variable ``FRETBURSTS_DIR`` containing the path.
