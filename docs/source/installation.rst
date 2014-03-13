Installation
============

FRETBursts is distributed as source code. To run FRETBursts you need to

-  install the dependencies
-  download FRETBursts sources from GitHub
-  execute the notebooks in the notebook dir

These 3 steps are described below.


Installing python
-----------------

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

Installing Git (optional)
-------------------------

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
----------------------

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


Downloading the data-samples
----------------------------

Data samples to run FRETBursts are provided separately from the source.
You can download the datasets **here (TODO)**.


Installing a compiler (optional)
--------------------------------

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