Installation
============

FRETBursts is currently distributed as source code.
The installation consist in installing a scientific python
distribution, downloading FRETBursts sources, and setting a folder
for the FRETBursts notebooks.

These 3 steps are described below.


Installing Python for scientific computing
------------------------------------------

On all the main platforms, the easiest way to install python and all
the scientific packages is using a python distribution like
`Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__ or
`Enthought Canopy <https://www.enthought.com/products/canopy/>`__.

By installing a python distribution all the python dependencies are
automatically met.

FRETBursts has been tested on Anaconda 1.9 or newer.

###Dependencies

If you prefer a manual installation, FRETBursts dependencies are:

 - Python 2.7
 - Numpy/Scipy (any version from 2013 on)
 - Matplotlib with qt (pyside) backend (1.3.x or greater)
 - IPython 1.x (2.x suggested)
 - PyTables 3.x (optional)
 - a modern browser (Chrome suggested)

For developing FRETBursts you should also install

 - sphinx 1.2.2 with the napoleon extension (sphinxcontrib-napoleon) to build the docs
 - pytest to execute the unit tests


Installing Git (optional)
-------------------------

FRETBursts uses `Git <http://git-scm.com/>`__ as revision control
system. Even if not necessary, we strongly recommend installing it because
FRETBursts notebooks can keep track of the software revision.
Furthermore, Git will make easier to download any future update.

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


Configure the IPython Notebook
------------------------------

To use the ipython notebook you have to launch a local notebook server from the folder
that contains the notebooks files.

From the terminal:

::
cd my_folder/FRETBursts_notebooks
ipython notebook

On windows, you can configure an IPython launcher by copying the IPython Notebook
icon and modifying the Start in filed. In detail: right click on the
*IPython Notebook icon* -> *Properties* and paste the notebook folder in the
*Start in* field.


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