Installation
============

FRETBursts can be installed as a standard python package or can be executed
from the source folder.

In both cases, installing a scientific python distribution is needed.


.. _package_install:

Quick: Package installation
---------------------------

If you are just starting with FRETBursts this is the preferred installation
method.

To install the last stable release type::

    pip install fretbursts==0.4rc4 --upgrade

To install the latest version from GitHub, type::

    pip install git+git://github.com/tritemio/FRETBursts.git

Alternatively you can clone FRETBursts git repository and run::

    python setup.py build
    python setup.py install

The optimized C extensions are installed in both cases.

After installation, FRETBurst can be imported with::

    from fretbursts import *

this will also import numpy (as `np`) and matplolib.pyplot (as `plt`).
Alternatively you can import FRETBursts in its own namespace with::

    import fretbursts as fb

For more info on running FRETBursts refer to the
`tutorials <https://github.com/tritemio/FRETBursts_notebooks>`_.


.. _source_install:

Quick: Source installation
--------------------------

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


Source installation description (long)
======================================

Windows: Source installation with SourceTree
---------------------------------------------

If you are not familiar with Python and Git and you are using MS Windows
follow these steps:

1. Install `Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__
(a scientific python distribution).

2. Install `SourceTree <http://www.sourcetreeapp.com/>`__ (GUI for Git).

3. Download FRETBursts using SourceTree (see below :ref:`download_fretbursts`).

4. Open FRETBursts folder, enter the `notebooks` folder and double-click
"Launch IPython Notebook Here.bat".

5. A new browser windows will open. Click on `FRETBursts Installation`
and run the notebook. Follow the instruction there to run the
tutorial notebooks.

For more info and customizations continue reading.


Installing Python for scientific computing
------------------------------------------

On all the main platforms, the easiest way to install python and all
the scientific packages is using a python distribution like
`Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__ or
`Enthought Canopy <https://www.enthought.com/products/canopy/>`__. Both
distributions include, in the free version, all the needed software (and much
more).

A pure-python dependency, lmfit, is not in Anaconda and is installed
during FRETBursts python package installation
or when running the Installation notebook.

FRETBursts has been tested on Anaconda 1.9 or newer.

Manual Dependencies
^^^^^^^^^^^^^^^^^^^

If you prefer a manual installation, FRETBursts dependencies are:

- Python 2.7
- Numpy/Scipy (any version from 2013 on)
- Matplotlib (1.3.x or greater) with QT4 backend (either PyQT4 or PySide).
- IPython 1.x (2.x recommended)
- PyTables 3.x (optional). To load/save the :ref:`HDF5 Ph-Data <hdf5-format>`.
- lmfit: (version 0.8 or higher) for powerful and user-friendly fitting
- Pandas (optional) currently used only by the fitting framework
- a modern browser (Chrome or Firefox suggested)

For developing FRETBursts you should also install

- sphinx 1.2.2 with the napoleon extension (sphinxcontrib-napoleon)
  to build this documentation.
- pytest to execute the unit tests.

Installing Git (optional)
-------------------------

FRETBursts uses `Git <http://git-scm.com/>`__ as revision control
system. Even if not necessary, we strongly recommend installing it because
FRETBursts notebooks can keep track of the FRETBursts software revision.
Furthermore, Git will make easy downloading future updates.

Unless you are familiar with Git it is preferable to install a graphical
interface like `SourceTree <http://www.sourcetreeapp.com/>`__.

On **Windows**, install SourceTree and, when asked, select the
single-user installation and choose to download the embedded Git.
Alternatively, for an independent system-wide Git installation,
download the windows binaries from the
`Git Homepage <http://git-scm.com/downloads>`__.

On **Mac OSX**, install SourceTree and configure it to use a system-wide
Git installation. Git can be installed system-wide using the
`homebrew <http://brew.sh/>`__ package manager.

On **Linux** Git is easily installed and usually comes with the **gitk**
graphical interface.


.. _download_fretbursts:

Obtaining FRETBursts sources
----------------------------

You can download a simple ZIP-ball containing FRETBursts by clicking on
**Download ZIP** on
`FRETBursts Homepage <https://github.com/tritemio/FRETBursts>`__ on GitHub.

However the preferred way is downloading FRETBursts through Git by
"cloning FRETBursts" (you will download the full history).
To clone the FRETBursts from the command line type::

    git clone https://github.com/tritemio/FRETBursts.git

When using SourceTree, click on *Clone/New* -> *Clone Repository* and paste
the `clone URL <https://github.com/tritemio/FRETBursts.git>`_
in *Source Path/URL*.


.. _install_notebook:

Configure FRETBursts to run from the source tree
------------------------------------------------

To run FRETBursts from the source folder (instead of installing the python
package) you first need to run the "FRETBursts Installation" notebook
that will create a configuration file (storing the sources path) and
install some dependencies.

To run the FRETBursts Installation notebook:

- On windows, click on "Launch IPython Notebook Server Here.bat" (inside the
  notebooks folder) and then click on "FRETBursts Installation".

- On the other platforms::

    cd notebook_folder
    ipython notebook

.. Note ::

    Once the configuration is done, you can load FRETBursts in any notebook
    by running `%run load_fretbursts`. Note that you need a copy of the
    `load_fretbursts.py <https://github.com/tritemio/FRETBursts/blob/master/notebooks/load_fretbursts.py>`_
    script in the notebook folder.


C compiler: manual installation
-------------------------------

Some core FRETBursts functions have a :ref:`cython version <fretbursts_cython>`
for higher execution speed. The cython functions require a C compiler that is
already installed when installing the Anaconda distribution.

The following paragraphs may be useful for users wanting to
manually install a C compiler.

On **Linux** the preferred compiler is GNU GCC, that is already installed (or
easily installed) in all the major distributions.
On **Windows**, the MS Visual Studio compiler is preferred. To install
it search on internet for the files VS2008ExpressWithSP1ENUX1504728.iso
and GRMSDKX\_EN\_DVD.iso.

On **Mac OSX** you should install the LLVM compiler included in Xcode.

