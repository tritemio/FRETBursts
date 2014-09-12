Installation
============

.. contents ::

FRETBursts is distributed as source code.
The installation consists in installing a scientific python
distribution, downloading FRETBursts sources, and running the Installation
notebook.



Windows Quick Installation steps
---------------------------------

If you are not familiar with Python and Git and you are using MS Windows
follows these steps:

1. Install `Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__
(a scientific python distribution).

2. Install `SourceTree <http://www.sourcetreeapp.com/>`__ (GUI for Git).

3. Download FRETBursts using SourceTree (see below :ref:`download_fretbursts`).

4. Open FRETBursts folder, enter the `notebooks` folder and double-click
"Launch IPython Notebook Here.bat".

5. A new browser windows will open. Click on `FRETBursts Installation`
and run the notebook. Follow the instruction there to run the
tutorials notebooks.

For more info and customizations continue reading.


Installing Python for scientific computing
------------------------------------------

On all the main platforms, the easiest way to install python and all
the scientific packages is using a python distribution like
`Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__ or
`Enthought Canopy <https://www.enthought.com/products/canopy/>`__. Both
distributions include, in the free version, all the needed software (and much
more).

Two dependencies, lmfit and pyside, are not installed by default by
Anaconda and will be installed by FRETBursts Installation notebook.

FRETBursts has been tested on Anaconda 1.9 or newer.

Manual Dependencies
^^^^^^^^^^^^^^^^^^^

If you prefer a manual installation, FRETBursts dependencies are:

- Python 2.7
- Numpy/Scipy (any version from 2013 on)
- Matplotlib with qt (pyside) backend (1.3.x or greater)
- IPython 1.x (2.x recommended)
- PyTables 3.x (optional). To load/save the :ref:`HDF5 Ph-Data <hdf5-format>`.
- lmfit: (version 0.8rc3 or higher) for poweful and user-friendly fitting
- Pandas (optional) currently used only by the fitting framework
- a modern browser (Chrome suggested)

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

Downloading FRETBursts
----------------------

You can download a simple ZIP-ball containing FRETBursts by clickind on
**Download ZIP** on
`FRETBursts Homepage <https://github.com/tritemio/FRETBursts>`__ on GitHub.


However the preferred way is downloading FRETBursts through Git (in other
words "cloning FRETBursts"). In this case copy the **clone URL** is::

    https://github.com/tritemio/FRETBursts.git

When using SourceTree, click on *Clone/New* -> *Clone Repository* and paste
the **clone URL** in *Source Path/URL*. You can choose where to put the
sources.

From the command line, type::

    git clone https://github.com/tritemio/FRETBursts.git


.. _install_fretbursts:

Installing FRETBursts
---------------------

Strictly speaking FRETBursts is not installed as it runs from the folder
where you download it. However some optional dependencies and a configuration
file is created by running the "FRETBursts Installation" notebook that you
find in the notebooks folder.


To run the FRETBursts Installation notebook:

- On windows, click on "Launch IPython Notebook Server Here.bat" (inside the
  notebooks folder)and then click on "FRETBursts Installation".

- On the other platforms::

    cd notebook_folder
    ipython notebook

.. Note ::

    Once the configuration is done, you can load FRETBursts in any notebook
    by running `run load_fretbursts`. Note that you need a copy of the
    `load_fretbursts.py <https://github.com/tritemio/FRETBursts/blob/master/notebooks/load_fretbursts.py>`_
    script in the notebook folder.


Running FRETBursts
------------------

We recommend starting by running the
`tutorial notebooks <https://github.com/tritemio/FRETBursts_notebooks>`__.
The easiest way to perform a new analysis is to modify (or copy) one of the
notebooks.

To run the FRETBursts notebooks,
`download <https://github.com/tritemio/FRETBursts_notebooks/archive/master.zip>`__
and decompress the ZIP-ball in a folder and launch an IPython Notebook server
**inside that folder**. For more details see
[IPython Notebook Startup Folder](http://fretbursts.readthedocs.org/installation.html#ipython-notebook-startup-folder) on FRETBursts documentation.

On the first run, the tutorial notebooks will automatically download
some public datasets of smFRET measurements that are provided for testing
demonstration.

These datasets are free to use for any purposes
(CC0 license). If you use these datasets please cite as:

* Ingargiola, Antonino; Chung, Sangyoon (2014): smFRET example datasets
  for the FRETBursts software. figshare.
  `DOI 10.6084/m9.figshare.1019906 <http://dx.doi.org/10.6084/m9.figshare.1019906>`_


.. _ipython_notebook_startup:

IPython Notebook startup folder
-------------------------------

To use the IPython Notebook you have to launch a local notebook server in
the folder containing the notebooks files (or in a parent folder).

On windows (Anaconda), you can copy and modify the IPython launcher you find in
the start menu. To change the
startup folder right click on the
*IPython Notebook icon* -> *Properties*, and set the new folder
in the *Start in* field.

On all the platforms, you can start IPython Notebook from the terminal
(cmd.exe on Windows) with::

    cd notebook_folder
    ipython notebook

.. Note ::

    The preferred browser is Chrome or Firefox. The use of MS Explorer is
    discouraged as its implementation of web standards is incomplete and not
    compliant.


Installing a compiler (obsolete)
--------------------------------

.. warning ::

    This paragraph is retained for historical reasons and because it may be useful
    for some user. However with recent versions of Anaconda the compiler is
    included and these steps are not necessary anymore.

Some core burst-search core functions can be optionally compiled to gain
significant execution speed. This process requires a compiler to be
installed.

On **Linux** the preferred compiler is GCC, that is easily available for
any distribution.

On **Windows**, the MS Visual Studio compiler is preferred. To install
it search on internet for the files VS2008ExpressWithSP1ENUX1504728.iso
and GRMSDKX\_EN\_DVD.iso.

On **Mac OSX** you should install the LLVM compiler included in Xcode.

*See also:*

* :doc:`cython`

