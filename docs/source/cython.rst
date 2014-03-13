Compiling Cython code
=====================

`Cython <http://cython.org/>`_ is a tool that, among other things, allows 
to translate python code into C code.
The C code can be then compiled into a dynamic library an transparently 
called from python like any other python library, but with the advantage
of a much higher execution speed.

For some core burst-search functions FRETBursts includes both a pure pyhton 
and a cython version. At import time, the code looks for the 
compiled version and, if not found, falls back to the pure python version.
Therefore, although the compiled cython version is completely optional,
it allows to gain significant execution speed in core functions that are 
potentially executed many times.

To compile the cython functions into machine code, enter the folder
`burstsearch` as type::

    python setup.py build_ext --inplace

This command requires that both `cython <http://cython.org/>`_ and a standard
compiler are installed.

Installing a compiler (optional)
--------------------------------

Linux
~~~~~

On **Linux** the preferred compiler is GCC, that is easily available for
any distribution.

Windows
~~~~~~~

On **Windows**, the MS Visual Studio v7 compiler is the preferred compiler. 
In principle, installing only MS Windows SDK v7.0 (GRMSDKX\_EN\_DVD.iso) should be sufficient.
However, sometime it maybe necessary to install the full MS Visual Studio 
Express 2008 (VS2008ExpressWithSP1ENUX1504728.iso). Both downloads are free.

After installation and a reboot, the cython version installed by python 
distributions like Anaconda should automatically find the compiler.

Mac OSX
~~~~~~~

On **Mac OSX** the LLVM compiler included in Xcode should be installed
(untested).