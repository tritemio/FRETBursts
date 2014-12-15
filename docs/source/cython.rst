.. _fretbursts_cython:

FRETBursts Cython extensions
============================

`Cython <http://cython.org/>`_ is a tool that, among other things, allows
to translate annotated python code into C code.
The C code can be then compiled into a dynamic library and transparently
called from python like any other python library, but with the advantage
of a much higher execution speed.

For some core burst-search functions FRETBursts includes both a pure pyhton
and a cython version. At import time, the code looks for the
compiled version and, if not found, falls back to the pure python version.
Therefore, although the compiled cython version is completely optional,
it allows to gain significant execution speed in core functions that are
potentially executed many times.

Usually the cython extensions are compiled during installation.
To manually build the extensions type::

    python setup.py build

from the FRETBursts source folder.
