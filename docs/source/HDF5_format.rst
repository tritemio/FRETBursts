.. _hdf5-advantages:

Why an HDF5-based smFRET file format
====================================

In this page we briefly introduce what the HDF5 format is and why it is
important for single-molecule FRET data.

What is HDF5?
-------------

`HDF5 <http://www.hdfgroup.org/HDF5/>`_ is standard and general-purposes
container-format for binary data (see also
`HDF on Wikipedia <http://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_).
The format can store any number of
multi-dimensional arrays with no size limit in a hierarchical fashion
(i.e. arrays can be put in folders and subfolders called groups).
Any dataset or folder can have metadata attached to it (for example a
description, a date, or an array of parameters).

The format is self-describing, so any HDF5 compatible application can read
the file content without knowing in advance the data-type (i.e. int32 or float)
or the byte layout (i.e. big-endian little-endian).

HDF5 supports transparent data compression using the zlib algorithm
or any third-party algorithm via plugins.

The format is an open standard supported by the non-profit organization
HDFGroup. Open-sources libraries to
read the format are available for all the main programming languages.

The HDF5 ecosystem
------------------

`Numerous organizations <http://www.hdfgroup.org/users.html>`_ use HDF5.
Just as an example, the native MATLAB format (.mat) is HDF5-based from
version 7.3 on.

Libraries to read the HDF5 format exist for the majority of
programming languages. Among the others, FORTRAN, C, C++, C#, Java, MATLAB,
Python, Mathematica, R have first-class support for the format.

LabView can read/write the format using
`h5labview <http://h5labview.sourceforge.net/>`_.

Origin natively support HDF5 from version 8.1.

Open-source and multi-platform viewers/editors are also available
(see `HDFView <http://www.hdfgroup.org/products/java/hdfview/index.html>`_ and
`ViTables <http://vitables.org/>`_).

Python, in particular, has 2 libraries that allow handling HDF5 files:

* `h5py <http://www.h5py.org/>`_
* `pytables <http://www.pytables.org/>`_

FRETBursts uses **pyTables**.

Why HDF5 and smFRET?
--------------------

Most of smFRET data around the world is acquired through a custom setup and
custom software. As a result the number of file formats is almost as large
as the number of existing setups.

A single, space-efficient and self-documenting file format like HDF5 is
highly preferable to the Babel of formats used today.

Numerous advantages can be easily envisioned:

* **Efficiency**: HDF5 is highly efficient both for space and speed. Libraries
  to interoperate with the format are broadly used and heavily tested.
  Scientists don't need to reinvent the wheel and can leverage the already
  available state-of-the art software technologies.

* **Long-term persistence:** in 5-10-20 years the data can be always read
  without relying on obscure, poorly document, (or in some case vendor
  specific) binary formats.

* **Easy interoperability:** a single format lowers the barriers for
  data-exchange and collaboration. A single format makes easier to compare
  the output of different analysis software, encourages reproducibility and
  foster collaboration between different groups.

HDF5 in FRETBursts
------------------

FRETBursts allows saving and loading smFRET data from and to
an HDF5-based file format called **Photon-HDF5**.

The **Photon-HDF5** is a pre-defined layout to be used with
smFRET and other data involving time-series of photon-data.

A description of the Photon-HDF5 format and its specifications can be found in
`Photon-HDF5 format <http://photon-hdf5.readthedocs.org/>`_.

For documentation on using the Photon-HDF5 format in *FRETBursts* see:

.. toctree::
    :maxdepth: 1

    HDF5_smFRET
