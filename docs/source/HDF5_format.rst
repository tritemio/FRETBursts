HDF5-based smFRET file format
=============================

FRETBursts allows to save and load confocal smFRET data from and to 
an `HDF5 file format <http://www.hdfgroup.org/HDF5/>`_ (see also
`HDF on Wikipedia <http://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_).

For documentation on HDF5-smFRET file format in *FRETBursts* see:

.. toctree::
    :maxdepth: 1

    HDF5_smFRET

For an introduction to HDF5 and why is importand for smFRET data read on.


What is HDF5?
-------------

HDF5 is standard and general-purposes container-format for binary data.
The format can store any number of 
multi-dimensional arrays with no size limit in a hierarchical fashion
(i.e. arrays can be put in folders and subfolders called groups).
Any dataset or folder can have metadata attached to it (for example a 
description, a date, or an ancillary array of parameters, etc...).

The format is self-describing, so any HDF5 compatible application can read
any array without needing to know the type (i.e. int32 or float) or the
byte layout (i.e. big-endian little-endian).

HDF5 supports transparent data compression using the zlib algorithm
or any third-party algorithm via plugins.

The format is an open standard supported by the non-profit organization 
HDFGroup. Open-sources libraries to
read the format are available for all the main programming languages.

The HDF5 ecosystem
------------------

`Numerous organizations <http://www.hdfgroup.org/users.html>`_ use HDF5. 
Just as an example, the native MATLAB format is HDF5-based from version 7.3 on.

Open-source libraries to read the HDF5 format exist for the majority of
programming languages. Among the others FORTRAN, C, C++, C#, Java, MATLAB,
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

Numerous advantages can be easily evisioned:

* **Efficiency**: HDF5 is highly efficient both for space and speed. Libraries
  to interoperate the format are broadly used and heavily tested. We scientist
  don't need to reinvent the wheel, we can (and should) use a state-of-the art 
  (software) technology already available.

* **Long-term persistence:** in 5-10-20 years the data can be always read 
  without relying on an opaque, poorly document, format.

* **Easy interoperability:** a single format lowers the barriers for 
  data-exchange and collaboration. A single format makes easier to compare 
  the output of different analysis software, encourages reproducibility and 
  foster collaboration between different groups.

