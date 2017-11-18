.. _running_fretbursts:

Running FRETBursts
==================

After installation, FRETBursts can be imported with::

    from fretbursts import *

that will also import numpy (as `np`) and matplolib.pyplot (as `plt`).
This is the syntax used throughout the tutorials.

Alternatively, you can import FRETBursts in its own namespace
(which is cleaner)::

    import fretbursts as fb

To get started with FRETBursts it is recommended that you download the
`FRETBursts notebooks <https://github.com/OpenSMFS/FRETBursts_notebooks/archive/master.zip>`__
that contains live tutorials ready to run and modify.

Why a notebook-based workflow
-----------------------------

`Jupyter Notebooks <http://jupyter.org/>`__ is the recommended
environment to perform interactive analysis with FRETBursts.

The `FRETBursts tutorials <https://github.com/OpenSMFS/FRETBursts_notebooks>`_
are Jupyter notebooks and, typically,
a new analysis is performed by copying and modifying an existing notebook.

The FRETBursts notebooks display and store the exact FRETBursts version
(including the revision) used in the execution. Saving the software revision
together with analysis commands and results allows long term reproducibility
and provides a lightweight approach for regression testing.

For more information on installing and first steps with Jupyter Notebook
see:

- `Jupyter/IPython Notebook Quick Start Guide <https://jupyter-notebook-beginner-guide.readthedocs.org>`__
