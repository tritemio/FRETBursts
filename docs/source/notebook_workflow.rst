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
`FRETBursts notebooks <https://github.com/tritemio/FRETBursts_notebooks/archive/master.zip>`__
that contains live tutorials ready to run and modify. Just copy the
FRETBursts notebooks in your *ipython notebook folder* and launch
the "IPython Notebook" application.

The following paragraphs explain the advantages of a notebook-based
workflow and provides more detailed informations on how to setup
the IPython Notebook application.

Why a notebook-based workflow
-----------------------------

`IPython Notebooks <http://ipython.org/notebook.html>`__ is the recommended
environment to perform interactive analysis with FRETBursts.

The `FRETBursts tutorials <https://github.com/tritemio/FRETBursts_notebooks>`_
are ipython notebooks and, typically,
a new analysis is performed by copying and modifying an existing notebook.

The FRETBursts notebooks display and store the exact FRETBursts version
(including the revision) used in the execution. Saving the software revision
together with analysis commands and results allows long term reproducibility
and provides a lightweight approach for regression testing.

.. _ipython_notebook_startup:

Starting the IPython Notebook server
------------------------------------

To use the IPython Notebook you need to keep the notebooks in the folder
(or subfolder) where the local notebook server is launched.

On Windows (Anaconda), to change the ipython notebook startup folder,
you can modify the "IPython launcher" found in the start menu. Just copy the
launcher and modify it with: right click on the *IPython Notebook icon* ->
*Properties*, and set the new folder in the *Start in* field.

On all the platforms, you can start IPython Notebook from the terminal
(cmd.exe on Windows) with::

    cd notebook_folder
    ipython notebook

.. Note ::

    The preferred browser is Chrome or Firefox. The use of MS Explorer is
    discouraged as its implementation of web standards is incomplete and not
    compliant.


Loading FRETBursts from source (development)
--------------------------------------------

When running FRETBursts directly
:ref:`from the source folder <source_install>`,
the import is performed by the ``load_fretbursts`` script with::

    %run load_fretbursts --sources

Note that, before running this script for the first time, you need to run
the "Installation notebook" (see :ref:`source_install`).
