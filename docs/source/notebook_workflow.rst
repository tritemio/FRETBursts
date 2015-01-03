Running FRETBursts
==================

Why a notebook-based workflow
-----------------------------

IPython Notebook is the recommended environment to perform interactive
analysis with FRETBursts.

Tutorials for **FRETBursts** are provided as
`IPython Notebooks <http://ipython.org/notebook.html>`__.

Typically, a new analysis is performed making a copy on an existing notebook
(used as a template) and applying all the needed modifications.

The FRETBursts notebooks display and store the exact
FRETBursts revision used during each execution. Saving the software revision
together with analysis commands and results allows long term reproducibility
and provides a lightweight approach for regression testing.

The following sections describe how to configure and use the IPython notebooks
to perform analysis with FRETBursts.


.. _ipython_notebook_startup:

Starting the IPython Notebook server
------------------------------------

To use the IPython Notebook you have to launch a local notebook server in
the folder containing the notebooks files (or in a parent folder).

On windows (Anaconda), you can copy and modify the IPython launcher you find in
the start menu. To change the startup folder right click on the
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


Loading FRETBursts from a notebook
----------------------------------

When running FRETBursts from a :ref:`package installation <package_install>`
you can import ``fretbursts`` and print citation information with::

    from fretbursts import *

This is the syntax used throughout the tutorials.

When running FRETBursts directly
:ref:`from the source folder <source_install>`,
the import is performed by the ``load_fretbursts`` script with::

    %run load_fretbursts

Note that before running this script for the first time you need to run the
:ref:`Installation notebook <install_notebook>`.

