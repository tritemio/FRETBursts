
Notebook-based workflow
=======================

IPython Notebook is the recommended environment to perform interactive
analysis with FRETBursts.

Tutorials for **FRETBursts** are provided as
<http://ipython.org/notebook.html>`IPython Notebooks`_.

Typically, a new analysis is performed making a copy on an existing notebook
(used as a template) and applying all the needed modifications.

The FRETBursts notebooks are "revision-control aware", meaning that the exact
FRETBursts revision used during each execution is stored (and displayed)
at load time. Saving the software revision together with analysis
commands and results allows long term reproducibility and provides a
lightweight approach for regression testing.

The following sections describe how to configure and use the IPython notebooks
to perform analysis with FRETBursts.


Loading FRETBursts from a notebook
----------------------------------

From inside a notebook, FRETBursts is loaded running a small script
(`load_fretbursts.py`) placed in the notebooks folder. For this reason,
a typical FRETBursts notebook always starts with the line:

::

    %run load_fretbursts

The script switches from notebook folder to FRETBursts source folder
and loads the software.

Before the first execution, you have to tell the script where the FRETBursts
source folder is. You can either paste the folder name in ``load_fretbursts.py`` or
set an environment variable named ``FRETBURSTS_DIR`` containing the path.
