
Notebook-based workflow
=======================


The user guide for **FRETBursts** is provided in the form of IPython
Notebooks. These notebooks are a sort of enhanced scripts that mix
(live) code, execution results and a rich HTML description in one single
document. Using an IPython Notebook a single document can contain the
analysis code, interleaved with descriptive content and analysis
results.

The preferred way to use the software is copy one of the provided
notebooks and execute and modify it to perform the desired analysis.

The FRETBursts is "revision-control aware", meaning that the exact
FRETBursts revision used during each execution is stored (and displayed)
at load time. Saving the software revision together with analysis
commands and results allows long term reproducibility and provides a
lightweight approach for regression testing.

Configure the notebook workflow
-------------------------------

A typical FRETBursts notebook starts with the line:

::

    %run load_fretbursts

This command switches from notebook folder to FRETBursts source folder
and loads the software.

We must tell the script where the FRETBursts folder is. You can either
paste the folder name in ``load_fretbursts.py`` or set an environment
variable ``FRETBURSTS_DIR`` containing the path.
