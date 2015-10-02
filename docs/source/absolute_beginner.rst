.. _absolute_beginner:

Getting started for the absolute python beginner
================================================

Before running FRETBursts you need to install a python distribution that
includes the Jupyter/IPython Notebook application.

You can find a quick guide for installing the software and running your first
notebook here:

- |jupyter_quick_guide|

.. |jupyter_quick_guide| raw:: html

   <a href="http://jupyter-notebook-beginner-guide.readthedocs.org/" target="_blank">Jupyter/IPython Notebook Quick Start Guide</a>

Be sure you can start Jupyter Notebook application and open a notebook before
moving to the next section.

Installing FRETBursts
---------------------

To install FRETBursts you need to type the following commands in a terminal
(i.e. ``cmd`` on Windows or ``Terminal`` on OSX)::

    conda config --add channels tritemio
    conda install fretbursts

The installation should take a few seconds. If you notice
any error please report it by opening a new issue on the
`FRETBursts GitHub page <https://github.com/tritemio/FRETBursts>`_.

Running FRETBursts tutorial notebook
------------------------------------

Download the ZIP file of
`FRETBursts notebooks <https://github.com/tritemio/FRETBursts_notebooks/archive/master.zip>`__
and extract it inside a folder accessible by the Jupyter Notebook App.

Next, in the new Jupyter Notebook Dashboard click on the folder containing
the FRETBursts notebooks and open the us-ALEX tutorial and follow the
instructions therein.

Remember, to run the notebooks step-by-step (one cell a time) keep pressing
*shift + enter*. To run the entire notebook in a single step click on menu
*Cell -> Run All*.

For more info how to run/edit a notebook see |run_notebook|.

.. |run_notebook| raw:: html

   <a href="http://jupyter-notebook-beginner-guide.readthedocs.org/en/latest/execute.html" target="_blank">Running the Jupyter Notebook</a>
