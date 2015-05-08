.. _absolute_beginner:

Getting started for the absolute python beginner
================================================

In this section you will find a step-by-step guide on how to install
FRETBursts and its dependencies, and on how to run the tutorial for
the first time. It assumes you have no previous experience with python.
Please refer to the section :ref:`installation` if you already know what
*conda* and *ipython notebook* are.

The step "zero" is installing a modern standard-compliant browser. Either
Mozilla Firefox or Google Chrome will work well. Please try to avoid
MS Explorer.

The step "one" is installing a scientific python distribution
that includes all the dependencies, in this guide we will use
the Anaconda distribution created by Continuum.

- Download `Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`_
  (free version), python 2.7, 64 bits.

- Install it using the default settings for a single user.

The next sections we describe how to configure the Jupyter Notebook application,
how to install FRETBursts, and how to run the tutorial using the Jupyter
Notebook application.

Jupyter Notebook
----------------

Jupyter Notebook (formerly IPython) is an application that runs inside your browser.
For our purposes, Jupyter Notebook is ran locally like a normal desktop
application and does not require accessing any remote server.

The notebook files (for example the FRETBursts tutorials) are text documents
produced by the Jupyter Notebook application that contains both code and
and rich text (paragraph, equations, figures, links, etc...).
Notebooks are both human-readable documents containing the analysis
description and the results (figures, tables, etc..) and executable documents
that can be run to perform a data analysis.

Launching Jupyter Notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~

When started, Jupyter Notebook can access only files within his start-up folder
(including any sub-folder). We need therefore to choose a folder
where we will place all the notebooks and set this folder as the
Jupyter Notebook start-up folder.

See below for platform-specific instructions on how to start Jupyter Notebook
in a specific folder.

On Windows
''''''''''

- Choose a folder that will contain all your notebooks including the
  FRETBursts notebook tutorial.

- Copy the IPython Notebook launcher from the menu to the desktop.

- Right click on the new launcher and change the "Start in" field by pasting
  the folder previously chosen.

- Double-click on the IPython Notebook launcher to start the Jupyter
  Notebook application. Note that this will open two new windows:
  a terminal (used only for error logging and for shut down) and a new
  browser window showing the Jupyter Notebook application.


On OSX
''''''

For simplicity we will start Jupyter Notebook in your home folder.
In this case all your notebooks and FRETBursts tutorial must be placed
in a sub-folder of your home folder (i.e. not in an external hard drive).

To launch Jupyter Notebook:

- Click on spotlight, type `terminal` to open a terminal window.

- Type ``ipython notebook``. This should open a new window in your
  default browser showing the content of your home folder.


Shutting down Jupyter Notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Briefly, to shut down Jupyter Notebook close the associated terminal.

More in details,
the Jupyter Notebook application appears in your browser at a default
address (*http://localhost:8888*).
Closing the browser will not close the Jupyter Notebook application.
You can reopen the previous address (just start typing *localhost*)
and the Jupyter Notebook application will be redisplayed.

You can run many copies of the Jupyter Notebook application and they will show
up at a similar address except that the number (which is the port number)
will increment for each new copy.

For simplicity we do not recommend running multiple copies of Jupyter Notebook.


Installing FRETBursts
---------------------

To install FRETBursts you need to type the following commands in a terminal.

.. note::

  On Windows, open a terminal windows by typing ``cmd`` in the start menu or
  by clicking on *Anaconda Command Prompt*.

.. note::

  On OSX, open a terminal windows by typing ``terminal`` in Spotlight in the
  top-right corner of your desktop.


In a terminal type (hit Enter after each line)::

    conda config --add channels tritemio
    conda install fretbursts

After a few seconds FRETBursts installation is complete. If you notice
any error please report it by opening a new issue on the
`FRETBursts GitHub page <https://github.com/tritemio/FRETBursts>`_.

Running FRETBursts tutorial notebook
------------------------------------

Start by downloading the ZIP file of the
`FRETBursts notebooks <https://github.com/tritemio/FRETBursts_notebooks/archive/master.zip>`__
and extracting it inside your notebook folder (or a sub-folder of it).

Then follow these steps:

- Launch the Jupyter Notebook like explained in the previous section.

- In the new Jupyter Notebook browser window click on the folder containing
  the FRETBursts notebooks and click on the us-ALEX tutorial. This will show
  the tutorial on a new browser tab.

- Click on the menu *Help -> User Interface Tour* to get acquitted with
  the Jupyter Notebook environment.

- You can run the notebooks step-by-step (one cell a time) by hitting
  *shift + enter*.

- You can run the whole notebook in a single step by clicking on the menu
  *Cell -> Run All*.

- To restart the notebook computations (the component performing the
  computation is called the kernel), click on the menu
  *Kernel -> Restart*.

.. note::

    Modifications to the notebooks are automatically saved every
    few minutes. It is suggested that you make a copy of the
    original tutorial (menu *File -> Make a copy ...*) and make
    modifications on the copy.

.. note::

    Closing the browser will not shut down the notebook computational kernel.
    The notebook can be reopened and it will be still running.
    To close a notebook and shut-down the kernel use the menu
    *File -> Close and Halt*. A this point the notebook is closed,
    Jupyter Notebook is still running and can open new notebooks.

.. warning::

    Please pay attention that if you open the **same** notebook on many
    tabs and do edits, the edits on different tabs can overwrite each other.
    To be safe, make sure you open each notebook in only one tab.
    If by mistake you open a notebook twice in two tabs, please close one tab.

Please refer to the `Jupyter Notebook documentation <http://ipython.org/notebook.html>`_
for more information on how to use the Jupyter Notebook environment.
