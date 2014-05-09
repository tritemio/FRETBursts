FRETBursts
==========

Overview
--------

**FRETBursts** is an open-source toolkit for analysis of single-molecule FRET
data acquired by single and multi-spot confocal systems.

FRETBursts aims to be a reference implementation for state-of-the-art
algorithms commonly used in smFRET burst analysis.

As input data, both single laser excitation and 2-laser Alternating Excitation
(ALEX) are supported.

Several background estimation and FRET efficiency fitting routines are
implemented. The burst search is an efficient version of the classical
sliding-window search and can be applied on different selection of timestamps
(all, d-only, a-only, etc...). A large number of post-search burst selection
criteria are included and ready to use (see `select_bursts.py`). Moreover,
defining a new burst selection criterium requires only a couple of lines of code.

A variety of preset plotting functions are already defined. Just to mention a few:
time-traces, scatter-plots of any burst data (i.e. E vs S, E vs size, etc...),
histograms (FRET, stoichiometry, inter-photon waiting times, 2D ALEX histogram, etc..),
kernel density estimations and much more (see `burst_plot.py`).
Thanks to the excellent [Matplotlib](http://matplotlib.org/) library,
FRETBursts can produce publication-quality plots out of the box.

FRETBursts can load one of the sample datasets (soon to be released) or any arbitrary
binary timestamps data, providing a suitable loader function. Writing a
loader function is extremely easy thanks to the
[binary loading capabilities](http://docs.scipy.org/doc/numpy/reference/routines.io.html)
of Numpy.

For bug reports please use the GitHub issue tracker. Also, fixes and/or enhancements
are welcome: just send a [pull request (PR)](https://help.github.com/articles/using-pull-requests).

For more info contact me at tritemio @ gmail.com.

Software Environment
--------------------

FRETBursts is written in the [python programming language](http://www.python.org/) using the standard
scientific stack of libraries (numpy, scipy, matplotlib).

Usage examples are provided as IPython notebooks.
[IPython Notebook](http://ipython.org/notebook.html) is an interactive web-based environment that allows
mixing rich text, math and graphics with (live) code, similarly to the Mathematica environment.
You can find a static HTML version of the notebooks below in section **[Usage examples](#usage-examples)**.

For a tutorial on using python for scientific computing:

* [Python Scientific Lecture Notes](http://scipy-lectures.github.io/)

Another useful resources for the IPython Notebook:

* [The IPython Notebook](http://ipython.org/ipython-doc/stable/interactive/notebook.html)
* [Notebook examples](http://nbviewer.ipython.org/github/ipython/ipython/blob/master/examples/Notebook/Index.ipynb)
* [A gallery of interesting IPython Notebooks](https://github.com/ipython/ipython/wiki/A-gallery-of-interesting-IPython-Notebooks)

##Installation

Briefly, the installation consist in installing a scientific python distribution,
downloading FRETBursts sources, and setting a folder for the FRETBursts notebooks.

FRETBursts is loaded running a small script (`load_fretbursts.py`) placed
in the notebooks folder. The first time you need to edit `load_fretbursts.py`
to specify where the FRETBursts source directory is on your system.

In the following you can find more detailed installation instructions
for different platforms.

###MS Windows

In order to run the code you need to install a scientific python
distribution like [Anaconda](https://store.continuum.io/cshop/anaconda/).
The free version of Anaconda includes all the needed dependencies.
Any other scientific python distribution (for example
[Enthought Canopy](https://www.enthought.com/products/canopy/))
will work as well.

Once a python distribution is installed, download the latest version
of [FRETBursts](https://github.com/tritemio/FRETBursts) from *GitHub*.
If new to git, we recomend to use the graphical application
[SourceTree](http://www.sourcetreeapp.com/) (selecting the option of
using the embedded git).

The most user friendly way to use FRETBursts is through an IPython Notebook.
The following paragraph shows how to configure it.

####Configuring IPython Notebook

When starting the IPython server, it will show a default folder for the notebooks.
You can create a launcher to start the IPython Notebook server on any local folder.

To create the launcher, right click on the
*IPython Notebook icon* -> *Properties* and paste
the notebook folder in the *Start in* field.

Now, on double click, a browser should pop up showing the list
of notebooks. Chrome browser is suggested.

###Linux and Mac OS X

On Linux or Mac OS X you can also use the [Anaconda](https://store.continuum.io/cshop/anaconda/) distribution.

Alternatively, these are the software dependencies (hint: on Mac OS X you can use MacPorts):

 - Python 2.7
 - Numpy/Scipy (any version from 2013 on)
 - Matplotlib with qt (pyside) backend (1.3.x or greater)
 - IPython 1.x (2.x suggested)
 - PyTables 3.x (optional)
 - a modern browser (Chrome suggested)


##Documentation

The FRETBursts documentation is hosted on ReadTheDocs:

* [FRETBursts Documentation]()

We also provide a list of IPython notebooks showing typical workflows 
for smFRET analysis and illustrates the FRETBursts functionality.

* [usALEX - Workflow](http://nbviewer.ipython.org/urls/raw.github.com/tritemio/FRETBursts/master/notebooks/usALEX%2520-%2520Workflow.ipynb)



##Development

The documentation is built using [Sphinx](http://sphinx-doc.org/) (1.2.2 or later) and
the [napoleon extension](https://pypi.python.org/pypi/sphinxcontrib-napoleon).
A notebook that builds the HTML docs can be found in `notebooks/dev/docs/`.

The unit tests are written with [pytest](http://pytest.org/latest/).
Notebooks that execute the unit tests can be found in `notebooks/dev/test/`.
In the same folder a notebook for regression testing is provided.


##Acknowledgements

This work was supported by NIH grants R01 GM069709 and R01 GM095904.

##License and Copyrights

FRETBursts - A bursts analysis toolkit for single and multi-spot smFRET data.

Copyright (C) 2014  Antonino Ingargiola - *tritemio @ gmail.com*

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 2, as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You can find a full copy of the license in the file LICENSE.txt
