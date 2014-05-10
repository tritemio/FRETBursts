FRETBursts
==========

Project description
-------------------

**[FRETBursts](https://github.com/tritemio/FRETBursts)** is an opensource
toolkit for analysis of timestamps series from confocal single-molecule FRET
(smFRET) experiments. In the spirit of
reproducible research, this software allows the authors and others to
reproduce previous research and to perform new one. FRETBursts is open to
public scrutiny and the authors are committed to promptly fix bugs
whenever they are discovered. We have a suite of unit tests

FRETBursts allows to analyze both [single-spot](http://dx.doi.org/10.1126/science.283.5408.1676)
and [multi-spot smFRET](http://dx.doi.org/10.1117/12.2003704) data.
Alternating laser excitation ([ALEX](http://dx.doi.org/10.1529/biophysj.104.054114))
scheme is supported as well.

Main analysis features includes:

- time-dependent background estimation (for example in 30s windows)
- sliding-window burst search with background-dependent threshold
- burst corrections: background, leakage (bleed-through), direct excitation,
gamma-factor
- per-bursts E and S calculation
- post-burst-search selection based on multiple criteria (for ex.:
burst size, burst width, E, S, ...). Defining a new burst selection
criterium requires only a couple of lines of code.
- fit routines for FRET efficiency (1 and 2-gaussians histogram fit,
MLE Poisson models, weighted least squares models,
weighted expectation maximization, etc...)

Moreover FRETBursts includes a large set of modular plot functions for
background, time-traces, rate-traces, E, S, ALEX histograms, weighted kernel
density estimation ([KDE](http://en.wikipedia.org/wiki/Kernel_density_estimation))
and more. Thanks to the excellent [Matplotlib](http://matplotlib.org/) library,
FRETBursts can produce publication-quality plots out of the box.

Motivations
-----------

This software aims to be a reference implementation for both established
and novel algorithms related to bursts analysis of smFRET data.

Despite the broad diffusion of smFRET experiments on freely diffusing
molecules, before FRETBursts, no complete smFRET burst analysis software is
freely available on internet. Each group have re-implemented the analysis
in-house with little or no code sharing. This situation is clearly sub-optimal
either because all the advances in the burst analysis are not readily
available to a wide public of smFRET users and because subtle differences in
implementation make the comparison of experiments performed by different
groups problematic.

We envision FRETBursts both as a state-of-the-art burst analysis package
for basic smFRET users, and as a benchmark for advanced users willing
to explore new algorithms or to compare alternative implementations.

Software Environment
--------------------
FRETBursts is written in the [python programming language](http://www.python.org/)
using the standard scientific-python stack of libraries (Numpy, Scipy, Matplotlib, IPython).

FRETBursts uses consolidated software engineering techniques (version control,
unit testing, regression testing) and a workflow based on IPython Notebooks
to ensure robustness and reproducibility of the results.

Tutorials, provided as IPython notebooks, can be edited and executed live.
[IPython Notebook](http://ipython.org/notebook.html) is an interactive
web-based environment that allows mixing rich text, math and graphics with
(live) code, similarly to the Mathematica environment.
You can find a static HTML version of the notebooks [below](#documentation).

###External resources

For a tutorial on using python for scientific computing:

* [Python Scientific Lecture Notes](http://scipy-lectures.github.io/)

Another useful resource for the IPython Notebook:

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

###Windows

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


##How to use

We provide a list of IPython notebooks showing typical workflows
for smFRET analysis and illustrating FRETBursts functionalities.
These notebooks can be executed locally using publically available datasets
(see below).

* [FRETBursts - usALEX Workflow](http://nbviewer.ipython.org/urls/raw.github.com/tritemio/FRETBursts_notebooks/master/notebooks/FRETBursts%2520-%2520usALEX%2520Workflow.ipynb)
* [FRETBursts - 8-spot smFRET analysis](http://nbviewer.ipython.org/urls/raw.github.com/tritemio/FRETBursts_notebooks/master/notebooks/FRETBursts%2520-%25208-spot%2520smFRET%2520analysis.ipynb)

The FRETBursts documentation is hosted on ReadTheDocs:

* [FRETBursts Documentation](http://fretbursts.readthedocs.org/)

We provide a public dataset [1] to testing and demonstration of the FRETBursts
functionalities. These data files are needed to run the tutorials.

[1] A. Ingargiola, S. Chung (2014): smFRET example datasets for the FRETBursts
software. [DOI 10.6084/m9.figshare.1019906](http://dx.doi.org/10.6084/m9.figshare.1019906)

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

Copyright (C) 2014 Antonino Ingargiola - *tritemio @ gmail.com*

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
version 2, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You can find a full copy of the license in the file LICENSE.txt

