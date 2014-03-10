
FRETBursts overview
===================

Description
-----------

**FRETBursts** is an opensource toolkit for analysis of timestamps
series from confocal single-molecule FRET (smFRET) experiments. In the
spirit of reproducible research, this software allows the authors and
others to reproduce previous research and to perform new one. FRETBursts
is open to public scrutiny and the authors are committed to promptly fix
bugs whenever they are discovered.

FRETBursts allows to analyze both
`single-spot <http://dx.doi.org/10.1126/science.283.5408.1676>`__ and
`multi-spot smFRET <http://dx.doi.org/10.1117/12.2003704>`__ data.
Alternating laser excitation
(`ALEX <http://dx.doi.org/10.1529/biophysj.104.054114>`__) scheme is
also supported.

Main analysis features includes:

-  time-dependent background estimation (for example in 30s windows)
-  sliding-window burst search with background-dependent threshold
-  burst corrections: background, leakage (bleed-through), direct
   excitation, gamma-factor
-  per-bursts E and S calculation
-  burst selections based on multiple criteria (for ex.: burst size,
   burst width, E, S, ...)
-  fit routines for FRET efficiency (1 and 2-gaussians histogram fit,
   MLE Poisson models, weighted least squares models, weighted
   expectation maximization, etc...)

Moreover FRETBursts includes a large set of modular plot functions for
background, time-traces, rate-traces, E, S, ALEX histogram and more.

Motivations
-----------

This software aims to be a reference implementation for both established
and novel algorithms related to bursts analysis of smFRET data.

Despite the broad diffusion of smFRET experiments on freely diffusing
molecules, before FRETBursts, no complete smFRET burst analysis software
was freely available on internet. Each group have re-implemented the
analysis in-house with little or no code sharing. This situation is
clearly sub-optimal either because all the advances in the burst
analysis are not readily available to a wide public of smFRET users and
because subtle differences in implementation make the comparison of
experiments performed by different groups problematic.

We envision FRETBursts both as a state-of-the-art burst analysis package
for basic smFRET users, and as a benchmark for advanced users willing to
explore new algorithms or to compare alternative implementations.

Implementation
--------------

**FRETBursts** is written in the Python programming language, and
depends on the standard scientific-python stack (Numpy, Scipy,
Matplotlib, IPython). All the software dependencies can be readily
fulfilled by installing one of the scientific python distributions
available free-of-charge for all the main platforms Windows, Mac OSX and
Linux.

Python is a modern, interpreted programming language with focus on
ease-of use and code readability. Python natively supports the main
programming paradigms: procedural, object-oriented and functional. The
python scientific stack is a comprehensive set of numerical libraries
mostly written in compiled languages (C and Fortran) providing the high
performances of a compiled language from an easy to use environment.
