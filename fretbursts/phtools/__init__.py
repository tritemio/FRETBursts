"""
This folder contains the core functions to manipulate timestamps,
including burst search and photon rates computations.
Additionally, data structures for storing and manipulating bursts data
are provided.

Burst search and photon counting functions (to count number of donor and acceptor
photons in each burts) are provided both as a pure python implementation and as
an optimized Cython (compiled) version. The cython version is usually 10 or 20
times faster. `burstlib.py` will load the Cython functions, falling back to the
pure python version if the compiled version is not found.
"""
