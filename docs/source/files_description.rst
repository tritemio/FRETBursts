
Description of the files
========================

Here a brief descriprion of the main FRETBursts files.


``burstlib.py``
---------------

.. automodule:: fretbursts.burstlib


``loader.py``
--------------

.. automodule:: fretbursts.loader


``select_bursts.py``
--------------------

See :mod:`fretbursts.select_bursts`.


``burst_plot.py``
-----------------

.. automodule:: fretbursts.burst_plot


``background.py``
-----------------

.. automodule:: fretbursts.background


``burstsearch`` (folder)
---------------------------

.. automodule:: fretbursts.burstsearch


``dataload`` (folder)
---------------------

.. automodule:: fretbursts.dataload


``fit`` (folder)
----------------

.. automodule:: fretbursts.fit

See :ref:`fit-section`.


.. _load_fretbursts:

``load_fretbursts.py``
----------------------

`load_fretbursts.py <https://github.com/tritemio/FRETBursts/blob/master/notebooks/load_fretbursts.py>`_
is a small script used to import **FRETBursts**
in an ipython notebook when you want to load fretbursts from the sources
folder (not from system installation).

The script is typically ran at the beginning of a notebook as::

    %run load_fretbursts --nogui --source

The script performs some basic operations:

-  Find the **FRETBursts** source folder (from a configuration file written by
   the installation notebook) and set two variables (``NOTEBOOK_DIR`` and
   ``FRETBURSTS_DIR``) to easily switch between the notebooks folder and
   the FRETBursts source folder.

-  Import **FRETBursts**, matplotlib (as ``plt``)  and some ipython functions
   (``display``, ``Math``, etc...)

-  If `git <http://git-scm.com/>`__ is found, it displays the current
   **FRETBursts** revision (and eventual files modified since last
   revision)

-  Enable inline plots (``%matplotlib inline``) and optionally enable the Qt
   GUI integration (to be able to launch open-file dialogs for example).

The script options are:

* ``--nogui`` do not load the Qt GUI
* ``--nompl`` do not load matplotlib
* ``--nostyle`` do not modify matplotlib default plot syle
* ``--source`` load FRETBursts from the source folder even if a system
  installation is found



