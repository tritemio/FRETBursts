
Description of the files
========================

Here a brief list of the main FRETBursts files.

``load_fretbursts.py``
----------------------

This is a small script in the notebook folder that is run in an IPython
Notebook to import **FRETBursts** (``run load_fretbursts``).

The script performs 3 basic operations:

-  Find the FRETBursts source folder (either from an environment
   variable or from a default value set in the script) and set two
   variables (``NOTEBOOK_DIR`` and ``FRETBURSTS_DIR``) to easily switch
   between the notebooks folder and the FRETBursts source folder.

-  Load **FRETBursts**.

-  If `git <http://git-scm.com/>`__ is found, it displays the current
   **FRETBursts** revision (and eventual files modified since last
   revision)

To quickly switch between the notebooks dir and the FRETBursts source
dir, use::

    %cd $NOTEBOOK_DIR

or::

    %cd $FRETBURSTS_DIR`



``burstlib.py``
---------------

.. automodule:: fretbursts.burstlib


``loader.py``
--------------

.. automodule:: fretbursts.loader


``burst_plot.py``
-----------------

.. automodule:: fretbursts.burst_plot


``background.py``
-----------------

.. automodule:: fretbursts.background


``select_bursts.py``
--------------------

.. automodule:: fretbursts.select_bursts


``burstsearch`` (folder)
---------------------------

.. automodule:: fretbursts.burstsearch


``dataload`` (folder)
---------------------

.. automodule:: fretbursts.dataload


``fit`` (folder)
----------------

.. automodule:: fretbursts.fit
