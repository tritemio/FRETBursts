
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


``fretburst.py``
----------------

.. automodule:: fretbursts


``burstlib.py``
---------------

.. automodule:: burstlib


``loaders.py``
--------------

.. automodule:: loaders


``burst_plot.py``
-----------------

.. automodule:: burst_plot


``background.py``
-----------------

.. automodule:: background


``select_bursts.py``
--------------------

.. automodule:: select_bursts


``burstsearch`` (folder)
---------------------------

.. automodule:: burstsearch


``dataload`` (folder)
---------------------

.. automodule:: dataload


``fit`` (folder)
----------------

.. automodule:: fit
