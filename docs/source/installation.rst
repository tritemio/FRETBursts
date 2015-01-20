.. _installation:

FRETBursts Installation
=======================

FRETBursts can be installed as a standard python package (either via PIP
or conda).

In any case, to easily install the different dependencies, you need to
install a **scientific python distribution**.

After installing FRETBursts, it is recommended that you download the
`FRETBursts notebooks <https://github.com/tritemio/FRETBursts_notebooks/archive/master.zip>`__
to get started. See :ref:`running_fretbursts`.


.. _package_install:

Package installation (PIP)
--------------------------

This is the preferred method if you are just starting with FRETBursts.

Start by installing `Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__
python distribution (free version).

Then you can install latest FRETBursts release with through PIP::

    pip install fretbursts --upgrade

Use the previous command every time you want to upgrade to the latest stable
release.

Package installation (conda)
----------------------------

If you prefer installing FRETBursts through `conda` (the Anaconda package
manager), you need to add the FRETBursts channel::

    conda config --add channels tritemio

and then install FRETBursts with::

    conda install fretbursts

Install from GitHub (development)
---------------------------------

You can install the latest development version from GitHub with::

    pip install git+git://github.com/tritemio/FRETBursts.git

Alternatively you can clone FRETBursts git repository and run::

    python setup.py build
    python setup.py install

The optimized C extensions are installed in both cases. Make sure that
the dependencies `lmfit` and `seaborn` have been installed. If not
just install them with::

    pip install lmfit seaborn

.. _source_install:

Source installation (development)
---------------------------------

To perform a "source installation", i.e. executing FRETBursts from the source
folder, download the code from GitHub, and execute the ipython notebook
`FRETBursts Installation`.

In the case of source installation, instead of the normal import, we need to
use a little helper script (:ref:`load_fretbursts.py <load_fretbursts>`) that
finds the sources and import FRETBursts.
Once the script `load_fretbursts.py` is copied in your notebook
folder you can load/import FRETBursts with::

    %run load_fretbursts --nogui -- sources

A copy of `load_fretbursts.py` can be found under `notebooks` in the
source tree.
