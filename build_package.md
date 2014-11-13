# Notes on building FREBursts package

## Sources cleanup

    python setup.py clean --all

## Build step

As a first step always do:

    python setup.py build

## Build distribution packages

**Source distribution**:

    python setup.py sdist

**Wheel distribution (platform specific)**

    python setup.py bdist_wheel

**Windows installer (GUI)**

    python setup.py bdist_wininst


# References

* [Python Packaging User Guide](https://python-packaging-user-guide.readthedocs.org/en/latest/)
* [Building and Distributing Packages with Setuptools](https://pythonhosted.org/setuptools/setuptools.html)


