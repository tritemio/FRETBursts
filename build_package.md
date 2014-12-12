# Developers note on package building

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
    
## Upload to PyPI

Once `.pypirc` is configured, run:
    
    twine upload dist/fretbursts-X.X.zip -ppassword
    
## Zenodo

To disable/enable the DOI creation for a release, login on Zenodo and enable/disable the repository in the GitHub section.


# References

* [Python Packaging User Guide](https://python-packaging-user-guide.readthedocs.org/en/latest/)
* [Building and Distributing Packages with Setuptools](https://pythonhosted.org/setuptools/setuptools.html)
* [PIP User Guide](http://pip.readthedocs.org/en/latest/user_guide.html)
* [setup.py vs requirements.txt](https://caremad.io/blog/setup-vs-requirement/)

