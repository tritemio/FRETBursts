Pre-release
-----------
- check that CI tests pass
- check that docs are built correctly

Release
-------
- git tag <version>
- python build sdist
- twine upload dist/FRETBursts-<version>.tar.gz

- cd notebooks
- python ../fretbursts/tests/nbrun.py
- copy notebooks to FRETBursts_notebooks, commit and push to github

- update fretbursts-feedstock recipe: version, sha256 and build number
- commit & push to personal fork, open PR for conda-forge/fretbursts-feedstock
- when build succeeds merge on conda-forge/fretbursts-feedstock

Update front-page
----------------
- Update README.md
- Regenerate github landing page using github settings
- cherry-pick commit to change fonts in gh-pages branch, push to GitHub.
