Release checklist
-----------------

To create a release, simply create a tag that starts with 'v' (e.g. 'v3.0.0')::

    git tag -a "v3.0.0"
    git push origin "v3.0.0"

The package will be automatically tested, released on GitHub and uploaded to PyPI.