Steps to Release scikit-ued
===========================

These are the steps to take to create a release of the package ``skued``:

1. Switch to the release branch

 .. code-block:: bash
 
    $ git checkout release
    $ git merge master --no-commit

2. Update the following files with the correct version.

    - ``skued/__init__.py``
    - ``CHANGELOG.rst``
    
3. Add and commit the changes.

 .. code-block:: bash
 
    $ git commit -a -m 'Release {VERSION}'
    $ git tag {VERSION} -m '{VERSION}'
    
4. Build the distributables and publish to release branch.

 .. code-block:: bash
 
    $ git push origin
    $ python setup.py sdist bdist_wheel upload
