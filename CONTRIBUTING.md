# How to contribute to scikit-ued

<!-- This contributing guide is inspired from scikit-image's [https://github.com/scikit-image/scikit-image]-->

Contents:

* [Development Process](#development-process)
* [Build environment setup](#build-environment-setup)
* [Guidelines](#guidelines)
* [Testing](#testing)
* [Writing documentation](#writing-documentation)
* [Bugs](#bugs)

## Development process

Here's the long and short of it:

1.  If you are a first-time contributor:
    * Go to <https://github.com/LaurentRDC/scikit-ued> and click the "fork" button to create your own copy of the project.

    * Clone the project to your local computer:

            git clone https://github.com/your-username/scikit-ued.git

    * Change the directory:

            cd scikit-ued

    * Add the upstream repository:

            git remote add upstream https://github.com/LaurentRDC/scikit-ued.git

    * Now, you have remote repositories named:

        * `upstream`, which refers to the `scikit-ued` repository
        * `origin`, which refers to your personal fork

2.  Develop your contribution:
    * Pull the latest changes from upstream:

            git checkout master
            git pull upstream master

    * Create a branch for the feature you want to work on. Since the branch name will appear in the merge message, use a sensible name such as 'my-new-feature':

            git checkout -b my-new-feature

    * Commit locally as you progress (`git add` and `git commit`)

3.  To submit your contribution:

    * Push your changes back to your fork on GitHub:

            git push origin my-new-feature

    * Enter your GitHub username and password (repeat contributors or advanced users can remove this step by [connecting to GitHub with SSH](https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh)).

    * Go to GitHub. The new branch will show up with a green Pull
        Request button - click it.

    * If you want, post on the [mailing
        list](https://mail.python.org/mailman3/lists/scikit-ued.python.org/)
        to explain your changes or to ask for review.


4.  Review process:

    * Reviewers (the other developers and interested community members) will write inline and/or general comments on your Pull Request (PR) to help you improve its implementation, documentation, and style.
    * To update your pull request, make your changes on your local repository and commit. As soon as those changes are pushed up (to the same branch as before) the pull request will update automatically.
    * Continuous integration is triggered after each Pull Request update to build the code, run unit tests, and other checks of your branch. The tests must pass before your PR can be merged. 
    * A pull request must be approved by a an admininistrator before merging.

5.  Document changes

    - If your change introduces any API modifications or fixes a known issue, please update `CHANGELOG.rst`.


### Divergence between `upstream master` and your feature branch

If GitHub indicates that the branch of your Pull Request can no longer be merged automatically, merge the master branch into yours:

    git fetch upstream master
    git merge upstream/master

If any conflicts occur, they need to be fixed before continuing. See which files are in conflict using:

    git status

Which displays a message like:

    Unmerged paths:
      (use "git add <file>..." to mark resolution)

      both modified:   file_with_conflict.txt

Inside the conflicted file, you'll find sections like these:

    <<<<<<< HEAD
    The way the text looks in your branch
    =======
    The way the text looks in the master branch
    >>>>>>> master

Choose one version of the text that should be kept, and delete the rest:

    The way the text looks in your branch

Now, add the fixed file:

    git add file_with_conflict.txt

Once you've fixed all merge conflicts, do:

    git commit

### Build environment setup

To create an appropriate development environment, you need to install the base requirements (`requirements.txt`) as well as extra, development requirements (`dev-requirements.txt`)

    pip install -r requirements
    pip install -r dev-requirements

## Guidelines

* All code should have tests.
* All code should be documented, to the same [standard](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard) as NumPy and SciPy.
* For new functionality, always add an example to the tutorials 
* No changes are ever committed without review and approval.

### Stylistic Guidelines

* Run the [`black`](https://black.readthedocs.io/en/stable/) code formatter.

* Use the following import conventions:

    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import ndimage as ndi

    # only in Cython code
    cimport numpy as cnp
    cnp.import_array()
    ```

* Refer to array dimensions as (plane), row, column, not as x, y, z. See `Coordinate conventions <numpy-images-coordinate-conventions>` in the user guide for more information.

* Use relative module imports, i.e. `from ..simulation import xyz` rather than `from skued.simulation import xyz`.

## Testing

`scikit-ued` has an extensive test suite that ensures that the code works as expected on many possible systems. `scikit-ued` uses the [`pytest`](https://docs.pytest.org/en/latest/) framework.

Run all tests using:

    pytest skued

To run the tests of a particular submodule (e.g. `skued/image`):

    pytest skued/image

## Writing documentation

Documentation is very important for any scientific software project. Documentation is located in the `docs` directory.

### Building docs

To build docs:

``` sh
python -m sphinx docs build
```

Then, all the HTML files will be generated in `build/sphinx/html/`.

### Testing documentation

Code snippets in the documentation are checked via `doctest`. To test the documentation snippets:

``` sh
python -m sphinx -b doctest docs build
```

## Bugs

Please [report bugs on GitHub](https://github.com/LaurentRDC/scikit-ued/issues).
