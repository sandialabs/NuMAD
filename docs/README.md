# Visit the [NuMAD website](http://NuMAD.github.io/NuMAD) for more information.

## NuMAD Documentation

### Compile Instructions

These instructions work for both Linux and Windows. For Windows, remember to
replace slashes (`/`) in paths with backslashes (`\ `).

#### Setup Sphinx (One Time Only)

1. Install [Anaconda Python](https://www.anaconda.com/distribution/).

2. Create the Sphinx environment:
   
   ```
   > conda create -c conda-forge -n _wssphinx git click colorama colorclass future pip sphinxcontrib-bibtex "sphinx_rtd_theme<1"
   > conda activate _wssphinx
   (_wssphinx) > pip install sphinxcontrib-matlabdomain sphinxext-remoteliteralinclude sphinx-multiversion
   (_wssphinx) > conda deactivate
   >
   ```

#### Testing the Current Branch

The documentation for the current branch can be built locally for inspection 
prior to publishing. They are built in the `docs/_build` directory. Note, 
unlike the final documentation, version tags and other branches will not be 
available. 

To test the current branch, use the following:

```
> conda activate _wssphinx
(_wssphinx) > cd path/to/NuMAD
(_wssphinx) > sphinx-build -b html docs docs/_build/html
(_wssphinx) > conda deactivate
>
```

NOTE: for a clean build ``(_wssphinx) > sphinx-build -b html -a docs docs/_build/html``

The front page of the docs can be accessed at 
`docs/_build/html/index.html`. 

#### Publishing Final Version Remotely

The NuMAD docs are rebuilt automatically following every merge commit made 
to the master or dev branch of the [sandialabs/NuMAD](
https://github.com/sandialabs/NuMAD) repository.


## Best Practices
  - Run spell check (not built into most text editors)

### Formatting Guidelines
  - use ``insert code`` to reference code
  - Title `####` with overline
  - Heading 1 `======`
  - Heading 2 `------`
  - Heading 3 `^^^^^^`
  - Heading 4 `""""""`
  - Heading 5 `++++++`
  - Use this style guide: https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html

### Terminology Guidelines
  - post-processing (not postprocessing)
  - pre-processing (not preprocessing)  
  - nonlinear (not non-linear)
  - MATLAB (not Matlab)
