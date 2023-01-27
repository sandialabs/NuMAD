.. _developer:

Developing pynumad
=====================

Installation
------------

To maintain a local installation, developers should use the following commands::
    
    git clone https://cee-gitlab.sandia.gov/ecamare/pynumad.git
    cd pynumad
    pip install -e .

Testing
-------
To run tests locally, run::

    pytest pynumad

at the root of the repository. Note that this requires the installation
of pytest.


Documentation
------------------

Building docs
^^^^^^^^^^^^^^^

To build docs locally, navigate to ``pynumad/docs`` and run::

    make html

After building, the static html files can be found in ``_build/html``.

Conventions
^^^^^^^^^^^^

Docstrings
^^^^^^^^^^^

The documentation for pynumad adheres to NumPy style docstrings. Not only does this
help to keep a consistent style, but it is also necessary for the API documentation
to be parsed and displayed correctly. For an example of what this should look like::

    def func(arg1, arg2):
    """Summary line.

    Extended description of function.

    Parameters
    ----------
    arg1 : int
        Description of arg1
    arg2 : str
        Description of arg2

    Returns
    -------
    bool
        Description of return value

    """
    return True

Additional examples can be found in the 
`napoleon documentation <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_.
The following boilerplate can be copy-pasted into the top of a function definition
to help get things started::

    """Summary line.

    Extended description of function.

    Parameters
    ----------

    Returns
    -------


    """