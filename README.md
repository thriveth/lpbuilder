lpbuilder
=========

Stand-alone version of the line profile builder of grism.
This is mostly to share a working version with the community that is easier to call and initialize than the full grism version. Development is still going on in grism, but of course suggestions etc. are welcome for this version too.

Usage:
------

For demonstration purposes, simply run the file as a script either calling `python profilebuilder.py` or using the `%run` magic from the IPython prompt, or whatever method one might prefer.

For a bit of more tinkering, import the module and check the docstrings for usage.

Dependencies:
-------------

* Python 2.7 (may work with earlier versions but not tested).
* Traits
* TraitsUI
* Chaco
* Pandas
* SciPy

The former three are all part of the Enthought Tool Suite which can be installed from PyPI as `ets`, the latter two are both part of the standard [SciPy stack][1].


[1]: http://scipy.org
