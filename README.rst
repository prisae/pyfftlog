.. image:: https://readthedocs.org/projects/pyfftlog/badge/?version=latest
   :target: https://pyfftlog.readthedocs.io/en/latest
   :alt: Documentation Status
.. image:: https://travis-ci.org/prisae/pyfftlog.svg?branch=master
   :target: https://travis-ci.org/prisae/pyfftlog
   :alt: Travis-CI
.. image:: https://coveralls.io/repos/github/prisae/pyfftlog/badge.svg?branch=master
   :target: https://coveralls.io/github/prisae/pyfftlog?branch=master
   :alt: Coveralls
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3830364.svg
   :target: https://doi.org/10.5281/zenodo.3830364
   :alt: Zenodo DOI


.. sphinx-inclusion-marker


`pyfftlog` - A python version of FFTLog
=======================================

This is a python version of the logarithmic FFT code *FFTLog* as presented in
Appendix B of `Hamilton (2000)
<http://dx.doi.org/10.1046/j.1365-8711.2000.03071.x>`_ and published at
`casa.colorado.edu/~ajsh/FFTLog <http://casa.colorado.edu/~ajsh/FFTLog>`_.

A simple `f2py`-wrapper (`fftlog`) can be found on `github.com/prisae/fftlog
<https://github.com/prisae/fftlog>`_.  Tests have shown that `fftlog` is a bit
faster than `pyfftlog`, but `pyfftlog` is easier to implement, as you only need
`NumPy` and `SciPy`, without the need to compile anything.

I hope that `FFTLog` will make it into `SciPy` in the future, which will make
this project redundant. (If you have the bandwidth and are willing to chip in
have a look at `SciPy PR #7310 <https://github.com/scipy/scipy/pull/7310>`_.)

Be aware that `pyfftlog` has not been tested extensively. It works fine for the
test from the original code, and my use case, which is `pyfftlog.fftl` with
`mu=0.5` (sine-transform), `q=0` (unbiased), `k=1`, `kropt=1`, and `tdir=1`
(forward). Please let me know if you encounter any issues.

- **Documentation**: https://pyfftlog.readthedocs.io
- **Source Code**: https://github.com/prisae/pyfftlog


Description of FFTLog from the FFTLog-Website
---------------------------------------------

FFTLog is a set of fortran subroutines that compute the fast Fourier or Hankel
(= Fourier-Bessel) transform of a periodic sequence of logarithmically spaced
points.

FFTLog can be regarded as a natural analogue to the standard Fast Fourier
Transform (FFT), in the sense that, just as the normal FFT gives the exact (to
machine precision) Fourier transform of a linearly spaced periodic sequence, so
also FFTLog gives the exact Fourier or Hankel transform, of arbitrary order m,
of a logarithmically spaced periodic sequence.

FFTLog shares with the normal FFT the problems of ringing (response to sudden
steps) and aliasing (periodic folding of frequencies), but under appropriate
circumstances FFTLog may approximate the results of a continuous Fourier or
Hankel transform.

The FFTLog algorithm was originally proposed by `Talman (1978)
<http://dx.doi.org/10.1016/0021-9991(78)90107-9>`_.

*For the full documentation, see* `casa.colorado.edu/~ajsh/FFTLog
<http://casa.colorado.edu/~ajsh/FFTLog>`_.


Installation
------------

You can install pyfftlog either via **conda**:

.. code-block:: console

   conda install -c conda-forge pyfftlog

or via **pip**:

.. code-block:: console

   pip install pyfftlog


License, Citation, and Credits
------------------------------

Released to the public domain under the `CC0 1.0 License
<http://creativecommons.org/publicdomain/zero/1.0>`_.

All releases have a Zenodo-DOI, which can be found on `10.5281/zenodo.3830364
<https://doi.org/10.5281/zenodo.3830364>`_.

Be kind and give credits by citing `Hamilton (2000)
<http://dx.doi.org/10.1046/j.1365-8711.2000.03071.x>`_. See the
`references-section
<https://pyfftlog.readthedocs.io/en/stable/references.html>`_ in the manual for
full references.
