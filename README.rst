`pyfftlog` - python version of FFTLog
=====================================

This is a python version of the logarithmic FFT code *FFTLog* as presented in
Appendix B of [Hamilton_2000]_ and published at `casa.colorado.edu/~ajsh/FFTLog
<http://casa.colorado.edu/~ajsh/FFTLog>`_.

A simple `f2py`-wrapper (`fftlog`) can be found on `github.com/prisae/fftlog
<https://github.com/prisae/fftlog>`_.  Tests have shown that `fftlog` is a bit
faster than `pyfftlog`, but `pyfftlog` is easier to implement, as you only need
`NumPy` and `SciPy`, without the need to compile anything.

I hope that `FFTLog` will make it into `SciPy` in the future, which will make
this project redundant.


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

The FFTLog algorithm was originally proposed by [Talman_1978]_.

*For the full documentation, see*
`casa.colorado.edu/~ajsh/FFTLog <http://casa.colorado.edu/~ajsh/FFTLog>`_.


References
----------

.. [Hamilton_2000] Hamilton, A. J. S., 2000, Uncorrelated modes of the
    non-linear power spectrum: Monthly Notices of the Royal Astronomical
    Society, 312, pages 257-284; DOI: `10.1046/j.1365-8711.2000.03071.x
    <http://dx.doi.org/10.1046/j.1365-8711.2000.03071.x>`_; Website of FFTLog:
    `casa.colorado.edu/~ajsh/FFTLog <http://casa.colorado.edu/~ajsh/FFTLog>`_.

.. [Talman_1978] Talman, J. D., 1978, Numerical Fourier and Bessel transforms
    in logarithmic variables: Journal of Computational Physics, 29, pages
    35-48; DOI: `10.1016/0021-9991(78)90107-9
    <http://dx.doi.org/10.1016/0021-9991(78)90107-9>`_.


License and Credits
-------------------

Released to the public domain under the `CC0 1.0 License
<http://creativecommons.org/publicdomain/zero/1.0>`_. Be kind and give credits
to [Hamilton_2000]_.
