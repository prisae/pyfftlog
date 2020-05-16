import numpy as np
from numpy.testing import assert_allclose

from pyfftlog import fhti, fht


def test_fftlog():
    # This is the example provided in `fftlogtest.f` of the original code. It
    # is the same test that is shown as an example in the gallery.
    # See the example for more details.

    # Parameters.
    logrmin = -4
    logrmax = 4
    n = 64
    mu = 0
    q = 0
    kr = 1
    kropt = 1
    tdir = 1
    logrc = (logrmin + logrmax)/2
    nc = (n + 1)/2.0
    dlogr = (logrmax - logrmin)/n
    dlnr = dlogr*np.log(10.0)

    # Calculate input function: $r^{\mu+1}\exp\left(-\frac{r^2}{2}\right)$.
    r = 10**(logrc + (np.arange(1, n+1) - nc)*dlogr)
    ar = r**(mu + 1)*np.exp(-r**2/2.0)

    # Initialize FFTLog transform - note fhti resets `kr`
    kr, xsave = fhti(n, mu, dlnr, q, kr, kropt)
    assert_allclose(kr, 0.953538967579)
    logkc = np.log10(kr) - logrc
    # rk = 10**(logrc - logkc)

    # Transform
    # ak = fftl(ar.copy(), xsave, rk, tdir)
    ak = fht(ar.copy(), xsave, tdir)

    # Calculate Output function: $k^{\mu+1}\exp\left(-\frac{k^2}{2}\right)$
    k = 10**(logkc + (np.arange(1, n+1) - nc)*dlogr)
    theo = k**(mu + 1)*np.exp(-k**2/2.0)

    # Check values
    assert_allclose(theo[theo > 1e-3], ak[theo > 1e-3], rtol=1e-8, atol=5e-5)
