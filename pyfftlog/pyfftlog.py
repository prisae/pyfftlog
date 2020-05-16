r"""

`pyfftlog` -- Python version of FFTLog
======================================

This is a Python version of the FFTLog Fortran code by Andrew Hamilton,
[Hami00]_.

The function :obj:`scipy.special.loggamma` replaces the file `cdgamma.f` in
the original code, and the functions :func:`scipy.fftpack.rfft` and
:func:`scipy.fftpack.irfft` replace the files `drffti.f`, `drfftf.f`, and
`drfftb.f` in the original code.

The original documentation has just been adjusted where necessary, and put into
a more pythonic format (e.g. using `Parameters` and `Returns` in the
documentation').

**What follows is the original documentation from the file `fftlog.f`:**

THE FFTLog CODE
---------------

FFTLog computes the discrete Fast Fourier Transform or Fast Hankel Transform
(of arbitrary real index) of a periodic logarithmic sequence.

- Version of 13 Mar 2000.
- For more information about FFTLog, see http://casa.colorado.edu/~ajsh/FFTLog.
- Andrew J S Hamilton March 1999.
- Refs: [Talm78]_.

FFTLog computes a discrete version of the Hankel Transform (= Fourier-Bessel
Transform) with a power law bias :math:`(k r)^q`

.. math::
    :label: ham1

    \tilde{a}(k) = \int^\infty_0 a(r) (k r)^q J_{\mu} (k r) k \,dr \, ,

.. math::
    :label: ham2

    a(r) = \int^\infty_0 \tilde{a}(k) (k r)^{-q} J_{\mu} (k r) r \,dk \, ,

where :math:`J_{\mu}` is the Bessel function of order :math:`\mu`. The index
:math:`\mu` may be any real number, positive or negative.

The input array :math:`a_j` is a periodic sequence of length :math:`n`,
uniformly logarithmically spaced with spacing :math:`dlnr`

.. math::
    :label: ham3

    a_j = a(r_j) \quad \text{at} \quad r_j = r_c \exp[(j-j_c) dlnr]

centred about the point :math:`r_c`. The central index :math:`j_c = (n+1)/2` is
1/2 integral if :math:`n` is even. Similarly, the output array
:math:`\tilde{a}_j` is a periodic sequence of length :math:`n`, also uniformly
logarithmically spaced with spacing :math:`dlnr`

.. math::
    :label: ham4

    \tilde{a}_j = \tilde{a}(k_j) \quad \text{at} \quad
    k_j = k_c \exp[(j-j_c) dlnr]

centred about the point :math:`k_c`.

The centre points :math:`r_c` and :math:`k_c` of the periodic intervals may be
chosen arbitrarily; but it would be normal to choose the product

.. math::
    :label: ham5

    kr = k_c r_c = k_j r_{(n+1-j)} = k_{(n+1-j)} r_j

to be about 1 (or 2, or pi, to taste).

The FFTLog algorithm is (see [Hami00]_):

1. FFT the input array :math:`a_j` to obtain the Fourier coefficients
   :math:`c_m` ;
2. Multiply :math:`c_m` by
   :math:`u_m = (kr)^{- i 2 m \pi/(n dlnr)} U_{\mu}[q + i 2 m \pi/(n dlnr)]`
   where :math:`U_{\mu}(x) = 2^x \Gamma[(\mu+1+x)/2] / \Gamma[(\mu+1-x)/2]` to
   obtain :math:`c_m u_m`;
3. FFT :math:`c_m u_m` back to obtain the discrete Hankel transform
   :math:`\tilde{a}_j`.

The Fourier sine and cosine transforms
``````````````````````````````````````

.. math::
    :label: ham6

    \tilde{A}(k) = \sqrt{2/\pi} \int^\infty_0 A(r) \sin(k r) \,dr \, ,

.. math::
    :label: ham7

    \tilde{A}(k) = \sqrt{2/\pi} \int^\infty_0 A(r) \cos(k r) \,dr \, ,

may be regarded as special cases of the Hankel transform with :math:`\mu = 1/2`
and :math:`-1/2` since

.. math::
    :label: ham8

    \sqrt{2/\pi} \sin(x) = \sqrt(x) J_{1/2} (x) \, ,

.. math::
    :label: ham9

    \sqrt{2/\pi} \cos(x) = \sqrt(x) J_{-1/2} (x) \, .


The Fourier transforms may be done by making the substitutions

.. math::
    :label: ham10

    A(r) = a(r) r^{q-1/2} \quad \text{and} \quad
    \tilde{A}(k) = \tilde{a}(k) k^{-q-1/2}

and Hankel transforming :math:`a(r)` with a power law bias :math:`(k r)^q`

.. math::
    :label: ham11

    \tilde{a}(k) = \int^\infty_0 a(r) (k r)^q J_{\pm 1/2} (k r) k \,dr \, .

Different choices of power law bias :math:`q` lead to different discrete
Fourier transforms of :math:`A(r)`, because the assumption of periodicity of
:math:`a(r) = A(r) r^{-q+(1/2)}` is different for different :math:`q`.

If :math:`A(r)` is a power law, :math:`A(r)` proportional to
:math:`r^{q-(1/2)}`, then applying a bias :math:`q` yields a discrete Fourier
transform :math:`\tilde{A}(k)` that is exactly equal to the continuous Fourier
transform, because then :math:`a(r)` is a constant, which is a periodic
function.

The Hankel transform
````````````````````

.. math::
    :label: ham12

    \tilde{A}(k) = \int^\infty_0 A(r) J_{\mu} (k r) k \,dr

may be done by making the substitutions

.. math::
    :label: ham13

    A(r) = a(r) r^q \quad \text{and} \quad \tilde{A}(k) = \tilde{a}(k) k^{-q}

and Hankel transforming :math:`a(r)` with a power law bias :math:`(k r)^q`

.. math::
    :label: ham14

    \tilde{a}(k) = \int^\infty_0 a(r) (k r)^q J_{\mu} (k r) k \,dr \, .

Different choices of power law bias :math:`q` lead to different discrete Hankel
transforms of :math:`A(r)`, because the assumption of periodicity of
:math:`a(r) = A(r) r^{-q}` is different for different :math:`q`.

If :math:`A(r)` is a power law, :math:`A(r)` proportional to :math:`r^q`, then
applying a bias :math:`q` yields a discrete Hankel transform
:math:`\tilde{A}(k)` that is exactly equal to the continuous Hankel transform,
because then :math:`a(r)` is a constant, which is a periodic function.

There are five routines:
````````````````````````
Comments in the subroutines contain further details.

1. **subroutine `fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)`**
   is an initialization routine.

2. **subroutine `fftl(n,a,rk,dir,wsave)`**
   computes the discrete Fourier sine or cosine transform of a logarithmically
   spaced periodic sequence. This is a driver routine that calls `fhtq`.

3. **subroutine `fht(n,a,dir,wsave)`**
   computes the discrete Hankel transform of a logarithmically spaced periodic
   sequence. This is a driver routine that calls `fhtq`.

4. **subroutine `fhtq(n,a,dir,wsave)`**
   computes the biased discrete Hankel transform of a logarithmically spaced
   periodic sequence. **This is the basic FFTLog routine.**

5. **real*8 function `krgood(mu,q,dlnr,kr)`**
   takes an input `kr` and returns the nearest low-ringing `kr`. This is an
   optional routine called by `fhti`.

**END of the original documentation from the file `fftlog.f`**

"""
import numpy as np
from scipy.special import loggamma
from scipy.fftpack import rfft, irfft

__all__ = ['fhti', 'fftl', 'fht', 'fhtq', 'krgood']


def fhti(n, mu, dlnr, q=0, kr=1, kropt=0):
    r"""Initialize the working array xsave used by fftl, fht, and fhtq.

    fhti initializes the working array xsave used by fftl, fht, and fhtq. fhti
    need be called once, whereafter fftl, fht, or fhtq may be called many
    times, as long as n, mu, q, dlnr, and kr remain unchanged. fhti should be
    called each time n, mu, q, dlnr, or kr is changed. The work array xsave
    should not be changed between calls to fftl, fht, or fhtq.

    Parameters
    ----------
    n : int
        Number of points in the array to be transformed; n may be any positive
        integer, but the FFT routines run fastest if n is a product of small
        primes 2, 3, 5.

    mu : float
        Index of J_mu in Hankel transform; mu may be any real number, positive
        or negative.

    dlnr : float
        Separation between natural log of points; dlnr may be positive or
        negative.

    q : float, optional
        Exponent of power law bias; q may be any real number, positive or
        negative. If in doubt, use q = 0, for which case the Hankel transform
        is orthogonal, i.e. self-inverse, provided also that, for n even, kr is
        low-ringing. Non-zero q may yield better approximations to the
        continuous Hankel transform for some functions.
        Defaults to 0 (unbiased).

    kr : float, optional
        k_c r_c where c is central point of array
        = k_j r_(n+1-j) = k_(n+1-j) r_j .
        Normally one would choose kr to be about 1 (default) (or 2, or pi, to
        taste).

    kropt : int, optional; {0, 1, 2, 3}
        - 0 to use input kr as is (default);
        - 1 to change kr to nearest low-ringing kr, quietly;
        - 2 to change kr to nearest low-ringing kr, verbosely;
        - 3 for option to change kr interactively.

    Returns
    -------
    kr : float, optional
        kr, adjusted depending on kropt.

    xsave : array
        Working array used by fftl, fht, and fhtq. Dimension:
        - for q = 0 (unbiased transform): n+3
        - for q != 0 (biased transform): 1.5*n+4
        If odd, last element is not needed.

    """

    # adjust kr
    if kropt == 0:    # keep kr as is
        pass
    elif kropt == 1:  # change kr to low-ringing kr quietly
        kr = krgood(mu, q, dlnr, kr)
    elif kropt == 2:  # change kr to low-ringing kr verbosely
        d = krgood(mu, q, dlnr, kr)
        if abs(kr/d - 1) >= 1e-15:
            kr = d
            print(" kr changed to ", kr)
    else:             # option to change kr to low-ringing kr interactively
        d = krgood(mu, q, dlnr, kr)
        if abs(kr/d-1.0) >= 1e-15:
            print(" change kr = ", kr)
            print(" to low-ringing kr = ", d)
            go = input("? [CR, y=yes, n=no, x=exit]: ")
            if go.lower() in ['', 'y']:
                kr = d
                print(" kr changed to ", kr)
            elif go.lower() == 'n':
                print(" kr left unchanged at ", kr)
            else:
                print("exit")
                return False

    # return if n is <= 0
    if n <= 0:
        return kr

    # The normal FFT is not initialized here as in the original FFTLog code, as
    # the `scipy.fftpack`-FFT routines `rfft` and `irfft` do that internally.
    # Therefore xsave in `pyfftlog` is 2*n+15 elements shorter, and named
    # xsave to not confuse it with xsave from the FFT.

    if q == 0:  # unbiased case (q = 0)
        ln2kr = np.log(2.0/kr)
        xp = (mu + 1)/2.0
        d = np.pi/(n*dlnr)

        m = np.arange(1, (n+1)/2)
        y = m*d  # y = m*pi/(n*dlnr)
        zp = loggamma(xp + 1j*y)
        arg = 2.0*(ln2kr*y + zp.imag)  # Argument of kr^(-2 i y) U_mu(2 i y)

        # Arange xsave: [q, dlnr, kr, cos, sin]
        xsave = np.empty(2*arg.size+3)
        xsave[0] = q
        xsave[1] = dlnr
        xsave[2] = kr
        xsave[3::2] = np.cos(arg)
        xsave[4::2] = np.sin(arg)

        # Altogether 3 + 2*(n/2) elements used for q = 0, which is n+3 for even
        # n, n+2 for odd n.

    else:       # biased case (q != 0)
        ln2 = np.log(2.0)
        ln2kr = np.log(2.0/kr)
        xp = (mu + 1 + q)/2.0
        xm = (mu + 1 - q)/2.0

        # first element of rest of xsave
        y = 0

        # case where xp or xm is a negative integer
        xpnegi = np.round(xp) == xp and xp <= 0
        xmnegi = np.round(xm) == xm and xm <= 0
        if xpnegi or xmnegi:

            # case where xp and xm are both negative integers
            # U_mu(q) = 2^q Gamma[xp]/Gamma[xm] is finite in this case
            if xpnegi and xmnegi:
                # Amplitude and Argument of U_mu(q)
                amp = np.exp(ln2*q)
                if xp > xm:
                    m = np.arange(1,  np.round(xp - xm)+1)
                    amp *= xm + m - 1
                elif xp < xm:
                    m = np.arange(1,  np.round(xm - xp)+1)
                    amp /= xp + m - 1
                arg = np.round(xp + xm)*np.pi

            else:  # one of xp or xm is a negative integer
                # Transformation is singular if xp is -ve integer, and inverse
                # transformation is singular if xm is -ve integer, but
                # transformation may be well-defined if sum_j a_j = 0, as may
                # well occur in physical cases.  Policy is to drop the
                # potentially infinite constant in the transform.

                if xpnegi:
                    print('fhti: (mu+1+q)/2 =', np.round(xp), 'is -ve integer',
                          ', yields singular transform:\ntransform will omit',
                          'additive constant that is generically infinite,',
                          '\nbut that may be finite or zero if the sum of the',
                          'elements of the input array a_j is zero.')
                else:
                    print('fhti: (mu+1-q)/2 =', np.round(xm), 'is -ve integer',
                          ', yields singular inverse transform:\n inverse',
                          'transform will omit additive constant that is',
                          'generically infinite,\nbut that may be finite or',
                          'zero if the sum of the elements of the input array',
                          'a_j is zero.')
                amp = 0
                arg = 0

        else:  # neither xp nor xm is a negative integer
            zp = loggamma(xp + 1j*y)
            zm = loggamma(xm + 1j*y)

            # Amplitude and Argument of U_mu(q)
            amp = np.exp(ln2*q + zp.real - zm.real)
            # note +Im(zm) to get conjugate value below real axis
            arg = zp.imag + zm.imag

        # first element: cos(arg) = ±1, sin(arg) = 0
        xsave1 = amp*np.cos(arg)

        # remaining elements of xsave
        d = np.pi/(n*dlnr)
        m = np.arange(1, (n+1)/2)
        y = m*d  # y = m pi/(n dlnr)
        zp = loggamma(xp + 1j*y)
        zm = loggamma(xm + 1j*y)
        # Amplitude and Argument of kr^(-2 i y) U_mu(q + 2 i y)
        amp = np.exp(ln2*q + zp.real - zm.real)
        arg = 2*ln2kr*y + zp.imag + zm.imag

        # Arrange xsave: [q, dlnr, kr, xsave1, cos, sin]
        xsave = np.empty(3*arg.size+4)
        xsave[0] = q
        xsave[1] = dlnr
        xsave[2] = kr
        xsave[3] = xsave1
        xsave[4::3] = amp
        xsave[5::3] = np.cos(arg)
        xsave[6::3] = np.sin(arg)

        # Altogether 3 + 3*(n/2)+1 elements used for q != 0, which is (3*n)/2+4
        # for even n, (3*n)/2+3 for odd n.  For even n, the very last element
        # of xsave [i.e. xsave(3*m+1)=sin(arg) for m=n/2] is not used within
        # FFTLog; if a low-ringing kr is used, this element should be zero.
        # The last element is computed in case somebody wants it.

    return kr, xsave


def fftl(a, xsave, rk=1, tdir=1):
    r"""Logarithmic fast Fourier transform FFTLog.

    This is a driver routine that calls :func:`fhtq`.

    `fftl` computes a discrete version of the Fourier sine (if mu = 1/2) or
    cosine (if mu = -1/2) transform

    .. math::

        \tilde{A}(k) = \sqrt{2/\pi} \int^\infty_0 A(r) \sin(k r) \,dr \, ,

    .. math::

        \tilde{A}(k) = \sqrt{2/\pi} \int^\infty_0 A(r) \cos(k r) \,dr \, ,

    by making the substitutions

    .. math::

        A(r) = a(r) r^{q-1/2} \quad \text{and} \quad
        \tilde{A}(k) = \tilde{a}(k) k^{-q-1/2}

    and applying a biased Hankel transform to :math:`a(r)`.

    The steps are:
    1. :math:`a(r) = A(r) r^[-dir (q-0.5)]`
    2. call `fhtq` to transform :math:`a(r) \rightarrow \tilde{a}(k)`
    3. :math:`\tilde{A}(k) = \tilde{a}(k) k^[-dir (q+0.5)]`

    `fhti` must be called before the first call to `fftl`, with `mu=1/2` for a
    sine transform, or `mu=-1/2` for a cosine transform.

    A call to `fftl` with `dir=1` followed by a call to `fftl` with `dir=-1`
    (and rk unchanged), or vice versa, leaves the array a unchanged.

    Parameters
    ----------
    a : array
        Array A(r) to transform: a(j) is A(r_j) at r_j = r_c exp[(j-jc) dlnr],
        where jc = (n+1)/2 = central index of array.

    xsave : array
        Working array set up by fhti.

    rk : float, optional
        r_c/k_c = r_j/k_j (a constant, the same constant for any j); rk is not
        (necessarily) the same quantity as kr.  rk is used only to multiply the
        output array by sqrt(rk)^dir, so if you want to do the normalization
        later, or you don't care about the normalization, you can set rk = 1.
        Defaults to 1.

    tdir : int, optional; {1, -1}
        -  1 for forward transform (default),
        - -1 for backward transform.

        A backward transform (dir = -1) is the same as a forward transform with
        q -> -q and rk -> 1/rk, for any kr if n is odd, for low-ringing kr if n
        is even.

    Returns
    -------
    a : array
        Transformed array Ã(k): a(j) is Ã(k_j) at k_j = k_c exp[(j-jc) dlnr].

    """
    fct = a.copy()
    q = xsave[0]
    dlnr = xsave[1]
    kr = xsave[2]

    # centre point of array
    jc = np.array((fct.size + 1)/2.0)
    j = np.arange(fct.size)+1

    # a(r) = A(r) (r/rc)^[-dir*(q-.5)]
    fct *= np.exp(-tdir*(q - 0.5)*(j - jc)*dlnr)

    # transform a(r) -> ã(k)
    fct = fhtq(fct, xsave, tdir)

    # Ã(k) = ã(k) k^[-dir*(q+.5)] rc^[-dir*(q-.5)]
    #      = ã(k) (k/kc)^[-dir*(q+.5)] (kc rc)^(-dir*q) (rc/kc)^(dir*.5)
    lnkr = np.log(kr)
    lnrk = np.log(rk)
    fct *= np.exp(-tdir*((q + 0.5)*(j - jc)*dlnr + q*lnkr - lnrk/2.0))

    return fct


def fht(a, xsave, tdir=1):
    r"""Fast Hankel transform FHT.

    This is a driver routine that calls :func:`fhtq`.

    `fht` computes a discrete version of the Hankel transform

    .. math::

        \tilde{A}(k) = \int^\infty_0 A(r) J_{\mu} (k r) k \,dr \,

    by making the substitutions

    .. math::

        A(r) = a(r) r^q \quad \text{and} \quad
        \tilde{A}(k) = \tilde{a}(k) k^{-q}

    and applying a biased Hankel transform to :math:`a(r)`.

    The steps are:
    1. :math:`a(r) = A(r) r^{-dir q}`
    2. call `fhtq` to transform :math:`a(r) \rightarrow \tilde{a}(k)`
    3. :math:`\tilde{A}(k) = \tilde{a}(k) k^{-dir q}`

    `fhti` must be called before the first call to `fht`.

    A call to `fht` with `dir=1` followed by a call to `fht` with `dir=-1`, or
    vice versa, leaves the array a unchanged.


    Parameters
    ----------
    a : array
        Array A(r) to transform: a(j) is A(r_j) at r_j = r_c exp[(j-jc) dlnr],
        where jc = (n+1)/2 = central index of array.

    xsave : array
        Working array set up by fhti.

    tdir : int, optional; {1, -1}
        -  1 for forward transform (default),
        - -1 for backward transform.

        A backward transform (dir = -1) is the same as a forward transform with
        q -> -q, for any kr if n is odd, for low-ringing kr if n is even.

    Returns
    -------
    a : array
        Transformed array Ã(k): a(j) is Ã(k_j) at k_j = k_c exp[(j-jc) dlnr].

    """
    fct = a.copy()
    q = xsave[0]
    dlnr = xsave[1]
    kr = xsave[2]

    # a(r) = A(r) (r/rc)^(-dir*q)
    if q != 0:
        #  centre point of array
        jc = np.array((fct.size + 1)/2.0)
        j = np.arange(fct.size)+1
        fct *= np.exp(-tdir*q*(j - jc)*dlnr)

    # transform a(r) -> ã(k)
    fct = fhtq(fct, xsave, tdir)

    # Ã(k) = ã(k) (k rc)^(-dir*q)
    #      = ã(k) (k/kc)^(-dir*q) (kc rc)^(-dir*q)
    if q != 0:
        lnkr = np.log(kr)
        fct *= np.exp(-tdir*q*((j - jc)*dlnr + lnkr))

    return fct


def fhtq(a, xsave, tdir=1):
    r"""Kernel routine of FFTLog.

    This is the basic FFTLog routine.

    `fhtq` computes a discrete version of the biased Hankel transform

    .. math::

        \tilde{a}(k) = \int^\infty_0 a(r) (k r)^q J_{\mu} (k r) k \,dr \, .

    `fhti` must be called before the first call to `fhtq`.

    A call to `fhtq` with `dir=1` followed by a call to `fhtq` with `dir=-1`,
    or vice versa, leaves the array a unchanged.

    Parameters
    ----------
    a : array
        Periodic array a(r) to transform: a(j) is a(r_j) at r_j = r_c
        exp[(j-jc) dlnr] where jc = (n+1)/2 = central index of array.

    xsave : array
        Working array set up by fhti.

    tdir : int, optional; {1, -1}
        -  1 for forward transform (default),
        - -1 for backward transform.

        A backward transform (dir = -1) is the same as a forward transform with
        q -> -q, for any kr if n is odd, for low-ringing kr if n is even.

    Returns
    -------
    a : array
        Transformed periodic array ã(k): a(j) is ã(k_j) at k_j = k_c exp[(j-jc)
        dlnr].

    """
    fct = a.copy()
    q = xsave[0]
    n = fct.size

    # normal FFT
    fct = rfft(fct)

    m = np.arange(1, n//2, dtype=int)  # index variable
    if q == 0:  # unbiased (q = 0) transform
        # multiply by (kr)^[- i 2 m pi/(n dlnr)] U_mu[i 2 m pi/(n dlnr)]
        ar = fct[2*m-1]
        ai = fct[2*m]
        fct[2*m-1] = ar*xsave[2*m+1] - ai*xsave[2*m+2]
        fct[2*m] = ar*xsave[2*m+2] + ai*xsave[2*m+1]
        # problem(2*m)atical last element, for even n
        if np.mod(n, 2) == 0:
            ar = xsave[-2]
            if (tdir == 1):  # forward transform: multiply by real part
                # Why? See http://casa.colorado.edu/~ajsh/FFTLog/index.html#ure
                fct[-1] *= ar
            elif (tdir == -1):  # backward transform: divide by real part
                # Real part ar can be zero for maximally bad choice of kr.
                # This is unlikely to happen by chance, but if it does, policy
                # is to let it happen.  For low-ringing kr, imaginary part ai
                # is zero by construction, and real part ar is guaranteed
                # nonzero.
                fct[-1] /= ar

    else:  # biased (q != 0) transform
        # multiply by (kr)^[- i 2 m pi/(n dlnr)] U_mu[q + i 2 m pi/(n dlnr)]
        # phase
        ar = fct[2*m-1]
        ai = fct[2*m]
        fct[2*m-1] = ar*xsave[3*m+2] - ai*xsave[3*m+3]
        fct[2*m] = ar*xsave[3*m+3] + ai*xsave[3*m+2]

        if tdir == 1:  # forward transform: multiply by amplitude
            fct[0] *= xsave[3]
            fct[2*m-1] *= xsave[3*m+1]
            fct[2*m] *= xsave[3*m+1]

        elif tdir == -1:  # backward transform: divide by amplitude
            # amplitude of m=0 element
            ar = xsave[3]
            if ar == 0:
                # Amplitude of m=0 element can be zero for some mu, q
                # combinations (singular inverse); policy is to drop
                # potentially infinite constant.
                fct[0] = 0
            else:
                fct[0] /= ar

            # remaining amplitudes should never be zero
            fct[2*m-1] /= xsave[3*m+1]
            fct[2*m] /= xsave[3*m+1]

        # problematical last element, for even n
        if np.mod(n, 2) == 0:
            m = int(n/2)
            ar = xsave[3*m+2]*xsave[3*m+1]
            if tdir == 1:  # forward transform: multiply by real part
                fct[-1] *= ar
            elif (tdir == -1):  # backward transform: divide by real part
                # Real part ar can be zero for maximally bad choice of kr.
                # This is unlikely to happen by chance, but if it does, policy
                # is to let it happen.  For low-ringing kr, imaginary part ai
                # is zero by construction, and real part ar is guaranteed
                # nonzero.
                fct[-1] /= ar

    # normal FFT back
    fct = irfft(fct)

    # reverse the array and at the same time undo the FFTs' multiplication by n
    # => Just reverse the array, the rest is already done in drfft.
    fct = fct[::-1]

    return fct


def krgood(mu, q, dlnr, kr):
    r"""Return optimal kr.

    Use of this routine is optional.

    Choosing kr so that

    .. math::

        (k r)^{- i pi/dlnr} U_{\mu}(q + i pi/dlnr)

    is real may reduce ringing of the discrete Hankel transform, because it
    makes the transition of this function across the period boundary smoother.

    Parameters
    ----------
    mu : float
        index of J_mu in Hankel transform; mu may be any real number, positive
        or negative.

    q : float
        exponent of power law bias; q may be any real number, positive or
        negative. If in doubt, use q = 0, for which case the Hankel transform
        is orthogonal, i.e. self-inverse, provided also that, for n even, kr is
        low-ringing. Non-zero q may yield better approximations to the
        continuous Hankel transform for some functions.

    dlnr : float
        separation between natural log of points; dlnr may be positive or
        negative.

    kr : float, optional
        k_c r_c where c is central point of array
        = k_j r_(n+1-j) = k_(n+1-j) r_j .
        Normally one would choose kr to be about 1 (default) (or 2, or pi, to
        taste).

    Returns
    -------
    krgood : float
        low-ringing value of kr nearest to input kr. ln(krgood) is always
        within dlnr/2 of ln(kr).

    """
    if dlnr == 0:
        return kr

    xp = (mu + 1.0 + q)/2.0
    xm = (mu + 1.0 - q)/2.0
    y = 1j*np.pi/(2.0*dlnr)
    zp = loggamma(xp + y)
    zm = loggamma(xm + y)

    # low-ringing condition is that following should be integral
    arg = np.log(2.0/kr)/dlnr + (zp.imag + zm.imag)/np.pi

    # return low-ringing kr
    return kr*np.exp((arg - np.round(arg))*dlnr)
