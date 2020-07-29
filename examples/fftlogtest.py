r"""
FFTLog-Test
===========

This example is a translation of `fftlogtest.f` from the Fortran package
`FFTLog`, which was presented in Appendix B of [Hami00]_ and published at
http://casa.colorado.edu/~ajsh/FFTLog. It serves as an example for the python
package `pyfftlog` (which is a Python version of `FFTLog`), in the same manner
as the original file `fftlogtest.f` serves as an example for Fortran package
`FFTLog`.

**What follows is the original documentation from the file `fftlogtest.f`:**

This is fftlogtest.f
--------------------

This is a simple test program to illustrate how `FFTLog` works. The test
transform is:

.. math::
    :label: hamtest1

    \int^\infty_0 r^{\mu+1} \exp\left(-\frac{r^2}{2} \right)\
    J_\mu(k, r)\ k\ {\rm d}r = k^{\mu+1} \exp\left(-\frac{k^2}{2}
    \right)


Disclaimer
'''''''''''

`FFTLog` does NOT claim to provide the most accurate possible solution of the
continuous transform (which is the stated aim of some other codes). Rather,
`FFTLog` claims to solve the exact discrete transform of a
logarithmically-spaced periodic sequence. If the periodic interval is wide
enough, the resolution high enough, and the function well enough behaved
outside the periodic interval, then `FFTLog` may yield a satisfactory
approximation to the continuous transform.

Observe:

1. How the result improves as the periodic interval is enlarged. With the
   normal FFT, one is not used to ranges orders of magnitude wide, but this is
   how `FFTLog` prefers it.
2. How the result improves as the resolution is increased. Because the function
   is rather smooth, modest resolution actually works quite well here.
3. That the central part of the transform is more reliable than the outer
   parts. Experience suggests that a good general strategy is to double the
   periodic interval over which the input function is defined, and then to
   discard the outer half of the transform.
4. That the best bias exponent seems to be :math:`q = 0`.
5. That for the critical index :math:`\mu = -1`, the result seems to be offset
   by a constant from the 'correct' answer.
6. That the result grows progressively worse as mu decreases below -1.

The analytic integral above fails for :math:`\mu \le -1`, but `FFTLog` still
returns answers.  Namely, `FFTLog` returns the analytic continuation of the
discrete transform.  Because of ambiguity in the path of integration around
poles, this analytic continuation is liable to differ, for :math:`\mu \le -1`,
by a constant from the 'correct' continuation given by the above equation.

`FFTLog` begins to have serious difficulties with aliasing as :math:`\mu`
decreases below :math:`-1`, because then :math:`r^{\mu+1} \exp(-r^2/2)` is far
from resembling a periodic function. You might have thought that it would help
to introduce a bias exponent :math:`q = \mu`, or perhaps :math:`q = \mu+1`, or
more, to make the function :math:`a(r) = A(r) r^{-q}` input to `fhtq` more
nearly periodic. In practice a nonzero :math:`q` makes things worse.

A symmetry argument lends support to the notion that the best exponent here
should be :math:`q = 0,` as empirically appears to be true. The symmetry
argument is that the function :math:`r^{\mu+1} \exp(-r^2/2)` happens to be the
same as its transform :math:`k^{\mu+1} \exp(-k^2/2)`. If the best bias exponent
were q in the forward transform, then the best exponent would be :math:`-q`
that in the backward transform; but the two transforms happen to be the same in
this case, suggesting :math:`q = -q`, hence :math:`q = 0`.

This example illustrates that you cannot always tell just by looking at a
function what the best bias exponent :math:`q` should be. You also have to look
at its transform.  The best exponent :math:`q` is, in a sense, the one that
makes both the function and its transform look most nearly periodic.

"""
import pyfftlog
import numpy as np
import matplotlib.pyplot as plt


###############################################################################
# Define the parameters you wish to use
# -------------------------------------
#
# The presets are the *Reasonable choices of parameters* from `fftlogtest.f`.

# Range of periodic interval
logrmin = -4
logrmax = 4

# Number of points (Max 4096)
n = 64

# Order mu of Bessel function
mu = 0

# Bias exponent: q = 0 is unbiased
q = 0

# Sensible approximate choice of k_c r_c
kr = 1

# Tell fhti to change kr to low-ringing value
# WARNING: kropt = 3 will fail, as interaction is not supported
kropt = 1

# Forward transform (changed from dir to tdir, as dir is a python fct)
tdir = 1


###############################################################################
# Computation related to the logarithmic spacing
# ----------------------------------------------

# Central point log10(r_c) of periodic interval
logrc = (logrmin + logrmax)/2

print(f"Central point of periodic interval at log10(r_c) = {logrc}")

# Central index (1/2 integral if n is even)
nc = (n + 1)/2.0

# Log-spacing of points
dlogr = (logrmax - logrmin)/n
dlnr = dlogr*np.log(10.0)


###############################################################################
# Compute input function: :math:`r^{\mu+1}\exp\left(-\frac{r^2}{2}\right)`
# ------------------------------------------------------------------------

r = 10**(logrc + (np.arange(1, n+1) - nc)*dlogr)
ar = r**(mu + 1)*np.exp(-r**2/2.0)


###############################################################################
# Initialize FFTLog transform - note fhti resets `kr`
# ---------------------------------------------------

kr, xsave = pyfftlog.fhti(n, mu, dlnr, q, kr, kropt)
print(f"pyfftlog.fhti: new kr = {kr}")

###############################################################################
# Call `pyfftlog.fht` (or `pyfftlog.fhtl`)
# ----------------------------------------

logkc = np.log10(kr) - logrc
print(f"Central point in k-space at log10(k_c) = {logkc}")

# rk = r_c/k_c
rk = 10**(logrc - logkc)

# Transform
# ak = pyfftlog.fftl(ar.copy(), xsave, rk, tdir)
ak = pyfftlog.fht(ar.copy(), xsave, tdir)

###############################################################################
# Compute Output function: :math:`k^{\mu+1}\exp\left(-\frac{k^2}{2}\right)`
# -------------------------------------------------------------------------

k = 10**(logkc + (np.arange(1, n+1) - nc)*dlogr)
theo = k**(mu + 1)*np.exp(-k**2/2.0)

###############################################################################
# Plot result
# -----------

plt.figure()

# Input
ax1 = plt.subplot(121)
plt.title(r'$r^{\mu+1}\ \exp(-r^2/2)$')
plt.xlabel('r')

plt.loglog(r, ar, 'k', lw=2)

plt.grid(axis='y', c='0.9')

# Transformed result
ax2 = plt.subplot(122, sharey=ax1)
plt.title(r'$k^{\mu+1} \exp(-k^2/2)$')
plt.xlabel('k')

plt.loglog(k, theo, 'k', lw=2, label='Theoretical')
plt.loglog(k, ak, 'r--', lw=2, label='FFTLog')

plt.legend()
plt.ylim([1e-8, 1e1])

ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.grid(axis='y', c='0.9')

# Switch off spines
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
plt.tight_layout(rect=[0, 0, 1, .9])

# Main title
plt.suptitle(r"$\int_0^\infty r^{\mu+1}\ \exp(-r^2/2)\ J_\mu(k,r)\ " +
             r"k\ {\rm d}r = k^{\mu+1} \exp(-k^2/2)$", y=0.98)

plt.show()

###############################################################################
# Print values
# ------------

print('           k                 a(k)       k^(mu+1) exp(-k^2/2)')
print('----------------------------------------------------------------')
for i in range(n):
    print(f"{k[i]:18.6e} {ak[i]:18.6e} {theo[i]:18.6e}")
