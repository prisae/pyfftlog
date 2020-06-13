r"""
Geophysical Electromagnetic modelling
=====================================

In this example we use `pyfftlog` to obtain time-domain EM data from
frequency-domain data and vice versa. We do this with using analytical
halfspace solution in both domains, and comparing the transformed responses to
the true result. The analytical halfspace solutions are computed using
`empymod` (see https://empymod.github.io).
"""
import empymod
import pyfftlog
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as iuSpline


###############################################################################
# Model and Survey parameters
# ---------------------------

# Acquisition parameters
freqs = np.logspace(-4, 4, 301)  # Frequencies (Hz)
times = np.logspace(-3, 3, 301)  # Times (s)
signal = 0                       # Impulse response

# Source and receiver
src = [0, 0, 100]     # At the origin, 100m below surface
rec = [6000, 0, 200]  # At an inline offset of 6 km, 200 m below surface

# Resistivity
depth = [0]      # Interface at z = 0, default for empymod.analytical
res = [2e14, 1]  # Horizontal resistivity [air, subsurface]
aniso = [1, 2]   # Anisotropy [air, subsurface]

# Collect parameters
analytical = {
    'src': src,
    'rec': rec,
    'res': res[1],
    'aniso': aniso[1],
    'solution': 'dhs',  # Diffusive half-space solution
    'verb': 2,
}

dipole = {
    'src': src,
    'rec': rec,
    'depth': depth,
    'res': res,
    'aniso': aniso,
    'ht': 'dlf',
    'verb': 2,
}


###############################################################################
# Analytical solutions :math:`f` and :math:`t` domains
# ----------------------------------------------------

# Frequency Domain
f_ana = empymod.analytical(**analytical, freqtime=freqs)

# Time Domain
t_ana = empymod.analytical(**analytical, freqtime=times, signal=signal)


###############################################################################
# FFTLog :: :math:`f \Rightarrow t`
# ---------------------------------
#
# `FFTLog` is directly implemented in `empymod`, but only for the
# frequency-to-time transformation.

ftl = empymod.dipole(**dipole, freqtime=times, signal=signal, ft='fftlog')


###############################################################################
# FFTLog :: :math:`t \Rightarrow f`
# ---------------------------------

# FFTLog parameters
pts_per_dec = 5   # Increase if not precise enough
add_dec = [0, 0]  # e.g. [-2, 2] to add 2 decades on each side
q = 0             # -1 - +1; can improve results

# Work with angular frequencies
wfreqs = 2*np.pi*freqs

# Calculate minimum and maximum required time
rmin = np.log10(1/wfreqs.max()) + add_dec[0]
rmax = np.log10(1/wfreqs.min()) + add_dec[1]
n = np.int(rmax - rmin)*pts_per_dec

# Pre-allocate output
f_fftlog = np.zeros(freqs.shape, dtype=complex)

# Loop over Sine, Cosine transform.
for mu in [0.5, -0.5]:

    # Central point log10(r_c) of periodic interval
    logrc = (rmin + rmax)/2

    # Central index (1/2 integral if n is even)
    nc = (n + 1)/2.

    # Log spacing of points
    dlogr = (rmax - rmin)/n
    dlnr = dlogr*np.log(10.)

    # Calculate required input x-values (time)
    time = 10**(logrc + (np.arange(1, n+1) - nc)*dlogr)

    # Initialize FFTLog
    kr, xsave = pyfftlog.fhti(n, mu, dlnr, q, kr=1, kropt=1)

    # Calculate fcalc with adjusted kr
    logkc = np.log10(kr) - logrc
    fcalc = 10**(logkc + (np.arange(1, n+1) - nc)*dlogr)

    # rk = r_c/k_r; adjust for Fourier transform scaling
    rk = 10**(logrc - logkc)*np.pi/2

    # Calculate required times with the analytical solution
    t_fftlog = empymod.analytical(**analytical, freqtime=time, signal=signal)

    # Carry out FFTLog
    f_resp = pyfftlog.fftl(t_fftlog.copy(), xsave, rk, 1)

    # Interpolate for required frequencies
    f_respspline = iuSpline(np.log(fcalc), f_resp)

    if mu > 0:
        f_fftlog += -1j*f_respspline(np.log(wfreqs))
    else:
        f_fftlog += f_respspline(np.log(wfreqs))


###############################################################################
# Compare them
# ------------

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 4))

# TIME DOMAIN
ax0.set_title(r'Frequency domain')
ax0.set_xlabel('Frequency (Hz)')
ax0.set_ylabel('Amplitude (V/m)')
ax0.semilogx(freqs, f_ana.real, 'k-', label='Analytical')
ax0.semilogx(freqs, f_ana.imag, 'k-')
ax0.semilogx(freqs, f_fftlog.real, 'C1--', label=r'FFTLog, $\mu=-0.5$')
ax0.semilogx(freqs, f_fftlog.imag, 'C2--', label=r'FFTLog, $\mu=0.5$')
ax0.legend(loc='best')
ax0.grid(which='both', c='.95')

# TIME DOMAIN
ax1.set_title(r'Time domain')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Amplitude (V/m)')
ax1.semilogx(times, t_ana, 'k', label='Analytical')
ax1.loglog(times, ftl, 'C3--', label=r'FFTLog, $\mu=0.5$')
ax1.legend(loc='best')
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
ax1.grid(which='both', c='.95')

fig.tight_layout()
fig.show()

###############################################################################
empymod.Report(pyfftlog)
