"""Analytic inversion: from a halo density profile, invert the virtual-space
coupling K(r) via the static field equations.

Generalized GL static (D=1, u=0, real scalar K): Psi = R(r) e^{i theta(r)}
  Real: nabla^2 R - R (theta')^2 + kappa (rho_target - R^2) R = 0
  Imag: 2 R' theta' + R theta'' + (2/r) R theta' + K(r) R = 0

Given R(r) = sqrt(rho_DM(r)):
  (theta')^2 = [nabla^2 R + kappa (rho_target - R^2) R] / R   (must be >= 0)
  K(r) = -[2 R' theta' + R theta'' + (2/r) R theta'] / R

Units: r in kpc, rho in Msun/kpc^3, kappa in kpc/Msun, K in 1/kpc^2.
"""
import numpy as np


def R_burkert(r, lrhoB, lr0):
    """Amplitude R(r) = sqrt(rho_Burkert(r))."""
    rho0B, r0 = 10 ** lrhoB, 10 ** lr0
    y = np.asarray(r, float) / r0
    rho = rho0B / ((1 + y) * (1 + y ** 2))
    return np.sqrt(np.maximum(rho, 1e-30))


def invert_K_from_R(r, R, kappa, rho_target):
    """Invert K(r) given radial grid r and amplitude R(r).

    Returns (K, theta_prime, valid) where valid flags (theta')^2 >= 0.
    """
    r = np.asarray(r, float)
    R = np.asarray(R, float)
    dr = r[1] - r[0]
    Rp = np.gradient(R, dr)
    Rpp = np.gradient(Rp, dr)
    lap_R = Rpp + 2 * Rp / r
    num = lap_R + kappa * (rho_target - R ** 2) * R
    t2 = num / R
    valid = bool(np.all(t2 >= -1e-12))
    thp = np.sqrt(np.maximum(t2, 0.0))  # positive branch
    thpp = np.gradient(thp, dr)
    K = -(2 * Rp * thp + R * thpp + (2 / r) * R * thp) / R
    return K, thp, valid


def lorentzian(r, K0, rc):
    return -K0 / (1 + (r / rc) ** 2)


def fit_lorentzian(r, K, rmin=0.05, rmax=8.0):
    """Fit K(r) = -K0 / (1 + (r/rc)^2) on [rmin, rmax]. Returns (K0, rc, R2)."""
    from scipy.optimize import curve_fit
    m = (r > rmin) & (r < rmax)
    popt, _ = curve_fit(lorentzian, r[m], K[m], p0=[200.0, 0.2], maxfev=20000)
    pred = lorentzian(r[m], *popt)
    ss_res = np.sum((K[m] - pred) ** 2)
    ss_tot = np.sum((K[m] - K[m].mean()) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(popt[0]), float(popt[1]), float(R2)
