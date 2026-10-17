"""Validity-boundary mechanism variables for SPARC K(r) inversion.

For each galaxy where the fixed-parameter inversion violates (theta')^2 >= 0,
compute 17 mechanism variables from the field equations:
  A. failure location (r_f, zone, negative fraction, depth)
  B. term decomposition at r_f (A vs B dominance, stabilizer sign)
  C. exact repair distances via LINEAR feasibility in kappa / rho_target
  D. galaxy/data features (Burkert params, Vmax, n_data, density ratio)

Self-contained: needs only numpy + scipy. Run from this folder:
    python mechanism.py --data ../data --out results

Output: results/mechanisms.json
"""
import argparse
import glob
import json
import os

import numpy as np
from scipy.optimize import least_squares

G = 4.3009e-6  # kpc (km/s)^2 / Msun
KAPPA0 = 10 ** -6.5
RHO0 = 10 ** 8.2
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


# ---- halo model ----
def R_burkert(r, lrhoB, lr0):
    rho0B, r0 = 10 ** lrhoB, 10 ** lr0
    y = np.asarray(r, float) / r0
    rho = rho0B / ((1 + y) * (1 + y ** 2))
    return np.sqrt(np.maximum(rho, 1e-30))


def v2_burkert(r_kpc, rho0B, r0):
    y = np.asarray(r_kpc, float) / r0
    s = np.linspace(0, 1, 400)
    out = np.zeros_like(y)
    for i, yy in enumerate(y):
        ss = s * yy
        integ = ss ** 2 / ((1 + ss) * (1 + ss ** 2))
        out[i] = _trapz(integ, ss)
    M = 4 * np.pi * rho0B * r0 ** 3 * out
    return G * M / np.maximum(r_kpc, 1e-9)


def load_rotmod(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append([float(x) for x in p[:6]])
    a = np.array(rows)
    m = (a[:, 1] > 0) & (a[:, 2] > 0) & np.isfinite(a).all(axis=1)
    return a[m]


def fit_burkert(r, v, e, vgas, vdisk, vbul):
    vbar2 = vgas * np.abs(vgas) + 0.5 * vdisk ** 2 + 0.7 * vbul ** 2

    def resid(p):
        lrhoB, lr0 = p
        vh2 = v2_burkert(r, 10 ** lrhoB, 10 ** lr0)
        vm = np.sqrt(np.maximum(vh2 + vbar2, 1e-9))
        return (v - vm) / e

    best = None
    for guess in ([7.0, 0.3], [8.0, 0.7], [6.5, 0.0]):
        sol = least_squares(resid, guess, bounds=([4.0, -1.0], [10.0, 1.5]),
                            max_nfev=4000)
        chi2 = float(np.sum(sol.fun ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, sol)
    return best[1].x


# ---- inversion helpers ----
def phi2_decomp(rg, R, kappa, rho_target):
    """(theta')^2 profile and A, C decomposition: N = A + kappa*C."""
    dr = rg[1] - rg[0]
    Rp = np.gradient(R, dr)
    Rpp = np.gradient(Rp, dr)
    A = Rpp + 2 * Rp / rg
    C = (rho_target - R ** 2) * R
    return (A + kappa * C) / R, A, C


def check_valid(rg, R, kappa=KAPPA0, rho_target=RHO0):
    t2, _, _ = phi2_decomp(rg, R, kappa, rho_target)
    return bool(np.all(t2 >= -1e-12))


def repair_distances(rg, R):
    """Exact repair distances via linear feasibility.

    kappa: need A(r)+k*C(r) >= 0 for all r -> interval [lo, hi].
    rho_target: need rho_t >= R^2 - A/(k*R) for all r -> always feasible.
    """
    _, A, C = phi2_decomp(rg, R, KAPPA0, RHO0)
    lo, hi = -np.inf, np.inf
    for a, c in zip(A, C):
        if c > 0:
            lo = max(lo, -a / c)
        elif c < 0:
            hi = min(hi, -a / c)
        elif a < 0:
            lo, hi = np.inf, -np.inf
            break
    kap_feasible = bool(lo <= hi) and bool(np.isfinite(lo))
    if kap_feasible and lo > 0:
        d_logkap = float(np.log10(max(float(lo), 1e-30) / KAPPA0))
    elif kap_feasible:
        d_logkap = 0.0
    else:
        d_logkap = float("inf")
    with np.errstate(divide="ignore", invalid="ignore"):
        need = R ** 2 - A / (KAPPA0 * R)
    rho_min = float(np.nanmax(need))
    d_logrho = float(np.log10(max(rho_min, 1e-30) / RHO0))
    return {"d_logkap": d_logkap, "kap_feasible": kap_feasible,
            "d_logrho": d_logrho}


def mechanism_variables(name, data_dir):
    d = load_rotmod(os.path.join(data_dir, f"{name}_rotmod.dat"))
    r, v, e = d[:, 0], d[:, 1], d[:, 2]
    lrhoB, lr0 = fit_burkert(r, v, e, d[:, 3], d[:, 4], d[:, 5])
    rho0B = 10 ** lrhoB
    rg = np.linspace(0.02, float(r.max()), 2000)
    R = R_burkert(rg, lrhoB, lr0)
    t2, A, C = phi2_decomp(rg, R, KAPPA0, RHO0)

    neg = t2 < 0
    idx = int(np.argmax(neg))
    rf = float(rg[idx])
    rmax = float(rg[-1])
    zone = "inner" if rf < 0.25 * rmax else ("mid" if rf < 0.6 * rmax else "outer")
    a_rf, c_rf = float(A[idx]), float(C[idx])
    B_rf = KAPPA0 * c_rf
    tot = a_rf + B_rf
    rep = repair_distances(rg, R)
    return {
        "name": name,
        "lrhoB": float(lrhoB), "lr0": float(lr0),
        "vmax": float(v.max()), "n_data": int(len(r)), "r_max": float(r.max()),
        "rf": rf, "rf_frac": float(rf / rmax), "zone": zone,
        "neg_frac": float(neg.mean()),
        "t2min": float(t2.min()),
        "log_abs_t2min": float(np.log10(max(abs(float(t2.min())), 1e-300))),
        "A_rf": a_rf, "B_rf": B_rf,
        "frac_A": float(a_rf / tot) if abs(tot) > 1e-300 else float("nan"),
        "B_neg": bool(B_rf < 0),
        "log_rho0_over_target": float(np.log10(rho0B / RHO0)),
        "d_logkap": rep["d_logkap"], "kap_feasible": rep["kap_feasible"],
        "d_logrho": rep["d_logrho"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="../data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    names = sorted(os.path.basename(f).replace("_rotmod.dat", "")
                   for f in glob.glob(os.path.join(args.data, "*_rotmod.dat")))
    print(f"{len(names)} galaxies.")
    boundary, datapoor = [], []
    for i, name in enumerate(names):
        d = load_rotmod(os.path.join(args.data, f"{name}_rotmod.dat"))
        if len(d) < 5:
            datapoor.append(name)
            continue
        r, v, e = d[:, 0], d[:, 1], d[:, 2]
        try:
            lrhoB, lr0 = fit_burkert(r, v, e, d[:, 3], d[:, 4], d[:, 5])
        except Exception as exc:
            print(f"[{i+1}] {name}: burkert fit failed ({exc})")
            continue
        rg = np.linspace(0.02, float(r.max()), 2000)
        if not check_valid(rg, R_burkert(rg, lrhoB, lr0)):
            boundary.append(name)
            print(f"[{i+1}/{len(names)}] {name}: beyond validity domain")

    print(f"\n{len(boundary)} beyond-domain, {len(datapoor)} data-poor: {datapoor}")
    mechs = []
    for j, name in enumerate(boundary):
        m = mechanism_variables(name, args.data)
        mechs.append(m)
        print(f"[{j+1}/{len(boundary)}] {name}: rf={m['rf']:.2f} "
              f"neg_frac={m['neg_frac']:.2f} B_neg={m['B_neg']} "
              f"d_logkap={m['d_logkap']:.2f} d_logrho={m['d_logrho']:.2f}")

    with open(os.path.join(args.out, "mechanisms.json"), "w") as f:
        json.dump({"datapoor": datapoor, "mechanisms": mechs}, f, indent=1)
    print(f"Wrote results/mechanisms.json ({len(mechs)} galaxies)")


if __name__ == "__main__":
    main()
