"""Full-sample validation: for every SPARC galaxy, fit a Burkert halo,
invert the coupling K(r), fit the Lorentzian, and collect (K0, rc).

Usage:
    python run_full_sample.py [--data DIR] [--out DIR] [--kappa 1e-6.5] [--rhotarget 1e8.2]

Data: SPARC *_rotmod.dat files (Lelli et al. 2016, AJ 152, 157).
Download from Zenodo record 16284118 and unpack so that DIR contains
the <NAME>_rotmod.dat files.

Output: k_results.json and K_universality_full.png in --out.
Expected runtime: ~15-30 min for the full SPARC sample on a laptop.
"""
import argparse
import glob
import json
import os

import numpy as np
from scipy.optimize import least_squares

from invert_K import R_burkert, invert_K_from_R, fit_lorentzian, lorentzian

# np.trapz was removed in NumPy 2.0 (renamed np.trapezoid)
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

G = 4.3009e-6  # kpc (km/s)^2 / Msun

KAPPA_DEFAULT = 10 ** -6.5
RHO_TARGET_DEFAULT = 10 ** 8.2


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
    return best[1].x  # (lrhoB, lr0)


def process_one(name, data_dir, kappa, rho_target):
    d = load_rotmod(os.path.join(data_dir, f"{name}_rotmod.dat"))
    if len(d) < 5:
        return {"name": name, "ok": False, "error": "too few data points"}
    r, v, e = d[:, 0], d[:, 1], d[:, 2]
    vgas, vdisk, vbul = d[:, 3], d[:, 4], d[:, 5]
    try:
        lrhoB, lr0 = fit_burkert(r, v, e, vgas, vdisk, vbul)
    except Exception as exc:
        return {"name": name, "ok": False, "error": f"burkert fit: {exc}"}
    rg = np.linspace(0.02, float(r.max()), 2000)
    R = R_burkert(rg, lrhoB, lr0)
    K, thp, valid = invert_K_from_R(rg, R, kappa, rho_target)
    if not valid:
        return {"name": name, "ok": False, "error": "(theta')^2 < 0"}
    try:
        K0, rc, R2 = fit_lorentzian(rg, K)
    except Exception as exc:
        return {"name": name, "ok": False, "error": f"lorentzian fit: {exc}"}
    return {"name": name, "ok": True, "n": len(r),
            "vmax": float(v.max()), "K0": K0, "rc": rc, "R2": R2}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data",
                    help="directory with <NAME>_rotmod.dat files")
    ap.add_argument("--out", default=".",
                    help="output directory for k_results.json and figure")
    ap.add_argument("--kappa", type=float, default=KAPPA_DEFAULT)
    ap.add_argument("--rhotarget", type=float, default=RHO_TARGET_DEFAULT)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.data, "*_rotmod.dat")))
    if not files:
        raise SystemExit(
            f"No *_rotmod.dat found in {args.data}. Download SPARC data "
            "(Zenodo 16284118) and unpack the rotmod files there.")
    names = [os.path.basename(f).replace("_rotmod.dat", "") for f in files]
    print(f"{len(names)} galaxies found.")

    results = []
    for i, name in enumerate(names):
        res = process_one(name, args.data, args.kappa, args.rhotarget)
        results.append(res)
        status = (f"K0={res['K0']:.0f} rc={res['rc']:.3f} R2={res['R2']:.3f}"
                  if res["ok"] else f"FAILED ({res['error']})")
        print(f"[{i + 1}/{len(names)}] {name:12s} {status}")

    ok = [r for r in results if r["ok"]]
    K0 = np.array([r["K0"] for r in ok])
    rc = np.array([r["rc"] for r in ok])
    print(f"\n{len(ok)}/{len(results)} succeeded.")
    print(f"K0 = {K0.mean():.1f} ± {K0.std():.1f} "
          f"({K0.std() / K0.mean() * 100:.1f}% scatter)")
    print(f"rc = {rc.mean():.4f} ± {rc.std():.4f} "
          f"({rc.std() / rc.mean() * 100:.1f}% scatter)")

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "k_results.json"), "w") as f:
        json.dump(results, f, indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    vm = np.array([r["vmax"] for r in ok])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].errorbar(vm, K0, yerr=K0.std(), fmt="ko", capsize=4, ms=6)
    axes[0].axhline(K0.mean(), color="r", ls="--",
                    label=f"Mean={K0.mean():.0f}±{K0.std():.0f}")
    axes[0].set_xlabel(r"$V_{\max}$ (km/s)")
    axes[0].set_ylabel(r"$K_0$ (kpc$^{-2}$)")
    axes[0].legend(); axes[0].grid(alpha=0.3)
    axes[0].set_title("Peak coupling strength")
    axes[1].errorbar(vm, rc, yerr=rc.std(), fmt="ko", capsize=4, ms=6)
    axes[1].axhline(rc.mean(), color="r", ls="--",
                    label=f"Mean={rc.mean():.3f}±{rc.std():.3f} kpc")
    axes[1].set_xlabel(r"$V_{\max}$ (km/s)")
    axes[1].set_ylabel(r"$r_c$ (kpc)")
    axes[1].legend(); axes[1].grid(alpha=0.3)
    axes[1].set_title("Coupling range")
    fig.suptitle(f"Universality of the inverted coupling ({len(ok)} SPARC galaxies)",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "K_universality_full.png"), dpi=150)
    print("Wrote k_results.json and K_universality_full.png")


if __name__ == "__main__":
    main()
