"""Sample validation: invert K(r) for 8 characteristic SPARC galaxies spanning
Vmax = 48-216 km/s, fit the Lorentzian, and plot the 2x4 panel.

Usage:
    python sample_validation.py [--data data] [--out figs]

Expected runtime: ~1-2 min on a laptop.
"""
import argparse
import os

import numpy as np

from invert_K import R_burkert, invert_K_from_R, fit_lorentzian, lorentzian
from full_validation import load_rotmod, fit_burkert

SAMPLE = ["DDO154", "D631-7", "NGC3741", "NGC6503",
          "NGC3198", "F583-1", "NGC0055", "NGC2903"]

KAPPA = 10 ** -6.5
RHO_TARGET = 10 ** 8.2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="figs")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.out, exist_ok=True)
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True, sharey=True)
    axes = axes.flatten()
    K0s, rcs = [], []
    r_fine = np.linspace(0.02, 8, 400)

    for ax, name in zip(axes, SAMPLE):
        d = load_rotmod(os.path.join(args.data, f"{name}_rotmod.dat"))
        r, v, e = d[:, 0], d[:, 1], d[:, 2]
        lrhoB, lr0 = fit_burkert(r, v, e, d[:, 3], d[:, 4], d[:, 5])
        rg = np.linspace(0.02, float(r.max()), 2000)
        R = R_burkert(rg, lrhoB, lr0)
        K, thp, valid = invert_K_from_R(rg, R, KAPPA, RHO_TARGET)
        assert valid, f"(theta')^2 < 0 for {name}"
        K0, rc, R2 = fit_lorentzian(rg, K)
        K0s.append(K0)
        rcs.append(rc)
        m = (rg > 0.05) & (rg < 8)
        ax.plot(rg[m], K[m], "k-", lw=1.2, alpha=0.7, label="Inverted")
        ax.plot(r_fine, lorentzian(r_fine, K0, rc), "r--", lw=1.5,
                label=f"$K_0$={K0:.0f}, $r_c$={rc:.3f}")
        ax.set_title(f"{name} ($V_{{\\max}}$={v.max():.0f} km/s)", fontsize=10)
        ax.set_xlim(0, 8)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ax.set_xlabel("r (kpc)")
        ax.set_ylabel(r"$K(r)$ (kpc$^{-2}$)")
        print(f"{name:10s} K0={K0:7.1f} rc={rc:.3f} R2={R2:.3f}")

    fig.suptitle(r"Inverted coupling $K(r)$ with Lorentzian fits: "
                 r"$K(r)=-K_0/[1+(r/r_c)^2]$", fontsize=13)
    fig.tight_layout()
    fp = os.path.join(args.out, "K_inversion_fits.png")
    fig.savefig(fp, dpi=150)
    print(f"saved {fp}")

    K0s = np.array(K0s)
    rcs = np.array(rcs)
    print(f"\nK0 = {K0s.mean():.1f} ± {K0s.std():.1f} "
          f"({K0s.std() / K0s.mean() * 100:.1f}% scatter)")
    print(f"rc = {rcs.mean():.4f} ± {rcs.std():.4f} "
          f"({rcs.std() / rcs.mean() * 100:.1f}% scatter)")


if __name__ == "__main__":
    main()
