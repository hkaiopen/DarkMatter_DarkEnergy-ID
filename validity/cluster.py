"""Cluster validity-boundary galaxies on mechanism variables.

No prior-cause labels are used: the data clusters itself on the 12
mechanism variables from mechanism.py, then clusters are profiled.
Model selection via GMM BIC + KMeans silhouette + bootstrap ARI.

Self-contained: needs scikit-learn + matplotlib. Run from this folder:
    python cluster.py
Reads results/mechanisms.json, writes results/clusters.png and
results/repair_plane.png, prints cluster profiles.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.decomposition import PCA

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

with open(os.path.join(RES, "mechanisms.json")) as f:
    data = json.load(f)
mechs = data["mechanisms"]
print(f"{len(mechs)} boundary galaxies; data-poor: {data['datapoor']}")

names = [m["name"] for m in mechs]
dKap = np.array([m["d_logkap"] if np.isfinite(m["d_logkap"]) else 3.0
                 for m in mechs])
X = np.column_stack([
    np.log10(np.array([m["rf"] for m in mechs]) + 1e-6),
    np.array([m["neg_frac"] for m in mechs]),
    np.array([m["log_abs_t2min"] for m in mechs]),
    np.log10(np.abs(np.array([m["frac_A"] for m in mechs])) + 1e-6),
    np.array([1.0 if m["B_neg"] else 0.0 for m in mechs]),
    dKap,
    np.array([1.0 if m["kap_feasible"] else 0.0 for m in mechs]),
    np.array([m["d_logrho"] for m in mechs]),
    np.array([m["lrhoB"] for m in mechs]),
    np.array([m["lr0"] for m in mechs]),
    np.log10(np.array([m["vmax"] for m in mechs])),
    np.array([m["log_rho0_over_target"] for m in mechs]),
])
Xs = StandardScaler().fit_transform(X)

print("\n=== model selection ===")
for k in (2, 3, 4):
    gm = GaussianMixture(n_components=k, n_init=10, random_state=0).fit(Xs)
    km = KMeans(n_clusters=k, n_init=20, random_state=0).fit(Xs)
    print(f"k={k}: GMM BIC={gm.bic(Xs):.1f} | "
          f"KMeans silhouette={silhouette_score(Xs, km.labels_):.3f}")

rng = np.random.default_rng(0)
print("\n=== bootstrap stability (ARI) ===")
for k in (2, 3):
    base = KMeans(n_clusters=k, n_init=20, random_state=0).fit(Xs).labels_
    aris = [adjusted_rand_score(
        base[idx], KMeans(n_clusters=k, n_init=10,
                          random_state=b).fit(Xs[idx]).labels_)
        for b in range(50)
        for idx in [rng.integers(0, len(Xs), len(Xs))]]
    print(f"k={k}: mean ARI={np.mean(aris):.3f} ± {np.std(aris):.3f}")

km = KMeans(n_clusters=2, n_init=20, random_state=0).fit(Xs)
gm = GaussianMixture(n_components=2, n_init=10, random_state=0).fit(Xs)
print(f"\nKMeans vs GMM agreement (ARI): "
      f"{adjusted_rand_score(km.labels_, gm.predict(Xs)):.3f}")
labels = km.labels_

print("\n=== cluster profiles ===")
for c in range(2):
    m = [mechs[i] for i in np.where(labels == c)[0]]
    dk = [x["d_logkap"] for x in m if np.isfinite(x["d_logkap"])]
    print(f"\nCluster {c} (n={len(m)}): "
          f"{'Type A: curvature-limited' if not any(x['B_neg'] for x in m) else 'Type B: stabilizer-overloaded'}")
    print(f"  B_neg: {sum(1 for x in m if x['B_neg'])}/{len(m)}")
    if dk:
        print(f"  d_logkap: {np.mean(dk):.2f} ± {np.std(dk):.2f}")
    else:
        print("  d_logkap: all infeasible")
    dr = [x["d_logrho"] for x in m]
    print(f"  d_logrho: {np.mean(dr):.2f} ± {np.std(dr):.2f}")
    print(f"  neg_frac: {np.mean([x['neg_frac'] for x in m]):.3f}")
    print(f"  <lrhoB>={np.mean([x['lrhoB'] for x in m]):.2f} "
          f"<lr0>={np.mean([x['lr0'] for x in m]):.2f} "
          f"<vmax>={np.mean([x['vmax'] for x in m]):.0f}")
    print(f"  members: {', '.join(x['name'] for x in m)}")

pc = PCA(n_components=2).fit_transform(Xs)
fig, ax = plt.subplots(figsize=(8, 6))
for c in range(2):
    sel = labels == c
    ax.scatter(pc[sel, 0], pc[sel, 1], s=60, label=f"Cluster {c} (n={sel.sum()})")
for i, nm in enumerate(names):
    ax.annotate(nm, (pc[i, 0], pc[i, 1]), fontsize=6, alpha=0.7)
ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(); ax.grid(alpha=0.3)
ax.set_title("Validity-boundary clustering (PCA)")
fig.tight_layout()
fig.savefig(os.path.join(RES, "clusters.png"), dpi=150)

bn = np.array([m["B_neg"] for m in mechs])
dr = np.array([m["d_logrho"] for m in mechs])
dk_plot = np.array([3.2 if not np.isfinite(x) else x for x in dKap])
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(dk_plot[~bn], dr[~bn], s=60, label="stabilizer OK (B>0)", alpha=0.8)
ax.scatter(dk_plot[bn], dr[bn], s=60, marker="s",
           label="stabilizer overloaded (B<0)", alpha=0.8)
for i, nm in enumerate(names):
    ax.annotate(nm, (dk_plot[i], dr[i]), fontsize=6, alpha=0.7)
ax.set_xlabel(r"$\Delta\log_{10}\kappa$ to repair (3.2 = infeasible)")
ax.set_ylabel(r"$\Delta\log_{10}\rho_{\rm target}$ to repair")
ax.legend(); ax.grid(alpha=0.3)
ax.set_title("Repair-distance plane")
fig.tight_layout()
fig.savefig(os.path.join(RES, "repair_plane.png"), dpi=150)
print("\nWrote results/clusters.png and results/repair_plane.png")
