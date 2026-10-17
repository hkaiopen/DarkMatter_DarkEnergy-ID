# Validity analysis: where the inversion holds, and why it stops

The K(r) inversion requires (θ')² ≥ 0 everywhere on the radial grid.
This folder studies the 46 SPARC galaxies where the fixed-parameter
inversion does not satisfy this condition — not as "failures", but as
probes of the method's **domain of validity**.

## Two mechanism types (data-driven)

Clustering on mechanism variables derived from the field equations
(no prior-cause labels used) yields two clean types (k=2; GMM and
KMeans agree exactly):

- **Type A — curvature-limited (26 galaxies).** The Burkert profile's
  nonzero central derivative drives a 1/r divergence in the curvature
  term A(r); a modest increase in κ (~3×) restores validity. Half are
  cured by moving the inner grid cutoff from 0.02 to 0.05 kpc.
  Interpretation: the fixed-κ stabilizing feedback is *marginally*
  insufficient against the profile's central curvature.

- **Type B — stabilizer-overloaded (19 galaxies).** Central density
  exceeds ρ_target, flipping the stabilizer sign (B<0); provably no κ
  increase can repair it (linear infeasibility) — ρ_target must be
  raised. Zero are cured by moving the inner cutoff.
  Interpretation: *projection overload* — the observed central density
  exceeds what the fixed-parameter virtual projection can stabilize.

## Contents

- `REPORT.md` — full analysis report (mechanism, clustering, r_min test,
  information-dynamics reading)
- `mechanism.py` — compute the 17 mechanism variables + exact repair
  distances (κ/ρ_target) via linear feasibility; self-contained
- `cluster.py` — GMM/KMeans clustering, model selection (BIC/silhouette),
  bootstrap stability; self-contained
- `results/` — `mechanisms.json` (data), `clusters.png`,
  `repair_plane.png` (figures)

## Run

```bash
pip install numpy scipy scikit-learn matplotlib
python mechanism.py --data ../data --out results
python cluster.py
```

## Key numbers

- 175 SPARC galaxies → 126 valid inversions at r_min=0.02 kpc
  (139 at r_min=0.05 kpc)
- 46 validity-boundary galaxies → 26 Type A + 19 Type B (+1 marginal)
- 4 galaxies have too few data points for any fit (listed in REPORT.md)
