# Validity-boundary analysis: 46 SPARC galaxies beyond the fixed-parameter inversion domain

## Core result (one sentence)

The 46 galaxies where the fixed-parameter K(r) inversion violates (θ')² ≥ 0
are not random outliers but **two mechanistically clean, data-coherent types**:
**Type A (26)** — curvature-limited, repaired by a modest κ increase;
**Type B (19)** — stabilizer-overloaded, unrepairable by κ, requiring higher
ρ_target. (+1 marginal boundary case; 4 galaxies have too few data points.)

## 1. Mathematical root of the boundary

The φ² numerator: N(r) = A(r) + κ·C(r), with
- A(r) = R'' + 2R'/r (geometric/curvature term)
- C(r) = (ρ_target − R²)·R (stabilizer coefficient)

**Key fact**: the Burkert profile has R'(0) = −√ρ₀/(2r₀) ≠ 0 (nonzero central
derivative — strictly speaking not a true flat core). Hence 2R'/r ~ −√ρ₀/(r₀·r)
diverges as 1/r toward −∞ at the center.

Consequences:
- For *any* Burkert halo, φ² < 0 sufficiently close to r = 0 — an intrinsic
  property of the Burkert parametrization, not of the physics.
- On the numerical grid (r_min = 0.02 kpc), validity is a **quantitative race**
  between the −1/r-divergent A(r) and the κ·C(r) stabilizer.
- 126 galaxies: stabilizer wins at r = 0.02. 46 galaxies: it does not.
- Resolution test (n=2000 vs n=8000): identical results — not a discretization artifact.

## 2. Data-driven clustering (k=2; GMM and KMeans agree exactly)

Twelve mechanism variables derived from the field equations (no prior-cause
labels), standardized, then clustered:

| Criterion | k=2 | k=3 | k=4 |
|-----------|-----|-----|-----|
| GMM BIC | **−1888.9** ✓ | −1823.9 | −1707.7 |
| KMeans silhouette | 0.474 | 0.548 | 0.445 |
| KMeans–GMM agreement (ARI) | **1.000** | — | — |
| Bootstrap ARI | 0.83 | 0.95 | — |

BIC selects k=2; the two algorithms return **identical partitions** —
the split is not an algorithmic artifact.

### Type A — curvature-limited (n=26)
- Stabilizer sign normal at failure point (B > 0) in 26/26
- κ-repairable: mean Δlog₁₀κ = +0.51 (~3× increase suffices)
- Small ρ_target adjustment: mean Δlog₁₀ρ_target = +0.21
- Negative φ² confined to 1.0% of radial range (thin central shell)
- Lower central densities: ⟨log₁₀ρ₀⟩ = 7.75

### Type B — stabilizer-overloaded (n=19)
- Stabilizer sign flipped at failure point (B < 0) in 19/19
  (central R² > ρ_target)
- **κ-unrepairable**: Δlog₁₀κ = ∞ (linear-feasibility proof — increasing κ
  worsens the C < 0 region)
- Requires ρ_target increase: mean Δlog₁₀ρ_target = +0.49
- Negative φ² spans 9.2% of radial range (substantially wider)
- Higher central densities: ⟨log₁₀ρ₀⟩ = 8.52

## 3. Inner-cutoff test: marginal vs robust (independent confirmation)

Relaxing the inner grid cutoff from r_min = 0.02 to 0.05 kpc:

| | Cured | Still out of domain |
|---|---|---|
| Type B (19) | **0** | 19 |
| Type A (26) | **13** | 13 |

- Type B: **zero cures** — robustly out of domain, grid-independent.
- Type A: half cured — 13 are marginal thin-shell cases, normal beyond r ≥ 0.05 kpc.

## 4. Information-dynamics reading

- **Type A** = stabilizing feedback *marginally insufficient*: the Burkert
  central 1/r singularity (a parametrization artifact) demands more phase
  curvature than fixed κ can coherently support. A small κ increase — or a
  genuinely flat-cored profile — restores oscillatory feedback.
  **Not a falsification; a parameter/representation issue.**
- **Type B** = *projection overload*: the observed central density exceeds
  the stability capacity of the fixed-parameter virtual projection. The phase
  is crushed to imaginary, feedback degrades from oscillatory to exponential.
  ρ_target (the virtual "stability capacity") must be raised; κ is useless.
  **The most physically informative boundary type.**

In one line: **valid inversion = virtual space can compile the halo into a
stably oscillating negative-feedback vortex; Type A = compilation resolution
marginally insufficient; Type B = the halo's density exceeds the virtual
stability capacity.**

## 5. Implications

- The 126/175 success rate is the quantitative result at r_min = 0.02 kpc
  (139/175 at r_min = 0.05 kpc).
- The 19 Type-B galaxies motivate galaxy-dependent (κᵢ, ρ_target,i) —
  the natural entry point for hierarchical Bayesian modeling.
- Replacing Burkert (R'(0) ≠ 0) with a genuinely flat-cored profile is
  predicted to cure most Type-A cases — a testable follow-up.
- Data-poor galaxies (too few points for any fit): D512-2, NGC6789,
  UGC00634, UGC07232 — set aside, not analyzed.

## Files

- `results/mechanisms.json` — 46 galaxies × 17 mechanism variables
- `results/clusters.png` — PCA cluster scatter
- `results/repair_plane.png` — repair-distance plane (Δlog₁₀κ vs
  Δlog₁₀ρ_target, colored by stabilizer sign)
- `mechanism.py` — mechanism variables + exact repair distances via linear
  feasibility (self-contained: needs only numpy, scipy)
- `cluster.py` — clustering + model selection + bootstrap (self-contained:
  needs scikit-learn, matplotlib)
