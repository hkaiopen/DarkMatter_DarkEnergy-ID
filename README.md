# Dark Matter & Dark Energy: Inverting the Coupling from SPARC

Inverting a universal virtual–real coupling $K(r)$ from SPARC galaxy rotation curves.

## Result

Across 126 of 175 SPARC galaxies, the inverted coupling is universal:

$$K(r) = -\frac{K_0}{1+(r/r_c)^2}$$

- $K_0 = 204 \pm 31\ \mathrm{kpc^{-2}}$ (15.3% scatter)
- $r_c = 0.145 \pm 0.012\ \mathrm{kpc}$ (8.4% scatter)
- Lorentzian fit $R^2 > 0.91$ for every successful galaxy

The core size $r_c$ matches the 8-galaxy pilot value (0.144 kpc) exactly —
universality confirmed at 16× the sample size.

New: $K_0$ and $r_c$ are anti-correlated (Pearson $r=-0.66$). The integrated
proxy $K_0 r_c^2$ (13.8% scatter) is more universal than $K_0$ alone,
suggesting a conserved vortex charge: cores vary from "peaked-narrow"
to "broad-shallow" at fixed total charge.

Interpretation: each dark matter halo is anchored by a virtual vortex — $r_c$ is the
virtual healing length, $K_0$ the topological charge. Dark matter (localized vortices)
and dark energy (homogeneous background) are two scales of one immersive projection.

![K(r) inversion](figs/K_inversion_fits.png)
![Universality](figs/K_universality.png)

## Usage

```bash
pip install -r requirements.txt
```

Sample validation (8 characteristic galaxies, ~1–2 min):

```bash
python src/sample_validation.py --data data --out figs
```

Full validation (all 175 SPARC galaxies, ~15–30 min):

```bash
python src/full_validation.py --data data --out results
```

## Layout

- `src/invert_K.py` — inversion core (static field equations → $K(r)$)
- `src/sample_validation.py` — 8-galaxy sample validation
- `src/full_validation.py` — full-sample universality stress test
- `data/` — SPARC `*_rotmod.dat` (175 galaxies)
- `figs/` — result figures
- `results/` — full-sample output (`k_results.json`, `K_universality_full.png`)

## Data

SPARC: Lelli, McGaugh & Schombert 2016, AJ 152, 157. Data via Zenodo 16284118.

## License

CC BY-NC-SA 4.0
