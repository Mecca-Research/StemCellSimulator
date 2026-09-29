# Models

Every model lives in `js/models/` as plain JavaScript (no rendering code),
uses the seeded RNG in `core/rng.js`, and is unit-tested in `tests/models/`.
Each section states what is an **established model** (with its reference) and
what is **phenomenological** (a qualitative choice made so the step can be
simulated). Parameters are dimensionless unless units are given.

## Numerics — `core/`

- `ode.js`: classical RK4; adaptive Dormand–Prince RK5(4) for stiff kinetics
  (protein assembly, MAPK); Gillespie's direct stochastic simulation algorithm;
  Hill functions.
- `linalg.js`: Jacobi eigen-decomposition, full PCA, and power-iteration PCA
  for large adjacency matrices.
- `spatial.js`: uniform spatial hash with exact integer bucket keys.

## 3D cell mechanics — `morpho/agents3d.js`

Centre-based, overdamped off-lattice model (Drasdo & Höhme 2005; Meineke,
Osborne et al. / Chaste):

```
η dxᵢ/dt = Σⱼ Fᵢⱼ + F_ext + √(2ηkT) ξ(t)
Fᵢⱼ = μ_rep s ln(1 + (d − s)/s)                 d < s   (repulsion)
Fᵢⱼ = μ_adh(τᵢ, τⱼ) (d − s) exp(−5 (d − s)/s)     s ≤ d < 1.25 s (adhesion)
```

Cells grow through G1–S–G2–M (phase fractions 0.45/0.30/0.19/0.06), divide
along an axis (planar for monolayers), and resolve mitosis into prophase …
cytokinesis for rendering. `dims: 2` models a monolayer of fixed height
(lateral radius ∝ √V). Contact inhibition arrests G1 cells with ≥ N neighbours.
Tests: relaxation to contact distance, adhesion, doubling once per cycle,
substrate confinement.

**Calibration from real data** (`#colony`): mechanical radius = ½ × measured
nearest-neighbour spacing of real nuclei (13.0 µm); each founder scaled by its
measured volume; monolayer height 11.6 µm; cell-cycle length T = t_M / MI with
MI = 0.062 measured on the CellProfiler mitosis field and t_M = 1 h
(adjustable) → T ≈ 16 h.

## First order

### Wnt — `pathways/wnt.js`
- `paperWnt`: **the paper's ODE verbatim** — dW/dt = k1 − k2W, dF/dt = k3W −
  k4F, dD/dt = k5F − k6D, dB/dt = k7D − k8B — with closed-form steady state
  W* = k1/k2, F* = k3W*/k4, D* = k5F*/k6, B* = k7D*/k8 and exact W(t).
- `canonicalWnt`: the paper's steps with the destruction complex and Axin2
  negative feedback in the spirit of Lee et al. 2003 (*PLoS Biol* 1:e10),
  reduced form (phenomenological rate constants).

### Notch — `pathways/notch.js`
- `collier`: lateral inhibition, Collier et al. 1996 (*J Theor Biol* 183:429):
  dNᵢ/dt = f(⟨D⟩ᵢ) − Nᵢ, dDᵢ/dt = v(g(Nᵢ) − Dᵢ), f(x) = xᵏ/(a + xᵏ),
  g(x) = 1/(1 + bxʰ); runs on any contact graph (the real 3D contacts in
  `#colony`). Test: isolated senders on a hexagonal sheet.
- `singleCell`: trans Delta binding → ADAM10 S2 cleavage → γ-secretase S3
  cleavage → NICD → nuclear import → Hes1 → repression of cis-Delta. Test:
  DAPT (γ-secretase = 0) abolishes NICD and Hes1.

### Hedgehog — `pathways/hedgehog.js`
Shh–PTCH1 binding, PTCH1 inhibition of SMO (Taipale et al. 2002), GLI
activator/repressor balance and PTCH1 negative feedback; French-flag readout
helper. See the file header for the references used.

## Second order
- `pathways/mapk.js`: Huang & Ferrell 1996 (*PNAS* 93:10078) cascade in
  Michaelis–Menten form; effective Hill coefficient nH = ln 81 / ln(EC90/EC10).
- `pathways/jakstat.js`: EPO receptor → JAK2 → STAT5 phosphorylation,
  dimerisation, nuclear cycling and SOCS/CIS feedback; total STAT conserved.
- `pathways/hematopoiesis.js`: GATA1–PU.1 cross-antagonism with
  self-activation (Huang et al. 2007, *Dev Biol* 305:695), Langevin decisions,
  the paper's lineage tree and stochastic population kinetics.

## Third order
- `pathways/myogenesis.js`: MRF hierarchy (Bentzinger, Wang & Rudnicki 2012) —
  Myf5/MyoD (bistable MyoD autoregulation, exact fold points), myogenin
  (repressed by growth factors), MRF4, MHC; myoblast culture with nematic
  alignment and volume/nuclei-conserving fusion. Kinetic forms are
  phenomenological Hill functions.
- `pathways/neurogenesis.js`: **the paper's Neurogenin ODEs verbatim**
  (analytic solution and steady state k_tx·k_tl/(k_dm·k_dp)); Sox2/Pax6 →
  Neurogenin → NeuroD with Notch/Hes1 lateral inhibition; guidance field =
  exact steady state of ∂C/∂t = D∇²C − λC + Σ Q δ (the paper's R(C) taken as
  linear decay plus point sources): C(r) = Σ Q e^(−r/ℓ)/(4πDr), ℓ = √(D/λ);
  growth-cone steering with gradient-sensing noise (Goodhill & Urbach).
- `pathways/hormone.js`: Goodwin (1965) loop with Griffith's oscillation
  condition (Hill n > 8 for equal clearance); Routh–Hurwitz Hopf threshold.

## Fourth order
- `morpho/reactionDiffusion.js`: grid and triangle-mesh Laplacians;
  Schnakenberg kinetics with linear Turing analysis (dispersion relation,
  unstable band; Murray, *Mathematical Biology II*); Gray–Scott (Pearson 1993).
- `morpho/gradient.js`: source–diffusion–degradation gradients,
  C(x) = C₀e^(−x/λ), λ = √(D/k); French-flag thresholds (Wolpert).
- `morpho/ecm.js`: collagen deposition and chemotaxis + haptotaxis
  (Anderson & Chaplain 1998) — the paper's "Cell Movement = Chemotaxis +
  Haptotaxis".
- `pathways/crosstalk.js`: Wnt–Notch crosstalk and combinatorial fate logic of
  the intestinal crypt (Fre et al. 2005; van Es et al. 2005).

## Epigenetics and the landscape — `pathways/epigenetics.js`
- The paper's reactions as mark fractions: dm/dt = k_M·DNMT·(1 − m) − k_DM·m,
  da/dt = k_A·HAT·(1 − a) − k_D·HDAC·a (closed-form steady states).
- Fate circuit (Huang et al. 2007): three attractors (multipotent + two fates)
  with strong self-activation, two after the bifurcation.
- Landscape U = −ln P_ss from Langevin trajectories (Wang, Xu & Wang 2008,
  *PNAS* 105:12271). Test: landscape minima coincide with the deterministic
  attractors.
- **Phenomenological coupling:** acetylation sets the noise (plasticity),
  methylation the repression strength (canalisation).

## Protein construction — `protein/assembly.js`
- **Hemoglobin**, the paper's scheme: mRNA → α/β → αβ → α₂β₂ → + 4 heme,
  plus AHSP chaperoning of free α (Kihm et al. 2002), precipitation of
  unpaired α (β-thalassaemia), HRI coupling of translation to heme (Chen 2007).
  Exact α-chain mass balance is tested.
- **Collagen I**, the paper's scheme: pro-α1/pro-α2 → hydroxylation +
  glycosylation (ascorbate-dependent prolyl-4-hydroxylase) → 2 α1 + 1 α2 →
  procollagen → secretion → ADAMTS2/BMP1 cleavage → tropocollagen → fibril
  (nucleation–elongation; Kadler et al. 1996). Scurvy (no ascorbate) and
  dermatosparaxis (no ADAMTS2) block fibrils; exact chain balance is tested.

## Organ level — `organ/nephron.js`
Nephron segments with textbook water-reabsorption fractions (proximal ≈65 %,
thin descending ≈15 %, thick ascending 0 %, distal ≈5 %, collecting duct ADH-
regulated); monotone filtrate profile, tested.

## Cellular map — `atlas/atlas.js`
164 cell types from the paper's lists in 17 lineages and 5 germ-layer groups,
with tissues, organs, organ systems and the paper's system interactions. Where
the paper's grouping differs from developmental origin (osteoclasts and
microglia are myeloid; Schwann cells, melanocytes and chromaffin cells come
from the neural crest) this is recorded in the data.

## Reaction network — `network/registry.js`
Aggregates every module's `REACTIONS`, builds the bipartite species–reaction
graph, embeds it by PCA of the (two-step diffused) adjacency matrix and a
spring refinement, and reports counts per order next to the paper's reported
cumulative estimates.

## Measurements on real microscopy — `tools/build_assets.py`

| Quantity | Value | Method |
|---|---|---|
| Nucleus volume (median, whole nuclei) | 756 µm³ | 3D Otsu + distance-transform h-maxima watershed |
| Cell volume (median, interior) | 3,188 µm³ | membrane-intensity watershed seeded by nuclei |
| Nucleus : cytoplasm volume | 0.39 | per cell |
| Monolayer height | 11.6 µm | z-extent of segmented cells |
| Nearest-neighbour spacing | 13.0 µm | nucleus centroids, in plane |
| Mitotic index | 0.062 (19 / 306) | 90th-percentile DNA intensity > median + 4 MAD |
| Nuclear-envelope targeting | E₀ 1.19 → E∞ 3.02, τ = 13.5 frames, R² 0.955 | rim/cytoplasm ratio, first-order fit |
| Erythrocyte profile vs Evans–Fung 3D model | r = 0.98, RMSE 0.063 | median radial absorbance of 59 isolated cells |
| Lymphocyte nucleus | 9.2 µm | RBC (7.82 µm) as internal ruler |
| Intestinal organoids | 91 detected, radius CV 0.33 | Canny + Hough circles |
