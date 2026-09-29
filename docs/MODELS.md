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
Quasi-steady Shh:PTCH1 binding f_b = Shh/(Shh + K_d) with ligand-induced
PTCH1 internalisation (Incardona et al. 2000); catalytic inhibition of SMO by
free PTCH1 (Taipale et al. 2002): dS/dt = k_Sa(1 − S) − (k_Si + k_Ptc(1 − f_b)P)S;
GLI processing into activator (SMO-driven) or GLI3 repressor; Shea–Ackers
promoter φ = (β + (A/K_A)²)/(1 + (A/K_A)² + (R/K_R)²) (Lai, Robertson &
Schaffer 2004; Saha & Schaffer 2006); Ptch1 as a GLI target gives negative
feedback and the temporal adaptation of Dessaud et al. 2007. Cyclopamine sets
k_Sa = 0. French-flag helpers with neural-tube domains (FoxA2 / Nkx2.2 /
Olig2 / Pax6-7). Parameters are from a small search, not a fit.

## Second order
- `pathways/mapk.js`: Huang & Ferrell 1996 (*PNAS* 93:10078) cascade with
  their concentrations (MAPKKK 3 nM, MAPKK and MAPK 1.2 µM, phosphatases) and
  K_m = 0.3 µM, k_cat = 150 s⁻¹, in competitive Michaelis–Menten form with exact
  per-tier steady states; effective Hill coefficient nH = ln 81 / ln(EC90/EC10)
  = 1.0 / 3.1 / 14 for MAPKKK / MAPKK / MAPK. HF96's mass-action model (with
  enzyme–substrate sequestration) gives 1.0 / 1.7 / 4.9: same ordering, less
  steep. Receptor (SCF/c-Kit), ERK nuclear cycling and Elk-1 → c-Fos are
  illustrative.
- `pathways/jakstat.js`: EPO receptor → JAK2 → STAT5 phosphorylation,
  dimerisation, nuclear import, dephosphorylation and export (structure of
  Swameye et al. 2003, *PNAS* 100:1028) with CIS/SOCS feedback (Vera et al.
  2008, *BMC Syst Biol* 2:38); total STAT5 conserved to 1e-9; illustrative
  parameters.
- `pathways/hematopoiesis.js`: GATA1–PU.1 cross-antagonism with
  self-activation (Huang et al. 2007, *Dev Biol* 305:695; three attractors at
  a = 1, pitchfork at a_c = 0.774), EPO / GM-CSF inputs (Chickarmane et al.
  2009), Langevin decisions; the paper's lineage tree with markers, sizes and
  stage durations; Gillespie compartment kinetics with niche-limited HSC
  self-renewal a(1 − H/K) (Schofield 1978; Marciniak-Czochra et al. 2009),
  H* = K(1 − 1/2a). Survival rules (EPO-dependent erythroid survival, SCF
  dependence of c-Kit⁺ progenitors) are phenomenological. Note: live imaging
  (Hoppe et al. 2016, *Nature* 535:299) indicates early lineage choice is not
  initiated by stochastic PU.1:GATA1 fluctuations; the toggle is presented as a
  model, not as settled mechanism.

## Third order
- `pathways/myogenesis.js` (see also Münsterberg 1995, Gustafsson 2002, Millay 2013): MRF hierarchy (Bentzinger, Wang & Rudnicki 2012) —
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
- `morpho/reactionDiffusion.js`: 2D/3D grid Laplacians (zero-flux or
  periodic) and a cotangent mesh Laplacian with mixed Voronoi areas (Meyer et
  al. 2003; < 1 % error against −l(l+1)/R² on a sphere); spherical harmonics
  and angular power spectra; Schnakenberg (1979), Gray–Scott (Pearson 1993
  presets) and Gierer–Meinhardt (1972) kinetics; linear Turing analysis
  (steady state, Jacobian, dispersion relation Re λ(k²), the four Turing
  conditions, unstable band, fastest mode; Murray, *Mathematical Biology II*),
  with sphere modes k² = l(l+1)/R² — on the organoid the predicted l* = 7
  matches the simulated dominant mode; explicit and semi-implicit (conjugate
  gradient) solvers.
- `morpho/gradient.js`: C = C₀e^(−x/λ), λ = √(D/k), finite-domain cosh form and
  numerical 1D/3D solvers; French-flag boundaries x = λ ln(C₀/T) (Wolpert
  1969); positional error σₓ = λ·CV (Gregor et al. 2007; Bollenbach et al.
  2008); accumulation time τ(x) = (1 + x/λ)/(2k) (Berezhkovskii et al. 2010);
  default D and k from the measured Dpp gradient (Kicheva et al. 2007); Shh/BMP
  neural-tube read-out.
- `morpho/ecm.js`: procollagen → tropocollagen → fibril (mass-conserving;
  nucleation–elongation after Kadler et al.), fibroblast persistent random walk
  (Dunn & Brown 1987) laying fibres along their motion (Canty & Kadler 2004),
  and `MigrationSim` v = χ/(1 + αC)∇C + ρ∇E + noise with MMP degradation
  (Anderson & Chaplain 1998) — the paper's "Cell Movement = Chemotaxis +
  Haptotaxis".
- `pathways/crosstalk.js`: Dishevelled ⊣ GSK-3β ⊣ β-catenin and ligand → NICD
  → Hes1 ⊣ Atoh1 → Dll1 per cell; crosstalk by Dishevelled–NICD sequestration
  (Axelrod et al. 1996) and Wnt-induced Jagged1 (Rodilla et al. 2009); fate
  from a thermodynamic two-site enhancer (Shea & Ackers 1985; Bintu et al.
  2005) giving stem / Paneth / absorptive / goblet fates (Fre et al. 2005;
  van Es et al. 2005). Reproduces the γ-secretase-inhibitor (all secretory)
  and constitutive-NICD (no secretory) phenotypes. Rates are dimensionless.

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
