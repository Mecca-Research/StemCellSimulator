# From the paper to the simulator

Where each part of *Modeling Multi-Order Chemical Reactions and Integrated
Cellular Pathways* (2024) is implemented. Scenes are opened with
`index.html#<scene>` (and `?view=<view>` where a scene has several views).

| Paper section | What the paper specifies | Implementation | Scene |
|---|---|---|---|
| Abstract / Introduction | Integrate multiple signalling pathways across orders of complexity, from stem cells to tissues and organs | Whole application; reactions tagged by order in every model's `REACTIONS` | all |
| Methodology 1–4 (data acquisition, batch graphs, combining) | Build graphs of reactions and merge them | `js/models/network/registry.js` (`allReactions`, `buildGraph`) | `#network` |
| Methodology 5 (PCA) | Dimensionality reduction by PCA | `embed()` — first three principal components of the adjacency matrix (`powerPCA`), spring-refined | `#network` |
| Methodology 6–7 (visualisation, validation) | Visualise in segments; cross-reference with existing data | 3D scenes; calibration against real microscopy (`tools/build_assets.py`, `calibration.json`) | `#gallery`, `#colony` |
| **First order** — Notch | Ligand (Delta/Jagged) binding → receptor cleavage → NICD → nucleus → transcription | `pathways/notch.js` (`singleCell`: S2/S3 cleavage, NICD import, Hes1; `collier`: lateral inhibition) | `#signaling`, `#colony` |
| First order — Wnt | Wnt + Frizzled/LRP5/6 → Dishevelled → GSK-3β inhibition → β-catenin stabilisation → nucleus → target genes | `pathways/wnt.js` (`paperWnt` = the paper's ODE verbatim; `canonicalWnt` with destruction complex and Axin2 feedback) | `#signaling` |
| First order — Hedgehog | Shh + Patched → Smoothened → GLI | `pathways/hedgehog.js` | `#signaling`, `#morphogenesis` |
| Mathematical models — Wnt ODE | dW/dt = k1 − k2W, dF/dt = k3W − k4F, dD/dt = k5F − k6D, dβ/dt = k7D − k8β | `paperWnt.rhs`, closed-form `steadyState`, tested against the exact solution | `#signaling` |
| **Second order** — MAPK/ERK | RTK activation → Raf → MEK → ERK → nuclear transcription factors | `pathways/mapk.js` (Huang–Ferrell cascade, ultrasensitivity) | `#hematopoiesis` |
| Second order — JAK/STAT | Cytokine binding → receptor dimerisation → JAK → STAT phosphorylation & nuclear translocation | `pathways/jakstat.js` (EPO-R/JAK2/STAT5 with SOCS feedback) | `#hematopoiesis` |
| Hematopoietic example | HSC → HSC + HSC (cytokines); HSC → CMP (GATA1, PU.1); CMP → (EPO) proerythroblast → erythroblast → reticulocyte → erythrocyte; CMP/CLP lineage graph | `pathways/hematopoiesis.js` (GATA1–PU.1 toggle, lineage tree, stochastic populations) | `#hematopoiesis` |
| **Third order** — Myogenesis | MyoD & Myf5 → myogenin & MRF4 → muscle-specific genes | `pathways/myogenesis.js` (MRF cascade, fusion) | `#differentiation?view=myogenesis` |
| Third order — Hormonal feedback | Secretion → receptor binding → negative feedback | `pathways/hormone.js` (Goodwin loop, Hopf threshold) | `#differentiation?view=hormone` |
| Neuronal example | Stem cell →(Sox2, Pax6)→ progenitor →(Neurogenin, NeuroD)→ neuron; Neurogenin mRNA/protein ODEs; neurite outgrowth ∂C/∂t = D∇²C + R(C); neurite outgrowth, synaptogenesis, axonal guidance | `pathways/neurogenesis.js` (paper's ODE verbatim + analytic solution; guidance field; growth cones) | `#differentiation?view=neurogenesis` |
| **Fourth order** — Crosstalk, convergence, combinatorial regulation | Wnt–Notch crosstalk, multiple transcription factors | `pathways/crosstalk.js` | `#morphogenesis`, `#signaling` |
| Fourth order — Morphogen gradients | Gradients of signalling molecules guide patterning | `morpho/gradient.js` (French flag) | `#morphogenesis` |
| Morphogenesis — Turing patterns | ∂u/∂t = Du∇²u + f(u,v), ∂v/∂t = Dv∇²v + g(u,v) | `morpho/reactionDiffusion.js` (Schnakenberg/Gray–Scott on grids and meshes, dispersion relation) | `#morphogenesis` |
| Morphogenesis — mechanics | Cell migration and adhesion; Cell Movement = Chemotaxis + Haptotaxis | `morpho/agents3d.js` (adhesion/repulsion mechanics), `morpho/ecm.js` (chemotaxis + haptotaxis) | `#colony`, `#morphogenesis` |
| Fourth order — ECM interactions | Cells interact with the ECM to form structured tissues | `morpho/ecm.js`, collagen kinetics | `#morphogenesis`, `#proteins?view=collagen` |
| Epigenetic modifications | DNA + DNMT → methylated DNA; histone + HAT → acetylated; + HDAC → histone | `pathways/epigenetics.js` (`marks` with closed-form steady state) | `#landscape` |
| Stem cell differentiation — gene expression regulation | Transcription factors, signalling, epigenetics | `pathways/epigenetics.js` (`circuit`, quasi-potential `landscape`) | `#landscape` |
| Constructing complex proteins — hemoglobin | DNA → mRNA → α/β globin → 2α + 2β → tetramer → + 4 heme | `protein/assembly.js` (`hemoglobin`, exact α balance) | `#proteins` |
| Constructing complex proteins — collagen | COL1A1/COL1A2 → pro-α chains → hydroxylation + glycosylation → 2 α1 + 1 α2 triple helix → secretion → propeptide cleavage → tropocollagen → fibril | `protein/assembly.js` (`collagen`, exact chain balance) | `#proteins?view=collagen` |
| Creating a cellular map | Specialised cell types by lineage and function (Sections 1–44) | `atlas/atlas.js` (164 cell types, 17 lineages, germ layers) | `#atlas` |
| Next order of complex structures | Tissues → organs → organ systems (tree) | `atlas/atlas.js` (`TISSUES`, `ORGANS`, `SYSTEMS`) | `#atlas?view=hierarchy` |
| Higher-level interactions | Cardiovascular–respiratory, nervous–muscular, endocrine–digestive, immune–lymphatic, urinary–cardiovascular, reproductive–endocrine, nervous–endocrine | `atlas/atlas.js` (`SYSTEM_LINKS`) | `#atlas?view=systems` |
| Organ example — kidney | Nephrons (podocytes, proximal/distal tubule, collecting duct), renal corpuscle | `organ/nephron.js` + real 3D kidney volume | `#kidney` |
| Order of chemical reactions (~500 / 1,051 / 2,051 / 3,021) | Reported cumulative counts per order | Shown next to the counts of explicitly modelled reactions (not conflated) | `#network` |

## Beyond the paper

- **Real 3D microscopy.** A 3D confocal stack of real cells is segmented into
  cells and nuclei, meshed, and used to seed the digital-twin colony
  (`#colony`).
- **Calibration.** Cell/nucleus volumes, monolayer height, spacing, mitotic
  index, nuclear-import kinetics, erythrocyte size and shape, and organoid size
  distribution are measured on real images (`assets/derived/calibration.json`).
- **Render modes.** Every scene can be viewed as a physical 3D rendering, as a
  confocal fluorescence image, or as a virtual H&E section.
