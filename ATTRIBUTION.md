# Attribution and licences

The simulator's source code is MIT-licensed (see `LICENSE`). Everything below
keeps its own licence.

## Reference microscopy redistributed in this repository

Every image or volume below is downloaded by `tools/fetch_micrographs.py` from
the URL pinned in `assets/micrographs/catalog.json`, verified against its
SHA-256, and turned into the files in `assets/derived/` by
`tools/build_assets.py`. Only the derived files are committed.

| Specimen | Credit | Licence | Used for |
|---|---|---|---|
| 3D confocal stack of a cell monolayer (membrane + DNA), `cells3d.tif` | Allen Institute for Cell Science, via the scikit-image data repository | CC0 1.0 | Digital-twin colony: 3D cell & nucleus meshes, volume rendering, calibration |
| Human cells in mitosis (DNA), CellProfiler *ExampleHuman* | David Root / Jason Moffat; Moffat et al. 2006, *Cell* 124:1283 | CC0 1.0 | Mitotic index → cell-cycle length; animated self-renewal |
| Normal epidermis and dermis, H&E 10× | Kilbad, Wikimedia Commons | Public domain | Virtual-H&E reference; animated keratinocyte differentiation; 3D epidermis |
| Mouse kidney tissue, 3D confocal | Genevieve Buckley, Monash Micro Imaging (2018), via scikit-image | CC0 1.0 | Organ-level volume next to the nephron model |
| Live imaging of nuclear-envelope targeting (`NPCsingleNucleus.tif`) | Andrea Boni & Jan Ellenberg; Boni et al. 2015, *J Cell Biol* 209:705 | CC0 1.0 | Fitted nuclear-translocation kinetics |
| Single cell, quantitative phase (`cell.png`) | Paul Müller et al. 2018, *Optics Express* 26:10729, via scikit-image | CC0 1.0 | Optical-thickness surface |
| Blood smears (images 1 and 2 + labels), Chula-RBC-12-Dataset | Chula PIC Lab; Naruenatthanaset et al. 2021, arXiv:2012.01321 | MIT | Erythrocyte internal ruler; absorbance profile vs 3D Evans–Fung model; lymphocyte size |
| Human intestinal organoids (bright field), OrgaQuant test data | Timothy Kassis; Kassis et al. 2019, *Sci Rep* 9:12479 | MIT | Organoid size distribution; animated lumen inflation |

## Further references (catalogued, not redistributed)

`assets/micrographs/references.json` lists 46 further openly licensed
micrographs for pluripotent colonies, mitosis, hematopoiesis, neurogenesis,
myogenesis, organoids, gastrulation, Turing and morphogen patterns, ECM and
Notch lateral-inhibition mosaics (Wikimedia Commons, eLife, GitHub-hosted
datasets). Their licence fields for Wikimedia Commons came from search-index
summaries and **must be re-confirmed on each file page** before reuse; three
EmbedSeg archives are CC BY-NC (non-commercial only). They are not downloaded
by the tools.

## Software

- [three.js](https://threejs.org) r186 and the addons in `vendor/three/` — MIT
  (`vendor/three/LICENSE`).
- Data pipeline: NumPy, SciPy, scikit-image, Pillow, tifffile (BSD/HPND
  licences; installed from PyPI, not vendored).
- Tests: Playwright (Apache-2.0, dev dependency).

## Source research

The four-order structure, pathways, reaction schemes, equations, cell lists and
body hierarchy implemented here come from *Modeling Multi-Order Chemical
Reactions and Integrated Cellular Pathways* (Startonix, 18 October 2024). The
established models used to make each step computable are cited in the source
files and in `docs/MODELS.md`.
