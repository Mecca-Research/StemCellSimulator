# Architecture

The simulator is a static web application (ES modules, no build step) plus a
Python data pipeline. Open `index.html` through any static server
(`npm run serve`) — WebGL 2 is required.

```
index.html            app shell (import map points "three" at vendor/three)
css/app.css           UI theme
js/main.js            scene registry, routing (#scene-id), transport, picking
js/engine/            rendering + UI toolkit shared by every scene
js/models/            pure simulation code (no DOM, no three.js) - unit tested in Node
js/scenes/            one module per 3D scene
assets/micrographs/   catalog.json (redistributed, SHA-pinned sources), references.json
assets/derived/       committed outputs of tools/build_assets.py
tools/                fetch_micrographs.py, build_assets.py (Python, scikit-image)
tests/models/         node:test unit tests for js/models
tests/tools/          pytest for the data pipeline
tests/e2e/            headless-Chromium smoke test of every scene
vendor/three/         three.js r186 (MIT) + the addons in use
```

## Scene contract

`js/scenes/index.js` lists scenes; each is loaded lazily and must default-export:

```js
export default {
  about: '<p>HTML shown at the top of the right panel</p>',
  paperRef: 'Which part of the paper this implements',
  async create(ctx) {
    // ctx.THREE, ctx.stage, ctx.root (THREE.Group - add everything here),
    // ctx.panel (control builder), ctx.assets (loaders), ctx.legend(items),
    // ctx.setStatus(text) (transport bar readout), ctx.query (URLSearchParams
    // from "#scene-id?view=x" deep links - use it to pick the initial view)
    return {
      update(dt, t) {},  // dt = seconds of animation * speed (0 while paused)
      reset() {},        // optional: restart the simulation in place
      onMode(mode) {},   // optional: 'physical' | 'confocal' | 'histology'
      dispose() {},      // optional: free textures etc. (meshes/materials under root are freed for you)
    };
  },
};
```

Rules: add objects only under `ctx.root`; keep simulation state in plain JS;
every chart/readout goes through `ctx.panel`; never touch the DOM outside the
panel; unit conversions (sim time per second) are chosen by the scene and shown
with `ctx.setStatus`.

## Engine API (js/engine)

| module | purpose |
|---|---|
| `materials.js` | `cellMaterial({ role, color, stain, useInst, opacity, rimPower, noise, noiseScale, emissive, gain, he, vertexColors, side })` - the uber-shader used for all biological surfaces. Roles: `membrane`, `cytoplasm`, `nucleus`, `reporter`, `solid`, `molecule`, `rbc`, `field`. Works with `InstancedMesh` (instance colours multiply `color`) and vertex colours. The three render modes (physical / confocal / virtual H&E) switch globally. `lineMaterial({ color, opacity, vertexColors })` for lines. |
| `instances.js` | `InstancePool(geometry, material, capacity, parent)`: `begin()`, `put(x,y,z,sx,sy,sz,color)`, `putOriented(pos, axis, sAlong, sPerp1, sPerp2, color)`, `putMatrix(m, color)`, `end()`. Grows automatically (the underlying `mesh` object is replaced - use a getter when registering it for picking). `idColor(id)` stable hues. |
| `geometry.js` | `unitSphere`, `erythrocyteGeometry` (Evans-Fung), `lobedNucleus`, `indentedNucleus`, `blobGeometry`, `foldedChainGeometry`, `collagenStrandCurve`, `decodeMesh`, `mergeGeometries`. |
| `volume.js` | `VolumeView(tex3d, extent, opts)` ray-marched volume (MIP / compositing / virtual H&E), `mosaicToVolume`, `mosaicToLabels`. |
| `stage.js` | `stage.frame(center, radius, dir)`, `stage.clipAxis` / `stage.clipRange` (section plane), `stage.setPickables([{ object, info(hit) => string|null }])`, `stage.setBloomScale(s)` (reset to 1 on every scene change), `stage.mode`, `stage.camera`, `stage.controls`. |
| `ui.js` | `panel.section(title)`, `slider({label,min,max,step,value,onChange,format,log})`, `toggle`, `select({label,options,value,onChange})`, `buttons([{label,onClick,primary}])`, `readouts(keys) -> {set(k,v)}`, `chart({title,series:[{name,color,dash}],yRange,logY,xLabel,capacity}) -> Chart` (`push(x, values)`, `set(xs, columns)`, `clear()`), `equation(text)`, `note(html)`, `table(headers, rows)`. |
| `assets.js` | `loadJSON(name)`, `loadBinary(name)`, `loadImageData(name)`, `loadImage(name)` for files in `assets/derived/`. |

## Models (js/models)

Pure functions/classes with deterministic RNG (`core/rng.js`). Every pathway
module exports a `REACTIONS` array
`{ id, order: 1..4, pathway, reactants: [], products: [], modifiers?: [], label }`
that the reaction-network scene aggregates.

- `core/`: `ode.js` (RK4, adaptive Dormand-Prince, Gillespie SSA, Hill functions),
  `rng.js`, `spatial.js` (spatial hash), `linalg.js` (Jacobi eigensolver, PCA).
- `morpho/agents3d.js`: centre-based 3D cell mechanics with cell cycle, mitosis
  stages, oriented division, monolayer (`dims: 2`) or volumetric growth.
- `pathways/*.js`: one module per signalling / differentiation programme.

## Render modes

- **3D (physical)**: wrap-lit translucent membranes, textured nuclei, bloom off.
- **Confocal**: emission only - membrane dyes glow at grazing angles (as in an
  optical section), DNA dyes fill nuclei with chromatin texture, additive
  blending, bloom.
- **H&E**: virtual hematoxylin & eosin by Beer-Lambert absorption with the
  Giacomelli et al. (2016) coefficients; the real confocal volumes are rendered
  through the same transform (DNA -> hematoxylin, membrane/cytoplasm -> eosin).
