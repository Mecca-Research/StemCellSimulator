// Real 3D stem-cell colony -> segmented cells -> growing 3D digital twin.
import * as THREE from 'three';
import { cellMaterial, lineMaterial } from '../engine/materials.js';
import { VolumeView, mosaicToVolume, mosaicToLabels } from '../engine/volume.js';
import { decodeMesh, unitSphere } from '../engine/geometry.js';
import { InstancePool, idColor } from '../engine/instances.js';
import { Tissue, Cell, PHASE_FRACTIONS, mitosisStage, doublingTime } from '../models/morpho/agents3d.js';
import { collier } from '../models/pathways/notch.js';
import { RNG } from '../models/core/rng.js';
import { fmt } from '../engine/chart.js';

const HOURS_PER_SECOND = 1.2;

export default {
  about: `
    <p>A real <b>3D confocal stack</b> of a cell monolayer (Allen Institute for Cell Science,
    CC0; 60 optical sections, membrane dye in red, DNA in green) is ray-marched directly.
    Every cell and nucleus was segmented in 3D by <code>tools/build_assets.py</code> and
    converted into a surface mesh. One cell is caught in <b>anaphase</b>.</p>
    <p>The <b>digital twin</b> starts from those exact cells (positions, volumes,
    nucleus shapes, contacts) and continues them forward in time with a 3D mechanical
    model: G1-S-G2-M cycle, mitotic rounding, oriented division, adhesion, contact
    inhibition and Delta-Notch lateral inhibition across real cell contacts.</p>`,
  paperRef: 'Paper: "Stage 1: Basic Stem Cell Reactions" → "Stem Cell Self-Renewal: HSC → HSC + HSC", and "Higher-Level Functional Interactions: cell-cell communication".',

  async create(ctx) {
    const { stage, root, panel, assets, legend, setStatus } = ctx;
    const [meta, cal, volImg, labImg, meshBuf] = await Promise.all([
      assets.loadJSON('hipsc_cells.json'),
      assets.loadJSON('calibration.json'),
      assets.loadImageData('hipsc_volume.jpg'),
      assets.loadImageData('hipsc_labels.png'),
      assets.loadBinary('hipsc_meshes.bin'),
    ]);
    const S = meta.summary;
    const [ex, ey, ez] = meta.extent_xyz_um;
    const layerBottom = Math.min(...meta.cells.map((c) => c.cell_centroid_um[2])) - S.cell_height_um_median / 2;

    // ------------------------------------------------------------ real data
    // image space (x = column, y = row, z = optical section) -> world (x, up, -y)
    const imgGroup = new THREE.Group();
    imgGroup.rotation.x = -Math.PI / 2;
    imgGroup.position.set(-ex / 2, 0, ey / 2);
    root.add(imgGroup);

    const volTex = mosaicToVolume(volImg, meta.mosaic, 2);
    const labels = mosaicToLabels(labImg, meta.mosaic);
    const volume = new VolumeView(volTex, [ex, ey, ez], { labels: labels.tex, c0: 0xff2a55, c1: 0x2cff7a, style: 0, steps: 200 });
    imgGroup.add(volume);

    const cellMeshes = [], nucMeshes = [];
    const realNucMat = cellMaterial({ role: 'nucleus', color: 0x8fb4ff });
    const mitoticNucMat = cellMaterial({ role: 'nucleus', color: 0xffe07a, emissive: 0.25 });
    const realMeshGroup = new THREE.Group();
    imgGroup.add(realMeshGroup);
    for (const c of meta.cells) {
      const cg = decodeMesh(meshBuf, meta, c.mesh.cell);
      const ng = decodeMesh(meshBuf, meta, c.mesh.nucleus);
      const cm = new THREE.Mesh(cg, cellMaterial({ role: 'membrane', color: idColor(c.id), opacity: 0.16 }));
      const nm = new THREE.Mesh(ng, c.mitotic ? mitoticNucMat : realNucMat);
      cm.userData.cell = nm.userData.cell = c;
      cm.renderOrder = 2;
      realMeshGroup.add(cm, nm);
      cellMeshes.push(cm);
      nucMeshes.push(nm);
    }

    // substrate (glass) grid under the monolayer
    const grid = new THREE.GridHelper(ex * 3, 30, 0x2a3a50, 0x16202c);
    grid.material = lineMaterial({ color: 0x2a3a50, opacity: 0.55 });
    grid.position.y = layerBottom;
    root.add(grid);

    // ------------------------------------------------------------ digital twin
    const twinGroup = new THREE.Group();
    root.add(twinGroup);
    const sphere = unitSphere(3);
    const bodyPool = new InstancePool(sphere, cellMaterial({ role: 'membrane', opacity: 0.2 }), 512, twinGroup);
    const nucPool = new InstancePool(sphere, cellMaterial({ role: 'nucleus', color: 0xa9c4ff, stain: 0x39ff88 }), 512, twinGroup);
    const repPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'reporter' }), 256, twinGroup);
    bodyPool.mesh.renderOrder = 3;

    const cellHeight = S.cell_height_um_median;

    const semi = S.nucleus_semi_axes_um_median;

    const params = {
      view: 'both', colorBy: 'founder', tM: 1.0, cycle: 1.0 / Math.max(cal.mitosis?.mitotic_index ?? 0.062, 1e-3),
      contact: true, noise: 0.08, adhesion: 1.4, maxCells: 700, notch: true, showVolume: true, showMeshes: true,
      volStyle: 0, winLo: 0.08, winHi: 0.85, labelTint: 0.0, growth: true,
    };

    let tissue, rng, history;
    // Mechanical radius from the measured in-plane spacing of real nuclei; each
    // cell is scaled by its own measured volume (cells cut by the field of view
    // get the median).
    const Rmech = S.nn_distance_um_median / 2;
    const cellR = (c) => (c.touches_border ? Rmech
      : Rmech * THREE.MathUtils.clamp(Math.sqrt(c.cell_volume_um3 / S.cell_volume_um3_median), 0.8, 1.25));
    const toWorld = ([x, y, z]) => [x - ex / 2, z, -(y - ey / 2)];

    function buildTwin() {
      rng = new RNG(20241018);
      tissue = new Tissue({
        dims: 2, noise: params.noise, muRep: 10, adhesion: () => params.adhesion, cut: 1.3,
        substrate: { y: layerBottom + cellHeight / 2, k: 6 }, rng,
        contactInhibition: params.contact ? 6 : null,
        onDivide: (m, a, b) => {
          for (const d of [a, b]) {
            d.state.N = (m.state.N ?? 0.5) + 0.02 * rng.normal();
            d.state.D = (m.state.D ?? 0.5) + 0.02 * rng.normal();
          }
        },
      });
      for (const c of meta.cells) {
        const p = toWorld(c.cell_centroid_um);
        p[1] = layerBottom + cellHeight / 2;
        // border cells are truncated by the field of view: use the median size
        const R = cellR(c);
        const cell = new Cell({
          p, R, cycle: params.cycle * (0.85 + 0.3 * rng.next()),
          state: { founder: c.id, N: 0.5 + 0.05 * rng.normal(), D: 0.5 + 0.05 * rng.normal(), real: true },
        });
        // desynchronise: place each real cell at a random point of its cycle and
        // back-compute its birth radius so the measured size is kept (R = R0 sqrt(1 + f))
        const fM = PHASE_FRACTIONS.M;
        const f = rng.next() * (1 - fM);
        let acc = 0;
        for (const ph of ['G1', 'S', 'G2']) {
          const L = PHASE_FRACTIONS[ph];
          if (f < acc + L) { cell.phase = ph; cell.phaseT = (f - acc) * cell.cycle; break; }
          acc += L;
        }
        cell.R0 = R / Math.sqrt(1 + f / (1 - fM));
        if (c.mitotic) {
          // the real anaphase cell: continue its division along the measured long axis
          // of its chromatin (eigenvector of smallest inertia = longest axis, zyx order)
          cell.R0 = R / Math.SQRT2;
          const ax = c.nucleus.axes_zyx[0];
          const w = toWorld([ax[2], ax[1], ax[0]]).map((v, k) => v - toWorld([0, 0, 0])[k]);
          w[1] = 0;
          const n = Math.hypot(...w) || 1;
          cell.phase = 'M';
          cell.phaseT = 0.6 * cell.cycle * PHASE_FRACTIONS.M;
          cell.axis = w.map((v) => v / n);
        }
        tissue.add(cell);
      }
      history = { t: [], n: [] };
      countChart?.clear();
      miChart?.clear();
    }

    // ------------------------------------------------------------ controls
    panel.section('View');
    panel.select({
      label: 'Show', value: params.view,
      options: [
        { value: 'both', label: 'Real volume + reconstruction + twin (side by side)' },
        { value: 'real', label: 'Real data only' },
        { value: 'twin', label: 'Digital twin only' },
      ],
      onChange: (v) => { params.view = v; layout(); },
    });
    panel.toggle({ label: 'Ray-marched confocal volume', value: params.showVolume, onChange: (v) => { params.showVolume = v; layout(); } });
    panel.toggle({ label: 'Segmented 3D cells & nuclei', value: params.showMeshes, onChange: (v) => { params.showMeshes = v; layout(); } });
    panel.select({
      label: 'Volume rendering', value: '0',
      options: [{ value: '0', label: 'Maximum-intensity projection' }, { value: '1', label: 'Emission-absorption compositing' }],
      onChange: (v) => { volume.material.uniforms.uStyle.value = +v; },
    });
    panel.slider({ label: 'Display window (low)', min: 0, max: 0.6, value: params.winLo, onChange: (v) => { params.winLo = v; volume.setWindow(params.winLo, params.winHi); } });
    panel.slider({ label: 'Display window (high)', min: 0.3, max: 1, value: params.winHi, onChange: (v) => { params.winHi = v; volume.setWindow(params.winLo, params.winHi); } });
    panel.slider({ label: 'Tint voxels by segmented cell', min: 0, max: 1, value: 0, onChange: (v) => { volume.material.uniforms.uLabelTint.value = v; } });
    volume.setWindow(params.winLo, params.winHi);

    panel.section('Digital twin');
    panel.select({
      label: 'Colour cells by', value: params.colorBy,
      options: [
        { value: 'founder', label: 'Founder (clonal lineage of each real cell)' },
        { value: 'phase', label: 'Cell-cycle phase' },
        { value: 'notch', label: 'Notch / Delta (lateral inhibition)' },
        { value: 'generation', label: 'Generation' },
      ],
      onChange: (v) => { params.colorBy = v; updateLegend(); },
    });
    const cycleReadout = panel.readouts(['Cell-cycle length T']);
    const tMs = panel.slider({
      label: 'Assumed M-phase duration t_M (h)', min: 0.5, max: 2, value: params.tM,
      onChange: (v) => {
        params.tM = v;
        params.cycle = v / Math.max(cal.mitosis?.mitotic_index ?? 0.062, 1e-3);
        for (const c of tissue.cells) c.cycle = params.cycle * (0.85 + 0.3 * rng.next());
        cycleReadout.set('Cell-cycle length T', `${params.cycle.toFixed(1)} h`);
      },
    });
    cycleReadout.set('Cell-cycle length T', `${params.cycle.toFixed(1)} h`);
    panel.note('T = t<sub>M</sub> / MI, with the mitotic index MI measured on the real CellProfiler mitosis field (see gallery).');
    panel.toggle({ label: 'Contact inhibition (G1 arrest at ≥ 6 neighbours)', value: params.contact, onChange: (v) => { params.contact = v; tissue.contactInhibition = v ? 6 : null; } });
    panel.toggle({ label: 'Delta–Notch lateral inhibition', value: params.notch, onChange: (v) => { params.notch = v; } });
    panel.slider({ label: 'Cell–cell adhesion μ_adh', min: 0, max: 4, value: params.adhesion, onChange: (v) => { params.adhesion = v; } });
    panel.slider({ label: 'Motility noise', min: 0, max: 1, value: params.noise, onChange: (v) => { params.noise = v; tissue.noise = v; } });
    panel.buttons([
      { label: 'Restart from real colony', primary: true, onClick: () => { buildTwin(); stage.simTime = 0; } },
    ]);

    panel.section('Measured on the real stack');
    panel.table(['quantity', 'value'], [
      ['cells segmented (3D)', S.n_cells],
      ['median nucleus volume (µm³)', S.nucleus_volume_um3_median],
      ['median cell volume (µm³)', S.cell_volume_um3_median],
      ['nucleus : cytoplasm volume', S.nc_ratio_median],
      ['monolayer height (µm)', S.cell_height_um_median],
      ['nucleus semi-axes (µm)', S.nucleus_semi_axes_um_median.map((v) => v.toFixed(1)).join(' × ')],
      ['3D contact neighbours', S.neighbours_mean_interior],
      ['mitotic index (2D field)', cal.mitosis?.mitotic_index ?? NaN],
    ]);

    panel.section('Growth');
    const countChart = panel.chart({
      title: 'Cell number (log) vs expected 2^(t/T)', logY: true, xLabel: 'hours',
      series: [{ name: 'twin', color: '#5ee0c1' }, { name: 'ideal', color: '#8a98ab', dash: [4, 3] }],
    });
    const miChart = panel.chart({
      title: 'Mitotic index: twin vs measured', xLabel: 'hours', yRange: [0, 0.2],
      series: [{ name: 'twin MI', color: '#ffb454' }, { name: 'measured', color: '#ff6fae', dash: [4, 3] }],
    });
    const live = panel.readouts(['cells', 'doubling time (fit)', 'Delta-high fraction', 'sender–sender contacts']);

    // ------------------------------------------------------------ rendering
    const phaseColor = { G1: new THREE.Color(0x4aa3ff), S: new THREE.Color(0x5ee0c1), G2: new THREE.Color(0xc38bff), M: new THREE.Color(0xffb454) };
    const tmpC = new THREE.Color();
    const white = new THREE.Color(0xffffff);
    let twinOffset = 0;

    function layout() {
      const v = params.view;
      const showReal = v !== 'twin';
      volume.visible = showReal && params.showVolume;
      realMeshGroup.visible = showReal && params.showMeshes;
      twinGroup.visible = v !== 'real';
      twinOffset = v === 'both' ? ex * 1.35 : 0;
      twinGroup.position.x = twinOffset;
      const cx = v === 'both' ? twinOffset / 2 : 0;
      stage.frame([cx, layerBottom + 6, 0], v === 'both' ? ex * 1.3 : ex * 0.8, new THREE.Vector3(0.1, 0.9, 0.75));
      applyMode(stage.mode);
      stage.clipAxis.set(0, 1, 0);
      stage.clipRange = [-cellHeight * 0.9, cellHeight * 0.9];
    }

    // confocal mode shows the real stack the way the microscope saw it (no mesh overlay)
    function applyMode(mode) {
      if (params.view !== 'twin') realMeshGroup.visible = params.showMeshes && !(mode === 'confocal' && params.showVolume);
      grid.visible = mode !== 'confocal';
    }

    function updateLegend() {
      const items = {
        founder: [{ color: '#7fd1ff', label: 'each real cell and its progeny keep one hue' }],
        phase: Object.entries(phaseColor).map(([k, c]) => ({ color: `#${c.getHexString()}`, label: k })),
        notch: [{ color: '#39ff88', label: 'Notch-high (receiver, progenitor)' }, { color: '#ff3b9d', label: 'Delta-high (sender)' }],
        generation: [{ color: '#4aa3ff', label: 'generation 0 (real)' }, { color: '#ffb454', label: 'later generations' }],
      }[params.colorBy];
      legend([
        { color: '#ff2a55', label: 'membrane dye (real volume)' },
        { color: '#2cff7a', label: 'DNA dye (real volume)' },
        ...items,
      ]);
    }

    function cellColor(c, out) {
      switch (params.colorBy) {
        case 'phase': return out.copy(phaseColor[c.phase]);
        case 'notch': {
          const d = THREE.MathUtils.clamp(c.state.D ?? 0.5, 0, 1);
          return out.setRGB(0.22 + 0.78 * d, 1 - 0.75 * d, 0.53 + 0.1 * d);
        }
        case 'generation': return out.setHSL(0.6 - Math.min(c.generation, 6) * 0.085, 0.75, 0.58);
        default: return idColor(c.state.founder ?? c.id, 0.62, 0.6, out);
      }
    }

    function drawTwin() {
      bodyPool.begin(); nucPool.begin(); repPool.begin();
      for (const c of tissue.cells) {
        const col = cellColor(c, tmpC);
        // confluent look: the drawn footprint slightly exceeds the mechanical radius
        const lat = c.R * 1.18, h = cellHeight / 2;
        const nucScale = 0.82 * Math.min(c.R / Rmech, 1.3);
        const m = c.mitosis;
        if (c.phase !== 'M') {
          bodyPool.put(c.p[0], c.p[1], c.p[2], lat, h, lat, col);
          // interphase nucleus: measured median semi-axes (long, short in-plane; vertical)
          nucPool.putOriented(c.p, [Math.cos(c.id), 0, Math.sin(c.id)],
            semi[0] * nucScale, semi[2], semi[1] * nucScale, white);
        } else {
          // mitotic rounding: flattened cell becomes a sphere of equal volume
          const vol = (4 / 3) * Math.PI * lat * lat * h;
          const rs = Math.cbrt((3 * vol) / (4 * Math.PI));
          const round = Math.min(m / 0.25, 1);
          const cy = c.p[1] + (rs - h) * round * 0.6;
          const stageName = mitosisStage(m);
          const sep = m < 0.55 ? 0 : Math.min((m - 0.55) / 0.35, 1); // anaphase -> cytokinesis
          if (sep < 0.35) {
            const sx = THREE.MathUtils.lerp(lat, rs, round) * (1 + 0.35 * sep);
            const sy = THREE.MathUtils.lerp(h, rs, round);
            bodyPool.putOriented([c.p[0], cy, c.p[2]], c.axis, sx, sy, THREE.MathUtils.lerp(lat, rs, round) * (1 - 0.15 * sep), col);
          } else {
            // cleavage furrow: two daughter lobes
            const r2 = rs / Math.cbrt(2) * 1.04;
            const d = rs * 0.62 * sep;
            for (const sg of [1, -1]) {
              bodyPool.putOriented([c.p[0] + sg * d * c.axis[0], cy, c.p[2] + sg * d * c.axis[2]], c.axis, r2, r2, r2, col);
            }
          }
          // chromatin: condensing nucleus -> metaphase plate -> two sister sets
          const bright = tmpC.setRGB(1, 0.95, 0.55);
          if (stageName === 'prophase' || stageName === 'prometaphase') {
            const k = 1 - 0.35 * Math.min(m / 0.32, 1);
            nucPool.put(c.p[0], cy, c.p[2], 4.2 * k, 3.6 * k, 4.2 * k, bright);
          } else if (stageName === 'metaphase') {
            nucPool.putOriented([c.p[0], cy, c.p[2]], c.axis, 0.9, 3.4, 3.4, bright);
          } else {
            const d = Math.min(1, (m - 0.55) / 0.3) * rs * 0.62;
            for (const sg of [1, -1]) {
              const telo = stageName === 'telophase' || stageName === 'cytokinesis';
              nucPool.putOriented([c.p[0] + sg * d * c.axis[0], cy, c.p[2] + sg * d * c.axis[2]], c.axis,
                telo ? 2.4 : 1.0, telo ? 2.6 : 3.0, telo ? 2.6 : 3.0, bright);
            }
          }
        }
        if (params.colorBy === 'notch' && (c.state.D ?? 0) > 0.6) {
          repPool.put(c.p[0], c.p[1] + cellHeight * 0.55, c.p[2], 1.2, 1.2, 1.2, tmpC.set(0xff3b9d));
        }
      }
      bodyPool.end(); nucPool.end(); repPool.end();
    }

    // Notch lateral inhibition over the twin's live contact graph
    function notchStep(dtH) {
      const cells = tissue.cells;
      const index = new Map(cells.map((c, i) => [c.id, i]));
      const s = { N: new Float64Array(cells.length), D: new Float64Array(cells.length) };
      cells.forEach((c, i) => { s.N[i] = c.state.N ?? 0.5; s.D[i] = c.state.D ?? 0.5; });
      const nbrs = cells.map((c) => (tissue.neighbours.get(c.id) ?? []).map((id) => index.get(id)).filter((j) => j !== undefined));
      const sub = Math.max(1, Math.ceil(dtH / 0.1));
      for (let k = 0; k < sub; k++) collier.step(s, nbrs, (dtH / sub) * 3, collier.defaults);
      cells.forEach((c, i) => { c.state.N = s.N[i]; c.state.D = s.D[i]; });
      return { s, nbrs };
    }

    // ------------------------------------------------------------ picking
    stage.setPickables([
      ...[...cellMeshes, ...nucMeshes].map((object) => ({
        object,
        info: () => {
          const c = object.userData.cell;
          if (!realMeshGroup.visible) return null;
          return `Real cell #${c.id}${c.mitotic ? ' · MITOTIC (anaphase)' : ''}\n` +
            `nucleus ${fmt(c.nucleus_volume_um3)} µm³ · cell ${fmt(c.cell_volume_um3)} µm³\n` +
            `N:C ${fmt(c.nc_ratio)} · height ${fmt(c.cell_height_um)} µm · ${c.neighbours.length} contacts` +
            (c.touches_border ? '\n(cut by field of view)' : '');
        },
      })),
      {
        get object() { return bodyPool.mesh; }, // the pool swaps meshes when it grows
        info: (hit) => {
          if (!twinGroup.visible) return null;
          const c = tissue.cells[bodyIndexToCell(hit.instanceId)];
          if (!c) return null;
          return `Twin cell ${c.id} · founder #${c.state.founder}\n${c.phase}${c.phase === 'M' ? ` (${mitosisStage(c.mitosis)})` : ''} · generation ${c.generation}\n` +
            `Notch ${fmt(c.state.N ?? 0)} · Delta ${fmt(c.state.D ?? 0)}`;
        },
      },
    ]);
    // body instances are 1 per cell except for late cytokinesis (2); map back
    function bodyIndexToCell(inst) {
      let k = 0;
      for (let i = 0; i < tissue.cells.length; i++) {
        const c = tissue.cells[i];
        const n = c.phase === 'M' && c.mitosis >= 0.55 + 0.35 * 0.35 ? 2 : 1;
        if (inst < k + n) return i;
        k += n;
      }
      return -1;
    }

    buildTwin();
    layout();
    updateLegend();
    drawTwin();

    let lastSample = -1;
    return {
      update(dt) {
        const dtH = dt * HOURS_PER_SECOND;
        if (dtH > 0 && twinGroup.visible) {
          const sub = Math.max(1, Math.ceil(dtH / 0.05));
          for (let k = 0; k < sub; k++) {
            tissue.growth = tissue.cells.length < params.maxCells;
            tissue.step(dtH / sub);
          }
          let notchStats = null;
          if (params.notch) notchStats = notchStep(dtH);
          const t = tissue.t;
          if (t - lastSample > 0.25 || lastSample < 0) {
            lastSample = t;
            history.t.push(t); history.n.push(tissue.cells.length);
            const n0 = meta.cells.length;
            countChart.push(t, [tissue.cells.length, n0 * 2 ** (t / params.cycle)]);
            const mi = tissue.count((c) => c.phase === 'M') / tissue.cells.length;
            miChart.push(t, [mi, cal.mitosis?.mitotic_index ?? NaN]);
            live.set('cells', tissue.cells.length);
            const recent = history.t.length > 8 ? -Math.min(history.t.length, 60) : 0;
            live.set('doubling time (fit)', `${fmt(doublingTime(history.t.slice(recent), history.n.slice(recent)))} h`);
            if (notchStats) {
              live.set('Delta-high fraction', collier.senderFraction(notchStats.s, 0.6));
              live.set('sender–sender contacts', collier.senderAdjacency(notchStats.s, notchStats.nbrs, 0.6));
            }
          }
        }
        drawTwin();
        setStatus(`t = ${tissue.t.toFixed(1)} h · ${tissue.cells.length} cells`);
      },
      reset() { buildTwin(); lastSample = -1; },
      onMode(mode) { applyMode(mode); },
      dispose() { volTex.dispose(); labels.tex.dispose(); },
    };
  },
};
