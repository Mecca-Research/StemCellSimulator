// Organ level: real 3D confocal kidney tissue next to a 3D nephron model.
import * as THREE from 'three';
import { cellMaterial } from '../engine/materials.js';
import { VolumeView, mosaicToVolume } from '../engine/volume.js';
import { unitSphere } from '../engine/geometry.js';
import { InstancePool } from '../engine/instances.js';
import { RNG } from '../models/core/rng.js';
import { NEPHRON_SEGMENTS, SEGMENT_BOUNDS, filtrateProfile } from '../models/organ/nephron.js';

export default {
  about: `
    <p>Left: a real <b>3D confocal volume of mouse kidney</b> (16 optical sections,
    three emission channels; G. Buckley, CC0) rendered by ray marching. Tubule
    cross-sections and two <b>glomeruli</b> (renal corpuscles) are visible.</p>
    <p>Right: a 3D <b>nephron model</b> built from the paper's organ map — renal
    corpuscle (glomerular capillary tuft, podocytes, mesangial cells, Bowman's
    capsule) → proximal tubule → loop of Henle → distal tubule → collecting duct —
    with epithelial cells arranged around each lumen and filtrate flowing through it.
    Segment reabsorption fractions are textbook physiology (≈65 % proximal, ≈25 %
    loop, ≈5 % distal, collecting duct regulated by ADH); the geometry is schematic.</p>`,
  paperRef: 'Paper: "Organs → Kidney: Nephrons (Podocytes, Proximal Tubule Cells, Distal Tubule Cells, Collecting Duct Cells), Renal Corpuscle (Glomerular Endothelial Cells, Mesangial Cells)" and "Urinary and Cardiovascular Systems".',

  async create(ctx) {
    const { stage, root, panel, assets, legend, setStatus } = ctx;
    const [meta, img] = await Promise.all([assets.loadJSON('kidney_volume.json'), assets.loadImageData('kidney_volume.jpg')]);
    const [ex, ey, ez] = meta.extent_xyz_um;

    // --------------------------------------------------------- real volume
    const volGroup = new THREE.Group();
    volGroup.rotation.x = -Math.PI / 2;
    volGroup.position.set(-ex - 40, 0, ey / 2);
    const volCenter = new THREE.Vector3(-ex / 2 - 40, 0, 0);
    root.add(volGroup);
    const tex = mosaicToVolume(img, meta.mosaic, 3);
    // true-ish emission colours: 605 nm red, 515 nm green, 450 nm blue
    const volume = new VolumeView(tex, [ex, ey, ez * 2.2], { c0: 0xff3030, c1: 0x30ff60, c2: 0x3a6bff, style: 0, steps: 160 });
    volume.setWindow(0.06, 0.9);
    volGroup.add(volume);
    panel.section('Real tissue volume');
    panel.select({
      label: 'Rendering', value: '0',
      options: [{ value: '0', label: 'Maximum-intensity projection' }, { value: '1', label: 'Emission-absorption compositing' }],
      onChange: (v) => { volume.material.uniforms.uStyle.value = +v; },
    });
    panel.slider({ label: 'Density (compositing)', min: 0.2, max: 4, value: 1, onChange: (v) => { volume.material.uniforms.uDensity.value = v; } });
    panel.slider({ label: 'Window low', min: 0, max: 0.5, value: 0.06, onChange: (v) => { const u = volume.material.uniforms; u.uWinLo.value.setScalar(v); } });
    panel.note(`Voxel ${meta.spacing_zyx_um[2]} × ${meta.spacing_zyx_um[1]} × ${meta.spacing_zyx_um[0]} µm; field ${ex.toFixed(0)} × ${ey.toFixed(0)} µm (z drawn 2.2× exaggerated so the 20 µm slab reads as 3D). The nephron model on the right is drawn at the same scale.`);

    // --------------------------------------------------------- nephron model
    // model units are scaled x2.4 so the glomerulus (~120 um) and tubules
    // (~40-60 um outer diameter) match the real tissue next to it
    const neph = new THREE.Group();
    neph.scale.setScalar(2.4);
    neph.position.set(40, 60, 0);
    root.add(neph);
    const rng = new RNG(7);

    // centre-line through the segments (schematic, compressed lengths)
    const pts = [
      [0, 40, 0], // urinary pole of Bowman's capsule
      [25, 48, 20], [45, 36, -10], [30, 22, -32], [55, 14, -20], [70, 26, 10], [92, 30, -8], // proximal convoluted
      [104, 10, 0], [106, -60, 4], [108, -125, 0], // descending limb
      [118, -140, 4], [130, -125, 0], // hairpin
      [132, -60, -4], [134, 5, -2], // thick ascending limb
      [120, 22, 24], [100, 32, 40], [80, 20, 52], [96, 6, 62], // distal convoluted (returns to the corpuscle: macula densa)
      [130, 0, 70], [160, -30, 72], [175, -130, 70], // collecting duct
    ];
    const curve = new THREE.CatmullRomCurve3(pts.map((p) => new THREE.Vector3(...p)), false, 'centripetal');
    const segs = NEPHRON_SEGMENTS;
    const bounds = SEGMENT_BOUNDS; // arc-length fractions per segment

    const lumenMat = cellMaterial({ role: 'cytoplasm', color: 0xbfd7ff, opacity: 0.12 });
    const lumen = new THREE.Mesh(new THREE.TubeGeometry(curve, 600, 7.5, 16, false), lumenMat);
    lumen.renderOrder = 4;
    neph.add(lumen);

    // ~3,700 epithelial cells: the nuclei sit inside translucent membranes, so
    // a coarser sphere (80 instead of 180 triangles) costs nothing visible
    const epi = new InstancePool(unitSphere(2), cellMaterial({ role: 'membrane', opacity: 0.55, noise: 0.35 }), 4096, neph);
    const epiNuc = new InstancePool(unitSphere(1), cellMaterial({ role: 'nucleus', color: 0xa9c4ff }), 4096, neph);
    const tmp = new THREE.Color();
    const up = new THREE.Vector3();
    const cells = [];
    const rings = 360;
    for (let i = 0; i < rings; i++) {
      const u = (i + 0.5) / rings;
      const si = bounds.findIndex((b, k) => u >= b && u < bounds[k + 1]);
      const seg = segs[Math.max(si, 0)];
      const p = curve.getPointAt(u), t = curve.getTangentAt(u);
      up.set(0, 1, 0); if (Math.abs(t.y) > 0.9) up.set(1, 0, 0);
      const n1 = new THREE.Vector3().crossVectors(t, up).normalize();
      const n2 = new THREE.Vector3().crossVectors(t, n1).normalize();
      const per = seg.cellsPerRing;
      for (let k = 0; k < per; k++) {
        const a = (k / per) * Math.PI * 2 + (i % 2) * (Math.PI / per);
        const dir = n1.clone().multiplyScalar(Math.cos(a)).add(n2.clone().multiplyScalar(Math.sin(a)));
        const r = seg.radius;
        cells.push({ p: p.clone().add(dir.clone().multiplyScalar(r)), dir, t: t.clone(), seg, u });
      }
    }
    const _m = new THREE.Matrix4(), _q = new THREE.Quaternion(), _s = new THREE.Vector3();
    function drawEpithelium() {
      epi.begin(); epiNuc.begin();
      for (const c of cells) {
        // cuboidal/columnar cell: radial height, tangential width
        _q.setFromUnitVectors(new THREE.Vector3(0, 1, 0), c.dir);
        _s.set(c.seg.width, c.seg.height, c.seg.width * 0.9);
        _m.compose(c.p, _q, _s);
        epi.putMatrix(_m, tmp.set(c.seg.color));
        _s.set(c.seg.width * 0.42, c.seg.height * 0.36, c.seg.width * 0.38);
        const np = c.p.clone().add(c.dir.clone().multiplyScalar(c.seg.height * (c.seg.basalNucleus ? 0.35 : 0.05)));
        _m.compose(np, _q, _s);
        epiNuc.putMatrix(_m);
      }
      epi.end(); epiNuc.end();
    }
    drawEpithelium();

    // renal corpuscle: capillary tuft + podocytes + Bowman's capsule
    const corpuscle = new THREE.Group();
    corpuscle.position.set(-6, 52, -6);
    neph.add(corpuscle);
    const capMat = cellMaterial({ role: 'solid', color: 0xd4495b, noise: 0.25, stain: 0xff5570 });
    for (let k = 0; k < 14; k++) {
      const loop = [];
      const c0 = new THREE.Vector3(rng.uniform(-8, 8), rng.uniform(-8, 8), rng.uniform(-8, 8));
      for (let j = 0; j < 9; j++) {
        const a = (j / 8) * Math.PI * 2;
        const v = new THREE.Vector3(Math.cos(a) * 9, Math.sin(a * 2) * 4, Math.sin(a) * 9).applyAxisAngle(new THREE.Vector3(rng.next(), rng.next(), rng.next()).normalize(), rng.uniform(0, Math.PI));
        loop.push(v.add(c0).clampLength(0, 17));
      }
      corpuscle.add(new THREE.Mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3(loop, true), 60, 1.6, 8, true), capMat));
    }
    const podo = new InstancePool(unitSphere(2), cellMaterial({ role: 'membrane', color: 0x9dffcf, opacity: 0.6 }), 256, corpuscle);
    const podoNuc = new InstancePool(unitSphere(1), cellMaterial({ role: 'nucleus', color: 0xa9c4ff }), 256, corpuscle);
    const feet = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter', color: 0x6fffc0 }), 2048, corpuscle);
    podo.begin(); podoNuc.begin(); feet.begin();
    const dir = [0, 0, 0];
    for (let k = 0; k < 60; k++) {
      rng.onSphere(dir);
      const r = 19;
      podo.put(dir[0] * r, dir[1] * r, dir[2] * r, 3.2, 3.2, 3.2, tmp.set(0x9dffcf));
      podoNuc.put(dir[0] * (r + 0.6), dir[1] * (r + 0.6), dir[2] * (r + 0.6), 1.8, 1.8, 1.8);
      for (let f = 0; f < 10; f++) { // interdigitating foot processes on the capillaries
        const d2 = [dir[0] + rng.normal() * 0.18, dir[1] + rng.normal() * 0.18, dir[2] + rng.normal() * 0.18];
        const n = Math.hypot(...d2);
        feet.put((d2[0] / n) * 17.2, (d2[1] / n) * 17.2, (d2[2] / n) * 17.2, 0.7, 0.7, 0.7);
      }
    }
    // mesangial cells in the tuft core
    for (let k = 0; k < 12; k++) {
      rng.onSphere(dir);
      podo.put(dir[0] * 6, dir[1] * 6, dir[2] * 6, 2.4, 2.4, 2.4, tmp.set(0xffc36e));
      podoNuc.put(dir[0] * 6, dir[1] * 6, dir[2] * 6, 1.4, 1.4, 1.4);
    }
    podo.end(); podoNuc.end(); feet.end();
    const capsule = new THREE.Mesh(new THREE.SphereGeometry(25, 48, 32, 0, Math.PI * 2, 0.35, Math.PI - 0.35),
      cellMaterial({ role: 'membrane', color: 0xcfd8ff, opacity: 0.1 }));
    capsule.rotation.x = Math.PI * 0.8;
    capsule.renderOrder = 5;
    corpuscle.add(capsule);

    // afferent / efferent arterioles
    const artMat = cellMaterial({ role: 'solid', color: 0xb8323f, stain: 0xff4466 });
    for (const [a, b] of [[[-60, 90, -30], [-10, 58, -8]], [[-8, 58, -2], [-55, 95, 20]]]) {
      const c = new THREE.CatmullRomCurve3([new THREE.Vector3(...a), new THREE.Vector3((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + 8, (a[2] + b[2]) / 2), new THREE.Vector3(...b)]);
      corpuscle.add(new THREE.Mesh(new THREE.TubeGeometry(c, 40, 3.2, 12), artMat).translateX(6).translateY(-52).translateZ(6));
    }

    // filtrate particles flowing along the tubule: survival ~ remaining fraction
    const flow = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter' }), 1024, neph);
    const params = { gfr: 1, adh: 0.5, showVolume: true, showModel: true, speed: 1 };
    const particles = [];
    let emitAcc = 0;
    let profile = filtrateProfile(params.adh);
    const remainChart = (() => {
      panel.section('Nephron model');
      panel.slider({ label: 'Filtration rate (relative GFR)', min: 0.2, max: 2, value: params.gfr, onChange: (v) => { params.gfr = v; } });
      panel.slider({ label: 'ADH (collecting-duct water permeability)', min: 0, max: 1, value: params.adh, onChange: (v) => { params.adh = v; profile = filtrateProfile(v); drawProfile(); } });
      panel.toggle({ label: 'Show real tissue volume', value: true, onChange: (v) => { volGroup.visible = v; } });
      panel.toggle({ label: 'Show nephron model', value: true, onChange: (v) => { neph.visible = v; } });
      return panel.chart({ title: 'Filtrate remaining along the nephron (fraction of filtered volume)', xLabel: 'position along nephron', yRange: [0, 1.05], series: [{ name: 'volume remaining', color: '#5ee0c1' }] });
    })();
    function drawProfile() {
      const xs = [], ys = [];
      for (let i = 0; i <= 100; i++) { xs.push(i / 100); ys.push(profile(i / 100)); }
      remainChart.set(xs, [ys]);
    }
    drawProfile();
    panel.table(['segment', 'cells (paper map)', 'reabsorbs'], segs.map((s) => [s.name, s.cells, s.reabsorbs]));

    legend([
      { color: '#ff3030', label: 'real volume: 605 nm channel' },
      { color: '#30ff60', label: 'real volume: 515 nm channel' },
      { color: '#3a6bff', label: 'real volume: 450 nm channel' },
      ...segs.map((s) => ({ color: `#${new THREE.Color(s.color).getHexString()}`, label: s.name })),
      { color: '#9dffcf', label: 'podocytes (foot processes)' },
      { color: '#ffc36e', label: 'mesangial cells' },
    ]);

    stage.frame([(volCenter.x + 460) / 2, -30, 40], 470, new THREE.Vector3(0.0, 0.8, 0.8));
    stage.clipAxis.set(0, 0, 1);

    stage.setPickables([
      {
        get object() { return epi.mesh; },
        info: (hit) => { const c = cells[hit.instanceId]; return c ? `${c.seg.name}\n${c.seg.cells}\nreabsorbs: ${c.seg.reabsorbs}` : null; },
      },
      { get object() { return podo.mesh; }, info: (hit) => (hit.instanceId < 60 ? 'Podocyte: foot processes wrap the glomerular capillaries (filtration slit diaphragm)' : 'Mesangial cell: structural support of the capillary tuft') },
    ]);

    let simT = 0;
    return {
      update(dt) {
        simT += dt;
        emitAcc += dt * 40 * params.gfr;
        while (emitAcc > 1 && particles.length < 900) { emitAcc -= 1; particles.push({ u: 0, off: [rng.normal() * 2.5, rng.normal() * 2.5], seed: rng.next() }); }
        flow.begin();
        for (let i = particles.length - 1; i >= 0; i--) {
          const q = particles[i];
          q.u += dt * 0.06 * params.gfr;
          // a particle is reabsorbed once its seed exceeds the remaining fraction
          if (q.u >= 1 || q.seed > profile(q.u)) { particles.splice(i, 1); continue; }
          const p = curve.getPointAt(Math.min(q.u, 1));
          flow.put(p.x + q.off[0], p.y + q.off[1], p.z, 0.9, 0.9, 0.9, tmp.setHSL(0.52 - 0.1 * q.u, 0.9, 0.6));
        }
        flow.end();
        corpuscle.rotation.y = Math.sin(simT * 0.2) * 0.05;
        setStatus(`filtrate particles in lumen: ${particles.length} · GFR ×${params.gfr.toFixed(2)}`);
      },
      dispose() { tex.dispose(); },
    };
  },
};
