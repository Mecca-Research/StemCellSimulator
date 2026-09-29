// Stage 3 (the paper's third order): differentiated cells and their first
// interactions, as three switchable 3D views:
//   (a) myogenesis  - MRF cascade in migrating, aligning myoblasts that fuse
//                      into striated, multinucleated myotubes
//   (b) neurogenesis - Sox2/Pax6 -> Neurogenin -> NeuroD with Delta-Notch
//                      lateral inhibition; axons steered by an analytic
//                      reaction-diffusion chemoattractant field; synapses
//   (c) hormonal feedback - Goodwin loop hypothalamus -> pituitary -> gland
import * as THREE from 'three';
import { cellMaterial, lineMaterial } from '../engine/materials.js';
import { unitSphere, mergeVertices } from '../engine/geometry.js';
import { InstancePool, idColor } from '../engine/instances.js';
import { fmt } from '../engine/chart.js';
import { RNG } from '../models/core/rng.js';
import { rk4Step, rk4Work, integrate } from '../models/core/ode.js';
import { Tissue, Cell } from '../models/morpho/agents3d.js';
import { mrf, MyoCulture, CULTURE_DEFAULTS, tubeCap, TUBE_FLAT } from '../models/pathways/myogenesis.js';
import { paperNgn, neuralGRN, guidance, NeuriteGrowth, OUTGROWTH_DEFAULTS } from '../models/pathways/neurogenesis.js';
import { goodwin, oscillationStats } from '../models/pathways/hormone.js';

const TAU = Math.PI * 2;
const clamp = (x, a, b) => (x < a ? a : x > b ? b : x);
const smooth = (a, b, x) => { const t = clamp((x - a) / (b - a), 0, 1); return t * t * (3 - 2 * t); };
const hex = (c) => `#${new THREE.Color(c).getHexString()}`;

const VIEWS = [
  { value: 'myogenesis', label: 'Myogenesis · MRF cascade & myoblast fusion' },
  { value: 'neurogenesis', label: 'Neurogenesis · Neurogenin, axon guidance, synapses' },
  { value: 'hormone', label: 'Hormonal feedback · Goodwin loop' },
];

export default {
  about: `
    <p>Three programmes of the paper's <b>third order</b> (differentiated cells and their
    first interactions), each as a live 3D model. Pick one under <i>View</i>.</p>
    <p><b>Myogenesis.</b> Myoblasts glide along their long axis on a dish, align nematically
    with their neighbours and run the MRF hierarchy in every nucleus: Wnt/Shh &rarr; Myf5,
    MyoD (with MyoD <i>positive autoregulation</i>, a bistable commitment switch) &rarr;
    myogenin (blocked by mitogens) &rarr; MRF4 &rarr; myosin heavy chain. Myogenin-high
    cells that touch an aligned partner <b>fuse</b> into multinucleated myotubes that
    lengthen (volume conserved), spread their nuclei and develop sarcomeric
    <b>striations</b> (2.2&nbsp;&micro;m period) as MHC accumulates. Fusion index = nuclei in
    myotubes / all nuclei, the standard culture metric.</p>
    <p><b>Neurogenesis.</b> A 3D cluster of Sox2/Pax6 progenitors runs the paper's
    Neurogenin equations inside a Delta&ndash;Notch lateral-inhibition network; cells whose
    NeuroD crosses threshold become neurons, migrate to the surface and grow an axon
    and dendrites. Growth cones read the <b>exact steady state</b> of
    &part;C/&part;t = D&nabla;&sup2;C &minus; &lambda;C + sources secreted by target cells
    (faint shells are iso-concentration surfaces), branch, and form synapses on contact.</p>
    <p><b>Hormonal feedback.</b> A Goodwin loop (hypothalamus &rarr; pituitary &rarr; gland,
    end hormone inhibits the first step). Hormone particles in the vessel are proportional
    to the ODE state; receptors on target cells light up with occupancy. The loop
    oscillates only when the Hill coefficient beats the Griffith/Hopf threshold (8 for
    equal clearance rates).</p>
    <p class="note"><b>Established vs phenomenological.</b> Established: the regulatory
    wiring (MRF hierarchy, Id/mitogen block, reserve cells; Sox2/Pax6/Hes1/Ngn2/NeuroD,
    proneural&ndash;Delta lateral inhibition), the paper's Neurogenin ODE and its analytic
    solution, the screened-Poisson guidance field and its gradient, Berg&ndash;Purcell
    gradient-sensing noise, the Goodwin equations and their Hopf criterion, mass-action
    receptor occupancy. Phenomenological (qualitative, not fitted): all Hill constants and
    rates of the MRF and progenitor networks, fusion hazard and myotube radius law,
    nuclear spreading, growth-cone steering gain/branching, neuronal migration, hormone
    rates (dimensionless axis, minutes), particle transport (a visual proxy, no delay).
    Scales: myogenesis and neurites in &micro;m; the endocrine axis is schematic.</p>`,
  paperRef: 'Paper: "Third Order — Differentiated Cells and Initial Interactions" → Myogenesis (MyoD/Myf5 → myogenin/MRF4 → muscle genes), Feedback loops in hormonal signalling, and the neuronal example Stem cell →(Sox2, Pax6)→ Neural progenitor →(Neurogenin, NeuroD)→ Neuron with d mRNA<sub>Ngn</sub>/dt, d Ngn/dt, ∂C/∂t = D∇²C + R(C) and "Cell Movement = Chemotaxis + Haptotaxis".',

  async create(ctx) {
    const { stage, root, panel, legend, setStatus } = ctx;
    const G = {
      sphere: unitSphere(3),
      sphereMid: unitSphere(2),
      sphereLo: unitSphere(1),
      // open cylinders / hemispheres whose local +x is the long axis (InstancePool.putOriented)
      tube: new THREE.CylinderGeometry(1, 1, 1, 28, 1, true).rotateZ(-Math.PI / 2),
      thin: new THREE.CylinderGeometry(1, 1, 1, 7, 1, true).rotateZ(-Math.PI / 2),
      hemi: new THREE.SphereGeometry(1, 28, 12, 0, TAU, 0, Math.PI / 2).rotateZ(-Math.PI / 2),
    };
    const shared = { stage, root, panel, legend, setStatus, G, mode: stage.mode };
    const factories = { myogenesis: makeMyogenesis, neurogenesis: makeNeurogenesis, hormone: makeHormone };
    const views = {};
    let current = null;

    function activate(id) {
      current = id;
      if (!views[id]) views[id] = factories[id](shared);
      for (const [k, v] of Object.entries(views)) v.group.visible = k === id;
      panel.clear();
      panel.section('View');
      panel.select({ label: 'Process', options: VIEWS, value: id, onChange: (v) => activate(v) });
      views[id].buildPanel(panel);
      views[id].enter();
      views[id].onMode?.(stage.mode, true);
    }
    const requested = ctx.query?.get('view');
    activate(factories[requested] ? requested : 'myogenesis');

    return {
      update(dt) { views[current].update(dt); },
      reset() { views[current].reset(); },
      onMode(mode) { shared.mode = mode; for (const [k, v] of Object.entries(views)) v.onMode?.(mode, k === current); },
    };
  },
};

// =====================================================================
// (a) Myogenesis
// =====================================================================
function makeMyogenesis(S) {
  const { stage, root, legend, setStatus, G } = S;
  const HOURS_PER_SECOND = 1.5;
  const group = new THREE.Group();
  root.add(group);

  const params = {
    start: 'committed', gf: 0.05, wnt: 0, shh: 0, vAuto: mrf.defaults.vAuto,
    align: CULTURE_DEFAULTS.align, v0: CULTURE_DEFAULTS.v0, kFuse: CULTURE_DEFAULTS.kFuse,
    grooves: false, colorBy: 'stage', n: 96,
  };
  let cult, lastSample = -1;

  // ---------------------------------------------------------------- scene
  const R = CULTURE_DEFAULTS.dishR;
  const polar = new THREE.PolarGridHelper(R + 6, 16, 8, 96);
  polar.material = lineMaterial({ color: 0x2a3a50, opacity: 0.5 });
  polar.position.y = -0.2;
  group.add(polar);
  const rim = new THREE.Mesh(new THREE.TorusGeometry(R + 8, 1.6, 10, 160), cellMaterial({ role: 'solid', color: 0x6d7f96, opacity: 1, he: [0.02, 0.1] }));
  rim.rotation.x = Math.PI / 2;
  group.add(rim);
  const grooveGeo = new THREE.BufferGeometry();
  {
    const pts = [];
    for (let z = -R; z <= R; z += 8) { const w = Math.sqrt(Math.max(R * R - z * z, 0)); pts.push(-w, -0.1, z, w, -0.1, z); }
    grooveGeo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
  }
  const grooveLines = new THREE.LineSegments(grooveGeo, lineMaterial({ color: 0x5a7390, opacity: 0.45 }));
  grooveLines.visible = false;
  group.add(grooveLines);

  const bodyMat = cellMaterial({ role: 'membrane', opacity: 0.3, useInst: 1 });
  const tubeMat = cellMaterial({ role: 'cytoplasm', opacity: 0.42, useInst: 1, he: [0.05, 0.95], rimPower: 1.8 });
  const nucMat = cellMaterial({ role: 'nucleus', color: 0xffffff, useInst: 1 });
  const bandMat = cellMaterial({ role: 'reporter', useInst: 1, he: [0.3, 1.4], gain: 0.9 });
  const flashMat = cellMaterial({ role: 'membrane', opacity: 0.05, useInst: 1, emissive: 0.9 });
  const body = new InstancePool(G.sphere, bodyMat, 160, group);
  const nuc = new InstancePool(G.sphereMid, nucMat, 256, group);
  const tubeB = new InstancePool(G.tube, tubeMat, 48, group);
  const tubeC = new InstancePool(G.hemi, tubeMat, 96, group);
  const bands = new InstancePool(G.tube, bandMat, 1024, group);
  const flash = new InstancePool(G.sphereLo, flashMat, 32, group);
  body.mesh.renderOrder = 3; tubeB.mesh.renderOrder = 3; tubeC.mesh.renderOrder = 3; flash.mesh.renderOrder = 4;

  // --------------------------------------------------------------- model
  function build() {
    cult = new MyoCulture({
      n: params.n, start: params.start, gf: params.gf, wnt: params.wnt, shh: params.shh, seed: 20240317,
      mrf: { vAuto: params.vAuto },
      culture: { align: params.align, v0: params.v0, kFuse: params.kFuse, grooves: params.grooves ? 1.2 : 0 },
    });
    lastSample = -1;
    mrfChart?.clear(); fiChart?.clear();
  }

  // --------------------------------------------------------------- colours
  const C_BASE = new THREE.Color(0x7688b8), C_MYOD = new THREE.Color(0x2fd3c9), C_MYOG = new THREE.Color(0xffc247), C_MHC = new THREE.Color(0xff4f86);
  const C_DAPI = new THREE.Color(0x4f7dff), C_MYOGN = new THREE.Color(0x4dff9a);
  const tmp = new THREE.Color(), tmp2 = new THREE.Color();
  function stageColor(s, out) {
    out.copy(C_BASE).lerp(C_MYOD, smooth(0.25, 0.7, s[1]));
    out.lerp(C_MYOG, smooth(0.15, 0.55, s[2]));
    return out.lerp(C_MHC, smooth(0.2, 0.7, s[4]));
  }
  function colorOf(s, th, out) {
    switch (params.colorBy) {
      case 'myod': return heat(s[1], out);
      case 'myog': return heat(s[2], out);
      case 'mhc': return heat(s[4], out);
      case 'orient': return out.setHSL((((th % Math.PI) + Math.PI) % Math.PI) / Math.PI, 0.8, 0.58);
      default: return stageColor(s, out);
    }
  }
  function heat(v, out) { return out.setHSL(0.66 - 0.66 * clamp(v, 0, 1), 0.85, 0.35 + 0.3 * clamp(v, 0, 1)); }
  function nucColor(s, out) {
    // confocal: DAPI, myogenin+ nuclei also carry the (green) myogenin immunostain
    if (S.mode === 'confocal') return out.copy(C_DAPI).lerp(C_MYOGN, 0.8 * smooth(0.25, 0.6, s[2]));
    if (S.mode === 'histology') return out.setRGB(1, 1, 1);
    return out.setRGB(0.86, 0.9, 1).lerp(tmp2.setRGB(1, 0.93, 0.82), smooth(0.25, 0.6, s[2]));
  }

  // --------------------------------------------------------------- drawing
  const pos = [0, 0, 0], dir = [1, 0, 0];
  function draw() {
    body.begin(); nuc.begin(); tubeB.begin(); tubeC.begin(); bands.begin(); flash.begin();
    for (const c of cult.cells) {
      dir[0] = Math.cos(c.th); dir[2] = Math.sin(c.th);
      pos[0] = c.x; pos[1] = c.H; pos[2] = c.z;
      body.putOriented(pos, dir, c.L / 2, c.H, c.W / 2, colorOf(c.s, c.th, tmp));
      pos[1] = c.H * 0.95;
      nuc.putOriented(pos, dir, 5.4, 1.9, 3.1, nucColor(c.s, tmp));
    }
    for (const t of cult.tubes) {
      const cap = tubeCap(t), rh = t.r * TUBE_FLAT, Lb = Math.max(t.len - 2 * cap, 0.01);
      const dx = Math.cos(t.th), dz = Math.sin(t.th);
      const col = colorOf(t.s, t.th, tmp);
      // tubes are more saturated than single cells: MHC-rich cytoplasm
      dir[0] = dx; dir[2] = dz;
      pos[0] = t.x; pos[1] = rh; pos[2] = t.z;
      tubeB.putOriented(pos, dir, Lb, rh, t.r, col);
      for (let sg = -1; sg <= 1; sg += 2) {
        pos[0] = t.x + sg * dx * Lb / 2; pos[2] = t.z + sg * dz * Lb / 2;
        dir[0] = sg * dx; dir[2] = sg * dz;
        tubeC.putOriented(pos, dir, cap, rh, t.r, col);
      }
      dir[0] = dx; dir[2] = dz;
      // nuclei: central chains that spread as the myotube matures
      const usable = Math.max(t.len / 2 - cap * 0.6, 1);
      for (const nu of t.nuclei) {
        const along = nu.u * usable, lat = nu.lat * t.r * 0.9;
        pos[0] = t.x + dx * along - dz * lat; pos[1] = rh * (1 + 0.35 * nu.h); pos[2] = t.z + dz * along + dx * lat;
        nuc.putOriented(pos, dir, 5.6, 2.1, 2.9, nucColor(t.s, tmp));
      }
      // sarcomeric striations (A-band spacing 2.2 um) appear with MHC, from the centre outwards
      const extent = smooth(0.3, 0.75, t.s[4]);
      if (extent > 0.02) {
        const half = (Lb / 2) * extent;
        const bandCol = bandColor(t.s[4], tmp2);
        for (let x = -half; x <= half; x += 2.2) {
          pos[0] = t.x + dx * x; pos[1] = rh; pos[2] = t.z + dz * x;
          bands.putOriented(pos, dir, 0.75, rh * 1.012, t.r * 1.012, bandCol);
        }
      }
    }
    for (const e of cult.events) {
      const a = (cult.t - e.t) / 1.5;
      const r = 6 + 22 * a;
      flash.put(e.x, 4, e.z, r, r * 0.4, r, tmp.setRGB(1, 0.92, 0.55).multiplyScalar(1 - a));
    }
    body.end(); nuc.end(); tubeB.end(); tubeC.end(); bands.end(); flash.end();
  }
  function bandColor(M, out) {
    if (S.mode === 'confocal') return out.setRGB(1, 0.25, 0.45).multiplyScalar(0.35 + 0.65 * smooth(0.3, 0.8, M));
    if (S.mode === 'histology') return out.setRGB(0.75, 0.55, 0.8);
    return out.setRGB(0.55, 0.12, 0.28);
  }

  // --------------------------------------------------------------- panel
  let mrfChart, fiChart, live, thrReadout;
  function buildPanel(panel) {
    panel.section('Culture');
    panel.select({
      label: 'Starting population', value: params.start,
      options: [
        { value: 'committed', label: 'Committed myoblasts (C2C12-like, MyoD on)' },
        { value: 'naive', label: 'Naive paraxial mesoderm (needs Wnt/Shh)' },
      ],
      onChange: (v) => {
        params.start = v;
        if (v === 'naive') { params.wnt = Math.max(params.wnt, 0.8); params.shh = Math.max(params.shh, 0.8); wntS.set(params.wnt); shhS.set(params.shh); }
        build();
      },
    });
    const gfS = panel.slider({ label: 'Mitogens / growth factors (serum)', min: 0, max: 1, value: params.gf, onChange: (v) => { params.gf = v; cult.env.gf = v; } });
    panel.buttons([
      { label: 'Growth medium', onClick: () => { params.gf = 1; cult.env.gf = 1; gfS.set(1); } },
      { label: 'Differentiation medium', primary: true, onClick: () => { params.gf = 0.05; cult.env.gf = 0.05; gfS.set(0.05); } },
    ]);
    const wntS = panel.slider({ label: 'Wnt (dorsal neural tube / ectoderm)', min: 0, max: 1, value: params.wnt, onChange: (v) => { params.wnt = v; cult.env.wnt = v; } });
    const shhS = panel.slider({ label: 'Shh (notochord / floor plate)', min: 0, max: 1, value: params.shh, onChange: (v) => { params.shh = v; cult.env.shh = v; } });
    panel.slider({
      label: 'MyoD autoregulation v', min: 0.5, max: 1.1, value: params.vAuto,
      onChange: (v) => { params.vAuto = v; cult.mp.vAuto = v; showThresholds(); },
    });
    thrReadout = panel.readouts(['switch turn-on input', 'switch turn-off input']);
    showThresholds();
    panel.note('Fold points of s = D − v·D⁴/(K⁴+D⁴): commitment needs input above the turn-on value; a negative turn-off value means the committed state is kept after Wnt/Shh are withdrawn.');

    panel.section('Mechanics & fusion');
    panel.slider({ label: 'Nematic alignment γ (1/h)', min: 0, max: 3, value: params.align, onChange: (v) => { params.align = v; cult.p.align = v; } });
    panel.slider({ label: 'Gliding speed (µm/h)', min: 0, max: 40, value: params.v0, onChange: (v) => { params.v0 = v; cult.p.v0 = v; } });
    panel.slider({ label: 'Fusion hazard (1/h)', min: 0, max: 1.5, value: params.kFuse, onChange: (v) => { params.kFuse = v; cult.p.kFuse = v; } });
    panel.toggle({ label: 'Micro-grooved substrate (contact guidance along x)', value: params.grooves, onChange: (v) => { params.grooves = v; cult.p.grooves = v ? 1.2 : 0; grooveLines.visible = v; } });
    panel.select({
      label: 'Colour cells by', value: params.colorBy,
      options: [
        { value: 'stage', label: 'MRF stage (MyoD → myogenin → MHC)' },
        { value: 'myod', label: 'MyoD level' },
        { value: 'myog', label: 'Myogenin level' },
        { value: 'mhc', label: 'MHC level' },
        { value: 'orient', label: 'Orientation (nematic director)' },
      ],
      onChange: (v) => { params.colorBy = v; updateLegend(); },
    });
    panel.buttons([{ label: 'Restart culture', primary: true, onClick: () => build() }]);

    panel.section('MRF cascade (per nucleus)');
    panel.equation('dMyoD/dt = k_D [ a_E·E(Wnt,Shh) + a_F·Myf5 + v·MyoD⁴/(K⁴+MyoD⁴) − MyoD ]');
    panel.equation('dMyog/dt = k_G [ MyoD²/(K²+MyoD²) · K_gf²/(K_gf²+GF²) − Myog ]');
    panel.equation('dMHC/dt = k_M [ h(Myog + MRF4) − MHC ]');
    panel.equation('dθᵢ/dt = γ Σⱼ sin 2(θⱼ − θᵢ) + τ_steric/ζ + √(2D_r) ξ');
    mrfChart = panel.chart({
      title: 'MRF levels (population mean over nuclei)', xLabel: 'hours', yRange: [0, 1],
      series: [
        { name: 'Myf5', color: '#8fa7ff' }, { name: 'MyoD', color: hex(C_MYOD) }, { name: 'myogenin', color: hex(C_MYOG) },
        { name: 'MRF4', color: '#ff9b5e' }, { name: 'MHC', color: hex(C_MHC) },
      ],
    });
    fiChart = panel.chart({
      title: 'Fusion index & nematic order', xLabel: 'hours', yRange: [0, 1],
      series: [{ name: 'fusion index', color: '#ffffff' }, { name: 'MHC+ nuclei', color: hex(C_MHC), dash: [4, 3] }, { name: 'nematic S', color: '#7fd1ff' }],
    });
    live = panel.readouts(['myoblasts', 'myotubes', 'nuclei in myotubes', 'largest myotube', 'fusion index', 'reserve (unfused) cells']);
    panel.note('Fusion index = nuclei in cells with ≥ 2 nuclei / all nuclei. C2C12 cultures typically reach 30–70 % after 3–5 days in differentiation medium; ~20–30 % "reserve cells" stay mononucleated (Yoshida et al. 1998).');
  }
  function showThresholds() {
    const th = mrf.thresholds({ ...mrf.defaults, vAuto: params.vAuto });
    if (!thrReadout) return;
    thrReadout.set('switch turn-on input', th.bistable ? fmt(th.on) : 'monostable');
    thrReadout.set('switch turn-off input', th.bistable ? (th.off < 0 ? `${fmt(th.off)} (irreversible)` : fmt(th.off)) : '—');
  }

  function updateLegend() {
    const common = [
      { color: '#dfe6ff', label: 'nuclei (myotube nuclei: central chains)' },
      { color: '#8c1f47', label: 'sarcomeric striations (2.2 µm period)' },
    ];
    const by = {
      stage: [
        { color: hex(C_BASE), label: 'uncommitted' }, { color: hex(C_MYOD), label: 'MyoD+ committed myoblast' },
        { color: hex(C_MYOG), label: 'myogenin+ (fusion-competent)' }, { color: hex(C_MHC), label: 'MHC+ myocyte / myotube' },
      ],
      myod: [{ color: '#2233aa', label: 'MyoD low' }, { color: '#ff5a50', label: 'MyoD high' }],
      myog: [{ color: '#2233aa', label: 'myogenin low' }, { color: '#ff5a50', label: 'myogenin high' }],
      mhc: [{ color: '#2233aa', label: 'MHC low' }, { color: '#ff5a50', label: 'MHC high' }],
      orient: [{ color: '#ff4d4d', label: 'hue = long-axis angle (0–180°)' }],
    }[params.colorBy];
    legend([...by, ...common, { color: '#fff0a0', label: 'fusion event' }]);
  }

  // --------------------------------------------------------------- picking
  const pickables = [
    {
      get object() { return body.mesh; },
      info: (hit) => {
        const c = cult.cells[hit.instanceId];
        if (!c) return null;
        return `Myoblast ${c.id}${c.reserve ? ' · reserve cell' : ''}\n${mrf.stage(c.s)}\n` +
          `Myf5 ${fmt(c.s[0])} · MyoD ${fmt(c.s[1])} · myogenin ${fmt(c.s[2])}\nMRF4 ${fmt(c.s[3])} · MHC ${fmt(c.s[4])}`;
      },
    },
    {
      get object() { return tubeB.mesh; },
      info: (hit) => {
        const t = cult.tubes[hit.instanceId];
        if (!t) return null;
        return `Myotube ${t.id} · ${t.nuclei.length} nuclei\nlength ${fmt(t.len)} µm · width ${fmt(2 * t.r)} µm · age ${fmt(t.age)} h\n` +
          `myogenin ${fmt(t.s[2])} · MRF4 ${fmt(t.s[3])} · MHC ${fmt(t.s[4])}${t.s[4] > 0.3 ? '\nstriated (sarcomeres assembling)' : ''}`;
      },
    },
  ];

  build();
  draw();
  const mean = new Float64Array(5);

  return {
    group,
    buildPanel,
    enter() {
      stage.frame([0, 0, 26], R * 1.02, new THREE.Vector3(0.1, 1, 0.78));
      stage.clipAxis.set(0, 1, 0);
      stage.clipRange = [-1, 16];
      stage.setPickables(pickables);
      updateLegend();
      lastSample = -1;
    },
    update(dt) {
      const dtH = dt * HOURS_PER_SECOND;
      if (dtH > 0) {
        const sub = Math.max(1, Math.ceil(dtH / 0.1));
        for (let k = 0; k < sub; k++) cult.step(dtH / sub);
      }
      if (cult.t - lastSample >= 0.5 || lastSample < 0) {
        lastSample = cult.t;
        cult.meanState(mean);
        mrfChart.push(cult.t, Array.from(mean));
        let nMHC = 0, nAll = 0;
        for (const c of cult.cells) { nAll++; if (c.s[4] > 0.5) nMHC++; }
        for (const t of cult.tubes) { nAll += t.nuclei.length; if (t.s[4] > 0.5) nMHC += t.nuclei.length; }
        fiChart.push(cult.t, [cult.fusionIndex(), nMHC / Math.max(nAll, 1), cult.nematicOrder()]);
        const inT = cult.tubes.reduce((s, t) => s + t.nuclei.length, 0);
        live.set('myoblasts', cult.cells.length);
        live.set('myotubes', cult.tubes.length);
        live.set('nuclei in myotubes', `${inT} (mean ${fmt(inT / Math.max(cult.tubes.length, 1))} / tube)`);
        live.set('largest myotube', cult.tubes.length ? `${Math.max(...cult.tubes.map((t) => t.nuclei.length))} nuclei` : '–');
        live.set('fusion index', `${(100 * cult.fusionIndex()).toFixed(1)} %`);
        live.set('reserve (unfused) cells', cult.cells.filter((c) => c.reserve).length);
      }
      draw();
      setStatus(`t = ${cult.t.toFixed(1)} h · ${cult.cells.length} myoblasts · ${cult.tubes.length} myotubes · fusion index ${(100 * cult.fusionIndex()).toFixed(0)} %`);
    },
    reset() { build(); },
    onMode(mode) {
      polar.visible = mode !== 'confocal';
      rim.visible = mode !== 'confocal';
      polar.material.opacity = mode === 'histology' ? 0.16 : 0.5;
      grooveLines.material.opacity = mode === 'histology' ? 0.2 : 0.45;
    },
  };
}

// =====================================================================
// (b) Neurogenesis
// =====================================================================
function makeNeurogenesis(S) {
  const { stage, root, legend, setStatus, G } = S;
  const HOURS_PER_SECOND = 0.8;
  const group = new THREE.Group();
  root.add(group);

  const params = {
    coupling: neuralGRN.defaults.coupling, notchExt: 0, drive: 1,
    ...paperNgn.defaults,
    D: guidance.defaults.D, lambda: guidance.defaults.lambda, Q: guidance.defaults.Q,
    kappa: OUTGROWTH_DEFAULTS.kappa, noise: OUTGROWTH_DEFAULTS.noise, branch: OUTGROWTH_DEFAULTS.branchRate, hapto: 0,
    showIso: true, colorBy: 'type', nCells: 110,
  };
  const N_TARGETS = 6;
  const R_CELL = 6.5;

  // targets: glowing chemoattractant-secreting cells around the cluster (golden-spiral directions)
  const targets = [];
  {
    const rng = new RNG(77);
    for (let i = 0; i < N_TARGETS; i++) {
      const y = 1 - (2 * (i + 0.5)) / N_TARGETS * 0.8 - 0.1;
      const r = Math.sqrt(1 - y * y), a = i * 2.39996 + 0.4;
      const d = 108 + 26 * rng.next();
      targets.push({ id: i, p: [Math.cos(a) * r * d, y * d * 0.8, Math.sin(a) * r * d], radius: 9, Q: params.Q, pulse: rng.next() * TAU });
    }
  }
  const gp = () => ({ ...guidance.defaults, D: params.D, lambda: params.lambda, Q: params.Q });

  // --------------------------------------------------------------- materials & pools
  // confocal gains < 1: ~100 overlapping cells and dense neurites would saturate additive blending
  const somaMat = cellMaterial({ role: 'membrane', opacity: 0.3, useInst: 1, gain: 0.28 });
  const nucMat = cellMaterial({ role: 'nucleus', color: 0xffffff, useInst: 1, gain: 0.32 });
  const targetMat = cellMaterial({ role: 'membrane', opacity: 0.35, useInst: 1, emissive: 0.35 });
  const targetCoreMat = cellMaterial({ role: 'reporter', useInst: 1, emissive: 0.8, he: [0.6, 0.3], gain: 0.8 });
  const neuriteMat = cellMaterial({ role: 'reporter', useInst: 1, he: [0.1, 0.8], emissive: 0.08, gain: 0.32 });
  const coneMat = cellMaterial({ role: 'reporter', useInst: 1, emissive: 0.7, he: [0.35, 0.6] });
  const flashMat = cellMaterial({ role: 'membrane', opacity: 0.06, useInst: 1, emissive: 1 });
  const soma = new InstancePool(G.sphere, somaMat, 160, group);
  const nuc = new InstancePool(G.sphereMid, nucMat, 160, group);
  const tBody = new InstancePool(G.sphere, targetMat, 8, group);
  const tCore = new InstancePool(G.sphereMid, targetCoreMat, 8, group);
  const segs = new InstancePool(G.thin, neuriteMat, 4096, group);
  const cones = new InstancePool(G.sphereMid, coneMat, 128, group);
  const filo = new InstancePool(G.thin, coneMat, 512, group);
  const boutons = new InstancePool(G.sphereMid, coneMat, 64, group);
  const flashes = new InstancePool(G.sphereLo, flashMat, 16, group);
  soma.mesh.renderOrder = 3; tBody.mesh.renderOrder = 3; flashes.mesh.renderOrder = 5;

  // iso-concentration surfaces of the summed field around each target (star-shaped from the source)
  const ISO_LEVELS = [{ f: 1, color: 0x5fd8ff, opacity: 0.01 }, { f: 0.6, color: 0x3d8bff, opacity: 0.006 }];
  const isoMeshes = [];
  const ico = new THREE.IcosahedronGeometry(1, 7);
  ico.deleteAttribute('normal'); ico.deleteAttribute('uv');
  const isoBase = mergeVertices(ico);
  ico.dispose();
  const isoDirs = isoBase.attributes.position.array.slice();
  for (const lv of ISO_LEVELS) {
    const mat = cellMaterial({ role: 'membrane', color: lv.color, stain: lv.color, opacity: lv.opacity, rimPower: 7, gain: 0.35 });
    for (const tg of targets) {
      const g = isoBase.clone();
      const m = new THREE.Mesh(g, mat);
      m.position.set(...tg.p);
      m.renderOrder = 6;
      m.userData = { level: lv, target: tg };
      group.add(m);
      isoMeshes.push(m);
    }
  }
  let isoDirty = true, isoLast = 0;
  function rebuildIso() {
    const p = gp();
    const src = targets.map((t) => ({ p: t.p, radius: t.radius, Q: p.Q }));
    const x = [0, 0, 0];
    for (const m of isoMeshes) {
      const level = m.userData.level.f * p.Kd; // nM: fraction of the receptor K_d
      const tg = m.userData.target, c = tg.p;
      const pos = m.geometry.attributes.position;
      const arr = pos.array;
      const f = (dx, dy, dz, r) => { x[0] = c[0] + dx * r; x[1] = c[1] + dy * r; x[2] = c[2] + dz * r; return guidance.concentration(x, src, p); };
      let clipped = 0;
      for (let i = 0; i < arr.length; i += 3) {
        const dx = isoDirs[i], dy = isoDirs[i + 1], dz = isoDirs[i + 2];
        // stay inside this source's Voronoi cell (the summed field is star-shaped there)
        let rMax = 400;
        for (const o of targets) {
          if (o === tg) continue;
          const ox = o.p[0] - c[0], oy = o.p[1] - c[1], oz = o.p[2] - c[2];
          const along = dx * ox + dy * oy + dz * oz;
          if (along > 0) rMax = Math.min(rMax, (0.5 * (ox * ox + oy * oy + oz * oz)) / along);
        }
        // march to the first crossing, then bisect inside the bracket
        let lo = tg.radius, hi = lo, r;
        if (f(dx, dy, dz, lo) < level) r = lo;
        else {
          while (hi < rMax && f(dx, dy, dz, hi) >= level) { lo = hi; hi = Math.min(hi + 6, rMax); if (hi === lo) break; }
          if (f(dx, dy, dz, hi) >= level) { r = rMax; clipped++; }
          else {
            for (let k = 0; k < 18; k++) { const mid = 0.5 * (lo + hi); if (f(dx, dy, dz, mid) > level) lo = mid; else hi = mid; }
            r = 0.5 * (lo + hi);
          }
        }
        arr[i] = dx * r; arr[i + 1] = dy * r; arr[i + 2] = dz * r;
      }
      pos.needsUpdate = true;
      m.geometry.computeVertexNormals();
      m.geometry.computeBoundingSphere();
      // a shell that reaches the Voronoi boundary has merged with its neighbours: not star-shaped, skip it
      m.visible = params.showIso && S.mode !== 'histology' && clipped < 0.01 * (arr.length / 3) && m.geometry.boundingSphere.radius > tg.radius * 1.2;
    }
    isoDirty = false;
  }

  // --------------------------------------------------------------- model state
  let rng, tis, Y, committed, growth, work, t, lastSample, segDrawn, idIndex;
  function grnParams() {
    return { ...neuralGRN.defaults, coupling: params.coupling, notchExt: params.notchExt, drive: params.drive, ktx: params.ktx, kdm: params.kdm, ktl: params.ktl, kdp: params.kdp };
  }
  let P = grnParams();
  function build() {
    rng = new RNG(4242);
    tis = new Tissue({
      dims: 3, growth: false, noise: 0.015, muRep: 10, adhesion: () => 2.2, rng,
      external: (c, out) => {
        const st = c.state;
        const r = Math.hypot(c.p[0], c.p[1], c.p[2]) || 1;
        if (st.committed && !st.anchored && t - st.tCommit > 0.8) {
          // newborn neurons migrate to the basal (outer) surface of the cluster: phenomenological
          const f = 6;
          out[0] += (f * c.p[0]) / r; out[1] += (f * c.p[1]) / r; out[2] += (f * c.p[2]) / r;
        } else if (!st.committed) {
          // progenitors stay compact (apical adhesion, phenomenological weak centring)
          out[0] -= 0.05 * c.p[0]; out[1] -= 0.05 * c.p[1]; out[2] -= 0.05 * c.p[2];
        }
      },
    });
    for (let i = 0; i < params.nCells; i++) {
      const d = rng.onSphere(), r = 36 * Math.cbrt(rng.next());
      tis.add(new Cell({ p: d.map((v) => v * r), R: R_CELL, state: { i, committed: false, anchored: false, sprouted: false, tCommit: 0 } }));
    }
    for (let k = 0; k < 240; k++) tis.step(0.05);
    tis.t = 0;
    idIndex = new Map(tis.cells.map((c, i) => [c.id, i]));
    Y = tis.cells.map(() => neuralGRN.initial(rng));
    committed = tis.cells.map(() => false);
    growth = new NeuriteGrowth({
      rng: new RNG(99),
      outgrowth: { kappa: params.kappa, noise: params.noise, branchRate: params.branch, hapto: params.hapto, maxAxonBranches: 3, maxAxonLen: 320, step: 2.5, maxSegments: 30000 },
      guidance: gp(),
    });
    work = neuralGRN.work();
    t = 0; lastSample = -1; segDrawn = 0;
    segs.begin(); segs.end();
    for (const ch of [netChart, fracChart, wireChart]) ch?.clear();
  }

  // --------------------------------------------------------------- colours
  const C_PROG = new THREE.Color(0x35d69e), C_NGN = new THREE.Color(0xffa63d), C_NEU = new THREE.Color(0xb57bff);
  const C_AXON = new THREE.Color(0xc9b4ff), C_DEND = new THREE.Color(0xff9fd0), C_TGT = new THREE.Color(0x49d2ff);
  const C_CONE = new THREE.Color(0xfff07a), C_BOUTON = new THREE.Color(0xffe4a0);
  // confocal: immunostain-like palette (Sox2 green, beta-III tubulin/TUJ1 red neurons and neurites)
  const F_PROG = new THREE.Color(0x2bff6a), F_NEU = new THREE.Color(0xff2f4f), F_AXON = new THREE.Color(0xff3552), F_DEND = new THREE.Color(0xff6a3a);
  const tmp = new THREE.Color(), tmp2 = new THREE.Color();
  function cellColor(i, out) {
    const y = Y[i];
    if (params.colorBy === 'hes') return out.setHSL(0.66 - 0.66 * clamp(y[6], 0, 1), 0.85, 0.5);
    if (params.colorBy === 'brainbow' && committed[i]) return idColor(tis.cells[i].id, 0.8, 0.62, out);
    const conf = S.mode === 'confocal';
    out.copy(conf ? F_PROG : C_PROG).lerp(C_NGN, smooth(0.15, 0.6, y[3]));
    if (committed[i]) out.lerp(conf ? F_NEU : C_NEU, smooth(0.5, 0.85, y[4]));
    return out;
  }
  function segColor(sg, out) {
    if (params.colorBy === 'brainbow') return idColor(sg.neuron, 0.8, 0.62, out);
    if (S.mode === 'confocal') return out.copy(sg.kind === 0 ? F_AXON : F_DEND);
    return out.copy(sg.kind === 0 ? C_AXON : C_DEND);
  }
  function nucColor(i, out) {
    if (S.mode === 'confocal') return out.setRGB(0.3, 0.45, 1).lerp(tmp2.setRGB(0.3, 1, 0.55), 0.7 * smooth(0.4, 0.8, Y[i][0]));
    if (S.mode === 'histology') return out.setRGB(1, 1, 1);
    return out.setRGB(0.88, 0.9, 1);
  }

  // --------------------------------------------------------------- simulation step
  function simStep(dt) {
    tis.step(dt);
    // pinned neurons (their neurites are anchored to them)
    for (const c of tis.cells) if (c.state.anchored) { c.p[0] = c.state.ax; c.p[1] = c.state.ay; c.p[2] = c.state.az; }
    // Delta-Notch contacts from the mechanical neighbour graph
    const nbrs = tis.cells.map((c) => (tis.neighbours.get(c.id) ?? []).map((id) => idIndex.get(id)));
    const fresh = neuralGRN.step(Y, committed, nbrs, dt, P, work);
    for (const i of fresh) { const st = tis.cells[i].state; st.committed = true; st.tCommit = t; }
    for (let i = 0; i < tis.cells.length; i++) {
      const c = tis.cells[i], st = c.state;
      if (!st.committed) continue;
      c.R += (5.2 - c.R) * Math.min(1, dt * 0.5); // neuronal somata are smaller
      const r = Math.hypot(c.p[0], c.p[1], c.p[2]);
      if (!st.anchored && (r > 40 || t - st.tCommit > 6) && t - st.tCommit > 1.5) {
        st.anchored = true; st.ax = c.p[0]; st.ay = c.p[1]; st.az = c.p[2];
      }
      if (st.anchored && !st.sprouted && t - st.tCommit > 2.2) {
        st.sprouted = true;
        growth.addNeuron(c.id, c.p, c.R, [c.p[0] / r, c.p[1] / r, c.p[2] / r]);
      }
    }
    growth.step(dt, sourceList);
    t += dt;
  }
  let sourceList = targets.map((tg) => ({ id: tg.id, p: tg.p, radius: tg.radius, Q: params.Q }));
  function applyGuidance() {
    growth.gp = gp();
    sourceList = targets.map((tg) => ({ id: tg.id, p: tg.p, radius: tg.radius, Q: params.Q }));
    isoDirty = true;
    lRead?.set('decay length ℓ = √(D/λ)', `${fmt(guidance.length(growth.gp))} µm`);
  }

  // --------------------------------------------------------------- drawing
  const dir = [1, 0, 0], pos = [0, 0, 0];
  function draw(time) {
    soma.begin(); nuc.begin(); tBody.begin(); tCore.begin(); cones.begin(); filo.begin(); boutons.begin(); flashes.begin();
    for (let i = 0; i < tis.cells.length; i++) {
      const c = tis.cells[i];
      soma.put(c.p[0], c.p[1], c.p[2], c.R * 1.05, c.R * 1.05, c.R * 1.05, cellColor(i, tmp));
      const k = committed[i] ? 0.62 : 0.66;
      nuc.put(c.p[0], c.p[1], c.p[2], c.R * k, c.R * k * 0.92, c.R * k, nucColor(i, tmp));
    }
    for (const tg of targets) {
      const pulse = 1 + 0.06 * Math.sin(time * 2.2 + tg.pulse);
      tBody.put(tg.p[0], tg.p[1], tg.p[2], tg.radius * pulse, tg.radius * pulse, tg.radius * pulse, C_TGT);
      tCore.put(tg.p[0], tg.p[1], tg.p[2], tg.radius * 0.55, tg.radius * 0.55, tg.radius * 0.55, tmp.copy(C_TGT).lerp(tmp2.setRGB(1, 1, 1), 0.45));
    }
    // neurite segments are immutable: only upload the new ones
    const all = growth.segments;
    if (segDrawn > all.length) segDrawn = 0;
    if (all.length > segDrawn) {
      segs.n = segDrawn;
      for (let s = segDrawn; s < all.length; s++) putSegment(all[s]);
      segs.end();
      const im = segs.mesh.instanceMatrix, ic = segs.mesh.instanceColor;
      im.clearUpdateRanges(); ic.clearUpdateRanges();
      im.addUpdateRange(segDrawn * 16, (all.length - segDrawn) * 16);
      ic.addUpdateRange(segDrawn * 3, (all.length - segDrawn) * 3);
      segDrawn = all.length;
    }
    // growth cones with flickering filopodia
    for (const tp of growth.tips) {
      if (!tp.active) continue;
      const big = tp.kind === 0;
      const s = big ? 1 : 0.6;
      cones.putOriented(tp.pos, tp.dir, 2.6 * s, 1.5 * s, 2.2 * s, C_CONE);
      if (!big) continue;
      const u = perp(tp.dir, _u), w = cross(tp.dir, u, _w);
      for (let f = 0; f < 4; f++) {
        const ph = tp.born * 13.1 + f * 1.7;
        // direction within a cone around the growth direction
        const a = 0.9 * Math.sin(ph) + 0.35 * Math.sin(time * 1.3 + ph);
        const b = 0.9 * Math.cos(ph * 1.3) + 0.35 * Math.cos(time * 1.1 + ph);
        dir[0] = tp.dir[0] + a * u[0] + b * w[0]; dir[1] = tp.dir[1] + a * u[1] + b * w[1]; dir[2] = tp.dir[2] + a * u[2] + b * w[2];
        const n = Math.hypot(dir[0], dir[1], dir[2]); dir[0] /= n; dir[1] /= n; dir[2] /= n;
        const L = 5 + 2.5 * Math.sin(time * 2 + ph);
        pos[0] = tp.pos[0] + dir[0] * L / 2; pos[1] = tp.pos[1] + dir[1] * L / 2; pos[2] = tp.pos[2] + dir[2] * L / 2;
        filo.putOriented(pos, dir, L, 0.22, 0.22, C_CONE);
      }
    }
    for (const sy of growth.synapses) {
      boutons.put(sy.pos[0], sy.pos[1], sy.pos[2], 2.3, 2.3, 2.3, C_BOUTON);
      const age = growth.t - sy.t;
      if (age < 1.2) {
        const a = age / 1.2, r = 3 + 16 * a;
        flashes.put(sy.pos[0], sy.pos[1], sy.pos[2], r, r, r, tmp.setRGB(1, 0.95, 0.6).multiplyScalar(1 - a));
      }
    }
    soma.end(); nuc.end(); tBody.end(); tCore.end(); cones.end(); filo.end(); boutons.end(); flashes.end();
  }
  const _mid = [0, 0, 0], _d = [0, 0, 0], _u = [0, 0, 0], _w = [0, 0, 0];
  function putSegment(sg) {
    _d[0] = sg.b[0] - sg.a[0]; _d[1] = sg.b[1] - sg.a[1]; _d[2] = sg.b[2] - sg.a[2];
    const L = Math.hypot(_d[0], _d[1], _d[2]) || 1e-3;
    _d[0] /= L; _d[1] /= L; _d[2] /= L;
    _mid[0] = (sg.a[0] + sg.b[0]) / 2; _mid[1] = (sg.a[1] + sg.b[1]) / 2; _mid[2] = (sg.a[2] + sg.b[2]) / 2;
    // drawn 1.15x thicker than the physical calibre so micron-thin neurites stay visible
    segs.putOriented(_mid, _d, L * 1.08, sg.r * 1.15, sg.r * 1.15, segColor(sg, tmp));
  }
  function redrawSegments() { segDrawn = 0; segs.begin(); segs.end(); }

  // --------------------------------------------------------------- panel
  let paperChart, netChart, fracChart, wireChart, live, lRead, errRead;
  function paperCurves() {
    const p = { ktx: params.ktx, kdm: params.kdm, ktl: params.ktl, kdp: params.kdp };
    const { t: ts, y } = integrate(paperNgn.rhs(p), [0, 0], 0, 6, 0.01, { every: 5 });
    const ex = ts.map((ti) => paperNgn.analytic(ti, p));
    let err = 0;
    y.forEach((v, i) => { for (let k = 0; k < 2; k++) err = Math.max(err, Math.abs(v[k] - ex[i][k]) / Math.max(Math.abs(ex[i][k]), 1e-6)); });
    paperChart.set(ts, [y.map((v) => v[0]), y.map((v) => v[1]), ex.map((v) => v[0]), ex.map((v) => v[1])]);
    const [ms, Ps] = paperNgn.steadyState(p);
    errRead.set('steady state m*, P*', `${fmt(ms)}, ${fmt(Ps)}`);
    errRead.set('max |RK4 − exact| / exact', err.toExponential(1));
  }
  function buildPanel(panel) {
    panel.section("Paper's Neurogenin equations");
    panel.equation('d mRNA_Ngn/dt = k_tx − k_dm · mRNA_Ngn');
    panel.equation('d Ngn/dt = k_tl · mRNA_Ngn − k_dp · Ngn');
    panel.equation('Ngn(t) = P* + A e^(−k_dm t) + (P₀ − P* − A) e^(−k_dp t),  A = k_tl (m₀ − m*)/(k_dp − k_dm)');
    const setRate = (k) => (v) => { params[k] = v; P = grnParams(); paperCurves(); };
    panel.slider({ label: 'k_tx transcription (1/h, max)', min: 0.2, max: 3, value: params.ktx, onChange: setRate('ktx') });
    panel.slider({ label: 'k_dm mRNA degradation (1/h)', min: 0.2, max: 3, value: params.kdm, onChange: setRate('kdm') });
    panel.slider({ label: 'k_tl translation (1/h)', min: 0.2, max: 4, value: params.ktl, onChange: setRate('ktl') });
    panel.slider({ label: 'k_dp protein degradation (1/h)', min: 0.2, max: 3, value: params.kdp, onChange: setRate('kdp') });
    paperChart = panel.chart({
      title: 'Constant k_tx from zero: RK4 (solid) vs analytic (dashed)', xLabel: 'hours',
      series: [
        { name: 'mRNA', color: '#ffd166' }, { name: 'Ngn', color: '#ff7b54' },
        { name: 'mRNA exact', color: '#ffffff', dash: [3, 4], width: 1.2 }, { name: 'Ngn exact', color: '#ffffff', dash: [3, 4], width: 1.2 },
      ],
    });
    errRead = panel.readouts(['steady state m*, P*', 'max |RK4 − exact| / exact']);
    paperCurves();

    panel.section('Progenitor network');
    panel.note('In each cell k_tx is modulated by Pax6 (activator) and Hes1 (repressor): k_tx,eff = k_tx · h⁺(Pax6) · h⁻(Hes1); Hes1 is driven by Delta on neighbouring cells (lateral inhibition).');
    panel.slider({ label: 'Delta–Notch coupling (0 = DAPT)', min: 0, max: 2, value: params.coupling, onChange: (v) => { params.coupling = v; P = grnParams(); } });
    panel.slider({ label: 'Exogenous Notch activation', min: 0, max: 0.12, step: 0.002, value: params.notchExt, onChange: (v) => { params.notchExt = v; P = grnParams(); } });
    panel.slider({ label: 'Neurogenic drive (× k_tx)', min: 0.5, max: 1.6, value: params.drive, onChange: (v) => { params.drive = v; P = grnParams(); } });
    panel.select({
      label: 'Colour cells by', value: params.colorBy,
      options: [
        { value: 'type', label: 'Fate (Sox2 progenitor → Ngn+ → NeuroD+ neuron)' },
        { value: 'hes', label: 'Hes1 (Notch activity)' },
        { value: 'brainbow', label: 'Brainbow (one hue per neuron)' },
      ],
      onChange: (v) => { params.colorBy = v; updateLegend(); redrawSegments(); },
    });
    netChart = panel.chart({
      title: 'Progenitors: mean Neurogenin mRNA / protein, Hes1, Sox2', xLabel: 'hours', yRange: [0, 1.4],
      series: [{ name: 'Ngn mRNA', color: '#ffd166' }, { name: 'Ngn', color: '#ff7b54' }, { name: 'Hes1', color: '#6fa8ff' }, { name: 'Sox2', color: hex(C_PROG) }],
    });
    fracChart = panel.chart({
      title: 'Fraction of cells', xLabel: 'hours', yRange: [0, 1],
      series: [{ name: 'neurons (NeuroD+)', color: hex(C_NEU) }, { name: 'Sox2+ progenitors', color: hex(C_PROG) }],
    });

    panel.section('Axon guidance: ∂C/∂t = D∇²C − λC + Σ Q δ(x − xₛ)');
    panel.note("The paper's reaction term R(C) is taken as first-order removal plus point secretion by the target cells, R(C) = −λC + Σ Q δ(x − xₛ). In unbounded 3D its steady state is the screened-Poisson (Yukawa) solution:");
    panel.equation('C(r) = Σₛ Q / (4πD rₛ) · exp(−rₛ/ℓ),   ℓ = √(D/λ)');
    panel.equation('∇C = −Σₛ Q e^(−rₛ/ℓ)(1 + rₛ/ℓ) / (4πD rₛ³) · (x − xₛ)');
    panel.equation('SNR = Δp √(N / 4p(1−p)),  p = C/(C+K_d),  Δp = K_d w|∇C|/(C+K_d)²');
    panel.equation('d ← norm( d + κ·SNR/(1+SNR)·∇C/|∇C| + η·h_ECM + σ ξ )');
    panel.slider({ label: 'Diffusion D (µm²/s)', min: 0.5, max: 50, log: true, value: params.D, onChange: (v) => { params.D = v; applyGuidance(); } });
    panel.slider({ label: 'Removal λ (1/h)', min: 0.1, max: 30, log: true, value: params.lambda, onChange: (v) => { params.lambda = v; applyGuidance(); } });
    panel.slider({ label: 'Secretion Q (molecules/s per target)', min: 100, max: 20000, log: true, value: params.Q, onChange: (v) => { params.Q = v; applyGuidance(); } });
    lRead = panel.readouts(['decay length ℓ = √(D/λ)']);
    panel.slider({ label: 'Chemotactic steering κ', min: 0, max: 1, value: params.kappa, onChange: (v) => { params.kappa = v; growth.p.kappa = v; } });
    panel.slider({ label: 'Haptotaxis η (laminin gradient, +x)', min: 0, max: 0.4, value: params.hapto, onChange: (v) => { params.hapto = v; growth.p.hapto = v; } });
    panel.slider({ label: 'Turning noise σ', min: 0, max: 0.7, value: params.noise, onChange: (v) => { params.noise = v; growth.p.noise = v; } });
    panel.slider({ label: 'Axon branching (1/h per growth cone)', min: 0, max: 0.6, value: params.branch, onChange: (v) => { params.branch = v; growth.p.branchRate = v; } });
    panel.toggle({ label: 'Show iso-concentration shells (C = K_d, 0.6 K_d)', value: params.showIso, onChange: (v) => { params.showIso = v; isoDirty = true; } });
    panel.buttons([{ label: 'Restart', primary: true, onClick: () => build() }]);
    wireChart = panel.chart({
      title: 'Wiring: neurite length (×100 µm) and synapses', xLabel: 'hours',
      series: [{ name: 'neurite length', color: hex(C_AXON) }, { name: 'synapses', color: hex(C_CONE) }],
    });
    live = panel.readouts(['neurons', 'synapses', 'active growth cones', 'mean SNR at axon growth cones', 'C at 50 µm from a target']);
    applyGuidance();
  }

  function updateLegend() {
    const fate = params.colorBy === 'hes'
      ? [{ color: '#2244dd', label: 'Hes1 low' }, { color: '#ff4040', label: 'Hes1 high (Notch on)' }]
      : params.colorBy === 'brainbow'
        ? [{ color: hex(C_PROG), label: 'progenitor' }, { color: '#ff7ad0', label: 'each neuron and its neurites: one hue' }]
        : S.mode === 'confocal'
          ? [{ color: hex(F_PROG), label: 'Sox2 (progenitor)' }, { color: hex(C_NGN), label: 'Neurogenin-high' }, { color: hex(F_NEU), label: 'βIII-tubulin (neuron)' },
            { color: hex(F_AXON), label: 'axon' }, { color: hex(F_DEND), label: 'dendrite' }, { color: '#4d73ff', label: 'DAPI (nuclei)' }]
          : [{ color: hex(C_PROG), label: 'Sox2/Pax6 progenitor' }, { color: hex(C_NGN), label: 'Neurogenin-high' }, { color: hex(C_NEU), label: 'NeuroD+ neuron' },
            { color: hex(C_AXON), label: 'axon' }, { color: hex(C_DEND), label: 'dendrite' }];
    legend([...fate,
      { color: hex(C_CONE), label: 'growth cone / synaptic bouton' },
      { color: hex(C_TGT), label: 'target cell secreting chemoattractant' },
      { color: '#5fd8ff', label: 'iso-concentration shells' }]);
  }

  const pickables = [
    {
      get object() { return soma.mesh; },
      info: (hit) => {
        const i = hit.instanceId, c = tis.cells[i];
        if (!c) return null;
        const y = Y[i];
        return `${committed[i] ? 'Neuron' : 'Neural progenitor'} ${c.id}\nSox2 ${fmt(y[0])} · Pax6 ${fmt(y[1])} · Hes1 ${fmt(y[6])}\n` +
          `Ngn mRNA ${fmt(y[2])} · Ngn ${fmt(y[3])} · NeuroD ${fmt(y[4])}\nDelta ${fmt(y[5])} · k_tx,eff ${fmt(neuralGRN.ktxEff(y, P))}/h`;
      },
    },
    {
      get object() { return tBody.mesh; },
      info: (hit) => {
        const tg = targets[hit.instanceId];
        if (!tg) return null;
        const p = gp();
        const n = growth.synapses.filter((s) => s.target === tg.id).length;
        return `Target cell ${tg.id + 1}\nsecretes Q = ${fmt(p.Q)} molecules/s\nC at surface ${fmt(guidance.concentration([tg.p[0] + tg.radius, tg.p[1], tg.p[2]], sourceList, p))} nM · ${n} synapses`;
      },
    },
    {
      get object() { return cones.mesh; },
      info: (hit) => {
        const act = growth.tips.filter((tp) => tp.active);
        const tp = act[hit.instanceId];
        if (!tp) return null;
        return `${tp.kind === 0 ? 'Axonal' : 'Dendritic'} growth cone (neuron ${tp.neuron}, order ${tp.order})\nlength ${fmt(tp.len)} µm` +
          (tp.kind === 0 ? `\nC ${fmt(tp.C)} nM · SNR ${fmt(tp.snr)}` : '');
      },
    },
  ];

  build();

  return {
    group,
    buildPanel,
    enter() {
      stage.frame([0, 0, 0], 150, new THREE.Vector3(0.35, 0.42, 1));
      stage.clipAxis.set(0, 0, 1);
      stage.setPickables(pickables);
      updateLegend();
      isoDirty = true;
    },
    update(dt, _t) {
      const dtH = dt * HOURS_PER_SECOND;
      if (dtH > 0) {
        const sub = Math.max(1, Math.ceil(dtH / 0.04));
        for (let k = 0; k < sub; k++) simStep(dtH / sub);
      }
      const now = performance.now();
      if (isoDirty && now - isoLast > 120) { isoLast = now; rebuildIso(); }
      if (t - lastSample >= 0.25 || lastSample < 0) {
        lastSample = t;
        let m = 0, pr = 0, h = 0, s = 0, np = 0, nN = 0, nS = 0;
        for (let i = 0; i < Y.length; i++) {
          if (committed[i]) { nN++; continue; }
          np++; m += Y[i][2]; pr += Y[i][3]; h += Y[i][6]; s += Y[i][0];
          if (Y[i][0] > 0.5) nS++;
        }
        np = Math.max(np, 1);
        netChart.push(t, [m / np, pr / np, h / np, s / np]);
        fracChart.push(t, [nN / Y.length, nS / Y.length]);
        wireChart.push(t, [growth.totalLength / 100, growth.synapses.length]);
        const act = growth.tips.filter((tp) => tp.active && tp.kind === 0);
        live.set('neurons', `${nN} / ${Y.length}`);
        live.set('synapses', growth.synapses.length);
        live.set('active growth cones', act.length);
        live.set('mean SNR at axon growth cones', act.length ? fmt(act.reduce((a, tp) => a + tp.snr, 0) / act.length) : '–');
        const tg = targets[0].p, p = gp();
        live.set('C at 50 µm from a target', `${fmt(guidance.concentration([tg[0] + 50, tg[1], tg[2]], sourceList, p))} nM`);
      }
      draw(performance.now() / 1000);
      let nN = 0;
      for (let i = 0; i < committed.length; i++) if (committed[i]) nN++;
      setStatus(`t = ${t.toFixed(1)} h · ${nN}/${Y.length} neurons · ${growth.synapses.length} synapses · neurites ${(growth.totalLength / 1000).toFixed(2)} mm`);
    },
    reset() { build(); },
    onMode(mode, active) { isoDirty = true; redrawSegments(); if (active) updateLegend(); },
  };
}

const _ax = [1, 0, 0], _ay = [0, 1, 0];
function perp(d, out) {
  const c = cross(d, Math.abs(d[0]) < 0.9 ? _ax : _ay, out);
  const n = Math.hypot(c[0], c[1], c[2]) || 1;
  c[0] /= n; c[1] /= n; c[2] /= n;
  return c;
}
function cross(a, b, out = [0, 0, 0]) {
  const x = a[1] * b[2] - a[2] * b[1], y = a[2] * b[0] - a[0] * b[2], z = a[0] * b[1] - a[1] * b[0];
  out[0] = x; out[1] = y; out[2] = z;
  return out;
}

// =====================================================================
// (c) Hormonal feedback
// =====================================================================
function makeHormone(S) {
  const { stage, root, legend, setStatus, G } = S;
  const MIN_PER_SECOND = 4;
  const group = new THREE.Group();
  root.add(group);
  const params = { n: 12, stim: 1, rate: 0.1, pulse: 3 };
  const p = () => ({ ...goodwin.defaults, n: params.n, stim: stimNow, b: params.rate, d: params.rate, f: params.rate });
  let stimNow = 1, pulseUntil = -1;

  // ------------------------------------------------------------ anatomy (schematic)
  const GL = {
    H: { name: 'hypothalamus-like (X)', c: [-120, 62, 0], n: 26, spread: [26, 16, 22], color: 0xa98bff },
    P: { name: 'pituitary-like (Y)', c: [-40, -8, 12], n: 30, spread: [20, 20, 20], color: 0x4fd3ff },
    T: { name: 'target gland (Z)', c: [88, -26, 0], n: 40, spread: [30, 20, 24], color: 0xffb347 },
    E: { name: 'target tissue (Z receptors)', c: [42, 92, -18], n: 24, spread: [28, 14, 20], color: 0xff7aa8 },
  };
  const ctrl = [
    GL.H.c, [-86, 30, 14], GL.P.c, [20, -40, 20], GL.T.c, [96, 40, 4], GL.E.c, [-40, 110, -12],
  ].map((v) => new THREE.Vector3(...v));
  const curve = new THREE.CatmullRomCurve3(ctrl, true, 'centripetal');
  const NS = 1600;
  const pts = curve.getSpacedPoints(NS);
  const frames = curve.computeFrenetFrames(NS, true);
  const uOf = (c) => {
    let best = 0, bd = Infinity;
    for (let i = 0; i < NS; i++) { const d = pts[i].distanceToSquared(new THREE.Vector3(...c)); if (d < bd) { bd = d; best = i; } }
    return best / NS;
  };
  const uH = uOf(GL.H.c), uP = uOf(GL.P.c), uT = uOf(GL.T.c), uE = uOf(GL.E.c);
  const span = (a, b) => ((b - a) % 1 + 1) % 1;
  const routes = [
    { k: 0, u0: uH, du: span(uH, uP), color: new THREE.Color(0xb48cff) },
    { k: 1, u0: uP, du: span(uP, uT), color: new THREE.Color(0x4fd3ff) },
    { k: 2, u0: uT, du: span(uT, uH), color: new THREE.Color(0xffb347) },
  ];
  const VR = 4.2;
  const vessel = new THREE.Mesh(new THREE.TubeGeometry(curve, 400, VR, 18, true),
    cellMaterial({ role: 'membrane', color: 0xff5a6e, stain: 0xff3b6b, opacity: 0.1, he: [0.02, 0.3] }));
  vessel.renderOrder = 2;
  group.add(vessel);

  const rngA = new RNG(515);
  const cells = [];
  const _q = new THREE.Vector3();
  for (const [key, gl] of Object.entries(GL)) {
    const list = [];
    let tries = 0;
    while (list.length < gl.n && tries < 5000) {
      tries++;
      const d = rngA.onSphere(), r = Math.cbrt(rngA.next());
      const q = [gl.c[0] + d[0] * r * gl.spread[0], gl.c[1] + d[1] * r * gl.spread[1], gl.c[2] + d[2] * r * gl.spread[2]];
      // keep the vessel lumen free
      let near = Infinity;
      _q.set(q[0], q[1], q[2]);
      for (let i = 0; i < NS; i += 8) near = Math.min(near, pts[i].distanceTo(_q));
      if (near < VR + 4.5) continue;
      if (list.some((o) => Math.hypot(o.p[0] - q[0], o.p[1] - q[1], o.p[2] - q[2]) < 10.5)) continue;
      list.push({ gland: key, p: q, R: 5.6 + 0.8 * rngA.next(), receptors: [] });
    }
    for (const c of list) {
      const nr = key === 'H' ? 6 : key === 'E' ? 9 : 0;
      for (let k = 0; k < nr; k++) c.receptors.push({ d: rngA.onSphere(), thr: rngA.next() });
      cells.push(c);
    }
  }

  const cellMat = cellMaterial({ role: 'membrane', opacity: 0.3, useInst: 1, gain: 0.55 });
  const nucMat = cellMaterial({ role: 'nucleus', color: 0xffffff, useInst: 1, gain: 0.45 });
  const partMat = cellMaterial({ role: 'molecule', useInst: 1, emissive: 0.35, he: [0.3, 0.6], gain: 0.8 });
  const recMat = cellMaterial({ role: 'reporter', useInst: 1, emissive: 0.25, he: [0.25, 0.8], gain: 0.7 });
  const cellPool = new InstancePool(G.sphere, cellMat, 128, group);
  const nucPool = new InstancePool(G.sphereMid, nucMat, 128, group);
  const partPool = new InstancePool(G.sphereLo, partMat, 512, group);
  const recPool = new InstancePool(G.sphereLo, recMat, 512, group);
  cellPool.mesh.renderOrder = 3;

  // phase portrait (X/X*, Y/Y*, Z/Z*) above the axis
  const PP = { o: new THREE.Vector3(-40, 2, -28), s: 18 };
  const axesGeo = new THREE.BufferGeometry();
  axesGeo.setAttribute('position', new THREE.Float32BufferAttribute([
    0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3,
  ].map((v) => v * PP.s), 3));
  axesGeo.setAttribute('color', new THREE.Float32BufferAttribute([
    ...routes[0].color.toArray(), ...routes[0].color.toArray(), ...routes[1].color.toArray(), ...routes[1].color.toArray(),
    ...routes[2].color.toArray(), ...routes[2].color.toArray()], 3));
  const axes = new THREE.LineSegments(axesGeo, lineMaterial({ color: 0xffffff, vertexColors: true, opacity: 0.9 }));
  axes.position.copy(PP.o);
  group.add(axes);
  const TRAJ = 1500;
  const trajPos = new Float32Array(TRAJ * 3);
  const trajGeo = new THREE.BufferGeometry();
  trajGeo.setAttribute('position', new THREE.BufferAttribute(trajPos, 3).setUsage(THREE.DynamicDrawUsage));
  trajGeo.setDrawRange(0, 0);
  const traj = new THREE.Line(trajGeo, lineMaterial({ color: 0xfff1c9, opacity: 0.85 }));
  traj.frustumCulled = false;
  traj.position.copy(PP.o);
  group.add(traj);
  const markMat = cellMaterial({ role: 'reporter', useInst: 1, emissive: 0.6 });
  const marks = new InstancePool(G.sphereMid, markMat, 4, group);
  let trajN = 0;

  // ------------------------------------------------------------ state
  let y, tMin, lastSample, work, parts, ss, hist;
  function build() {
    const pp = p();
    ss = goodwin.steadyState(pp);
    // start away from the steady state (as after a perturbation)
    y = Float64Array.from([ss[0] * 0.2, ss[1] * 0.3, ss[2] * 0.4, ss[3] * 0.4]);
    tMin = 0; lastSample = -1; pulseUntil = -1;
    work = rk4Work(4);
    parts = routes.map(() => []);
    hist = { t: [], z: [] };
    trajN = 0; trajGeo.setDrawRange(0, 0);
    levelChart?.clear(); occChart?.clear();
  }

  // ------------------------------------------------------------ panel
  let levelChart, occChart, stab, live;
  function showStability() {
    const pp = p();
    const st = goodwin.stability(pp);
    ss = goodwin.steadyState(pp);
    stab?.set('loop gain G = n·u/(1+u)', fmt(st.gain));
    stab?.set('Hopf threshold G_c', fmt(st.gainCrit));
    stab?.set('critical Hill n', fmt(goodwin.criticalN(pp)));
    stab?.set('linear prediction', st.stable ? 'damped (stable steady state)' : `limit cycle (T ≈ ${fmt(st.periodHopf)} min at onset)`);
  }
  function buildPanel(panel) {
    panel.section('Goodwin negative-feedback loop');
    panel.equation('dX/dt = s·a / (1 + (Z/K)ⁿ) − b·X      (hypothalamus)');
    panel.equation('dY/dt = c·X − d·Y                    (pituitary)');
    panel.equation('dZ/dt = e·Y − f·Z                    (target gland)');
    panel.equation('dO/dt = k_on·Z·(1 − O) − k_off·O     (target receptors)');
    panel.equation('Hopf: G = n u/(1+u) > (b+d+f)(bd+bf+df)/(bdf) − 1  (= 8 if b = d = f)');
    const nS = panel.slider({ label: 'Feedback cooperativity (Hill n)', min: 1, max: 20, step: 0.5, value: params.n, onChange: (v) => { params.n = v; showStability(); } });
    const setN = (v) => { params.n = v; nS.set(v); showStability(); };
    panel.slider({ label: 'Stimulus s (secretion drive)', min: 0.2, max: 3, value: params.stim, onChange: (v) => { params.stim = v; showStability(); } });
    panel.slider({ label: 'Clearance rate b = d = f (1/min)', min: 0.03, max: 0.3, value: params.rate, onChange: (v) => { params.rate = v; showStability(); } });
    panel.buttons([
      { label: 'Acute stimulus (×3, 15 min)', primary: true, onClick: () => { pulseUntil = tMin + 15; } },
      { label: 'n = 4 (damped)', onClick: () => setN(4) },
      { label: 'n = 12 (oscillating)', onClick: () => setN(12) },
    ]);
    stab = panel.readouts(['loop gain G = n·u/(1+u)', 'Hopf threshold G_c', 'critical Hill n', 'linear prediction']);
    levelChart = panel.chart({
      title: 'Hormone levels relative to steady state (X/X*, Y/Y*, Z/Z*)', xLabel: 'minutes', capacity: 600,
      series: [{ name: 'X releasing', color: hex(routes[0].color) }, { name: 'Y tropic', color: hex(routes[1].color) }, { name: 'Z end', color: hex(routes[2].color) }],
    });
    occChart = panel.chart({
      title: 'Receptor occupancy', xLabel: 'minutes', yRange: [0, 1], capacity: 600,
      series: [{ name: 'feedback receptors (hypothalamus)', color: hex(GL.H.color) }, { name: 'target-cell receptors O', color: hex(GL.E.color) }],
    });
    live = panel.readouts(['measured behaviour', 'period (peak to peak)', 'relative amplitude of Z']);
    panel.note('Receptor occupancy by mass action (K_d = k_off/k_on). The feedback term 1/(1 + (Z/K)ⁿ) is the unoccupied fraction of a cooperative feedback receptor. Rates are dimensionless and illustrative (minutes), not fitted.');
    showStability();
  }

  function updateLegend() {
    legend([
      { color: hex(GL.H.color), label: 'hypothalamus-like cells: secrete X' },
      { color: hex(GL.P.color), label: 'pituitary-like cells: secrete Y' },
      { color: hex(GL.T.color), label: 'target gland: secretes Z' },
      { color: hex(GL.E.color), label: 'target tissue with Z receptors' },
      { color: '#ffd35a', label: 'receptor occupied (lit) / free (dim)' },
      { color: '#fff1c9', label: 'phase portrait (X, Y, Z)/steady state' },
    ]);
  }

  // ------------------------------------------------------------ drawing
  const tmp = new THREE.Color(), tmp2 = new THREE.Color();
  const glandCol = Object.fromEntries(Object.entries(GL).map(([k, g]) => [k, new THREE.Color(g.color)]));
  const REC_OFF = new THREE.Color(0x3a4658), REC_ON = new THREE.Color(0xffd35a);
  function lut(u, out) {
    const f = (((u % 1) + 1) % 1) * NS;
    const i = Math.floor(f), a = f - i;
    const p0 = pts[i], p1 = pts[(i + 1) % NS];
    out[0] = p0.x + (p1.x - p0.x) * a; out[1] = p0.y + (p1.y - p0.y) * a; out[2] = p0.z + (p1.z - p0.z) * a;
    return i;
  }
  const P3 = [0, 0, 0];
  function draw(pp, time) {
    const Zk = Math.pow(Math.max(y[2], 0) / pp.K, pp.n);
    const secX = (stimNow * pp.a) / (1 + Zk), secX0 = pp.b * ss[0];
    const rate = { H: secX / secX0, P: y[0] / ss[0], T: y[1] / ss[1], E: 1 };
    const fbOcc = goodwin.feedbackOccupancy(y[2], pp);
    cellPool.begin(); nucPool.begin(); recPool.begin(); partPool.begin(); marks.begin();
    for (const c of cells) {
      const act = clamp(rate[c.gland], 0, 2.5);
      const col = tmp.copy(glandCol[c.gland]).multiplyScalar(0.45 + 0.4 * act);
      cellPool.put(c.p[0], c.p[1], c.p[2], c.R, c.R * 0.92, c.R, col);
      nucPool.put(c.p[0], c.p[1], c.p[2], c.R * 0.55, c.R * 0.5, c.R * 0.55,
        S.mode === 'confocal' ? tmp2.setRGB(0.3, 0.45, 1) : S.mode === 'histology' ? tmp2.setRGB(1, 1, 1) : tmp2.setRGB(0.88, 0.9, 1));
      const occ = c.gland === 'H' ? fbOcc : y[3];
      for (const r of c.receptors) {
        const on = smooth(r.thr - 0.06, r.thr + 0.06, occ);
        const rr = c.R * 1.03;
        const s = 0.9 + 0.5 * on;
        recPool.put(c.p[0] + r.d[0] * rr, c.p[1] + r.d[1] * rr * 0.92, c.p[2] + r.d[2] * rr, s, s, s, tmp2.copy(REC_OFF).lerp(REC_ON, on));
      }
    }
    for (const rt of routes) {
      for (const q of parts[rt.k]) {
        const i = lut(rt.u0 + q.s * rt.du, P3);
        const N = frames.normals[i], B = frames.binormals[i];
        const ox = Math.cos(q.a) * q.r, oy = Math.sin(q.a) * q.r;
        partPool.put(P3[0] + N.x * ox + B.x * oy, P3[1] + N.y * ox + B.y * oy, P3[2] + N.z * ox + B.z * oy, 1.5, 1.5, 1.5, rt.color);
      }
    }
    // current state + steady state on the phase portrait
    const cx = (y[0] / ss[0]) * PP.s, cy = (y[1] / ss[1]) * PP.s, cz = (y[2] / ss[2]) * PP.s;
    marks.put(PP.o.x + Math.min(cx, 3.2 * PP.s), PP.o.y + Math.min(cy, 3.2 * PP.s), PP.o.z + Math.min(cz, 3.2 * PP.s), 2.4, 2.4, 2.4, tmp.setRGB(1, 0.95, 0.75));
    marks.put(PP.o.x + PP.s, PP.o.y + PP.s, PP.o.z + PP.s, 1.6, 1.6, 1.6, tmp.setRGB(0.6, 0.65, 0.7));
    cellPool.end(); nucPool.end(); recPool.end(); partPool.end(); marks.end();
  }

  // particle counts proportional to the ODE state (visual proxy)
  const NMAX = 150;
  function updateParticles(dtS) {
    for (const rt of routes) {
      const arr = parts[rt.k];
      const target = Math.round(NMAX * clamp(y[rt.k] / (2.2 * ss[rt.k]), 0, 1));
      const speed = 1 / (3.2 + 2.5 * rt.du); // routes traversed in a few seconds
      for (const q of arr) { q.s += dtS * speed * q.v; if (q.s > 1) q.s -= 1; }
      let add = Math.min(target - arr.length, 6);
      while (add-- > 0) arr.push({ s: rngA.next() * 0.04, a: rngA.next() * TAU, r: Math.sqrt(rngA.next()) * (VR - 1.6), v: 0.8 + 0.4 * rngA.next() });
      let remove = Math.min(arr.length - target, 6);
      while (remove-- > 0) arr.splice(rngA.int(arr.length), 1);
    }
  }

  const pickables = [
    {
      get object() { return cellPool.mesh; },
      info: (hit) => {
        const c = cells[hit.instanceId];
        if (!c) return null;
        const pp = p();
        const g = GL[c.gland];
        let extra = '';
        if (c.gland === 'H') extra = `\nX secretion ${fmt((stimNow * pp.a) / (1 + Math.pow(y[2] / pp.K, pp.n)))} /min · feedback occupancy ${(100 * goodwin.feedbackOccupancy(y[2], pp)).toFixed(0)} %`;
        if (c.gland === 'P') extra = `\nY secretion c·X = ${fmt(pp.c * y[0])} /min`;
        if (c.gland === 'T') extra = `\nZ secretion e·Y = ${fmt(pp.e * y[1])} /min`;
        if (c.gland === 'E') extra = `\nZ-receptor occupancy O = ${(100 * y[3]).toFixed(0)} %`;
        return `${g.name}${extra}\nX ${fmt(y[0])} · Y ${fmt(y[1])} · Z ${fmt(y[2])}`;
      },
    },
  ];

  build();

  return {
    group,
    buildPanel,
    enter() {
      stage.frame([-12, 30, 0], 138, new THREE.Vector3(0.18, 0.3, 1));
      stage.clipAxis.set(0, 0, 1);
      stage.setPickables(pickables);
      updateLegend();
    },
    update(dt) {
      const dtM = dt * MIN_PER_SECOND;
      if (dtM > 0) {
        const sub = Math.max(1, Math.ceil(dtM / 0.05));
        const h = dtM / sub;
        for (let k = 0; k < sub; k++) {
          stimNow = params.stim * (tMin < pulseUntil ? params.pulse : 1);
          rk4Step(goodwin.rhs(p()), tMin, y, h, work);
          tMin += h;
        }
        stimNow = params.stim * (tMin < pulseUntil ? params.pulse : 1);
      }
      const pp = p();
      ss = goodwin.steadyState({ ...pp, stim: params.stim });
      updateParticles(dt);
      if (tMin - lastSample >= 0.5 || lastSample < 0) {
        lastSample = tMin;
        levelChart.push(tMin, [y[0] / ss[0], y[1] / ss[1], y[2] / ss[2]]);
        occChart.push(tMin, [goodwin.feedbackOccupancy(y[2], pp), y[3]]);
        hist.t.push(tMin); hist.z.push(y[2]);
        if (hist.t.length > 800) { hist.t.shift(); hist.z.shift(); }
        // phase-portrait trajectory
        if (trajN >= TRAJ) { trajPos.copyWithin(0, 3); trajN = TRAJ - 1; }
        trajPos[trajN * 3] = Math.min(y[0] / ss[0], 3.2) * PP.s;
        trajPos[trajN * 3 + 1] = Math.min(y[1] / ss[1], 3.2) * PP.s;
        trajPos[trajN * 3 + 2] = Math.min(y[2] / ss[2], 3.2) * PP.s;
        trajN++;
        trajGeo.setDrawRange(0, trajN);
        trajGeo.attributes.position.needsUpdate = true;
        if (hist.t.length > 40) {
          const o = oscillationStats(hist.t, hist.z, { from: tMin - 200 });
          live.set('measured behaviour', tMin < 150 ? 'transient…' : o.sustained ? 'sustained oscillation' : o.relAmplitude < 0.01 ? 'at steady state (homeostasis)' : 'damped oscillation');
          live.set('period (peak to peak)', Number.isFinite(o.period) ? `${fmt(o.period)} min` : '–');
          live.set('relative amplitude of Z', fmt(o.relAmplitude));
        }
      }
      draw(pp, performance.now() / 1000);
      const st = goodwin.stability({ ...pp, stim: params.stim });
      setStatus(`t = ${tMin.toFixed(0)} min · n = ${params.n} · ${st.stable ? `damped (G = ${fmt(st.gain)} < ${fmt(st.gainCrit)})` : `oscillatory (G = ${fmt(st.gain)} > ${fmt(st.gainCrit)})`}${tMin < pulseUntil ? ' · stimulus pulse' : ''}`);
    },
    reset() { build(); },
    onMode() {},
  };
}
