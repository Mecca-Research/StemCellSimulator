// Waddington's epigenetic landscape computed from a gene-circuit model
// (U = -ln P), with DNMT / HAT / HDAC chromatin marks from the paper.
import * as THREE from 'three';
import { cellMaterial } from '../engine/materials.js';
import { unitSphere } from '../engine/geometry.js';
import { InstancePool } from '../engine/instances.js';
import { RNG } from '../models/core/rng.js';
import { rk4Step, rk4Work } from '../models/core/ode.js';
import { marks, circuit, landscape, sampleU } from '../models/pathways/epigenetics.js';

const MAXX = 2.4;

function viridis(t, out) {
  // compact viridis approximation
  const c = [[0.267, 0.005, 0.329], [0.229, 0.322, 0.546], [0.128, 0.567, 0.551], [0.369, 0.789, 0.383], [0.993, 0.906, 0.144]];
  const x = Math.min(Math.max(t, 0), 1) * 4, i = Math.min(Math.floor(x), 3), f = x - i;
  return out.setRGB(c[i][0] + (c[i + 1][0] - c[i][0]) * f, c[i][1] + (c[i + 1][1] - c[i][1]) * f, c[i][2] + (c[i + 1][2] - c[i][2]) * f);
}

export default {
  about: `
    <p><b>Waddington's landscape, computed.</b> Two master regulators that activate
    themselves and repress each other (the GATA1–PU.1 motif; Huang et al. 2007) define
    the fate dynamics. Stochastic trajectories give the steady-state probability P, and
    the landscape is its <b>quasi-potential U = −ln P</b> (Wang, Xu &amp; Wang 2008):
    valleys are fates, ridges are barriers.</p>
    <p>As development proceeds, self-activation weakens and the single multipotent valley
    splits into two committed valleys — a bifurcation the cells (marbles) roll through.
    The paper's chromatin reactions (DNA + DNMT → methylated DNA; histone + HAT ⇄ HDAC)
    set the marks; their coupling to the landscape (acetylation → plasticity/noise,
    methylation → deeper canalisation) is <b>phenomenological</b>.</p>`,
  paperRef: 'Paper: "Epigenetic Modifications: DNA Methylation and Histone Modification" and "Gene Expression Regulation: controlled by transcription factors, signaling molecules, and epigenetic modifications".',

  async create(ctx) {
    const { stage, root, panel, legend, setStatus, query } = ctx;
    const rng = new RNG(99);
    stage.setBloomScale(0.35);
    const mp = { ...marks.defaults };
    let chrom = marks.steadyState(mp); // [methylation, acetylation]
    const params = { view: query?.get('view') === 'state' ? 'state' : 'waddington', devSpeed: 0.06, cells: 160, tau: 0.5, noiseMul: 1 };

    // ------------------------------------------------ Waddington surface (tau x fate)
    const NT = 28, ND = 90;
    const W = 120, Ddepth = 150;
    const surfGeo = new THREE.PlaneGeometry(W, Ddepth, ND - 1, NT - 1);
    surfGeo.rotateX(-Math.PI / 2);
    surfGeo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(ND * NT * 3), 3));
    const surf = new THREE.Mesh(surfGeo, cellMaterial({ role: 'field', vertexColors: true, noise: 0, side: THREE.DoubleSide }));
    root.add(surf);
    let U2 = null; // Float32Array NT x ND rows
    const uScale = 7, drop = 70;

    function computeWaddington() {
      const hist = new Float64Array(ND);
      U2 = new Float32Array(NT * ND);
      const out = [0, 0];
      for (let r = 0; r < NT; r++) {
        const tau = r / (NT - 1);
        const p = circuit.params(tau, chrom[0], chrom[1]);
        p.sigma *= params.noiseMul;
        hist.fill(0);
        const lr = new RNG(1000 + r);
        for (let w = 0; w < 260; w++) {
          const s = [0.5 + lr.uniform(-0.2, 0.2), 0.5 + lr.uniform(-0.2, 0.2)];
          for (let k = 0; k < 520; k++) {
            circuit.step(p, s, 0.03, lr, out);
            if (k < 80) continue;
            const d = (s[0] - s[1]) / MAXX; // -1..1
            hist[Math.min(ND - 1, Math.max(0, Math.floor(((d + 1) / 2) * ND)))] += 1;
          }
        }
        // smooth and convert to -ln P
        let tot = 0; const sm = new Float64Array(ND);
        for (let i = 0; i < ND; i++) { let v = 0, wsum = 0; for (let k = -2; k <= 2; k++) { const j = i + k; if (j >= 0 && j < ND) { const wk = 3 - Math.abs(k); v += wk * hist[j]; wsum += wk; } } sm[i] = v / wsum; tot += sm[i]; }
        let lo = Infinity;
        for (let i = 0; i < ND; i++) { const u = -Math.log(sm[i] / tot + 0.3 / (ND * 20)); U2[r * ND + i] = u; lo = Math.min(lo, u); }
        for (let i = 0; i < ND; i++) U2[r * ND + i] = Math.min(U2[r * ND + i] - lo, 6);
      }
      const pos = surfGeo.attributes.position, col = surfGeo.attributes.color, c = new THREE.Color();
      for (let r = 0; r < NT; r++) for (let i = 0; i < ND; i++) {
        const k = r * ND + i;
        const u = U2[k];
        // rows run from far (tau = 0, high) to near (tau = 1, low)
        pos.setY(k, u * uScale - (r / (NT - 1)) * drop + drop / 2);
        viridis(1 - u / 6, c);
        col.setXYZ(k, c.r, c.g, c.b);
      }
      pos.needsUpdate = true; col.needsUpdate = true;
      surfGeo.computeVertexNormals();
    }
    function heightAt(tau, d) {
      const fr = Math.min(Math.max(tau, 0), 1) * (NT - 1), fi = Math.min(Math.max((d / MAXX + 1) / 2, 0), 1) * (ND - 1);
      const r = Math.min(Math.floor(fr), NT - 2), i = Math.min(Math.floor(fi), ND - 2), u = fr - r, v = fi - i;
      const g = (rr, ii) => U2[rr * ND + ii];
      const U = (g(r, i) * (1 - v) + g(r, i + 1) * v) * (1 - u) + (g(r + 1, i) * (1 - v) + g(r + 1, i + 1) * v) * u;
      return U * uScale - Math.min(Math.max(tau, 0), 1) * drop + drop / 2;
    }

    // ------------------------------------------------ state-space surface (x, y) at tau
    const NS = 64;
    const ssGeo = new THREE.PlaneGeometry(W, W, NS - 1, NS - 1);
    ssGeo.rotateX(-Math.PI / 2);
    ssGeo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(NS * NS * 3), 3));
    const ss = new THREE.Mesh(ssGeo, cellMaterial({ role: 'field', vertexColors: true, noise: 0, side: THREE.DoubleSide }));
    root.add(ss);
    let L = null;
    function computeState() {
      const p = circuit.params(params.tau, chrom[0], chrom[1]);
      p.sigma *= params.noiseMul;
      L = landscape(p, { N: NS, max: MAXX, walkers: 400, steps: 700, burn: 120, dt: 0.02, seed: 5 });
      const pos = ssGeo.attributes.position, col = ssGeo.attributes.color, c = new THREE.Color();
      for (let j = 0; j < NS; j++) for (let i = 0; i < NS; i++) {
        // plane vertex order: rows from -z (far) to +z; map row j -> y index (NS-1-j)
        const k = j * NS + i;
        const u = Math.min(L.U[(NS - 1 - j) * NS + i], 7);
        pos.setY(k, u * 5);
        viridis(1 - u / 7, c);
        col.setXYZ(k, c.r, c.g, c.b);
      }
      pos.needsUpdate = true; col.needsUpdate = true;
      ssGeo.computeVertexNormals();
    }
    const ssHeight = (x, y) => Math.min(sampleU(L, x, y), 7) * 5;

    // ------------------------------------------------ cells (marbles)
    const pool = new InstancePool(unitSphere(3), cellMaterial({ role: 'reporter', color: 0xffffff }), 512, root);
    const cells = [];
    function spawn(n) {
      cells.length = 0;
      for (let i = 0; i < n; i++) cells.push({ x: 0.55 + rng.normal() * 0.05, y: 0.55 + rng.normal() * 0.05, tau: rng.uniform(0, 0.15), trail: [] });
    }
    spawn(params.cells);
    const colA = new THREE.Color(0xff5a6e), colB = new THREE.Color(0x5aa8ff), colM = new THREE.Color(0xf2f2f2), tmp = new THREE.Color();
    const fateOf = (c) => (c.x - c.y > 0.45 ? 'A' : c.y - c.x > 0.45 ? 'B' : 'M');

    // ------------------------------------------------ controls
    panel.section('View');
    panel.select({
      label: 'Landscape', value: params.view,
      options: [{ value: 'waddington', label: 'Waddington: developmental time × fate' }, { value: 'state', label: 'State space: GATA1 × PU.1 at a fixed stage' }],
      onChange: (v) => { params.view = v; layout(); },
    });
    const tauSlider = panel.slider({ label: 'Developmental stage τ (state-space view)', min: 0, max: 1, value: params.tau, onChange: (v) => { params.tau = v; scheduleRecompute(); } });
    panel.section('Chromatin (paper reactions)');
    panel.slider({ label: 'DNMT (DNA methyltransferase)', min: 0, max: 1.5, value: mp.DNMT, onChange: (v) => { mp.DNMT = v; } });
    panel.slider({ label: 'HAT (histone acetyltransferase)', min: 0, max: 1.5, value: mp.HAT, onChange: (v) => { mp.HAT = v; } });
    panel.slider({ label: 'HDAC (histone deacetylase)', min: 0, max: 1.5, value: mp.HDAC, onChange: (v) => { mp.HDAC = v; } });
    panel.equation('dm/dt = k_M·DNMT·(1−m) − k_DM·m\nda/dt = k_A·HAT·(1−a) − k_D·HDAC·a\nnoise σ = 0.04 + 0.22·a     (plasticity)\nrepression b ∝ 0.7 + 0.6·m (canalisation)');
    const marksOut = panel.readouts(['methylation m', 'acetylation a']);
    panel.buttons([
      { label: 'Recompute landscape', onClick: () => recompute() },
      { label: 'Reprogram (OSKM-like reset)', primary: true, onClick: () => { for (const c of cells) { c.tau = 0; c.x = 0.55 + rng.normal() * 0.05; c.y = 0.55 + rng.normal() * 0.05; } } },
    ]);
    panel.note('Reprogramming sends every cell back to the top valley — the landscape analogue of making induced pluripotent stem cells (like the hiPSCs in the colony scene).');
    panel.section('Cells');
    panel.slider({ label: 'Speed of development', min: 0, max: 0.2, value: params.devSpeed, onChange: (v) => { params.devSpeed = v; } });
    panel.slider({ label: 'Noise multiplier', min: 0.2, max: 3, value: 1, onChange: (v) => { params.noiseMul = v; scheduleRecompute(); } });
    const fateChart = panel.chart({ title: 'Fate of the cell population', yRange: [0, 1], xLabel: 't', series: [
      { name: 'multipotent', color: '#f2f2f2' }, { name: 'fate A (GATA1-high)', color: '#ff5a6e' }, { name: 'fate B (PU.1-high)', color: '#5aa8ff' }] });
    const bif = panel.chart({ title: 'Bifurcation: attractors (GATA1 − PU.1) vs stage τ', xLabel: 'τ', series: [
      { name: 'fate A', color: '#ff5a6e', points: true, width: 0.01 }, { name: 'multipotent', color: '#f2f2f2', points: true, width: 0.01 }, { name: 'fate B', color: '#5aa8ff', points: true, width: 0.01 }] });
    function drawBifurcation() {
      const xs = [], a = [], m = [], b = [];
      for (let i = 0; i <= 24; i++) {
        const tau = i / 24;
        const att = circuit.attractors(circuit.params(tau, chrom[0], chrom[1]), { grid: 6, T: 40 });
        xs.push(tau);
        const ds = att.map(([x, y]) => x - y);
        a.push(Math.max(...ds)); b.push(Math.min(...ds));
        const mid = ds.find((d) => Math.abs(d) < 0.2);
        m.push(mid === undefined ? NaN : mid);
      }
      bif.set(xs, [a, m, b]);
    }

    legend([{ color: '#f2f2f2', label: 'uncommitted cell' }, { color: '#ff5a6e', label: 'committed to fate A (GATA1-high, erythroid)' }, { color: '#5aa8ff', label: 'committed to fate B (PU.1-high, myeloid)' }, { color: '#fde725', label: 'low U (valley)' }, { color: '#440154', label: 'high U (ridge)' }]);

    let timer = null;
    function scheduleRecompute() { clearTimeout(timer); timer = setTimeout(recompute, 250); }
    function recompute() {
      if (params.view === 'waddington') computeWaddington(); else computeState();
      drawBifurcation();
    }
    function layout() {
      const wad = params.view === 'waddington';
      surf.visible = wad; ss.visible = !wad;
      tauSlider.input.disabled = wad;
      recompute();
      // classic Waddington view: stand at the low (committed) end looking up-slope
      if (wad) stage.frame([0, 0, -10], 105, new THREE.Vector3(0.45, 0.85, 1.0));
      else stage.frame([0, 10, 0], 90, new THREE.Vector3(0.6, 1.0, 0.9));
      stage.clipAxis.set(1, 0, 0);
    }
    layout();

    const work = rk4Work(2);
    const drift = [0, 0];
    let t = 0, lastPush = -1, lastMarks = chrom.slice();
    return {
      update(dt) {
        t += dt;
        // chromatin marks relax toward the slider-defined steady state
        const y = Float64Array.from(chrom);
        rk4Step(marks.rhs(mp), t, y, dt, work);
        chrom = [y[0], y[1]];
        marksOut.set('methylation m', chrom[0]);
        marksOut.set('acetylation a', chrom[1]);
        if (Math.abs(chrom[0] - lastMarks[0]) + Math.abs(chrom[1] - lastMarks[1]) > 0.04) { lastMarks = chrom.slice(); scheduleRecompute(); }

        const wad = params.view === 'waddington';
        pool.begin();
        const counts = { A: 0, B: 0, M: 0 };
        for (const c of cells) {
          const tau = wad ? c.tau : params.tau;
          const p = circuit.params(tau, chrom[0], chrom[1]);
          p.sigma *= params.noiseMul;
          const w = [c.x, c.y];
          for (let k = 0; k < 4; k++) circuit.step(p, w, dt * 0.9, rng, drift);
          c.x = Math.min(w[0], MAXX); c.y = Math.min(w[1], MAXX);
          if (wad) c.tau = Math.min(c.tau + dt * params.devSpeed, 1);
          const f = fateOf(c); counts[f]++;
          const col = f === 'A' ? colA : f === 'B' ? colB : colM;
          if (wad) {
            const d = c.x - c.y;
            const px = (d / MAXX) * (W / 2), pz = -Ddepth / 2 + c.tau * Ddepth;
            pool.put(px, heightAt(c.tau, d) + 1.6, pz, 1.6, 1.6, 1.6, tmp.copy(col));
          } else if (L) {
            const px = (c.x / MAXX - 0.5) * W, pz = (0.5 - c.y / MAXX) * W;
            pool.put(px, ssHeight(c.x, c.y) + 1.4, pz, 1.4, 1.4, 1.4, tmp.copy(col));
          }
        }
        pool.end();
        if (t - lastPush > 0.25) {
          lastPush = t;
          const n = cells.length;
          fateChart.push(t, [counts.M / n, counts.A / n, counts.B / n]);
        }
        setStatus(`t = ${t.toFixed(1)} · ${counts.M} uncommitted · ${counts.A} fate A · ${counts.B} fate B`);
      },
      reset() { spawn(params.cells); fateChart.clear(); t = 0; },
      dispose() { clearTimeout(timer); },
    };
  },
};
