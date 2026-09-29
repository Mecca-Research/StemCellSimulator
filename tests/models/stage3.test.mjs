import test from 'node:test';
import assert from 'node:assert/strict';
import { integrate } from '../../js/models/core/ode.js';
import { RNG } from '../../js/models/core/rng.js';
import { Tissue, Cell } from '../../js/models/morpho/agents3d.js';
import { mrf, MyoCulture, segSeg2D, myoblastVolume, REACTIONS as MYO } from '../../js/models/pathways/myogenesis.js';
import { paperNgn, neuralGRN, guidance, NeuriteGrowth, REACTIONS as NEURO } from '../../js/models/pathways/neurogenesis.js';
import { goodwin, oscillationStats, REACTIONS as HORM } from '../../js/models/pathways/hormone.js';

// ------------------------------------------------------------ neurogenin ODE
test("paper's Neurogenin ODE matches its analytic solution (< 1e-6 relative)", () => {
  for (const p of [paperNgn.defaults, { ktx: 0.7, kdm: 0.9, ktl: 3, kdp: 0.9 }]) { // includes k_dm = k_dp
    for (const [m0, P0] of [[0, 0], [2.5, 0.3]]) {
      const { t, y } = integrate(paperNgn.rhs(p), [m0, P0], 0, 12, 0.005, { every: 100 });
      t.forEach((ti, i) => {
        const [ma, Pa] = paperNgn.analytic(ti, p, m0, P0);
        assert.ok(Math.abs(y[i][0] - ma) <= 1e-6 * Math.max(Math.abs(ma), 1e-3), `mRNA t=${ti}: ${y[i][0]} vs ${ma}`);
        assert.ok(Math.abs(y[i][1] - Pa) <= 1e-6 * Math.max(Math.abs(Pa), 1e-3), `protein t=${ti}: ${y[i][1]} vs ${Pa}`);
      });
    }
  }
});

test("paper's Neurogenin ODE: steady state k_tx k_tl / (k_dm k_dp)", () => {
  const p = { ktx: 1.7, kdm: 0.8, ktl: 2.4, kdp: 1.3 };
  const [ms, Ps] = paperNgn.steadyState(p);
  assert.ok(Math.abs(ms - p.ktx / p.kdm) < 1e-12);
  assert.ok(Math.abs(Ps - (p.ktx * p.ktl) / (p.kdm * p.kdp)) < 1e-12);
  const { y } = integrate(paperNgn.rhs(p), [0, 0], 0, 40, 0.01);
  assert.ok(Math.abs(y.at(-1)[1] - Ps) / Ps < 1e-8, `numerical ${y.at(-1)[1]} vs ${Ps}`);
  const [, Pa] = paperNgn.analytic(1e3, p);
  assert.ok(Math.abs(Pa - Ps) < 1e-12);
});

// ------------------------------------------------------------ MRF cascade
test('MyoD autoregulation is a bistable switch with hysteresis (turn-on > turn-off)', () => {
  // reversible variant: both fold points at physical inputs
  const pRev = { ...mrf.defaults, vAuto: 0.7 };
  const analytic = mrf.thresholds(pRev);
  assert.ok(analytic.bistable);
  assert.ok(analytic.on > analytic.off && analytic.off > 0, `on ${analytic.on} off ${analytic.off}`);
  const sweep = mrf.hysteresisSweep(pRev);
  assert.ok(sweep.on > sweep.off, `sweep on ${sweep.on} off ${sweep.off}`);
  assert.ok(Math.abs(sweep.on - analytic.on) < 0.01, 'sweep jump-up agrees with the lower fold point');
  assert.ok(Math.abs(sweep.off - analytic.off) < 0.01, 'sweep jump-down agrees with the upper fold point');
  // default: commitment is irreversible (the MyoD-high state survives zero input)
  const def = mrf.thresholds();
  assert.ok(def.on > 0 && def.off < 0, `default on ${def.on} off ${def.off}`);
  const high = mrf.relaxMyoD(0, 1);
  const low = mrf.relaxMyoD(0, 0);
  assert.ok(high > 0.5 && low < 0.01, 'same input, two stable states (memory of commitment)');
});

test('Wnt + Shh commit naive mesoderm to MyoD, which persists after the signals are removed', () => {
  const env = { wnt: 0.9, shh: 0.9, gf: 1 };
  const f = mrf.rhs(mrf.defaults, env);
  const y = Float64Array.from([0, 0, 0, 0, 0]);
  const r1 = integrate(f, y, 0, 60, 0.05).y.at(-1);
  assert.ok(r1[0] > 0.3 && r1[1] > 0.6, `Myf5 ${r1[0]} MyoD ${r1[1]}`);
  env.wnt = 0; env.shh = 0;
  const r2 = integrate(f, r1, 0, 100, 0.05).y.at(-1);
  assert.ok(r2[0] < 0.05, 'Myf5 follows the signal');
  assert.ok(r2[1] > 0.6, 'MyoD stays on (autoregulation)');
  const weak = integrate(mrf.rhs(mrf.defaults, { wnt: 0.25, shh: 0.25, gf: 1 }), [0, 0, 0, 0, 0], 0, 100, 0.05).y.at(-1);
  assert.ok(weak[1] < 0.2, `sub-threshold signalling does not commit (MyoD ${weak[1]})`);
});

test('growth factors keep myogenin low; withdrawal raises myogenin, MRF4 and MHC', () => {
  const run = (gf) => integrate(mrf.rhs(mrf.defaults, { wnt: 0, shh: 0, gf }), [0.3, 0.85, 0, 0, 0], 0, 96, 0.05).y.at(-1);
  const gm = run(1), dm = run(0.05);
  assert.ok(gm[2] < 0.1 && gm[4] < 0.02, `growth medium: myogenin ${gm[2]} MHC ${gm[4]}`);
  assert.ok(dm[2] > 0.5 && dm[3] > 0.4 && dm[4] > 0.5, `differentiation medium: myogenin ${dm[2]} MRF4 ${dm[3]} MHC ${dm[4]}`);
  // hierarchy: MyoD precedes myogenin, myogenin precedes MHC
  const { t, y } = integrate(mrf.rhs(mrf.defaults, { wnt: 0, shh: 0, gf: 0.05 }), [0.3, 0.85, 0, 0, 0], 0, 96, 0.05);
  const half = (k) => t[y.findIndex((v) => v[k] > 0.5 * y.at(-1)[k])];
  assert.ok(half(2) < half(3) && half(3) < half(4), `t1/2 myogenin ${half(2)} < MRF4 ${half(3)} < MHC ${half(4)}`);
});

test('segment-segment distance (2D)', () => {
  const o = {};
  segSeg2D(0, 0, 10, 0, 5, 3, 5, 8, o);
  assert.ok(Math.abs(o.d - 3) < 1e-12 && Math.abs(o.ax - 5) < 1e-12);
  segSeg2D(0, 0, 10, 0, 12, 0, 20, 0, o); // collinear end-to-end
  assert.ok(Math.abs(o.d - 2) < 1e-12);
  segSeg2D(0, 0, 10, 10, 0, 10, 10, 0, o); // crossing
  assert.ok(o.d < 1e-12);
});

test('myoblast culture: fusion conserves nuclei and volume; growth factor blocks fusion', () => {
  const dm = new MyoCulture({ n: 70, gf: 0.05, seed: 4 });
  const gm = new MyoCulture({ n: 70, gf: 1, seed: 4 });
  const nuc0 = dm.totalNuclei(), vol0 = dm.totalVolume();
  let divVol = 0;
  const origDivide = dm.divide.bind(dm);
  dm.divide = (c) => { const d = origDivide(c); divVol += myoblastVolume(d); return d; };
  for (let k = 0; k < 720; k++) { dm.step(0.1); gm.step(0.1); } // 72 h
  assert.equal(dm.totalNuclei(), nuc0 + dm.divisions, 'nuclei are conserved through fusion');
  assert.ok(Math.abs(dm.totalVolume() - vol0 - divVol) / vol0 < 1e-9, 'cytoplasmic volume is conserved through fusion');
  assert.ok(dm.fusions > 10 && dm.tubes.length > 3, `fusions ${dm.fusions}`);
  const fi = dm.fusionIndex();
  assert.ok(fi > 0.4 && fi < 0.95, `fusion index in differentiation medium ${fi}`);
  assert.ok(dm.tubes.some((t) => t.nuclei.length >= 4), 'multinucleated myotubes form');
  assert.equal(gm.fusionIndex(), 0, 'no fusion in growth medium');
  assert.ok(gm.cells.length > 70, 'myoblasts proliferate in growth medium');
  // myotubes are longer than myoblasts, and their state has progressed to MHC
  const m = dm.meanState();
  assert.ok(m[4] > 0.3, `mean MHC ${m[4]}`);
  for (const t of dm.tubes) assert.ok(t.len > dm.p.L * 0.9);
});

// ------------------------------------------------------------ neural network
function cluster(seed = 5, n = 90) {
  const rng = new RNG(seed);
  const tis = new Tissue({ noise: 0, growth: false, adhesion: () => 2, rng, muRep: 10 });
  for (let i = 0; i < n; i++) { const d = rng.onSphere(); const r = 36 * Math.cbrt(rng.next()); tis.add(new Cell({ p: d.map((v) => v * r), R: 6.5 })); }
  for (let k = 0; k < 300; k++) tis.step(0.05);
  const idx = new Map(tis.cells.map((c, i) => [c.id, i]));
  return tis.cells.map((c) => (tis.neighbours.get(c.id) ?? []).map((id) => idx.get(id)));
}

function runGRN(nbrs, over, T = 24) {
  const rng = new RNG(3);
  const p = { ...neuralGRN.defaults, ...over };
  const Y = nbrs.map(() => neuralGRN.initial(rng));
  const com = nbrs.map(() => false), when = nbrs.map(() => -1);
  const work = neuralGRN.work();
  for (let t = 0; t < T; t += 0.05) {
    for (const i of neuralGRN.step(Y, com, nbrs, 0.05, p, work)) when[i] = t;
  }
  return { com, when, Y };
}

test('Delta-Notch lateral inhibition selects isolated neurons; DAPT or exogenous Notch flip the outcome', () => {
  const nbrs = cluster();
  const n = nbrs.length;
  const normal = runGRN(nbrs, {});
  const nNeurons = normal.com.filter(Boolean).length;
  assert.ok(nNeurons > 0.05 * n && nNeurons < 0.6 * n, `neurons ${nNeurons}/${n}`);
  // salt and pepper: neurons committing together are not neighbours
  let bad = 0;
  for (let i = 0; i < n; i++) if (normal.com[i]) for (const j of nbrs[i]) if (normal.com[j] && Math.abs(normal.when[i] - normal.when[j]) < 1) bad++;
  assert.equal(bad, 0, 'no simultaneous neuron-neuron neighbours');
  // neurons have lost Sox2; progenitors keep it
  const soxN = normal.Y.filter((_, i) => normal.com[i]).reduce((s, y) => s + y[0], 0) / nNeurons;
  const soxP = normal.Y.filter((_, i) => !normal.com[i]).reduce((s, y) => s + y[0], 0) / (n - nNeurons);
  assert.ok(soxN < soxP, `Sox2 neurons ${soxN} < progenitors ${soxP}`);
  const dapt = runGRN(nbrs, { coupling: 0 });
  assert.equal(dapt.com.filter(Boolean).length, n, 'gamma-secretase block: every progenitor differentiates');
  const notch = runGRN(nbrs, { notchExt: 0.15 });
  assert.equal(notch.com.filter(Boolean).length, 0, 'constitutive Notch keeps progenitors');
});

// ------------------------------------------------------------ guidance field
test('guidance field: analytic gradient matches finite differences and points to the source', () => {
  const p = guidance.defaults;
  const sources = [{ p: [0, 0, 0], radius: 6 }, { p: [80, 40, -30], radius: 6, Q: 2500 }];
  const rng = new RNG(8);
  const g = [0, 0, 0];
  for (let k = 0; k < 40; k++) {
    const x = [rng.uniform(-150, 200), rng.uniform(-120, 150), rng.uniform(-150, 120)];
    if (sources.some((s) => Math.hypot(x[0] - s.p[0], x[1] - s.p[1], x[2] - s.p[2]) < 10)) continue;
    const C = guidance.gradient(x, sources, p, g);
    assert.ok(Math.abs(C - guidance.concentration(x, sources, p)) < 1e-12 * Math.max(C, 1));
    const h = 1e-3;
    for (let a = 0; a < 3; a++) {
      const xp = x.slice(), xm = x.slice();
      xp[a] += h; xm[a] -= h;
      const fd = (guidance.concentration(xp, sources, p) - guidance.concentration(xm, sources, p)) / (2 * h);
      assert.ok(Math.abs(fd - g[a]) <= 1e-6 * Math.hypot(...g) + 1e-12, `axis ${a}: fd ${fd} analytic ${g[a]}`);
    }
  }
  // single source: gradient is anti-parallel to (x - x_s) and C(r) = Q exp(-r/l)/(4 pi D r)
  const one = [{ p: [10, -5, 3], radius: 1 }];
  const x = [70, 20, -40];
  guidance.gradient(x, one, p, g);
  const d = x.map((v, i) => one[0].p[i] - v);
  const cos = (g[0] * d[0] + g[1] * d[1] + g[2] * d[2]) / (Math.hypot(...g) * Math.hypot(...d));
  assert.ok(cos > 1 - 1e-12, `cos ${cos}`);
  const r = Math.hypot(...d), l = guidance.length(p);
  const Cexp = (p.Q / (4 * Math.PI * p.D * r)) * Math.exp(-r / l) / 0.6022;
  assert.ok(Math.abs(guidance.concentration(x, one, p) - Cexp) / Cexp < 1e-12);
  // l = sqrt(D / lambda) with D in um^2/s and lambda in 1/h
  assert.ok(Math.abs(l - Math.sqrt((p.D * 3600) / p.lambda)) < 1e-12);
  // steady state satisfies D lap C - lambda C = 0 away from the source (radial Laplacian, finite differences)
  const Cr = (rr) => guidance.concentration([one[0].p[0] + rr, one[0].p[1], one[0].p[2]], one, p);
  const rr = 60, hh = 0.01;
  const lap = (Cr(rr + hh) - 2 * Cr(rr) + Cr(rr - hh)) / (hh * hh) + (2 / rr) * (Cr(rr + hh) - Cr(rr - hh)) / (2 * hh);
  assert.ok(Math.abs(p.D * 3600 * lap - p.lambda * Cr(rr)) < 1e-4 * p.lambda * Cr(rr), 'reaction-diffusion residual');
});

test('gradient-sensing SNR peaks near C = Kd and falls with distance', () => {
  const p = guidance.defaults;
  const src = [{ p: [0, 0, 0], radius: 5 }];
  const g = [0, 0, 0];
  const snrAt = (r) => { const C = guidance.gradient([r, 0, 0], src, p, g); return guidance.snr(C, Math.hypot(...g), p); };
  assert.ok(snrAt(40) > snrAt(120) && snrAt(120) > snrAt(250));
  // for a fixed relative gradient, the occupancy-difference signal is maximal at C = Kd
  const rel = 0.02 / p.width; // 2 % across the growth cone
  const s = (C) => guidance.snr(C, C * rel, p);
  assert.ok(s(p.Kd) > s(p.Kd * 10) && s(p.Kd) > s(p.Kd / 10));
});

test('growth cones steered by the field reach the target and form a synapse; without chemotaxis they do not', () => {
  const src = [{ p: [0, 0, 0], radius: 8, id: 'target' }];
  const trial = (kappa, seed) => {
    const ng = new NeuriteGrowth({ seed, outgrowth: { kappa, branchRate: 0, dendrites: 0 } });
    ng.neurons.set(1, { axonBranches: 0, synapses: 0 });
    ng.tips.push(ng.tip(1, 0, [140, 0, 0], [0, 1, 0], 1, 0, 600));
    for (let t = 0; t < 10 && !ng.synapses.length; t += 0.05) ng.step(0.05, src);
    return ng;
  };
  let guided = 0, blind = 0;
  for (let seed = 1; seed <= 8; seed++) {
    const a = trial(0.3, seed);
    if (a.synapses.length) {
      guided++;
      const s = a.synapses[0];
      assert.ok(Math.abs(Math.hypot(...s.pos) - 8) < 1e-9, 'bouton sits on the target surface');
      assert.ok(Math.abs(a.totalLength - a.segments.length * a.p.step) < 1e-9, 'total length = sum of segments');
    }
    if (trial(0, seed).synapses.length) blind++;
  }
  assert.equal(guided, 8);
  assert.ok(blind <= 1, `unguided hits ${blind}`);
});

// ------------------------------------------------------------ Goodwin
test('Goodwin oscillator: n = 12 sustains oscillation, n = 4 is damped (Griffith n > 8)', () => {
  const sim = (n) => {
    const p = { ...goodwin.defaults, n };
    const { t, y } = integrate(goodwin.rhs(p), [0, 0, 0, 0], 0, 3000, 0.05, { every: 10 });
    return { p, t, Z: y.map((v) => v[2]), y };
  };
  const hi = sim(12), lo = sim(4);
  const oh = oscillationStats(hi.t, hi.Z, { from: 1500 });
  const ol = oscillationStats(lo.t, lo.Z, { from: 1500 });
  assert.ok(oh.sustained && oh.relAmplitude > 0.2 && oh.decay > 0.98, `n=12 ${JSON.stringify({ a: oh.relAmplitude, d: oh.decay })}`);
  assert.ok(!ol.sustained && ol.relAmplitude < 1e-3, `n=4 rel. amplitude ${ol.relAmplitude}`);
  // n = 4 converges to the analytic steady state
  const ss = goodwin.steadyState(lo.p);
  lo.y.at(-1).forEach((v, i) => assert.ok(Math.abs(v - ss[i]) / ss[i] < 1e-4, `${goodwin.species[i]} ${v} vs ${ss[i]}`));
  // linear stability agrees: loop gain vs Hopf threshold 8 for equal rates
  const s12 = goodwin.stability(hi.p), s4 = goodwin.stability(lo.p);
  assert.equal(s12.gainCrit, 8);
  assert.ok(!s12.stable && s4.stable);
  const nc = goodwin.criticalN(goodwin.defaults);
  assert.ok(nc > 8 && nc < 8.2, `critical n ${nc}`);
  // near the Hopf point the period approaches 2 pi / omega (omega^2 = bd + bf + df)
  const near = sim(nc + 0.05);
  const on = oscillationStats(near.t, near.Z, { from: 2000 });
  assert.ok(Math.abs(on.period - s12.periodHopf) / s12.periodHopf < 0.05, `period ${on.period} vs ${s12.periodHopf}`);
  // receptor occupancy follows mass action
  const O = ss[3], Z = ss[2];
  assert.ok(Math.abs(O - Z / (Z + lo.p.koff / lo.p.kon)) < 1e-12);
});

test('Goodwin: unequal degradation rates raise the Hopf threshold above 8 (secant condition)', () => {
  const p = { ...goodwin.defaults, b: 0.02, d: 0.1, f: 0.5 };
  const { gainCrit } = goodwin.stability(p);
  assert.ok(gainCrit > 8, `Gc ${gainCrit}`); // equal rates are the most favourable case
  // so the same n = 12 loop that oscillates with equal rates is stable here
  assert.ok(goodwin.stability({ ...p, n: 12 }).stable);
  const { t, y } = integrate(goodwin.rhs({ ...p, n: 12 }), [0, 0, 0, 0], 0, 3000, 0.05, { every: 10 });
  assert.ok(!oscillationStats(t, y.map((v) => v[2]), { from: 2000 }).sustained);
});

// ------------------------------------------------------------ reaction lists
test('stage-3 REACTIONS are well formed and order 3', () => {
  const all = [...MYO, ...NEURO, ...HORM];
  const ids = new Set();
  for (const r of all) {
    assert.equal(r.order, 3);
    assert.ok(typeof r.id === 'string' && !ids.has(r.id), r.id);
    ids.add(r.id);
    assert.ok(Array.isArray(r.reactants) && Array.isArray(r.products) && r.label && r.pathway);
    assert.ok(r.reactants.length + r.products.length > 0, r.id);
    if (r.modifiers) assert.ok(Array.isArray(r.modifiers));
  }
  assert.ok(MYO.length >= 8 && NEURO.length >= 8 && HORM.length >= 6);
});
