import test from 'node:test';
import assert from 'node:assert/strict';
import { RNG } from '../../js/models/core/rng.js';
import { rk45 } from '../../js/models/core/ode.js';
import * as mapk from '../../js/models/pathways/mapk.js';
import * as js from '../../js/models/pathways/jakstat.js';
import * as hema from '../../js/models/pathways/hematopoiesis.js';

const { toggle } = hema;

// ------------------------------------------------------------------ MAPK

test('MAPK cascade is ultrasensitive and sensitivity grows down the cascade (Huang & Ferrell 1996)', () => {
  const h = mapk.cascadeHill();
  // a single Michaelis-Menten cycle far from saturation is hyperbolic: nH = 1
  assert.ok(Math.abs(h.MKKK.nH - 1) < 0.05, `MAPKKK nH ${h.MKKK.nH}`);
  assert.ok(h.MKK.nH > h.MKKK.nH + 0.5, `MAPKK nH ${h.MKK.nH}`);
  assert.ok(h.MAPK.nH > 2, `MAPK nH ${h.MAPK.nH}`);
  assert.ok(h.MAPK.nH > h.MKK.nH, `MAPK ${h.MAPK.nH} vs MAPKK ${h.MKK.nH}`);
  // amplification: the last tier switches at far lower input than the first tier's half-activation
  assert.ok(h.MAPK.EC50 < 0.2 * h.MKKK.EC50, `EC50 ${h.MAPK.EC50} vs ${h.MKKK.EC50}`);
});

test('MAPK steady state: conservation of every tier and zero net flux at each step', () => {
  const p = mapk.HF96;
  for (const E1 of [1e-6, 2e-5, 2.5e-5, 1e-4, 1e-2]) {
    const s = mapk.cascadeSteadyState(E1, p);
    assert.ok(Math.abs(s.MKK + s.MKKP + s.MKKPP - p.MKKtot) < 1e-12);
    assert.ok(Math.abs(s.MAPK + s.MAPKP + s.MAPKPP - p.MAPKtot) < 1e-12);
    assert.ok(s.MKKKa >= 0 && s.MKKKa <= p.MKKKtot);
    // first MAPK step: kinase flux (competitive MM) equals phosphatase flux
    const K = p.Km;
    const kin = (p.kcat * s.MKKPP * s.MAPK) / K / (1 + (s.MAPK + s.MAPKP) / K);
    const pho = (p.kcat * p.MAPKPase * s.MAPKP) / K / (1 + (s.MAPKP + s.MAPKPP) / K);
    assert.ok(Math.abs(kin - pho) <= 1e-8 * Math.max(kin, 1e-12), `E1 ${E1}: ${kin} vs ${pho}`);
  }
  // responses are monotone in the stimulus
  const d = mapk.doseResponse([1e-6, 1e-5, 2e-5, 3e-5, 1e-4, 1e-3]);
  for (const k of ['MKKK', 'MKK', 'MAPK']) for (let i = 1; i < d[k].length; i++) assert.ok(d[k][i] >= d[k][i - 1] - 1e-12);
});

test('MAPK ODE (nuclear shuttle off) converges to the algebraic steady state', () => {
  const p = { ...mapk.model.defaults, L: 1, kin: 0 };
  const r = rk45(mapk.model.rhs(p), mapk.model.initial(), 0, 3000, { rtol: 1e-7, atol: 1e-11, h0: 0.01 });
  const ss = mapk.cascadeSteadyState(mapk.inputFromLigand(1, p), p);
  const y = r.y;
  const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-9);
  assert.ok(rel(y[2], ss.MKKKa) < 1e-3, `Raf* ${y[2]} vs ${ss.MKKKa}`);
  assert.ok(rel(y[4], ss.MKKPP) < 1e-3, `MEK-PP ${y[4]} vs ${ss.MKKPP}`);
  assert.ok(rel(y[6], ss.MAPKPP) < 1e-3, `ERK-PP ${y[6]} vs ${ss.MAPKPP}`);
  // receptor steady state matches the closed form
  assert.ok(Math.abs((2 * y[1]) / p.Rtot - mapk.receptorActivation(1, p)) < 1e-4);
});

test('MAPK full model with nuclear translocation conserves receptor and ERK and induces c-Fos', () => {
  const p = { ...mapk.model.defaults, L: 1 };
  const f = mapk.model.rhs(p);
  let ok = true;
  const r = rk45(f, mapk.model.initial(), 0, 1500, {
    rtol: 1e-6, atol: 1e-10, h0: 0.01,
    onStep: (t, y) => {
      const tot = mapk.model.totals(y, p);
      if (tot.RTK > p.Rtot + 1e-9 || tot.ERKphosOrNuclear > p.MAPKtot + 1e-9) ok = false;
      if (y.some((v) => v < -1e-9)) ok = false;
    },
  });
  assert.ok(ok, 'conservation bounds or positivity violated');
  // conservation: the implicit (eliminated) pools, differentiated from their own fluxes,
  // cancel the summed derivatives of the explicit species of the same protein
  const dy = new Float64Array(11);
  f(0, r.y, dy);
  const y = r.y, K = p.Km;
  const R = p.Rtot - y[0] - 2 * y[1];
  assert.ok(Math.abs(-(p.kon * p.L * R - p.koff * y[0]) + dy[0] + 2 * dy[1]) < 1e-12, 'receptor');
  const M0 = p.MKKtot - y[3] - y[4];
  const dM0 = -(p.kcat * y[2] * M0) / K / (1 + (M0 + y[3]) / K) + (p.kcat * p.MKKPase * y[3]) / K / (1 + (y[3] + y[4]) / K);
  assert.ok(Math.abs(dM0 + dy[3] + dy[4]) < 1e-9, 'MEK');
  const E0 = p.MAPKtot - y[5] - y[6] - y[7] - y[8];
  const dE0 = -(p.kcat * y[4] * E0) / K / (1 + (E0 + y[5]) / K) + (p.kcat * p.MAPKPase * y[5]) / K / (1 + (y[5] + y[6]) / K) + p.kex * y[8];
  assert.ok(Math.abs(dE0 + dy[5] + dy[6] + dy[7] + dy[8]) < 1e-9, 'ERK');
  assert.ok(r.y[7] > 0.05 * p.MAPKtot, 'ERK-PP accumulates in the nucleus');
  assert.ok(r.y[10] > 0 && dy[10] > -1e-6, 'c-Fos induced');
  // no ligand -> no signalling
  const r0 = rk45(mapk.model.rhs({ ...p, L: 0 }), mapk.model.initial(), 0, 500, { h0: 0.01 });
  assert.ok(r0.y.every((v) => Math.abs(v) < 1e-12));
});

test('growth-factor dose response of ERK is switch-like', () => {
  assert.ok(mapk.nuclearERK(0.05) < 0.05);
  assert.ok(mapk.nuclearERK(3) > 0.9);
  const h = mapk.effectiveHill((L) => mapk.nuclearERK(L) / mapk.nuclearERK(1000), 1e-4, 1000);
  assert.ok(h.nH > 2, `SCF -> ERK nH ${h.nH}`);
});

// ------------------------------------------------------------------ JAK/STAT

test('JAK2/STAT5 conserves total STAT5 through phosphorylation, dimerisation and nuclear cycling', () => {
  for (const EPO of [0.1, 1, 10]) {
    const p = { ...js.jakstat.defaults, EPO };
    const { y } = js.timeCourse(p, 180, { every: 2 });
    for (const s of y) assert.ok(Math.abs(js.jakstat.totalSTAT(s) - p.STATtot) < 1e-9, `total ${js.jakstat.totalSTAT(s)}`);
    for (const s of y) assert.ok(s.every((v) => v > -1e-12));
  }
});

test('CIS/SOCS negative feedback lowers peak and sustained nuclear pSTAT5', () => {
  const p = js.jakstat.defaults;
  const withFb = js.responseSummary(p);
  const noFb = js.responseSummary({ ...p, fb: 0 });
  assert.ok(withFb.peak < 0.9 * noFb.peak, `peak ${withFb.peak} vs ${noFb.peak}`);
  assert.ok(withFb.late < 0.5 * noFb.late, `late ${withFb.late} vs ${noFb.late}`);
  // adaptation: with feedback the response overshoots and relaxes below its peak
  assert.ok(withFb.late < 0.5 * withFb.peak);
  // nuclear accumulation lags phosphorylation (Swameye et al.: minutes, not seconds)
  assert.ok(withFb.tPeak > 5 && withFb.tPeak < 40, `tPeak ${withFb.tPeak}`);
  // no EPO -> no signal; drive increases with EPO
  assert.equal(js.stat5Drive(0), 0);
  let prev = 0;
  for (const E of [0.01, 0.1, 1, 10]) { const d = js.stat5Drive(E); assert.ok(d > prev); prev = d; }
});

// ------------------------------------------------------------------ GATA1-PU.1 toggle

test('toggle without self-activation is a bistable switch (two attractors + saddle)', () => {
  const p = { ...toggle.defaults, a1: 0, a2: 0 };
  const fps = toggle.fixedPoints(p);
  const stable = fps.filter((f) => f.stable), saddles = fps.filter((f) => f.type === 'saddle');
  assert.equal(stable.length, 2);
  assert.equal(saddles.length, 1);
  // analytic saddle at x = theta by symmetry
  assert.ok(Math.abs(saddles[0].g - 0.5) < 1e-6 && Math.abs(saddles[0].u - 0.5) < 1e-6);
  // attractors are mirror images: (x_hi, x_lo) and (x_lo, x_hi)
  assert.ok(Math.abs(stable[0].g - stable[1].u) < 1e-6 && Math.abs(stable[0].u - stable[1].g) < 1e-6);
});

test('strong self-activation adds a central primed attractor (Huang et al. 2007)', () => {
  const p = { ...toggle.defaults, a1: 1, a2: 1 };
  const stable = toggle.fixedPoints(p).filter((f) => f.stable);
  assert.equal(stable.length, 3);
  // with a = b = 1, x = 1 solves the symmetric equation exactly
  assert.ok(stable.some((f) => Math.abs(f.g - 1) < 1e-6 && Math.abs(f.u - 1) < 1e-6));
  const ac = toggle.criticalSelfActivation();
  assert.ok(ac > 0 && ac < 1, `a_c ${ac}`);
  // just below a_c the centre is gone
  const below = toggle.fixedPoints({ ...toggle.defaults, a1: ac - 0.05, a2: ac - 0.05 }).filter((f) => f.stable);
  assert.equal(below.length, 2);
});

test('Langevin step reduces to Euler without noise and noise-free cells stay in their attractor', () => {
  const p = { ...toggle.defaults, a1: 0, a2: 0 };
  const s = { g: 0.8, u: 0.3 };
  const r = toggle.rates(0.8, 0.3, p);
  toggle.langevinStep(s, 0.01, p, 0, new RNG(1));
  assert.ok(Math.abs(s.g - (0.8 + 0.01 * r[0])) < 1e-12 && Math.abs(s.u - (0.3 + 0.01 * r[1])) < 1e-12);
  const hi = toggle.fixedPoints(p).find((f) => f.stable && f.g > f.u);
  const c = { g: hi.g + 0.05, u: hi.u + 0.05 };
  for (let i = 0; i < 2000; i++) toggle.langevinStep(c, 0.01, p, 0, new RNG(1));
  assert.ok(Math.abs(c.g - hi.g) < 1e-3 && Math.abs(c.u - hi.u) < 1e-3);
});

test('EPO (GATA1 input) shifts the commitment probability toward the GATA1-high fate', () => {
  const opts = { n: 400, sigma: 0.1 };
  const p0 = hema.fateProbability({ ...toggle.defaults }, { ...opts, rng: new RNG(7) });
  const pE = hema.fateProbability({ ...toggle.defaults, uG: 0.06 }, { ...opts, rng: new RNG(7) });
  const pM = hema.fateProbability({ ...toggle.defaults, uP: 0.06 }, { ...opts, rng: new RNG(7) });
  assert.ok(Math.abs(p0 - 0.5) < 0.08, `unbiased ${p0}`);
  assert.ok(pE > p0 + 0.2, `EPO ${pE} vs ${p0}`);
  assert.ok(pM < p0 - 0.2, `GM-CSF ${pM} vs ${p0}`);
  // every simulated cell is committed at the end of the protocol (far from the diagonal)
  const s = hema.commitCell({ g: 1, u: 1 }, toggle.defaults, { rng: new RNG(3), sigma: 0.1 });
  assert.ok(Math.abs(s.g - s.u) > 0.5);
});

// ------------------------------------------------------------------ lineage

test('stage rates reproduce amplification and residence of the deterministic stages', () => {
  const eb = hema.stageRates(hema.LINEAGE.EB);
  assert.equal(eb.A, 8);
  assert.ok(Math.abs(eb.exit / (eb.exit - eb.div) - 8) < 1e-12);
  // one proerythroblast -> 16 reticulocytes
  const pro = hema.stageRates(hema.LINEAGE.ProEB);
  assert.equal(pro.A * eb.A, 16);
  // every transition points to a node of the tree
  for (const tr of hema.TRANSITIONS) {
    assert.ok(hema.LINEAGE[tr.from], tr.id);
    for (const t of tr.to) assert.ok(t === 'Pyrenocyte' || hema.LINEAGE[t], `${tr.id} -> ${t}`);
  }
});

test('niche-limited self-renewal gives a finite HSC steady state matching the mean field', () => {
  const env = { selfRenewal: 0.8, K: 200 };
  const mf = hema.meanFieldSteadyState({ ...hema.SSA_DEFAULTS, ...env });
  assert.ok(Math.abs(mf.HSC - 200 * (1 - 1 / 1.6)) < 1e-9);
  // start far below the steady state
  const ssa = new hema.LineageSSA({ rng: new RNG(2024), env, init: { HSC: 10 } });
  const ts = [];
  let sumH = 0, sumN = 0, sumRBC = 0, n = 0;
  for (let d = 0; d < 700; d++) {
    ssa.step(1);
    ts.push(ssa.count('HSC'));
    if (d >= 300) { sumH += ssa.count('HSC'); sumN += ssa.count('Neu'); n++; }
    if (d >= 500) sumRBC += ssa.count('RBC');
  }
  const H = sumH / n;
  assert.ok(Math.abs(H - mf.HSC) / mf.HSC < 0.1, `HSC ${H} vs ${mf.HSC}`);
  assert.ok(Math.max(...ts) < 200 && Math.min(...ts.slice(100)) > 0);
  assert.ok(Math.abs(sumN / n - mf.Neu) / mf.Neu < 0.15, `Neu ${sumN / n} vs ${mf.Neu}`);
  // RBCs (120-day lifespan) have filled up after > 4 lifespans
  assert.ok(Math.abs(sumRBC / 200 - mf.RBC) / mf.RBC < 0.15, `RBC ${sumRBC / 200} vs ${mf.RBC}`);
});

test('self-renewal below one half exhausts the stem-cell pool', () => {
  const ssa = new hema.LineageSSA({ rng: new RNG(5), env: { selfRenewal: 0.3 }, init: { HSC: 40 } });
  for (let d = 0; d < 300; d++) ssa.step(1);
  assert.equal(ssa.count('HSC'), 0);
  assert.equal(hema.meanFieldSteadyState({ ...hema.SSA_DEFAULTS, selfRenewal: 0.3 }).HSC, 0);
});

test('EPO withdrawal collapses erythroid output but not myeloid output', () => {
  const hi = hema.meanFieldSteadyState({ ...hema.SSA_DEFAULTS, stat5: 0.9 });
  const lo = hema.meanFieldSteadyState({ ...hema.SSA_DEFAULTS, stat5: 0.02 });
  assert.ok(lo.Retic < 0.5 * hi.Retic, `${lo.Retic} vs ${hi.Retic}`);
  assert.ok(Math.abs(lo.Neu - hi.Neu) < 1e-9);
});

test('SCF/c-Kit withdrawal (no ERK) cuts progenitor output but leaves the HSC pool intact', () => {
  const on = hema.meanFieldSteadyState({ ...hema.SSA_DEFAULTS, erk: 1 });
  const off = hema.meanFieldSteadyState({ ...hema.SSA_DEFAULTS, erk: 0 });
  assert.equal(on.HSC, off.HSC);
  assert.ok(off.Neu < 0.3 * on.Neu, `Neu ${off.Neu} vs ${on.Neu}`);
  assert.ok(off.RBC < 0.3 * on.RBC, `RBC ${off.RBC} vs ${on.RBC}`);
  // the SCF -> ERK map is switch-like, so the loss happens over a narrow SCF range
  assert.ok(hema.apoptosisRate('GMP', { ...hema.SSA_DEFAULTS, erk: mapk.nuclearERK(3) }) < 0.1);
  assert.ok(hema.apoptosisRate('GMP', { ...hema.SSA_DEFAULTS, erk: mapk.nuclearERK(0.1) }) > 0.5);
});

test('every stage-2 module exports well-formed order-2 REACTIONS', () => {
  for (const [name, R] of [['mapk', mapk.REACTIONS], ['jakstat', js.REACTIONS], ['hematopoiesis', hema.REACTIONS]]) {
    assert.ok(R.length >= 10, name);
    const ids = new Set();
    for (const r of R) {
      assert.equal(r.order, 2, r.id);
      assert.ok(typeof r.pathway === 'string' && typeof r.label === 'string' && r.label.length > 3, r.id);
      assert.ok(Array.isArray(r.reactants) && Array.isArray(r.products), r.id);
      assert.ok(r.reactants.length + r.products.length > 0, r.id);
      for (const s of [...r.reactants, ...r.products, ...(r.modifiers ?? [])]) assert.ok(typeof s === 'string' && s.length > 0, r.id);
      assert.ok(!ids.has(r.id), `duplicate ${r.id}`);
      ids.add(r.id);
    }
  }
});
