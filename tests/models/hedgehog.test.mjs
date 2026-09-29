import test from 'node:test';
import assert from 'node:assert/strict';
import { integrate } from '../../js/models/core/ode.js';
import {
  hedgehog, frenchFlag, frenchFlagField, NEURAL_TUBE_THRESHOLDS, NEURAL_TUBE_DOMAINS, REACTIONS,
} from '../../js/models/pathways/hedgehog.js';

const P = hedgehog.defaults;
const DOSES = [0, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100];

test('Shh dose-response is monotone: GliA and target rise, GliR and free PTCH1 fall', () => {
  const d = hedgehog.doseResponse(P, DOSES);
  for (let i = 1; i < d.length; i++) {
    assert.ok(d[i].A > d[i - 1].A, `GliA not increasing at Shh ${d[i].shh}`);
    assert.ok(d[i].readout > d[i - 1].readout, `read-out not increasing at Shh ${d[i].shh}`);
    assert.ok(d[i].R < d[i - 1].R, `GliR not decreasing at Shh ${d[i].shh}`);
    assert.ok(d[i].S > d[i - 1].S, `SMO activity not increasing at Shh ${d[i].shh}`);
    const free = (r) => (1 - hedgehog.boundFraction(r.shh, P.Kd)) * r.P;
    assert.ok(free(d[i]) < free(d[i - 1]), `free PTCH1 not decreasing at Shh ${d[i].shh}`);
  }
  // switch-like: off without ligand, near-saturated at 100 Kd
  assert.ok(d[0].readout < 0.05, `basal read-out ${d[0].readout}`);
  assert.ok(d.at(-1).readout > 0.85, `saturated read-out ${d.at(-1).readout}`);
  // Ptch1 is a GLI target: total receptor rises at high dose even though free PTCH1 falls
  assert.ok(d.at(-1).P > d[0].P);
});

test('steady state is unique (monostable): high and low initial conditions agree', () => {
  for (const shh of [0, 1, 3, 30]) {
    const lo = hedgehog.steadyState({ ...P, shh });
    const hi = hedgehog.steadyState({ ...P, shh }, { y0: [0.1, 0.95, 0.1, 3, 0, 4] });
    lo.forEach((v, i) => assert.ok(Math.abs(v - hi[i]) < 1e-4, `Shh ${shh} ${hedgehog.species[i]}: ${v} vs ${hi[i]}`));
  }
});

test('Shh binding (QSS) and catalytic PTCH1 inhibition match their analytic forms', () => {
  assert.equal(hedgehog.boundFraction(1, 1), 0.5);
  assert.equal(hedgehog.boundFraction(0, 1), 0);
  assert.ok(Math.abs(hedgehog.boundFraction(9, 1) - 0.9) < 1e-12);
  // at steady state dS/dt = 0 gives S* = kSa / (kSa + kSi + kPtc (1 - f_b) P*)
  for (const shh of [0, 0.5, 2, 20]) {
    const q = { ...P, shh };
    const y = hedgehog.steadyState(q);
    const S = hedgehog.smoSteady(y[0], q);
    assert.ok(Math.abs(S - y[1]) < 1e-6, `Shh ${shh}: S ${y[1]} vs analytic ${S}`);
  }
  // unbound PTCH1 (no Shh): SMO activity is low
  assert.ok(hedgehog.steadyState(P)[1] < 0.1);
});

test('negative feedback: Ptch1 induction lowers GliA versus a no-feedback cell', () => {
  for (const shh of [1, 3, 10]) {
    const fb = hedgehog.steadyState({ ...P, shh });
    const nofb = hedgehog.steadyState({ ...P, shh, sP: 0 });
    assert.ok(fb[0] > nofb[0], `feedback raises PTCH1 at Shh ${shh}`);
    assert.ok(fb[3] < nofb[3], `GliA with feedback ${fb[3]} vs without ${nofb[3]} (Shh ${shh})`);
  }
  // feedback shifts the dose-response to higher Shh (desensitisation)
  const half = (pp) => {
    const d = hedgehog.doseResponse(pp, [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]);
    const top = hedgehog.doseResponse(pp, [1000])[0].readout;
    return d.find((r) => r.readout >= top / 2)?.shh ?? Infinity;
  };
  assert.ok(half(P) > half({ ...P, sP: 0 }));
});

test('temporal adaptation (Dessaud et al. 2007): GLI activity overshoots then declines with feedback only', () => {
  const run = (sP) => {
    const q = { ...P, sP };
    const y0 = hedgehog.steadyState({ ...q, shh: 0 });
    const { y } = integrate(hedgehog.rhs({ ...q, shh: 3 }), y0, 0, 200, 0.01, { every: 10 });
    const r = y.map((s) => hedgehog.promoter(s[3], s[4], q));
    return { peak: Math.max(...r), end: r.at(-1) };
  };
  const fb = run(P.sP), nofb = run(0);
  assert.ok(fb.peak > 1.15 * fb.end, `feedback: peak ${fb.peak} end ${fb.end}`);
  assert.ok(nofb.peak < 1.01 * nofb.end, `no feedback: peak ${nofb.peak} end ${nofb.end}`);
});

test('cyclopamine (SMO activation blocked) abolishes GliA and the target even at high Shh', () => {
  const on = hedgehog.steadyState({ ...P, shh: 10 });
  const cyc = hedgehog.steadyState({ ...P, shh: 10, kSa: 0 });
  assert.ok(cyc[1] < 1e-9, `SMO* ${cyc[1]}`);
  assert.ok(cyc[3] < 0.02 * on[3], `GliA ${cyc[3]} vs ${on[3]}`);
  assert.ok(cyc[5] < 0.02 * on[5], `target ${cyc[5]} vs ${on[5]}`);
  // Gli3 repressor dominates instead
  assert.ok(cyc[4] > on[4]);
});

test('French flag: thresholds map a graded signal to ordered, contiguous domains', () => {
  assert.equal(frenchFlag(0.1, [0.3, 0.6]), 0);
  assert.equal(frenchFlag(0.3, [0.3, 0.6]), 1);
  assert.equal(frenchFlag(0.59, [0.3, 0.6]), 1);
  assert.equal(frenchFlag(0.9, [0.3, 0.6]), 2);
  assert.equal(NEURAL_TUBE_DOMAINS.length, NEURAL_TUBE_THRESHOLDS.length + 1);
  // cells along an exponential Shh gradient (ventral -> dorsal) read out by the full model
  const xs = Array.from({ length: 40 }, (_, i) => i / 39);
  const shh = xs.map((x) => 60 * Math.exp(-x / 0.15));
  const read = hedgehog.doseResponse(P, shh).map((r) => r.readout);
  const dom = frenchFlagField(read);
  for (let i = 1; i < dom.length; i++) assert.ok(dom[i] <= dom[i - 1], 'domain index must fall with distance');
  // every neural-tube domain appears, each as one contiguous block
  for (let k = 0; k < NEURAL_TUBE_DOMAINS.length; k++) {
    const idx = [...dom].map((d, i) => (d === k ? i : -1)).filter((i) => i >= 0);
    assert.ok(idx.length > 0, `domain ${NEURAL_TUBE_DOMAINS[k].name} missing`);
    assert.equal(idx.at(-1) - idx[0] + 1, idx.length, `domain ${k} not contiguous`);
  }
  assert.equal(dom[0], NEURAL_TUBE_DOMAINS.length - 1, 'floor plate at the source');
});

test('REACTIONS are well-formed order-1 Hedgehog steps with unique ids', () => {
  const ids = new Set();
  for (const r of REACTIONS) {
    assert.equal(r.order, 1);
    assert.equal(r.pathway, 'Hedgehog');
    assert.ok(Array.isArray(r.reactants) && Array.isArray(r.products));
    assert.ok(r.label && !ids.has(r.id));
    ids.add(r.id);
  }
  const names = new Set(REACTIONS.flatMap((r) => [...r.reactants, ...r.products, ...(r.modifiers ?? [])]));
  for (const s of ['Shh', 'PTCH1', 'SMO*', 'GliA', 'GLI3R', 'GSK-3b', 'beta-TrCP']) assert.ok(names.has(s), s);
});
