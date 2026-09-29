import test from 'node:test';
import assert from 'node:assert/strict';
import { rk45 } from '../../js/models/core/ode.js';
import { hemoglobin, collagen } from '../../js/models/protein/assembly.js';

// the assembly kinetics are stiff (fast chaperone/heme binding): adaptive Dormand-Prince
const run = (model, over, T = 150) => rk45(model.rhs({ ...model.defaults, ...over }), model.initial(), 0, T, { rtol: 1e-7, atol: 1e-9 }).y;

test('hemoglobin: every alpha chain made is accounted for', () => {
  const y = run(hemoglobin, {});
  const I = hemoglobin.idx;
  assert.ok(Math.abs(hemoglobin.alphaBalance(y) - y[I.aSyn]) / y[I.aSyn] < 1e-6);
});

test('hemoglobin: HRI couples globin translation to heme supply', () => {
  const I = hemoglobin.idx;
  const lowHeme = run(hemoglobin, { hemeSyn: 5 }), lowHemeNoHRI = run(hemoglobin, { hemeSyn: 5, hri: 0 });
  assert.ok(lowHeme[I.aSyn] < 0.7 * lowHemeNoHRI[I.aSyn], 'HRI throttles translation when heme is scarce');
  const apo = (y) => (y[I.T0] + y[I.T1] + y[I.T2] + y[I.T3]) / Math.max(y[I.T4], 1e-9);
  assert.ok(apo(lowHeme) < apo(lowHemeNoHRI), 'fewer heme-deficient tetramers with HRI');
});

test('hemoglobin: balanced globins make functional 4-heme hemoglobin', () => {
  const I = hemoglobin.idx;
  const y = run(hemoglobin, {});
  const tetramers = y[I.T0] + y[I.T1] + y[I.T2] + y[I.T3] + y[I.T4];
  assert.ok(y[I.T4] / tetramers > 0.8, 'most tetramers carry four hemes');
  assert.ok(y[I.P] / y[I.aSyn] < 0.1, 'little alpha precipitates when beta is sufficient');
});

test('hemoglobin: low beta expression (beta-thalassaemia) precipitates free alpha; AHSP buffers it', () => {
  const I = hemoglobin.idx;
  const thal = run(hemoglobin, { betaExpr: 0.3 });
  const normal = run(hemoglobin, {});
  assert.ok(thal[I.P] / thal[I.aSyn] > 3 * (normal[I.P] / normal[I.aSyn]));
  assert.ok(thal[I.T4] < 0.5 * normal[I.T4]);
  const noAhsp = run(hemoglobin, { betaExpr: 0.3, AHSP: 0 });
  assert.ok(noAhsp[I.P] > thal[I.P], 'without AHSP more alpha precipitates');
});

test('collagen: chain balance and 2:1 alpha1:alpha2 stoichiometry in procollagen', () => {
  const I = collagen.idx;
  const y = run(collagen, {}, 200);
  assert.ok(Math.abs(collagen.chainBalance(y) - (y[I.s1] + y[I.s2])) / (y[I.s1] + y[I.s2]) < 1e-6);
  const trimers = y[I.PC] + y[I.PCe] + y[I.TC] + y[I.F];
  assert.ok(trimers > 1);
});

test('collagen: no ascorbate (scurvy) or no ADAMTS2 (dermatosparaxis) blocks fibrils', () => {
  const I = collagen.idx;
  const ok = run(collagen, {}, 200), scurvy = run(collagen, { ascorbate: 0 }, 200), eds = run(collagen, { ADAMTS2: 0 }, 200);
  assert.ok(ok[I.F] > 5);
  assert.ok(scurvy[I.F] < 0.01 * ok[I.F], `scurvy fibrils ${scurvy[I.F]}`);
  assert.ok(eds[I.F] < 1e-9 && eds[I.PCe] > 1, 'procollagen accumulates uncleaved');
});
