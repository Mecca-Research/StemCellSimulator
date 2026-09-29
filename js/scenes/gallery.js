// Real micrographs, the measurements taken from them, and "animated
// variations": the real pixels of each specimen driven by the models, next to
// 3D reconstructions/models of the same process.
import * as THREE from 'three';
import { cellMaterial, lineMaterial, releaseTree } from '../engine/materials.js';
import { unitSphere, erythrocyteGeometry } from '../engine/geometry.js';
import { InstancePool } from '../engine/instances.js';
import { Panel, el } from '../engine/ui.js';
import { RNG } from '../models/core/rng.js';
import { fmt } from '../engine/chart.js';

const CATALOG_URL = new URL('../../assets/micrographs/catalog.json', import.meta.url);
const REFS_URL = new URL('../../assets/micrographs/references.json', import.meta.url);

// ---------------------------------------------------------------- sprites
// Instanced quads that each show a rectangle of a real micrograph through a
// soft circular mask: real nuclei/cells become movable agents.
const spriteVert = /* glsl */`
attribute vec4 aRect;
attribute float aAlpha;
varying vec2 vUv;
varying vec4 vRect;
varying float vA;
void main() {
  vUv = uv; vRect = aRect; vA = aAlpha;
  gl_Position = projectionMatrix * modelViewMatrix * instanceMatrix * vec4(position, 1.0);
}`;
const spriteFrag = /* glsl */`
uniform sampler2D map;
uniform float uSoft;
uniform vec3 uTint;
varying vec2 vUv;
varying vec4 vRect;
varying float vA;
void main() {
  float r = length(vUv - 0.5) * 2.0;
  float m = 1.0 - smoothstep(1.0 - uSoft, 1.0, r);
  vec4 c = texture2D(map, mix(vRect.xy, vRect.zw, vUv));
  gl_FragColor = vec4(c.rgb * uTint, m * vA);
  if (gl_FragColor.a < 0.01) discard;
  #include <colorspace_fragment>
}`;

class SpriteLayer {
  constructor(texture, capacity, { additive = false, soft = 0.35 } = {}) {
    const geo = new THREE.PlaneGeometry(1, 1);
    geo.rotateX(-Math.PI / 2);
    this.rect = new THREE.InstancedBufferAttribute(new Float32Array(capacity * 4), 4);
    this.alpha = new THREE.InstancedBufferAttribute(new Float32Array(capacity), 1);
    geo.setAttribute('aRect', this.rect);
    geo.setAttribute('aAlpha', this.alpha);
    this.material = new THREE.ShaderMaterial({
      vertexShader: spriteVert, fragmentShader: spriteFrag, transparent: true, depthWrite: false,
      blending: additive ? THREE.AdditiveBlending : THREE.NormalBlending,
      uniforms: { map: { value: texture }, uSoft: { value: soft }, uTint: { value: new THREE.Color(1, 1, 1) } },
    });
    this.mesh = new THREE.InstancedMesh(geo, this.material, capacity);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = 6;
    this.capacity = capacity;
    this.n = 0;
    this.m = new THREE.Matrix4();
    this.q = new THREE.Quaternion();
    this.s = new THREE.Vector3();
    this.p = new THREE.Vector3();
    this.yAxis = new THREE.Vector3(0, 1, 0);
  }
  begin() { this.n = 0; }
  /** rect = [u0, v0, u1, v1] in texture coordinates */
  put(x, y, z, sx, sz, angle, rect, alpha = 1) {
    if (this.n >= this.capacity) return;
    const i = this.n++;
    this.q.setFromAxisAngle(this.yAxis, angle);
    this.m.compose(this.p.set(x, y, z), this.q, this.s.set(sx, 1, sz));
    this.mesh.setMatrixAt(i, this.m);
    this.rect.setXYZW(i, rect[0], rect[1], rect[2], rect[3]);
    this.alpha.setX(i, alpha);
  }
  end() {
    this.mesh.count = this.n;
    this.mesh.instanceMatrix.needsUpdate = true;
    this.rect.needsUpdate = true;
    this.alpha.needsUpdate = true;
  }
}

function texFrom(img, { srgb = true, nearest = false } = {}) {
  const t = new THREE.Texture(img);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  if (nearest) t.minFilter = t.magFilter = THREE.NearestFilter;
  t.anisotropy = 4;
  t.needsUpdate = true;
  return t;
}

function imagePlane(texture, w, h, opacity = 1) {
  const g = new THREE.PlaneGeometry(w, h);
  g.rotateX(-Math.PI / 2);
  const m = new THREE.MeshBasicMaterial({ map: texture, transparent: opacity < 1, opacity, toneMapped: false, depthWrite: false });
  const mesh = new THREE.Mesh(g, m);
  mesh.renderOrder = 1;
  return mesh;
}

/** texture rect around pixel (x, y) with half-size s for an image W x H (y down). */
function rectAt(x, y, s, W, H) {
  return [(x - s) / W, 1 - (y + s) / H, (x + s) / W, 1 - (y - s) / H];
}

// ---------------------------------------------------------------- views
const VIEWS = {
  mitosis: {
    label: 'Human cells in mitosis → self-renewal (DNA stain)',
    async build({ group, sub, assets, stage, legend }) {
      const [meta, img, cal] = await Promise.all([assets.loadJSON('mitosis.json'), assets.loadImage('mitosis.png'), assets.loadJSON('calibration.json')]);
      const tex = texFrom(img);
      const [W, H] = meta.size;
      const k = 0.5; // world units per pixel
      const plane = imagePlane(tex, W * k, H * k, 0.0);
      group.add(plane);
      const rng = new RNG(5);
      const inter = meta.nuclei.filter((n) => !n.mitotic);
      const mito = meta.nuclei.filter((n) => n.mitotic);
      // every agent borrows the real pixels of an interphase or a mitotic nucleus
      const agents = meta.nuclei.map((n) => ({
        x: n.x, y: n.y, r: n.r, look: n, mitotic: n.mitotic, t: n.mitotic ? rng.uniform(0, 0.6) : 0,
        cycle: rng.uniform(0, 1), ang: n.orientation, vx: 0, vy: 0,
      }));
      const layer = new SpriteLayer(tex, 4000, { additive: true, soft: 0.45 });
      group.add(layer.mesh);
      const T = 1 / Math.max(cal.mitosis.mitotic_index, 1e-3); // cycle in units of t_M
      const params = { original: 0, speed: 1 };
      sub.section('Animated variation');
      sub.note(`Every nucleus below is a crop of the real image. Interphase nuclei drift; nuclei entering M take the appearance of a real condensed mitotic nucleus, separate along their spindle axis and re-form two interphase nuclei. Division rate: cycle T = t<sub>M</sub>/MI = ${T.toFixed(1)} t<sub>M</sub> from the measured mitotic index (${fmt(cal.mitosis.mitotic_index)}; ${cal.mitosis.n_mitotic} of ${cal.mitosis.n_nuclei} nuclei, rule: 90th-percentile DNA intensity > median + 4 MAD).`);
      sub.slider({ label: 'Overlay original image', min: 0, max: 1, value: 0, onChange: (v) => { plane.material.opacity = v; plane.material.transparent = true; } });
      const chart = sub.chart({ title: 'Nuclei and mitotic index vs time', xLabel: 't / t_M', series: [{ name: 'nuclei ÷ 100', color: '#5ee0c1' }, { name: 'mitotic index', color: '#ffb454' }, { name: 'measured MI', color: '#ff6fae', dash: [4, 3] }] });
      legend([{ color: '#ffffff', label: 'real interphase nuclei (crops)' }, { color: '#ffe07a', label: 'real mitotic figures (crops)' }]);
      stage.frame([0, 0, 0], W * k * 0.6, new THREE.Vector3(0, 1, 0.55));
      let t = 0, lastPush = -1;
      return {
        update(dt) {
          const dT = dt * 0.6; // 1 s = 0.6 t_M
          t += dT;
          const born = [];
          for (const a of agents) {
            if (a.dead) continue;
            if (!a.mitotic) {
              a.cycle += dT / T;
              if (a.cycle >= 1 && agents.length + born.length < 1800) { a.mitotic = true; a.t = 0; a.look = rng.pick(mito); a.ang = rng.uniform(0, Math.PI); }
              a.vx += rng.normal() * 0.4 * dT - a.vx * 0.5 * dT; a.vy += rng.normal() * 0.4 * dT - a.vy * 0.5 * dT;
              a.x += a.vx; a.y += a.vy;
            } else {
              a.t += dT;
              if (a.t >= 1) { // telophase done: two daughters
                a.dead = true;
                const d = a.r * 1.3;
                for (const sg of [1, -1]) {
                  born.push({ x: a.x + sg * d * Math.cos(a.ang + Math.PI / 2), y: a.y + sg * d * Math.sin(a.ang + Math.PI / 2), r: a.r * 0.85,
                    look: rng.pick(inter), mitotic: false, t: 0, cycle: rng.uniform(0, 0.1), ang: rng.uniform(0, 6.28), vx: 0, vy: 0 });
                }
              }
            }
          }
          if (born.length) { for (let i = agents.length - 1; i >= 0; i--) if (agents[i].dead) agents.splice(i, 1); agents.push(...born); }
          layer.begin();
          for (const a of agents) {
            const s = a.look.r * 1.7;
            let sep = 0;
            if (a.mitotic && a.t > 0.5) sep = (a.t - 0.5) / 0.5; // anaphase -> telophase
            const rect = rectAt(a.look.x, a.look.y, s, W, H);
            const sx = (a.x - W / 2) * k, sz = (a.y - H / 2) * k;
            if (sep > 0) {
              const d = sep * a.r * 1.1 * k;
              const ca = Math.cos(a.ang + Math.PI / 2), sa = Math.sin(a.ang + Math.PI / 2);
              for (const sg of [1, -1]) layer.put(sx + sg * d * ca, 0.1, sz + sg * d * sa, s * 2 * k * 0.8, s * 2 * k * 0.8, a.ang, rect, 1);
            } else layer.put(sx, 0.1, sz, s * 2 * k, s * 2 * k, a.mitotic ? a.ang : 0, rect, 1);
          }
          layer.end();
          if (t - lastPush > 0.05) {
            lastPush = t;
            const mi = agents.filter((a) => a.mitotic).length / agents.length;
            chart.push(t, [agents.length / 100, mi, cal.mitosis.mitotic_index]);
          }
          return `t = ${t.toFixed(2)} t_M · ${agents.length} nuclei`;
        },
        dispose() { tex.dispose(); },
      };
    },
  },

  npc: {
    label: 'Live imaging: protein moving to the nuclear envelope',
    async build({ group, sub, assets, stage, legend }) {
      const [kin, img] = await Promise.all([assets.loadJSON('npc_kinetics.json'), assets.loadImage('npc_timelapse.png')]);
      const tex = texFrom(img, { srgb: false });
      const [tw, th] = kin.tile;
      const cols = kin.cols;
      const rows = Math.ceil(kin.frames / cols);
      const k = 0.35;
      const frameMat = new THREE.ShaderMaterial({
        uniforms: { map: { value: tex }, uA: { value: new THREE.Vector4() }, uB: { value: new THREE.Vector4() }, uMix: { value: 0 }, uHalf: { value: 0 } },
        vertexShader: 'varying vec2 vUv; void main(){ vUv = uv; gl_Position = projectionMatrix * modelViewMatrix * vec4(position,1.0); }',
        fragmentShader: `uniform sampler2D map; uniform vec4 uA; uniform vec4 uB; uniform float uMix; varying vec2 vUv;
          void main(){ float a = texture2D(map, mix(uA.xy, uA.zw, vUv)).r; float b = texture2D(map, mix(uB.xy, uB.zw, vUv)).r;
          float v = mix(a, b, uMix); bool dna = vUv.x < 0.5; vec3 c = dna ? vec3(0.2, 0.55, 1.0) * v : vec3(0.25, 1.0, 0.45) * v;
          gl_FragColor = vec4(c * 1.4, 1.0); }`,
      });
      const g = new THREE.PlaneGeometry(tw * k, th * k); g.rotateX(-Math.PI / 2);
      const movie = new THREE.Mesh(g, frameMat);
      movie.position.set(-tw * k * 0.35, 0, 0);
      group.add(movie);
      const rectOf = (f) => { const c = f % cols, r = Math.floor(f / cols); return new THREE.Vector4(c / cols, 1 - (r + 1) / rows, (c + 1) / cols, 1 - r / rows); };
      // 3D nucleus model beside the movie: envelope enrichment follows the fitted kinetics
      const nuc = new THREE.Mesh(unitSphere(4), cellMaterial({ role: 'nucleus', color: 0x86a8ff, stain: 0x3a7bff }));
      nuc.scale.set(18, 14, 16);
      nuc.position.set(tw * k * 0.45, 16, 0);
      group.add(nuc);
      const envPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter' }), 2000, group);
      const cytoPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'reporter' }), 2000, group);
      const rng = new RNG(4);
      const envPts = Array.from({ length: 1400 }, () => rng.onSphere([0, 0, 0]));
      const cytoPts = Array.from({ length: 900 }, () => { const d = rng.onSphere([0, 0, 0]); const r = rng.uniform(1.25, 2.1); return d.map((v) => v * r); });
      const E = (t) => kin.fit.E_inf - (kin.fit.E_inf - kin.fit.E0) * Math.exp(-t / kin.fit.tau_frames);
      sub.section('Measured kinetics');
      sub.note(`Real 15-frame movie (Boni et al. 2015, CC0): left DNA, right the envelope-targeted protein. The rim/cytoplasm intensity ratio was measured on every frame and fitted with E(t) = E∞ − (E∞ − E₀)e^(−t/τ): E₀ = ${kin.fit.E0}, E∞ = ${kin.fit.E_inf}, τ = ${kin.fit.tau_frames} frames, R² = ${kin.fit.r2}. The 3D nucleus redistributes the same protein between cytoplasm and envelope with that fitted curve — the calibration used for nuclear translocation in the Stage-1 pathway scene.`);
      const chart = sub.chart({ title: 'Envelope enrichment: measured (dots) vs fit', xLabel: 'frame', series: [{ name: 'measured', color: '#5ee0c1', points: true, width: 0.01 }, { name: 'fit', color: '#ffb454' }] });
      const xs = kin.rim_over_cytoplasm.map((_, i) => i);
      chart.set(xs, [kin.rim_over_cytoplasm, xs.map(E)]);
      legend([{ color: '#3a8cff', label: 'DNA (real)' }, { color: '#40ff73', label: 'nuclear-envelope protein (real)' }, { color: '#9dff6a', label: 'model: envelope-bound' }, { color: '#ffd166', label: 'model: cytoplasmic' }]);
      stage.frame([0, 0, 0], tw * k * 0.75, new THREE.Vector3(0, 1, 0.8));
      let t = 0;
      const col = new THREE.Color();
      return {
        update(dt) {
          t = (t + dt * 2.2) % (kin.frames - 1 + 2);
          const f = Math.min(t, kin.frames - 1);
          const f0 = Math.floor(f), f1 = Math.min(f0 + 1, kin.frames - 1);
          frameMat.uniforms.uA.value.copy(rectOf(f0));
          frameMat.uniforms.uB.value.copy(rectOf(f1));
          frameMat.uniforms.uMix.value = f - f0;
          const frac = (E(f) - 1) / (kin.fit.E_inf - 1); // bound fraction proxy 0..1
          const nEnv = Math.round(envPts.length * THREE.MathUtils.clamp(frac, 0.05, 1));
          const nCyt = Math.round(cytoPts.length * THREE.MathUtils.clamp(1 - frac * 0.8, 0.1, 1));
          envPool.begin();
          for (let i = 0; i < nEnv; i++) { const d = envPts[i]; envPool.put(nuc.position.x + d[0] * 18.3, nuc.position.y + d[1] * 14.3, d[2] * 16.3, 0.45, 0.45, 0.45, col.set(0x9dff6a)); }
          envPool.end();
          cytoPool.begin();
          for (let i = 0; i < nCyt; i++) { const d = cytoPts[i]; cytoPool.put(nuc.position.x + d[0] * 18, nuc.position.y + d[1] * 14, d[2] * 16, 0.35, 0.35, 0.35, col.set(0xffd166)); }
          cytoPool.end();
          nuc.rotation.y += dt * 0.15;
          return `frame ${f.toFixed(1)} / ${kin.frames - 1} · envelope enrichment ${E(f).toFixed(2)}`;
        },
        dispose() { tex.dispose(); },
      };
    },
  },

  skin: {
    label: 'Skin H&E → 3D stratified epidermis (keratinocyte differentiation)',
    async build({ group, sub, assets, stage, legend }) {
      const [meta, plateImg, img] = await Promise.all([assets.loadJSON('skin.json'), assets.loadImage('skin_plate.jpg'), assets.loadImage('skin.jpg')]);
      const plateTex = texFrom(plateImg), tex = texFrom(img);
      const [W, H] = meta.size;
      const k = 0.18;
      const plate = imagePlane(plateTex, W * k, H * k);
      plate.position.x = -W * k * 0.55;
      group.add(plate);
      // surface line (top of epidermis) sampled every 16 px
      const top = meta.epidermis_top_px_every16;
      const topAt = (x) => { const i = Math.min(Math.max(Math.round(x / 16), 0), top.length - 1); return top[i] > 0 ? top[i] : 200; };
      const rng = new RNG(8);
      const layer = new SpriteLayer(tex, 3000, { soft: 0.5 });
      layer.mesh.position.x = plate.position.x;
      group.add(layer.mesh);
      const nuclei = meta.nuclei.map((n) => {
        const depth = n.y - topAt(n.x);
        return { ...n, x0: n.x, y0: n.y, epi: depth > -5 && depth < 170, depth, phase: rng.next(), look: n };
      });
      const basal = nuclei.filter((n) => n.epi && n.depth > 60);
      const params = { turnover: 1 };
      sub.section('Animated variation (real nuclei)');
      sub.note('Real hematoxylin-stained nuclei are cut out of the slide (the background plate has them removed). Epidermal nuclei move from the basal layer toward the surface, flatten and fade as keratinocytes cornify; basal cells divide to replace them. Dermal fibroblast nuclei stay nearly still. The speed is illustrative (human epidermal turnover takes ~4 weeks).');
      sub.slider({ label: 'Turnover speed', min: 0.2, max: 3, value: 1, onChange: (v) => { params.turnover = v; } });

      // ---- 3D stratified epidermis model to the right, same palette in H&E mode
      const model = new THREE.Group();
      model.position.set(W * k * 0.62, 0, 0);
      group.add(model);
      const sx = 120, sz = 70;
      const cellPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'cytoplasm', color: 0xf3b7c9, opacity: 0.55, he: [0.03, 0.5] }), 3000, model);
      const nucPool = new InstancePool(unitSphere(2), cellMaterial({ role: 'nucleus', color: 0x8f7bd9 }), 3000, model);
      const dermis = new THREE.Mesh(new THREE.BoxGeometry(sx, 24, sz), cellMaterial({ role: 'solid', color: 0xf0a3bd, he: [0.05, 0.42] }));
      dermis.position.y = -12;
      model.add(dermis);
      const fibers = [];
      for (let i = 0; i < 160; i++) { const p = [rng.uniform(-sx / 2, sx / 2), rng.uniform(-22, -2), rng.uniform(-sz / 2, sz / 2)]; const a = rng.uniform(0, Math.PI); fibers.push(p[0] - Math.cos(a) * 6, p[1], p[2] - Math.sin(a) * 6, p[0] + Math.cos(a) * 6, p[1] + rng.normal(), p[2] + Math.sin(a) * 6); }
      const fg = new THREE.BufferGeometry(); fg.setAttribute('position', new THREE.Float32BufferAttribute(fibers, 3));
      model.add(new THREE.LineSegments(fg, lineMaterial({ color: 0xffd7e4, opacity: 0.6 })));
      // Epidermal proliferative units: columns of keratinocytes over one basal
      // stem/progenitor cell; each basal division pushes the column up by one
      // cell, so cells traverse basal -> spinous -> granular -> cornified layers.
      const columns = [];
      const colW = 7.2, perCol = 9, spacing = 5.2;
      for (let ix = -sx / 2 + colW / 2; ix < sx / 2; ix += colW) {
        for (let iz = -sz / 2 + colW / 2; iz < sz / 2; iz += colW * 0.87) {
          columns.push({ x: ix + (Math.round(iz / colW) % 2) * colW * 0.5, z: iz, phase: rng.next(), rate: rng.uniform(0.8, 1.2) });
        }
      }
      const kc = columns;
      const layerOf = (h) => (h < 4 ? 'basal' : h < 24 ? 'spinous' : h < 34 ? 'granular' : 'cornified');
      const counts = sub.chart({ title: '3D model: cells per epidermal layer', xLabel: 'time', series: [
        { name: 'basal', color: '#8f7bd9' }, { name: 'spinous', color: '#f3b7c9' }, { name: 'granular', color: '#ffb454' }, { name: 'cornified', color: '#ff6fae' }] });
      legend([{ color: '#8f7bd9', label: 'nuclei (hematoxylin)' }, { color: '#f3b7c9', label: 'keratinocyte cytoplasm (eosin)' }, { color: '#ff6fae', label: 'cornified (anucleate) squames' }]);
      sub.note('Tip: switch to <b>H&amp;E</b> and turn on <b>section</b> in the transport bar to cut the 3D epidermis like a histology slide.');
      stage.frame([-8, 10, 0], 205, new THREE.Vector3(0, 0.95, 0.75));
      stage.clipAxis.set(0, 0, 1);
      stage.clipRange = [-40, 40];
      let t = 0, lastPush = -1;
      const col = new THREE.Color();
      return {
        update(dt) {
          t += dt * params.turnover;
          // real-pixel variation
          layer.begin();
          for (const n of nuclei) {
            let x = n.x0, y = n.y0, a = 1, flat = 1;
            if (n.epi) {
              const span = Math.max(n.depth + 10, 20);
              const u = ((t * 0.05 + n.phase) % 1);
              y = n.y0 - u * span; // toward the surface
              const f = u;
              flat = 1 + f * 1.4;
              a = f > 0.8 ? (1 - f) / 0.2 : Math.min(1, f / 0.05);
            } else {
              x += Math.sin(t * 0.3 + n.phase * 20) * 0.6;
            }
            const s = n.r * 2.1;
            layer.put((x - W / 2) * k, 0.2, (y - H / 2) * k, s * 2 * k * flat, s * 2 * k / flat, 0, rectAt(n.look.x, n.look.y, s, W, H), a);
          }
          layer.end();
          // 3D model: columns advance by one cell per basal division
          cellPool.begin(); nucPool.begin();
          const cnt = { basal: 0, spinous: 0, granular: 0, cornified: 0 };
          for (const c of kc) {
            c.phase += dt * params.turnover * 0.12 * c.rate;
            const f = c.phase % 1;
            for (let k = 0; k < perCol; k++) {
              const h = (k + f) * spacing - spacing * 0.5;
              if (h < 0) continue;
              const L = layerOf(h); cnt[L]++;
              const flatten = h < 4 ? 0 : Math.min((h - 4) / 30, 1);
              // basal: columnar; spinous: polyhedral; granular/cornified: flattened squames
              const w = h < 4 ? 2.6 : 3.4 + flatten * 3.6, hh = h < 4 ? 3.4 : 2.9 * (1 - flatten * 0.78);
              const shed = h > (perCol - 1.3) * spacing ? 1 - (h - (perCol - 1.3) * spacing) / (spacing * 0.8) : 1;
              if (shed <= 0) continue;
              cellPool.put(c.x, h + hh, c.z, w, hh * shed, w * 0.93, col.set(L === 'cornified' ? 0xff9fc3 : L === 'granular' ? 0xf6b3c6 : L === 'basal' ? 0xe99ab4 : 0xf3b7c9));
              if (L !== 'cornified') nucPool.put(c.x, h + hh, c.z, 1.7 * (1 - flatten * 0.45), 1.5 * (1 - flatten * 0.7), 1.7 * (1 - flatten * 0.45));
            }
          }
          cellPool.end(); nucPool.end();
          if (t - lastPush > 0.5) { lastPush = t; counts.push(t, [cnt.basal, cnt.spinous, cnt.granular, cnt.cornified]); }
          return `t = ${t.toFixed(1)} · switch to H&E mode to compare the 3D model with the real slide`;
        },
        dispose() { plateTex.dispose(); tex.dispose(); },
      };
    },
  },

  rbc: {
    label: 'Blood smear → 3D biconcave erythrocytes',
    async build({ group, sub, assets, stage, legend }) {
      const [rbc, img] = await Promise.all([assets.loadJSON('rbc.json'), assets.loadImage('rbc_normal.jpg')]);
      const tex = texFrom(img);
      const W = 640, H = 480, k = 0.25;
      const bgMat = new THREE.MeshBasicMaterial({ color: new THREE.Color().setRGB(0.55, 0.64, 0.68, THREE.SRGBColorSpace), toneMapped: false });
      const bg = new THREE.Mesh(new THREE.PlaneGeometry(W * k, H * k).rotateX(-Math.PI / 2), bgMat);
      bg.position.x = -W * k * 0.6;
      group.add(bg);
      const layer = new SpriteLayer(tex, 400, { soft: 0.3 });
      layer.mesh.position.x = bg.position.x;
      group.add(layer.mesh);
      // each sprite shows one of the real isolated erythrocytes measured by the pipeline
      const rng = new RNG(12);
      const dpx = rbc.normal.diameter_px_median;
      const crops = rbc.normal.centres_px;
      const cells = Array.from({ length: 150 }, () => ({ x: rng.uniform(0, W), y: rng.uniform(0, H), crop: rng.pick(crops), rot: rng.uniform(0, 6.28) }));
      // 3D erythrocytes in shear flow on the right
      const model = new THREE.Group();
      model.position.x = W * k * 0.55;
      group.add(model);
      const rbcGeo = erythrocyteGeometry(7.82, 48);
      const pool = new InstancePool(rbcGeo, cellMaterial({ role: 'rbc', color: 0xd9363e, stain: 0xff5533, noise: 0.15 }), 200, model);
      const flowCells = Array.from({ length: 120 }, () => ({ p: [rng.uniform(-40, 40), rng.uniform(-18, 18), rng.uniform(-25, 25)], q: new THREE.Quaternion().setFromEuler(new THREE.Euler(rng.uniform(0, 6), rng.uniform(0, 6), 0)), w: new THREE.Vector3(rng.normal(), rng.normal(), rng.normal()).normalize() }));
      const big = new THREE.Mesh(erythrocyteGeometry(7.82, 96), cellMaterial({ role: 'rbc', color: 0xd9363e }));
      big.scale.setScalar(3.2);
      big.position.set(0, 34, 0);
      model.add(big);
      sub.section('Measurement');
      sub.note(`The erythrocyte is its own ruler: median diameter ${fmt(rbc.normal.diameter_px_median)} px ≡ 7.82 µm (${fmt(rbc.normal.um_per_px)} µm/px, ${rbc.normal.n_cells} isolated normal cells). The measured absorbance profile across a cell (central pallor at ${fmt(rbc.normal.pallor_min_over_max)} of the rim maximum) is compared with the Evans–Fung thickness of the 3D model: correlation ${fmt(rbc.normal.profile_correlation)}, RMSE ${fmt(rbc.normal.profile_rmse)}. The lymphocyte in a second field measures ${fmt(rbc.lymphocyte.lymphocyte_nucleus_diameter_um)} µm (nucleus).`);
      const ch = sub.chart({ title: 'Radial profile: real absorbance vs 3D model thickness', xLabel: 'r / R', yRange: [0, 1.1], series: [{ name: 'measured (smear)', color: '#ff6fae', points: true }, { name: 'Evans–Fung 3D model', color: '#5ee0c1' }] });
      ch.set(rbc.normal.r_over_R, [rbc.normal.absorbance_profile, rbc.normal.evans_fung_thickness]);
      legend([{ color: '#d9363e', label: '3D erythrocytes (Evans–Fung biconcave disc)' }, { color: '#a8b8c2', label: 'real erythrocytes (crops) drifting in the smear plane' }]);
      stage.frame([-20, 10, 0], 165, new THREE.Vector3(0, 0.9, 0.9));
      const m = new THREE.Matrix4(), s = new THREE.Vector3(1, 1, 1), dq = new THREE.Quaternion(), pp = new THREE.Vector3();
      let t = 0;
      return {
        update(dt) {
          t += dt;
          layer.begin();
          for (const c of cells) {
            c.x = (c.x + dt * 8 * (0.6 + 0.4 * Math.sin(c.y / 60))) % W;
            c.rot += dt * 0.1;
            const r = dpx * 0.62;
            layer.put((c.x - W / 2) * k, 0.1, (c.y - H / 2) * k, r * 2 * k, r * 2 * k, c.rot, rectAt(c.crop[0], c.crop[1], r, W, H), 1);
          }
          layer.end();
          pool.begin();
          for (const c of flowCells) {
            c.p[0] += dt * (4 + c.p[1] * 0.35); // simple shear flow
            if (c.p[0] > 45) c.p[0] = -45;
            dq.setFromAxisAngle(c.w, dt * (0.4 + Math.abs(c.p[1]) * 0.05)); // tumbling
            c.q.multiply(dq);
            m.compose(pp.set(...c.p), c.q, s);
            pool.putMatrix(m);
          }
          pool.end();
          big.rotation.x = Math.sin(t * 0.4) * 0.9;
          big.rotation.z = t * 0.2;
          return `3D red cells tumbling in shear flow · diameter 7.82 µm`;
        },
        dispose() { tex.dispose(); },
      };
    },
  },

  organoids: {
    label: 'Intestinal organoids → 3D epithelial cysts',
    async build({ group, sub, assets, stage, legend }) {
      const [meta, img, plateImg] = await Promise.all([assets.loadJSON('organoids.json'), assets.loadImage('organoids.jpg'), assets.loadImage('organoids_plate.jpg')]);
      const tex = texFrom(img), plateTex = texFrom(plateImg);
      const [W, H] = meta.size;
      const k = 0.35;
      const plate = imagePlane(plateTex, W * k, H * k);
      group.add(plate);
      const layer = new SpriteLayer(tex, 400, { soft: 0.25 });
      group.add(layer.mesh);
      const rng = new RNG(21);
      const orgs = meta.circles.map((c) => ({ ...c, g: rng.uniform(0.6, 1.4), z: rng.uniform(4, 40) }));
      // 3D reconstruction hypothesis: each detected organoid as a hollow epithelial shell at a random depth in the gel
      const cellPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'membrane', color: 0xffd9a8, opacity: 0.6 }), 20000, group);
      const nucPool = new InstancePool(unitSphere(1), cellMaterial({ role: 'nucleus', color: 0x9db4ff }), 20000, group);
      const dirs = Array.from({ length: 220 }, () => rng.onSphere([0, 0, 0]));
      const params = { growth: 1, show3d: true };
      sub.section('Measured');
      sub.table(['quantity', 'value'], [['organoids detected (Hough)', meta.n_organoids], ['median radius (px)', meta.radius_px_median], ['radius CV', meta.radius_cv], ['p90 / median radius', meta.radius_p90_over_median]]);
      sub.note('Animated variation: each detected organoid (real crop) inflates its lumen with logistic growth; the 3D shells show the same organoids as epithelial cysts at an assumed depth in the gel (the image itself is 2D, so depth is illustrative).');
      sub.slider({ label: 'Growth speed', min: 0, max: 3, value: 1, onChange: (v) => { params.growth = v; } });
      sub.toggle({ label: 'Show 3D epithelial shells', value: true, onChange: (v) => { params.show3d = v; } });
      legend([{ color: '#ffd9a8', label: 'epithelial cells of each cyst (3D)' }, { color: '#9db4ff', label: 'nuclei' }]);
      stage.frame([0, 10, 0], W * k * 0.6, new THREE.Vector3(0, 1, 0.8));
      let t = 0;
      return {
        update(dt) {
          t += dt * params.growth;
          layer.begin();
          cellPool.begin(); nucPool.begin();
          for (const o of orgs) {
            const scale = 1 + 0.6 / (1 + Math.exp(-(t * 0.25 * o.g - 2))) - 0.6 / (1 + Math.exp(2));
            const r = o.r * scale;
            const x = (o.x - W / 2) * k, z = (o.y - H / 2) * k;
            layer.put(x, 0.2, z, r * 2.3 * k, r * 2.3 * k, 0, rectAt(o.x, o.y, o.r * 1.15, W, H), params.show3d ? 0.35 : 1);
            if (params.show3d && o.r > 12) {
              const R = r * k, n = Math.min(dirs.length, Math.round(40 * (r / 12) ** 2));
              for (let i = 0; i < n; i++) {
                const d = dirs[i];
                cellPool.put(x + d[0] * R, o.z + d[1] * R, z + d[2] * R, R * 0.2, R * 0.2, R * 0.2);
                nucPool.put(x + d[0] * R * 0.98, o.z + d[1] * R * 0.98, z + d[2] * R * 0.98, R * 0.09, R * 0.09, R * 0.09);
              }
            }
          }
          layer.end(); cellPool.end(); nucPool.end();
          return `lumen inflation · t = ${t.toFixed(1)}`;
        },
        dispose() { tex.dispose(); plateTex.dispose(); },
      };
    },
  },

  stack: {
    label: 'hiPSC confocal stack · z-scan through 60 optical sections',
    async build({ group, sub, assets, stage, legend }) {
      const [meta, img] = await Promise.all([assets.loadJSON('hipsc_cells.json'), assets.loadImage('hipsc_volume.jpg')]);
      const tex = texFrom(img, { srgb: false });
      const { cols, slices } = meta.mosaic;
      const rows = Math.ceil(slices / cols);
      const [ex, ey, ez] = meta.extent_xyz_um;
      const mk = (z, bright) => {
        const c = z % cols, r = Math.floor(z / cols);
        const mat = new THREE.ShaderMaterial({
          transparent: true, depthWrite: false, blending: THREE.AdditiveBlending, side: THREE.DoubleSide,
          uniforms: { map: { value: tex }, uRect: { value: new THREE.Vector4(c / cols, 1 - (r + 1) / rows, (c + 1) / cols, 1 - r / rows) }, uGain: { value: bright } },
          vertexShader: 'varying vec2 vUv; void main(){ vUv = uv; gl_Position = projectionMatrix * modelViewMatrix * vec4(position,1.0); }',
          fragmentShader: `uniform sampler2D map; uniform vec4 uRect; uniform float uGain; varying vec2 vUv;
            void main(){ vec3 s = texture2D(map, mix(uRect.xy, uRect.zw, vUv)).rgb;
            vec3 c = vec3(1.0, 0.16, 0.33) * s.r + vec3(0.17, 1.0, 0.48) * s.g; gl_FragColor = vec4(c * uGain, 1.0); }`,
        });
        const g = new THREE.PlaneGeometry(ex, ey); g.rotateX(-Math.PI / 2);
        const m = new THREE.Mesh(g, mat);
        m.position.y = z * meta.spacing_zyx_um[0] * 2;
        return m;
      };
      const ghosts = [];
      for (let z = 0; z < slices; z += 3) { const m = mk(z, 0.06); ghosts.push(m); group.add(m); }
      const current = mk(0, 1.0);
      group.add(current);
      sub.note('Each plane is one real optical section (membrane red, DNA green); the bright plane scans through the stack while faint planes show every third section (z drawn 2× for clarity).');
      legend([{ color: '#ff2a55', label: 'membrane dye' }, { color: '#2cff7a', label: 'DNA dye' }]);
      stage.frame([0, 20, 0], ex * 0.75, new THREE.Vector3(0.3, 0.7, 1));
      let t = 0;
      return {
        update(dt) {
          t += dt * 8;
          const z = Math.floor(t) % slices;
          const c = z % cols, r = Math.floor(z / cols);
          current.material.uniforms.uRect.value.set(c / cols, 1 - (r + 1) / rows, (c + 1) / cols, 1 - r / rows);
          current.position.y = z * meta.spacing_zyx_um[0] * 2;
          return `optical section ${z + 1} / ${slices} · z = ${(z * meta.spacing_zyx_um[0]).toFixed(1)} µm`;
        },
        dispose() { tex.dispose(); },
      };
    },
  },

  qpi: {
    label: 'Single cell, quantitative phase → 3D optical-thickness surface',
    async build({ group, sub, assets, stage, legend }) {
      const data = await assets.loadImageData('cell_qpi.png');
      const W = data.width, H = data.height, step = 3;
      const nx = Math.floor(W / step), ny = Math.floor(H / step);
      const g = new THREE.PlaneGeometry(W * 0.107, H * 0.107, nx - 1, ny - 1);
      g.rotateX(-Math.PI / 2);
      const pos = g.attributes.position;
      const colors = new Float32Array(pos.count * 3);
      const c = new THREE.Color();
      let lo = 255, hi = 0;
      for (let i = 0; i < data.data.length; i += 4) { lo = Math.min(lo, data.data[i]); hi = Math.max(hi, data.data[i]); }
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
        const v = (data.data[((j * step) * W + i * step) * 4] - lo) / (hi - lo);
        const idx = j * nx + i;
        pos.setY(idx, v * 9);
        c.setHSL(0.66 - 0.66 * v, 0.85, 0.35 + 0.3 * v);
        colors.set([c.r, c.g, c.b], idx * 3);
      }
      g.setAttribute('color', new THREE.BufferAttribute(colors, 3));
      g.computeVertexNormals();
      const surf = new THREE.Mesh(g, cellMaterial({ role: 'field', vertexColors: true, noise: 0 }));
      group.add(surf);
      sub.note('Quantitative phase φ ∝ (n<sub>cell</sub> − n<sub>medium</sub>)·thickness, so the phase map of this suspended cell (CC0, Müller et al. 2018) is shown as a 3D optical-thickness surface (height ∝ φ, arbitrary scale). Pixel size 0.107 µm.');
      legend([{ color: '#3355ff', label: 'low phase' }, { color: '#ff5533', label: 'high phase (thicker / denser)' }]);
      stage.frame([0, 3, 0], W * 0.107 * 0.6, new THREE.Vector3(0.4, 0.8, 1));
      return { update() { return 'drag to orbit the phase surface'; } };
    },
  },
};

export default {
  about: `
    <p>The reference microscopy used to build and calibrate the models, each shown
    with an <b>animated variation</b>: the real pixels of the specimen (nuclei, red cells,
    organoids cut out of the image) are set in motion by the corresponding model, and a
    3D model of the same process sits alongside for comparison.</p>
    <p>All redistributed images are CC0, public domain or MIT (see the table below and
    <code>ATTRIBUTION.md</code>). A further 46 openly-licensed references per process are
    catalogued in <code>assets/micrographs/references.json</code>.</p>`,
  paperRef: 'Paper: "Validation: cross-referencing results with existing data" and "Visualizing the Cellular Map".',

  async create(ctx) {
    const { root, panel, stage, assets, legend, setStatus } = ctx;
    stage.setBloomScale(0.25); // bright slides would otherwise bloom
    let group = new THREE.Group();
    root.add(group);
    const params = { view: VIEWS[ctx.query?.get('view')] ? ctx.query.get('view') : 'mitosis' };
    panel.section('Specimen');
    panel.select({
      label: 'Real micrograph', value: params.view,
      options: Object.entries(VIEWS).map(([value, v]) => ({ value, label: v.label })),
      onChange: (v) => { params.view = v; show(v); },
    });
    const subRoot = el('div');
    panel.current.append(subRoot);
    const sub = new Panel(subRoot);

    panel.section('Sources & licences (redistributed)');
    const [catalog, refs] = await Promise.all([
      fetch(CATALOG_URL).then((r) => r.json()),
      fetch(REFS_URL).then((r) => r.json()).catch(() => ({ entries: [] })),
    ]);
    panel.table(['specimen', 'licence'], catalog.entries.filter((e) => !e.id.endsWith('labels')).map((e) => [e.title, e.license]));
    panel.section(`Further references (${refs.entries.length}, not redistributed)`);
    const list = el('div', { class: 'note' });
    const byProc = {};
    for (const r of refs.entries) for (const p of r.process) (byProc[p] ??= []).push(r);
    for (const [proc, items] of Object.entries(byProc)) {
      list.append(el('div', { html: `<b>${proc.replace(/_/g, ' ')}</b>` }));
      for (const r of items.slice(0, 6)) list.append(el('div', { html: `· <a href="${r.file_page}" target="_blank" rel="noopener">${r.title.slice(0, 70)}</a> <span class="tag">${r.license}</span>` }));
    }
    panel.current.append(list);

    let active = null, token = 0;
    async function show(key) {
      const my = ++token;
      active?.dispose?.();
      active = null;
      root.remove(group);
      releaseTree(group);
      group = new THREE.Group();
      root.add(group);
      sub.clear();
      const inst = await VIEWS[key].build({ group, sub, assets, stage, legend });
      if (my !== token) { inst.dispose?.(); return; }
      active = inst;
    }
    await show(params.view);

    return {
      update(dt) {
        if (!active) return;
        const status = active.update(dt);
        if (status) setStatus(status);
        sub.drawCharts();
      },
      dispose() { active?.dispose?.(); },
    };
  },
};
