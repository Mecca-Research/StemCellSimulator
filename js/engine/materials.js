// Shared material system. Every biological surface uses one uber-shader
// ("cell material") with three render modes that the whole app switches
// together:
//
//   0 physical   wrap-lit translucent membranes, subsurface-ish nuclei
//   1 confocal   emission only: membrane dyes light up at grazing angles (as a
//                membrane stain does in an optical section), DNA dyes fill
//                the nucleus with chromatin texture, reporters glow
//   2 histology  virtual H&E: Beer-Lambert absorption of hematoxylin and eosin
//                (Giacomelli et al. 2016, PLoS ONE 11:e0159337 coefficients)
//
// Materials register themselves so a mode switch only flips uniforms and
// blending state; scenes never rebuild meshes.
import * as THREE from 'three';

export const MODES = { physical: 0, confocal: 1, histology: 2 };
const ROLE_IDS = { membrane: 0, cytoplasm: 1, nucleus: 2, reporter: 3, solid: 4, molecule: 5, rbc: 6, field: 7 };

// uniforms shared by reference across all materials
export const shared = {
  uMode: { value: 0 },
  uTime: { value: 0 },
  uLightDir: { value: new THREE.Vector3(0.45, 0.8, 0.55).normalize() },
};

const registry = new Set();
let currentMode = 'physical';

const vertexShader = /* glsl */`
#include <common>
#include <color_pars_vertex>
#include <clipping_planes_pars_vertex>
uniform float uNoiseScale;
varying vec3 vNormalV;
varying vec3 vViewPos;
varying vec3 vObjPos;
void main() {
  #include <color_vertex>
  vec3 transformed = position;
  vec3 objectNormal = normal;
  vec3 seed = vec3(0.0);
  #ifdef USE_INSTANCING
    mat3 im3 = mat3(instanceMatrix);
    objectNormal = inverse(transpose(im3)) * objectNormal;
    transformed = (instanceMatrix * vec4(transformed, 1.0)).xyz;
    seed = instanceMatrix[3].xyz * 0.37;
  #endif
  vec4 mvPosition = modelViewMatrix * vec4(transformed, 1.0);
  vNormalV = normalize(normalMatrix * objectNormal);
  vViewPos = -mvPosition.xyz;
  vObjPos = position * uNoiseScale + seed;
  gl_Position = projectionMatrix * mvPosition;
  #include <clipping_planes_vertex>
}`;

const fragmentShader = /* glsl */`
#include <common>
#include <color_pars_fragment>
#include <clipping_planes_pars_fragment>
uniform int uMode;
uniform int uRole;
uniform vec3 uColor;
uniform vec3 uStain;
uniform float uUseInst;
uniform float uOpacity;
uniform float uRimPower;
uniform float uNoise;
uniform float uEmissive;
uniform float uGain;
uniform vec2 uHE;
uniform vec3 uLightDir;
uniform float uTime;
varying vec3 vNormalV;
varying vec3 vViewPos;
varying vec3 vObjPos;

float hash3(vec3 p) { p = fract(p * 0.3183099 + 0.1); p *= 17.0; return fract(p.x * p.y * p.z * (p.x + p.y + p.z)); }
float vnoise(vec3 x) {
  vec3 i = floor(x); vec3 f = fract(x); f = f * f * (3.0 - 2.0 * f);
  return mix(mix(mix(hash3(i + vec3(0,0,0)), hash3(i + vec3(1,0,0)), f.x),
                 mix(hash3(i + vec3(0,1,0)), hash3(i + vec3(1,1,0)), f.x), f.y),
             mix(mix(hash3(i + vec3(0,0,1)), hash3(i + vec3(1,0,1)), f.x),
                 mix(hash3(i + vec3(0,1,1)), hash3(i + vec3(1,1,1)), f.x), f.y), f.z);
}
float fbm(vec3 p) { float a = 0.5, s = 0.0; for (int i = 0; i < 4; i++) { s += a * vnoise(p); p *= 2.03; a *= 0.5; } return s; }

void main() {
  #include <clipping_planes_fragment>
  vec3 N = normalize(vNormalV);
  if (!gl_FrontFacing) N = -N;
  vec3 V = normalize(vViewPos);
  float ndv = clamp(abs(dot(N, V)), 0.0, 1.0);
  float rim = pow(1.0 - ndv, uRimPower);
  vec3 base = uColor;
  #if defined( USE_COLOR ) || defined( USE_INSTANCING_COLOR )
    base *= vColor.rgb;
  #endif
  float n = fbm(vObjPos);
  bool membraneLike = (uRole == 0 || uRole == 1);
  vec3 col; float alpha;

  if (uMode == 0) {
    vec3 L = normalize(uLightDir);
    float wrap = membraneLike ? 0.6 : 0.3;
    float diff = clamp((dot(N, L) + wrap) / (1.0 + wrap), 0.0, 1.0);
    vec3 H = normalize(L + V);
    float spec = pow(max(dot(N, H), 0.0), uRole == 5 ? 64.0 : 40.0) * (uRole == 2 ? 0.18 : 0.4);
    float hemi = 0.5 + 0.5 * N.y;
    vec3 amb = mix(vec3(0.07, 0.08, 0.11), vec3(0.26, 0.28, 0.33), hemi);
    vec3 albedo = base * (1.0 - uNoise * 0.45 + uNoise * 0.6 * n);
    col = albedo * (amb + diff * vec3(1.0, 0.96, 0.9)) + spec;
    col += base * rim * (membraneLike ? 0.9 : 0.35);
    col += base * uEmissive;
    alpha = membraneLike ? clamp(uOpacity + rim * 0.75, 0.0, 1.0) : uOpacity;
  } else if (uMode == 1) {
    vec3 dye = mix(uStain, base, uUseInst);
    float s;
    if (membraneLike) s = rim * 1.15 + 0.015;
    else if (uRole == 2) s = (0.16 + 0.62 * n * n) * (0.45 + 0.55 * ndv);
    else if (uRole == 3 || uRole == 5) s = 0.75 + 0.35 * ndv;
    else if (uRole == 6) s = 0.1 + rim * 0.4;
    else if (uRole == 7) s = 0.9;
    else s = 0.06 + rim * 0.25;
    col = dye * s * uGain;
    alpha = 1.0;
  } else {
    // virtual H&E: transmitted light through hematoxylin (dH) and eosin (dE)
    vec3 bH = vec3(0.860, 1.000, 0.300);
    vec3 bE = vec3(0.050, 1.000, 0.544);
    float dH = uHE.x, dE = uHE.y;
    if (uRole == 2) { dH *= 0.75 + 0.7 * n; }
    else { dE *= 0.8 + 0.45 * n; }
    col = exp(-(bH * dH + bE * dE) * 2.3);
    float shade = 0.82 + 0.18 * clamp(dot(N, normalize(uLightDir)), 0.0, 1.0);
    col *= shade;
    if (uUseInst > 0.5) col = mix(col, col * base, 0.5);
    alpha = membraneLike ? clamp(0.45 + rim * 0.5, 0.0, 1.0) : 1.0;
  }
  gl_FragColor = vec4(col, alpha);
  #include <tonemapping_fragment>
  #include <colorspace_fragment>
}`;

const ROLE_DEFAULTS = {
  membrane: { stain: 0xff3b6b, opacity: 0.22, rimPower: 2.2, noise: 0.15, he: [0.04, 0.45] },
  cytoplasm: { stain: 0xff3b6b, opacity: 0.35, rimPower: 1.6, noise: 0.35, he: [0.05, 0.6] },
  nucleus: { stain: 0x39ff88, opacity: 1.0, rimPower: 2.0, noise: 0.8, he: [1.05, 0.15] },
  reporter: { stain: 0xffffff, opacity: 1.0, rimPower: 2.0, noise: 0.1, he: [0.1, 0.4], useInst: 1 },
  solid: { stain: 0x3355aa, opacity: 1.0, rimPower: 2.5, noise: 0.3, he: [0.08, 0.45] },
  molecule: { stain: 0xffffff, opacity: 1.0, rimPower: 2.5, noise: 0.05, he: [0.2, 0.4], useInst: 1 },
  rbc: { stain: 0xff5533, opacity: 1.0, rimPower: 2.0, noise: 0.1, he: [0.0, 1.35] },
  field: { stain: 0xffffff, opacity: 1.0, rimPower: 2.0, noise: 0.0, he: [0.3, 0.3], useInst: 1 },
};

/**
 * Create a registered cell material.
 * @param {object} o
 * @param {keyof ROLE_DEFAULTS} o.role
 * @param {number|THREE.Color} [o.color]    albedo / tint (multiplied with instance or vertex colors)
 * @param {number} [o.stain]                confocal dye colour (ignored when useInst = 1)
 * @param {number} [o.useInst]              1 = use instance/vertex colour as dye colour in confocal
 * @param {boolean} [o.vertexColors]
 */
export function cellMaterial(o = {}) {
  const role = o.role ?? 'membrane';
  const d = ROLE_DEFAULTS[role];
  const mat = new THREE.ShaderMaterial({
    vertexShader,
    fragmentShader,
    clipping: true,
    vertexColors: !!o.vertexColors,
    side: o.side ?? (role === 'membrane' || role === 'cytoplasm' ? THREE.DoubleSide : THREE.FrontSide),
    uniforms: {
      uMode: shared.uMode,
      uTime: shared.uTime,
      uLightDir: shared.uLightDir,
      uRole: { value: ROLE_IDS[role] },
      uColor: { value: new THREE.Color(o.color ?? 0xffffff) },
      uStain: { value: new THREE.Color(o.stain ?? d.stain) },
      uUseInst: { value: o.useInst ?? d.useInst ?? 0 },
      uOpacity: { value: o.opacity ?? d.opacity },
      uRimPower: { value: o.rimPower ?? d.rimPower },
      uNoise: { value: o.noise ?? d.noise },
      uNoiseScale: { value: o.noiseScale ?? 0.45 },
      uEmissive: { value: o.emissive ?? 0 },
      uGain: { value: o.gain ?? 1 },
      uHE: { value: new THREE.Vector2(...(o.he ?? d.he)) },
    },
  });
  mat.userData.role = role;
  mat.userData.physicalOpacity = mat.uniforms.uOpacity.value;
  registry.add(mat);
  applyModeTo(mat, currentMode);
  return mat;
}

/** Registered line material that stays legible on every background. */
export function lineMaterial(o = {}) {
  const mat = new THREE.LineBasicMaterial({
    color: o.color ?? 0xffffff, vertexColors: !!o.vertexColors, transparent: true,
    opacity: o.opacity ?? 0.8, depthWrite: false,
  });
  mat.userData.role = 'line';
  mat.userData.baseColor = new THREE.Color(o.color ?? 0xffffff);
  mat.userData.baseOpacity = o.opacity ?? 0.8;
  registry.add(mat);
  applyModeTo(mat, currentMode);
  return mat;
}

function applyModeTo(mat, mode) {
  const role = mat.userData.role;
  if (role === 'line') {
    mat.blending = mode === 'confocal' ? THREE.AdditiveBlending : THREE.NormalBlending;
    mat.color.copy(mat.userData.baseColor);
    if (mode === 'histology') mat.color.multiplyScalar(0.55);
    mat.opacity = mat.userData.baseOpacity;
    mat.needsUpdate = true;
    return;
  }
  const membraneLike = role === 'membrane' || role === 'cytoplasm';
  const baseOpacity = mat.userData.physicalOpacity;
  if (mode === 'confocal') {
    mat.transparent = true;
    mat.blending = THREE.AdditiveBlending;
    mat.depthWrite = false;
  } else if (mode === 'histology') {
    mat.transparent = membraneLike;
    mat.blending = THREE.NormalBlending;
    mat.depthWrite = !membraneLike;
  } else {
    const translucent = membraneLike || baseOpacity < 1;
    mat.transparent = translucent;
    mat.blending = THREE.NormalBlending;
    mat.depthWrite = !translucent;
  }
  mat.needsUpdate = true;
}

export function setMode(mode) {
  currentMode = mode;
  shared.uMode.value = MODES[mode];
  for (const m of registry) applyModeTo(m, mode);
}

export function getMode() { return currentMode; }

export function releaseMaterial(mat) {
  registry.delete(mat);
  mat.dispose();
}

/** Dispose every registered material that belongs to a scene root. */
export function releaseTree(root) {
  root.traverse((obj) => {
    if (obj.geometry) obj.geometry.dispose();
    const mats = Array.isArray(obj.material) ? obj.material : obj.material ? [obj.material] : [];
    for (const m of mats) {
      if (registry.has(m)) releaseMaterial(m);
      else m.dispose?.();
    }
  });
}
