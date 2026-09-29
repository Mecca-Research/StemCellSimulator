// Ray-marched volume renderer for real confocal z-stacks.
// Supports maximum-intensity projection and front-to-back emission/absorption
// compositing, per-channel transfer windows, a segmentation-label texture for
// highlighting reconstructed cells, the global section plane, and a virtual
// H&E mode (Beer-Lambert with hematoxylin <- DNA, eosin <- membrane/cytoplasm).
import * as THREE from 'three';
import { shared } from './materials.js';

const vert = /* glsl */`
varying vec3 vPos;
varying vec3 vCamObj;
void main() {
  vPos = position;
  vCamObj = (inverse(modelMatrix) * vec4(cameraPosition, 1.0)).xyz;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const frag = /* glsl */`
precision highp float;
precision highp sampler3D;
uniform sampler3D uVol;
uniform sampler3D uLab;
uniform int uHasLab;
uniform vec3 uExtent;
uniform vec3 uC0; uniform vec3 uC1; uniform vec3 uC2;
uniform vec3 uWinLo; uniform vec3 uWinHi;
uniform int uStyle;           // 0 MIP, 1 composite
uniform float uDensity;
uniform int uSteps;
uniform float uHighlight;     // label id to highlight (0 none)
uniform float uLabelTint;     // tint strength for labelled cells
uniform int uMode;
uniform vec4 uClip;           // world-space plane (n, d)
uniform int uClipOn;
uniform float uOpacity;
uniform mat4 modelMatrix;   // three.js only predeclares it for vertex shaders
varying vec3 vPos;
varying vec3 vCamObj;

vec2 hitBox(vec3 o, vec3 d) {
  vec3 inv = 1.0 / d;
  vec3 t0 = (vec3(0.0) - o) * inv, t1 = (uExtent - o) * inv;
  vec3 tmin = min(t0, t1), tmax = max(t0, t1);
  return vec2(max(max(tmin.x, tmin.y), tmin.z), min(min(tmax.x, tmax.y), tmax.z));
}
vec3 win(vec3 v) { return clamp((v - uWinLo) / max(uWinHi - uWinLo, vec3(1e-4)), 0.0, 1.0); }
vec3 hue(float id) { return 0.55 + 0.45 * cos(6.2831 * (fract(id * 0.61803) + vec3(0.0, 0.33, 0.67))); }

void main() {
  vec3 dir = normalize(vPos - vCamObj);
  vec2 t = hitBox(vCamObj, dir);
  t.x = max(t.x, 0.0);
  if (t.x >= t.y) discard;
  float len = t.y - t.x;
  float dt = len / float(uSteps);
  // jitter start to hide slicing artefacts
  float j = fract(sin(dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
  vec3 p = vCamObj + dir * (t.x + dt * j);
  vec3 acc = vec3(0.0); float aacc = 0.0; vec3 mip = vec3(0.0);
  vec3 odH = vec3(0.0); float nIn = 0.0;
  float stepScale = dt / 0.5;
  for (int i = 0; i < 512; i++) {
    if (i >= uSteps) break;
    vec3 uvw = p / uExtent;
    bool inside = true;
    if (uClipOn == 1) {
      vec3 wp = (modelMatrix * vec4(p, 1.0)).xyz;
      inside = dot(uClip.xyz, wp) + uClip.w >= 0.0;
    }
    if (inside) {
      vec3 v = win(texture(uVol, uvw).rgb);
      float lab = 0.0;
      if (uHasLab == 1) lab = texture(uLab, uvw).r * 255.0;
      float cid = mod(lab, 128.0);
      if (uMode == 2) {
        // virtual H&E: stain densities averaged along the ray, i.e. the look
        // of a thin physical section regardless of view angle or slab depth
        odH += vec3(0.860, 1.0, 0.300) * v.g * 1.6 + vec3(0.050, 1.0, 0.544) * (0.25 + v.r) * 0.6;
        nIn += 1.0;
      } else {
        vec3 c = uC0 * v.r + uC1 * v.g + uC2 * v.b;
        if (uHasLab == 1 && cid > 0.5) {
          vec3 tint = hue(cid);
          float hl = (uHighlight > 0.5 && abs(cid - uHighlight) < 0.5) ? 1.0 : 0.0;
          c = mix(c, c * tint * 1.6 + tint * 0.05 * max(v.r, v.g), uLabelTint);
          c += hl * vec3(0.9, 0.9, 0.3) * (0.15 + v.g);
        }
        if (uStyle == 0) {
          mip = max(mip, c);
        } else {
          float a = clamp(max(max(v.r, v.g), v.b) * uDensity * stepScale * 0.12, 0.0, 1.0);
          acc += (1.0 - aacc) * c * a;
          aacc += (1.0 - aacc) * a;
          if (aacc > 0.985) break;
        }
      }
    }
    p += dir * dt;
  }
  vec4 outc;
  if (uMode == 2) {
    vec3 T = exp(-(odH / max(nIn, 1.0)) * uDensity * 1.6);
    outc = vec4(T, 1.0 - min(min(T.r, T.g), T.b) * 0.15);
    outc.rgb = T; outc.a = uOpacity;
  } else if (uStyle == 0) {
    outc = vec4(mip, max(max(mip.r, mip.g), mip.b) * uOpacity);
  } else {
    outc = vec4(acc, aacc * uOpacity);
  }
  gl_FragColor = outc;
  #include <colorspace_fragment>
}`;

/**
 * Build a 3D texture from a slice mosaic image.
 * @param {ImageData} img mosaic image (tiles laid out row-major, `cols` per row)
 */
export function mosaicToVolume(img, { cols, slices, tile: [w, h] }, channels = 3) {
  const data = new Uint8Array(w * h * slices * 4);
  const src = img.data;
  const W = img.width;
  for (let z = 0; z < slices; z++) {
    const ox = (z % cols) * w, oy = Math.floor(z / cols) * h;
    for (let y = 0; y < h; y++) {
      let si = ((oy + y) * W + ox) * 4;
      let di = (z * w * h + y * w) * 4;
      for (let x = 0; x < w; x++, si += 4, di += 4) {
        data[di] = src[si];
        data[di + 1] = channels > 1 ? src[si + 1] : 0;
        data[di + 2] = channels > 2 ? src[si + 2] : 0;
        data[di + 3] = 255;
      }
    }
  }
  const tex = new THREE.Data3DTexture(data, w, h, slices);
  tex.format = THREE.RGBAFormat;
  tex.type = THREE.UnsignedByteType;
  tex.minFilter = tex.magFilter = THREE.LinearFilter;
  tex.unpackAlignment = 1;
  tex.needsUpdate = true;
  return tex;
}

/** Single-channel label volume (nearest sampling so ids never blend). */
export function mosaicToLabels(img, { cols, slices, tile: [w, h] }) {
  const data = new Uint8Array(w * h * slices);
  const src = img.data;
  const W = img.width;
  for (let z = 0; z < slices; z++) {
    const ox = (z % cols) * w, oy = Math.floor(z / cols) * h;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) data[z * w * h + y * w + x] = src[((oy + y) * W + ox + x) * 4];
    }
  }
  const tex = new THREE.Data3DTexture(data, w, h, slices);
  tex.format = THREE.RedFormat;
  tex.type = THREE.UnsignedByteType;
  tex.minFilter = tex.magFilter = THREE.NearestFilter;
  tex.unpackAlignment = 1;
  tex.needsUpdate = true;
  return { tex, data, dims: [w, h, slices] };
}

export class VolumeView extends THREE.Mesh {
  /**
   * @param {THREE.Data3DTexture} vol
   * @param {number[]} extent xyz extent in scene units (um)
   */
  constructor(vol, extent, opts = {}) {
    const [ex, ey, ez] = extent;
    const geo = new THREE.BoxGeometry(ex, ey, ez);
    geo.translate(ex / 2, ey / 2, ez / 2);
    const labTex = opts.labels ?? new THREE.Data3DTexture(new Uint8Array(1), 1, 1, 1);
    if (!opts.labels) { labTex.format = THREE.RedFormat; labTex.needsUpdate = true; }
    const mat = new THREE.ShaderMaterial({
      vertexShader: vert,
      fragmentShader: frag,
      side: THREE.BackSide,
      transparent: true,
      depthWrite: false,
      uniforms: {
        uVol: { value: vol },
        uLab: { value: labTex },
        uHasLab: { value: opts.labels ? 1 : 0 },
        uExtent: { value: new THREE.Vector3(ex, ey, ez) },
        uC0: { value: new THREE.Color(opts.c0 ?? 0xff2a55) },
        uC1: { value: new THREE.Color(opts.c1 ?? 0x2cff7a) },
        uC2: { value: new THREE.Color(opts.c2 ?? 0x3a7bff) },
        uWinLo: { value: new THREE.Vector3(0.08, 0.08, 0.08) },
        uWinHi: { value: new THREE.Vector3(0.9, 0.9, 0.9) },
        uStyle: { value: opts.style ?? 0 },
        uDensity: { value: opts.density ?? 1.0 },
        uSteps: { value: opts.steps ?? 220 },
        uHighlight: { value: 0 },
        uLabelTint: { value: 0 },
        uMode: shared.uMode,
        uClip: { value: new THREE.Vector4(0, 0, 1, 0) },
        uClipOn: { value: 0 },
        uOpacity: { value: 1 },
      },
    });
    super(geo, mat);
    this.extent = extent;
    this.renderOrder = 10;
    this.frustumCulled = false;
  }

  /** Mirror the renderer's global section plane into the shader. */
  syncClip(plane, on) {
    const u = this.material.uniforms;
    u.uClipOn.value = on ? 1 : 0;
    // three.js keeps points with distanceToPoint >= 0; the shader uses the same rule
    if (on) u.uClip.value.set(plane.normal.x, plane.normal.y, plane.normal.z, plane.constant);
  }

  setWindow(lo, hi) {
    this.material.uniforms.uWinLo.value.setScalar(lo);
    this.material.uniforms.uWinHi.value.setScalar(hi);
  }
}
