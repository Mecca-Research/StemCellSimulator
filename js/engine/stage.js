// Renderer, camera, controls, post-processing, global section plane, picking
// and the fixed-timestep simulation loop shared by all scenes.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { shared, setMode, getMode } from './materials.js';

const BACKGROUNDS = {
  physical: new THREE.Color(0x05070b),
  confocal: new THREE.Color(0x000000),
  histology: new THREE.Color(0xf4ecf1),
};
const BLOOM = { physical: 0.12, confocal: 0.5, histology: 0.0 };
const BLOOM_THRESHOLD = { physical: 0.6, confocal: 0.28, histology: 1 };

export class Stage {
  constructor(canvas) {
    this.canvas = canvas;
    const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true, powerPreference: 'high-performance' });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.05;
    renderer.localClippingEnabled = true;
    this.renderer = renderer;

    this.scene = new THREE.Scene();
    this.scene.background = BACKGROUNDS.physical.clone();
    this.camera = new THREE.PerspectiveCamera(42, 1, 0.05, 20000);
    this.camera.position.set(0, 40, 120);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;

    // three.js lights for the few built-in materials (lines, points, standard)
    const hemi = new THREE.HemisphereLight(0xdde8ff, 0x1a1410, 0.9);
    const dir = new THREE.DirectionalLight(0xffffff, 1.6);
    dir.position.set(40, 80, 60);
    this.scene.add(hemi, dir);

    this.root = new THREE.Group();
    this.scene.add(this.root);

    this.clipPlane = new THREE.Plane(new THREE.Vector3(0, 0, -1), 0);
    this.clipOn = false;
    this.clipAxis = new THREE.Vector3(0, 0, 1); // scenes may point this elsewhere
    this.clipRange = [-50, 50];

    const size = new THREE.Vector2(1, 1);
    this.composer = new EffectComposer(renderer);
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.bloom = new UnrealBloomPass(size, BLOOM.physical, 0.5, BLOOM_THRESHOLD.physical);
    this.composer.addPass(this.bloom);
    this.composer.addPass(new OutputPass());

    this.bloomScale = 1; // scenes showing bright transmitted-light images turn bloom down
    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this.pickables = [];
    this.modeListeners = new Set();

    this.running = true;
    this.speed = 1;
    this.simTime = 0;
    this.frames = 0;
    this.clock = new THREE.Clock();
    this.onFrame = null;
    // adaptive resolution: keep interaction smooth on weak GPUs / software GL
    this.maxPixelRatio = Math.min(window.devicePixelRatio || 1, 2);
    this.pixelRatio = this.maxPixelRatio;
    this._frameAcc = 0;
    this._frameN = 0;

    this._resize = () => this.resize();
    window.addEventListener('resize', this._resize);
    new ResizeObserver(this._resize).observe(canvas);
    this.resize();
  }

  resize() {
    const w = this.canvas.clientWidth || 1, h = this.canvas.clientHeight || 1;
    this.renderer.setSize(w, h, false);
    this.composer.setSize(w, h);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }

  setMode(mode) {
    setMode(mode);
    this.scene.background = BACKGROUNDS[mode].clone();
    this.bloom.strength = BLOOM[mode] * this.bloomScale;
    this.bloom.threshold = BLOOM_THRESHOLD[mode];
    this.renderer.toneMapping = mode === 'histology' ? THREE.NoToneMapping : THREE.ACESFilmicToneMapping;
    for (const fn of this.modeListeners) fn(mode);
  }

  get mode() { return getMode(); }

  setBloomScale(s) {
    this.bloomScale = s;
    this.bloom.strength = BLOOM[getMode()] * s;
  }

  /** Frame the camera on a bounding sphere. */
  frame(center, radius, dir = new THREE.Vector3(0.35, 0.55, 1)) {
    const c = new THREE.Vector3(...(Array.isArray(center) ? center : [center.x, center.y, center.z]));
    const dist = radius / Math.sin(THREE.MathUtils.degToRad(this.camera.fov / 2)) * 1.05;
    this.camera.position.copy(c).add(dir.clone().normalize().multiplyScalar(dist));
    this.camera.near = Math.max(radius / 200, 0.01);
    this.camera.far = dist * 20 + radius * 10;
    this.camera.updateProjectionMatrix();
    this.controls.target.copy(c);
    this.controls.update();
    this.clipRange = [-radius, radius];
    this.clipCenter = c;
  }

  /** Configure the global section plane: axis in world space, t in [-1, 1]. */
  setClip(on, t = 0) {
    this.clipOn = on;
    const c = this.clipCenter ?? new THREE.Vector3();
    const [a, b] = this.clipRange;
    const offset = a + (b - a) * (t + 1) / 2;
    const axis = this.clipAxis.clone().normalize();
    // keep points p with dot(axis, p - p0) <= 0, i.e. everything below the section
    const p0 = c.clone().add(axis.clone().multiplyScalar(offset));
    this.clipPlane.setFromNormalAndCoplanarPoint(axis.negate(), p0);
    this.renderer.clippingPlanes = on ? [this.clipPlane] : [];
    this.root.traverse((o) => o.syncClip?.(this.clipPlane, on));
  }

  setPickables(list) { this.pickables = list; }

  pick(clientX, clientY) {
    if (!this.pickables.length) return null;
    const r = this.canvas.getBoundingClientRect();
    this.pointer.set(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    this.raycaster.setFromCamera(this.pointer, this.camera);
    const hits = this.raycaster.intersectObjects(this.pickables.map((p) => p.object), false);
    for (const h of hits) {
      if (this.clipOn && this.clipPlane.distanceToPoint(h.point) < 0) continue;
      const entry = this.pickables.find((p) => p.object === h.object);
      if (entry) {
        const info = entry.info(h);
        if (info) return info;
      }
    }
    return null;
  }

  start() {
    const tick = () => {
      this._raf = requestAnimationFrame(tick);
      const raw = this.clock.getDelta();
      const dt = Math.min(raw, 1 / 20);
      this.adaptResolution(raw);
      shared.uTime.value += dt;
      if (this.onFrame) this.onFrame(this.running ? dt * this.speed : 0, dt);
      this.controls.update();
      if (this.bloom.strength > 0.001) this.composer.render();
      else this.renderer.render(this.scene, this.camera);
      this.frames++;
    };
    tick();
  }

  adaptResolution(raw) {
    this._frameAcc += raw;
    if (++this._frameN < 20) return;
    const mean = this._frameAcc / this._frameN;
    this._frameAcc = 0; this._frameN = 0;
    let pr = this.pixelRatio;
    if (mean > 1 / 24 && pr > 0.35) pr = Math.max(0.35, pr * 0.8);
    else if (mean < 1 / 50 && pr < this.maxPixelRatio) pr = Math.min(this.maxPixelRatio, pr * 1.15);
    if (pr !== this.pixelRatio) {
      this.pixelRatio = pr;
      this.renderer.setPixelRatio(pr);
      this.composer.setPixelRatio?.(pr);
      this.resize();
    }
  }

  screenshot() {
    return this.canvas.toDataURL('image/png');
  }
}
