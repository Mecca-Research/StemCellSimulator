// Asset loaders with a small in-memory cache.
const cache = new Map();

function once(key, fn) {
  if (!cache.has(key)) cache.set(key, fn().catch((e) => { cache.delete(key); throw e; }));
  return cache.get(key);
}

export const ASSET_BASE = new URL('../../assets/derived/', import.meta.url);

export function assetURL(name) {
  return new URL(name, ASSET_BASE).href;
}

export function loadJSON(name) {
  return once(`json:${name}`, async () => {
    const r = await fetch(assetURL(name));
    if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
    return r.json();
  });
}

export function loadBinary(name) {
  return once(`bin:${name}`, async () => {
    const r = await fetch(assetURL(name));
    if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
    return r.arrayBuffer();
  });
}

/** Decode an image to raw RGBA pixels without colour management. */
export function loadImageData(name) {
  return once(`img:${name}`, async () => {
    const r = await fetch(assetURL(name));
    if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
    const blob = await r.blob();
    const bmp = await createImageBitmap(blob, { colorSpaceConversion: 'none', premultiplyAlpha: 'none' });
    const canvas = typeof OffscreenCanvas !== 'undefined'
      ? new OffscreenCanvas(bmp.width, bmp.height)
      : Object.assign(document.createElement('canvas'), { width: bmp.width, height: bmp.height });
    const ctx = canvas.getContext('2d', { willReadFrequently: true, colorSpace: 'srgb' });
    ctx.drawImage(bmp, 0, 0);
    const data = ctx.getImageData(0, 0, bmp.width, bmp.height);
    bmp.close?.();
    return data;
  });
}

/** Load an image as a THREE texture source (HTMLImageElement). */
export function loadImage(name) {
  return once(`el:${name}`, () => new Promise((resolve, reject) => {
    const im = new Image();
    im.onload = () => resolve(im);
    im.onerror = () => reject(new Error(`image ${name} failed`));
    im.src = assetURL(name);
  }));
}
