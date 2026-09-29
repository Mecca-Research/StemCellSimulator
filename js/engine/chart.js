// Lightweight streaming line chart on a 2D canvas (no dependencies).
export class Chart {
  /**
   * @param {HTMLCanvasElement} canvas
   * @param {object} o
   * @param {{name:string,color:string,dash?:number[]}[]} o.series
   * @param {number} [o.capacity=400] samples kept per series
   * @param {[number,number]|null} [o.yRange] fixed y range (null = auto)
   * @param {boolean} [o.logY]
   * @param {string} [o.xLabel] [o.yLabel]
   */
  constructor(canvas, o) {
    this.c = canvas;
    this.ctx = canvas.getContext('2d');
    this.series = o.series.map((s) => ({ ...s, data: [] }));
    this.capacity = o.capacity ?? 400;
    this.yRange = o.yRange ?? null;
    this.logY = !!o.logY;
    this.xLabel = o.xLabel ?? '';
    this.yLabel = o.yLabel ?? '';
    this.markers = [];
    this.xs = [];
    this.dirty = true;
  }

  push(x, values) {
    this.xs.push(x);
    values.forEach((v, i) => this.series[i]?.data.push(v));
    if (this.xs.length > this.capacity) {
      this.xs.shift();
      for (const s of this.series) s.data.shift();
    }
    this.dirty = true;
  }

  /** Replace all data (for static curves such as dispersion relations). */
  set(xs, columns) {
    this.xs = xs.slice();
    columns.forEach((col, i) => { if (this.series[i]) this.series[i].data = col.slice(); });
    this.dirty = true;
  }

  clear() {
    this.xs = [];
    for (const s of this.series) s.data = [];
    this.markers = [];
    this.dirty = true;
  }

  draw() {
    if (!this.dirty) return;
    this.dirty = false;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const W = this.c.clientWidth, H = this.c.clientHeight;
    if (!W || !H) return;
    if (this.c.width !== Math.round(W * dpr) || this.c.height !== Math.round(H * dpr)) {
      this.c.width = Math.round(W * dpr);
      this.c.height = Math.round(H * dpr);
    }
    const g = this.ctx;
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    const padL = 38, padR = 8, padT = 8, padB = 30;
    const pw = W - padL - padR, ph = H - padT - padB;
    const xs = this.xs;
    if (!xs.length) { this.legend(g, padL, H); return; }
    const x0 = xs[0], x1 = xs[xs.length - 1] === x0 ? x0 + 1 : xs[xs.length - 1];
    let y0, y1;
    const tf = (v) => (this.logY ? Math.log10(Math.max(v, 1e-6)) : v);
    if (this.yRange) [y0, y1] = this.yRange.map(tf);
    else {
      y0 = Infinity; y1 = -Infinity;
      for (const s of this.series) for (const v of s.data) { if (Number.isFinite(v)) { const t = tf(v); if (t < y0) y0 = t; if (t > y1) y1 = t; } }
      if (!Number.isFinite(y0)) { y0 = 0; y1 = 1; }
      if (y1 - y0 < 1e-9) { y1 += 0.5; y0 -= 0.5; }
      const pad = (y1 - y0) * 0.08; y0 -= pad; y1 += pad;
    }
    const X = (x) => padL + ((x - x0) / (x1 - x0)) * pw;
    const Y = (v) => padT + (1 - (tf(v) - y0) / (y1 - y0)) * ph;

    // grid + axes
    g.strokeStyle = '#1d2733'; g.lineWidth = 1; g.fillStyle = '#7d8ca0'; g.font = '10px ui-monospace, monospace';
    for (let i = 0; i <= 4; i++) {
      const v = y0 + ((y1 - y0) * i) / 4;
      const py = padT + (1 - i / 4) * ph;
      g.beginPath(); g.moveTo(padL, py); g.lineTo(padL + pw, py); g.stroke();
      const lab = this.logY ? `1e${v.toFixed(1)}` : fmt(v);
      g.textAlign = 'right'; g.fillText(lab, padL - 4, py + 3);
    }
    g.textAlign = 'center';
    for (let i = 0; i <= 4; i++) {
      const v = x0 + ((x1 - x0) * i) / 4;
      g.fillText(fmt(v), padL + (pw * i) / 4, padT + ph + 12);
    }
    if (this.xLabel) { g.textAlign = 'right'; g.fillText(this.xLabel, padL + pw, H - 3); }

    // zero line
    if (!this.logY && y0 < 0 && y1 > 0) {
      g.strokeStyle = '#3a4a61'; g.beginPath(); g.moveTo(padL, Y(0)); g.lineTo(padL + pw, Y(0)); g.stroke();
    }
    for (const m of this.markers) {
      g.strokeStyle = m.color ?? '#ffb454'; g.setLineDash([3, 3]);
      g.beginPath(); g.moveTo(X(m.x), padT); g.lineTo(X(m.x), padT + ph); g.stroke(); g.setLineDash([]);
    }
    g.save();
    g.beginPath(); g.rect(padL, padT, pw, ph); g.clip();
    for (const s of this.series) {
      g.strokeStyle = s.color; g.lineWidth = s.width ?? 1.6; g.setLineDash(s.dash ?? []);
      g.beginPath();
      let started = false;
      for (let i = 0; i < xs.length; i++) {
        const v = s.data[i];
        if (!Number.isFinite(v)) { started = false; continue; }
        const px = X(xs[i]), py = Y(v);
        if (!started) { g.moveTo(px, py); started = true; } else g.lineTo(px, py);
      }
      g.stroke();
      if (s.points) {
        g.fillStyle = s.color;
        for (let i = 0; i < xs.length; i++) if (Number.isFinite(s.data[i])) { g.beginPath(); g.arc(X(xs[i]), Y(s.data[i]), 2.4, 0, 6.283); g.fill(); }
      }
    }
    g.setLineDash([]);
    g.restore();
    this.legend(g, padL, H);
  }

  legend(g, padL, H) {
    let lx = padL;
    g.font = '10.5px system-ui, sans-serif'; g.textAlign = 'left';
    for (const s of this.series) {
      if (!s.name) continue;
      g.fillStyle = s.color; g.fillRect(lx, H - 12, 9, 3);
      g.fillStyle = '#b8c4d3'; g.fillText(s.name, lx + 12, H - 8);
      lx += 16 + g.measureText(s.name).width + 8;
    }
  }
}

export function fmt(v) {
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 1e4 || a < 1e-2) return v.toExponential(1);
  if (a >= 100) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  return v.toFixed(2);
}
