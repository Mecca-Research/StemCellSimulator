// Declarative control-panel builder used by every scene.
import { Chart, fmt } from './chart.js';

function el(tag, attrs = {}, children = []) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') e.className = v;
    else if (k === 'text') e.textContent = v;
    else if (k === 'html') e.innerHTML = v;
    else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
    else e.setAttribute(k, v);
  }
  for (const c of children) e.append(c);
  return e;
}

export class Panel {
  constructor(root) {
    this.root = root;
    this.charts = [];
    this.current = root;
  }

  clear() {
    this.root.replaceChildren();
    this.charts = [];
    this.current = this.root;
  }

  section(title) {
    const s = el('div', { class: 'ctl-section' }, [el('h3', { text: title })]);
    this.root.append(s);
    this.current = s;
    return this;
  }

  slider({ label, min, max, step = 0.01, value, onChange, format = fmt, log = false }) {
    const toSlider = (v) => (log ? Math.log10(v) : v);
    const fromSlider = (s) => (log ? 10 ** s : s);
    const val = el('span', { class: 'val', text: format(value) });
    const input = el('input', {
      type: 'range', min: toSlider(min), max: toSlider(max), step: log ? (Math.log10(max) - Math.log10(min)) / 200 : step,
      value: toSlider(value), 'aria-label': label,
    });
    input.addEventListener('input', () => {
      const v = fromSlider(parseFloat(input.value));
      val.textContent = format(v);
      onChange?.(v);
    });
    this.current.append(el('div', { class: 'ctl' }, [el('div', { class: 'ctl-row' }, [el('span', { text: label }), val]), input]));
    return {
      input,
      set: (v) => { input.value = toSlider(v); val.textContent = format(v); },
    };
  }

  toggle({ label, value = false, onChange }) {
    const input = el('input', { type: 'checkbox' });
    input.checked = value;
    input.addEventListener('change', () => onChange?.(input.checked));
    this.current.append(el('label', { class: 'ctl ctl-toggle' }, [input, el('span', { text: label })]));
    return { input, set: (v) => { input.checked = v; } };
  }

  select({ label, options, value, onChange }) {
    const sel = el('select', { 'aria-label': label });
    for (const o of options) {
      const opt = el('option', { value: o.value ?? o, text: o.label ?? o });
      sel.append(opt);
    }
    if (value !== undefined) sel.value = value;
    sel.addEventListener('change', () => onChange?.(sel.value));
    this.current.append(el('div', { class: 'ctl' }, [el('div', { class: 'ctl-row' }, [el('span', { text: label })]), sel]));
    return { input: sel, set: (v) => { sel.value = v; } };
  }

  buttons(list) {
    const row = el('div', { class: 'btn-row' });
    const out = [];
    for (const b of list) {
      const btn = el('button', { class: `btn${b.primary ? ' primary' : ''}`, text: b.label, onclick: b.onClick });
      row.append(btn);
      out.push(btn);
    }
    this.current.append(row);
    return out;
  }

  /** Key/value readouts; returns {set(key, value)} */
  readouts(keys) {
    const grid = el('div', { class: 'readouts' });
    const cells = {};
    for (const k of keys) {
      const v = el('span', { class: 'v', text: '–' });
      grid.append(el('span', { class: 'k', text: k }), v);
      cells[k] = v;
    }
    this.current.append(grid);
    return {
      set: (k, value) => { if (cells[k]) cells[k].textContent = typeof value === 'number' ? fmt(value) : value; },
    };
  }

  chart({ title, height = 150, ...opts }) {
    const canvas = el('canvas', { class: 'chart' });
    canvas.style.height = `${height}px`;
    const wrap = el('div', { class: 'chart-wrap' }, [title ? el('div', { class: 'chart-title', text: title }) : '', canvas]);
    this.current.append(wrap);
    const ch = new Chart(canvas, opts);
    this.charts.push(ch);
    return ch;
  }

  equation(text) {
    this.current.append(el('div', { class: 'eq', text }));
  }

  note(html) {
    this.current.append(el('p', { class: 'note', html }));
  }

  html(html) {
    const d = el('div', { html });
    this.current.append(d);
    return d;
  }

  table(headers, rows) {
    const t = el('table', { class: 'data' });
    t.append(el('tr', {}, headers.map((h) => el('th', { text: h }))));
    for (const r of rows) {
      t.append(el('tr', {}, r.map((c) => el('td', { class: typeof c === 'number' ? 'num' : '', text: typeof c === 'number' ? fmt(c) : c }))));
    }
    this.current.append(t);
    return t;
  }

  drawCharts() {
    for (const c of this.charts) c.draw();
  }
}

export function setLegend(root, items) {
  root.replaceChildren(...items.map((it) => el('div', { class: 'legend-row' }, [
    el('span', { class: 'legend-sw', style: `background:${it.color}` }),
    el('span', { text: it.label }),
  ])));
}

export { el };
