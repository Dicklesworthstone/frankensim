// Rendering and measurement only. No constitutive or simulation equations here.
export function format(value, digits = 3) {
  if (value == null || !Number.isFinite(value)) return '—';
  if (value === 0) return '0';
  if (Math.abs(value) < 0.001 || Math.abs(value) >= 1e5) return value.toExponential(2);
  return Number(value.toPrecision(digits)).toString();
}
export function nearest(times, x) {
  let lo = 0, hi = times.length - 1;
  while (lo < hi) {const mid = (lo + hi) >>> 1; if (times[mid] < x) lo = mid + 1; else hi = mid;}
  return lo > 0 && Math.abs(times[lo - 1] - x) < Math.abs(times[lo] - x) ? lo - 1 : lo;
}
export function spectrum(values, rate) {
  if (values.length < 4) return null;
  let n = 1; while (n < values.length) n *= 2;
  const re = new Float64Array(n), im = new Float64Array(n); let weight = 0;
  for (let i = 0; i < values.length; i++) {
    const w = 0.5 * (1 - Math.cos(2 * Math.PI * i / (values.length - 1)));
    re[i] = values[i] * w; weight += w;
  }
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
    if (i < j) {[re[i], re[j]] = [re[j], re[i]];}
  }
  for (let size = 2; size <= n; size *= 2) {
    const a = -2 * Math.PI / size, wr = Math.cos(a), wi = Math.sin(a);
    for (let begin = 0; begin < n; begin += size) {
      let ur = 1, ui = 0;
      for (let k = 0; k < size / 2; k++) {
        const i = begin + k, j = i + size / 2;
        const vr = re[j] * ur - im[j] * ui, vi = re[j] * ui + im[j] * ur;
        re[j] = re[i] - vr; im[j] = im[i] - vi; re[i] += vr; im[i] += vi;
        [ur, ui] = [ur * wr - ui * wi, ur * wi + ui * wr];
      }
    }
  }
  const x = [], y = [];
  for (let k = 0; k <= n / 2; k++) {
    x.push(k * rate / n); y.push(Math.hypot(re[k], im[k]) / weight * (k === 0 || k === n / 2 ? 1 : 2));
  }
  return {x, y};
}
function surface(canvas) {
  const box = canvas.getBoundingClientRect(), ratio = Math.min(window.devicePixelRatio || 1, 2);
  const width = Math.max(1, box.width), height = Math.max(1, box.height);
  const w = Math.round(width * ratio), h = Math.round(height * ratio);
  if (canvas.width !== w || canvas.height !== h) {canvas.width = w; canvas.height = h;}
  const ctx = canvas.getContext('2d'); ctx.setTransform(ratio, 0, 0, ratio, 0, 0); ctx.clearRect(0, 0, width, height);
  const css = getComputedStyle(document.documentElement);
  const color = key => css.getPropertyValue('--' + key).trim();
  return {ctx, width, height, color};
}
export function palette() {
  const css = getComputedStyle(document.documentElement);
  return [css.getPropertyValue('--cold').trim(), css.getPropertyValue('--hot').trim(),
    '#b59be9', '#9ac78d', '#d987b2', '#cab674'];
}
// Keep each occupied pixel's actual min/max samples rather than stride-dropping peaks.
export function envelopeIndices(x, y, xmin, xmax, width) {
  const points = []; let bucket = -1, imin = 0, imax = 0;
  const flush = () => {if (bucket >= 0) {points.push(Math.min(imin, imax)); if (imin !== imax) points.push(Math.max(imin, imax));}};
  for (let i = 0; i < x.length; i++) {
    if (x[i] < xmin || x[i] > xmax) continue;
    const next = Math.floor((x[i] - xmin) / (xmax - xmin) * width);
    if (next !== bucket) {flush(); bucket = next; imin = i; imax = i;}
    else {if (y[i] < y[imin]) imin = i; if (y[i] > y[imax]) imax = i;}
  }
  flush(); return points;
}
export function drawPlot(canvas, lines, {range = [0, 1], cursor = null, xunit = 's', yunit = '', positive = false} = {}) {
  const {ctx: c, width: w, height: h, color} = surface(canvas);
  const pad = {l: 52, r: 20, t: 16, b: 25}, pw = Math.max(1, w - pad.l - pad.r), ph = Math.max(1, h - pad.t - pad.b);
  let min = Infinity, max = -Infinity;
  for (const line of lines) for (let i = 0; i < line.x.length; i++) {
    if (line.x[i] >= range[0] && line.x[i] <= range[1]) {min = Math.min(min, line.y[i]); max = Math.max(max, line.y[i]);}
  }
  if (!Number.isFinite(min)) {min = -1; max = 1;}
  if (positive) min = 0;
  if (max === min) {const span = Math.abs(max) * .1 || 1; min -= positive ? 0 : span; max += span;}
  const margin = (max - min) * .1; max += margin; if (!positive) min -= margin;
  const X = x => pad.l + (x - range[0]) / (range[1] - range[0]) * pw;
  const Y = y => pad.t + ph - (y - min) / (max - min) * ph;
  c.font = '9px ui-monospace, monospace'; c.lineWidth = 1; c.strokeStyle = color('line-soft');
  c.fillStyle = color('dim');
  for (let i = 0; i <= 4; i++) {
    const x = pad.l + pw * i / 4, y = pad.t + ph * i / 4;
    c.beginPath(); c.moveTo(x, pad.t); c.lineTo(x, pad.t + ph); c.stroke();
    c.beginPath(); c.moveTo(pad.l, y); c.lineTo(pad.l + pw, y); c.stroke();
    c.textAlign = 'right'; c.fillText(format(max - (max - min) * i / 4), pad.l - 8, y + 3);
    c.textAlign = 'center'; const v = range[0] + (range[1] - range[0]) * i / 4;
    c.fillText(format(xunit === 'ms' ? v * 1000 : v), x, h - 8);
  }
  c.textAlign = 'left'; c.fillText(yunit, 4, 9); c.textAlign = 'right'; c.fillText(xunit, w - 3, 9);
  c.save(); c.beginPath(); c.rect(pad.l, pad.t, pw, ph); c.clip();
  for (const line of lines) {
    c.strokeStyle = line.color; c.lineWidth = 1.25; c.setLineDash(line.dashed ? [4, 4] : []); c.beginPath();
    const points = envelopeIndices(line.x, line.y, range[0], range[1], Math.ceil(pw));
    points.forEach((i, n) => {const x = X(line.x[i]), y = Y(line.y[i]); if (n === 0) c.moveTo(x, y); else c.lineTo(x, y);});
    c.stroke();
  }
  c.setLineDash([]);
  if (cursor != null && cursor >= range[0] && cursor <= range[1]) {
    c.strokeStyle = color('muted'); c.setLineDash([3, 4]); c.beginPath(); c.moveTo(X(cursor), pad.t); c.lineTo(X(cursor), pad.t + ph); c.stroke();
  }
  c.restore();
  return {start: pad.l, width: pw, range};
}
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const mix = (a, b, t) => a.map((x, i) => x * (1 - t) + b[i] * t);
const rgb = (v, light = 1) => `rgb(${v.map(x => Math.round(clamp(x * light, 0, 255))).join(',')})`;
export class SectionView {
  constructor(canvas, onSelect) {
    this.canvas = canvas; this.onSelect = onSelect; this.yaw = -0.48; this.pitch = 0.28; this.zoom = 1; this.options = {}; this.faces = [];
    let drag = null;
    canvas.addEventListener('pointerdown', e => {drag = {id: e.pointerId, x: e.clientX, y: e.clientY, distance: 0}; canvas.setPointerCapture(e.pointerId);});
    canvas.addEventListener('pointermove', e => {
      if (!drag || drag.id !== e.pointerId) return;
      const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.distance += Math.abs(dx) + Math.abs(dy);
      this.yaw += dx * .006; this.pitch = clamp(this.pitch + dy * .006, -1.1, 1.1); drag.x = e.clientX; drag.y = e.clientY; this.invalidate();
    });
    canvas.addEventListener('pointerup', e => {
      if (!drag || drag.id !== e.pointerId) return;
      if (drag.distance < 5) {
        const r = canvas.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
        for (const face of [...this.faces].reverse()) {
          let inside = false;
          for (let i = 0, j = face.p.length - 1; i < face.p.length; j = i++) {
            const a = face.p[i], b = face.p[j];
            if ((a[1] > y) !== (b[1] > y) && x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0]) inside = !inside;
          }
          if (inside) {this.onSelect(face.layer); break;}
        }
      }
      drag = null;
    });
    canvas.addEventListener('pointercancel', () => {drag = null;});
    canvas.addEventListener('wheel', e => {e.preventDefault(); this.zoomBy(Math.exp(-e.deltaY * .001));}, {passive: false});
    canvas.addEventListener('keydown', e => {
      if (!['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown', '+', '-', '=', '0'].includes(e.key)) return;
      e.preventDefault();
      if (e.key === '0') this.reset();
      else if (['+', '=', '-'].includes(e.key)) this.zoomBy(e.key === '-' ? .9 : 1.1);
      else {this.yaw += e.key === 'ArrowLeft' ? -.08 : e.key === 'ArrowRight' ? .08 : 0;
        this.pitch = clamp(this.pitch + (e.key === 'ArrowUp' ? -.08 : e.key === 'ArrowDown' ? .08 : 0), -1.1, 1.1); this.invalidate();}
    });
    this.observer = new ResizeObserver(() => this.invalidate()); this.observer.observe(canvas); this.invalidate();
  }
  reset() {this.yaw = -.48; this.pitch = .28; this.zoom = 1; this.invalidate();}
  zoomBy(f) {this.zoom = clamp(this.zoom * f, .5, 1.8); this.invalidate();}
  update(options) {this.options = options; this.invalidate();}
  invalidate() {
    if (this.pending || document.hidden) return;
    this.pending = requestAnimationFrame(() => {this.pending = null; this.draw();});
  }
  draw() {
    const {ctx: c, width: w, height: h, color} = surface(this.canvas);
    const {snapshot, state = 'hot', field = 'temperature', selected = 0, origin = 'Preview', baseline} = this.options;
    const s = snapshot?.section; const temperatures = s?.[state + '_temperature_k'];
    const stress = state === 'hot' ? s?.hot_stress_pa : null;
    const values = field === 'stress' ? stress : temperatures;
    const upper = Math.max(300.001, ...(s?.hot_temperature_k || []), ...(s?.cold_temperature_k || []),
      ...(baseline?.section?.hot_temperature_k || []));
    const lower = Math.min(300, ...(s?.cold_temperature_k || []), ...(s?.hot_temperature_k || []));
    const stressMax = Math.max(1, ...((s?.hot_stress_pa || []).map(Math.abs)));
    const colors = values ? values.map(v => field === 'stress'
      ? mix([83, 171, 197], [246, 149, 92], clamp((v / stressMax + 1) / 2, 0, 1))
      : mix([68, 163, 191], [246, 158, 99], clamp((v - lower) / (upper - lower), 0, 1)))
      : (snapshot ? [[80, 110, 130], [105, 119, 131]] : [[72, 167, 193], [153, 116, 88]]);
    const scale = Math.min(w / 7.8, h / 3.9) * this.zoom, cx = w * .49, cy = h * .47;
    const rotate = ([x, y, z]) => {
      const a = x * Math.cos(this.yaw) + z * Math.sin(this.yaw), b = -x * Math.sin(this.yaw) + z * Math.cos(this.yaw);
      return [a, y * Math.cos(this.pitch) - b * Math.sin(this.pitch), y * Math.sin(this.pitch) + b * Math.cos(this.pitch)];
    };
    const project = p => {const [x, y, z] = rotate(p), k = 7 / (7 - z); return [cx + x * scale * k, cy - y * scale * k, z];};
    // A diagrammatic ground grid; never an FE mesh or extra spatial observations.
    c.strokeStyle = color('line-soft'); c.lineWidth = .6;
    for (let i = -8; i <= 8; i++) {
      const a = project([i * .6, -1.1, -4]), b = project([i * .6, -1.1, 4]);
      c.beginPath(); c.moveTo(a[0], a[1]); c.lineTo(b[0], b[1]); c.stroke();
      const u = project([-5, -1.1, i * .5]), v = project([5, -1.1, i * .5]);
      c.beginPath(); c.moveTo(u[0], u[1]); c.lineTo(v[0], v[1]); c.stroke();
    }
    const faces = [], radius = s ? s.radii_m[0] / s.radii_m[1] : .5;
    const P = (x, r, a) => [x, r * Math.cos(a), r * Math.sin(a)];
    const add = (pts, layer, light = 1, cap = false) => {
      const p = pts.map(project); faces.push({p, layer, cap, light, z: p.reduce((a, b) => a + b[2], 0) / p.length});
    };
    const length = 4.05, segments = 64;
    for (let i = 0; i < segments; i++) {
      const a = i / segments * 2 * Math.PI, b = (i + 1) / segments * 2 * Math.PI;
      const light = .6 + .34 * Math.max(0, Math.cos(a - .6));
      for (let j = 0; j < 10; j++) {
        const x = -length / 2 + length * j / 10, y = x + length / 10;
        add([P(x, radius, a), P(y, radius, a), P(y, radius, b), P(x, radius, b)], 0, light);
      }
      for (const x of [-length / 2, length / 2]) add([[x, 0, 0], P(x, radius, a), P(x, radius, b)], 0, .9, true);
      // Remove a front/upper wedge to expose the core and material interface.
      if (a >= .1 && a < 1.47) continue;
      for (let j = 0; j < 10; j++) {
        const x = -length / 2 + length * j / 10, y = x + length / 10;
        add([P(x, 1, a), P(y, 1, a), P(y, 1, b), P(x, 1, b)], 1, light);
      }
      for (const x of [-length / 2, length / 2]) add([P(x, radius, a), P(x, 1, a), P(x, 1, b), P(x, radius, b)], 1, 1.02, true);
    }
    for (const a of [Math.PI / 32, 15 * Math.PI / 32]) add([P(-length / 2, radius, a), P(length / 2, radius, a),
      P(length / 2, 1, a), P(-length / 2, 1, a)], 1, 1.1);
    faces.sort((a, b) => a.z - b.z); this.faces = faces;
    for (const f of faces) {
      c.beginPath(); f.p.forEach((p, i) => i ? c.lineTo(p[0], p[1]) : c.moveTo(p[0], p[1])); c.closePath();
      c.fillStyle = rgb(colors[f.layer], f.light); c.fill();
      c.strokeStyle = f.cap && f.layer === selected ? '#d9eff333' : '#06131b25'; c.lineWidth = .45; c.stroke();
    }
    // Dimension labels report the actual specimen, despite diagrammatic aspect ratio.
    const da = project([-length / 2, -1.35, 0]), db = project([length / 2, -1.35, 0]);
    c.strokeStyle = color('dim'); c.lineWidth = .65; c.beginPath(); c.moveTo(da[0], da[1]); c.lineTo(db[0], db[1]); c.stroke();
    for (const p of [da, db]) {c.beginPath(); c.moveTo(p[0], p[1] - 4); c.lineTo(p[0], p[1] + 4); c.stroke();}
    c.font = '9px ui-monospace, monospace'; c.textAlign = 'center'; c.fillStyle = color('muted');
    c.fillText(s ? `${format(s.length_m)} m · fixed reference geometry` : 'REFERENCE GEOMETRY · NOT A RESULT', (da[0] + db[0]) / 2, (da[1] + db[1]) / 2 + 16);
    c.font = '9px ui-monospace, monospace'; c.textAlign = 'left'; c.fillStyle = color('dim');
    c.fillText(snapshot ? `${origin.toUpperCase()} · ${state.toUpperCase()} SNAPSHOT` : 'GEOMETRY PREVIEW · NO SOLVER DATA', 20, 16);
    if (values) {
      const gx = w - 140, gy = 18, gradient = c.createLinearGradient(gx, gy, gx + 100, gy);
      gradient.addColorStop(0, rgb([68, 163, 191])); gradient.addColorStop(1, rgb([246, 158, 99]));
      c.fillStyle = gradient; c.fillRect(gx, gy, 100, 3); c.fillStyle = color('muted'); c.font = '8px ui-monospace, monospace';
      c.textAlign = 'left'; c.fillText(field === 'stress' ? `${format(-stressMax)} Pa` : `${format(lower, 5)} K`, gx, gy + 14);
      c.textAlign = 'right'; c.fillText(field === 'stress' ? `${format(stressMax)} Pa` : `${format(upper, 5)} K`, gx + 100, gy + 14);
    }
    c.textAlign = 'left';
  }
}
