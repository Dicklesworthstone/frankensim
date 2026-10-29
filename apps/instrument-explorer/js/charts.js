// SVG charts drawn from the built data. Colours come from CSS classes so light and dark themes both work.
import { s, h, fmt } from './app.js';

export const SECTION_COLORS = {
  Reeds: '#4c6ef5',
  'Flutes and pipes': '#15aabf',
  Brass: '#f59f00',
  'Bowed strings': '#e64980',
  'Plucked strings': '#be4bdb',
  Keyboard: '#7950f2',
  Mallets: '#12b886',
  Drums: '#fd7e14',
  Voice: '#82c91e',
  'Several families': '#868e96',
  'Shared physics and tools': '#adb5bd',
};

const title = text => s('title', null, text);

/**
 * One measure against its ruler: the shaded band is the spread between real recordings (0 to the largest real-vs-real
 * distance), grey ticks are the individual real-vs-real distances, the dashed line the same-note (strict) maximum, and the
 * dot our render's distance from the recording.
 */
export function rulerBar(m) {
  const W = 240;
  const H = 22;
  const values = [m.gap, m.ruler_max, m.ruler_strict_max, ...m.ruler.map(r => r.value)].filter(v => Number.isFinite(v));
  const max = Math.max(...values.map(Math.abs), 1e-9) * 1.08;
  const x = v => 6 + ((W - 12) * Math.abs(v)) / max;
  const kids = [s('line.track', { x1: 6, x2: W - 6, y1: H / 2, y2: H / 2 })];
  if (Number.isFinite(m.ruler_max)) {
    kids.push(s('rect.spread', { x: 6, y: 4, width: x(m.ruler_max) - 6, height: H - 8, rx: 2 },
      title(`real spread: up to ${fmt.num(m.ruler_max, 3)} ${m.unit}`)));
  }
  for (const r of m.ruler) {
    kids.push(s('line.rtick', { x1: x(r.value), x2: x(r.value), y1: 5, y2: H - 5 }, title(`${r.label}: ${fmt.num(r.value, 3)}`)));
  }
  if (Number.isFinite(m.ruler_strict_max) && m.ruler_strict_max !== m.ruler_max) {
    kids.push(s('line.strict', { x1: x(m.ruler_strict_max), x2: x(m.ruler_strict_max), y1: 2, y2: H - 2 },
      title(`same note only: up to ${fmt.num(m.ruler_strict_max, 3)}`)));
  }
  if (Number.isFinite(m.gap)) {
    kids.push(s(`circle.gap.${m.inside ? 'in' : m.inside === false ? 'out' : 'na'}`, { cx: x(m.gap), cy: H / 2, r: 5 },
      title(`ours: ${fmt.num(m.gap, 3)} ${m.unit}`)));
  }
  return s('svg.ruler', { viewBox: `0 0 ${W} ${H}`, width: W, height: H }, kids);
}

function niceMax(v) {
  if (v <= 0) return 1;
  const p = 10 ** Math.floor(Math.log10(v));
  return [1, 2, 2.5, 5, 10].map(k => k * p).find(k => k >= v);
}

/** Stacked daily bars. rows: [{label, by: {key: n}}]; keys in stacking order. */
export function stackedBars(rows, keys, colors, { height = 220, yLabel = '' } = {}) {
  const W = 900;
  const H = height;
  const L = 40;
  const B = 34;
  const T = 10;
  const totals = rows.map(r => keys.reduce((a, k) => a + (r.by[k] || 0), 0));
  const ymax = niceMax(Math.max(...totals, 1));
  const bw = (W - L - 10) / Math.max(rows.length, 1);
  const y = v => T + (H - T - B) * (1 - v / ymax);
  const kids = [];
  // counts: whole-number ticks (a 25 top gives 0, 5 … 25, not 6.25 and 18.75)
  const p = 10 ** Math.floor(Math.log10(ymax / 4));
  const step = [1, 2, 2.5, 5, 10].map(k => k * p).find(st => Number.isInteger(st) && ymax / st <= 5) || Math.max(1, ymax / 4);
  for (let v = 0; v <= ymax + 1e-9; v += step) {
    kids.push(s('line.grid', { x1: L, x2: W - 10, y1: y(v), y2: y(v) }), s('text.axis', { x: L - 6, y: y(v) + 4, 'text-anchor': 'end' }, String(v)));
  }
  rows.forEach((r, i) => {
    let acc = 0;
    const x0 = L + i * bw + bw * 0.12;
    for (const k of keys) {
      const n = r.by[k] || 0;
      if (!n) continue;
      kids.push(s('rect', { x: x0, y: y(acc + n), width: bw * 0.76, height: y(acc) - y(acc + n), fill: colors[k] || '#999' },
        title(`${r.label}: ${n} ${k}`)));
      acc += n;
    }
    const every = Math.ceil(rows.length / 12);
    if (i % every === 0) kids.push(s('text.axis', { x: x0 + bw * 0.38, y: H - B + 16, 'text-anchor': 'middle' }, fmt.date(r.label)));
  });
  if (yLabel) kids.push(s('text.axis', { x: L + 4, y: T + 4 }, yLabel));
  return s('svg.chart', { viewBox: `0 0 ${W} ${H}`, preserveAspectRatio: 'xMidYMid meet' }, kids);
}

export function legend(keys, colors) {
  return h('div.legend', null, keys.map(k => h('span', null, h('i', { style: { background: colors[k] } }), k)));
}

/**
 * Inside fraction over grading snapshots, one line per note. series: [{name, points: [{t: Date, y, label}], muted?, color?}];
 * a muted series is a thin grey line with no end label (its points keep their titles), for when there are too many to colour.
 */
export function historyChart(series, { height = 200 } = {}) {
  // a phone draws it narrower, so its text is not scaled down to unreadable; tall enough for every end label 12 apart
  const W = (window.innerWidth || 1200) < 760 ? 420 : 900;
  const H = Math.max(height, series.filter(x => !x.muted).length * 12 + 16);
  const L = 44;
  const B = 30;
  const T = 10;
  const R = 110;
  const all = series.flatMap(x => x.points);
  if (!all.length) return null;
  let t0 = Math.min(...all.map(p => p.t));
  let t1 = Math.max(...all.map(p => p.t));
  if (t1 - t0 < 3.6e6) {
    t0 -= 3.6e6;
    t1 += 3.6e6;
  }
  const x = t => L + ((W - L - R) * (t - t0)) / (t1 - t0);
  const y = v => T + (H - T - B) * (1 - v);
  const kids = [];
  for (const v of [0, 0.25, 0.5, 0.75, 1]) {
    kids.push(s('line.grid', { x1: L, x2: W - R, y1: y(v), y2: y(v) }), s('text.axis', { x: L - 6, y: y(v) + 4, 'text-anchor': 'end' }, `${v * 100}%`));
  }
  const days = Math.ceil((t1 - t0) / 864e5);
  // a date label needs about 48 units, so a long span labels every few days
  const every = Math.max(1, Math.ceil((days * 48) / (W - L - R)));
  for (let d = 0; d <= days; d += every) {
    // local midnights: fmt.date labels in local time, as the grading times in the cards are
    const t = new Date(t0);
    t.setHours(0, 0, 0, 0);
    t.setDate(t.getDate() + d);
    if (t < t0 || t > t1) continue;
    kids.push(s('line.grid.v', { x1: x(t), x2: x(t), y1: T, y2: H - B }), s('text.axis', { x: x(t), y: H - B + 16, 'text-anchor': 'middle' }, fmt.date(t.toISOString())));
  }
  const palette = ['#4c6ef5', '#e64980', '#12b886', '#f59f00', '#7950f2', '#15aabf', '#fd7e14', '#82c91e', '#be4bdb', '#868e96'];
  const ends = [];
  let k = 0;
  // muted series first, so the coloured lines draw over them
  for (const ser of [...series.filter(x => x.muted), ...series.filter(x => !x.muted)]) {
    const pts = ser.points.slice().sort((a, b) => a.t - b.t);
    // a step: a grading's score holds until the next grading replaces it, so nothing slopes between them
    const line = pts.map((p, i) => (i ? `${x(p.t)},${y(pts[i - 1].y)} ` : '') + `${x(p.t)},${y(p.y)}`).join(' ');
    if (ser.muted) {
      if (pts.length > 1) kids.push(s('polyline.faint', { points: line }));
      for (const p of pts) kids.push(s('circle.faint', { cx: x(p.t), cy: y(p.y), r: 2.5 }, title(`${ser.name} ${p.label}`)));
      continue;
    }
    // as a style, so a CSS variable works as well as a hex colour
    const c = ser.color || palette[k++ % palette.length];
    if (pts.length > 1) kids.push(s('polyline', { points: line, fill: 'none', 'stroke-width': 2, style: `stroke: ${c}` }));
    for (const p of pts) kids.push(s('circle', { cx: x(p.t), cy: y(p.y), r: 3.5, style: `fill: ${c}` }, title(`${ser.name} ${p.label}`)));
    ends.push({ name: ser.name, c, y: y(pts[pts.length - 1].y) + 4 });
  }
  // End labels in their line's colour (as a style: the .axis rule would override a fill attribute), pushed at least
  // 12 px apart so lines that finish at the same value stay readable.
  ends.sort((a, b) => a.y - b.y);
  for (let i = 1; i < ends.length; i++) ends[i].y = Math.max(ends[i].y, ends[i - 1].y + 12);
  const over = ends.length ? ends[ends.length - 1].y - (H - 4) : 0;
  for (const e of ends) kids.push(s('text.axis', { x: W - R + 8, y: over > 0 ? e.y - over : e.y, style: `fill: ${e.c}` }, e.name));
  return s('svg.chart', { viewBox: `0 0 ${W} ${H}` }, kids);
}

/** Each run's identified body modes against the pre-registered bands. */
export function bodyModes(e) {
  const runs = e.runs;
  const W = 900;
  const L = 150;
  const R = 20;
  const T = 30;
  const rowH = 18;
  const H = T + rowH * runs.length + 30;
  const found = runs.flatMap(r => Object.values(r.modes).map(m => m.f_hz)).filter(Number.isFinite);
  const f0 = Math.min(200, 50 * Math.floor(Math.min(...found, 200) / 50 - 0.2));
  const f1 = Math.max(720, 50 * Math.ceil(Math.max(...found, 720) / 50 + 0.2));
  const x = f => L + ((W - L - R) * (f - f0)) / (f1 - f0);
  const kids = [];
  for (const b of e.bands) {
    kids.push(s('rect.band', { x: x(b.lo), y: T - 6, width: x(b.hi) - x(b.lo), height: rowH * runs.length + 6 }, title(`${b.id}: ${b.lo}-${b.hi} Hz`)));
    kids.push(s('text.bandlabel', { x: (x(b.lo) + x(b.hi)) / 2, y: T - 12, 'text-anchor': 'middle' }, b.id.toUpperCase()));
  }
  for (let f = f0 + 50; f < f1; f += 50) {
    kids.push(s('line.grid.v', { x1: x(f), x2: x(f), y1: T - 4, y2: T + rowH * runs.length }), s('text.axis', { x: x(f), y: H - 10, 'text-anchor': 'middle' }, `${f}`));
  }
  kids.push(s('text.axis', { x: W - R, y: H - 10, 'text-anchor': 'end' }, 'Hz'));
  const shapes = { a0: 'circle', cbr: 'rect', 'b1-': 'diamond', 'b1+': 'triangle' };
  runs.forEach((run, i) => {
    const cy = T + rowH * i + rowH / 2;
    const superseded = run.run.startsWith('cracked');
    kids.push(s(`text.axis${superseded ? '.dim' : ''}`, { x: L - 8, y: cy + 4, 'text-anchor': 'end' }, run.run));
    kids.push(s('line.grid', { x1: L, x2: W - R, y1: cy, y2: cy }));
    for (const [id, m] of Object.entries(run.modes)) {
      if (!Number.isFinite(m.f_hz)) continue;
      const band = e.bands.find(b => b.id === id);
      const hit = band && m.f_hz >= band.lo && m.f_hz <= band.hi;
      const cls = `mode.${hit ? 'hit' : 'miss'}${superseded ? '.dim' : ''}`;
      const cx = x(Math.min(Math.max(m.f_hz, f0), f1));
      const tip = title(`${run.run} ${id.toUpperCase()}: ${m.f_hz} Hz, ${m.outcome}. ${m.note}`);
      const shape = shapes[id] || 'circle';
      if (shape === 'circle') kids.push(s(`circle.${cls}`, { cx, cy, r: 5 }, tip));
      else if (shape === 'rect') kids.push(s(`rect.${cls}`, { x: cx - 4.5, y: cy - 4.5, width: 9, height: 9 }, tip));
      else if (shape === 'diamond') kids.push(s(`path.${cls}`, { d: `M${cx},${cy - 6}L${cx + 6},${cy}L${cx},${cy + 6}L${cx - 6},${cy}Z` }, tip));
      else kids.push(s(`path.${cls}`, { d: `M${cx},${cy - 6}L${cx + 6},${cy + 5}L${cx - 6},${cy + 5}Z` }, tip));
    }
  });
  return s('svg.chart', { viewBox: `0 0 ${W} ${H}` }, kids);
}

const REGIME = { playable: 'Helmholtz (playable)', raucous: 'raucous', 'surface-or-dead': 'surface sound or dead', null: 'not classified' };

/** Bow-force sweep: one strip per string (or bow position), bow force on a log axis, colour = regime. */
export function bowedSweep(e) {
  const rows = e.rows.filter(r => Number.isFinite(r.force_n));
  const keyOf = r => `${fmt.num(r.f1_hz, 0)} Hz, β ${fmt.num(r.beta, 3)}`;
  const keys = [...new Set(rows.map(keyOf))].sort((a, b) => parseFloat(a) - parseFloat(b));
  const W = 900;
  const L = 150;
  const R = 20;
  const T = 10;
  const rowH = 22;
  const H = T + rowH * keys.length + 34;
  const fmin = Math.min(...rows.map(r => r.force_n));
  const fmax = Math.max(...rows.map(r => r.force_n));
  const x = f => L + ((W - L - R) * Math.log(f / fmin)) / Math.log(fmax / fmin || 10);
  const kids = [];
  for (let p = Math.floor(Math.log10(fmin)); p <= Math.ceil(Math.log10(fmax)); p++) {
    for (const k of [1, 2, 5]) {
      const f = k * 10 ** p;
      if (f < fmin * 0.999 || f > fmax * 1.001) continue;
      kids.push(s('line.grid.v', { x1: x(f), x2: x(f), y1: T, y2: T + rowH * keys.length }), s('text.axis', { x: x(f), y: H - 16, 'text-anchor': 'middle' }, `${f}`));
    }
  }
  kids.push(s('text.axis', { x: W - R, y: H - 2, 'text-anchor': 'end' }, 'bow force [N]'));
  keys.forEach((k, i) => {
    const cy = T + rowH * i + rowH / 2;
    kids.push(s('text.axis', { x: L - 8, y: cy + 4, 'text-anchor': 'end' }, k));
    for (const r of rows.filter(q => keyOf(q) === k)) {
      kids.push(s(`circle.regime.r-${r.cls || 'none'}`, { cx: x(r.force_n), cy, r: 5 },
        title(`${fmt.num(r.force_n, 3)} N, ${REGIME[r.cls] || r.cls}; slip fraction ${fmt.num(r.slip, 3)}; spectral peak ${fmt.num(r.peak_hz, 1)} Hz`)));
    }
  });
  return h('div', null, s('svg.chart', { viewBox: `0 0 ${W} ${H}` }, kids),
    h('div.legend', null, ['playable', 'raucous', 'surface-or-dead', 'none'].map(c => h('span', null, h(`i.regime.r-${c}`), REGIME[c === 'none' ? null : c]))));
}

/** Drum decay per octave band: the render (filled; an arrow when it is a lower bound) against the recordings (grey). */
export function bandDecays(e) {
  const W = 900;
  const H = 260;
  const L = 56;
  const B = 40;
  const T = 12;
  const R = 12;
  const bands = e.bands_hz;
  const pieces = e.pieces;
  const allMs = pieces.flatMap(p => [...p.render_decays, ...p.references.flatMap(r => r.decays)].map(d => d.decay_ms)).filter(v => v > 0);
  const lo = 10;
  const hi = niceMax(Math.max(...allMs) * 1.2);
  const y = v => T + ((H - T - B) * (1 - Math.log(Math.max(v, lo) / lo) / Math.log(hi / lo)));
  const groupW = (W - L - R) / bands.length;
  const pieceW = groupW / (pieces.length + 1);
  const kids = [];
  for (const v of [10, 30, 100, 300, 1000, 3000, 10000]) {
    if (v > hi) continue;
    kids.push(s('line.grid', { x1: L, x2: W - R, y1: y(v), y2: y(v) }), s('text.axis', { x: L - 6, y: y(v) + 4, 'text-anchor': 'end' }, v >= 1000 ? `${v / 1000} s` : `${v} ms`));
  }
  const palette = ['#fd7e14', '#e64980', '#4c6ef5', '#12b886', '#7950f2'];
  bands.forEach((fc, bi) => {
    const gx = L + groupW * bi;
    kids.push(s('text.axis', { x: gx + groupW / 2, y: H - B + 16, 'text-anchor': 'middle' }, fc >= 1000 ? `${fc / 1000}k` : String(fc)));
    pieces.forEach((p, pi) => {
      const cx = gx + pieceW * (pi + 1);
      for (const r of p.references) {
        const d = r.decays.find(q => q.fc_hz === fc);
        if (d && d.decay_ms > 0 && d.mark !== 'nm') kids.push(s('circle.refdot', { cx, cy: y(d.decay_ms), r: 2.5 }, title(`${p.piece} recording ${r.id}: ${d.decay_ms} ms`)));
      }
      const d = p.render_decays.find(q => q.fc_hz === fc);
      if (d && d.decay_ms > 0) {
        const c = palette[pi % palette.length];
        const lower = d.mark === '>';
        const cy = y(d.decay_ms);
        kids.push(lower
          ? s('path', { d: `M${cx},${cy - 8}L${cx + 5},${cy}L${cx - 5},${cy}Z`, fill: c }, title(`${p.piece} render: more than ${d.decay_ms} ms (still ringing at the end of the clip)`))
          : s('circle', { cx, cy, r: 4.5, fill: c }, title(`${p.piece} render: ${d.decay_ms} ms`)));
      }
    });
  });
  kids.push(s('text.axis', { x: W - R, y: H - 4, 'text-anchor': 'end' }, 'octave band centre [Hz]'));
  return h('div', null, s('svg.chart', { viewBox: `0 0 ${W} ${H}` }, kids),
    h('div.legend', null, pieces.map((p, i) => h('span', null, h('i', { style: { background: palette[i % palette.length] } }), `${p.piece} (render)`)),
      h('span', null, h('i.refdot'), 'recordings of the same kind of drum')));
}
