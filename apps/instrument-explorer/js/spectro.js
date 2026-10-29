// Time-frequency maps. The grader's images (bare webp: x = time, 4 px per 10 ms frame; y = log frequency from f_lo at the
// bottom to f_hi at the top) get frequency and time ticks and a playhead overlaid. Clips the grader never saw (references,
// work-in-progress renders) get the same kind of map computed here: Hann STFT, log-frequency rows, 80 dB under the clip's
// own peak, magma colours.
import { h } from './app.js';

const TICKS_HZ = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];

function hzLabel(f) {
  return f >= 1000 ? `${f / 1000}k` : String(f);
}

/** How far a measured note is from its target: in tune (10 cents), near (50 cents), out, or no target. */
/** A note whose waveform has no stable period (autocorrelation peak under 0.5) has no pitch to compare: it counts
 *  as a miss, and its measured "f0" is not drawn or quoted. */
export function unpitched(n) {
  return Boolean(n && n.target_hz && n.periodicity !== undefined && n.periodicity !== null && n.periodicity < 0.5);
}

export function pitchClass(n) {
  if (n && n.refused) return 'refused';  // the solver refused this note: its slot is silent, its message kept
  if (unpitched(n)) return 'out';
  if (!n || n.cents === null || n.cents === undefined) return 'none';
  const c = Math.abs(n.cents);
  return c <= 10 ? 'in' : c <= 50 ? 'near' : 'out';
}

const SVG_NS = 'http://www.w3.org/2000/svg';

/** The notes a phrase asks for, drawn on its map: a dashed line at each target pitch (coloured by how far the render
 *  lands from it) and, when the render is more than 50 cents away, a solid line where it actually sounds. */
function noteOverlay(notes, { fLo, fHi, duration }) {
  const svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('class', 'overlay');
  svg.setAttribute('viewBox', '0 0 1000 1000');
  svg.setAttribute('preserveAspectRatio', 'none');
  svg.setAttribute('aria-hidden', 'true');
  const y = f => 1000 * (1 - Math.log(f / fLo) / Math.log(fHi / fLo));
  const line = (cls, x0, x1, yy) => {
    const l = document.createElementNS(SVG_NS, 'line');
    l.setAttribute('class', cls);
    l.setAttribute('x1', x0);
    l.setAttribute('x2', x1);
    l.setAttribute('y1', yy);
    l.setAttribute('y2', yy);
    svg.append(l);
  };
  for (const n of notes) {
    const x0 = (1000 * n.at) / duration;
    const x1 = (1000 * (n.at + n.dur)) / duration;
    line('sep', x0, x0, 0);
    const sep = svg.lastChild;
    sep.setAttribute('y2', 1000);
    const cls = pitchClass(n);
    if (n.target_hz && n.target_hz > fLo && n.target_hz < fHi) line(`tgt c-${cls}`, x0 + 2, x1 - 2, y(n.target_hz));
    if (n.f0_hz && n.f0_hz > fLo && n.f0_hz < fHi && !unpitched(n) && (cls === 'out' || !n.target_hz)) line(`got c-${n.target_hz ? 'out' : 'none'}`, x0 + 2, x1 - 2, y(n.f0_hz));
  }
  return svg;
}

/** Wrap a map element (img or canvas) with ticks and a playhead that follows `player`; a click seeks.
 *  `notes` (optional) draws the phrase's target pitches; `side` ('sim' | 'real') marks the label. */
export function mapFrame(mapEl, { fLo, fHi, duration, player, label, notes, side, tall, seekAndPlay = true }) {
  const labelled = Boolean(notes && notes.length && notes.length <= 40);
  const ticks = h('div.ticks');
  for (const f of TICKS_HZ) {
    if (f <= fLo * 1.05 || f >= fHi * 0.95) continue;
    const y = 100 * (1 - Math.log(f / fLo) / Math.log(fHi / fLo));
    if ((labelled && y < 9) || y > 92) continue; // note names run along the top, time ticks along the bottom
    ticks.append(h('span.ftick', { style: { top: `${y}%` } }, hzLabel(f)));
  }
  const step = duration > 12 ? 2 : duration > 6 ? 1 : 0.5;
  for (let t = step, i = 1; t < duration - 0.05; t += step, i++) {
    // every other time tick drops out on a narrow map
    ticks.append(h(`span.ttick${i % 2 ? '.alt' : ''}`, { style: { left: `${(100 * t) / duration}%` } }, `${t.toFixed(t % 1 ? 1 : 0)} s`));
  }
  const head = h('div.playhead');
  const overlay = notes && notes.length ? noteOverlay(notes, { fLo, fHi, duration }) : null;
  const labels = labelled
    ? h('div.notelabels', null, notes.map(n => h('span', { style: { left: `${(100 * n.at) / duration}%` } }, n.note.replace(/^lip /, ''))))
    : null;
  // a note the solver refused leaves its slot silent: say so on the simulation's map (never on the recording's)
  const refused = side !== 'real' && notes ? notes.filter(n => n.refused).map(n => h('div.refusedslot',
    { style: { left: `${(100 * n.at) / duration}%`, width: `${(100 * n.dur) / duration}%` } }, h('span', null, 'refused'))) : [];
  const tag = label ? h(`span.maplabel${side ? `.${side}` : ''}`, null, label) : null;
  const frame = h(`div.mapframe${tall ? '.tall' : ''}`, null, mapEl, overlay, refused, labels, ticks, head, tag);
  if (typeof ResizeObserver === 'function') {
    // Note names need ~30 px per note; on a narrow map they would run into each other (the note strip still lists them).
    const shortest = labels ? Math.min(...notes.map(n => n.dur)) / duration : 1;
    new ResizeObserver(([e]) => {
      frame.classList.toggle('crowded', Boolean(labels) && e.contentRect.width * shortest < 30);
      frame.classList.toggle('narrow', e.contentRect.width < 560);
      if (tag && labels) {
        // the map's name sits top right, in the row of note names: drop it a row when a name would run under it
        frame.classList.remove('tagdown');
        const t = tag.getBoundingClientRect();
        frame.classList.toggle('tagdown', !frame.classList.contains('crowded') && [...labels.children].some(s => {
          const r = s.getBoundingClientRect();
          return r.right > t.left - 4 && r.left < t.right;
        }));
      }
    }).observe(frame);
  }
  if (player) {
    const move = p => {
      head.style.left = `${(100 * Math.min(p.position, duration)) / duration}%`;
      head.style.opacity = p.playing || p.position > 0 ? '1' : '0';
    };
    player.on(move, frame);
    frame.addEventListener('click', e => {
      const r = frame.getBoundingClientRect();
      const t = ((e.clientX - r.left) / r.width) * duration;
      player.load().then(() => {
        player.seek(t);
        if (seekAndPlay) player.play();
      }, () => {});
    });
    frame.tabIndex = 0;
    frame.setAttribute('role', 'slider');
    frame.setAttribute('aria-label', `${label || 'Spectrogram'}: arrow keys move the playhead, Enter plays`);
    frame.setAttribute('aria-valuemin', '0');
    frame.setAttribute('aria-valuemax', String(duration));
    frame.addEventListener('keydown', e => {
      if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') {
        e.preventDefault();
        player.seek(player.position + (e.key === 'ArrowRight' ? 1 : -1) * (e.shiftKey ? 5 : 0.5));
      } else if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        player.toggle();
      }
    });
    frame.title = 'Click to play from here';
  }
  return frame;
}

// Magma, 9 stops (matplotlib), interpolated linearly.
const MAGMA = [
  [0, 0, 4], [28, 16, 68], [79, 18, 123], [129, 37, 129], [181, 54, 122],
  [229, 80, 100], [251, 135, 97], [254, 194, 135], [252, 253, 191],
];

function magma(v) {
  const x = Math.min(Math.max(v, 0), 1) * (MAGMA.length - 1);
  const i = Math.min(Math.floor(x), MAGMA.length - 2);
  const f = x - i;
  const a = MAGMA[i];
  const b = MAGMA[i + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}

/** In-place radix-2 FFT of (re, im), length a power of two. */
function fft(re, im) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * Math.PI) / len;
    const wr = Math.cos(ang);
    const wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1;
      let ci = 0;
      for (let k = 0; k < len / 2; k++) {
        const a = i + k;
        const b = a + len / 2;
        const tr = re[b] * cr - im[b] * ci;
        const ti = re[b] * ci + im[b] * cr;
        re[b] = re[a] - tr;
        im[b] = im[a] - ti;
        re[a] += tr;
        im[a] += ti;
        const nr = cr * wr - ci * wi;
        ci = cr * wi + ci * wr;
        cr = nr;
      }
    }
  }
}

/** A log-frequency magma spectrogram of an AudioBuffer on a canvas; returns {canvas, fLo, fHi, duration}. */
export function computeMap(buffer, { width = 1000, rows = 300, fLo = 40, fHi = 16000, spanDb = 80 } = {}) {
  const sr = buffer.sampleRate;
  const n = buffer.length;
  const mono = new Float32Array(n);
  for (let c = 0; c < buffer.numberOfChannels; c++) {
    const d = buffer.getChannelData(c);
    for (let i = 0; i < n; i++) mono[i] += d[i] / buffer.numberOfChannels;
  }
  fHi = Math.min(fHi, 0.48 * sr);
  const size = sr > 60000 ? 8192 : 4096;
  const cols = Math.max(1, Math.min(width, Math.floor(n / 64)));
  const hop = n / cols;
  const win = new Float32Array(size).map((_, i) => 0.5 - 0.5 * Math.cos((2 * Math.PI * i) / (size - 1)));
  // Each output row reads the bin at its log-spaced frequency, interpolated between neighbouring bins.
  const rowBin = new Float32Array(rows).map((_, r) => (fLo * (fHi / fLo) ** (1 - r / (rows - 1)) * size) / sr);
  const db = new Float32Array(cols * rows);
  const re = new Float64Array(size);
  const im = new Float64Array(size);
  let peak = -Infinity;
  for (let c = 0; c < cols; c++) {
    const centre = Math.round((c + 0.5) * hop);
    re.fill(0);
    im.fill(0);
    for (let i = 0; i < size; i++) {
      const k = centre - size / 2 + i;
      if (k >= 0 && k < n) re[i] = mono[k] * win[i];
    }
    fft(re, im);
    for (let r = 0; r < rows; r++) {
      const b = rowBin[r];
      const i = Math.floor(b);
      const f = b - i;
      const p0 = re[i] * re[i] + im[i] * im[i];
      const p1 = re[i + 1] * re[i + 1] + im[i + 1] * im[i + 1];
      const v = 10 * Math.log10(p0 * (1 - f) + p1 * f + 1e-30);
      db[r * cols + c] = v;
      if (v > peak) peak = v;
    }
  }
  const canvas = document.createElement('canvas');
  canvas.width = cols;
  canvas.height = rows;
  const g = canvas.getContext('2d');
  const img = g.createImageData(cols, rows);
  for (let i = 0; i < cols * rows; i++) {
    const [r, gg, b] = magma((db[i] - peak + spanDb) / spanDb);
    img.data[4 * i] = r;
    img.data[4 * i + 1] = gg;
    img.data[4 * i + 2] = b;
    img.data[4 * i + 3] = 255;
  }
  g.putImageData(img, 0, 0);
  canvas.className = 'computed';
  return { canvas, fLo, fHi, duration: n / sr };
}
