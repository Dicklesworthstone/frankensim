// Time-frequency maps. The grader's images (bare webp: x = time, 4 px per 10 ms frame; y = log frequency from f_lo at the
// bottom to f_hi at the top) get frequency and time ticks and a playhead overlaid. Clips the grader never saw (references,
// work-in-progress renders) get the same kind of map computed here: Hann STFT, log-frequency rows, 80 dB under the clip's
// own peak, magma colours.
import { h } from './app.js';

const TICKS_HZ = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];

function hzLabel(f) {
  return f >= 1000 ? `${f / 1000}k` : String(f);
}

/** Wrap a map element (img or canvas) with ticks and a playhead that follows `player`; a click seeks. */
export function mapFrame(mapEl, { fLo, fHi, duration, player, label }) {
  const ticks = h('div.ticks');
  for (const f of TICKS_HZ) {
    if (f <= fLo * 1.05 || f >= fHi * 0.95) continue;
    const y = 100 * (1 - Math.log(f / fLo) / Math.log(fHi / fLo));
    ticks.append(h('span.ftick', { style: { top: `${y}%` } }, hzLabel(f)));
  }
  const step = duration > 6 ? 1 : 0.5;
  for (let t = step; t < duration - 0.05; t += step) {
    ticks.append(h('span.ttick', { style: { left: `${(100 * t) / duration}%` } }, `${t.toFixed(t % 1 ? 1 : 0)} s`));
  }
  const head = h('div.playhead');
  const frame = h('div.mapframe', null, mapEl, ticks, head, label ? h('span.maplabel', null, label) : null);
  if (player) {
    const move = p => {
      head.style.left = `${(100 * Math.min(p.position, duration)) / duration}%`;
      head.style.opacity = p.playing || p.position > 0 ? '1' : '0';
    };
    player.on(move);
    frame.addEventListener('click', e => {
      const r = frame.getBoundingClientRect();
      player.seek(((e.clientX - r.left) / r.width) * duration);
      player.play();
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
