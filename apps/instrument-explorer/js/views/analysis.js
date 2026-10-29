// One phrase, everything about it on one page: the simulation and the real recording side by side (maps, then every
// measure as a lane on the same time axis), a panel that follows the playhead note by note, the per-note table, the
// long-term spectrum, and how both sounds were made and measured. Everything moves with the player.
import { h, s, fmt, link, getJSON, DATA, empty } from '../app.js';
import { controls } from '../audio.js';
import { mapFrame, pitchClass, unpitched } from '../spectro.js';
import { phrasePlayer, KIND_LABEL, pitchTally } from './phrase.js';

let titleText = 'Phrase';
export function title() {
  return titleText;
}

// Each lane: which frame series, how to scale it, and what it means in a sentence.
const LANES = [
  { key: 'cents', label: 'Pitch', unit: '¢ from the note asked for', lo: -100, hi: 100, ticks: [-100, -50, 0, 50, 100], clip: true,
    help: 'How far the sounding pitch is from the note asked for, every 10 ms. Triangles at the edge mean it is off the chart (a different note).' },
  { key: 'level', label: 'Loudness', unit: 'dB (A-weighted)', lo: -60, hi: 0, ticks: [-60, -40, -20, 0],
    help: 'A-weighted level under each clip\'s own peak: the shape of every attack, swell and decay.' },
  { key: 'centroid', label: 'Brightness', unit: 'Hz (spectral centroid)', lo: 100, hi: 8000, log: true, ticks: [100, 300, 1000, 3000, 8000],
    help: 'The centre of mass of the spectrum. Higher is brighter.' },
  { key: 'hnr', label: 'Harmonicity', unit: 'dB (harmonics over noise)', lo: -5, hi: 45, ticks: [0, 15, 30, 45],
    help: 'How periodic the sound is. Real bowed and blown notes carry bow, breath and turbulence noise; a synthesiser-like tone sits far above them.' },
  { key: 'rough', label: 'Noise and roughness', unit: 'share of energy off the harmonics', lo: 0, hi: 1, ticks: [0, 0.5, 1],
    help: 'The fraction of the energy that is not at a harmonic of the note: breath and bow noise, buzz, and, at a note\'s start, scratch or distortion.' },
  { key: 'flux', label: 'Onsets', unit: 'spectral flux', lo: 0, hi: null, ticks: [],
    help: 'How fast the spectrum changes. Peaks mark attacks and transients; a real player\'s articulation shows here.' },
];

const BANDS = {
  // display bands for the tiles: close / some way off / far off. They describe the size of a difference; they are not
  // the grader's rulers (those come from real-vs-real spread and live on the instrument pages).
  timbre_distance_db: [3, 6],
  partials_diff_db: [4, 8],
  harmonicity_diff_db: [3, 8],
  noise_floor_diff_db: [3, 8],
  brightness_ratio: [0.15, 0.4],
  attack_ratio: [0.25, 1],
  onset_roughness_diff: [0.08, 0.2],
  speaks_ratio: [0.5, 2],
  wobble_diff_cents: [1.5, 4],
};

function median(v) {
  const a = v.filter(x => x !== null && x !== undefined).sort((p, q) => p - q);
  return a.length ? a[Math.floor(a.length / 2)] : null;
}

function band(key, v) {
  const b = BANDS[key];
  if (!b || v === null || v === undefined) return 'none';
  const m = key.endsWith('_ratio') ? Math.abs(Math.log(v)) : Math.abs(v);
  const lim = key.endsWith('_ratio') ? b.map(x => Math.log(1 + x)) : b;
  return m <= lim[0] ? 'in' : m <= lim[1] ? 'near' : 'out';
}

function signed(v, d = 1, unit = '') {
  if (v === null || v === undefined) return '–';
  const r = Number(v.toFixed(d)); // the sign follows the shown digits, so -0.04 reads ±0.0, not −0.0
  return `${r > 0 ? '+' : r < 0 ? '−' : '±'}${Math.abs(r).toFixed(d)}${unit}`;
}

function tiles(p, an) {
  const out = [];
  const t = pitchTally(p);
  const total = t.in + t.near + t.out;
  if (total) {
    out.push({ cls: t.in === total ? 'in' : t.out ? 'out' : 'near', big: `${t.in}/${total}`, label: 'notes in tune',
      sub: t.in === total ? 'every note within 10 cents'
        : [t.near ? `${t.near} within 50 cents` : null, t.out - t.refused ? `${t.out - t.refused} more than 50 cents off` : null,
          t.refused ? `${t.refused} refused by the solver` : null].filter(Boolean).join(', ') });
  }
  const sm = an.summary || {};
  if (sm.timbre_distance_db !== undefined) {
    out.push({ cls: band('timbre_distance_db', sm.timbre_distance_db), big: `${sm.timbre_distance_db.toFixed(1)} dB`, label: 'timbre distance',
      sub: 'median gap between note spectra (third-octave bands)' });
  }
  if (sm.partials_diff_db !== undefined) {
    out.push({ cls: band('partials_diff_db', sm.partials_diff_db), big: `${sm.partials_diff_db.toFixed(1)} dB`, label: 'harmonic levels',
      sub: 'mean gap over the first 12 harmonics' });
  }
  if (sm.harmonicity_diff_db !== undefined) {
    const v = sm.harmonicity_diff_db;
    out.push({ cls: band('harmonicity_diff_db', v), big: signed(v, 0, ' dB'), label: 'harmonicity',
      sub: v > 3 ? 'cleaner than the real instrument: too pure' : v < -3 ? 'noisier than the real instrument' : 'as periodic as the real instrument' });
  }
  if (sm.noise_floor_diff_db !== undefined) {
    const v = sm.noise_floor_diff_db;
    out.push({ cls: band('noise_floor_diff_db', v), big: signed(v, 0, ' dB'), label: 'noise between harmonics',
      sub: v < -3 ? 'too clean: less noise than the real one' : v > 3 ? 'more noise than the real one' : 'about the real noise floor' });
  }
  if (sm.brightness_ratio !== undefined) {
    const v = sm.brightness_ratio;
    out.push({ cls: band('brightness_ratio', v), big: `×${v.toFixed(2)}`, label: 'brightness',
      sub: v > 1.15 ? 'brighter than the recording' : v < 0.87 ? 'duller than the recording' : 'as bright as the recording' });
  }
  if (sm.attack_ratio !== undefined) {
    const v = sm.attack_ratio;
    out.push({ cls: band('attack_ratio', v), big: `×${v.toFixed(2)}`, label: 'attack time',
      sub: v > 1.25 ? 'slower to speak' : v < 0.8 ? 'faster, more abrupt starts' : 'speaks as fast as the real one' });
  }
  if (sm.speaks_ratio !== undefined) {
    const v = sm.speaks_ratio;
    const sim = median(an.notes.map(n => n.sim.speaks_ms));
    const ref = median(an.notes.map(n => (n.ref ? n.ref.speaks_ms : null)));
    out.push({ cls: band('speaks_ratio', v), big: `${sim} ms`, label: 'speaks in',
      sub: `real: ${ref} ms; ${v > 1.5 ? 'scrapes or wavers before the note settles' : v < 0.67 ? 'settles quicker than a player' : 'settles like a player'}` });
  }
  if (sm.onset_roughness_diff !== undefined) {
    const v = sm.onset_roughness_diff;
    out.push({ cls: band('onset_roughness_diff', v), big: signed(v * 100, 0, ' pt'), label: 'rough starts',
      sub: v > 0.08 ? 'scratchier or more distorted note starts' : v < -0.08 ? 'cleaner starts than the real one' : 'starts like the real one' });
  }
  if (sm.wobble_diff_cents !== undefined) {
    const v = sm.wobble_diff_cents;
    out.push({ cls: band('wobble_diff_cents', v), big: signed(v, 1, '¢'), label: 'wobble',
      sub: v < -1.5 ? 'steadier than any player: machine-like' : v > 4 ? 'less steady than the real one' : 'as steady as a player' });
  }
  return out;
}

function val(arr, i) {
  return arr && i >= 0 && i < arr.length ? arr[i] : null;
}

/** One measure lane: a canvas with the simulation (brass) and the recording (verdigris) on the shared time axis. */
function lane(spec, fr, notes, duration, hop) {
  const canvas = h('canvas', { 'aria-hidden': 'true' });
  // a tap shows the explanation in place: title tooltips never appear on touch screens
  const text = h('p.lanetext', { hidden: true }, spec.help);
  const help = h('button.lanehelp', { type: 'button', 'aria-label': `What is ${spec.label}?`, 'aria-expanded': 'false', title: spec.help }, '?');
  help.addEventListener('click', e => {
    e.stopPropagation();
    text.hidden = !text.hidden;
    help.setAttribute('aria-expanded', String(!text.hidden));
  });
  const wrap = h('div.lane', null, h('div.lanehead', null, h('b', null, spec.label), h('span', null, spec.unit), help), text, canvas);
  const sim = fr.sim[spec.key];
  const ref = fr.ref ? fr.ref[spec.key] : null;
  let hi = spec.hi;
  if (hi === null) {
    const all = [...(sim || []), ...(ref || [])].filter(v => v !== null).sort((a, b) => a - b);
    hi = all.length ? Math.max(1e-6, all[Math.floor(all.length * 0.995)] * 1.1) : 1;
  }
  const yOf = (v, H) => {
    let u = spec.log ? Math.log(v / spec.lo) / Math.log(hi / spec.lo) : (v - spec.lo) / (hi - spec.lo);
    return H - 4 - Math.min(Math.max(u, 0), 1) * (H - 8);
  };
  let active = 0;
  const draw = () => {
    const W = canvas.clientWidth;
    const H = canvas.clientHeight;
    if (!W || !H) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    const g = canvas.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    const css = getComputedStyle(document.documentElement);
    g.strokeStyle = css.getPropertyValue('--line').trim();
    g.fillStyle = css.getPropertyValue('--muted').trim();
    g.font = '10px ' + css.getPropertyValue('--mono');
    g.lineWidth = 1;
    let lastLabel = -Infinity;
    for (const tv of [...spec.ticks].sort((a, b) => yOf(b, H) - yOf(a, H))) { // bottom to top
      const y = Math.round(yOf(tv, H)) + 0.5;
      g.globalAlpha = tv === 0 && spec.key === 'cents' ? 0.9 : 0.45;
      g.beginPath();
      g.moveTo(0, y);
      g.lineTo(W, y);
      g.stroke();
      // kept inside the canvas (the top tick's label would be cut in half), and skipped when it would sit on the last
      const ly = Math.min(Math.max(y - 2, 10), H - 2);
      if (lastLabel - ly < 11) continue;
      lastLabel = ly;
      g.globalAlpha = 1;
      g.fillText(spec.log ? (tv >= 1000 ? `${tv / 1000}k` : String(tv)) : String(tv), W - 30, ly);
    }
    g.globalAlpha = 0.35;
    for (const n of notes) {
      const x = Math.round((n.at / duration) * W) + 0.5;
      g.beginPath();
      g.moveTo(x, 0);
      g.lineTo(x, H);
      g.stroke();
    }
    g.globalAlpha = 1;
    const series = [[ref, css.getPropertyValue('--real').trim(), 1], [sim, css.getPropertyValue('--sim').trim(), 0]];
    for (const [arr, color, side] of series) {
      if (!arr) continue;
      g.strokeStyle = color;
      g.fillStyle = color;
      g.lineWidth = side === active ? 2 : 1.2;
      g.globalAlpha = side === active ? 1 : 0.7;
      g.beginPath();
      let pen = false;
      for (let i = 0; i < arr.length; i++) {
        const v = arr[i];
        const t = i * hop;
        if (t > duration) break;
        if (v === null) { pen = false; continue; }
        const x = (t / duration) * W;
        if (spec.clip && (v > hi || v < spec.lo)) {
          pen = false;
          if (i % 6 === 0) { // a small triangle at the edge: off the chart
            const yy = v > hi ? 3 : H - 3;
            g.fillRect(x - 1, yy - 1, 2.5, 2.5);
          }
          continue;
        }
        const y = yOf(v, H);
        if (pen) g.lineTo(x, y);
        else g.moveTo(x, y);
        pen = true;
      }
      g.stroke();
    }
    g.globalAlpha = 1;
  };
  return { el: wrap, draw, setActive(i) { active = i; draw(); }, at(i) { return [val(sim, i), val(ref, i)]; }, spec };
}

function fmtLane(spec, v) {
  if (v === null || v === undefined) return '–';
  if (spec.key === 'cents') return Math.abs(v) >= 600 ? `${v > 0 ? '+' : '−'}${Math.round(Math.abs(v) / 100)} st` : `${v > 0 ? '+' : v < 0 ? '−' : '±'}${Math.abs(Math.round(v))}¢`;
  if (spec.key === 'centroid') return `${Math.round(v)} Hz`;
  if (spec.key === 'rough') return `${Math.round(v * 100)} %`;
  if (spec.key === 'flux') return v.toFixed(2);
  return `${v.toFixed(1)} dB`;
}

/** Grouped bars of the first harmonics' levels: simulation vs recording, for one note. */
function partialBars(simDb, refDb) {
  const n = Math.max((simDb || []).length, (refDb || []).length, 1);
  const W = 320;
  const H = 120;
  const bw = W / n;
  const y = v => Math.min(H, Math.max(0, (-v / 60) * H));
  const svg = s('svg.partials', { viewBox: `0 0 ${W} ${H + 16}`, role: 'img', 'aria-label': 'Harmonic levels of this note: simulation and recording' });
  for (const db of [0, -20, -40, -60]) {
    svg.append(s('line.grid', { x1: 0, x2: W, y1: y(db), y2: y(db) }), s('text.ax', { x: W - 2, y: y(db) - 2, 'text-anchor': 'end' }, `${db}`));
  }
  for (let k = 0; k < n; k++) {
    const x0 = k * bw;
    if (simDb && simDb[k] !== undefined) svg.append(s('rect.sim', { x: x0 + bw * 0.12, y: y(simDb[k]), width: bw * 0.36, height: H - y(simDb[k]) }));
    if (refDb && refDb[k] !== undefined) svg.append(s('rect.real', { x: x0 + bw * 0.52, y: y(refDb[k]), width: bw * 0.36, height: H - y(refDb[k]) }));
    if (n <= 16) svg.append(s('text.ax', { x: x0 + bw / 2, y: H + 12, 'text-anchor': 'middle' }, String(k + 1)));
  }
  return svg;
}

/** Long-term average spectrum, third-octave bands, both tracks. */
function ltasChart(lt) {
  // drawn at about the width it is shown at, so the axis text stays readable on a phone
  const W = (window.innerWidth || 1200) < 760 ? 380 : 760;
  const H = 220;
  const pad = { l: 36, r: 10, t: 10, b: 26 };
  const n = lt.bands_hz.length;
  const x = i => pad.l + (i / (n - 1)) * (W - pad.l - pad.r);
  const y = v => pad.t + Math.min(1, Math.max(0, -v / 60)) * (H - pad.t - pad.b);
  const svg = s('svg.ltas', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'Long-term average spectrum of the simulation and the recording' });
  for (const db of [0, -20, -40, -60]) {
    svg.append(s('line.grid', { x1: pad.l, x2: W - pad.r, y1: y(db), y2: y(db) }), s('text.ax', { x: pad.l - 4, y: y(db) + 3, 'text-anchor': 'end' }, `${db}`));
  }
  lt.bands_hz.forEach((f, i) => {
    if ([31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000].some(c => Math.abs(f / c - 1) < 0.03)) {
      svg.append(s('text.ax', { x: x(i), y: H - 8, 'text-anchor': 'middle' }, f >= 1000 ? `${Math.round(f / 1000)}k` : String(Math.round(f))));
    }
  });
  const path = arr => arr.map((v, i) => `${i ? 'L' : 'M'}${x(i).toFixed(1)},${y(v).toFixed(1)}`).join('');
  if (lt.ref) svg.append(s('path.real', { d: path(lt.ref) }));
  svg.append(s('path.sim', { d: path(lt.sim) }));
  return svg;
}

const NOTE_COLS = [
  ['cents', 'Pitch', v => (v === null || v === undefined ? '–' : Math.abs(v) >= 600 ? `${v > 0 ? '+' : '−'}${Math.round(Math.abs(v) / 100)} st` : signed(v, 0, '¢'))],
  ['attack_ms', 'Attack', v => (v === null || v === undefined ? '–' : `${v} ms`)],
  ['speaks_ms', 'Speaks in', (v, m) => (v === null || v === undefined ? '–' : m && m.speaks === false ? `>${v} ms` : `${v} ms`)],
  ['centroid_hz', 'Brightness', v => (v === null || v === undefined ? '–' : `${Math.round(v)} Hz`)],
  ['hnr_db', 'Harmonicity', v => (v === null || v === undefined ? '–' : `${v.toFixed(0)} dB`)],
  ['noise_floor_db', 'Noise floor', v => (v === null || v === undefined ? '–' : `${v.toFixed(0)} dB`)],
  ['onset_roughness', 'Rough start', v => (v === null || v === undefined ? '–' : `${Math.round(v * 100)} %`)],
  ['wobble_cents', 'Wobble', v => (v === null || v === undefined ? '–' : `${v.toFixed(1)}¢`)],
];

function noteTable(p, an, player) {
  const notes = p.notes || [];
  if (!notes.length) return null;
  const hasRef = an.notes.some(n => n.ref);
  const head = h('tr', null, h('th', null, 'Note'), NOTE_COLS.map(([, label]) => h('th', { colspan: hasRef ? 2 : 1 }, label)),
    hasRef ? h('th', null, 'Timbre gap') : null, h('th', null, ''));
  const sub = hasRef ? h('tr.sub', null, h('th'), NOTE_COLS.map(() => [h('th.sim', null, 'sim'), h('th.real', null, 'real')]), h('th'), h('th')) : null;
  const span = 1 + NOTE_COLS.length * (hasRef ? 2 : 1) + (hasRef ? 1 : 0);
  const rows = an.notes.map(row => {
    const n = notes[row.i];
    if (n.refused) {  // a silent slot: the solver's message stands in for the measurements
      return h('tr.n-refused', { 'data-i': String(row.i) }, h('td.note', null, n.note),
        h('td.refused', { colspan: span }, h('b', null, 'refused by the solver'), ` (not played): ${n.refused}`));
    }
    const play = h('button.mini', { type: 'button', 'aria-label': `Play ${n.note}` }, '▶');
    play.addEventListener('click', () => player.load().then(() => { player.seek(n.at); player.play(); }, () => {}));
    const tr = h(`tr.n-${pitchClass({ cents: row.sim.cents })}`, { 'data-i': String(row.i) }, h('td.note', null, n.note),
      NOTE_COLS.map(([k, , f]) => [h('td.sim', null, f(row.sim[k], row.sim)), hasRef ? h('td.real', null, f(row.ref ? row.ref[k] : null, row.ref)) : null]),
      hasRef ? h(`td.gap.c-${band('timbre_distance_db', row.timbre_distance_db)}`, null, row.timbre_distance_db === null || row.timbre_distance_db === undefined ? '–' : `${row.timbre_distance_db.toFixed(1)} dB`) : null,
      h('td', null, play));
    return tr;
  });
  return h('div.tablewrap', null, h('table.notetable', null, h('thead', null, head, sub), h('tbody', null, rows)));
}

const METHODS = [
  ['Pitch', 'Every 10 ms, the autocorrelation of a window holding three periods of the lowest pitch (85 ms at 50 Hz), divided by the window\'s own autocorrelation (Boersma\'s method, as in Praat), searched from its first zero crossing. Rumble below 25 Hz is filtered out first. A frame cannot always tell a period from twice it, so its octave is settled toward the note\'s own pitch measured over its whole steady part, never toward the note asked for (that would pull a wrong note nearer the target). The whole-note pitch checks its candidates against the spectrum: a period with no line at its own frequency whose lines all fall at multiples of k is k times too long, and one whose in-between harmonics sit 20 dB under the rest is a subharmonic. A note whose autocorrelation peak stays under 0.5 is marked "no steady pitch". Per note: the median over its steady part (after the first 120 ms).'],
  ['Loudness', 'A-weighted RMS level of each 43 ms window, in dB under the clip\'s own peak, so two clips at different gains compare by shape.'],
  ['Brightness', 'Spectral centroid above 50 Hz, power-weighted.'],
  ['Harmonicity', 'Harmonics-to-noise ratio from the peak r of the window-corrected autocorrelation used for pitch: 10 log10(r / (1 − r)).'],
  ['Noise and roughness', 'Share of the window\'s energy farther than 3 % (or 1.5 bins) from every harmonic of its own pitch, or of every note of a chord; where the frame has no reliable pitch, of the note\'s own steady pitch. Rough start: the same over a note\'s first 120 ms. Below about 70 Hz, and on low chords, the harmonics are closer together than this window can separate, so there is nothing between them to measure and no value is shown.'],
  ['Speaks in', 'How long a note takes to settle: from its start until the sound repeats at its own period (a chord: at the shortest period all its notes share) with a normalised correlation above 0.9, and keeps doing so for 50 ms. A scraping or multiply-slipping bow, a reed that has not locked and a breathy flue transient are all aperiodic. ">" marks a note that never settles before it ends.'],
  ['Noise floor', 'The median spectrum level midway between harmonics, under the strongest harmonic, over a note\'s steady part.'],
  ['Attack', 'Time from 10 % to 90 % of a note\'s peak amplitude (−20 dB to −1 dB).'],
  ['Wobble', 'Standard deviation of the pitch (cents) over the steady part after removing a straight-line drift, over the frames within 300 cents of the note\'s median (a stray mistracked frame is not vibrato): vibrato plus jitter. Real players are never perfectly steady; a note with fewer than half its frames that close has none.'],
  ['Harmonic levels', 'The first 16 harmonics of each note, from a long FFT of its steady part, in dB under the strongest.'],
  ['Timbre gap', 'Mean absolute difference of the two notes\' third-octave spectra after matching their overall level.'],
  ['Long-term spectrum', 'Third-octave average over the whole phrase, under its strongest band.'],
];

export async function render(main, r) {
  const [idx, phr] = await Promise.all([getJSON('index.json'), getJSON('phrases.json')]);
  const p = (phr.phrases || []).find(x => x.id === r.id);
  if (!p) {
    main.append(h('div.error', null, h('h2', null, 'No such phrase'), h('p', null, h('a', { href: '#/' }, 'Back to all phrases'))));
    return;
  }
  const fam = (idx.families || []).find(f => f.id === p.family) || { title: p.family };
  titleText = `${p.title} · ${fam.title}`;
  const an = p.analysis ? await getJSON(p.analysis).catch(() => null) : null;
  const player = phrasePlayer(p, fam.title);
  player.load().catch(() => {});
  const real = p.real || p.ref;
  const duration = p.duration_s;
  const notes = p.notes || [];

  const head = h('header.anahead', null,
    h('nav.crumbs', { 'aria-label': 'Breadcrumb' }, h('a', { href: '#/' }, 'Listen'), h('span', null, '›'),
      h('a', { href: link('f', p.family) }, fam.title), h('span', null, '›'), h('span', null, KIND_LABEL[p.kind] || p.kind)),
    h('h1', null, p.title),
    h('p.lede', null, real
      ? ['The simulation and a real recording, measured the same way and drawn on one time axis. Press play: every lane, the maps and the note panel follow the playhead.',
        h('span.kbdonly', null, ' Keys 1 and 2 switch what you hear without losing your place.')]
      : 'The simulation measured every 10 ms. There is no real recording of this passage to set beside it.'));

  const tileList = an ? tiles(p, an) : [];
  // two even rows on a wide screen (9 tiles: 5 + 4), rather than one full row and a stray pair
  const tileRow = an ? h('div.tiles', { style: `--tilecols: ${Math.ceil(tileList.length / 2)}` }, tileList.map(t => h(`div.tile.c-${t.cls}`, null, h('b', null, t.big), h('span.tl', null, t.label), h('span.ts', null, t.sub)))) : null;

  const transport = h('div.anatransport', null, controls(player, { names: real ? ['Simulation', p.real ? 'Recording' : 'Real notes'] : undefined, keys: Boolean(real) }));

  const fig = { fLo: p.f_lo_hz, fHi: p.f_hi_hz, duration, player };
  const maps = h('div.anamaps', null,
    mapFrame(h('img', { src: DATA + p.map, alt: `Spectrogram of the simulation: ${p.title}` }), { ...fig, notes, label: 'Simulation', side: 'sim' }));
  if (real && real.map) {
    maps.append(mapFrame(h('img', { src: DATA + real.map, alt: `Spectrogram of the recording: ${p.title}` }),
      { ...fig, notes: p.ref ? notes : null, label: p.real ? 'Recording' : 'Real notes', side: 'real' }));
  }

  // Lanes and the hover/playhead readout.
  let lanes = [];
  let stack = null;
  let hover = null;
  let readout = null;
  if (an) {
    const hop = an.hop_s;
    // a lane with nothing to draw (pitch on a chord, a drum groove) is left out, not shown empty
    const drawn = arr => Array.isArray(arr) && arr.some(v => v !== null);
    lanes = LANES.filter(L => drawn(an.frames.sim[L.key]) || drawn(an.frames.ref && an.frames.ref[L.key]))
      .map(L => lane(L, an.frames, notes, duration, hop));
    const head2 = h('div.playhead.big');
    hover = h('div.hoverline');
    readout = h('div.readout', { 'aria-live': 'off' });
    stack = h('div.lanestack', null, lanes.map(l => l.el), hover, head2, readout);
    const idxAt = t => Math.round(t / hop);
    const showReadout = (t, xFrac) => {
      const i = idxAt(t);
      const n = notes.find(nn => t >= nn.at && t < nn.at + nn.dur);
      readout.replaceChildren(h('div.rt', null, `${fmt.time(t)}${n ? ` · ${n.note}` : ''}`),
        ...lanes.map(l => {
          const [a, b] = l.at(i);
          return h('div.rr', null, h('span', null, l.spec.label), h('b.sim', null, fmtLane(l.spec, a)), real ? h('b.real', null, fmtLane(l.spec, b)) : null);
        }));
      // beside the pointer where it fits, else on its other side, else as far right as the lanes allow (a phone's
      // lanes are barely wider than the readout)
      const sw = stack.clientWidth;
      const w = readout.offsetWidth;
      const x = xFrac * sw;
      readout.style.left = `${x + 12 + w <= sw ? x + 12 : x - 12 - w >= 0 ? x - 12 - w : Math.max(4, sw - w - 4)}px`;
    };
    stack.addEventListener('pointermove', e => {
      const rc = stack.getBoundingClientRect();
      const fx = (e.clientX - rc.left) / rc.width;
      hover.style.left = `${fx * 100}%`;
      hover.style.opacity = '1';
      readout.style.opacity = '1';
      showReadout(fx * duration, fx);
    });
    stack.addEventListener('pointerleave', () => {
      hover.style.opacity = '0';
      readout.style.opacity = '0';
    });
    stack.addEventListener('click', e => {
      const rc = stack.getBoundingClientRect();
      const t = ((e.clientX - rc.left) / rc.width) * duration;
      player.load().then(() => { player.seek(t); player.play(); }, () => {});
    });
    player.on(q => {
      head2.style.left = `${(100 * Math.min(q.position, duration)) / duration}%`;
      head2.style.opacity = q.playing || q.position > 0 ? '1' : '0';
    }, stack);
  }

  // The panel that follows the playhead.
  const nowNote = h('div.nownote');
  const nowBody = h('div.nowbody');
  const aside = h('aside.nowpanel', { 'aria-label': 'The note under the playhead' }, h('h4', null, 'At the playhead'), nowNote, nowBody);
  let lastI = -2;
  let tableRows = [];
  // Stacked under the lanes (narrow screens) the panel is in the page's flow: a note's rows replacing "Press play" would
  // shove everything below it down mid-playback and back up at the end. It keeps the tallest height it has shown.
  const holdHeight = () => {
    if (!aside.isConnected || getComputedStyle(aside).position === 'sticky') return;
    if (aside.offsetHeight > (parseFloat(aside.style.minHeight) || 0)) aside.style.minHeight = `${aside.offsetHeight}px`;
  };
  const showNote = i => {
    if (i === lastI) return;
    lastI = i;
    for (const tr of tableRows) tr.classList.toggle('now', an && i >= 0 && an.notes[i] && tr.dataset.i === String(an.notes[i].i));
    if (!an || i < 0 || !an.notes[i]) {
      nowNote.replaceChildren(h('b', null, notes.length ? 'Press play' : 'Whole phrase'),
        h('span', null, notes.length ? `${fmt.plural(notes.length, p.kind === 'chords' ? 'chord' : 'note')}; the panel follows each one.` : 'One continuous render.'));
      nowBody.replaceChildren(an && an.summary && Object.keys(an.summary).length
        ? h('p.small.muted', null, 'The tiles above summarise the per-note differences; the table below lists every note.')
        : h('p.small.muted', null, `${matchMedia('(hover: none)').matches ? 'Touch' : 'Hover'} the lanes to read every measure at any moment.`));
      holdHeight();
      return;
    }
    const row = an.notes[i];
    const n = notes[row.i];
    if (n.refused) {
      nowNote.replaceChildren(h('b', null, n.note), h('span', null, [n.string ? `${n.string} string` : null,
        n.target_hz ? `asked ${fmt.num(n.target_hz, 1)} Hz` : null, 'silent'].filter(Boolean).join(' · ')));
      nowBody.replaceChildren(h('p.small', null, 'The solver refused this note part way through its render, rather than return a wrong answer: ', h('code', null, n.refused), '.'));
      holdHeight();
      return;
    }
    nowNote.replaceChildren(h('b', null, n.note), h('span', null, [n.string ? `${n.string} string` : null, n.target_hz ? `asked ${fmt.num(n.target_hz, 1)} Hz` : null,
      unpitched(n) ? 'no steady pitch' : n.f0_hz ? `sounds ${fmt.num(n.f0_hz, 1)} Hz` : null].filter(Boolean).join(' · ')));
    const line = (label, k, f) => h('div.nl', null, h('span', null, label), h('b.sim', null, f(row.sim[k], row.sim)), real ? h('b.real', null, f(row.ref ? row.ref[k] : null, row.ref)) : null);
    // replaceChildren takes nodes only (an array would be written out as text), and the narrow layout wants the rows
    // as direct children, so the list is flattened rather than wrapped.
    nowBody.replaceChildren(...[
      h('div.nlhead', null, h('span'), h('span.sim', null, 'simulation'), real ? h('span.real', null, 'recording') : null),
      NOTE_COLS.map(([k, label, f]) => line(label, k, f)),
      row.timbre_distance_db !== undefined && row.timbre_distance_db !== null ? h('div.nl.gap', null, h('span', null, 'Timbre gap'), h(`b.c-${band('timbre_distance_db', row.timbre_distance_db)}`, null, `${row.timbre_distance_db.toFixed(1)} dB`)) : null,
      row.sim.partials_db ? [h('h4', null, 'Harmonic levels (dB)'), partialBars(row.sim.partials_db, row.ref ? row.ref.partials_db : null)] : null,
    ].flat().filter(Boolean));
    holdHeight();
  };
  showNote(-1);
  // Once the page is on screen, size the panel for a full note (one with harmonic levels) before anything plays.
  if (an && an.notes.length) {
    requestAnimationFrame(() => {
      if (!aside.isConnected || lastI !== -1) return;
      showNote(Math.max(0, an.notes.findIndex(rw => rw.sim.partials_db)));
      showNote(-1);
    });
  }
  player.on(q => {
    const t = q.position;
    const i = q.playing || t > 0 ? an ? an.notes.findIndex(rw => t >= notes[rw.i].at && t < notes[rw.i].at + notes[rw.i].dur) : -1 : -1;
    if (i >= 0 || !q.playing) showNote(i >= 0 ? i : (q.position > 0 ? lastI : -1));
  }, aside);
  // the brass/verdigris emphasis follows what is being heard
  let shownSide = -1;
  player.on(q => {
    if (q.active === shownSide) return;
    shownSide = q.active;
    lanes.forEach(l => l.setActive(q.active));
  }, aside);

  const grid = h('div.anagrid', null, h('div.anamain', null, maps, stack, h('div.legendrow', null,
    h('span', null, h('span.sw.sim'), 'simulation'), real ? h('span', null, h('span.sw.real'), p.real ? 'recording' : 'real notes') : null,
    h('span.muted', null, matchMedia('(hover: none)').matches ? 'touch to read values · tap to play from there' : 'hover to read values · click to play from there'))), aside);

  const sections = [];
  const table = an ? noteTable(p, an, player) : null;
  if (table) tableRows = [...table.querySelectorAll('tbody tr')];
  if (table) sections.push(h('section.anasec', null, h('h2', null, 'Every note, measured'), h('p.muted.small', null, 'Each note\'s steady part (after its first 120 ms), simulation beside the recording. Click ▶ to hear it.'), table));
  if (an && an.ltas) {
    sections.push(h('section.anasec', null, h('h2', null, 'Where the energy sits'),
      h('p.muted.small', null, 'Long-term average spectrum of the whole phrase in third-octave bands, each under its own strongest band.'),
      ltasChart(an.ltas), an.quality ? h('p.small.muted', null,
        `Crest factor ${an.quality.sim.crest_db} dB (simulation)${an.quality.ref ? `, ${an.quality.ref.crest_db} dB (recording)` : ''}; clipped samples: ${an.quality.sim.clipped_samples}.`,
        an.envelope_correlation !== undefined && an.envelope_correlation !== null ? ` Loudness contours correlate at r = ${an.envelope_correlation.toFixed(2)}.` : '') : null));
  }
  const provenance = [h('div', null, h('h4', null, 'How the simulation was made'), h('p', null, p.how || ''))];
  if (p.ref) {
    const credits = new Map();
    for (const src of p.ref.sources || []) credits.set(`${src.licence}|${src.credit}`, src);
    provenance.push(h('div', null, h('h4', null, 'How the real phrase was made'), h('p', null, p.ref.how),
      h('ul.credits', null, [...credits.values()].map(c => h('li', null, c.credit || c.source, c.licence && c.licence !== c.credit ? h('span.muted', null, ` · ${c.licence}`) : null)))));
  } else if (p.real) {
    provenance.push(h('div', null, h('h4', null, 'The recording'), h('p', null, p.credit || 'A real recording of the same passage.')));
  }
  sections.push(h('section.anasec.howbox', null, provenance));
  sections.push(h('section.anasec', null, h('h2', null, 'How each measure is computed'),
    h('p.muted.small', null, 'Both tracks go through the same code (tools/phrases/analyze.py). The tile colours describe the size of a difference only: close, some way off, far off. The grader\'s rulers, built from the spread between real recordings, are on each instrument\'s page.'),
    h('dl.methods', null, METHODS.map(([k, v]) => [h('dt', null, k), h('dd', null, v)]))));

  main.append(h('div.anapage', null, head, tileRow, transport, grid, sections));
  if (!an) main.append(empty('The measurements for this phrase have not been computed yet.'));

  // Lanes draw once laid out, and again whenever their width changes.
  if (lanes.length && typeof ResizeObserver === 'function') {
    let pending = 0;
    new ResizeObserver(() => {
      cancelAnimationFrame(pending);
      pending = requestAnimationFrame(() => lanes.forEach(l => l.draw()));
    }).observe(stack);
  } else {
    requestAnimationFrame(() => lanes.forEach(l => l.draw()));
  }
}
