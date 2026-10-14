// Soundboard resolution for the live piano (shared by worker.js and lane-worker.js).
//
// The number of soundboard coordinates is the audible lever of the live piano. Measured against
// the studio render of the same Strandbar score (grand_piano, 128 coordinates), octave balance
// error / 2-8 kHz share at the attacks: 32 -> 2.95 dB / -6.2 dB, 48 -> 2.66 / -3.4,
// 64 -> 1.28 / -0.7, 96 -> 0.77 / -1.7, 128 -> 0.36 / -1.3 (engine-piano README of the patch).
// The 32-coordinate board is embedded in the wasm; larger bridge-informed reductions of the same
// (0, 2400] Hz FE slice (examples/gen_board.rs, spec N,N/2,200,400,800,1600,2300) are fetched
// from ./boards/ only when this device can compute them fast enough.
export const BOARD_SIZES = [128, 96, 64, 48];
export const boardUrl = (n) => new URL(`./boards/model-d-board-2400-r${n}.bin`, import.meta.url);

// Cost of the same score relative to the embedded 32-coordinate board (V8, measured on the
// Strandbar part with the Woodbury bridge solve). Only used to pick which larger board is worth
// fetching and measuring; the choice itself is always a measurement on this device.
export const RELATIVE_COST = { 32: 1, 48: 1.12, 64: 1.3, 96: 1.75, 128: 2.35 };

// Speed probe on a throwaway engine: twelve keys struck together (the Strandbar part's densest
// moment: 12 sounding keys, ~2400 dynamic string modes), 512 frames to settle, then 48 blocks of
// 128 frames timed. Returns audio seconds per wall second.
export const PROBE_KEYS = [38, 45, 50, 54, 56, 57, 59, 61, 63, 64, 66, 73];
export function probeRealtime(engine) {
  for (const k of PROBE_KEYS) engine.noteOn(k, 3.5);
  engine.render(512);
  const t0 = performance.now();
  for (let i = 0; i < 48; i++) engine.render(128);
  return (48 * 128 * 1000) / 48000 / (performance.now() - t0);
}

// The lanes' probe: what the piano costs while music plays. The chord above is held and never damped,
// so it ran 1.4-1.8x faster than real scores at the same dynamic modes (Giant Steps, Character Select
// replayed in V8); the dampers of released keys still ringing and the strikes that keep the cap
// binding are a large part of the cost. Here that cluster is struck and released at once (it rings
// on its dampers), and a comping pattern keeps striking: each beat (10240 frames, ~280 bpm) a bass
// key and a four-key chord, two upper keys half a beat later, both released at 3/4 beat. After 4096
// frames the next 12288 are timed. On 32 and 64 coordinates it came within 3-8 % of Giant Steps'
// sustained speed. Pass a fresh engine (the strikes would otherwise land on ringing keys).
const COMP_BASS = [36, 41, 43, 38, 45, 40];
const COMP_CHORDS = [[55, 59, 62, 66], [57, 60, 64, 67], [54, 58, 61, 65], [56, 60, 63, 69], [59, 62, 66, 71]];
const COMP_UPPER = [[74, 78], [76, 79], [73, 77]];
export function probeSustained(engine, { beat = 10240, settle = 4096, timed = 12288 } = {}) {
  const half = Math.round(beat / 256) * 128, release = Math.round((beat * 0.75) / 128) * 128;
  let k = 0, chord = null, bass = null, upper = null, t0 = 0;
  for (let f = 0; f < settle + timed; f += 128) {
    if (f === 0) for (const x of PROBE_KEYS) engine.noteOn(x, 3);
    if (f === 2048) for (const x of PROBE_KEYS) engine.noteOff(x);
    const p = f % beat;
    if (p === 0) {
      if (bass !== null) engine.noteOff(bass);
      bass = COMP_BASS[k % COMP_BASS.length];
      engine.noteOn(bass, 3);
      chord = COMP_CHORDS[k % COMP_CHORDS.length];
      for (const x of chord) engine.noteOn(x, 2.6);
      k += 1;
    } else if (p === half) {
      upper = COMP_UPPER[k % COMP_UPPER.length];
      for (const x of upper) engine.noteOn(x, 2.4);
    } else if (p === release) {
      for (const x of chord) engine.noteOff(x);
      if (upper) for (const x of upper) engine.noteOff(x);
      upper = null;
    }
    if (f === settle) t0 = performance.now();
    engine.render(128);
  }
  return (timed * 1000) / 48000 / (performance.now() - t0);
}

// Share of the probe's cost that does not shrink with the mode cap (board, radiation, per-block work):
// cost(cap) / cost(full) = a + (1 - a) cap / full. Measured on the probe in V8: 0.36-0.39 on 32
// coordinates, 0.44-0.49 on 64; 96 and 128 extrapolated. Only a first guess: the cap is re-measured.
export const FIXED_COST = { 32: 0.38, 48: 0.42, 64: 0.46, 96: 0.52, 128: 0.58 };
// the mode cap at which a probe measured at `realtime` (cap `from`, cost share `a` fixed) reaches `target`
export function capFor(realtime, target, full, a, from = full) {
  const share = (realtime / target) * (a + ((1 - a) * from) / full); // cost(cap)/cost(full) wanted
  return Math.floor((full * (share - a)) / (1 - a));
}

// The lanes' choice of board AND mode cap, for `target` headroom on the sustained probe. The board is
// the audible lever and the cap barely is: against a render at 128 coordinates and the full cap
// (Giant Steps, Character Select; mean |dB| change of 1/3-octave band levels), 32 coordinates are off
// by 2.15-2.47 dB at any cap from 1600 up, 48 by 1.71-1.77, 64 by 1.24-1.38 at caps 2400 down to 1200
// (and 64 at cap 1200 ran faster than 32 at 2400), 96 by 0.71. So a larger board wins while it keeps
// at least `minShare` of the full cap; otherwise the embedded board, with the cap lowered as far as
// it needs (not below `floor`; floor >= the full cap keeps the full cap). A lowered cap is measured
// once more on a fresh engine; if that is still more than 10 % short, it is corrected once from the
// two points. make(bytes|null, coords, cap?) -> engine with that mode cap (per-key budget in
// proportion; no cap = the engine's own). `full` defaults to the first engine's maxDynamicModes.
// Returns {coordinates, bytes, cap, realtime}: realtime = the probe's speed at that cap (measured,
// or the two-point estimate after a correction).
export async function chooseQuality({ make, measure = probeSustained, target, full = 0, sizes = BOARD_SIZES, minShare = 0.5, floor = 400, log = () => {} }) {
  const timed = (bytes, coords, cap) => {
    const engine = make(bytes, coords, cap);
    try {
      if (!full) full = JSON.parse(engine.info()).maxDynamicModes;
      return measure(engine);
    } finally {
      engine.free();
    }
  };
  const fit = (bytes, coords, realtime) => {
    if (realtime >= target || floor >= full) return { cap: full, realtime };
    const a = FIXED_COST[coords || 32] ?? 0.5;
    const cap = Math.max(floor, capFor(realtime, target, full, a));
    const measured = timed(bytes, coords, cap);
    log(`board ${coords || 32}, cap ${cap}: ${measured.toFixed(2)}x`);
    if (measured >= 0.9 * target || cap <= floor) return { cap, realtime: measured };
    // still short: the fixed share from the two points (cost affine in the cap), clamped against
    // timing noise, and at most a 25 % step, as the governor's
    const ratio = realtime / measured; // cost(cap) / cost(full)
    const fixed = Math.min(0.65, Math.max(0.25, (ratio - cap / full) / (1 - cap / full)));
    const lower = Math.max(floor, Math.ceil(0.75 * cap), Math.min(cap - 1, capFor(measured, target, full, fixed, cap)));
    const at = (c) => fixed + ((1 - fixed) * c) / full;
    return { cap: lower, realtime: (measured * at(cap)) / at(lower) };
  };
  const realtime32 = timed(null, 0, full || undefined);
  log(`board 32: ${realtime32.toFixed(2)}x`);
  let tries = 0;
  for (const n of sizes) {
    if (tries >= 2) break;
    const a = FIXED_COST[n] ?? 0.5;
    // even at minShare of the cap this board would not reach the target: do not fetch it
    if (realtime32 / (RELATIVE_COST[n] || n / 32) < target * (a + (1 - a) * minShare)) continue;
    tries += 1;
    try {
      const bytes = await fetchBoard(n);
      const realtime = timed(bytes, n, full);
      log(`board ${n}: ${realtime.toFixed(2)}x`);
      if (realtime >= target) return { coordinates: n, bytes, cap: full, realtime };
      if (capFor(realtime, target, full, a) >= minShare * full) return { coordinates: n, bytes, ...fit(bytes, n, realtime) };
    } catch (error) {
      log(`board ${n}: ${error && error.message ? error.message : error}`);
    }
  }
  return { coordinates: 0, bytes: null, ...fit(null, 0, realtime32) };
}

export async function fetchBoard(n) {
  const response = await fetch(boardUrl(n));
  if (!response.ok) throw new Error(`board r${n}: HTTP ${response.status}`);
  return new Uint8Array(await response.arrayBuffer());
}

// make(bytes|null, coords) -> engine; measure(engine) -> real-time factor of a probe.
// Returns {engine, coordinates, realtime}: the largest board whose probe runs at least `target`
// times faster than real time. Start-up cost: the embedded board plus at most two fetched ones.
export async function chooseBoard({ make, measure = probeRealtime, target, sizes = BOARD_SIZES, log = () => {} }) {
  const base = make(null, 0);
  const realtime32 = measure(base);
  log(`board 32: ${realtime32.toFixed(2)}x`);
  let tries = 0;
  for (const n of sizes) {
    if (tries >= 2) break;
    if (realtime32 / (RELATIVE_COST[n] || n / 32) < target) continue;
    tries += 1;
    let engine = null;
    try {
      const bytes = await fetchBoard(n);
      engine = make(bytes, n);
      const realtime = measure(engine);
      log(`board ${n}: ${realtime.toFixed(2)}x`);
      if (realtime >= target) {
        base.free();
        return { engine, coordinates: n, realtime, bytes };
      }
    } catch (error) {
      log(`board ${n}: ${error && error.message ? error.message : error}`);
    }
    engine?.free();
  }
  return { engine: base, coordinates: 0, realtime: realtime32, bytes: null };
}
