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
// Strandbar part with the Woodbury bridge solve; the sustained probe at the full cap gives 1.09-1.12,
// 1.27-1.32, 1.67-1.82, 2.33-2.47). Only used to pick which larger board is worth fetching and at
// which cap to time it; the choice itself is always a measurement on this device.
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
// coordinates, 0.43-0.50 on 64, 0.48-0.49 on 48, 0.54-0.64 on 96, 0.70-0.73 on 128 (caps 700-1800,
// medians of 4-6 interleaved rounds). Only a first guess: the cap is re-measured.
export const FIXED_COST = { 32: 0.38, 48: 0.48, 64: 0.46, 96: 0.58, 128: 0.7 };
// the mode cap at which a probe measured at `realtime` (cap `from`, cost share `a` fixed) reaches `target`
export function capFor(realtime, target, full, a, from = full) {
  const share = (realtime / target) * (a + ((1 - a) * from) / full); // cost(cap)/cost(full) wanted
  return Math.floor((full * (share - a)) / (1 - a));
}

// The lowest starting mode cap at which each fetched board still sounds better than the next smaller
// board at any cap. Dense player rolls (Giant Steps, Character Select, Bird blues; 20 s each) rendered
// against 128 coordinates at the full cap 2400, mean |dB| change of 1/3-octave band levels:
//   128: 0.16-0.20 at cap 1600, 0.26-0.37 at 1200
//    96: 0.71-0.72 at 2400, 0.76-0.88 at 1200, 0.77-0.92 at 1000, 0.85-1.16 at 800
//    64: 1.20-1.38 at 2400, 1.23-1.39 at 1200, 1.23-1.37 at 1000, 1.28-1.54 at 800
//    48: 1.71-1.75 at 2400, 1.73-1.78 at 1000, 1.75-1.88 at 800;   32: 2.15-2.47 at 2400, 2.12-2.42 at 1000
// At its floor each board beats the next smaller one in every song by 0.34 dB or more (96 at cap 800
// by only 0.15 on Giant Steps). At these caps the governor can still buy back about 1.2x by lowering
// the cap to its own floor, although the board-sized work does not shrink with the cap.
export const MIN_CAP = { 128: 1200, 96: 1000, 64: 1000, 48: 1000 };

// The lanes' choice of board AND mode cap, for `target` headroom on the sustained probe. The board is
// the audible lever and the cap barely is above the floors of MIN_CAP. So the boards are ranked by
// fidelity (`sizes`, largest first), each planned at the cap the embedded board's probe predicts for
// it (RELATIVE_COST, FIXED_COST). A board whose predicted cap is below its floor is not fetched; the
// best one that keeps its floor is fetched and timed once, at that cap, and taken if it reaches
// `target` (at a lowered cap, within 10 %). If it is short, its cap is corrected once along its cost
// line (at most a 25 % step, not below its floor) and taken if the estimate is within 10 % of `target`;
// if a 25 % step is not enough but its own timing still puts the cap at or above its floor, it is timed
// once more there (as the fit below). The plans follow the measured / predicted speed: after a
// shortfall the next smaller board is planned, after a faster run than predicted a better board that
// now keeps its floor is tried too and taken if it passes. At most two boards are fetched (about 0.3 s to build and time 64, 1.2 s 96).
// Without one, the embedded board, with the cap lowered as far as it needs (not below `floor`; floor >=
// the full cap keeps the full cap): a lowered cap is measured once more on a fresh engine; if that is
// still more than 10 % short, it is corrected once from the two points. make(bytes|null, coords, cap?)
// -> engine with that mode cap (per-key budget in proportion; no cap = the engine's own). `full`
// defaults to the first engine's maxDynamicModes. Returns {coordinates, bytes, cap, realtime}:
// realtime = the probe's speed at that cap (measured, or the estimate after a correction).
export async function chooseQuality({ make, measure = probeSustained, target, full = 0, sizes = BOARD_SIZES, minCap = MIN_CAP, floor = 400, log = () => {} }) {
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
  // A board's plan: the cap the probes so far predict for `target` (drift = measured / predicted speed
  // of the last board timed: load that changed since the first probe, or a first probe that ran before
  // the wasm was fully optimized), and its floor.
  let drift = 1;
  const plan = (n) => {
    const a = FIXED_COST[n] ?? 0.5;
    const share = (c) => a + ((1 - a) * c) / full; // cost(cap) / cost(full)
    const atFull = (drift * realtime32) / (RELATIVE_COST[n] || n / 32);
    const cap = atFull >= target ? full : Math.min(full, capFor(atFull, target, full, a));
    return { a, share, atFull, cap, lowest: Math.min(full, Math.max(floor, minCap[n] ?? full / 2)) };
  };
  let pick = null, pickAt = sizes.length, tries = 0;
  const tried = new Set();
  while (tries < 2) {
    // the best-sounding board not yet timed that is predicted to keep its floor (and beats the pick)
    const i = sizes.findIndex((n, k) => k < pickAt && !tried.has(n) && plan(n).cap >= plan(n).lowest);
    if (i < 0) break;
    const n = sizes[i];
    const { a, share, atFull, cap, lowest } = plan(n);
    tried.add(n);
    tries += 1;
    try {
      const bytes = await fetchBoard(n);
      const realtime = timed(bytes, n, cap);
      log(`board ${n}, cap ${cap}: ${realtime.toFixed(2)}x (predicted ${(atFull / share(cap)).toFixed(2)})`);
      drift *= (realtime * share(cap)) / atFull;
      if (realtime >= target || (cap < full && realtime >= 0.9 * target)) {
        // taken; a faster run than predicted may bring a better board within reach (the loop's next pass)
        pick = { coordinates: n, bytes, cap, realtime };
        pickAt = i;
        continue;
      }
      const wanted = capFor(realtime, target, full, a, cap); // where this board's own timing puts the cap
      const lower = Math.max(lowest, Math.ceil(0.75 * cap), Math.min(cap - 1, wanted));
      const estimate = (realtime * share(cap)) / share(lower);
      if (lower < cap && estimate >= 0.9 * target) {
        log(`board ${n}, cap ${lower}: ${estimate.toFixed(2)}x (estimated)`);
        pick = { coordinates: n, bytes, cap: lower, realtime: estimate };
        pickAt = i;
      } else if (wanted >= lowest && wanted < cap) {
        // further short than a 25 % step covers (a slow moment or a pessimistic first probe): time it
        // once more where its own timing puts the cap, as the embedded board's fit does
        const again = timed(bytes, n, wanted);
        log(`board ${n}, cap ${wanted}: ${again.toFixed(2)}x`);
        drift *= (again * share(wanted)) / (realtime * share(cap));
        if (again >= 0.9 * target) {
          pick = { coordinates: n, bytes, cap: wanted, realtime: again };
          pickAt = i;
        }
      }
    } catch (error) {
      log(`board ${n}: ${error && error.message ? error.message : error}`);
    }
  }
  return pick || { coordinates: 0, bytes: null, ...fit(null, 0, realtime32) };
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
