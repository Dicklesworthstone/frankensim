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
