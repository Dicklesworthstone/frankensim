// Loudness and peak of each live kit piece's standard stroke on every grid a page may pick, measured
// with the kit wasm itself, for public/live/kit/levels.js (grid-independent output level).
// usage: node tools/kit_levels.mjs KIT_DIR > KIT_DIR/engine_levels.json   (engine-kit/build-wasm.sh)
//
// The standard stroke is the live kit's: src/live/kitStage.ts hitCentre() at the piece's default
// stick speed (the player's STICK). Loudness is the largest RMS over a 400 ms window (20 ms hop) of
// the first second at unit gain [Pa]; peak is the largest |sample| [Pa]. The table records the wasm's
// SHA-256, and the Worker uses it only beside that exact wasm.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

const dir = process.argv[2];
if (!dir) { console.error('usage: node tools/kit_levels.mjs KIT_DIR'); process.exit(2); }
const wasm = readFileSync(`${dir}/engine_bg.wasm`);
const mod = await import(pathToFileURL(`${dir}/engine.js`).href);
mod.initSync({ module: wasm });

// every grid a page asks for (worker default 2.5-8, Strandbar from 2, the player's restarts up to 10)
const SCALES = [2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10];
const REFERENCE = 5;
const SPEED = { kick: 3.0, snare: 4.5, 'snare-off': 4.5, 'tom-12': 3.5, 'tom-14': 3.5, 'hihat-closed': 2.5, 'hihat-open': 2.5, crash: 6.0, ride: 2.5 };
const spot = (name, R) => (name === 'kick' ? [0, 0.02] : name.startsWith('hihat') || name === 'crash' ? [0.85 * R, 0] : name === 'ride' ? [0.6 * R, 0] : [0.35 * R, 0.1 * R]);
const FRAMES = 48000, WINDOW = 19200, HOP = 960;

const pieces = {};
for (const name of JSON.parse(mod.pieces())) {
  if (!(name in SPEED)) continue;
  const rows = {};
  for (const scale of SCALES) {
    let p;
    try { p = new mod.KitPiece(name, scale); } catch { continue; } // grid too fine for 48 kHz
    const [x, y] = spot(name, p.radius());
    p.strike(x, y, SPEED[name]);
    const cs = new Float64Array(FRAMES + 1);
    let peak = 0;
    for (let k = 0; k < FRAMES; k += 128) {
      const out = p.render(128, 1);
      for (let i = 0; i < out.length; i++) { const v = out[i]; cs[k + i + 1] = cs[k + i] + v * v; peak = Math.max(peak, Math.abs(v)); }
    }
    p.free();
    let level = 0;
    for (let i = 0; i + WINDOW <= FRAMES; i += HOP) level = Math.max(level, Math.sqrt((cs[i + WINDOW] - cs[i]) / WINDOW));
    rows[scale] = { level: Number(level.toPrecision(6)), peak: Number(peak.toPrecision(6)) };
  }
  pieces[name] = rows;
  console.error(`${name}: ${Object.keys(rows).length} grids`);
}
console.log(JSON.stringify({
  about: 'Standard-stroke loudness (largest 400 ms RMS, Pa) and peak (Pa) per live kit piece and grid scale; see public/live/kit/levels.js',
  wasmSha256: createHash('sha256').update(wasm).digest('hex'),
  reference: REFERENCE,
  windowMs: (WINDOW * 1000) / 48000,
  pieces,
}, null, 1));
