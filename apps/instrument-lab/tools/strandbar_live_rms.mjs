// RMS of the live engines playing one Strandbar loop, for level matching against the studio render.
// usage: node tools/strandbar_live_rms.mjs EVENTS.json PUBLIC_DIR [KIT_GRID_SCALE] [LANES_PREFIX]
// EVENTS.json: {loopFrames, tailFrames, piano:[{offset,on?,off?,speed?}], kick:[{offset,x,y,speed}], hihat:[...]}
// Prints {piano, kick, hihat} RMS in each lane's own units (piano: engine output; kit: Pa).
// With LANES_PREFIX, each lane's samples are also written to LANES_PREFIX-<lane>.f32 (little-endian f32),
// so the caller can measure the peak of the lanes summed at their gains.
import { readFileSync, writeFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

const [evPath, pub, scaleArg, lanesPrefix] = process.argv.slice(2);
const ev = JSON.parse(readFileSync(evPath, 'utf8'));
const scale = Number(scaleArg) || 5;
const load = async (dir) => {
  const mod = await import(pathToFileURL(`${pub}/live/${dir}/engine.js`).href);
  mod.initSync({ module: readFileSync(`${pub}/live/${dir}/engine_bg.wasm`) });
  return mod;
};

function play(render, apply, events, frames) {
  const out = new Float32Array(frames);
  const q = [...events].sort((a, b) => a.offset - b.offset);
  let at = 0, k = 0;
  while (at < frames) {
    while (k < q.length && q[k].offset <= at) apply(q[k++]);
    const n = Math.min(128, frames - at, k < q.length ? q[k].offset - at : Infinity);
    out.set(render(n), at);
    at += n;
  }
  return out;
}
const rms = (x, from, to) => {
  let s = 0;
  for (let i = from; i < to; i++) s += x[i] * x[i];
  return Math.sqrt(s / Math.max(1, to - from));
};

const frames = ev.loopFrames + ev.tailFrames;
const result = { scale };
const keep = (id, y) => lanesPrefix && writeFileSync(`${lanesPrefix}-${id}.f32`, Buffer.from(y.buffer, y.byteOffset, y.byteLength));
const kit = await load('kit');
for (const [id, piece] of [['kick', 'kick'], ['hihat', 'hihat-closed']]) {
  const p = new kit.KitPiece(piece, scale);
  const y = play((n) => p.render(n, 1), (e) => p.strike(e.x, e.y, e.speed), ev[id], frames);
  result[id] = rms(y, 0, frames);
  keep(id, y);
  p.free();
}
const pianoMod = await load('piano');
const pe = new pianoMod.PianoEngine(48000, 0, 0, 0, 0);
const t0 = performance.now();
const y = play((n) => pe.render(n), (e) => (e.on !== undefined ? pe.noteOn(e.on, e.speed) : pe.noteOff(e.off)), ev.piano, frames);
result.piano = rms(y, 0, frames);
keep('piano', y);
result.pianoRealtime = (frames / 48000) / ((performance.now() - t0) / 1000);
result.pianoInfo = JSON.parse(pe.info());
console.log(JSON.stringify(result));
