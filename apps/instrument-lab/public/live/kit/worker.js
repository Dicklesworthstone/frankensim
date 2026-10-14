// Live drum-kit piece worker (ES module Worker): owns one fs-vkplate piece (engine.js) and
// feeds one lane of fs-live-mix. Samples are microphone pressure [Pa] times a fixed gain, which the
// page sets for the 5x grid; on the grid picked here it is scaled so that a stroke is as loud as on
// 5x (levels.js, engine_levels.json beside the wasm).
//
// Main thread -> worker
//   {type:'init', piece, port, gain?, aheadFrames?, chunk?, scales?, minRealtime?, fieldHz?, fieldFrames?, bursts?}
//       chunk        most frames per audio message (default 128; lane.js)
//       fieldHz      fields per second of wall time (0: no fields)
//       fieldFrames  instead: one field every this many AUDIO frames (even spacing however far
//                    ahead the lane renders; for a score heard a second after it is computed)
//       bursts       false: no slow-motion burst after each stroke (default true when fields are on)
//       scales       grid coarsening factors to try, finest first (default [4, 5, 6, 7, 8])
//       minRealtime  the finest scale that runs at least this much faster than real time wins (1.7)
//   {type:'strike', x, y, speed, frame?}   x, y metres from the centre; speed m/s
//   {type:'sequence', startFrame, loopFrames, events:[{offset, x, y, speed}]}
//   {type:'score', events:[{frame, x, y, speed}]}   a fixed score, played once
//   {type:'stopScore'}                      drop the rest of the score (ringing continues)
//   {type:'resync', epoch}                  render for the clock as it is now (after a mix flush)
//   {type:'stop'}                           end the sequence (ringing continues)
//   {type:'gain', gain}                     a new page gain (scaled by the same levelFix)
// Worker -> main thread
//   {type:'ready', info, xy:Float32Array, h, radius, scale, realtime, levelFix, levels}
//       levelFix  the factor on the page's gain for this grid (1 on 5x); levels: 'ok', or 'none' when
//                 engine_levels.json is missing or was measured with another wasm (gain unscaled)
//   {type:'burst', frame, stepMs, count, v:Float32Array}    the first milliseconds after a stroke: `count`
//                                                          surface velocities `stepMs` apart, back to back
//   {type:'field', frame, epoch, v:Float32Array, kinetic, force}   ~fieldHz while the piece moves; force:
//                                                          the largest stick force since the last field
//   {type:'stats', realtime, renderedFrames, skippedFrames, lateEvents, asleep}
//   {type:'error', message}
import init, { KitPiece } from './engine.js';
import { createLane } from '../lane.js';
import { levelFix } from './levels.js';

let piece = null;
let lane = null;
let gain = 0.03; // = the page's gain times fix
let fix = 1;
let fieldEvery = 1000 / 30;
let fieldFrames = 0; // > 0: fields paced by audio frames instead of wall time
let fieldDue = 0;
let bursts = true;
let lastField = 0;
let lastStats = 0;
let forceMax = 0;
// right after a stroke: the surface every BURST_STEP samples for BURST_N snapshots (a slow-motion replay)
const BURST_STEP = 12;
const BURST_N = 32;
let burst = null;

const post = (m, t) => self.postMessage(m, t ?? []);

function bench(name, scale) {
  const p = new KitPiece(name, scale);
  p.strike(0.3 * p.radius(), 0, 3);
  p.render(256, 1); // first touch of every buffer
  const t0 = performance.now();
  for (let i = 0; i < 40; i++) p.render(128, 1);
  const ms = performance.now() - t0;
  p.free();
  return (40 * 128 * 1000) / 48000 / ms;
}

// The level table, only if it was measured with this very wasm (a stale table would scale wrongly).
async function loadLevels(wasmBytes) {
  try {
    const r = await fetch(new URL('./engine_levels.json', import.meta.url));
    if (!r.ok) return null;
    const table = await r.json();
    if (!globalThis.crypto || !crypto.subtle) return table; // plain-http dev server: cannot check
    const h = new Uint8Array(await crypto.subtle.digest('SHA-256', await wasmBytes));
    return Array.from(h, (b) => b.toString(16).padStart(2, '0')).join('') === table.wasmSha256 ? table : null;
  } catch {
    return null;
  }
}

function choose(name, scales, minRealtime) {
  let fallback = null;
  for (const s of scales) {
    try {
      const rt = bench(name, s);
      if (rt >= minRealtime) return { scale: s, realtime: rt };
      fallback = { scale: s, realtime: rt };
    } catch {
      // Scale requires oversampling at 48 kHz; try the next coarser grid.
    }
  }
  return fallback || { scale: 5, realtime: bench(name, 5) };
}

self.onmessage = async (event) => {
  const m = event.data || {};
  try {
    if (m.type === 'init') {
      const wasm = await fetch(new URL('./engine_bg.wasm', import.meta.url));
      const levels = loadLevels(wasm.clone().arrayBuffer());
      await init({ module_or_path: wasm });
      const table = await levels;
      const pageGain = Number.isFinite(m.gain) ? m.gain : gain;
      if (Number.isFinite(m.fieldHz)) fieldEvery = m.fieldHz > 0 ? 1000 / m.fieldHz : Infinity;
      if (Number(m.fieldFrames) > 0) { fieldFrames = Math.max(128, Math.round(Number(m.fieldFrames))); fieldEvery = Infinity; }
      if (m.bursts === false) bursts = false;
      const pick = choose(m.piece, m.scales || [2.5, 3, 3.5, 4, 5, 6, 7, 8], m.minRealtime || 1.7);
      piece = new KitPiece(m.piece, pick.scale);
      fix = levelFix(table, m.piece, pick.scale, pageGain);
      gain = pageGain * fix;
      const engine = {
        idle: () => piece.asleep(),
        apply: (e) => {
          piece.strike(e.x, e.y, e.speed);
          if (bursts && (fieldEvery !== Infinity || fieldFrames > 0)) burst = { frames: [], since: 0, frame: lane ? lane.frame : 0 };
          if (fieldFrames > 0) fieldDue = 0; // show the stroke's first moment
        },
        render: (n) => {
          let out;
          if (burst) {
            // step in short pieces so the first milliseconds can be watched in slow motion
            out = new Float32Array(n);
            let at = 0;
            while (at < n) {
              const k = Math.min(BURST_STEP - burst.since, n - at);
              out.set(piece.render(k, gain), at);
              at += k;
              burst.since += k;
              if (burst.since === BURST_STEP) {
                burst.since = 0;
                burst.frames.push(piece.velocity());
                if (burst.frames.length >= BURST_N) {
                  const all = new Float32Array(burst.frames.length * burst.frames[0].length);
                  burst.frames.forEach((f, i) => all.set(f, i * f.length));
                  post({ type: 'burst', frame: burst.frame, stepMs: (BURST_STEP * 1000) / 48000, count: burst.frames.length, v: all }, [all.buffer]);
                  burst = null;
                  out.set(piece.render(n - at, gain), at);
                  break;
                }
              }
            }
          } else out = piece.render(n, gain);
          forceMax = Math.max(forceMax, piece.force());
          const now = performance.now();
          let fieldNow = false;
          if (fieldFrames > 0) {
            fieldDue -= n;
            if (fieldDue <= 0) { fieldDue += fieldFrames; if (fieldDue <= 0) fieldDue = fieldFrames; fieldNow = !piece.asleep(); }
          } else if (!piece.asleep() && now - lastField >= fieldEvery) {
            lastField = now;
            fieldNow = true;
          }
          if (fieldNow) {
            const v = piece.velocity();
            post({ type: 'field', frame: lane ? lane.frame + n : 0, epoch: lane ? lane.epoch : 0, v, kinetic: piece.kinetic(), force: forceMax }, [v.buffer]);
            forceMax = 0;
          }
          if (now - lastStats >= 500 && lane) {
            lastStats = now;
            const s = lane.stats;
            post({ type: 'stats', realtime: s.wallMs > 0 ? s.audioMs / s.wallMs : 0, renderedFrames: s.renderedFrames, skippedFrames: s.skippedFrames, lateEvents: s.lateEvents, asleep: piece.asleep() });
          }
          return out;
        },
      };
      lane = createLane(m.port, engine, { aheadFrames: m.aheadFrames || 1024, chunk: m.chunk });
      const xy = piece.surfaceXY();
      post({ type: 'ready', info: JSON.parse(piece.info()), xy, h: piece.surfaceH(), radius: piece.radius(), scale: pick.scale, realtime: pick.realtime, levelFix: fix, levels: table ? 'ok' : 'none' }, [xy.buffer]);
    } else if (!piece || !lane) {
      return;
    } else if (m.type === 'strike') {
      lane.event({ frame: Number.isFinite(m.frame) ? m.frame : undefined, x: m.x, y: m.y, speed: m.speed });
    } else if (m.type === 'sequence') {
      lane.sequence(m.startFrame, m.loopFrames, m.events);
    } else if (m.type === 'score') {
      lane.score(Array.isArray(m.events) ? m.events : []);
    } else if (m.type === 'stopScore') {
      lane.stopScore();
    } else if (m.type === 'resync') {
      lane.resync(m.epoch);
    } else if (m.type === 'stop') {
      lane.stopSequence();
    } else if (m.type === 'gain') {
      gain = (Number(m.gain) || 0) * fix;
    }
  } catch (e) {
    post({ type: 'error', message: e && e.message ? String(e.message) : String(e) });
  }
};
