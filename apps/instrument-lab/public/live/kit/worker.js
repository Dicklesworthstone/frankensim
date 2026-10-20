// Live drum-kit piece worker (ES module Worker): owns one fs-vkplate piece (engine.js) and
// feeds one lane of fs-live-mix. Samples are microphone pressure [Pa] times a fixed gain, which the
// page sets for the 5x grid; on the grid picked here it is scaled so that a stroke is as loud as on
// 5x (levels.js, engine_levels.json beside the wasm).
//
// Main thread -> worker
//   {type:'init', piece, port, gasRecord?, gain?, aheadFrames?, chunk?, scales?, minRealtime?, benchMs?, fieldHz?, fieldFrames?, bursts?}
//       gasRecord    shared material record: `T_K P_Pa RH`, or
//                    `mixture T_K P_Pa mole|mass COUNT SPECIES FRACTION ...`
//                    Held for every grid trial, clone and fallback construction.
//       chunk        most frames per audio message (default 128; lane.js)
//       fieldHz      fields per second of wall time (0: no fields)
//       fieldFrames  instead: one field every this many AUDIO frames (even spacing however far
//                    ahead the lane renders; for a score heard a second after it is computed)
//       bursts       false: no slow-motion burst after each stroke (default true when fields are on)
//       scales       grid coarsening factors to try, finest first (default [2.5, 3, 3.5, 4, 5, 6, 7, 8])
//       minRealtime  the finest scale that runs at least this much faster than real time wins (1.7);
//                    timed from the coarse end (choose)
//       benchMs      each grid is timed over this much of one stroke (default 112 ms). A piece that costs
//                    more once it rings than at its attack needs longer: the snare's wires touch the head
//                    from ~85 ms after a stroke, and its steps then cost ~1.6x as much
//   {type:'strike', x, y, speed, frame?}   x, y metres from the centre; speed m/s
//   {type:'sequence', startFrame, loopFrames, events:[{offset, x, y, speed}]}
//   {type:'score', events:[{frame, x, y, speed}]}   a fixed score, played once
//   {type:'stopScore'}                      drop the rest of the score (ringing continues)
//   {type:'resync', epoch}                  render for the clock as it is now (after a mix flush)
//   {type:'stop'}                           end the sequence (ringing continues)
//   {type:'gain', gain}                     a new page gain (scaled by the same levelFix)
// Worker -> main thread
//   {type:'ready', info, gasRecord?, xy:Float32Array, h, radius, scale, realtime, levelFix, levels}
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

// How fast the grid of `p` (a piece as `new` left it) runs: one stroke, then 128-frame blocks timed from 256
// frames after it (the first touch of every buffer) until `benchMs` after it (at least 40 blocks: the default
// 112 ms). The timing stops as soon as the window can no longer reach `minRealtime`, since the rest could only
// lower its speed. `p` is freed.
function bench(p, benchMs, minRealtime) {
  try {
    p.strike(0.3 * p.radius(), 0, 3);
    p.render(256, 1);
    const blocks = Math.max(40, Math.round((benchMs * 48 - 256) / 128));
    const most = (blocks * 128) / 48 / minRealtime; // wall ms a passing window may take
    const t0 = performance.now();
    let n = 0, ms = 0;
    while (n < blocks && ms <= most) {
      p.render(128, 1);
      n++;
      ms = performance.now() - t0;
    }
    return (n * 128) / 48 / ms;
  } finally {
    p.free();
  }
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

// A copy of `p` in its present state, or null when this engine cannot copy. A build serves public/live under
// live/<content hash>/ (vite.config.ts), so the glue and the wasm come from one build; this guard stays for a
// pairing that can still differ (a dev server mid-rebuild, a host that rewrites caching): a wasm from before
// copies has no kitpiece_copy, so copy() throws, and an engine.js from before copies has no copy() at all. The
// chosen grid is then built again to play, as it was before copies: the same samples, one more settle.
function copyOf(p) {
  try {
    return p.copy();
  } catch {
    return null;
  }
}

// The material owner admits the record; the browser never substitutes another
// gas after a refused constructor or when a deployed WASM predates withGas.
function makePiece(name, scale, gasRecord) {
  if (gasRecord === undefined) return new KitPiece(name, scale);
  if (typeof gasRecord !== 'string' || !gasRecord.trim()) {
    throw new Error('A nonempty atmosphere record is required.');
  }
  if (typeof KitPiece.withGas !== 'function') {
    throw new Error('This kit engine needs rebuilding with atmosphere support; the selected gas was not applied.');
  }
  return KitPiece.withGas(name, scale, gasRecord);
}

// The finest of `scales` (listed finest first) that runs `minRealtime` times faster than real time; if none
// does, the coarsest. The grids are timed from the coarse end, stepping finer until one falls short. A finer
// grid has more nodes and never runs faster, so this is the grid a finest-first search picks. But no grid finer
// than the first one that falls short is built: building is the cost of start-up (the snare takes 0.1 s at 6x,
// 0.5 s at 2.5x).
// Each grid is built once: before its build is struck and timed, a copy is kept, and the chosen grid's copy is
// the piece that plays (`piece`; null when the engine cannot copy, see copyOf: the caller builds the grid again).
// A copy renders exactly what a fresh `new` would, and takes under a millisecond, while a snare's or hi-hat's
// build is almost all its settle (the wire preload, the pedal: 0.05-0.9 s). The timings run on builds, not
// copies: a copy runs 1-3 % faster (its memory is laid out in one piece), so it would read high.
// The first timing in a Worker runs while V8 is still optimizing the engine code, and can read slow. A
// shortfall there is timed once more, on a new build, before it ends the search.
let timed = 0;
function choose(name, scales, minRealtime, benchMs, gasRecord) {
  let pick = null;
  for (let i = scales.length - 1; i >= 0; i--) {
    let kept = null, rt;
    try {
      const built = makePiece(name, scales[i], gasRecord);
      kept = copyOf(built);
      rt = bench(built, benchMs, minRealtime);
      if (timed++ === 0 && rt < minRealtime) rt = bench(makePiece(name, scales[i], gasRecord), benchMs, minRealtime);
    } catch {
      kept?.free();
      continue; // refused (the grid needs oversampling at 48 kHz) or failed: skip it, as a finest-first search would
    }
    if (rt < minRealtime) {
      if (!pick) return { scale: scales[i], realtime: rt, piece: kept };
      kept?.free();
      return pick;
    }
    pick?.piece?.free();
    pick = { scale: scales[i], realtime: rt, piece: kept };
  }
  if (pick) return pick;
  const built = makePiece(name, 5, gasRecord);
  const kept = copyOf(built);
  try {
    return { scale: 5, realtime: bench(built, benchMs, minRealtime), piece: kept };
  } catch (e) {
    kept?.free();
    throw e;
  }
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
      const gasRecord = m.gasRecord;
      const pick = choose(m.piece, m.scales || [2.5, 3, 3.5, 4, 5, 6, 7, 8], m.minRealtime || 1.7, Number(m.benchMs) > 0 ? Number(m.benchMs) : 112, gasRecord);
      piece = pick.piece || makePiece(m.piece, pick.scale, gasRecord);
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
      post({ type: 'ready', info: JSON.parse(piece.info()), gasRecord, xy, h: piece.surfaceH(), radius: piece.radius(), scale: pick.scale, realtime: pick.realtime, levelFix: fix, levels: table ? 'ok' : 'none' }, [xy.buffer]);
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
