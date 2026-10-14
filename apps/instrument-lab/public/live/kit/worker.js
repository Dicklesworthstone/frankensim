// Live drum-kit piece worker (ES module Worker): owns one fs-vkplate piece (engine.js) and
// feeds one lane of fs-live-mix. Samples are microphone pressure [Pa] times a fixed gain.
//
// Main thread -> worker
//   {type:'init', piece, port, gain?, aheadFrames?, scales?, minRealtime?, fieldHz?}   fieldHz 0: no fields
//       scales       grid coarsening factors to try, finest first (default [4, 5, 6, 7, 8])
//       minRealtime  the finest scale that runs at least this much faster than real time wins (1.7)
//   {type:'strike', x, y, speed, frame?}   x, y metres from the centre; speed m/s
//   {type:'sequence', startFrame, loopFrames, events:[{offset, x, y, speed}]}
//   {type:'stop'}                           end the sequence (ringing continues)
//   {type:'gain', gain}
// Worker -> main thread
//   {type:'ready', info, xy:Float32Array, h, radius, scale, realtime}
//   {type:'burst', frame, stepMs, count, v:Float32Array}    the first milliseconds after a stroke: `count`
//                                                          surface velocities `stepMs` apart, back to back
//   {type:'field', frame, v:Float32Array, kinetic, force}   ~fieldHz while the piece moves; force: the
//                                                          largest stick force since the last field
//   {type:'stats', realtime, renderedFrames, skippedFrames, lateEvents, asleep}
//   {type:'error', message}
import init, { KitPiece } from './engine.js';
import { createLane } from '../lane.js';

let piece = null;
let lane = null;
let gain = 0.03;
let fieldEvery = 1000 / 30;
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

function choose(name, scales, minRealtime) {
  let fallback = null;
  for (const s of scales) {
    const rt = bench(name, s);
    if (rt >= minRealtime) return { scale: s, realtime: rt };
    fallback = { scale: s, realtime: rt };
  }
  return fallback;
}

self.onmessage = async (event) => {
  const m = event.data || {};
  try {
    if (m.type === 'init') {
      await init({ module_or_path: new URL('./engine_bg.wasm', import.meta.url) });
      gain = Number.isFinite(m.gain) ? m.gain : gain;
      if (Number.isFinite(m.fieldHz)) fieldEvery = m.fieldHz > 0 ? 1000 / m.fieldHz : Infinity;
      const pick = choose(m.piece, m.scales || [4, 5, 6, 7, 8], m.minRealtime || 1.7);
      piece = new KitPiece(m.piece, pick.scale);
      const engine = {
        idle: () => piece.asleep(),
        apply: (e) => { piece.strike(e.x, e.y, e.speed); if (fieldEvery !== Infinity) burst = { frames: [], since: 0, frame: lane ? lane.frame : 0 }; },
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
          if (!piece.asleep() && now - lastField >= fieldEvery) {
            lastField = now;
            const v = piece.velocity();
            post({ type: 'field', frame: lane ? lane.frame + n : 0, v, kinetic: piece.kinetic(), force: forceMax }, [v.buffer]);
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
      lane = createLane(m.port, engine, { aheadFrames: m.aheadFrames || 1024 });
      const xy = piece.surfaceXY();
      post({ type: 'ready', info: JSON.parse(piece.info()), xy, h: piece.surfaceH(), radius: piece.radius(), scale: pick.scale, realtime: pick.realtime }, [xy.buffer]);
    } else if (!piece || !lane) {
      return;
    } else if (m.type === 'strike') {
      lane.event({ frame: Number.isFinite(m.frame) ? m.frame : undefined, x: m.x, y: m.y, speed: m.speed });
    } else if (m.type === 'sequence') {
      lane.sequence(m.startFrame, m.loopFrames, m.events);
    } else if (m.type === 'stop') {
      lane.stopSequence();
    } else if (m.type === 'gain') {
      gain = Number(m.gain) || 0;
    }
  } catch (e) {
    post({ type: 'error', message: e && e.message ? String(e.message) : String(e) });
  }
};
