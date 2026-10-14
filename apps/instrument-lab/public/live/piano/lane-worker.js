// Live piano on a lane of fs-live-mix (ES module Worker): the same engine as worker.js, but
// every block is stamped with its context frame and key events can be scheduled on exact
// frames, so the piano can play in time with live drum pieces computed in other Workers.
//
// Main thread -> worker
//   {type:'init', port, options?, aheadFrames?, gain?}
//       gain     one fixed factor on every sample (level matching against another rendering; default 1)
//       options  {boardCoordinates?, courseModeBudget?, maxDynamicModes?, fullScalePa?}  (0 = engine default)
//                without boardCoordinates the worker measures itself: the default, then 24, then 16
//       minRealtime  the speed the measurement asks for (default 1.5)
//   {type:'note', midi, speed, frame?} | {type:'off', midi, frame?} | {type:'sustain', value, frame?}
//   {type:'sequence', startFrame, loopFrames, events:[{offset, on?, off?, speed?}]}
//   {type:'stop'}   end the sequence and release every key
// Worker -> main thread
//   {type:'ready', info, boardCoordinates, realtime} | {type:'stats', realtime, renderedFrames, skippedFrames, lateEvents, voices} | {type:'error', message}
import init, { PianoEngine } from './engine.js';
import { createLane } from '../lane.js';

let engine = null;
let lane = null;
let lastStats = 0;
let gain = 1;
const post = (m) => self.postMessage(m);

self.onmessage = async (event) => {
  const m = event.data || {};
  try {
    if (m.type === 'init') {
      await init({ module_or_path: new URL('./engine_bg.wasm', import.meta.url) });
      const o = m.options || {};
      gain = Number.isFinite(m.gain) ? m.gain : 1;
      const make = (coords) => new PianoEngine(48000, coords, Number(o.courseModeBudget) || 0, Number(o.maxDynamicModes) || 0, Number(o.fullScalePa) || 0);
      // a five-note stab with its bass: keep the most soundboard coordinates this device can
      // compute at least `minRealtime` times faster than real time
      const tries = o.boardCoordinates ? [Number(o.boardCoordinates)] : [0, 24, 16];
      let realtime = 0;
      for (const coords of tries) {
        engine?.free();
        engine = make(coords);
        for (const k of [38, 57, 61, 64, 66, 69]) engine.noteOn(k, 3.5);
        engine.render(512);
        const t0 = performance.now();
        for (let i = 0; i < 48; i++) engine.render(128);
        realtime = (48 * 128 * 1000) / 48000 / (performance.now() - t0);
        engine.free();
        engine = make(coords);
        if (realtime >= (m.minRealtime || 1.5)) break;
      }
      const adapter = {
        idle: () => false, // the soundboard and every sympathetic string keep moving
        apply: (e) => {
          if (e.on !== undefined) engine.noteOn(e.on, e.speed);
          else if (e.off !== undefined) engine.noteOff(e.off);
          else if (e.sustain !== undefined) engine.sustain(e.sustain);
          else if (e.allOff) engine.allNotesOff();
        },
        render: (n) => {
          const out = engine.render(n);
          if (gain !== 1) for (let i = 0; i < out.length; i++) out[i] *= gain;
          const now = performance.now();
          if (now - lastStats >= 500 && lane) {
            lastStats = now;
            const s = lane.stats;
            const st = engine.stats();
            post({ type: 'stats', realtime: s.wallMs > 0 ? s.audioMs / s.wallMs : 0, renderedFrames: s.renderedFrames, skippedFrames: s.skippedFrames, lateEvents: s.lateEvents, voices: st[1], peakPa: st[4] });
          }
          return out;
        },
      };
      lane = createLane(m.port, adapter, { aheadFrames: m.aheadFrames || 2048 });
      const info = JSON.parse(engine.info());
      post({ type: 'ready', info, boardCoordinates: info.boardCoordinates, realtime });
    } else if (!engine || !lane) {
      return;
    } else if (m.type === 'note') {
      lane.event({ frame: Number.isFinite(m.frame) ? m.frame : undefined, on: m.midi, speed: m.speed });
    } else if (m.type === 'off') {
      lane.event({ frame: Number.isFinite(m.frame) ? m.frame : undefined, off: m.midi });
    } else if (m.type === 'sustain') {
      lane.event({ frame: Number.isFinite(m.frame) ? m.frame : undefined, sustain: m.value });
    } else if (m.type === 'sequence') {
      lane.sequence(m.startFrame, m.loopFrames, m.events);
    } else if (m.type === 'stop') {
      lane.stopSequence();
      lane.event({ allOff: true });
    }
  } catch (e) {
    post({ type: 'error', message: e && e.message ? String(e.message) : String(e) });
  }
};
