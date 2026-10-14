// Live piano on a lane of fs-live-mix (ES module Worker): the same engine as worker.js, but
// every block is stamped with its context frame and key events can be scheduled on exact
// frames, so the piano can play in time with live drum pieces computed in other Workers.
//
// Main thread -> worker
//   {type:'init', port, options?, aheadFrames?, gain?, minRealtime?, snapshotHz?}
//       gain     one fixed factor on every sample (level matching against another rendering; default 1)
//       options  {boardCoordinates?, courseModeBudget?, maxDynamicModes?, fullScalePa?}  (0 = engine default)
//                boardCoordinates: 32 = the embedded board, 48/64/96/128 = that board file from
//                ./boards/, 1..16 = that many exact low modes. Without it the worker measures
//                itself and keeps the LARGEST board it computes `minRealtime` times faster than
//                real time (boards.js): 32 first, then at most two larger ones, else 16.
//       minRealtime  the speed the measurement asks for (default 1.35)
//       snapshotHz   display snapshots per second of AUDIO (default 0: none). With it, 'ready' also
//                    carries the board geometry and the key geometry.
//   {type:'note', midi, speed, frame?} | {type:'off', midi, frame?} | {type:'sustain', value, frame?}
//   {type:'sequence', startFrame, loopFrames, events:[{offset, on?, off?, speed?}]}
//   {type:'score', events:[{frame, on?, off?, speed?, sustain?}]}   a fixed score, played once
//   {type:'stopScore'}   drop the rest of the score and release every key and the pedal
//   {type:'resync', epoch}   render for the clock as it is now (after the mix flushed this lane)
//   {type:'stop'}   end the sequence and release every key
// Worker -> main thread
//   {type:'ready', info, boardCoordinates, realtime, board?:{xy, shapes, modeHz}, keyGeometry?}
//   {type:'stats', realtime, renderedFrames, skippedFrames, lateEvents, voices, peakPa,
//    dynamicModes, noteOns, refusedNoteOns, forcedRetirements}
//   {type:'snapshot', frame, epoch, data}   engine.snapshot() for the state at context frame `frame`
//                                    (layout: parseSnapshot in client.js); `epoch` as on the audio
//   {type:'error', message}
import init, { PianoEngine } from './engine.js';
import { createLane } from '../lane.js';
import { chooseBoard, fetchBoard, probeRealtime } from './boards.js';

let engine = null;
let lane = null;
let lastStats = 0;
let gain = 1;
let snapEvery = 0;   // audio frames between snapshots (0: none)
let snapDue = 0;
const post = (m, t) => self.postMessage(m, t ?? []);
const FALLBACK_REALTIME = 0.75;

self.onmessage = async (event) => {
  const m = event.data || {};
  try {
    if (m.type === 'init') {
      await init({ module_or_path: new URL('./engine_bg.wasm', import.meta.url) });
      const o = m.options || {};
      gain = Number.isFinite(m.gain) ? m.gain : 1;
      snapEvery = Number(m.snapshotHz) > 0 ? Math.max(256, Math.round(48000 / Math.min(120, Number(m.snapshotHz)))) : 0;
      const limits = [Number(o.courseModeBudget) || 0, Number(o.maxDynamicModes) || 0, Number(o.fullScalePa) || 0];
      // bytes: a fetched board file (PianoEngine.withBoard); otherwise the embedded board or its
      // exact low-mode truncation (coords 1..16; 0/32 = all 32 embedded coordinates)
      const make = (bytes, coords) => bytes ? PianoEngine.withBoard(48000, bytes, ...limits)
        : new PianoEngine(48000, coords > 16 ? 0 : coords, ...limits);
      const measure = probeRealtime; // twelve keys struck together, timed on a throwaway engine
      const target = m.minRealtime || 1.35;
      const asked = Number(o.boardCoordinates) || 0;
      let bytes = null;
      let coords = asked;
      let realtime = 0;
      if (asked > 32) {
        bytes = await fetchBoard(asked);
      } else if (!asked) {
        const pick = await chooseBoard({ make, measure, target });
        pick.engine.free();
        ({ bytes, coordinates: coords, realtime } = pick);
        // the old last resort, 16 exact low modes, only when even the embedded board cannot keep
        // up with the probe (it is about twice the old five-note stab's load)
        if (!bytes && realtime < FALLBACK_REALTIME) { coords = 16; realtime = 0; }
      }
      if (!realtime) {
        const probe = make(bytes, coords);
        realtime = measure(probe);
        probe.free();
      }
      engine = make(bytes, coords);
      let active = false;
      let quietFrames = 0;
      const adapter = {
        idle: () => !active, // at rest before the first strike or once the soundboard has died away
        apply: (e) => {
          if (e.on !== undefined) { active = true; quietFrames = 0; engine.noteOn(e.on, e.speed); }
          else if (e.off !== undefined) engine.noteOff(e.off);
          else if (e.sustain !== undefined) engine.sustain(e.sustain);
          else if (e.allOff) engine.allNotesOff();
        },
        render: (n) => {
          const out = engine.render(n);
          let peak = 0;
          for (let i = 0; i < out.length; i++) {
            const v = out[i] * gain;
            out[i] = v;
            const a = Math.abs(v);
            if (a > peak) peak = a;
          }
          if (peak < 1e-4) {
            quietFrames += out.length;
            if (quietFrames >= 14400 && engine.stats()[1] === 0) active = false;
          } else {
            quietFrames = 0;
          }
          if (snapEvery > 0 && lane) {
            snapDue -= n;
            if (snapDue <= 0) {
              snapDue += snapEvery;
              if (snapDue <= 0) snapDue = snapEvery;
              const data = engine.snapshot();
              post({ type: 'snapshot', frame: lane.frame + n, epoch: lane.epoch, data }, [data.buffer]);
            }
          }
          const now = performance.now();
          if (now - lastStats >= 500 && lane) {
            lastStats = now;
            const s = lane.stats;
            const st = engine.stats();
            post({ type: 'stats', realtime: s.wallMs > 0 ? s.audioMs / s.wallMs : 0, renderedFrames: s.renderedFrames, skippedFrames: s.skippedFrames, lateEvents: s.lateEvents, voices: st[1], peakPa: st[4], dynamicModes: st[2], noteOns: st[10], refusedNoteOns: st[11], forcedRetirements: st[6] });
          }
          return out;
        },
      };
      lane = createLane(m.port, adapter, { aheadFrames: m.aheadFrames || 2048 });
      const info = JSON.parse(engine.info());
      if (snapEvery > 0) {
        const xy = engine.boardSurfaceXY();
        const shapes = engine.boardSurfaceShapes();
        const modeHz = engine.boardModeHz();
        const keyGeometry = engine.keyGeometry();
        post({ type: 'ready', info, boardCoordinates: info.boardCoordinates, realtime, board: { xy, shapes, modeHz }, keyGeometry },
          [xy.buffer, shapes.buffer, modeHz.buffer, keyGeometry.buffer]);
      } else {
        post({ type: 'ready', info, boardCoordinates: info.boardCoordinates, realtime });
      }
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
    } else if (m.type === 'score') {
      lane.score(Array.isArray(m.events) ? m.events : []);
    } else if (m.type === 'stopScore') {
      lane.stopScore();
      lane.event({ allOff: true });
    } else if (m.type === 'resync') {
      lane.resync(m.epoch);
    } else if (m.type === 'stop') {
      lane.stopSequence();
      lane.event({ allOff: true });
    }
  } catch (e) {
    post({ type: 'error', message: e && e.message ? String(e.message) : String(e) });
  }
};
