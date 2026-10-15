// Live piano on a lane of fs-live-mix (ES module Worker): the same engine as worker.js, but
// every block is stamped with its context frame and key events can be scheduled on exact
// frames, so the piano can play in time with live drum pieces computed in other Workers.
//
// Main thread -> worker
//   {type:'init', port, options?, aheadFrames?, chunk?, gain?, minRealtime?, snapshotHz?}
//       chunk    most frames per audio message (default 128; lane.js)
//       gain     one fixed factor on every sample (level matching against another rendering; default 1)
//       options  {boardCoordinates?, courseModeBudget?, maxDynamicModes?, fullScalePa?}  (0 = engine default)
//                boardCoordinates: 32 = the embedded board, 48/64/96/128 = that board file from
//                ./boards/, 1..16 = that many exact low modes. Without it the worker times a
//                sustained workload (boards.js probeSustained) and picks board AND starting mode
//                cap together (chooseQuality): the largest board that runs `minRealtime` times
//                faster than real time at a cap no lower than its floor (boards.js MIN_CAP; 32
//                first, then the best board predicted to keep its floor, timed at the predicted
//                cap; at most two larger ones), else the embedded board with the cap lowered
//                until it does; 16 only when even the floor cap is expected below real time.
//                A fixed board gets its cap
//                the same way; an explicit maxDynamicModes is kept.
//       minRealtime  the speed the measurement asks for (default 1.3)
//       snapshotHz   display snapshots per second of AUDIO (default 0: none). With it, 'ready' also
//                    carries the board geometry and the key geometry.
//       governor     false turns the load governor off (default on; see createGovernor)
//   {type:'note', midi, speed, frame?} | {type:'off', midi, frame?} | {type:'sustain', value, frame?}
//   {type:'sequence', startFrame, loopFrames, events:[{offset, on?, off?, speed?}]}
//   {type:'score', events:[{frame, on?, off?, speed?, sustain?}]}   a fixed score, played once
//   {type:'stopScore'}   drop the rest of the score and release every key and the pedal
//   {type:'resync', epoch}   render for the clock as it is now (after the mix flushed this lane)
//   {type:'stop'}   end the sequence and release every key
// Worker -> main thread
//   {type:'ready', info, boardCoordinates, realtime, modeCap, board?:{xy, shapes, modeHz}, keyGeometry?}
//       realtime  the probe's speed at the chosen board and cap; modeCap  the starting mode cap
//   {type:'stats', realtime, renderedFrames, skippedFrames, lateEvents, voices, peakPa,
//    dynamicModes, noteOns, refusedNoteOns, forcedRetirements,
//    modeCap, fullModeCap, courseModeBudget, governorDown, governorUp}
//       modeCap  the dynamic string mode cap the governor grants now (fullModeCap = the asked one)
//   {type:'snapshot', frame, epoch, data}   engine.snapshot() for the state at context frame `frame`
//                                    (layout: parseSnapshot in client.js); `epoch` as on the audio
//   {type:'error', message}
import init, { PianoEngine } from './engine.js';
import { createLane } from '../lane.js';
import { BOARD_SIZES, chooseQuality, fetchBoard, probeSustained } from './boards.js';

let engine = null;
let lane = null;
let lastStats = 0;
let gain = 1;
let snapEvery = 0;   // audio frames between snapshots (0: none)
let snapDue = 0;
const post = (m, t) => self.postMessage(m, t ?? []);

// Load governor. The piano's cost follows the dynamic string modes sounding at once (about
// linearly), and a dense passage can need more than this device computes in real time. When the
// look-ahead runs down while the measured speed is below real time (or the lane had to skip), the
// mode cap for strikes from then on is lowered to the modes sounding now times the measured speed
// over GOV.target, so the quietest keys yield at the next strikes. The per-key mode budget follows
// the cap in proportion, so about as many keys keep sounding and each loses upper partials instead:
// against the full-quality render that changed the spectrogram about a third as much as retiring
// quiet keys at the same cap, and the partial band (it dulls every key) is never touched. The cap
// returns 10 % at a time only once the look-ahead is full again and the speed measured at the cap
// leaves GOV.raiseHeadroom of room at the raised cap, so a passage settles instead of see-sawing.
// A device that keeps up never meets the lowering condition, so its piano is unchanged.
const GOV = {
  everyMs: 400, minFrames: 4096,      // judge windows of at least this much wall time and audio
  lowCushion: 0.85, slow: 1.0,        // lower after two windows with look-ahead < 85 % and speed < real time ...
  skipSlow: 1.1,                      // ... or at once when the lane skipped and speed < 1.1
  target: 1.25, minFactor: 0.75, maxFactor: 0.92, coolMs: 800,
  fullCushion: 0.85, holdMs: 3000, raise: 1.1, raiseHeadroom: 1.35,
  floorModes: 400, floorBudget: 40,
};
// per-key string mode budget for a mode cap: in proportion to the cap
// The startup choice never starts below this cap: offline (V8, board 32) cap 700 ran Giant Steps and
// Character Select 2.65-2.86x faster than real time and cap 450 2.76x, while 450 moved the band levels
// 1.67-1.76 dB from the full cap against 0.59-0.86 dB at 700. Below it only the governor goes.
const START_FLOOR_MODES = 700;
const budgetFor = (full, cap) => Math.max(GOV.floorBudget, Math.min(full.courseModeBudget, Math.round((full.courseModeBudget * cap) / full.maxDynamicModes)));
function createGovernor(eng, full, aheadFrames, startCap = full.maxDynamicModes) {
  const g = { cap: startCap, budget: budgetFor(full, startCap), down: 0, up: 0, at: -1, changedAt: -1e9, audioMs: 0, wallMs: 0, skipped: 0, rendered: 0, cushion: 0, slowWindows: 0, modes: 0 };
  const set = (cap) => {
    const budget = budgetFor(full, cap);
    if (!eng.setQuality(budget, full.partialBandHz, cap)) return false;
    g.cap = cap;
    g.budget = budget;
    return true;
  };
  // called after each block: `frames` is how far ahead of the clock the lane has rendered
  g.check = (stats, frames) => {
    const t = performance.now();
    if (g.at >= 0 && t - g.at < GOV.everyMs) return;
    const audioMs = stats.audioMs - g.audioMs;
    const skipped = stats.skippedFrames - g.skipped;
    if (g.at >= 0 && audioMs * 48 < GOV.minFrames && !skipped) return;
    const wallMs = stats.wallMs - g.wallMs;
    const rendered = stats.renderedFrames - g.rendered;
    const span = t - g.at;
    const cushionBefore = g.cushion;
    const fresh = g.at < 0 || span > 4 * GOV.everyMs; // first window, or rendering paused in between
    const modesNow = eng.stats()[2];
    const modes = fresh ? modesNow : Math.max(modesNow, g.modes); // sounding modes, at either end of the window
    Object.assign(g, { at: t, audioMs: stats.audioMs, wallMs: stats.wallMs, skipped: stats.skippedFrames, rendered: stats.renderedFrames, cushion: frames, modes: modesNow });
    if (fresh || !(wallMs > 0)) { g.slowWindows = 0; return; }
    let speed = audioMs / wallMs; // speed of the engine calls themselves
    // behind at both ends of the window, the lane rendered without waiting: frames produced per
    // wall frame is the true speed, including time this Worker was not scheduled at all
    if (frames < 0.9 * aheadFrames && cushionBefore < 0.9 * aheadFrames) speed = Math.min(speed, rendered / (span * 48));
    const cushion = frames / aheadFrames;
    g.slowWindows = cushion < GOV.lowCushion && speed < GOV.slow ? g.slowWindows + 1 : 0;
    if (g.slowWindows >= 2 || (skipped > 0 && speed < GOV.skipSlow)) {
      if (t - g.changedAt < GOV.coolMs || g.cap <= GOV.floorModes) return;
      const factor = Math.min(GOV.maxFactor, Math.max(GOV.minFactor, speed / GOV.target));
      const cap = Math.max(GOV.floorModes, Math.min(g.cap - 1, Math.round(Math.min(g.cap, modes) * factor)));
      if (set(cap)) { g.down += 1; g.changedAt = t; g.slowWindows = 0; }
    } else if (g.cap < full.maxDynamicModes && cushion >= GOV.fullCushion && !skipped && t - g.changedAt >= GOV.holdMs) {
      const cap = Math.min(full.maxDynamicModes, Math.ceil(g.cap * GOV.raise));
      if ((speed * modes) / cap >= GOV.raiseHeadroom && set(cap)) { g.up += 1; g.changedAt = t; }
    }
  };
  return g;
}

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
      // exact low-mode truncation (coords 1..16; 0/32 = all 32 embedded coordinates). cap: a lower
      // dynamic mode cap, the per-key budget in proportion (as the governor sets it)
      let full = null; // the quality asked for (read from the first engine made)
      // a new engine's limits for a lower cap (the same call on an unplayed copy gives it a new engine's limits)
      const atCap = (e, cap) => {
        if (cap && cap < full.maxDynamicModes) e.setQuality(budgetFor(full, cap), full.partialBandHz, cap);
        return e;
      };
      const make = (bytes, coords, cap) => {
        const e = bytes ? PianoEngine.withBoard(48000, bytes, ...limits) : new PianoEngine(48000, coords > 16 ? 0 : coords, ...limits);
        if (!full) {
          const i = JSON.parse(e.info());
          full = { courseModeBudget: i.courseModeBudget, partialBandHz: i.partialBandHz, maxDynamicModes: i.maxDynamicModes };
        }
        return atCap(e, cap);
      };
      // Board and starting mode cap are chosen together by timing a sustained, damped, cap-bound
      // workload (boards.js), for `minRealtime` headroom, so the governor need not catch the first
      // bars. A fixed board keeps its board and still gets its cap; an explicit maxDynamicModes is
      // kept as it is.
      const target = m.minRealtime || 1.3;
      const asked = Number(o.boardCoordinates) || 0;
      let bytes = asked > 32 ? await fetchBoard(asked) : null;
      let coords = asked;
      const pick = await chooseQuality({ make: asked ? (_b, _c, k) => make(bytes, coords, k) : make, target, sizes: asked ? [] : BOARD_SIZES,
        floor: Number(o.maxDynamicModes) > 0 ? Infinity : START_FLOOR_MODES });
      if (!asked) ({ bytes, coordinates: coords } = pick);
      let { cap, realtime } = pick;
      // the old last resort, 16 exact low modes, only when even the embedded board at the floor cap
      // is expected below real time
      if (!asked && !bytes && cap <= START_FLOOR_MODES && realtime < 1) {
        pick.engine?.free();
        pick.engine = null;
        coords = 16;
        cap = full.maxDynamicModes;
        const probe = make(null, 16);
        realtime = probeSustained(probe);
        probe.free();
      }
      // the chosen board was built once, for its timing: its copy from before the timing plays
      // (bit for bit what a new engine at this cap renders)
      engine = pick.engine ? atCap(pick.engine, cap) : make(bytes, coords, cap);
      const info = JSON.parse(engine.info());
      const aheadFrames = m.aheadFrames || 2048;
      const gov = m.governor === false ? null : createGovernor(engine, full, aheadFrames, cap);
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
          if (gov && lane) gov.check(lane.stats, lane.frame + n - lane.now());
          const now = performance.now();
          if (now - lastStats >= 500 && lane) {
            lastStats = now;
            const s = lane.stats;
            const st = engine.stats();
            post({ type: 'stats', realtime: s.wallMs > 0 ? s.audioMs / s.wallMs : 0, renderedFrames: s.renderedFrames, skippedFrames: s.skippedFrames, lateEvents: s.lateEvents, voices: st[1], peakPa: st[4], dynamicModes: st[2], noteOns: st[10], refusedNoteOns: st[11], forcedRetirements: st[6],
              modeCap: gov ? gov.cap : full.maxDynamicModes, fullModeCap: full.maxDynamicModes, courseModeBudget: gov ? gov.budget : full.courseModeBudget, governorDown: gov ? gov.down : 0, governorUp: gov ? gov.up : 0 });
          }
          return out;
        },
      };
      lane = createLane(m.port, adapter, { aheadFrames, chunk: m.chunk });
      if (snapEvery > 0) {
        const xy = engine.boardSurfaceXY();
        const shapes = engine.boardSurfaceShapes();
        const modeHz = engine.boardModeHz();
        const keyGeometry = engine.keyGeometry();
        post({ type: 'ready', info, boardCoordinates: info.boardCoordinates, realtime, modeCap: cap, board: { xy, shapes, modeHz }, keyGeometry },
          [xy.buffer, shapes.buffer, modeHz.buffer, keyGeometry.buffer]);
      } else {
        post({ type: 'ready', info, boardCoordinates: info.boardCoordinates, realtime, modeCap: cap });
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
