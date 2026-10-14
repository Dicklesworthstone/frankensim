// Live piano engine worker (ES module Worker).
//
// Owns the WebAssembly engine (FrankenSim's reduced real-time piano), renders
// ahead in 128-frame chunks to keep a short FIFO in the shared `fs-live-sink`
// AudioWorklet, and reports state. It never post-processes audio: a sample is
// the pressure [Pa] at the declared Rayleigh receiver divided by the declared
// fixed full-scale pressure. No gain riding, no limiter, no fade.
//
// Main thread -> worker
//   {type:'init', sampleRate, port, options?, wasmUrl?}
//       port     MessagePort adopted by the sink ({type:'connect', port} on the node's port)
//       options  {boardCoordinates?, courseModeBudget?, maxDynamicModes?, fullScalePa?,
//                 targetMs?, snapshotHz?, governor?}   (0/undefined = engine default; targetMs default 18)
//   {type:'cmd', name, args?, id?}
//       noteOn(midi, speedMS)   post-escapement hammer velocity in m/s, (0, 8]
//       noteOff(midi)
//       sustain(value)          pedal travel 0..1
//       allNotesOff()
//       setQuality({courseModeBudget?, partialBandHz?, maxDynamicModes?})
//       setTargetMs(ms) | setSnapshotHz(hz) | setGovernor(bool)
//       partials(midi)          -> {type:'reply', id, name, value: Float64Array}
//       bench({keys?, seconds?, speedMS?, sustain?})  offline speed test on a
//                               SEPARATE engine instance -> {type:'bench', id, results}
//       dispose()
// Worker -> main thread
//   {type:'ready', info, board:{xy, shapes, modeHz}, keyGeometry, initMs}
//       keyGeometry  Float64Array, per key: midi, length_m, strike_fraction, unison, f1_hz
//   {type:'snapshot', data:Float32Array, renderedFrames, aheadMs}   (layout: see parseSnapshot in client.js)
//   {type:'stats', realtimeFactor, load, worstChunkMs, queuedMs, targetMs, underruns, underrunFrames,
//    voices, dynamicModes, dynamicStrings, peakPa, clippedSamples, forcedRetirements, retirements,
//    contactFailures, refusedNoteOns, rebuilds, renderedFrames, playedFrames, quality, governorSteps}
//   {type:'reply', id, name, value} | {type:'bench', id, results}
//   {type:'error', code:'refused'|'engine'|'worker_error', message, id?, command?}
import init, { PianoEngine } from './engine.js';

const CHUNK = 128;

let engine = null;
let sinkPort = null;
let sampleRate = 48000;
let info = null;
let targetFrames = Math.round(0.018 * 48000);
let snapshotIntervalMs = 1000 / 30;
let timer = 0;
let sentTotal = 0;
let rendered = 0;
let level = { queuedFrames: 0, receivedFrames: 0, playedFrames: 0, underruns: 0, underrunFrames: 0, at: 0 };
let lastSnapshotAt = 0;
let lastStatsAt = 0;
let windowWallMs = 0;
let windowAudioMs = 0;
let worstChunkMs = 0;
let engineOptions = null;

// Quality the user asked for, and what the load governor currently grants.
let asked = { courseModeBudget: 0, partialBandHz: 0, maxDynamicModes: 0 };
let granted = { courseModeBudget: 0, partialBandHz: 0, maxDynamicModes: 0 };
let governor = true;
let governorSteps = 0;
let loadAverage = 0;
let lastGovernorAt = 0;

const post = (message, transfer) => self.postMessage(message, transfer ?? []);
// Error codes: 'refused' (a command the engine declined), 'engine' (a refusal
// raised inside the mechanics, e.g. a strike withdrawn), 'worker_error'.
const fail = (error, extra = {}) => post({ type: 'error', code: 'worker_error', message: error && error.message ? String(error.message) : String(error), ...extra });

function estimateQueue(now) {
  const inFlight = sentTotal - level.receivedFrames;
  const consumed = level.at > 0 ? ((now - level.at) * sampleRate) / 1000 : 0;
  return Math.max(0, level.queuedFrames + inFlight - consumed);
}

function renderChunk() {
  const t0 = performance.now();
  const samples = engine.render(CHUNK);
  const elapsed = performance.now() - t0;
  const chunkMs = (CHUNK * 1000) / sampleRate;
  windowWallMs += elapsed;
  windowAudioMs += chunkMs;
  if (elapsed > worstChunkMs) worstChunkMs = elapsed;
  loadAverage += 0.02 * (elapsed / chunkMs - loadAverage);
  sinkPort.postMessage({ type: 'audio', channels: [samples] }, [samples.buffer]);
  sentTotal += CHUNK;
  rendered += CHUNK;
}

// Load governor: when the smoothed render time exceeds 80 % of real time the
// DYNAMIC MODE BUDGET for subsequently struck keys is lowered (reported in
// stats.quality); it recovers below 45 %. It never touches the audio signal.
function govern(now) {
  if (!governor || now - lastGovernorAt < 400) return;
  lastGovernorAt = now;
  const floorModes = 400;
  if (loadAverage > 0.8 && granted.maxDynamicModes > floorModes) {
    granted.maxDynamicModes = Math.max(floorModes, Math.floor(granted.maxDynamicModes * 0.8));
    granted.courseModeBudget = Math.max(60, Math.floor(granted.courseModeBudget * 0.85));
    governorSteps += 1;
    engine.setQuality(granted.courseModeBudget, granted.partialBandHz, granted.maxDynamicModes);
  } else if (loadAverage < 0.45 && (granted.maxDynamicModes < asked.maxDynamicModes || granted.courseModeBudget < asked.courseModeBudget)) {
    granted.maxDynamicModes = Math.min(asked.maxDynamicModes, Math.ceil(granted.maxDynamicModes * 1.1));
    granted.courseModeBudget = Math.min(asked.courseModeBudget, Math.ceil(granted.courseModeBudget * 1.08));
    engine.setQuality(granted.courseModeBudget, granted.partialBandHz, granted.maxDynamicModes);
  }
}

function pump() {
  if (!engine || !sinkPort) return;
  let now = performance.now();
  let queue = estimateQueue(now);
  // Bounded work per turn, so commands are still served when the engine is
  // slower than real time (the sink then reports honest underruns).
  let budget = Math.ceil(targetFrames / CHUNK) + 2;
  while (queue < targetFrames && budget-- > 0) {
    renderChunk();
    queue += CHUNK;
  }
  now = performance.now();
  const error = engine.lastError();
  if (error) post({ type: 'error', code: 'engine', message: error });
  govern(now);
  if (snapshotIntervalMs > 0 && now - lastSnapshotAt >= snapshotIntervalMs) {
    lastSnapshotAt = now;
    const data = engine.snapshot();
    post({ type: 'snapshot', data, renderedFrames: rendered, aheadMs: (estimateQueue(now) * 1000) / sampleRate }, [data.buffer]);
  }
  if (now - lastStatsAt >= 250) {
    lastStatsAt = now;
    const s = engine.stats();
    post({
      type: 'stats',
      realtimeFactor: windowWallMs > 0 ? windowAudioMs / windowWallMs : null,
      load: loadAverage,
      worstChunkMs,
      queuedMs: (estimateQueue(now) * 1000) / sampleRate,
      targetMs: (targetFrames * 1000) / sampleRate,
      underruns: level.underruns,
      underrunFrames: level.underrunFrames,
      voices: s[1],
      dynamicModes: s[2],
      dynamicStrings: s[3],
      peakPa: s[4],
      clippedSamples: s[5],
      forcedRetirements: s[6],
      retirements: s[7],
      contactFailures: s[8],
      rebuilds: s[9],
      noteOns: s[10],
      refusedNoteOns: s[11],
      renderedFrames: rendered,
      playedFrames: level.playedFrames,
      quality: { ...granted },
      governorSteps,
    });
    windowWallMs = 0;
    windowAudioMs = 0;
    worstChunkMs = 0;
  }
}

function onSinkMessage(event) {
  const message = event.data;
  if (message && message.type === 'level') {
    level = { ...message, at: performance.now() };
    pump();
  }
}

function makeEngine(options) {
  return new PianoEngine(
    sampleRate,
    Number(options.boardCoordinates) || 0,
    Number(options.courseModeBudget) || 0,
    Number(options.maxDynamicModes) || 0,
    Number(options.fullScalePa) || 0,
  );
}

// Offline speed test on a separate instance: real-time factor versus the
// number of simultaneously sounding keys. Nothing is sent to the sink.
function bench(args) {
  const keys = Array.isArray(args.keys) && args.keys.length ? args.keys
    : [36, 43, 48, 52, 55, 60, 64, 67, 72, 76, 79, 84, 31, 40, 57, 62, 24, 29, 33, 38, 45, 50, 53, 69, 74, 81, 88, 91];
  const counts = Array.isArray(args.counts) && args.counts.length ? args.counts : [1, 2, 4, 6, 8, 10, 12, 16];
  const seconds = Number(args.seconds) > 0 ? Number(args.seconds) : 1.5;
  const speed = Number(args.speedMS) > 0 ? Number(args.speedMS) : 2.5;
  const results = [];
  for (const count of counts) {
    const probe = makeEngine({ ...engineOptions, maxDynamicModes: 1e9 });
    probe.sustain(args.sustain === undefined ? 1 : Number(args.sustain));
    for (const key of keys.slice(0, count)) probe.noteOn(key, speed);
    const chunks = Math.round((seconds * sampleRate) / CHUNK);
    let worst = 0;
    let peak = 0;
    const t0 = performance.now();
    for (let c = 0; c < chunks; c++) {
      const tc = performance.now();
      const out = probe.render(CHUNK);
      worst = Math.max(worst, performance.now() - tc);
      for (let i = 0; i < out.length; i++) { const a = Math.abs(out[i]); if (a > peak) peak = a; }
    }
    const wallMs = performance.now() - t0;
    const s = probe.stats();
    results.push({
      keys: count,
      dynamicStrings: s[3],
      dynamicModes: s[2],
      soundingKeys: s[1],
      realtimeFactor: (chunks * CHUNK * 1000) / sampleRate / wallMs,
      worstChunkMs: worst,
      chunkMs: (CHUNK * 1000) / sampleRate,
      peakFullScale: peak,
      peakPa: s[4],
      error: probe.lastError() ?? null,
    });
    probe.free();
  }
  return results;
}

function command(message) {
  const { name, id } = message;
  const args = Array.isArray(message.args) ? message.args : message.args === undefined ? [] : [message.args];
  if (!engine) { fail('engine is not initialised', { id, command: name }); return; }
  switch (name) {
    case 'noteOn': {
      if (!engine.noteOn(Number(args[0]), Number(args[1]))) fail(`noteOn(${args[0]}, ${args[1]}) refused: key must be 21..108 and hammer velocity in (0, 8] m/s`, { id, command: name, code: 'refused' });
      break;
    }
    case 'noteOff': engine.noteOff(Number(args[0])); break;
    case 'sustain': if (!engine.sustain(Number(args[0]))) fail('sustain must lie in [0,1]', { id, command: name, code: 'refused' }); break;
    case 'allNotesOff': engine.allNotesOff(); break;
    case 'setQuality': {
      const q = args[0] || {};
      asked = {
        courseModeBudget: Number(q.courseModeBudget) || asked.courseModeBudget,
        partialBandHz: Number(q.partialBandHz) || asked.partialBandHz,
        maxDynamicModes: Number(q.maxDynamicModes) || asked.maxDynamicModes,
      };
      granted = { ...asked };
      if (!engine.setQuality(asked.courseModeBudget, asked.partialBandHz, asked.maxDynamicModes)) fail('invalid quality limits', { id, command: name, code: 'refused' });
      break;
    }
    case 'setTargetMs': {
      const ms = Number(args[0]);
      if (ms >= 5 && ms <= 500) targetFrames = Math.round((ms * sampleRate) / 1000); else fail('targetMs must lie in 5..500', { id, command: name, code: 'refused' });
      break;
    }
    case 'setSnapshotHz': {
      const hz = Number(args[0]);
      snapshotIntervalMs = hz > 0 ? 1000 / Math.min(hz, 120) : 0;
      break;
    }
    case 'setGovernor': governor = Boolean(args[0]); break;
    case 'partials': {
      const value = engine.partialsHz(Number(args[0]));
      post({ type: 'reply', id, name, value }, [value.buffer]);
      return;
    }
    case 'bench': {
      try { post({ type: 'bench', id, results: bench(args[0] || {}) }); } catch (error) { fail(error, { id, command: name }); }
      return;
    }
    case 'dispose': {
      if (timer) clearInterval(timer);
      timer = 0;
      engine.free();
      engine = null;
      if (sinkPort) sinkPort.onmessage = null;
      post({ type: 'reply', id, name, value: true });
      self.close();
      return;
    }
    default:
      fail(`unknown command ${name}`, { id, command: name });
      return;
  }
  // Serve the strike right away when the queue has room.
  pump();
}

async function start(message) {
  const t0 = performance.now();
  sampleRate = Number(message.sampleRate) || 48000;
  engineOptions = message.options || {};
  if (engineOptions.targetMs >= 5) targetFrames = Math.round((engineOptions.targetMs * sampleRate) / 1000);
  else targetFrames = Math.round(0.018 * sampleRate);
  if (engineOptions.snapshotHz !== undefined) snapshotIntervalMs = engineOptions.snapshotHz > 0 ? 1000 / engineOptions.snapshotHz : 0;
  if (engineOptions.governor !== undefined) governor = Boolean(engineOptions.governor);
  await init(message.wasmUrl ? { module_or_path: message.wasmUrl } : undefined);
  engine = makeEngine(engineOptions);
  info = JSON.parse(engine.info());
  asked = { courseModeBudget: info.courseModeBudget, partialBandHz: info.partialBandHz, maxDynamicModes: info.maxDynamicModes };
  granted = { ...asked };
  const xy = engine.boardSurfaceXY();
  const shapes = engine.boardSurfaceShapes();
  const modeHz = engine.boardModeHz();
  const keyGeometry = engine.keyGeometry();
  sinkPort = message.port || null;
  if (sinkPort) {
    sinkPort.onmessage = onSinkMessage;
    if (sinkPort.start) sinkPort.start();
  }
  post({ type: 'ready', info, board: { xy, shapes, modeHz }, keyGeometry, initMs: performance.now() - t0 },
    [xy.buffer, shapes.buffer, modeHz.buffer, keyGeometry.buffer]);
  // Prime the FIFO with one buffer of the whole target (silence computed by
  // the engine itself), so playback never starts on an empty queue.
  if (sinkPort) {
    const prime = Math.min(8192, Math.ceil(targetFrames / CHUNK) * CHUNK);
    const samples = engine.render(prime);
    sinkPort.postMessage({ type: 'audio', channels: [samples] }, [samples.buffer]);
    sentTotal += prime;
    rendered += prime;
  }
  // The sink's level reports drive the pump; this slow timer only covers a
  // suspended or not-yet-started audio clock. No busy spin.
  timer = setInterval(pump, 8);
  pump();
}

self.onmessage = (event) => {
  const message = event.data;
  if (!message || typeof message !== 'object') return;
  try {
    if (message.type === 'init') {
      start(message).catch((error) => fail(error, { command: 'init' }));
    } else if (message.type === 'cmd') {
      command(message);
    }
  } catch (error) {
    fail(error, { id: message.id, command: message.name });
  }
};
