// Live drum engine worker (ES module Worker).
//
// Owns the WebAssembly engine (FrankenSim's own drum mechanics), renders ahead
// in small chunks to keep a short FIFO in the `fs-live-sink` AudioWorklet, and
// reports state. It never post-processes audio: samples are pressure [Pa] at a
// declared microphone divided by the declared fixed full-scale pressure.
//
// Main thread -> worker
//   {type:'init', sampleRate, port, params?, autoCalibrate?, calibration?, targetMs?,
//    adaptiveTarget?, snapshotHz?, monitor?, engine?: 'auto'|'baseline'}
//   {type:'cmd', id?, name, args?}   names: configure | strike | setSubsteps | calibrate | bench |
//                                    reset | setMonitor | setTargetMs | setSnapshotHz | info | dispose
// Worker -> main thread
//   {type:'ready', info, mesh:{nodes, triangles}, modeShapes:[batter, resonant], buildMs,
//    eigenSolved, rebuilt, calibration?}
//   {type:'snapshot', state:Float64Array, displacements:Float32Array, renderedFrames, aheadMs}
//   {type:'stats', realtimeFactor, worstChunkMs, queuedMs, targetMs, targetRaises, underruns, underrunFrames,
//    renderedFrames, playedFrames, substeps, mechanicalRateHz, monitor, engineVariant, halted}
//   {type:'calibration', results:[{substeps, realtimeFactor, measured, ...}], chosen, minFactor, sustained}
//   {type:'strike', id?, outcome:'queued'|'deferred', x, y, speed}
//   {type:'done', id, name}            acknowledgement of a command that has no other reply
//   {type:'error', id?, command?, code, message}
// Two builds of the same engine, bit-identical in every result:
//   engine_fma.*  hardware fused multiply-add through WebAssembly relaxed-simd
//   engine.*      baseline: runs everywhere, software fused multiply-add (slower)
let engine = null;
let engineVariant = null;
async function loadEngine(message) {
  const order = message.engine === 'baseline' ? ['engine'] : ['engine_fma', 'engine'];
  let lastError = null;
  for (const name of order) {
    try {
      const module = await import(`./${name}.js`);
      await module.default();
      engine = module;
      engineVariant = name === 'engine_fma' ? 'relaxed-simd' : 'baseline';
      return;
    } catch (error) {
      // A browser without relaxed-simd rejects the first module at compile
      // time; the baseline build is then the honest (slower) choice.
      lastError = error;
    }
  }
  throw Object.assign(new Error(`could not load the WebAssembly engine: ${lastError && lastError.message}`), { code: 'wasm_load' });
}

const OUTPUT_RATE = 48000;
const CHUNK = 128;
const MAX_EXTRAPOLATED_FRAMES = 512;

let drum = null;
let sinkPort = null;
let sampleRate = OUTPUT_RATE;
let monitor = 'batter';
// 12 ms queued: measured 0 underruns over 10 s soaks in headless Chrome at
// 20, 12 and 8 ms. If a device does underrun, the target is raised 4 ms at a
// time (never above 40 ms) and reported; it is never lowered silently.
let targetFrames = Math.round(0.012 * OUTPUT_RATE);
let adaptiveTarget = true;
let targetRaises = 0;
let underrunsAtLastStats = 0;
let snapshotIntervalMs = 1000 / 30;
let timer = 0;
let sentTotal = 0;      // frames ever sent to the sink (matches its cumulative receivedFrames)
let rendered = 0;       // frames rendered since the last (re)build
let level = { queuedFrames: 0, receivedFrames: 0, playedFrames: 0, underruns: 0, underrunFrames: 0, at: 0 };
let epoch = 0;
let lastSnapshotAt = 0;
let lastStatsAt = 0;
let windowWallMs = 0;
let windowAudioMs = 0;
let worstChunkMs = 0;
let halted = null;
let info = null;

const post = (message, transfer) => self.postMessage(message, transfer ?? []);

function fail(error, extra = {}) {
  const code = error && error.code ? String(error.code) : 'worker_error';
  const message = error && error.message ? String(error.message) : String(error);
  post({ type: 'error', code, message, ...extra });
}

function estimateQueue(now) {
  const inFlight = sentTotal - level.receivedFrames;
  // Extrapolate playback only across a couple of report intervals. If the
  // sink stops reporting (suspended AudioContext, closed tab audio) the
  // estimate stops falling, so the engine pauses instead of rendering an
  // unbounded backlog that would later play as stale, late audio.
  const sinceReport = level.at > 0 ? ((now - level.at) * sampleRate) / 1000 : 0;
  const consumed = Math.min(sinceReport, MAX_EXTRAPOLATED_FRAMES);
  return Math.max(0, level.queuedFrames + inFlight - consumed);
}

function renderChunk() {
  const t0 = performance.now();
  const planar = drum.render(CHUNK);
  const elapsed = performance.now() - t0;
  windowWallMs += elapsed;
  windowAudioMs += (CHUNK * 1000) / OUTPUT_RATE;
  if (elapsed > worstChunkMs) worstChunkMs = elapsed;
  const refusal = drum.take_refusal();
  if (refusal) {
    const parsed = JSON.parse(refusal);
    halted = parsed;
    post({ type: 'error', code: parsed.code, message: `${parsed.message} (the voice is halted and silent; send reset or configure)` });
  }
  let channels;
  if (monitor === 'split') {
    channels = [planar.slice(0, CHUNK), planar.slice(CHUNK, 2 * CHUNK)];
  } else if (monitor === 'resonant') {
    channels = [planar.slice(CHUNK, 2 * CHUNK)];
  } else if (monitor === 'sum') {
    // Declared electrical sum of the two microphone signals, unit gains.
    const mix = planar.slice(0, CHUNK);
    for (let i = 0; i < CHUNK; i++) mix[i] += planar[CHUNK + i];
    channels = [mix];
  } else {
    channels = [planar.slice(0, CHUNK)];
  }
  sinkPort.postMessage({ type: 'audio', channels }, channels.map((c) => c.buffer));
  sentTotal += CHUNK;
  rendered += CHUNK;
}

function pump() {
  if (!drum || !sinkPort) return;
  let now = performance.now();
  let queue = estimateQueue(now);
  // Bounded work per turn so commands are still served when the engine is
  // slower than real time (the sink then reports honest underruns).
  let budget = Math.ceil(targetFrames / CHUNK) + 2;
  while (queue < targetFrames && budget-- > 0) {
    renderChunk();
    queue += CHUNK;
  }
  now = performance.now();
  if (now - lastSnapshotAt >= snapshotIntervalMs) {
    lastSnapshotAt = now;
    const state = drum.snapshot();
    const displacements = drum.displacements();
    post({ type: 'snapshot', state, displacements, renderedFrames: rendered, aheadMs: (estimateQueue(now) * 1000) / sampleRate },
      [state.buffer, displacements.buffer]);
  }
  if (now - lastStatsAt >= 250) {
    lastStatsAt = now;
    if (adaptiveTarget && level.underruns > underrunsAtLastStats && targetFrames < 0.04 * OUTPUT_RATE) {
      targetFrames = Math.min(Math.round(0.04 * OUTPUT_RATE), targetFrames + Math.round(0.004 * OUTPUT_RATE));
      targetRaises += 1;
    }
    underrunsAtLastStats = level.underruns;
    post({
      type: 'stats',
      realtimeFactor: windowWallMs > 0 ? windowAudioMs / windowWallMs : null,
      worstChunkMs,
      queuedMs: (estimateQueue(now) * 1000) / sampleRate,
      targetMs: (targetFrames * 1000) / sampleRate,
      targetRaises,
      underruns: level.underruns,
      underrunFrames: level.underrunFrames,
      renderedFrames: rendered,
      playedFrames: level.playedFrames,
      substeps: info ? info.substeps : null,
      mechanicalRateHz: info ? info.mechanical_rate_hz : null,
      monitor,
      engineVariant,
      halted,
    });
    windowWallMs = 0;
    windowAudioMs = 0;
    worstChunkMs = 0;
  }
}

function onSinkMessage(event) {
  const message = event.data;
  if (message && message.type === 'level') {
    // Levels from before our last flush describe audio that no longer exists.
    // A newer epoch means another producer flushed the sink: adopt it.
    if (message.epoch < epoch) return;
    epoch = message.epoch;
    level = { ...message, at: performance.now() };
    pump();
  }
}

// The old instrument stops (queued audio is dropped) and the new one starts
// silent. Nothing is cross-faded: a fade would be invented sound.
function restartStream() {
  epoch += 1;
  sinkPort.postMessage({ type: 'flush', epoch, primeFrames: Math.max(0, targetFrames - CHUNK) });
  // Until the sink confirms the flush, assume everything already delivered was dropped.
  level = { queuedFrames: 0, receivedFrames: sentTotal, playedFrames: 0, underruns: 0, underrunFrames: 0, at: 0 };
  rendered = 0;
  underrunsAtLastStats = 0;
  halted = null;
  windowWallMs = 0;
  windowAudioMs = 0;
  worstChunkMs = 0;
}

function publishReady(extra) {
  info = JSON.parse(drum.info());
  const nodes = drum.mesh_nodes();
  const triangles = drum.mesh_triangles();
  const batter = drum.mode_shapes(0);
  const resonant = drum.mode_shapes(1);
  post({ type: 'ready', info, mesh: { nodes, triangles }, modeShapes: [batter, resonant], ...extra },
    [nodes.buffer, triangles.buffer, batter.buffer, resonant.buffer]);
}

// Measure a scratch voice of the current instrument at each candidate rate.
// The scenario is a 20 Hz roll (a stroke every 50 ms): the stick is in the
// system for roughly a sixth of the time, which is harder than single strokes.
function measure(args = {}) {
  const candidates = (args.candidates ?? [16, 8, 4, 2, 1]).filter((s) => [1, 2, 4, 8, 16].includes(s));
  const minFactor = Number.isFinite(args.minFactor) ? args.minFactor : 1.6;
  const frames = Math.max(1200, Math.min(48000, Math.floor(args.frames ?? 9600)));
  const interval = Math.max(0, Math.floor(args.strikeIntervalFrames ?? 2400));
  const radius = info ? info.params.radius_m : 0.17;
  const x = args.x ?? 0.35 * radius;
  const y = args.y ?? 0.0;
  const speed = args.speed ?? 2.0;
  const results = [];
  // A short probe at the cheapest admissible rate bounds how long the real
  // measurements can take on a slow device.
  const admissible = candidates.filter((s) => !info || info.params.band_high_hz / (OUTPUT_RATE * s) < 0.45);
  if (admissible.length === 0) return { results: [], chosen: null, minFactor, sustained: false };
  const cheapest = Math.min(...admissible);
  let perStepMs = null;
  try {
    const probeFrames = 2400;
    const t0 = performance.now();
    const r = drum.bench(cheapest, probeFrames, interval > 0 ? Math.min(interval, 1200) : 0, x, y, speed);
    perStepMs = (performance.now() - t0) / r[1];
  } catch (error) {
    return { results: [{ substeps: cheapest, refused: String(error.message ?? error) }], chosen: null, minFactor, sustained: false };
  }
  for (const substeps of admissible.slice().sort((a, b) => b - a)) {
    const predicted = 1000 / OUTPUT_RATE / (perStepMs * substeps);
    if (predicted < 0.5 * minFactor) {
      results.push({ substeps, measured: false, predictedRealtimeFactor: predicted });
      continue;
    }
    try {
      const t0 = performance.now();
      const r = drum.bench(substeps, frames, interval, x, y, speed);
      const wallMs = performance.now() - t0;
      results.push({
        substeps, measured: true, realtimeFactor: (frames * 1000) / OUTPUT_RATE / wallMs, wallMs,
        frames: r[0], mechanicalSteps: r[1], stepsWithStick: r[2], stepsInContact: r[3], strikes: r[4],
        peakPa: r[5], balanceResidualJ: r[6], nsPerStep: (wallMs * 1e6) / r[1],
      });
    } catch (error) {
      results.push({ substeps, measured: false, refused: String(error.message ?? error) });
    }
  }
  const sustaining = results.filter((r) => r.measured && r.realtimeFactor >= minFactor).map((r) => r.substeps);
  const chosen = sustaining.length > 0 ? Math.max(...sustaining) : cheapest;
  return { results, chosen, minFactor, sustained: sustaining.length > 0, scenario: { frames, strikeIntervalFrames: interval, x, y, speed } };
}

function calibrate(args) {
  const outcome = measure(args);
  if (outcome.chosen && info && outcome.chosen !== info.substeps) {
    drum.rebuild(JSON.stringify({ substeps: outcome.chosen }));
    info = JSON.parse(drum.info());
  }
  return outcome;
}

async function start(message) {
  sampleRate = message.sampleRate;
  sinkPort = message.port;
  if (sampleRate !== OUTPUT_RATE) {
    throw Object.assign(new Error(`the live engine renders ${OUTPUT_RATE} Hz only; create the AudioContext with sampleRate ${OUTPUT_RATE} (got ${sampleRate}). No resampling is applied.`), { code: 'sample_rate' });
  }
  if (typeof message.monitor === 'string') monitor = message.monitor;
  if (Number.isFinite(message.targetMs)) targetFrames = Math.round((message.targetMs * OUTPUT_RATE) / 1000);
  if (message.adaptiveTarget === false) adaptiveTarget = false;
  if (Number.isFinite(message.snapshotHz) && message.snapshotHz > 0) snapshotIntervalMs = 1000 / message.snapshotHz;
  await loadEngine(message);
  const t0 = performance.now();
  drum = new engine.LiveDrum(JSON.stringify(message.params ?? {}));
  const buildMs = performance.now() - t0;
  info = JSON.parse(drum.info());
  let calibration = null;
  if (message.autoCalibrate !== false) {
    calibration = calibrate(message.calibration ?? {});
    post({ type: 'calibration', ...calibration });
  }
  sinkPort.onmessage = onSinkMessage;
  restartStream();
  publishReady({ buildMs, eigenSolved: true, rebuilt: false, calibration, defaults: JSON.parse(engine.default_params()), scope: JSON.parse(engine.scope()), engineVariant });
  timer = setInterval(pump, 4);
  pump();
}

function command(message) {
  const { id, name } = message;
  const args = message.args ?? {};
  switch (name) {
    case 'configure': {
      const t0 = performance.now();
      const eigenSolved = drum.rebuild(JSON.stringify(args.params ?? args));
      const buildMs = performance.now() - t0;
      restartStream();
      publishReady({ id, buildMs, eigenSolved, rebuilt: true });
      break;
    }
    case 'setSubsteps': {
      const t0 = performance.now();
      const eigenSolved = drum.rebuild(JSON.stringify({ substeps: args.substeps }));
      const buildMs = performance.now() - t0;
      restartStream();
      publishReady({ id, buildMs, eigenSolved, rebuilt: true });
      break;
    }
    case 'strike': {
      const outcome = drum.strike(args.x, args.y, args.speed);
      post({ type: 'strike', id, outcome, x: args.x, y: args.y, speed: args.speed });
      pump();
      break;
    }
    case 'calibrate': {
      const outcome = calibrate(args);
      restartStream();
      post({ type: 'calibration', id, ...outcome });
      publishReady({ id, buildMs: 0, eigenSolved: false, rebuilt: true, calibration: outcome });
      break;
    }
    case 'bench': {
      // Blocks this worker while it measures; the sink will report underruns.
      post({ type: 'calibration', id, applied: false, ...measure(args) });
      restartStream();
      break;
    }
    case 'reset':
      drum.reset();
      restartStream();
      post({ type: 'done', id, name });
      break;
    case 'setMonitor':
      if (!['batter', 'resonant', 'split', 'sum'].includes(args.mode)) {
        throw Object.assign(new Error('monitor must be batter, resonant, split or sum'), { code: 'bad_parameter' });
      }
      monitor = args.mode;
      post({ type: 'done', id, name });
      break;
    case 'setTargetMs':
      if (!(args.ms >= 3 && args.ms <= 500)) {
        throw Object.assign(new Error('targetMs must lie in 3..500'), { code: 'bad_parameter' });
      }
      targetFrames = Math.round((args.ms * OUTPUT_RATE) / 1000);
      post({ type: 'done', id, name });
      break;
    case 'setSnapshotHz':
      if (!(args.hz > 0 && args.hz <= 240)) {
        throw Object.assign(new Error('snapshotHz must lie in (0, 240]'), { code: 'bad_parameter' });
      }
      snapshotIntervalMs = 1000 / args.hz;
      post({ type: 'done', id, name });
      break;
    case 'info':
      publishReady({ id, buildMs: 0, eigenSolved: false, rebuilt: false });
      break;
    case 'dispose':
      clearInterval(timer);
      if (sinkPort) sinkPort.postMessage({ type: 'flush' });
      if (drum) drum.free();
      drum = null;
      post({ type: 'done', id, name });
      self.close();
      break;
    default:
      throw Object.assign(new Error(`unknown command "${name}"`), { code: 'unknown_command' });
  }
}

self.onmessage = (event) => {
  const message = event.data;
  if (!message || typeof message !== 'object') return;
  if (message.type === 'init') {
    start(message).catch((error) => fail(error, { command: 'init' }));
  } else if (message.type === 'cmd') {
    if (!drum) {
      fail(Object.assign(new Error('engine is not initialised'), { code: 'not_ready' }), { id: message.id, command: message.name });
      return;
    }
    try {
      command(message);
    } catch (error) {
      fail(error, { id: message.id, command: message.name });
    }
  }
};
