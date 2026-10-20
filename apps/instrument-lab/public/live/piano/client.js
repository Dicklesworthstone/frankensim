// Framework-free client of the live piano engine (plain ES module; the typed
// twin lives in src/live/pianoClient.ts). It wires
//   AudioContext -> AudioWorkletNode('fs-live-sink')  <-- MessagePort -->  piano/worker.js (wasm)
// and exposes note/pedal commands plus snapshot/stats callbacks.
//
// Units: strikes are POST-ESCAPEMENT hammer velocities in m/s; audio samples
// are pressure [Pa] / info.fullScalePa at the declared Rayleigh receiver.

/** MIDI velocity 127 maps to this hammer speed; linear and UNCALIBRATED (grand_piano/midi.rs). */
export const MIDI_MAX_VELOCITY_M_S = 4.5;
export const midiVelocityToMS = (velocity) => (MIDI_MAX_VELOCITY_M_S * Math.max(0, Math.min(127, velocity))) / 127;

const KEY_HEADER = 9;

/** Decode the flat Float32Array snapshot posted by the worker. */
export function parseSnapshot(data) {
  const keyCount = data[0];
  const stations = data[1];
  const boardModes = data[2];
  // The bank adds one signed exchange value after the original 18 header values.
  const boardStart = data.length - boardModes - keyCount * (KEY_HEADER + stations);
  if (boardStart !== 18 && boardStart !== 19) throw new RangeError('unsupported piano snapshot layout; rebuild the piano WASM and client together');
  const keys = [];
  let offset = boardStart + boardModes;
  for (let k = 0; k < keyCount; k++) {
    keys.push({
      midi: data[offset],
      held: data[offset + 1] > 0,
      damped: data[offset + 2] > 0,
      hammerPositionM: data[offset + 3], // 0 = undeflected string plane; rest at -0.020 m
      hammerForceN: data[offset + 4],
      relativeEnergyJ: data[offset + 5],
      partials: data[offset + 6],
      strings: data[offset + 7],
      bridgeDisplacementM: data[offset + 8],
      // Transverse displacement [m] of the key's first string at x/L = (j + 0.5)/stations.
      stringShapeM: data.subarray(offset + KEY_HEADER, offset + KEY_HEADER + stations),
    });
    offset += KEY_HEADER + stations;
  }
  return {
    sustain: data[3],
    ledger: {
      inputJ: data[4],
      storedJ: data[5],
      dissipatedJ: data[6],
      residualJ: data[7], // input + bank exchange - stored - dissipated - retired: what the books do not explain
      stringLossJ: data[8],
      boardLossJ: data[9],
      damperLossJ: data[10],
      hammerLossJ: data[11], // felt (crush + relaxation) and shank
      catchLossJ: data[12],
      retiredJ: data[13], // relative string energy discarded when keys went back to the static image
      couplingJ: boardStart === 19 ? data[18] : 0, // signed numerical work of the bank's explicit exchange
    },
    dynamicModes: data[14],
    dynamicStrings: data[15],
    peakPa: data[16],
    forcedRetirements: data[17],
    // Bare soundboard modal coordinates [m sqrt(kg)]; with board.shapes they give w(x, y).
    boardCoordinates: data.subarray(boardStart, boardStart + boardModes),
    keys,
    raw: data,
  };
}

export class PianoClient {
  /**
   * @param {AudioContext} audioContext  running or suspended context (48 kHz recommended)
   * @param {object} [options] {baseUrl='/live/', destination, boardCoordinates, courseModeBudget,
   *   maxDynamicModes, fullScalePa, targetMs, snapshotHz, governor, sinkUrl, workerUrl}
   */
  static async create(audioContext, options = {}) {
    const base = new URL(options.baseUrl ?? '/live/', globalThis.location?.href);
    const sinkUrl = options.sinkUrl ?? new URL('sink-worklet.js', base).href;
    const workerUrl = options.workerUrl ?? new URL('piano/worker.js', base).href;
    await audioContext.audioWorklet.addModule(sinkUrl);
    const node = new AudioWorkletNode(audioContext, 'fs-live-sink', { numberOfInputs: 0, numberOfOutputs: 1, outputChannelCount: [2] });
    node.connect(options.destination ?? audioContext.destination);
    const channel = new MessageChannel();
    node.port.postMessage({ type: 'connect', port: channel.port2 }, [channel.port2]);
    const worker = new Worker(workerUrl, { type: 'module' });
    const client = new PianoClient(audioContext, node, worker);
    const { baseUrl, destination, sinkUrl: _s, workerUrl: _w, ...engineOptions } = options;
    await new Promise((resolve, reject) => {
      client._ready = { resolve, reject };
      worker.postMessage({ type: 'init', sampleRate: audioContext.sampleRate, port: channel.port1, options: engineOptions }, [channel.port1]);
    });
    return client;
  }

  constructor(audioContext, node, worker) {
    this.audioContext = audioContext;
    this.node = node;
    this.worker = worker;
    this.info = null;
    this.board = null;
    this.initMs = 0;
    this.lastStats = null;
    this._snapshot = new Set();
    this._stats = new Set();
    this._error = new Set();
    this._pending = new Map();
    this._nextId = 1;
    this._ready = null;
    worker.onmessage = (event) => this._receive(event.data);
    worker.onerror = (event) => this._fail(new Error(event.message || 'piano worker failed to load'));
  }

  _fail(error) {
    if (this._ready) { this._ready.reject(error); this._ready = null; }
    for (const listener of this._error) listener(error.message);
  }

  _receive(message) {
    switch (message.type) {
      case 'ready':
        this.info = message.info;
        this.board = message.board;
        this.initMs = message.initMs;
        if (this._ready) { this._ready.resolve(); this._ready = null; }
        break;
      case 'snapshot':
        if (this._snapshot.size) {
          const parsed = parseSnapshot(message.data);
          parsed.renderedFrames = message.renderedFrames;
          parsed.aheadMs = message.aheadMs;
          for (const listener of this._snapshot) listener(parsed);
        }
        break;
      case 'stats':
        this.lastStats = message;
        for (const listener of this._stats) listener(message);
        break;
      case 'reply':
      case 'bench': {
        const pending = this._pending.get(message.id);
        if (pending) { this._pending.delete(message.id); pending.resolve(message.type === 'bench' ? message.results : message.value); }
        break;
      }
      case 'error': {
        const pending = message.id !== undefined ? this._pending.get(message.id) : null;
        if (pending) { this._pending.delete(message.id); pending.reject(new Error(message.message)); }
        if (this._ready && message.command === 'init') { this._ready.reject(new Error(message.message)); this._ready = null; }
        for (const listener of this._error) listener(message.message);
        break;
      }
      default:
        break;
    }
  }

  _send(name, ...args) { this.worker.postMessage({ type: 'cmd', name, args }); }
  _request(name, ...args) {
    const id = this._nextId++;
    return new Promise((resolve, reject) => {
      this._pending.set(id, { resolve, reject });
      this.worker.postMessage({ type: 'cmd', name, args, id });
    });
  }

  /** Strike a key (MIDI 21..108) with a post-escapement hammer velocity in m/s, (0, 8]. */
  noteOn(midi, speedMS) { this._send('noteOn', midi, speedMS); }
  /** Strike with a MIDI velocity through the reference's uncalibrated linear map (127 -> 4.5 m/s). */
  noteOnMidi(midi, velocity) { if (velocity > 0) this._send('noteOn', midi, midiVelocityToMS(velocity)); else this._send('noteOff', midi); }
  noteOff(midi) { this._send('noteOff', midi); }
  /** Sustain pedal travel 0..1 (damper drag scales with (1 - travel)^2). */
  sustain(value) { this._send('sustain', Math.max(0, Math.min(1, Number(value)))); }
  allNotesOff() { this._send('allNotesOff'); }
  /** Limits for keys struck afterwards: {courseModeBudget, partialBandHz, maxDynamicModes}. */
  setQuality(quality) { this._send('setQuality', quality); }
  setTargetMs(ms) { this._send('setTargetMs', ms); }
  setSnapshotHz(hz) { this._send('setSnapshotHz', hz); }
  setGovernor(on) { this._send('setGovernor', Boolean(on)); }
  /** Stiff-string partial frequencies [Hz] the live model retains for a key. */
  partials(midi) { return this._request('partials', midi); }
  /** Offline speed test on a separate engine instance: real-time factor vs sounding keys. */
  bench(options = {}) { return this._request('bench', options); }

  onSnapshot(listener) { this._snapshot.add(listener); return () => this._snapshot.delete(listener); }
  onStats(listener) { this._stats.add(listener); return () => this._stats.delete(listener); }
  onError(listener) { this._error.add(listener); return () => this._error.delete(listener); }

  /** Estimated key-to-ear latency [ms]: queue + receiver travel/filter + output device. */
  latencyEstimateMs() {
    const queued = this.lastStats ? this.lastStats.queuedMs : 0;
    const receiver = this.info ? this.info.microphoneLatencyMs : 0;
    const output = ((this.audioContext.baseLatency || 0) + (this.audioContext.outputLatency || 0)) * 1000;
    return { queuedMs: queued, receiverMs: receiver, outputMs: output, totalMs: queued + receiver + output };
  }

  dispose() {
    try { this.worker.postMessage({ type: 'cmd', name: 'dispose' }); } catch { /* already gone */ }
    setTimeout(() => this.worker.terminate(), 200);
    try { this.node.port.postMessage({ type: 'flush' }); } catch { /* closed */ }
    try { this.node.disconnect(); } catch { /* already disconnected */ }
    this._snapshot.clear();
    this._stats.clear();
    this._error.clear();
  }
}
