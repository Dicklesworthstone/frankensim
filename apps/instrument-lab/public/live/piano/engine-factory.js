// Shared by standalone and scheduled-lane playback, including their timing probes.
// This chooses an existing WASM constructor, never a JavaScript audio resampler.
export function createPianoEngine(Engine, sampleRate, options = {}, bytes = null) {
  const substeps = options.substeps === undefined ? 1 : Number(options.substeps);
  if (![1, 2, 4].includes(substeps) || !Number.isInteger(sampleRate)
      || sampleRate < 8000 || sampleRate * substeps > 192000) {
    throw new RangeError('piano substeps must be 1, 2 or 4 with a mechanical rate at most 192000 Hz');
  }
  const limits = [Number(options.courseModeBudget) || 0, Number(options.maxDynamicModes) || 0, Number(options.fullScalePa) || 0];
  const coords = Number(options.boardCoordinates) || 0;
  if (options.scaleCsv !== undefined || options.stretchTuning !== undefined || options.gas !== undefined) {
    let scale = '';
    if (options.scaleCsv !== undefined) {
      if (typeof options.scaleCsv !== 'string' || !options.scaleCsv.trim()
          || new TextEncoder().encode(options.scaleCsv).length > 256 * 1024) {
        throw new RangeError('scaleCsv must be nonempty physical SI CSV text of at most 256 KiB');
      }
      scale = options.scaleCsv;
    }
    const tuning = tuningValues(options.stretchTuning);
    if (options.gas !== undefined) {
      if (typeof options.gas !== 'string' || !options.gas.trim()
          || new TextEncoder().encode(options.gas).length > 8192) {
        throw new RangeError('gas must be a nonempty material atmosphere declaration of at most 8 KiB');
      }
      if (typeof Engine.withAtmosphere !== 'function') {
        throw new Error('selected atmosphere requires rebuilding the piano WASM with withAtmosphere support');
      }
      return Engine.withAtmosphere(sampleRate, bytes ?? new Uint8Array(),
        bytes ? 0 : (coords > 16 ? 0 : coords), ...limits, substeps, scale, tuning, options.gas);
    }
    if (typeof Engine.withPreparation !== 'function') {
      throw new Error('selected physical scale preparation requires rebuilding the piano WASM with withPreparation support');
    }
    return Engine.withPreparation(sampleRate, bytes ?? new Uint8Array(),
      bytes ? 0 : (coords > 16 ? 0 : coords), ...limits, substeps, scale, tuning);
  }
  // Preserve the original calls at 1x, including compatibility with deployed
  // WASM. An explicitly higher clock MUST NOT silently fall back on old code.
  if (substeps === 1) {
    return bytes ? Engine.withBoard(sampleRate, bytes, ...limits)
      : new Engine(sampleRate, coords > 16 ? 0 : coords, ...limits);
  }
  const method = bytes ? 'withBoardOversampling' : 'withOversampling';
  if (typeof Engine[method] !== 'function') {
    throw new Error('selected piano oversampling requires rebuilding the piano WASM with PianoAudio support');
  }
  return Engine[method](sampleRate, bytes ?? (coords > 16 ? 0 : coords), ...limits, substeps);
}

// Input admission only. The Rust owner alone computes all tensions and partials.
function tuningValues(value) {
  if (value === undefined) return new Float64Array();
  const fields = ['concertAHz', 'bassOctaveWeights', 'trebleOctaveWeights', 'maxPartialHz'];
  if (!value || typeof value !== 'object' || Array.isArray(value)
      || Object.keys(value).some((key) => !fields.includes(key))) {
    throw new TypeError('stretchTuning needs explicit concertAHz and optional bass/treble weights and maxPartialHz');
  }
  if (!Number.isFinite(value.concertAHz) || value.concertAHz < 430 || value.concertAHz > 450) {
    throw new RangeError('stretchTuning.concertAHz must be finite and in 430..450 Hz');
  }
  const weights = (input) => {
    if (input === undefined) return [1, 0, 0, 0];
    if (!(Array.isArray(input) || ArrayBuffer.isView(input)) || input.length !== 4) {
      throw new TypeError('octave weights need four numeric entries for 2:1, 4:2, 6:3 and 8:4');
    }
    const entries = Array.from(input); // holes cannot pass Array.every unnoticed
    const scale = Math.max(...entries);
    if (entries.some((x) => !Number.isFinite(x) || x < 0) || !(scale > 0)
        || entries.some((x) => x > 0 && x / scale === 0)) {
      throw new RangeError('octave weights must be finite, nonnegative and representable, with a positive entry');
    }
    return entries;
  };
  let band = 0;
  if (value.maxPartialHz !== undefined) {
    if (!Number.isFinite(value.maxPartialHz) || value.maxPartialHz <= 0) {
      throw new RangeError('stretchTuning.maxPartialHz must be finite and positive');
    }
    band = value.maxPartialHz;
  }
  return new Float64Array([value.concertAHz, ...weights(value.bassOctaveWeights),
    ...weights(value.trebleOctaveWeights), band]);
}
