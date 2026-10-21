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
  if (options.scaleCsv !== undefined || options.stretchTuning !== undefined || options.gas !== undefined
      || options.stringMaterials !== undefined) {
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
    }
    if (options.stringMaterials !== undefined) {
      const materials = stringMaterialValues(options.stringMaterials);
      if (typeof Engine.withStringMaterials !== 'function') {
        throw new Error('selected string materials require rebuilding the piano WASM with withStringMaterials support');
      }
      return Engine.withStringMaterials(sampleRate, bytes ?? new Uint8Array(),
        bytes ? 0 : (coords > 16 ? 0 : coords), ...limits, substeps, scale, tuning,
        options.gas ?? '293.15 101325 0', ...materials);
    }
    if (options.gas !== undefined) {
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

// Frame unchanged canonical packs, including byte-offset views. Rust owns
// source admission, the concentric-section integrals and complete scale coverage.
function stringMaterialValues(value) {
  const fields = ['packs', 'geometry', 'temperatureAxis', 'temperatureK',
    'viscosityProperty', 'frequencyAxis', 'frequencyHz'];
  if (!value || typeof value !== 'object' || Array.isArray(value)
      || Object.keys(value).some((key) => !fields.includes(key))) {
    throw new TypeError('stringMaterials needs explicit packs, geometry, temperature and frequency axes/values, and viscosityProperty');
  }
  const { packs, geometry, temperatureAxis, temperatureK, viscosityProperty, frequencyAxis, frequencyHz } = value;
  if (!Array.isArray(packs) || packs.length < 1 || packs.length > 128) {
    throw new RangeError('stringMaterials.packs must contain 1..128 nonempty canonical material packs');
  }
  let total = 0;
  for (const pack of packs) {
    if (!(pack instanceof Uint8Array || pack instanceof ArrayBuffer) || !pack.byteLength) {
      throw new TypeError('each stringMaterials pack must be a nonempty Uint8Array or ArrayBuffer');
    }
    total += pack.byteLength;
    if (total > 16 * 1024 * 1024) throw new RangeError('combined stringMaterials packs must fit 16 MiB');
  }
  const encoder = new TextEncoder();
  if (typeof geometry !== 'string' || !geometry.trim() || encoder.encode(geometry).length > 256 * 1024) {
    throw new RangeError('stringMaterials.geometry must be nonempty CSV text of at most 256 KiB');
  }
  for (const name of ['temperatureAxis', 'frequencyAxis', 'viscosityProperty']) {
    if (typeof value[name] !== 'string' || !value[name] || value[name].trim() !== value[name]
        || encoder.encode(value[name]).length > 128) {
      throw new TypeError(`stringMaterials.${name} must be an explicit trimmed name of at most 128 UTF-8 bytes`);
    }
  }
  if (temperatureAxis === frequencyAxis) throw new TypeError('stringMaterials temperature and frequency axes must be distinct');
  if (!Number.isFinite(temperatureK) || temperatureK <= 0 || !Number.isFinite(frequencyHz) || frequencyHz <= 0) {
    throw new RangeError('stringMaterials temperatureK and frequencyHz must be explicit finite positive numbers');
  }
  const bytes = new Uint8Array(total), lengths = new Uint32Array(packs.length);
  let offset = 0;
  for (let i = 0; i < packs.length; i++) {
    const pack = packs[i] instanceof Uint8Array ? packs[i] : new Uint8Array(packs[i]);
    bytes.set(pack, offset);
    lengths[i] = pack.byteLength;
    offset += pack.byteLength;
  }
  return [bytes, lengths, geometry, temperatureAxis, temperatureK, viscosityProperty, frequencyAxis, frequencyHz];
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
