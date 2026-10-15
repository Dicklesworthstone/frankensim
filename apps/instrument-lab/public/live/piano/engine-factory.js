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
