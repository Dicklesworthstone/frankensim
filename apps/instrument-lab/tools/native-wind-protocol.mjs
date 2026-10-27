// Transport limits only. music_render owns the complete source grammar and physics.
export const MAX_SOURCE_BYTES = 1024 * 1024;
export const MAX_OUTPUT_BYTES = 16 * 1024 * 1024;
export const MAX_LOG_BYTES = 256 * 1024;
export const DEFAULT_TIMEOUT_MS = 60000;
// Supported editor families. The native --self-contained loader policy, rather
// than the schema name, refuses any external material/artifact callbacks.
export const SUPPORTED_SCHEMAS = ['frankensim-reed-performance-v1', 'frankensim-flue-dg-performance-v1'];

export class NativeWindError extends Error {
  constructor(status, message) { super(message); this.status = status; }
}

export function validateSource(bytes) {
  if (!bytes.length) throw new NativeWindError(400, 'Choose or enter a wind performance file.');
  if (bytes.length > MAX_SOURCE_BYTES) throw new NativeWindError(413, 'Performance exceeds the 1 MiB input limit.');
  let text;
  try { text = new TextDecoder('utf-8', {fatal: true, ignoreBOM: true}).decode(bytes); }
  catch { throw new NativeWindError(400, 'Performance must contain valid UTF-8 bytes.'); }
  if (!SUPPORTED_SCHEMAS.includes(text.split(/\r?\n/, 1)[0]))
    throw new NativeWindError(400, 'This text-only renderer accepts self-contained reed and flue-DG performances. Pack-backed sources require the native CLI.');
  // Do not trim, normalize newlines, parse records, or re-encode the admitted bytes.
  return bytes;
}

export function validateWave(bytes, limit = MAX_OUTPUT_BYTES) {
  const bad = () => { throw new NativeWindError(502, 'Native renderer did not produce a complete PCM16 WAV.'); };
  if (bytes.length > limit) throw new NativeWindError(502, 'Native WAV exceeds the output limit.');
  if (bytes.length < 44 || bytes.toString('ascii', 0, 4) !== 'RIFF'
      || bytes.toString('ascii', 8, 12) !== 'WAVE' || bytes.readUInt32LE(4) + 8 !== bytes.length) bad();
  let pos = 12, format = false, data = false, alignment = 0, dataLength = 0;
  while (pos + 8 <= bytes.length) {
    const name = bytes.toString('ascii', pos, pos + 4), size = bytes.readUInt32LE(pos + 4);
    const end = pos + 8 + size;
    if (end > bytes.length) bad();
    if (name === 'fmt ') {
      if (format || size < 16 || bytes.readUInt16LE(pos + 8) !== 1
          || bytes.readUInt16LE(pos + 10) !== 1 || bytes.readUInt32LE(pos + 12) !== 48000
          || bytes.readUInt16LE(pos + 22) !== 16) bad();
      alignment = bytes.readUInt16LE(pos + 20);
      if (alignment !== 2 || bytes.readUInt32LE(pos + 16) !== 96000) bad();
      format = true;
    } else if (name === 'data') {
      if (data || !size) bad();
      data = true; dataLength = size;
    }
    pos = end + (size & 1);
  }
  if (pos !== bytes.length || !format || !data || dataLength % alignment) bad();
  return bytes;
}
