// Transport limits only. music_render owns the complete source grammar and physics.
export const MAX_SOURCE_BYTES = 1024 * 1024;
export const MAX_ATTACHMENTS = 32;
export const MAX_ATTACHMENT_FILE_BYTES = 8 * 1024 * 1024;
export const MAX_ATTACHMENT_BYTES = 32 * 1024 * 1024;
// Base64 source and attachments, plus a bounded path/JSON envelope.
export const MAX_REQUEST_BYTES = 4 * Math.ceil((MAX_SOURCE_BYTES + MAX_ATTACHMENT_BYTES) / 3) + 65536;
export const MAX_OUTPUT_BYTES = 16 * 1024 * 1024;
export const MAX_LOG_BYTES = 256 * 1024;
export const DEFAULT_TIMEOUT_MS = 60000;
// The configured native CLI owns the grammar and all material/artifact decoding.
// --self-contained plus the explicit attached-source list confines native callbacks.
export const SUPPORTED_SCHEMAS = ['frankensim-reed-performance-v1', 'frankensim-flue-dg-performance-v1',
  'frankensim-material-modal-dg-valve-v1', 'frankensim-material-modal-valve-v1',
  'frankensim-plate-valve-performance-v1', 'frankensim-flue-performance-v1'];

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
    throw new NativeWindError(400, 'Choose a native reed, flue or material-valve performance supported by music_render wind.');
  // Do not trim, normalize newlines, parse records, or re-encode the admitted bytes.
  return bytes;
}

function decodedBase64(value, limit) {
  if (typeof value !== 'string' || value.length > 4 * Math.ceil(limit / 3)
      || value.length % 4 !== 0 || !/^[A-Za-z0-9+/]*={0,2}$/.test(value))
    throw new NativeWindError(400, 'Source and attachments require canonical base64 bytes.');
  const bytes = Buffer.from(value, 'base64');
  if (bytes.length > limit || bytes.toString('base64') !== value)
    throw new NativeWindError(400, 'Source or attachment encoding exceeds its byte limit or is invalid.');
  return bytes;
}

/** Transport envelope only: no physics records or material constants are interpreted. */
export function validateInput(bytes, bundled = false) {
  if (!bundled) return {source:validateSource(bytes), attachments:[], decimate:false};
  if (bytes.length > MAX_REQUEST_BYTES) throw new NativeWindError(413, 'Experiment exceeds the complete request limit.');
  let packet;
  try { packet = JSON.parse(new TextDecoder('utf-8', {fatal:true}).decode(bytes)); }
  catch { throw new NativeWindError(400, 'Experiment must be UTF-8 JSON.'); }
  if (packet?.format !== 'frankensim-native-wind-input-v1' || !Array.isArray(packet.attachments)
      || packet.attachments.length > MAX_ATTACHMENTS)
    throw new NativeWindError(400, 'Experiment requires its version and at most 32 source attachments.');
  if (Object.keys(packet).some(key => !['format', 'source', 'attachments', 'decimate'].includes(key))
      || (packet.decimate !== undefined && typeof packet.decimate !== 'boolean'))
    throw new NativeWindError(400, 'Experiment contains unsupported fields or an invalid native decimation option.');
  const source = validateSource(decodedBase64(packet.source, MAX_SOURCE_BYTES));
  const names = new Set(['input.performance', 'output.wav', 'output.provenance.json']);
  let total = 0;
  const attachments = packet.attachments.map(file => {
    if (!file || typeof file !== 'object' || Object.keys(file).some(key => !['path', 'bytes'].includes(key)))
      throw new NativeWindError(400, 'Source attachments require only their path and exact bytes.');
    const path = file?.path;
    if (typeof path !== 'string' || Buffer.byteLength(path, 'utf8') > 512 || !path
        || path.startsWith('-') || path.includes('\\') || path.includes(':') || /[\u0000-\u0020\u007f]/.test(path)
        || path.split('/').some(part => !part || part === '.' || part === '..' || part.startsWith('.')))
      throw new NativeWindError(400, 'Attachment paths must be explicit relative paths without whitespace or traversal.');
    if (names.has(path) || [...names].some(name => path.startsWith(name + '/') || name.startsWith(path + '/')))
      throw new NativeWindError(400, 'Attachment paths collide with another input or the native output.');
    names.add(path);
    const content = decodedBase64(file.bytes, MAX_ATTACHMENT_FILE_BYTES);
    if (!content.length) throw new NativeWindError(400, 'Source attachments must be nonempty.');
    total += content.length;
    if (total > MAX_ATTACHMENT_BYTES) throw new NativeWindError(413, 'Attachments exceed the 32 MiB combined limit.');
    return {path, bytes:content};
  });
  return {source, attachments, decimate:packet.decimate ?? false};
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
