// The existing native request packet is also the portable experiment file.
// Browser editing and the loopback service share these byte/path bounds; only
// the Rust performer interprets physical records or material coefficients.
export const MAX_SOURCE_BYTES = 1024 * 1024;
export const MAX_ATTACHMENTS = 32;
export const MAX_ATTACHMENT_FILE_BYTES = 8 * 1024 * 1024;
export const MAX_ATTACHMENT_BYTES = 32 * 1024 * 1024;
export const MAX_REQUEST_BYTES = 4 * Math.ceil((MAX_SOURCE_BYTES + MAX_ATTACHMENT_BYTES) / 3) + 65536;
export const MAX_HISTORY_BYTES = 16 * 1024 * 1024;
export const FINITE_GAS_STRING_RESULT_FORMAT = 'frankensim-finite-gas-string-result-v1';
export const SUPPORTED_SCHEMAS = ['frankensim-reed-performance-v1', 'frankensim-flue-dg-performance-v1',
  'frankensim-material-modal-dg-valve-v1', 'frankensim-material-modal-valve-v1',
  'frankensim-plate-valve-performance-v1', 'frankensim-flue-performance-v1'];
export const STRING_SCHEMAS = ['frankensim-material-strings-v1', 'frankensim-material-strings-v2',
  'frankensim-material-strings-v3'];
export const RESOLVED_SCHEMAS = ['frankensim-resolved-dg-performance-v1'];
export const FINITE_GAS_STRING_SCHEMAS = ['frankensim-finite-gas-string-v1'];

export class NativeWindError extends Error {
  constructor(status, message) { super(message); this.status = status; }
}

const utf8Bytes = value => new TextEncoder().encode(value).byteLength;

export function validateSource(bytes, renderer = 'wind') {
  if (!['wind', 'strings', 'resolved', 'finite-gas-string'].includes(renderer)) throw new NativeWindError(400, 'Choose a supported native physical performer.');
  if (!(bytes instanceof Uint8Array) || !bytes.length) throw new NativeWindError(400, 'Choose or enter a performance file.');
  if (bytes.length > MAX_SOURCE_BYTES) throw new NativeWindError(413, 'Performance exceeds the 1 MiB input limit.');
  let text;
  try { text = new TextDecoder('utf-8', {fatal:true, ignoreBOM:true}).decode(bytes); }
  catch { throw new NativeWindError(400, 'Performance must contain valid UTF-8 bytes.'); }
  const schemas = renderer === 'strings' ? STRING_SCHEMAS : renderer === 'resolved' ? RESOLVED_SCHEMAS
    : renderer === 'finite-gas-string' ? FINITE_GAS_STRING_SCHEMAS : SUPPORTED_SCHEMAS;
  if (!schemas.includes(text.split(/\r?\n/, 1)[0]))
    throw new NativeWindError(400, renderer === 'strings'
      ? 'Choose an original material-string performance supported by material_strings.'
      : renderer === 'resolved'
        ? 'Choose an original resolved-flow performance supported by music_render resolved.'
        : renderer === 'finite-gas-string'
          ? 'Choose an original heated string scenario supported by finite_gas_string.'
        : 'Choose a native reed, flue or material-valve performance supported by music_render wind.');
  // Do not trim, normalize newlines, parse records, or re-encode admitted bytes.
  return bytes;
}

function base64(bytes) {
  let binary = '';
  for (let start = 0; start < bytes.length; start += 32768)
    binary += String.fromCharCode(...bytes.subarray(start, start + 32768));
  return btoa(binary);
}

function decodedBase64(value, limit) {
  if (typeof value !== 'string' || value.length > 4 * Math.ceil(limit / 3)
      || value.length % 4 !== 0 || !/^[A-Za-z0-9+/]*={0,2}$/.test(value))
    throw new NativeWindError(400, 'Source and attachments require canonical base64 bytes.');
  let bytes;
  try {
    const binary = atob(value);
    bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  }
  catch { throw new NativeWindError(400, 'Source or attachment encoding is invalid.'); }
  if (bytes.length > limit || base64(bytes) !== value)
    throw new NativeWindError(400, 'Source or attachment encoding exceeds its byte limit or is invalid.');
  return bytes;
}

function validateGasSources(value, attachedNames) {
  if (value === undefined) return undefined;
  if (!value || typeof value !== 'object' || Array.isArray(value)
      || Object.keys(value).some(key => !['components', 'waterSpecies', 'phase'].includes(key))
      || !Array.isArray(value.components) || !value.components.length || value.components.length > 128)
    throw new NativeWindError(400, 'Gas sources require 1..128 explicit species/model/transport associations.');
  if (!['frozen', 'liquid-water'].includes(value.phase))
    throw new NativeWindError(400, 'Choose an explicit frozen or liquid-water gas phase policy.');
  if (value.waterSpecies !== null && (typeof value.waterSpecies !== 'string'
      || !value.waterSpecies || value.waterSpecies === 'none' || utf8Bytes(value.waterSpecies) > 512
      || /[\u0000-\u0020\u007f]/.test(value.waterSpecies)))
    throw new NativeWindError(400, 'Gas water identity must be a canonical species ID or null for none.');
  if (value.phase === 'liquid-water' && value.waterSpecies === null)
    throw new NativeWindError(400, 'The liquid-water policy requires an explicit water species ID.');
  const fields = ['speciesPath', 'nasa9Path', 'transportPath'];
  const components = value.components.map(component => {
    if (!component || typeof component !== 'object' || Array.isArray(component)
        || Object.keys(component).some(key => !fields.includes(key))
        || fields.some(key => typeof component[key] !== 'string' || !attachedNames.has(component[key])))
      throw new NativeWindError(400, 'Every gas component needs three exact paths to attached source files.');
    return {speciesPath:component.speciesPath, nasa9Path:component.nasa9Path, transportPath:component.transportPath};
  });
  return {components, waterSpecies:value.waterSpecies, phase:value.phase};
}

function validateStringProfiles(value, renderer, attachedNames) {
  if (value === undefined) return undefined;
  const keys = ['bodyPath', 'stopPath', 'contactThermalPath'];
  if (renderer !== 'strings' || !value || typeof value !== 'object' || Array.isArray(value)
      || Object.keys(value).some(key => !keys.includes(key))
      || Object.values(value).some(path => typeof path !== 'string' || !attachedNames.has(path)))
    throw new NativeWindError(400, 'String profiles must name attached body, stop or contact-thermal source files.');
  return Object.fromEntries(keys.filter(key => value[key] !== undefined).map(key => [key, value[key]]));
}

function validateCoatingStrainGuard(value, renderer, profiles) {
  if (value === undefined) return undefined;
  const keys = ['maximumBytes', 'maximumWork'];
  if (renderer !== 'strings' || !profiles?.bodyPath || !value || typeof value !== 'object'
      || Array.isArray(value) || Object.keys(value).some(key => !keys.includes(key))
      || keys.some(key => !Number.isSafeInteger(value[key]) || value[key] <= 0))
    throw new NativeWindError(400, 'Coating strain checks require a material-string body profile and positive safe-integer byte/work limits.');
  return {maximumBytes:value.maximumBytes, maximumWork:value.maximumWork};
}

function validateFiniteGasString(value, renderer, attachedNames, gasSources, decimate) {
  if (renderer !== 'finite-gas-string') {
    if (value !== undefined) throw new NativeWindError(400, 'Finite-gas string sources require their own performer.');
    return undefined;
  }
  const keys = ['solidPath', 'gasRecordPath'];
  if (!value || typeof value !== 'object' || Array.isArray(value)
      || Object.keys(value).some(key => !keys.includes(key))
      || keys.some(key => typeof value[key] !== 'string' || !attachedNames.has(value[key])))
    throw new NativeWindError(400, 'Finite-gas strings require explicit attached solid-pack and initial-gas-record paths.');
  if (!gasSources || gasSources.phase !== 'frozen')
    throw new NativeWindError(400, 'Finite-gas strings require attached gas sources with the frozen phase policy.');
  if (decimate) throw new NativeWindError(400, 'Finite-gas strings require an authored 48 kHz clock; native decimation is unavailable.');
  return {solidPath:value.solidPath, gasRecordPath:value.gasRecordPath};
}

function admitAttachmentPath(path, names) {
  if (typeof path !== 'string' || utf8Bytes(path) > 512 || !path
      || path.startsWith('-') || path.includes('\\') || path.includes(':') || /[\u0000-\u0020\u007f]/.test(path)
      || path.split('/').some(part => !part || part === '.' || part === '..' || part.startsWith('.')))
    throw new NativeWindError(400, 'Attachment paths must be explicit relative paths without whitespace or traversal.');
  if (names.has(path) || [...names].some(name => path.startsWith(name + '/') || name.startsWith(path + '/')))
    throw new NativeWindError(400, 'Attachment paths collide with another input or the native output.');
  names.add(path);
  return path;
}

function validateMaterialManifests(value, attachments, names) {
  if (value === undefined) return undefined;
  if (!Array.isArray(value) || !value.length || attachments.length + value.length > MAX_ATTACHMENTS)
    throw new NativeWindError(400, 'Material manifests and attached sources together must fit the source-file limit.');
  const original = new Set(attachments.map(file => file.path));
  return value.map(binding => {
    if (!binding || typeof binding !== 'object' || Array.isArray(binding)
        || Object.keys(binding).some(key => !['packPath', 'manifestPath', 'materialState'].includes(key))
        || !original.has(binding.manifestPath))
      throw new NativeWindError(400, 'Each material compilation needs an original attached manifest and a new pack path.');
    const packPath = admitAttachmentPath(binding.packPath, names);
    const state = binding.materialState;
    if (state !== undefined && (!state || typeof state !== 'object' || Array.isArray(state)
        || Object.keys(state).some(key => !['chemistry', 'phase', 'process', 'revision'].includes(key))
        || ['chemistry', 'phase', 'process'].some(key => typeof state[key] !== 'string'
          || !state[key].trim() || utf8Bytes(state[key]) > 512 || /[\u0000-\u001f\u007f]/.test(state[key]))
        || !Number.isInteger(state.revision) || state.revision < 0 || state.revision > 4294967295))
      throw new NativeWindError(400, 'A named material card needs explicit chemistry, phase, process and an unsigned 32-bit revision.');
    return {packPath, manifestPath:binding.manifestPath,
      ...(state === undefined ? {} : {materialState:{chemistry:state.chemistry, phase:state.phase,
        process:state.process, revision:state.revision}})};
  });
}

/** Decode an attached text source for editing; binary packs stay byte attachments. */
export function editableAttachmentText(bytes) {
  if (!(bytes instanceof Uint8Array) || !bytes.length)
    throw new NativeWindError(400, 'Choose a nonempty text source to edit.');
  if (bytes.length > MAX_ATTACHMENT_FILE_BYTES)
    throw new NativeWindError(413, 'Text source exceeds the 8 MiB attachment limit.');
  let text;
  try { text = new TextDecoder('utf-8', {fatal:true, ignoreBOM:true}).decode(bytes); }
  catch { throw new NativeWindError(400, 'This attachment is not UTF-8 text. Use Replace for binary source packs.'); }
  if (/[\u0000-\u0008\u000b\u000c\u000e-\u001f\u007f]/.test(text))
    throw new NativeWindError(400, 'This attachment contains binary controls. Use Replace for binary source packs.');
  return text;
}

/** Preserve original bytes when a textarea has only normalized its line endings. */
export function editedAttachmentBytes(original, text) {
  const before = editableAttachmentText(original);
  if (typeof text !== 'string')
    throw new NativeWindError(400, 'The edited source must be text.');
  if (text === before || text === before.replace(/\r\n?/g, '\n')) return original;
  const bytes = new TextEncoder().encode(text);
  editableAttachmentText(bytes);
  return bytes;
}

/** Transport envelope only: no physics records or material constants are interpreted. */
export function validateInput(bytes, bundled = false) {
  if (!bundled) return {source:validateSource(bytes), attachments:[], decimate:false};
  if (!(bytes instanceof Uint8Array)) throw new NativeWindError(400, 'Experiment must be UTF-8 JSON bytes.');
  if (bytes.length > MAX_REQUEST_BYTES) throw new NativeWindError(413, 'Experiment exceeds the complete request limit.');
  let packet;
  try { packet = JSON.parse(new TextDecoder('utf-8', {fatal:true}).decode(bytes)); }
  catch { throw new NativeWindError(400, 'Experiment must be UTF-8 JSON.'); }
  if (packet?.format !== 'frankensim-native-wind-input-v1' || !Array.isArray(packet.attachments)
      || packet.attachments.length > MAX_ATTACHMENTS)
    throw new NativeWindError(400, 'Experiment requires its version and at most 32 source attachments.');
  if (Object.keys(packet).some(key => !['format', 'source', 'attachments', 'decimate', 'gasSources', 'renderer', 'stringProfiles', 'coatingStrainGuard', 'finiteGasString', 'materialManifests'].includes(key))
      || (packet.decimate !== undefined && typeof packet.decimate !== 'boolean'))
    throw new NativeWindError(400, 'Experiment contains unsupported fields or an invalid native decimation option.');
  const renderer = packet.renderer === undefined ? 'wind' : packet.renderer;
  const source = validateSource(decodedBase64(packet.source, MAX_SOURCE_BYTES), renderer);
  const names = new Set(['input.performance', 'output.wav', 'output.provenance.json']);
  if (renderer === 'finite-gas-string') names.add('output.history.csv');
  let total = 0;
  const attachments = packet.attachments.map(file => {
    if (!file || typeof file !== 'object' || Object.keys(file).some(key => !['path', 'bytes'].includes(key)))
      throw new NativeWindError(400, 'Source attachments require only their path and exact bytes.');
    const path = admitAttachmentPath(file?.path, names);
    const content = decodedBase64(file.bytes, MAX_ATTACHMENT_FILE_BYTES);
    if (!content.length) throw new NativeWindError(400, 'Source attachments must be nonempty.');
    total += content.length;
    if (total > MAX_ATTACHMENT_BYTES) throw new NativeWindError(413, 'Attachments exceed the 32 MiB combined limit.');
    return {path, bytes:content};
  });
  const materialManifests = validateMaterialManifests(packet.materialManifests, attachments, names);
  const attachedNames = new Set([...attachments.map(file => file.path),
    ...(materialManifests ?? []).map(binding => binding.packPath)]);
  const gasSources = validateGasSources(packet.gasSources, attachedNames);
  const stringProfiles = validateStringProfiles(packet.stringProfiles, renderer, attachedNames);
  const coatingStrainGuard = validateCoatingStrainGuard(packet.coatingStrainGuard, renderer, stringProfiles);
  const finiteGasString = validateFiniteGasString(packet.finiteGasString, renderer, attachedNames, gasSources, packet.decimate);
  return {source, attachments, decimate:packet.decimate ?? false, ...(renderer !== 'wind' ? {renderer} : {}),
    ...(gasSources ? {gasSources} : {}), ...(stringProfiles ? {stringProfiles} : {}),
    ...(coatingStrainGuard ? {coatingStrainGuard} : {}),
    ...(finiteGasString ? {finiteGasString} : {}),
    ...(materialManifests ? {materialManifests} : {})};
}

/** Same versioned packet for HTTP rendering and files. No token, executable or output path is stored. */
export function serializeNativeExperiment(source, attachments, options = {}) {
  return JSON.stringify({format:'frankensim-native-wind-input-v1',
    source:base64(new TextEncoder().encode(source)),
    attachments:attachments.map(file => ({path:file.path, bytes:base64(file.bytes)})),
    decimate:options.decimate ?? false,
    ...(options.renderer !== undefined ? {renderer:options.renderer} : {}),
    ...(options.gasSources !== undefined ? {gasSources:options.gasSources} : {}),
    ...(options.stringProfiles !== undefined ? {stringProfiles:options.stringProfiles} : {}),
    ...(options.coatingStrainGuard !== undefined ? {coatingStrainGuard:options.coatingStrainGuard} : {}),
    ...(options.finiteGasString !== undefined ? {finiteGasString:options.finiteGasString} : {}),
    ...(options.materialManifests !== undefined ? {materialManifests:options.materialManifests} : {})});
}
