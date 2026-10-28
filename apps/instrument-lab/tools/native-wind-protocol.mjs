// Transport limits only. The selected native performer owns the source grammar and physics.
import {inspectNativeWave, NATIVE_WAVE_LIMIT} from '../src/native-wave.mjs';
import {NativeWindError} from '../src/native-experiment.mjs';
export {MAX_SOURCE_BYTES, MAX_ATTACHMENTS, MAX_ATTACHMENT_FILE_BYTES, MAX_ATTACHMENT_BYTES,
  MAX_REQUEST_BYTES, SUPPORTED_SCHEMAS, STRING_SCHEMAS, NativeWindError,
  validateSource, validateInput} from '../src/native-experiment.mjs';

export const MAX_OUTPUT_BYTES = NATIVE_WAVE_LIMIT;
export const MAX_LOG_BYTES = 256 * 1024;
export const DEFAULT_TIMEOUT_MS = 60000;

export function validateWave(bytes, limit = MAX_OUTPUT_BYTES) {
  if (bytes.length > limit) throw new NativeWindError(502, 'Native WAV exceeds the output limit.');
  try { inspectNativeWave(bytes, limit); }
  catch (error) { throw new NativeWindError(502, error.message); }
  return bytes;
}
