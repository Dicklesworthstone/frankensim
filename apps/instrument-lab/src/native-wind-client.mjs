// Transport only. Source parsing, material laws and synthesis belong to the selected native performer.
import {NATIVE_WAVE_LIMIT, PRESSURE_RECEIPT_HEADER, verifyPressureReceipt} from './native-wave.mjs';
import {MAX_HISTORY_BYTES, FINITE_GAS_STRING_RESULT_FORMAT,
  serializeNativeExperiment, validateInput} from './native-experiment.mjs';

const STATUS_URL = "/api/native-wind/status";
const RENDER_URL = "/api/native-wind/render";

function positiveInteger(value) {
  return Number.isSafeInteger(value) && value > 0;
}

export async function nativeWindStatus(signal, request = globalThis.fetch) {
  const response = await request(STATUS_URL, {
    signal, cache: "no-store", credentials: "same-origin",
    headers: { Accept: "application/json" },
  });
  if (!response.ok || !response.headers.get("content-type")?.includes("application/json")) {
    throw new Error("The native renderer is unavailable on this host.");
  }
  const status = await response.json();
  if (typeof status?.available !== "boolean" || typeof status?.busy !== "boolean"
    || !positiveInteger(status.maximumSourceBytes) || !positiveInteger(status.maximumOutputBytes)
    || !positiveInteger(status.timeoutMs)
    || (status.available && (typeof status.token !== "string" || status.token.length === 0))) {
    throw new Error("The native renderer returned an invalid status.");
  }
  return status;
}

async function refusal(response) {
  if (response.headers.get("content-type")?.includes("application/json")) {
    try {
      const body = await response.json();
      if (typeof body?.error === "string" && body.error.length > 0) return body.error;
    } catch { /* Keep the HTTP refusal when a proxy returns an unreadable response. */ }
  }
  return `The native renderer refused this request (HTTP ${response.status}).`;
}

async function boundedResponse(response, maximumBytes, label, type) {
  const declared = response.headers.get("content-length");
  if (declared !== null && (!/^\d+$/.test(declared) || Number(declared) > maximumBytes)) {
    await response.body?.cancel();
    throw new Error(`The native ${label} exceeds the renderer's output limit.`);
  }
  if (!response.body) throw new Error(`The native renderer returned no ${label}.`);
  const reader = response.body.getReader();
  const chunks = [];
  let bytes = 0;
  try {
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      bytes += value.byteLength;
      if (bytes > maximumBytes) {
        await reader.cancel();
        throw new Error(`The native ${label} exceeds the renderer's output limit.`);
      }
      chunks.push(value);
    }
  } finally {
    reader.releaseLock();
  }
  return new Blob(chunks, {type});
}

async function boundedWave(response, maximumBytes) {
  const wave = await boundedResponse(response, maximumBytes, 'WAV', 'audio/wav');
  const header = new DataView(await wave.slice(0, 12).arrayBuffer());
  if (wave.size < 44 || header.getUint32(0, false) !== 0x52494646
    || header.getUint32(8, false) !== 0x57415645 || header.getUint32(4, true) !== wave.size - 8) {
    throw new Error("The native renderer returned an invalid WAV.");
  }
  return wave;
}

function resultBytes(value, limit, label) {
  if (typeof value !== 'string' || !value.length || value.length > 4 * Math.ceil(limit / 3)
      || value.length % 4 !== 0 || !/^[A-Za-z0-9+/]*={0,2}$/.test(value))
    throw new Error(`Native ${label} is invalid or exceeds its output limit.`);
  let binary;
  try { binary = atob(value); }
  catch { throw new Error(`Native ${label} is not valid base64.`); }
  if (binary.length > limit || btoa(binary) !== value)
    throw new Error(`Native ${label} is invalid or exceeds its output limit.`);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return bytes;
}

async function finiteGasStringResult(response, status) {
  const limit = 4 * Math.ceil((status.maximumOutputBytes + status.maximumHistoryBytes) / 3) + 65536;
  const encoded = await boundedResponse(response, limit, 'result', 'application/json');
  let result;
  try { result = JSON.parse(new TextDecoder('utf-8', {fatal:true}).decode(await encoded.arrayBuffer())); }
  catch { throw new Error('The native finite-gas string result is not valid UTF-8 JSON.'); }
  if (!result || result.format !== FINITE_GAS_STRING_RESULT_FORMAT
      || Object.keys(result).some(key => !['format', 'wave', 'history', 'receipt'].includes(key)))
    throw new Error('The native finite-gas string result has an unsupported format.');
  const wave = resultBytes(result.wave, status.maximumOutputBytes, 'WAV');
  const history = resultBytes(result.history, status.maximumHistoryBytes, 'history');
  try { new TextDecoder('utf-8', {fatal:true}).decode(history); }
  catch { throw new Error('The native history is not valid UTF-8 CSV.'); }
  const receipt = await verifyPressureReceipt(result.receipt, wave, 'finite-gas-string');
  return {wave:new Blob([wave], {type:'audio/wav'}), receipt,
    history:new Blob([history], {type:'text/csv;charset=utf-8'})};
}

export async function renderNativeWind(source, status, signal, request = globalThis.fetch) {
  return (await sendNativeWind(source, [], status, {}, signal, request)).wave;
}

/** Exact attached source files; the native loader remains the material authority. */
export async function renderNativeWindBundle(source, attachments, status, options = {}, signal, request = globalThis.fetch) {
  return (await sendNativeWind(source, attachments, status, options, signal, request)).wave;
}

/** Original WAV plus verified native pressure calibration, when the executable emits it. */
export async function renderNativePerformance(source, attachments, status, options = {}, signal, request = globalThis.fetch) {
  return sendNativeWind(source, attachments, status, options, signal, request);
}

async function sendNativeWind(source, attachments, status, options, signal, request) {
  const renderer = options.renderer === undefined ? 'wind' : options.renderer;
  if (!['wind', 'strings', 'resolved', 'finite-gas-string'].includes(renderer)) throw new Error("Choose a supported native physical performer.");
  if (renderer === 'strings' && status.renderers?.strings !== true) {
    throw new Error("The material-string renderer is unavailable. Configure material_strings on the local server.");
  }
  if (renderer === 'wind' && status.renderers?.wind === false) {
    throw new Error("The wind renderer is unavailable. Configure music_render on the local server.");
  }
  if (renderer === 'resolved' && status.renderers?.resolved !== true) {
    throw new Error("The resolved-flow renderer is unavailable. Update music_render and the local server.");
  }
  if (renderer === 'finite-gas-string') {
    if (status.renderers?.['finite-gas-string'] !== true)
      throw new Error('The finite-gas string renderer is unavailable. Configure finite_gas_string on the local server.');
    if (!positiveInteger(status.maximumHistoryBytes) || status.maximumHistoryBytes > MAX_HISTORY_BYTES
        || !positiveInteger(status.maximumOutputBytes) || status.maximumOutputBytes > NATIVE_WAVE_LIMIT)
      throw new Error('The finite-gas string renderer returned invalid history or WAV output limits.');
  }
  if (options.stringProfiles !== undefined && (renderer !== 'strings' || !options.stringProfiles
      || typeof options.stringProfiles !== 'object' || Array.isArray(options.stringProfiles))) {
    throw new Error("Choose explicit attached profiles for the material-string renderer.");
  }
  if (options.coatingStrainGuard !== undefined) {
    if (renderer !== 'strings') throw new Error('Coating strain checks require the material-string performer.');
    if (status.supportsCoatingStrainGuard !== true)
      throw new Error('This local renderer cannot retain coating strain checks. Update the local server.');
  }
  if (options.materialManifests !== undefined && status.supportsMaterialManifests !== true)
    throw new Error('This local renderer cannot compile material manifests. Configure the xtask material compiler on the local server.');
  if (options.decimate !== undefined && typeof options.decimate !== 'boolean') {
    throw new Error("Choose whether to use native decimation explicitly.");
  }
  if (options.gasSources !== undefined) {
    if (!options.gasSources || typeof options.gasSources !== "object"
        || !Array.isArray(options.gasSources.components) || !options.gasSources.components.length) {
      throw new Error("Gas sources need explicit species, NASA-9 and transport-pack associations.");
    }
    if (status.supportsGasSources !== true) {
      throw new Error("This local renderer does not support attached gas sources. Update the local server.");
    }
  }
  if (!status.available) throw new Error("The native renderer is unavailable on this host.");
  if (typeof status.token !== "string" || status.token.length === 0) {
    throw new Error("Check the native renderer connection before rendering.");
  }
  if (typeof source !== "string" || source.trim().length === 0) {
    throw new Error("Enter a performance description before rendering.");
  }
  const bytes = new TextEncoder().encode(source);
  if (bytes.byteLength > status.maximumSourceBytes) {
    throw new Error("The performance description exceeds the renderer's source limit.");
  }
  let body = bytes, contentType = "text/plain;charset=utf-8";
  if (attachments.length || options.decimate || options.gasSources !== undefined || renderer !== 'wind'
      || options.stringProfiles !== undefined || options.coatingStrainGuard !== undefined || options.finiteGasString !== undefined || options.materialManifests !== undefined) {
    if (![status.maximumAttachments, status.maximumAttachmentFileBytes,
      status.maximumAttachmentBytes, status.maximumRequestBytes].every(positiveInteger)) {
      throw new Error("This local renderer does not support attached source files. Update the local server.");
    }
    const compiledCount = Array.isArray(options.materialManifests) ? options.materialManifests.length : 0;
    if (attachments.length + compiledCount > status.maximumAttachments) throw new Error("Too many source attachments and compiled packs.");
    let total = 0;
    for (const file of attachments) {
      if (typeof file.path !== "string" || !(file.bytes instanceof Uint8Array) || !file.bytes.byteLength
        || file.bytes.byteLength > status.maximumAttachmentFileBytes) {
        throw new Error("Each attachment needs its source path and nonempty bytes within the per-file limit.");
      }
      total += file.bytes.byteLength;
      if (total > status.maximumAttachmentBytes) throw new Error("Source attachments exceed the combined byte limit.");
    }
    body = serializeNativeExperiment(source, attachments, options);
    contentType = "application/json;charset=utf-8";
    if (new TextEncoder().encode(body).byteLength > status.maximumRequestBytes) {
      throw new Error("The complete experiment exceeds the renderer's request limit.");
    }
    if (renderer === 'finite-gas-string' || options.finiteGasString !== undefined || options.coatingStrainGuard !== undefined || options.materialManifests !== undefined)
      validateInput(new TextEncoder().encode(body), true);
  }
  const response = await request(RENDER_URL, {
    method: "POST", signal, cache: "no-store", credentials: "same-origin",
    headers: { "Content-Type": contentType, Accept: "audio/wav, application/json",
      "X-Native-Wind-Token": status.token },
    body,
  });
  if (!response.ok) throw new Error(await refusal(response));
  if (renderer === 'finite-gas-string') {
    if (response.headers.get('content-type')?.split(';')[0]?.trim() !== 'application/json') {
      await response.body?.cancel();
      throw new Error('The native finite-gas string renderer returned an unexpected response.');
    }
    return finiteGasStringResult(response, status);
  }
  if (response.headers.get("content-type")?.split(";")[0]?.trim() !== "audio/wav") {
    await response.body?.cancel();
    throw new Error("The native renderer returned an unexpected response instead of audio.");
  }
  const wave = await boundedWave(response, status.maximumOutputBytes);
  const header = response.headers.get(PRESSURE_RECEIPT_HEADER);
  let receipt = null;
  if (header !== null) {
    let value;
    try {
      if (header.length > 2048) throw new Error('oversized receipt');
      value = JSON.parse(header);
    } catch { throw new Error('The native renderer returned invalid pressure calibration.'); }
    receipt = await verifyPressureReceipt(value, new Uint8Array(await wave.arrayBuffer()), renderer);
  }
  return {wave, receipt};
}
