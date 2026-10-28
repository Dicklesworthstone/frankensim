// Transport only. Source parsing, material laws and synthesis belong to the selected native performer.
import {PRESSURE_RECEIPT_HEADER, verifyPressureReceipt} from './native-wave.mjs';
import {serializeNativeExperiment} from './native-experiment.mjs';

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

async function boundedWave(response, maximumBytes) {
  const declared = response.headers.get("content-length");
  if (declared !== null && (!/^\d+$/.test(declared) || Number(declared) > maximumBytes)) {
    await response.body?.cancel();
    throw new Error("The native WAV exceeds the renderer's output limit.");
  }
  if (!response.body) throw new Error("The native renderer returned no WAV.");
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
        throw new Error("The native WAV exceeds the renderer's output limit.");
      }
      chunks.push(value);
    }
  } finally {
    reader.releaseLock();
  }
  const wave = new Blob(chunks, { type: "audio/wav" });
  const header = new DataView(await wave.slice(0, 12).arrayBuffer());
  if (bytes < 44 || header.getUint32(0, false) !== 0x52494646
    || header.getUint32(8, false) !== 0x57415645 || header.getUint32(4, true) !== bytes - 8) {
    throw new Error("The native renderer returned an invalid WAV.");
  }
  return wave;
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
  if (!['wind', 'strings', 'resolved'].includes(renderer)) throw new Error("Choose the wind, material-string or resolved-flow renderer.");
  if (renderer === 'strings' && status.renderers?.strings !== true) {
    throw new Error("The material-string renderer is unavailable. Configure material_strings on the local server.");
  }
  if (renderer === 'wind' && status.renderers?.wind === false) {
    throw new Error("The wind renderer is unavailable. Configure music_render on the local server.");
  }
  if (renderer === 'resolved' && status.renderers?.resolved !== true) {
    throw new Error("The resolved-flow renderer is unavailable. Update music_render and the local server.");
  }
  if (options.stringProfiles !== undefined && (renderer !== 'strings' || !options.stringProfiles
      || typeof options.stringProfiles !== 'object' || Array.isArray(options.stringProfiles))) {
    throw new Error("Choose explicit attached profiles for the material-string renderer.");
  }
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
  if (attachments.length || options.decimate || options.gasSources !== undefined || renderer !== 'wind' || options.stringProfiles !== undefined) {
    if (![status.maximumAttachments, status.maximumAttachmentFileBytes,
      status.maximumAttachmentBytes, status.maximumRequestBytes].every(positiveInteger)) {
      throw new Error("This local renderer does not support attached source files. Update the local server.");
    }
    if (attachments.length > status.maximumAttachments) throw new Error("Too many source attachments.");
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
  }
  const response = await request(RENDER_URL, {
    method: "POST", signal, cache: "no-store", credentials: "same-origin",
    headers: { "Content-Type": contentType, Accept: "audio/wav, application/json",
      "X-Native-Wind-Token": status.token },
    body,
  });
  if (!response.ok) throw new Error(await refusal(response));
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
