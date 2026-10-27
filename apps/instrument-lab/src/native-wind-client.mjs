// Transport only. All source parsing, material laws and synthesis belong to music_render.
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
  return sendNativeWind(source, [], status, {}, signal, request);
}

/** Exact attached source files; the native loader remains the material authority. */
export async function renderNativeWindBundle(source, attachments, status, options = {}, signal, request = globalThis.fetch) {
  return sendNativeWind(source, attachments, status, options, signal, request);
}

async function sendNativeWind(source, attachments, status, options, signal, request) {
  if (options.decimate !== undefined && typeof options.decimate !== 'boolean') {
    throw new Error("Choose whether to use native decimation explicitly.");
  }
  if (options.gasSources !== undefined) {
    if (!options.gasSources || typeof options.gasSources !== "object"
        || !Array.isArray(options.gasSources.components) || !options.gasSources.components.length) {
      throw new Error("Gas sources need explicit species, NASA-9 and transport-claim associations.");
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
  if (attachments.length || options.decimate || options.gasSources !== undefined) {
    if (![status.maximumAttachments, status.maximumAttachmentFileBytes,
      status.maximumAttachmentBytes, status.maximumRequestBytes].every(positiveInteger)) {
      throw new Error("This local renderer does not support attached source files. Update the local server.");
    }
    if (attachments.length > status.maximumAttachments) throw new Error("Too many source attachments.");
    let total = 0;
    const files = attachments.map(file => {
      if (typeof file.path !== "string" || !(file.bytes instanceof Uint8Array) || !file.bytes.byteLength
        || file.bytes.byteLength > status.maximumAttachmentFileBytes) {
        throw new Error("Each attachment needs its source path and nonempty bytes within the per-file limit.");
      }
      total += file.bytes.byteLength;
      if (total > status.maximumAttachmentBytes) throw new Error("Source attachments exceed the combined byte limit.");
      return {path:file.path, bytes:base64(file.bytes)};
    });
    body = JSON.stringify({format:"frankensim-native-wind-input-v1", source:base64(bytes), attachments:files, decimate:options.decimate ?? false,
      ...(options.gasSources !== undefined ? {gasSources:options.gasSources} : {})});
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
  return boundedWave(response, status.maximumOutputBytes);
}

// Keep encoding bounded without spreading a large binary array on the call stack.
function base64(bytes) {
  let binary = "";
  for (let start = 0; start < bytes.length; start += 32768) {
    binary += String.fromCharCode(...bytes.subarray(start, start + 32768));
  }
  return btoa(binary);
}
