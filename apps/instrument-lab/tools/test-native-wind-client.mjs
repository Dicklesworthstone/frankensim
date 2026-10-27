import assert from "node:assert/strict";
import test from "node:test";
import { nativeWindStatus, renderNativeWind } from "../src/native-wind-client.mjs";

const status = { available: true, busy: false, token: "test-session",
  maximumSourceBytes: 1_048_576, maximumOutputBytes: 16_777_216, timeoutMs: 60_000 };
const json = (value, code = 200) => new Response(JSON.stringify(value), {
  status: code, headers: { "Content-Type": "application/json" },
});
function wave() {
  const bytes = new Uint8Array(48);
  const view = new DataView(bytes.buffer);
  bytes.set(new TextEncoder().encode("RIFF"), 0);
  view.setUint32(4, bytes.length - 8, true);
  bytes.set(new TextEncoder().encode("WAVEfmt "), 8);
  view.setUint32(16, 16, true);
  view.setUint16(20, 1, true);
  view.setUint16(22, 1, true);
  view.setUint32(24, 48_000, true);
  view.setUint32(28, 96_000, true);
  view.setUint16(32, 2, true);
  view.setUint16(34, 16, true);
  bytes.set(new TextEncoder().encode("data"), 36);
  view.setUint32(40, 4, true);
  return bytes;
}
const audio = (bytes = wave(), headers = {}) => new Response(bytes, {
  headers: { "Content-Type": "audio/wav", ...headers },
});

test("status uses the same-origin uncached endpoint and keeps the server's limits/token", async () => {
  const controller = new AbortController();
  const result = await nativeWindStatus(controller.signal, async (url, options) => {
    assert.equal(url, "/api/native-wind/status");
    assert.equal(options.cache, "no-store");
    assert.equal(options.credentials, "same-origin");
    assert.equal(options.signal, controller.signal);
    return json(status);
  });
  assert.deepEqual(result, status);
});

test("a static host's HTML fallback is unavailable, not a native service", async () => {
  await assert.rejects(nativeWindStatus(undefined, async () => new Response("<!doctype html>", {
    headers: { "Content-Type": "text/html" },
  })), /unavailable/);
  await assert.rejects(nativeWindStatus(undefined, async () => json({ available: true })), /invalid status/);
  await assert.rejects(nativeWindStatus(undefined, async () => json({ ...status, token: "" })), /invalid status/);
  const unavailable = { ...status, available: false };
  delete unavailable.token;
  assert.equal((await nativeWindStatus(undefined, async () => json(unavailable))).available, false);
});

test("source reaches the native endpoint as exact UTF-8, with its session token", async () => {
  const source = "frankensim-reed-performance-v1\r\n# µ · brass\ntrack\tblow\n";
  const expected = wave();
  const result = await renderNativeWind(source, status, undefined, async (url, options) => {
    assert.equal(url, "/api/native-wind/render");
    assert.equal(options.method, "POST");
    assert.equal(options.credentials, "same-origin");
    assert.equal(options.headers["X-Native-Wind-Token"], "test-session");
    assert.equal(options.headers["Content-Type"], "text/plain;charset=utf-8");
    assert.deepEqual(options.body, new TextEncoder().encode(source));
    return audio(expected);
  });
  assert.equal(result.type, "audio/wav");
  assert.deepEqual(new Uint8Array(await result.arrayBuffer()), expected);
});

test("the source limit counts UTF-8 bytes and does not send truncated or empty descriptions", async () => {
  const never = async () => { assert.fail("invalid source must not be sent"); };
  await assert.rejects(renderNativeWind("µµµ", { ...status, maximumSourceBytes: 5 }, undefined, never), /source limit/);
  await assert.rejects(renderNativeWind(" \n\t", status, undefined, never), /Enter a performance/);
  await assert.rejects(renderNativeWind("source", { ...status, token: "" }, undefined, never), /connection/);
});

test("native source refusals retain the actual diagnostic and never become preview audio", async () => {
  const diagnostic = "FS-DUCT-BAD-PARAMETER: duct and exterior need a common static pressure";
  await assert.rejects(renderNativeWind("source", status, undefined, async () => json({ error: diagnostic }, 400)),
    { message: diagnostic });
  await assert.rejects(renderNativeWind("source", status, undefined, async () => json({ error: "Renderer busy" }, 409)),
    { message: "Renderer busy" });
  await assert.rejects(renderNativeWind("source", status, undefined, async () => new Response("proxy error", { status: 502 })),
    /HTTP 502/);
});

test("success requires the native WAV media type and a complete RIFF container", async () => {
  await assert.rejects(renderNativeWind("source", status, undefined, async () => new Response(wave(), {
    headers: { "Content-Type": "text/html" },
  })), /unexpected response/);
  const invalid = wave();
  invalid[0] = 0;
  await assert.rejects(renderNativeWind("source", status, undefined, async () => audio(invalid)), /invalid WAV/);
  await assert.rejects(renderNativeWind("source", status, undefined, async () => audio(new Uint8Array(2))), /invalid WAV/);
  const truncated = wave().slice(0, 46);
  await assert.rejects(renderNativeWind("source", status, undefined, async () => audio(truncated)), /invalid WAV/);
});

test("both advertised and streamed oversized outputs are refused before preview allocation", async () => {
  await assert.rejects(renderNativeWind("source", { ...status, maximumOutputBytes: 47 }, undefined,
    async () => audio(wave(), { "Content-Length": "48" })), /output limit/);
  let cancelled = false;
  const body = new ReadableStream({
    start(controller) { controller.enqueue(new Uint8Array(65)); },
    cancel() { cancelled = true; },
  });
  await assert.rejects(renderNativeWind("source", { ...status, maximumOutputBytes: 64 }, undefined,
    async () => new Response(body, { headers: { "Content-Type": "audio/wav" } })), /output limit/);
  assert.equal(cancelled, true);
});

test("cancellation reaches the in-flight native request", async () => {
  const controller = new AbortController();
  const pending = renderNativeWind("source", status, controller.signal, async (_url, options) => {
    assert.equal(options.signal, controller.signal);
    return new Promise((_resolve, reject) => options.signal.addEventListener("abort", () => reject(options.signal.reason), { once: true }));
  });
  controller.abort();
  await assert.rejects(pending, { name: "AbortError" });
});
