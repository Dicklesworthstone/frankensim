import assert from "node:assert/strict";
import test from "node:test";
import { readFile } from "node:fs/promises";
import { nativeWindStatus, renderNativeWind, renderNativeWindBundle } from "../src/native-wind-client.mjs";
import { validateInput } from "./native-wind-protocol.mjs";

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

test("binary source bundles retain exact material bytes, paths, performance and WAV amplitude", async () => {
  const extended = {...status, maximumAttachments:32, maximumAttachmentFileBytes:8*1024*1024,
    maximumAttachmentBytes:32*1024*1024, maximumRequestBytes:48*1024*1024};
  const source = "frankensim-material-modal-dg-valve-v1\n# µ native data\n";
  const bytes = Uint8Array.from({length:130001}, (_, i) => i % 256);
  const expected = wave(); expected[44] = 0x9a; expected[45] = 0x37;
  const result = await renderNativeWindBundle(source, [{path:"packs/reed.fsmcdpk", bytes}], extended, {}, undefined, async (url, options) => {
    assert.equal(url, "/api/native-wind/render");
    assert.equal(options.headers["Content-Type"], "application/json;charset=utf-8");
    assert.equal(options.headers["X-Native-Wind-Token"], status.token);
    const packet = JSON.parse(options.body);
    assert.equal(packet.format, "frankensim-native-wind-input-v1");
    assert.deepEqual(Buffer.from(packet.source, "base64"), Buffer.from(source));
    assert.equal(packet.attachments[0].path, "packs/reed.fsmcdpk");
    assert.deepEqual(Buffer.from(packet.attachments[0].bytes, "base64"), Buffer.from(bytes));
    return audio(expected);
  });
  assert.deepEqual(new Uint8Array(await result.arrayBuffer()), expected);
});

test("attachment limits and older servers refuse without discarding a selected material", async () => {
  const attachment = {path:"wall.fsmcdpk", bytes:new Uint8Array([1, 2, 3])};
  const never = async () => assert.fail("invalid attachment request must not be sent");
  await assert.rejects(renderNativeWindBundle("source", [attachment], status, {}, undefined, never), /does not support/);
  const extended = {...status, maximumAttachments:1, maximumAttachmentFileBytes:4, maximumAttachmentBytes:4,
    maximumRequestBytes:1024};
  await assert.rejects(renderNativeWindBundle("source", [attachment, attachment], extended, {}, undefined, never), /Too many/);
  await assert.rejects(renderNativeWindBundle("source", [{...attachment, bytes:new Uint8Array(5)}], extended, {}, undefined, never), /per-file limit/);
  await assert.rejects(renderNativeWindBundle("source", [attachment], {...extended, maximumAttachmentBytes:2}, {}, undefined, never), /combined byte limit/);
  await assert.rejects(renderNativeWindBundle("source", [attachment], {...extended, maximumRequestBytes:8}, {}, undefined, never), /request limit/);
});

test("higher mechanical clocks request the original native decimator without modifying source or WAV", async () => {
  const connected = {...status, maximumAttachments:32, maximumAttachmentFileBytes:8388608,
    maximumAttachmentBytes:33554432, maximumRequestBytes:48000000};
  const source = "frankensim-flue-dg-performance-v1\noriginal high-clock source";
  await renderNativeWindBundle(source, [], connected, {decimate:true}, undefined, async (_url, options) => {
    const packet = JSON.parse(options.body);
    assert.equal(packet.decimate, true);
    assert.deepEqual(packet.attachments, []);
    assert.equal(Buffer.from(packet.source, "base64").toString(), source);
    return audio();
  });
});

test("canonical DG flute and reed atmosphere records reach the native packet unchanged", async () => {
  const connected = {...status, maximumAttachments:32, maximumAttachmentFileBytes:8388608,
    maximumAttachmentBytes:33554432, maximumRequestBytes:48000000};
  for (const path of ["flue-dg-baffled.performance", "flue-dg-molecular.performance", "material_modal_valve/paired-dg-chimney.performance"]) {
    const bytes = await readFile(new URL(`../../../crates/fs-couple/examples/${path}`, import.meta.url));
    const source = bytes.toString("utf8");
    if (path === "flue-dg-molecular.performance") {
      assert.match(source, /^dg_molecular iso9613-terrestrial-air$/m);
      assert.match(source, /^observer_gas molecular \d+ \d+ \d+ iso9613-terrestrial-air$/m);
    } else assert.match(source, /^observer_gas stokes-kirchhoff \d+ \d+ \d+$/m);
    assert.match(source, /^(?:clock 96000|audio 96000 96000 \d+ \S+)$/m);
    await renderNativeWindBundle(source, [], connected, {decimate:true}, undefined, async (_url, options) => {
      const packet = validateInput(Buffer.from(options.body), true);
      assert.deepEqual(Buffer.from(packet.source), bytes);
      assert.equal(packet.renderer ?? "wind", "wind");
      assert.equal(packet.decimate, true);
      return audio();
    });
  }
});

test("gas source associations accompany exact source and pack bytes without browser coefficients", async () => {
  const connected = {...status, supportsGasSources:true, maximumAttachments:32, maximumAttachmentFileBytes:8388608,
    maximumAttachmentBytes:33554432, maximumRequestBytes:48000000};
  const gasSources = {components:[
    {speciesPath:"gas/custom.species", nasa9Path:"gas/custom.nasa9", transportPath:"gas/custom.transport"},
    {speciesPath:"gas/second.species", nasa9Path:"gas/second.nasa9", transportPath:"gas/second.transport"},
  ], waterSpecies:null, phase:"frozen"};
  const attachments = gasSources.components.flatMap(component => Object.values(component))
    .map((path, index) => ({path, bytes:new Uint8Array([0, 255, index, 13, 10])}));
  const source = "frankensim-flue-performance-v1\r\n# authored thermochemical record remains verbatim\r\n";
  await renderNativeWindBundle(source, attachments, connected, {gasSources}, undefined, async (_url, options) => {
    const packet = JSON.parse(options.body);
    assert.deepEqual(packet.gasSources, gasSources);
    assert.equal(Buffer.from(packet.source, "base64").toString(), source);
    for (let index = 0; index < attachments.length; index++) {
      assert.equal(packet.attachments[index].path, attachments[index].path);
      assert.deepEqual(Buffer.from(packet.attachments[index].bytes, "base64"), Buffer.from(attachments[index].bytes));
    }
    return audio();
  });
});

test("selected gas sources cannot silently disappear on an older server or malformed option", async () => {
  const connected = {...status, maximumAttachments:32, maximumAttachmentFileBytes:8388608,
    maximumAttachmentBytes:33554432, maximumRequestBytes:48000000};
  const gasSources = {components:[{speciesPath:"s", nasa9Path:"n", transportPath:"t"}], waterSpecies:null, phase:"frozen"};
  const never = async () => assert.fail("unsupported gas source request must not be sent");
  await assert.rejects(renderNativeWindBundle("source", [], connected, {gasSources}, undefined, never), /does not support attached gas sources/);
  for (const gasSources of [null, {}, {components:[]}]) {
    await assert.rejects(renderNativeWindBundle("source", [], {...connected, supportsGasSources:true},
      {gasSources}, undefined, never), /Gas sources need explicit/);
  }
});

test("material strings preserve their performer, source, profiles, gas associations and native PCM", async () => {
  const connected = {...status, renderers:{wind:false, strings:true}, supportsGasSources:true,
    maximumAttachments:32, maximumAttachmentFileBytes:8388608,
    maximumAttachmentBytes:33554432, maximumRequestBytes:48000000};
  const source = "frankensim-material-strings-v3\r\n# µ original string material data\r\n";
  const stringProfiles = {bodyPath:"profiles/body.profile", stopPath:"profiles/fingers.profile",
    contactThermalPath:"profiles/contact-heat.profile"};
  const gasSources = {components:[{speciesPath:"gas/custom.species", nasa9Path:"gas/custom.nasa9",
    transportPath:"gas/custom.transport"}], waterSpecies:null, phase:"frozen"};
  const attachments = [...Object.values(stringProfiles), ...Object.values(gasSources.components[0]), "materials/wire.fsmcdpk"]
    .map((path, index) => ({path, bytes:new Uint8Array([0, 255, index, 13, 10])}));
  const expected = wave(); expected[44] = 0x21; expected[45] = 0xa7;
  const result = await renderNativeWindBundle(source, attachments, connected,
    {renderer:"strings", stringProfiles, gasSources}, undefined, async (url, options) => {
      assert.equal(url, "/api/native-wind/render");
      assert.equal(options.headers["X-Native-Wind-Token"], status.token);
      assert.equal(options.headers["Content-Type"], "application/json;charset=utf-8");
      const packet = JSON.parse(options.body);
      assert.equal(packet.renderer, "strings");
      assert.equal(packet.decimate, false);
      assert.deepEqual(packet.stringProfiles, stringProfiles);
      assert.deepEqual(packet.gasSources, gasSources);
      assert.deepEqual(Buffer.from(packet.source, "base64"), Buffer.from(source));
      assert.deepEqual(packet.attachments.map(file => ({path:file.path, bytes:Buffer.from(file.bytes, "base64")})),
        attachments.map(file => ({path:file.path, bytes:Buffer.from(file.bytes)})));
      return audio(expected);
    });
  assert.deepEqual(new Uint8Array(await result.arrayBuffer()), expected);
  // Even a source with no external files needs an explicit performer; plain text selects wind.
  await renderNativeWindBundle(source, [], connected, {renderer:"strings", decimate:true}, undefined,
    async (_url, options) => {
      assert.equal(options.headers["Content-Type"], "application/json;charset=utf-8");
      const packet = JSON.parse(options.body);
      assert.equal(packet.renderer, "strings");
      assert.equal(packet.decimate, true);
      assert.deepEqual(packet.attachments, []);
      assert.equal(Buffer.from(packet.source, "base64").toString(), source);
      return audio();
    });
});

test("a selected performer or string profile cannot fall back to another native engine", async () => {
  const connected = {...status, renderers:{wind:true, strings:true}, maximumAttachments:32,
    maximumAttachmentFileBytes:8388608, maximumAttachmentBytes:33554432, maximumRequestBytes:48000000};
  const never = async () => assert.fail("unsupported performer must be refused before sending");
  for (const renderers of [undefined, {wind:true, strings:false}]) {
    await assert.rejects(renderNativeWindBundle("source", [], {...connected, renderers},
      {renderer:"strings"}, undefined, never), /material-string renderer is unavailable/);
  }
  await assert.rejects(renderNativeWind("source", {...connected, renderers:{wind:false, strings:true}},
    undefined, never), /wind renderer is unavailable/);
  for (const renderer of [null, "piano"]) {
    await assert.rejects(renderNativeWindBundle("source", [], connected, {renderer}, undefined, never),
      /Choose the wind, material-string or resolved-flow renderer/);
  }
  await assert.rejects(renderNativeWindBundle("source", [], connected,
    {stringProfiles:{bodyPath:"body.profile"}}, undefined, never), /material-string renderer/);
  for (const stringProfiles of [null, [], "body.profile"]) {
    await assert.rejects(renderNativeWindBundle("source", [], connected,
      {renderer:"strings", stringProfiles}, undefined, never), /explicit attached profiles/);
  }
  await assert.rejects(renderNativeWindBundle("source", [], connected,
    {renderer:"strings", decimate:1}, undefined, never), /decimation explicitly/);
});

test("resolved sources always select their command, including experiments without attachments", async () => {
  const connected = {...status, renderers:{wind:true, strings:false, resolved:true},
    maximumAttachments:32, maximumAttachmentFileBytes:8388608,
    maximumAttachmentBytes:33554432, maximumRequestBytes:48000000};
  const source = (await readFile(new URL("../../../crates/fs-couple/examples/resolved-dg.performance", import.meta.url))).toString("utf8");
  for (const decimate of [false, true]) {
    await renderNativeWindBundle(source, [], connected, {renderer:"resolved", decimate}, undefined,
      async (_url, options) => {
        assert.equal(options.headers["Content-Type"], "application/json;charset=utf-8");
        const packet = validateInput(new TextEncoder().encode(options.body), true);
        assert.equal(packet.renderer, "resolved");
        assert.equal(packet.decimate, decimate);
        assert.deepEqual(packet.attachments, []);
        assert.deepEqual(packet.source, new TextEncoder().encode(source));
        return audio();
      });
  }
  const never = async () => assert.fail("a missing resolved route must not send a wind request");
  for (const renderers of [undefined, {wind:true, strings:true}, {wind:true, strings:false, resolved:false}]) {
    await assert.rejects(renderNativeWindBundle(source, [], {...connected, renderers},
      {renderer:"resolved"}, undefined, never), /resolved-flow renderer is unavailable/);
  }
});
