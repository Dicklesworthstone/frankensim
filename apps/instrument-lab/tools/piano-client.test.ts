// bun test tools/piano-client.test.ts
// Real client lifecycle and snapshot parsing; controlled Web Audio/Worker boundaries.
// These checks do not qualify audio fidelity.
import { afterEach, test } from "bun:test";
import assert from "node:assert/strict";
const globals = { Worker: globalThis.Worker, AudioWorkletNode: globalThis.AudioWorkletNode };
afterEach(() => Object.assign(globalThis, globals));
function deferred() {
  let resolve!: () => void;
  let reject!: (error: Error) => void;
  const promise = new Promise<void>((yes, no) => { resolve = yes; reject = no; });
  return { promise, resolve, reject };
}
const turn = () => new Promise<void>(resolve => setTimeout(resolve, 0));

for (const header of [18, 19]) {
test(`piano snapshot readers preserve board, key geometry and bank exchange: header ${header}`, async () => {
  const { parseSnapshot } = await import("../public/live/piano/client.js");
  const { parsePianoSnapshot } = await import("../src/live/pianoClient");
  const data = new Float32Array(header + 2 + 2 * (9 + 3));
  data.set([2, 3, 2, 0.75, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 20, 6, 2.5, 1]);
  if (header === 19) data[18] = -0.125;
  data.set([0.01, 0.02], header);
  data.set([57, 1, 0, -0.005, 12, 0.125, 8, 3, 0.001, 0.1, 0.2, 0.3], header + 2);
  data.set([72, 0, 1, -0.02, 0, 0.25, 4, 2, 0.002, 0.4, 0.5, 0.6], header + 14);
  const plain = parseSnapshot(data), typed = parsePianoSnapshot(data);
  assert.deepEqual(plain.keys.map((k: { midi: number }) => k.midi), [57, 72]);
  assert.deepEqual(typed.keys.map(k => k.midi), [57, 72]);
  assert.deepEqual(plain.keys, typed.keys);
  assert.deepEqual(plain.boardCoordinates, data.subarray(header, header + 2));
  assert.deepEqual(typed.boardCoordinates, plain.boardCoordinates);
  assert.deepEqual(plain.ledger, typed.ledger);
  assert.equal(typed.ledger.couplingJ, header === 19 ? -0.125 : 0);
  const unsupported = new Float32Array(data.length + 2);
  unsupported.set(data);
  assert.throws(() => parseSnapshot(unsupported), /snapshot layout/);
  assert.throws(() => parsePianoSnapshot(unsupported), /snapshot layout/);
});
}

test("a native piano worker failure rejects requests, stops audio and permits a new client", async () => {
  const workers: WorkerBoundary[] = [];
  const nodes: NodeBoundary[] = [];
  class WorkerBoundary {
    onmessage: ((event: { data: object }) => void) | null = null;
    onerror: ((event: ErrorEvent) => void) | null = null;
    terminated = false;
    messages: { type: string; name?: string }[] = [];
    port: MessagePort | null = null;
    constructor() { workers.push(this); }
    postMessage(message: { type: string; name?: string; port?: MessagePort }) {
      this.messages.push(message);
      if (message.type === "init") {
        this.port = message.port!;
        queueMicrotask(() => this.onmessage?.({ data: {
          type: "ready", info: {}, board: {}, initMs: 1,
        } }));
      }
    }
    terminate() { this.terminated = true; this.port?.close(); }
  }
  class NodeBoundary {
    disconnected = false;
    disposed = false;
    producers: MessagePort[] = [];
    port = { onmessage: null, postMessage: (m: { type: string; port?: MessagePort }) => {
      if (m.port) this.producers.push(m.port);
      if (m.type === "dispose") { this.disposed = true; this.producers.forEach(p => p.close()); }
    } };
    constructor() { nodes.push(this); }
    connect() {}
    disconnect() { this.disconnected = true; }
  }
  Object.assign(globalThis, { Worker: WorkerBoundary, AudioWorkletNode: NodeBoundary });
  const context = {
    sampleRate: 48000, state: "running", destination: {},
    audioWorklet: { async addModule() {} },
  } as unknown as AudioContext;
  const { LivePianoClient } = await import("../src/live/pianoClient.ts");
  const client = await LivePianoClient.create(context);
  try {
    const errors: Error[] = [];
    client.onError(error => errors.push(error));
    const pending = client.partials(60);
    void pending.catch(() => {});
    const fail = workers[0]!.onerror!;
    const event = { message: "render trapped", preventDefault() {} } as ErrorEvent;
    fail(event);
    assert.equal(errors.length, 1, "a post-ready native error must reach the client listener");
    assert.match(errors[0]!.message, /render trapped/);
    await assert.rejects(pending, /render trapped/);
    assert.ok(workers[0]!.terminated, "a failed timer must not keep rendering or throwing");
    assert.ok(nodes[0]!.disposed && nodes[0]!.disconnected);
    const count = workers[0]!.messages.length;
    client.noteOn(60, 2);
    assert.equal(workers[0]!.messages.length, count);
    fail(event);
    assert.equal(errors.length, 1, "an in-flight error cannot report the same failure again");
    const retry = await LivePianoClient.create(context);
    retry.noteOn(60, 2);
    assert.equal(workers[1]!.messages.at(-1)?.name, "noteOn");
    retry.dispose();
    Object.assign(context, { state: "suspended", async resume() { throw new Error("output device unavailable"); } });
    await assert.rejects(LivePianoClient.create(context), /output device unavailable/);
    assert.ok(workers[2]!.terminated, "a failed final resume must clean up the newly built worker");
    assert.ok(nodes[2]!.disposed && nodes[2]!.disconnected);

    const resuming = deferred();
    const resumed = deferred();
    Object.assign(context, { state: "suspended", async resume() {
      resuming.resolve();
      await resumed.promise;
      Object.assign(context, { state: "running" });
    } });
    const starting = LivePianoClient.create(context);
    // Dispose a mistakenly returned client too, so this regression leaves no audio resources behind.
    const refused = assert.rejects(starting.then(created => { created.dispose(); }), /render trapped during resume/);
    await resuming.promise;
    workers[3]!.onerror!({ message: "render trapped during resume", preventDefault() {} } as ErrorEvent);
    await turn();
    const stoppedDuringResume = workers[3]!.terminated && nodes[3]!.disposed && nodes[3]!.disconnected;
    resumed.resolve();
    await refused;
    assert.ok(stoppedDuringResume, "a failed worker must be released even while the output device is still suspended");
    const resumedRetry = await LivePianoClient.create(context);
    resumedRetry.noteOn(64, 2);
    assert.equal(workers[4]!.messages.at(-1)?.name, "noteOn");
    resumedRetry.dispose();
  } finally {
    client.dispose();
    workers.forEach(worker => worker.terminate());
    nodes.forEach(node => node.producers.forEach(port => port.close()));
  }
});


import { readFileSync } from "node:fs";
import { runInNewContext } from "node:vm";
test("the real live sink preserves samples and stops processing when disposed", () => {
  function port() {
    return { closed: false, onmessage: null as ((event: { data: unknown }) => void) | null,
      postMessage() {}, close() { this.closed = true; } };
  }
  interface Sink {
    queuedFrames: number;
    receive(message: object): void;
    process(inputs: unknown[], outputs: Float32Array[][]): boolean;
  }
  let Processor!: new () => Sink;
  const control = port(), producer = port();
  // Run the actual worklet source with only the browser host boundary supplied.
  runInNewContext(readFileSync(new URL("../public/live/sink-worklet.js", import.meta.url), "utf8"), {
    AudioWorkletProcessor: class { port = control; },
    registerProcessor: (_name: string, ctor: typeof Processor) => { Processor = ctor; },
    currentFrame: 0, currentTime: 0, sampleRate: 48000,
  });
  const sink = new Processor();
  sink.receive({ type: "connect", port: producer });
  const samples = new Float32Array([0.25, -0.5, 0.125, 0]);
  sink.receive({ type: "audio", channels: [samples] });
  const out = new Float32Array(4);
  assert.equal(sink.process([], [[out]]), true);
  assert.deepEqual(out, samples);
  sink.receive({ type: "audio", channels: [samples] });
  assert.equal(sink.queuedFrames, 4);
  sink.receive({ type: "dispose" });
  assert.equal(sink.process([], [[new Float32Array(4)]]), false);
  assert.equal(sink.queuedFrames, 0);
  assert.ok(control.closed && producer.closed);
  assert.equal(control.onmessage, null);
  assert.equal(producer.onmessage, null);
  sink.receive({ type: "audio", channels: [samples] });
  assert.equal(sink.queuedFrames, 0, "in-flight messages cannot revive a disposed sink");
});
