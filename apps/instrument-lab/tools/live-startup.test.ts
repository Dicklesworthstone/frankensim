import assert from "node:assert/strict";
import { test } from "bun:test";
import { loadWorkletModule } from "../src/live/liveAudio";
import { KitLane, PianoLane } from "../src/live/mix";

test("kit lane startup and replacement forward the caller's full atmosphere record", async () => {
  const previousWorker = globalThis.Worker;
  const workers: WorkerBoundary[] = [], ports: MessagePort[] = [];
  class WorkerBoundary {
    onmessage: ((event: { data: object }) => void) | null = null;
    onerror: ((event: ErrorEvent) => void) | null = null;
    messages: { type: string; gasRecord?: string; port?: MessagePort }[] = [];
    constructor() { workers.push(this); }
    postMessage(message: { type: string; gasRecord?: string; port?: MessagePort }) {
      this.messages.push(message); if (message.port) ports.push(message.port);
    }
    terminate() {}
  }
  Object.assign(globalThis, { Worker: WorkerBoundary });
  const mix = { attach(_id: string, port: MessagePort) { ports.push(port); }, dropLane() {} } as unknown as ConstructorParameters<typeof KitLane>[0];
  const lanes: KitLane[] = [];
  try {
    for (const gasRecord of ["293.15 101325 0.6", "mixture 300 80000 mass 2 N2 0.8 CO2 0.2"]) {
      const lane = new KitLane(mix, "kick", { gain: 0.2, gasRecord }); lanes.push(lane);
      const worker = workers.at(-1)!;
      assert.equal(worker.messages[0]!.gasRecord, gasRecord);
      worker.onmessage!({ data: { type: "ready", gasRecord } });
      assert.equal((await lane.ready).gasRecord, gasRecord);
      lane.dispose();
    }
  } finally {
    lanes.forEach(lane => lane.dispose()); ports.forEach(port => port.close());
    Object.assign(globalThis, { Worker: previousWorker });
  }
});

test("failed worklet URLs change on retry; concurrent and successful loads are shared", async () => {
  const urls: string[] = [];
  const context = { audioWorklet: { async addModule(url: string) {
    urls.push(url);
    // Chrome retains failed module URLs within an audio context.
    if (!url.includes("retry=")) throw new Error("temporary download failure");
  } } } as unknown as AudioContext;
  for (const url of ["/live/sink-worklet.js", "/live/mix-worklet.js?version=1"]) {
    const first = loadWorkletModule(context, url);
    assert.equal(loadWorkletModule(context, url), first);
    await assert.rejects(first, /temporary download failure/);
    const retried = loadWorkletModule(context, url);
    assert.equal(loadWorkletModule(context, url), retried);
    await retried;
    assert.equal(loadWorkletModule(context, url), retried);
  }
  assert.deepEqual(urls, ["/live/sink-worklet.js", "/live/sink-worklet.js?retry=1",
    "/live/mix-worklet.js?version=1", "/live/mix-worklet.js?version=1&retry=1"]);
});

for (const failure of ["startup-message", "startup-native", "attach", "running-message", "running-native"]) {
  test(`failed lane releases its worker, settles readiness and permits replacement: ${failure}`, async () => {
    const previousWorker = globalThis.Worker;
    const workers: WorkerBoundary[] = [];
    const attached: MessagePort[] = [];
    const dropped: string[] = [];
    class WorkerBoundary {
      onmessage: ((event: { data: object }) => void) | null = null;
      onerror: ((event: ErrorEvent) => void) | null = null;
      terminated = false;
      constructor() { workers.push(this); }
      postMessage() {}
      terminate() { this.terminated = true; }
    }
    Object.assign(globalThis, { Worker: WorkerBoundary });
    let refuseAttach = failure === "attach";
    const mix = {
      attach(_id: string, port: MessagePort) {
        if (refuseAttach) throw new Error("mixer port unavailable");
        attached.push(port);
      },
      dropLane(id: string) { dropped.push(id); attached.forEach(p => p.close()); },
    } as unknown as ConstructorParameters<typeof KitLane>[0];
    const lane = new KitLane(mix, "kick", { gain: 0.2 });
    try {
      const reports: Error[] = [];
      lane.onFailure = e => reports.push(e);
      const worker = workers[0]!;
      const message = worker.onmessage!;
      const native = worker.onerror!;
      const running = failure.startsWith("running-");
      if (running) { message({ data: { type: "ready" } }); await lane.ready; }
      if (failure === "attach") {
        assert.doesNotThrow(() => message({ data: { type: "ready" } }));
      } else if (failure.endsWith("native")) {
        let prevented = false;
        native({ message: "render trapped", preventDefault() { prevented = true; } } as ErrorEvent);
        assert.ok(prevented);
      } else message({ data: { type: "error", message: "render trapped" } });
      if (!running) await assert.rejects(lane.ready, failure === "attach" ? /mixer port unavailable/ : /render trapped/);
      assert.ok(worker.terminated);
      assert.equal(worker.onmessage, null);
      assert.equal(worker.onerror, null);
      assert.equal(reports.length, running ? 1 : 0);
      assert.deepEqual(dropped, running ? ["kick"] : [], "unattached failures must not drop another lane with the same ID");
      message({ data: { type: "ready" } });
      assert.equal(attached.length, running ? 1 : 0, "late readiness cannot revive a failed lane");
      refuseAttach = false;
      const replacement = new KitLane(mix, "kick", { gain: 0.2 });
      workers[1]!.onmessage!({ data: { type: "ready" } });
      await replacement.ready;
      replacement.dispose();
      assert.ok(workers[1]!.terminated);
    } finally {
      lane.dispose();
      workers.forEach(w => w.terminate());
      attached.forEach(p => p.close());
      Object.assign(globalThis, { Worker: previousWorker });
    }
  });
}

for (const kind of ["kit", "piano"] as const) for (const failure of ["construct", "init"] as const) {
  test(`${kind} lane releases resources when worker ${failure} throws synchronously`, async () => {
    const previousWorker = globalThis.Worker, previousChannel = globalThis.MessageChannel;
    const ports: { closeCount: number; close(): void }[] = [];
    const workers: WorkerBoundary[] = [];
    const refusal = new Error(`worker ${failure} refused`);
    let fail = true;
    class WorkerBoundary {
      onmessage: ((event: { data: object }) => void) | null = null;
      onerror: ((event: ErrorEvent) => void) | null = null;
      terminated = false;
      constructor() {
        if (fail && failure === "construct") throw refusal;
        workers.push(this);
      }
      postMessage() { if (fail && failure === "init") throw refusal; }
      terminate() { this.terminated = true; }
    }
    const makePort = () => {
      const port = { closeCount: 0, close() { this.closeCount++; } };
      ports.push(port);
      return port;
    };
    Object.assign(globalThis, { Worker: WorkerBoundary, MessageChannel: class {
      port1 = makePort();
      port2 = makePort();
    } });
    const mix = { attach() {}, dropLane() {} } as unknown as ConstructorParameters<typeof KitLane>[0];
    const create = () => kind === "kit" ? new KitLane(mix, "kick", { gain: 0.2 }) : new PianoLane(mix);
    let replacement: KitLane | PianoLane | undefined;
    try {
      assert.throws(create, error => error === refusal, "preserve the original startup error");
      assert.equal(ports.length, 2);
      assert.ok(ports.every(p => p.closeCount === 1), "both untransferred ports must close");
      assert.ok(workers.every(w => w.terminated), "failed init must terminate the worker");
      assert.ok(workers.every(w => w.onmessage === null && w.onerror === null));
      fail = false;
      replacement = create();
      workers.at(-1)!.onmessage!({ data: { type: "ready" } });
      await replacement.ready;
    } finally {
      replacement?.dispose();
      Object.assign(globalThis, { Worker: previousWorker, MessageChannel: previousChannel });
    }
  });
}
