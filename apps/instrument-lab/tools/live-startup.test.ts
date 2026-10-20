import assert from "node:assert/strict";
import { test } from "bun:test";
import { loadWorkletModule } from "../src/live/liveAudio";
import { KitLane } from "../src/live/mix";

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
