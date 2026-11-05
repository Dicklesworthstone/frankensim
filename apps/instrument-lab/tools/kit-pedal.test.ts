import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { runInNewContext } from "node:vm";
import { test } from "node:test";

// Run the production worker and audio-frame scheduler. Only WASM physics and
// browser ports are substituted; native rod/contact trajectories have Rust tests.
async function kit(options: Record<string, unknown> = {}, prepare = true) {
  let wall = 0, benchmarking = true, corrections = 0;
  const tasks: (() => void)[] = [], messages: any[] = [], audio: any[] = [];
  const instances: Piece[] = [];
  class Piece {
    pedal = 0;
    hits = 0;
    active = false;
    freed = false;
    inputs: unknown[] = [];
    changes: number[][] = [];
    name: string;
    scale: number;
    constructor(name: string, scale: number) { this.name = name; this.scale = scale; instances.push(this); }
    static withGas(...args: [string, number, ...unknown[]]) {
      const p = new Piece(args[0], args[1]); p.inputs = args.slice(2); return p;
    }
    static withMaterials = Piece.withGas;
    static withGasSources = Piece.withGas;
    prepareHiHatPedal(force: number) { this.pedal = force; }
    copy() {
      const p = new Piece(this.name, this.scale);
      p.pedal = this.pedal; p.inputs = this.inputs; return p;
    }
    setPedalForce(force: number, seconds: number) {
      this.changes.push([force, seconds]); this.pedal = force; this.active = true;
    }
    strike() { this.hits++; this.active = true; }
    render(n: number, gain: number) {
      if (benchmarking) wall++;
      return new Float32Array(n).fill((this.pedal + this.hits * 100) * gain);
    }
    asleep() { return !this.active; }
    free() { this.freed = true; }
    radius() { return 0.2; }
    surfaceXY() { return new Float32Array([0, 0]); }
    surfaceH() { return 0.01; }
    info() { return JSON.stringify({ piece: this.name }); }
    force() { return 0; }
    kinetic() { return 0; }
    velocity() { return new Float32Array(1); }
  }
  if (!prepare) Object.defineProperty(Piece.prototype, "prepareHiHatPedal", { value: undefined });
  const port = { onmessage: null as any, postMessage(m: any) { audio.push(m); } };
  const self = { onmessage: null as any, postMessage(m: any) { messages.push(m); } };
  class Channel {
    port1 = { onmessage: null as any };
    port2 = { postMessage: () => tasks.push(() => this.port1.onmessage({ data: 0 })) };
  }
  const lane = readFileSync(new URL("../public/live/lane.js", import.meta.url), "utf8")
    .replace("export function createLane", "function createLane");
  const worker = readFileSync(new URL("../public/live/kit/worker.js", import.meta.url), "utf8")
    .replace(/^import .*;\n/gm, "")
    .replaceAll("import.meta.url", '"https://example.invalid/live/kit/worker.js"');
  runInNewContext(`${lane}\n${worker}`, {
    self, KitPiece: Piece, init: async () => {}, MessageChannel: Channel, setInterval: () => 0,
    performance: { now: () => wall }, URL, Float32Array, Uint8Array,
    fetch: async () => ({ ok: true, clone() { return this; },
      arrayBuffer: async () => new ArrayBuffer(0), json: async () => ({ measured: true }) }),
    gasSourcePacket: (_record: unknown, sources: unknown) => sources,
    levelFix: () => { corrections++; return 7; },
  });
  const send = (m: Record<string, unknown>) => self.onmessage({ data: m });
  await send({ type: "init", piece: "hihat", port, gain: 1, scales: [7],
    minRealtime: 1, benchMs: 40, bursts: false, fieldHz: 0,
    hihatPedal: { initialForceN: 80 }, ...options });
  if (prepare) assert.ok(messages.some(m => m.type === "ready"), JSON.stringify(messages));
  benchmarking = false;
  const live = instances.find(p => !p.freed);
  return { send, messages, instances, live, corrections,
    clock() {
      port.onmessage({ data: { type: "clock", frame: 0, sampleRate: 48000 } });
      let turns = 0;
      while (tasks.length) { assert.ok(turns++ < 30, "scheduler must yield at its horizon"); tasks.shift()!(); }
    },
    sample(frame: number) {
      const block = audio.find(m => m.frame <= frame && m.frame + m.channels[0].length > frame);
      assert.ok(block, `no audio at frame ${frame}`);
      return block.channels[0][frame - block.frame];
    },
  };
}

for (const mode of ["direct", "score", "sequence"]) {
  test(`${mode}: pedal moves the ringing instrument at its audio frame without a stick strike`, async () => {
    const k = await kit();
    const events = [
      { frame: 512, x: 0.05, y: 0, speed: 2 },
      { frame: 641, pedalForceN: 0, transitionSeconds: 0.12 },
      { frame: 700, pedalForceN: 80, transitionSeconds: 0.12 },
      { frame: 700, x: 0.05, y: 0, speed: 2 },
    ];
    if (mode === "score") await k.send({ type: "score", events });
    else if (mode === "sequence") await k.send({ type: "sequence", startFrame: 512,
      loopFrames: 2048, events: events.map(({ frame, ...e }) => ({ offset: frame - 512, ...e })) });
    else for (const e of events) await k.send("pedalForceN" in e
      ? { type: "pedal", frame: e.frame, forceN: e.pedalForceN, transitionSeconds: e.transitionSeconds }
      : { type: "strike", ...e });
    assert.equal(k.live!.hits, 0, "events wait for the audio clock");
    k.clock();
    assert.equal(k.live!.hits, 2, "foot movement cannot manufacture a stick impact");
    assert.deepEqual(k.live!.changes, [[0, 0.12], [80, 0.12]]);
    assert.equal(k.sample(512), 180);
    assert.equal(k.sample(640), 180, "already rendered ringing survives until the pedal frame");
    assert.equal(k.sample(641), 100);
    assert.equal(k.sample(699), 100);
    assert.equal(k.sample(700), 280, "same-frame pedal target precedes the next stroke");
  });
}

test("initial pedal settlement preserves custom sources and precedes the benchmark clone", async () => {
  const gas = "mixture 300 80000 mole 2 N2 0.8 CO2 0.2";
  const k = await kit({ gasRecord: gas, gasSources: ["species", "thermo", "transport"] });
  assert.equal(k.live!.pedal, 80);
  assert.equal(k.live!.hits, 0, "the trial stroke belongs only to the discarded benchmark piece");
  assert.deepEqual(k.live!.inputs, [gas, "species", "thermo", "transport"]);
  const ready = k.messages.find(m => m.type === "ready");
  assert.equal(ready.hihatPedal.initialForceN, 80);
  assert.equal(ready.levelFix, 1);
  assert.equal(ready.levels, "none", "fixed open/closed calibration cannot normalize a moving pedal");
  assert.equal(k.corrections, 0);
});

test("old WASM fails explicitly and releases trials instead of dropping initial pedal force", async () => {
  const k = await kit({}, false);
  assert.ok(k.messages.some(m => m.type === "error" && /initial hi-hat pedal preparation/.test(m.message)));
  assert.ok(!k.messages.some(m => m.type === "ready"));
  assert.ok(k.instances.length > 0 && k.instances.every(p => p.freed));
});
