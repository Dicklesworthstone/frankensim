import { test } from "bun:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";

interface Sink {
  receive(message: object): void;
  process(inputs: Float32Array[][], outputs: Float32Array[][]): boolean;
  queuedFrames: number;
  playedFrames: number;
  underruns: number;
  underrunFrames: number;
}

function makeSink(): Sink {
  let SinkClass!: new () => Sink;
  const source = readFileSync(new URL("../public/live/sink-worklet.js", import.meta.url), "utf8");
  // Exercise the real processor with only its browser host boundary supplied.
  new Function("AudioWorkletProcessor", "registerProcessor", "currentFrame", "currentTime", "sampleRate", source)(
    class { port = { onmessage: null, postMessage() {}, close() {} }; },
    (_name: string, processor: new () => Sink) => { SinkClass = processor; }, 0, 0, 48000,
  );
  return new SinkClass();
}

function quantum(sink: Sink, channels = 2): Float32Array[] {
  const output = Array.from({ length: channels }, () => new Float32Array(128));
  assert.equal(sink.process([], [output]), true);
  return output;
}

test("whole mono chunks retain every float bit on all output channels without temporary views", () => {
  const sink = makeSink();
  const samples = Float32Array.from({ length: 128 }, (_, i) => (i - 64) / 128);
  samples[0] = -0;
  samples[1] = 2 ** -149;
  samples.subarray = () => { throw new Error("whole chunk allocated a view"); };
  sink.receive({ type: "audio", channels: [samples] });
  for (const channel of quantum(sink, 3)) assert.deepEqual(new Uint32Array(channel.buffer), new Uint32Array(samples.buffer));
  assert.equal(sink.queuedFrames, 0);
  assert.equal(sink.playedFrames, 128);
  assert.equal(sink.underruns, 0);
});

test("stereo chunks preserve channel mapping including extra output channels", () => {
  const sink = makeSink();
  const left = new Float32Array(128).fill(0.25), right = new Float32Array(128).fill(-0.5);
  sink.receive({ type: "audio", channels: [left, right] });
  const output = quantum(sink, 3);
  assert.deepEqual(output, [left, right, right]);
});

test("partial chunks retain offsets, ordering and honest missing-frame counts", () => {
  const sink = makeSink();
  const a = Float32Array.from({ length: 200 }, (_, i) => i / 256);
  const b = Float32Array.from({ length: 33 }, (_, i) => -i / 256);
  sink.receive({ type: "audio", channels: [a] });
  sink.receive({ type: "audio", channels: [b] });
  const first = quantum(sink, 1)[0]!;
  const second = quantum(sink, 1)[0]!;
  assert.deepEqual(first, a.slice(0, 128));
  const expected = new Float32Array(128);
  expected.set(a.slice(128)); expected.set(b, 72);
  assert.deepEqual(second, expected);
  assert.equal(sink.playedFrames, 233);
  assert.equal(sink.queuedFrames, 0);
  assert.equal(sink.underruns, 1);
  assert.equal(sink.underrunFrames, 23);
});

test("whole source after a partial source copies into the correct output offset", () => {
  const sink = makeSink();
  const a = new Float32Array(17).fill(-0.25), b = new Float32Array(111).fill(0.75);
  sink.receive({ type: "audio", channels: [a] });
  sink.receive({ type: "audio", channels: [b] });
  const expected = new Float32Array(128); expected.set(a); expected.set(b, 17);
  assert.deepEqual(quantum(sink, 1)[0], expected);
  assert.equal(sink.underruns, 0);
});
