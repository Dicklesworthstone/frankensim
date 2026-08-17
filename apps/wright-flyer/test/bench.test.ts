// Harness unit battery (bead wf-root-guzez.1.6.1, E0.6a).
import assert from "node:assert/strict";
import { test } from "node:test";
import { bench, percentileOf } from "../src/bench/harness.ts";

test("percentileOf: nearest-rank on sorted input, empty-safe", () => {
  assert.equal(percentileOf([1, 2, 3, 4], 0.5), 3);
  assert.equal(percentileOf([1, 2, 3, 4], 0.99), 4);
  assert.equal(percentileOf([7], 0.01), 7);
  assert.equal(percentileOf([], 0.99), 0);
});

test("bench: runs warmup+samples, reports per-op scaling and positive rates", () => {
  let calls = 0;
  const r = bench("probe", 100, () => {
    calls += 1;
    let acc = 0;
    for (let i = 0; i < 1000; i += 1) {
      acc += Math.sqrt(i);
    }
    if (acc < 0) {
      throw new Error("unreachable");
    }
  }, 20, 5);
  assert.equal(calls, 25, "warmup runs are executed then discarded");
  assert.equal(r.samples, 20);
  assert.equal(r.batchSize, 100);
  assert.ok(r.p50_us > 0 && r.p99_us >= r.p50_us, "ordered positive percentiles");
  assert.ok(r.opsPerSec > 0);
});
