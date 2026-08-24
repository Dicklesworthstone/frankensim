// E0.6b unit battery for the shared bench kernels (src/bench/kernels.ts).
// Run: node --test test/benchKernels.test.ts
// Pins the cross-host contract: fixed kernel order, sane percentile
// fields, typed NO-DATA coverage, and node/SAB availability semantics.

import test from "node:test";
import assert from "node:assert/strict";
import { runBenchSuite, standingNoData, KERNEL_NAMES } from "../src/bench/kernels.ts";

const SAB_AVAILABLE = typeof SharedArrayBuffer !== "undefined";

function assertSaneRow(row) {
  assert.equal(typeof row.name, "string");
  assert.ok(row.batchSize > 0, `${row.name} batchSize`);
  assert.ok(row.samples > 0, `${row.name} samples`);
  assert.ok(row.p50_us > 0, `${row.name} p50 positive`);
  assert.ok(row.p50_us <= row.p95_us, `${row.name} p50 <= p95`);
  assert.ok(row.p95_us <= row.p99_us, `${row.name} p95 <= p99`);
  assert.ok(row.opsPerSec > 0, `${row.name} opsPerSec positive`);
}

test("suite returns the fixed kernel order for this context", () => {
  const { rows } = runBenchSuite();
  const expected = SAB_AVAILABLE ? [...KERNEL_NAMES] : KERNEL_NAMES.filter((n) => n !== "seqlock-publish-256f64");
  assert.deepEqual(
    rows.map((r) => r.name),
    expected,
    "kernel order is part of the cross-host artifact contract",
  );
});

test("every row carries sane monotone percentiles", () => {
  const { rows } = runBenchSuite();
  for (const row of rows) {
    assertSaneRow(row);
  }
});

test("standing NO-DATA rows are always present; seqlock/GPU NO-DATA typed in headless context", () => {
  const { noData } = runBenchSuite();
  const names = new Set(noData.map((n) => n.name));
  for (const standing of standingNoData()) {
    assert.ok(names.has(standing.name), `missing standing NO-DATA ${standing.name}`);
    assert.equal(typeof standing.reason, "string");
  }
  if (!SAB_AVAILABLE) {
    assert.ok(names.has("seqlock-publish-256f64"), "fallback origin must type its seqlock gap");
  } else {
    assert.ok(!names.has("seqlock-publish-256f64"), "measured kernel must not double-report as NO-DATA");
  }
  if (typeof WebGL2RenderingContext === "undefined") {
    assert.ok(names.has("float32-gpu-upload"), "headless context must type its GPU upload gap");
  }
});

test("two consecutive runs keep identical ordering (artifact alignment)", () => {
  const a = runBenchSuite();
  const b = runBenchSuite();
  assert.deepEqual(
    a.rows.map((r) => r.name),
    b.rows.map((r) => r.name),
  );
});
