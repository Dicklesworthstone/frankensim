// E8.1 battery (bead wf-root-guzez.9.1): sweep-engine laws. CRN
// verified per member (the SAME seed at every design point — the
// seed matrix's rows are constant); RunSpecId cache hits exactly
// (second sweep calls the runner ZERO times; a model-version bump
// misses); the QoS-subordination law EXECUTED (closed gate = zero
// runner calls = sweeps cannot add deadline misses, V-14); progress
// streams; CSV deterministic; caps at cap AND cap+1.
// Repro: node --test test/sweeps.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_AXIS_POINTS,
  MAX_DESIGN_POINTS,
  MAX_ENSEMBLE,
  exportCsv,
  makeSweepEngine,
  makeSweepGrid,
  runSpecId,
  type DesignPoint,
} from "../src/sweeps.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-sweeps","case":"${kase}",${payload}}`);
}

function grid2d(): DesignPoint[] {
  const g = makeSweepGrid(
    { name: "headwind", values: [9, 10, 11] },
    { name: "rho", values: [1.25, 1.294] },
  );
  assert.ok(g.ok);
  return g.ok ? g.value : [];
}

test("CRN: every member uses the SAME seed at every design point", () => {
  const points = grid2d();
  const engine = makeSweepEngine(points, [101, 202, 303], "wf-model-v1", new Map());
  assert.ok(engine.ok);
  if (!engine.ok) return;
  const e = engine.value;
  const open = () => true;
  const runner = (p: DesignPoint, seed: number) => seed * 1000 + p.index;
  while (e.step(open, runner)) {
    /* drain */
  }
  assert.equal(e.records.length, 6 * 3);
  // The seed matrix: rows (members) constant across points.
  for (const member of [0, 1, 2]) {
    const seeds = new Set(
      e.records.filter((r) => r.member === member).map((r) => r.seed),
    );
    assert.equal(seeds.size, 1, `member ${member} must ride ONE realization`);
  }
  // And distinct members ride distinct realizations.
  const distinct = new Set(e.records.map((r) => r.seed));
  assert.equal(distinct.size, 3);
  jlog("crn", `"members":3,"points":6`);
});

test("RunSpecId cache: exact hits, model-version bump misses", () => {
  const points = grid2d();
  const cache = new Map<string, number>();
  let calls = 0;
  const runner = (p: DesignPoint, seed: number) => {
    calls += 1;
    return seed + p.index;
  };
  const open = () => true;
  const first = makeSweepEngine(points, [7], "wf-model-v1", cache);
  assert.ok(first.ok);
  if (first.ok) {
    while (first.value.step(open, runner)) {
      /* drain */
    }
  }
  assert.equal(calls, 6, "first sweep computes");
  // Second sweep over the SAME specs: zero runner calls, all cached.
  const second = makeSweepEngine(points, [7], "wf-model-v1", cache);
  assert.ok(second.ok);
  if (second.ok) {
    while (second.value.step(open, runner)) {
      /* drain */
    }
    assert.ok(second.value.records.every((r) => r.fromCache));
  }
  assert.equal(calls, 6, "cache hit by RunSpecId — zero recomputes");
  // A model-version bump is a full miss (the key is load-bearing).
  const bumped = makeSweepEngine(points, [7], "wf-model-v2", cache);
  assert.ok(bumped.ok);
  if (bumped.ok) {
    while (bumped.value.step(open, runner)) {
      /* drain */
    }
  }
  assert.equal(calls, 12, "version bump recomputes");
  // The id itself is order-stable over config keys.
  const p = points[0];
  if (p !== undefined) {
    assert.equal(
      runSpecId(p, 7, "wf-model-v1"),
      runSpecId({ index: 0, config: { rho: 1.25, headwind: 9 } }, 7, "wf-model-v1"),
    );
  }
  jlog("cache", `"calls_after_three_sweeps":${calls}`);
});

test("V-14: a closed QoS gate dispatches NOTHING (sweeps cannot add deadline misses)", () => {
  const points = grid2d();
  const engine = makeSweepEngine(points, [7], "wf-model-v1", new Map());
  assert.ok(engine.ok);
  if (!engine.ok) return;
  const e = engine.value;
  let calls = 0;
  const runner = () => {
    calls += 1;
    return 0;
  };
  // Gate closed: no unit runs, ever.
  for (let i = 0; i < 100; i += 1) {
    assert.equal(e.step(() => false, runner), false);
  }
  assert.equal(calls, 0, "closed gate = ZERO sweep work");
  assert.equal(e.progress().completed, 0);
  // Gate reopens: the sweep resumes exactly where it paused.
  let budget = 2;
  const throttled = () => budget > 0;
  while (e.step(throttled, () => {
    calls += 1;
    budget -= 1;
    return 0;
  })) {
    /* drain until throttle */
  }
  assert.equal(calls, 2, "throttled gate admits exactly its headroom");
  assert.equal(e.progress().completed, 2);
  jlog("qos", `"closed_gate_calls":0,"throttled_calls":2`);
});

test("progress streaming + deterministic CSV", () => {
  const points = grid2d();
  const engine = makeSweepEngine(points, [11, 22], "wf-model-v1", new Map());
  assert.ok(engine.ok);
  if (!engine.ok) return;
  const e = engine.value;
  assert.deepEqual(e.progress(), { completed: 0, total: 12 });
  const runner = (p: DesignPoint, seed: number) => seed + p.index * 0.5;
  const open = () => true;
  let seen = 0;
  while (e.step(open, runner)) {
    seen += 1;
    assert.equal(e.progress().completed, seen, "progress streams per unit");
  }
  const csv = exportCsv(e, points);
  const lines = csv.split("\n");
  assert.equal(lines[0], "point,headwind,rho,member,seed,value,from_cache");
  assert.equal(lines.length, 13, "header + 12 rows");
  // Deterministic: a re-run of an identical engine exports the same
  // bytes.
  const twin = makeSweepEngine(points, [11, 22], "wf-model-v1", new Map());
  assert.ok(twin.ok);
  if (twin.ok) {
    while (twin.value.step(open, runner)) {
      /* drain */
    }
    assert.equal(exportCsv(twin.value, points), csv, "bit-identical CSV");
  }
  jlog("csv", `"rows":12`);
});

test("caps at cap AND cap+1", () => {
  const axis = (n: number) => ({
    name: "a",
    values: Array.from({ length: n }, (_, i) => i),
  });
  assert.ok(makeSweepGrid(axis(MAX_AXIS_POINTS)).ok, "AT axis cap");
  const overAxis = makeSweepGrid(axis(MAX_AXIS_POINTS + 1));
  assert.ok(!overAxis.ok && overAxis.refusal.code === "sweep-axis-invalid");
  // Total cap: 32x32 = 1024 admits; 33x32 refuses.
  assert.ok(makeSweepGrid(axis(32), { name: "b", values: axis(32).values }).ok);
  const overTotal = makeSweepGrid(axis(33), { name: "b", values: axis(32).values });
  assert.ok(!overTotal.ok && overTotal.refusal.code === "sweep-grid-too-large");
  // Ensemble caps + duplicate-seed refusal + unnamed model.
  const pts = grid2d();
  const seeds = (n: number) => Array.from({ length: n }, (_, i) => i);
  assert.ok(makeSweepEngine(pts, seeds(MAX_ENSEMBLE), "m", new Map()).ok, "AT ensemble cap");
  const overE = makeSweepEngine(pts, seeds(MAX_ENSEMBLE + 1), "m", new Map());
  assert.ok(!overE.ok && overE.refusal.code === "sweep-ensemble-invalid");
  const dup = makeSweepEngine(pts, [1, 1], "m", new Map());
  assert.ok(!dup.ok && dup.refusal.code === "sweep-ensemble-invalid");
  const unnamed = makeSweepEngine(pts, [1], "  ", new Map());
  assert.ok(!unnamed.ok && unnamed.refusal.code === "sweep-model-version-missing");
  jlog("caps", `"axis":${MAX_AXIS_POINTS},"total":${MAX_DESIGN_POINTS}`);
});
