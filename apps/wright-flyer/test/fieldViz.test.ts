// E7.2 battery (bead wf-root-guzez.8.3): field-viz logic laws.
// Glyph budget at cap AND cap+1 (typed refusal, never silent
// downsampling); per-item glyph filtering (invalid + core excluded);
// divergence overlay teaching toggle matches the API duals with an
// identical mask law; RK4 streamline determinism + the solid-rotation
// radius oracle + domain exit; wake age-fade monotone
// presentation-only; legend HONESTY falsifier ("total flow" never
// with a supported force-coupled omission); probe binding caps.
// Repro: node --test test/fieldViz.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  GRAD_NORM_FLOOR,
  MAX_GLYPHS,
  MAX_PROBES,
  MAX_SEEDS,
  MAX_STEPS,
  WAKE_FADE_PRESENTATION_ONLY,
  bindProbes,
  buildGlyphInstances,
  divergenceOverlay,
  integrateStreamlines,
  legendConfig,
  wakeAgeFade,
  type FieldArrays,
} from "../src/fieldViz.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-fieldviz","case":"${kase}",${payload}}`);
}

function makeField(n: number, overrides?: Partial<Record<string, unknown>>): FieldArrays {
  const points = new Float64Array(3 * n);
  const u = new Float64Array(3 * n);
  const divAnalytic = new Float64Array(n);
  const divFd = new Float64Array(n);
  const gradNorm = new Float64Array(n);
  const validity = new Uint8Array(n).fill(1);
  const singularityCore = new Uint8Array(n);
  for (let i = 0; i < n; i += 1) {
    points[3 * i] = i;
    u[3 * i] = 8;
    u[3 * i + 2] = 0.5;
    divAnalytic[i] = 1e-12;
    divFd[i] = 2e-9;
    gradNorm[i] = 0.3;
  }
  return {
    n,
    points,
    u,
    divAnalytic,
    divFd,
    gradNorm,
    validity,
    singularityCore,
    omittedComponents: ["physical-wake"],
    forceCoupledSupported: ["bound-circulation"],
    ...(overrides as object),
  } as FieldArrays;
}

test("glyph budget: AT the cap admits, one more refuses typed", () => {
  const at = buildGlyphInstances(makeField(MAX_GLYPHS));
  assert.ok(at.ok);
  assert.equal(at.ok && at.value.count, MAX_GLYPHS);
  const over = buildGlyphInstances(makeField(MAX_GLYPHS + 1));
  assert.ok(!over.ok);
  assert.equal(!over.ok && over.refusal.code, "glyph-count-exceeded");
  jlog("glyph-cap", `"cap":${MAX_GLYPHS}`);
});

test("glyph filtering: invalid and core points excluded per item", () => {
  const f = makeField(6);
  f.validity[1] = 0;
  f.singularityCore[3] = 1;
  const r = buildGlyphInstances(f);
  assert.ok(r.ok);
  if (r.ok) {
    assert.equal(r.value.count, 4);
    // Positions of kept glyphs are exactly points 0,2,4,5 (x = index).
    const xs = [r.value.positions[0], r.value.positions[3], r.value.positions[6], r.value.positions[9]];
    assert.deepEqual(xs, [0, 2, 4, 5]);
    // Directions are unit vectors of u.
    const mag = Math.hypot(8, 0, 0.5);
    assert.ok(Math.abs(r.value.directions[0] - 8 / mag) < 1e-15);
    assert.ok(Math.abs(r.value.magnitudes[0] - mag) < 1e-15);
  }
  jlog("glyph-filter", `"kept":4`);
});

test("divergence overlay: the teaching toggle matches the duals, mask law identical", () => {
  const f = makeField(5);
  f.gradNorm[1] = GRAD_NORM_FLOOR / 2; // under floor -> masked
  f.singularityCore[2] = 1; // core -> masked
  f.validity[3] = 0; // invalid -> masked
  const an = divergenceOverlay(f, "analytic");
  const fd = divergenceOverlay(f, "finite-difference");
  for (let i = 0; i < 5; i += 1) {
    assert.equal(an.absolute[i], Math.abs(f.divAnalytic[i]), `analytic abs ${i}`);
    assert.equal(fd.absolute[i], Math.abs(f.divFd[i]), `fd abs ${i}`);
    assert.equal(an.masked[i], fd.masked[i], `mask parity ${i}`);
  }
  assert.deepEqual(Array.from(an.masked), [0, 1, 1, 1, 0]);
  assert.ok(Number.isNaN(an.normalized[1]) && Number.isNaN(fd.normalized[1]));
  assert.ok(Math.abs(an.normalized[0] - Math.abs(f.divAnalytic[0]) / 0.3) < 1e-18);
  // The duals genuinely differ (the toggle is not vacuous).
  assert.notEqual(an.absolute[0], fd.absolute[0]);
  jlog("divergence-toggle", `"masked":[0,1,1,1,0]`);
});

test("streamlines: uniform field is exact, solid rotation preserves radius, deterministic", () => {
  // Uniform field: after k steps the line is exactly seed + k*step*u.
  const uniform = integrateStreamlines(() => [2, 0, 0] as const, [[0, 0, 0]], 0.1, 10);
  assert.ok(uniform.ok);
  if (uniform.ok) {
    const pts = uniform.value[0].points;
    assert.equal(pts.length, 3 * 11);
    assert.ok(Math.abs(pts[30] - 2.0) < 1e-12, `x after 10 steps: ${pts[30]}`);
  }
  // Solid rotation about z: RK4 must hold the radius to ~1e-8 over a
  // quarter turn (the classical conservation oracle).
  const rot = (p: readonly [number, number, number]) => [-p[1], p[0], 0] as const;
  const quarter = integrateStreamlines(rot, [[1, 0, 0]], 0.001, 1571);
  assert.ok(quarter.ok);
  if (quarter.ok) {
    const pts = quarter.value[0].points;
    const last = pts.length - 3;
    const r = Math.hypot(pts[last], pts[last + 1]);
    assert.ok(Math.abs(r - 1) < 1e-8, `radius drift ${r - 1}`);
  }
  // Determinism: bitwise identical polylines.
  const a = integrateStreamlines(rot, [[1, 0, 0]], 0.01, 200);
  const b = integrateStreamlines(rot, [[1, 0, 0]], 0.01, 200);
  assert.ok(a.ok && b.ok);
  if (a.ok && b.ok) assert.deepEqual(a.value[0].points, b.value[0].points);
  // Domain exit: sampler null past x=1 -> ended left-domain, partial line.
  const bounded = integrateStreamlines(
    (p) => (p[0] > 1 ? null : ([2, 0, 0] as const)),
    [[0, 0, 0]],
    0.1,
    100,
  );
  assert.ok(bounded.ok);
  if (bounded.ok) {
    assert.equal(bounded.value[0].ended, "left-domain");
    assert.ok(bounded.value[0].points.length < 3 * 101);
  }
  jlog("streamlines", `"radius_oracle":"1e-8"`);
});

test("streamline caps at cap AND cap+1", () => {
  const s = () => [1, 0, 0] as const;
  const seeds: Array<readonly [number, number, number]> = Array.from(
    { length: MAX_SEEDS },
    (_, i) => [i, 0, 0] as const,
  );
  assert.ok(integrateStreamlines(s, seeds, 0.1, 1).ok, "AT seed cap");
  const overSeeds = [...seeds, [0, 0, 0] as const];
  const r1 = integrateStreamlines(s, overSeeds, 0.1, 1);
  assert.ok(!r1.ok && r1.refusal.code === "streamline-seeds-invalid");
  assert.ok(integrateStreamlines(s, [[0, 0, 0]], 0.1, MAX_STEPS).ok, "AT step cap");
  const r2 = integrateStreamlines(s, [[0, 0, 0]], 0.1, MAX_STEPS + 1);
  assert.ok(!r2.ok && r2.refusal.code === "streamline-steps-invalid");
  const r3 = integrateStreamlines(s, [[0, 0, 0]], 0, 10);
  assert.ok(!r3.ok && r3.refusal.code === "streamline-step-invalid");
  jlog("streamline-caps", `"seeds":${MAX_SEEDS},"steps":${MAX_STEPS}`);
});

test("wake age-fade: monotone, clamped, presentation-only", () => {
  assert.equal(wakeAgeFade(0, 100), 1);
  assert.equal(wakeAgeFade(100, 100), 0);
  assert.equal(wakeAgeFade(150, 100), 0);
  assert.equal(wakeAgeFade(-5, 100), 1);
  let prev = 1;
  for (let a = 0; a <= 100; a += 5) {
    const v = wakeAgeFade(a, 100);
    assert.ok(v <= prev + 1e-15, `monotone at ${a}`);
    prev = v;
  }
  assert.equal(wakeAgeFade(10, 0), 0, "degenerate maxAge is opaque-none");
  assert.match(WAKE_FADE_PRESENTATION_ONLY, /presentation-only/);
  jlog("age-fade", `"monotone":true`);
});

test("legend honesty: 'total flow' never with a supported force-coupled omission", () => {
  // Omitting a SUPPORTED force-coupled component: label demoted and
  // the omissions are named.
  const dishonestCandidate = makeField(2, {
    omittedComponents: ["bound-circulation", "visualization-only"],
    forceCoupledSupported: ["bound-circulation"],
  });
  const demoted = legendConfig("velocity", dishonestCandidate);
  assert.equal(demoted.label, "selected components");
  assert.ok(demoted.omittedComponents.includes("bound-circulation"));
  assert.equal(demoted.alwaysVisible, true);
  // Only unsupported/vis-only omissions: the claim is honest.
  const honest = makeField(2, {
    omittedComponents: ["physical-wake", "visualization-only"],
    forceCoupledSupported: ["bound-circulation"],
  });
  assert.equal(legendConfig("velocity", honest).label, "total flow");
  // Non-velocity legends never claim totality.
  assert.equal(legendConfig("divergence-normalized", honest).label, "eps_div (normalized, masked)");
  jlog("legend-honesty", `"falsifier":"executed"`);
});

test("probe binding: nearest point, deterministic tie-break, cap", () => {
  const f = makeField(10);
  const r = bindProbes(f, [[4.4, 0, 0], [0, 0, 0]]);
  assert.ok(r.ok);
  if (r.ok) {
    assert.equal(r.value[0].pointIndex, 4);
    assert.equal(r.value[1].pointIndex, 0);
    assert.equal(r.value[0].valid, true);
  }
  const many = Array.from({ length: MAX_PROBES + 1 }, () => [0, 0, 0] as const);
  const over = bindProbes(f, many);
  assert.ok(!over.ok && over.refusal.code === "probe-count-exceeded");
  assert.ok(bindProbes(f, many.slice(0, MAX_PROBES)).ok, "AT the cap admits");
  jlog("probes", `"cap":${MAX_PROBES}`);
});
