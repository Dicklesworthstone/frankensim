// E7.3 battery (bead wf-root-guzez.8.4): the passthrough LAW.
// Overlay values equal sim-plane state bit-for-bit (per-strip
// Object.is oracle, NaN-safe, negative-zero-strict); the
// RECOMPUTATION FALSIFIER feeds a state whose components do NOT sum
// to its net and requires the gnomon to show the state's net
// verbatim (a recomputing renderer would show the sum); doctored
// overlays are caught by the exported audit; probe charts match the
// logged series bit-for-bit with caps at cap AND cap+1.
// Repro: node --test test/forceOverlay.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_CHART_SAMPLES,
  MAX_STRIPS,
  buildForceOverlay,
  buildProbeChart,
  firstChartDivergence,
  firstOverlayDivergence,
  type StripLoadsState,
} from "../src/forceOverlay.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-forceoverlay","case":"${kase}",${payload}}`);
}

function makeState(nStrips: number): StripLoadsState {
  const positions = new Float64Array(3 * nStrips);
  const forces = new Float64Array(3 * nStrips);
  for (let i = 0; i < nStrips; i += 1) {
    positions[3 * i] = i * 0.5;
    positions[3 * i + 2] = 2.5;
    forces[3 * i] = -12.5 - i;
    forces[3 * i + 1] = i % 2 === 0 ? 0 : -0; // negative zero matters
    forces[3 * i + 2] = 160.25 + i * 0.125;
  }
  return {
    nStrips,
    positions,
    forces,
    thrustN: [180.5, 0, 3.25],
    thrustAt: [-1.2, 0, 2.4],
    weightN: [0, 0, -3315.75],
    cgAt: [0.1, 0, 1.9],
    // DELIBERATELY not the sum of the components: the sim plane owns
    // this number; a recomputing overlay would disagree with it.
    netN: [17.125, -0.5, 42.0],
  };
}

test("overlay is bit-for-bit the state (incl. negative zero), audit clean", () => {
  const state = makeState(8);
  const r = buildForceOverlay(state);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.equal(firstOverlayDivergence(r.value, state), -1);
  // Per-strip verbatim, negative zero preserved (Object.is strict).
  assert.ok(Object.is(r.value.stripVectors[4], -0));
  // A copy, not a view: later state mutation cannot retint it.
  state.forces[0] = 999;
  assert.equal(r.value.stripVectors[0], -12.5);
  jlog("verbatim", `"divergence":-1`);
});

test("RECOMPUTATION FALSIFIER: gnomon shows the state's net, never a renderer sum", () => {
  const state = makeState(4);
  const r = buildForceOverlay(state);
  assert.ok(r.ok);
  if (!r.ok) return;
  // The components deliberately do NOT sum to netN…
  let sum = 0;
  for (let i = 0; i < state.nStrips; i += 1) sum += state.forces[3 * i] ?? 0;
  sum += state.thrustN[0] + state.weightN[0];
  assert.notEqual(sum, state.netN[0], "fixture must be inconsistent");
  // …and the overlay must show netN VERBATIM anyway.
  assert.ok(Object.is(r.value.net[0], 17.125));
  // A doctored overlay (the recomputation) is CAUGHT by the audit.
  const doctored = { ...r.value, net: [sum, r.value.net[1], r.value.net[2]] as const };
  assert.notEqual(firstOverlayDivergence(doctored, state), -1);
  jlog("recompute-falsifier", `"renderer_sum":${sum},"state_net":17.125`);
});

test("doctored strip vector is caught at its exact index", () => {
  const state = makeState(6);
  const r = buildForceOverlay(state);
  assert.ok(r.ok);
  if (!r.ok) return;
  const vecs = Float64Array.from(r.value.stripVectors);
  vecs[7] = vecs[7]! * (1 + 1e-16) + 1e-300; // one ulp-ish nudge
  const doctored = { ...r.value, stripVectors: vecs };
  assert.equal(firstOverlayDivergence(doctored, state), 7, "first divergence localized");
  jlog("doctored-strip", `"index":7`);
});

test("strip caps at cap AND cap+1, mismatched arrays refuse", () => {
  assert.ok(buildForceOverlay(makeState(MAX_STRIPS)).ok, "AT cap");
  const over = buildForceOverlay(makeState(MAX_STRIPS + 1));
  assert.ok(!over.ok && over.refusal.code === "strip-count-invalid");
  const bad = { ...makeState(4), forces: new Float64Array(11) };
  const r = buildForceOverlay(bad);
  assert.ok(!r.ok && r.refusal.code === "strip-arrays-mismatched");
  jlog("caps", `"max_strips":${MAX_STRIPS}`);
});

test("probe charts match the logged series bit-for-bit", () => {
  const ticks = Float64Array.from([10, 11, 12, 13]);
  const values = Float64Array.from([0.5, -0, Number.MIN_VALUE, 42]);
  const r = buildProbeChart("hinge moment", 1.5, ticks, values);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.equal(firstChartDivergence(r.value, ticks, values), -1);
  // A resampled/smoothed chart is CAUGHT.
  const smoothed = {
    ...r.value,
    values: Float64Array.from([0.5, 0.25, Number.MIN_VALUE, 42]),
  };
  assert.equal(firstChartDivergence(smoothed, ticks, values), 1);
  assert.equal(r.value.referenceHeightM, 1.5, "reference height declared");
  jlog("chart", `"divergence":-1`);
});

test("probe chart refusals: caps, order, reference height", () => {
  const mono = (n: number) => Float64Array.from({ length: n }, (_, i) => i);
  assert.ok(buildProbeChart("x", 0, mono(MAX_CHART_SAMPLES), mono(MAX_CHART_SAMPLES)).ok);
  const over = buildProbeChart("x", 0, mono(MAX_CHART_SAMPLES + 1), mono(MAX_CHART_SAMPLES + 1));
  assert.ok(!over.ok && over.refusal.code === "probe-series-length-invalid");
  const unordered = buildProbeChart(
    "x",
    0,
    Float64Array.from([1, 3, 2]),
    Float64Array.from([0, 0, 0]),
  );
  assert.ok(!unordered.ok && unordered.refusal.code === "probe-series-unordered");
  const mismatched = buildProbeChart("x", 0, mono(3), mono(4));
  assert.ok(!mismatched.ok && mismatched.refusal.code === "probe-series-mismatched");
  const badRef = buildProbeChart("x", Number.NaN, mono(2), mono(2));
  assert.ok(!badRef.ok && badRef.refusal.code === "probe-reference-invalid");
  jlog("chart-refusals", `"max_samples":${MAX_CHART_SAMPLES}`);
});
