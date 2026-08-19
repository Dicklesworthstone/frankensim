// E2.2 pose battery (bead wf-root-guzez.3.2): the rig's pure core —
// slaving law, clamp-and-report at the stops (cap AND beyond), counter-
// rotation, the ±25% schematic-preview flag at the boundary and one ulp
// past, and finite-input refusals.
// Repro: node --test test/pose.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";
import {
  CANARD_TRAVEL_DEG,
  REFERENCE_DIMS,
  RUDDER_SLAVING,
  SCHEMATIC_PREVIEW_FRACTION,
  WARP_LIMIT_DEG,
  computePose,
} from "../src/airframe/pose.ts";

const base = {
  canardDeg: 0,
  warpDeg: 0,
  rudderDeg: 0,
  coupled: true,
  propAngleRad: 0,
};

test("slaved rudder follows the flown 2.5 ratio; uncoupled uses the command", () => {
  const p = computePose({ ...base, warpDeg: 4 });
  assert.ok(Math.abs(p.rudderRad - ((RUDDER_SLAVING * 4 * Math.PI) / 180)) < 1e-12);
  const free = computePose({ ...base, coupled: false, warpDeg: 4, rudderDeg: -3 });
  assert.ok(Math.abs(free.rudderRad - ((-3 * Math.PI) / 180)) < 1e-12);
  assert.equal(p.clamped, false);
  console.log(JSON.stringify({ suite: "wf-pose", case: "slaving", rudderRad: p.rudderRad }));
});

test("stops clamp AND report at the cap and beyond (never silent)", () => {
  const at = computePose({ ...base, canardDeg: CANARD_TRAVEL_DEG });
  assert.equal(at.clamped, false, "at the stop is admitted");
  const over = computePose({ ...base, canardDeg: CANARD_TRAVEL_DEG + 0.0001 });
  assert.equal(over.clamped, true, "beyond the stop must REPORT");
  assert.ok(Math.abs(over.canardRad - (CANARD_TRAVEL_DEG * Math.PI) / 180) < 1e-12);
  // Warp stop drives the slaved rudder through ITS stop consistently.
  const warpOver = computePose({ ...base, warpDeg: WARP_LIMIT_DEG + 1 });
  assert.equal(warpOver.clamped, true);
  assert.ok(warpOver.warpTipRad <= (WARP_LIMIT_DEG * Math.PI) / 180 + 1e-12);
});

test("props counter-rotate (crossed chain) and cradle tracks warp", () => {
  const p = computePose({ ...base, propAngleRad: 1.25, warpDeg: WARP_LIMIT_DEG });
  assert.equal(p.leftPropRad, 1.25);
  assert.equal(p.rightPropRad, -1.25);
  assert.ok(Math.abs(p.cradleOffsetM - 0.12) < 1e-12, "full warp = full cradle throw");
});

test("schematic-preview flags at the ±25% boundary and one step past", () => {
  const atBoundary = {
    ...REFERENCE_DIMS,
    span_m: REFERENCE_DIMS.span_m * (1 + SCHEMATIC_PREVIEW_FRACTION),
  };
  assert.equal(computePose(base, atBoundary).schematicPreview, false, "AT the cap: not flagged");
  const past = { ...REFERENCE_DIMS, span_m: REFERENCE_DIMS.span_m * (1 + SCHEMATIC_PREVIEW_FRACTION) * (1 + 1e-9) };
  assert.equal(computePose(base, past).schematicPreview, true, "past the cap: flagged");
  const shrunk = { ...REFERENCE_DIMS, rudder_area_m2: REFERENCE_DIMS.rudder_area_m2 * 0.7 };
  assert.equal(computePose(base, shrunk).schematicPreview, true, "shrink flags too");
});

test("non-finite controls refuse", () => {
  assert.throws(() => computePose({ ...base, warpDeg: Number.NaN }), /finite/);
  assert.throws(() => computePose({ ...base, propAngleRad: Infinity }), /finite/);
});
