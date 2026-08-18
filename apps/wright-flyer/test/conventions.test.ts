// E1.4 conventions battery (bead wf-root-guzez.2.4): machine checks over the
// frozen frame-conventions artifact — handedness computed from the declared
// unit vectors (not trusted as prose), cross-artifact id closure against
// geometry-conventions-v1 and flyer-reference, and the re-expression rule's
// required fields.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

const load = (name: string): any =>
  JSON.parse(readFileSync(new URL(`../../../data/wright-flyer/${name}`, import.meta.url), "utf8"));

const frame = load("frame-conventions-v1.json");
const geom = load("geometry-conventions-v1.json");
const reference = load("flyer-reference.json");

test("artifact is frozen, versioned, and carries no_claims", () => {
  assert.equal(frame.schema, "org.frankensim.wright-flyer.frame-conventions.v1");
  assert.equal(frame.status, "frozen-conventions");
  assert.ok(String(frame.no_claims).length > 20);
});

test("body axes are right-handed BY COMPUTATION from the declared vectors", () => {
  const [x, y, z] = [frame.body_axes.x_unit, frame.body_axes.y_unit, frame.body_axes.z_unit] as [
    number[],
    number[],
    number[],
  ];
  const cross = [
    x[1]! * y[2]! - x[2]! * y[1]!,
    x[2]! * y[0]! - x[0]! * y[2]!,
    x[0]! * y[1]! - x[1]! * y[0]!,
  ];
  assert.deepEqual(cross, z, "x cross y must equal z (right-handed FRD)");
  console.log(JSON.stringify({ suite: "wf-conventions", case: "handedness", cross }));
});

test("every frozen block carries a versioned id", () => {
  for (const key of [
    "body_axes",
    "moment_signs",
    "angle_conventions",
    "control_signs",
    "wind_reference",
    "uncertainty_doctrine",
    "dataset_reexpression_rule",
  ]) {
    assert.match(frame[key].id, /-v\d+$/, `${key} id must be versioned`);
  }
});

test("control signs are moment-defined and cover all three 1903 controls", () => {
  const cs = frame.control_signs;
  assert.match(cs.canard_dc, /POSITIVE pitch/);
  assert.match(cs.warp_dw, /POSITIVE roll/);
  assert.match(cs.rudder_dr, /POSITIVE yaw/);
  assert.match(cs.rudder_dr, /SLAVED/);
  assert.match(cs.doctrine, /COMMANDED MOMENT/);
});

test("wind reference disaggregates both 1903 instruments with the unrecorded-height rule", () => {
  const wr = frame.wind_reference;
  assert.ok(wr.instruments["wrights-hand-richard"]);
  assert.ok(wr.instruments["government-weather-bureau"]);
  assert.match(wr.unrecorded_rule, /null/);
  assert.deepEqual(wr.schema_fields, [
    "instrument",
    "height_m_or_null",
    "averaging_interval_s_or_null",
    "provenance",
  ]);
});

test("cross-artifact closure: geometry sibling frozen; reference labels resolve", () => {
  assert.equal(geom.status, "frozen-conventions");
  assert.ok(geom.aspect_ratio.AR_plane && geom.aspect_ratio.AR_system, "AR pair frozen in sibling");
  // flyer-reference convention labels must cite the frozen artifacts where
  // they lean on them (per-value oracle, not totals).
  for (const id of ["wingspan_m", "aspect_ratio_plane", "aspect_ratio_system", "wing_area_both_m2"]) {
    assert.match(
      reference.values[id].convention,
      /geometry-conventions-v1/,
      `${id} must cite the frozen geometry conventions`,
    );
  }
  assert.match(reference.values.wind_dec17_mps.convention, /WindReference/);
  assert.match(frame.control_signs.rudder_dr, /2\.5/, "slaving ratio consistent with flyer-reference");
});

test("re-expression rule names its required fields and the not-yet-ingested tables", () => {
  const rule = frame.dataset_reexpression_rule;
  assert.match(rule.rule, /axes_id/);
  assert.match(rule.rule, /typed refusal/);
  assert.match(rule.reexpressed_now, /not yet in-repo/, "honesty: tables outstanding, not claimed");
});
