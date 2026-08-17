// E1.2 reference battery (bead wf-root-guzez.2.2): per-value oracles over
// flyer-reference.json — unit-conversion arithmetic re-derived independently,
// AR/gross-mass consistency, cross-file dossier linkage, and the demotion
// doctrine (unverifiable values must be tunable-with-provenance or
// absent-by-verification, never silently kept as verified).

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

interface RefValue {
  value: number | string;
  original: string;
  status: string;
  sources: string[];
  dossier_record: string;
  convention: string;
  variants: { value_original: string; source: string; note?: string }[];
}

const load = (name: string): unknown =>
  JSON.parse(readFileSync(new URL(`../../../data/wright-flyer/${name}`, import.meta.url), "utf8"));

const doc = load("flyer-reference.json") as { schema: string; values: Record<string, RefValue>; no_claims: string };
const dossier = load("source-dossier-v1.json") as { records: { id: string }[] };
const dossierIds = new Set(dossier.records.map((r) => r.id));

const STATUSES = new Set([
  "verified",
  "verified-range",
  "tunable-with-provenance",
  "absent-by-verification",
]);

const num = (id: string): number => {
  const v = doc.values[id]!.value;
  assert.equal(typeof v, "number", `${id} not numeric`);
  return v as number;
};

test("every value carries status, sources, convention, and a real dossier link", () => {
  assert.equal(doc.schema, "org.frankensim.wright-flyer.flyer-reference.v1");
  for (const [id, v] of Object.entries(doc.values)) {
    assert.ok(STATUSES.has(v.status), `${id}: bad status ${v.status}`);
    assert.ok(v.sources.length >= 1, `${id}: no sources`);
    assert.ok(v.convention.length > 0, `${id}: no convention label`);
    assert.ok(v.original.length > 0, `${id}: no original-units record`);
    assert.ok(dossierIds.has(v.dossier_record), `${id}: dossier record ${v.dossier_record} not in source-dossier-v1`);
    console.log(JSON.stringify({ suite: "wf-reference", case: "value", id, status: v.status }));
  }
  assert.ok(Object.keys(doc.values).length >= 30, "expected the full §3 sweep");
});

test("unit conversions re-derived independently (ft/lb/hp -> SI)", () => {
  const FT = 0.3048;
  const LB = 0.45359237;
  const close = (id: string, expect: number, tolFrac = 0.002): void => {
    const got = num(id);
    assert.ok(
      Math.abs(got - expect) / expect < tolFrac,
      `${id}: ${got} vs derived ${expect.toFixed(4)}`,
    );
  };
  close("wingspan_m", (40 + 4 / 12) * FT);
  close("chord_m", 6.5 * FT);
  close("wing_area_both_m2", 510 * FT * FT);
  close("empty_weight_kg", 605 * LB);
  close("gross_weight_kg", 750 * LB);
  close("prop_diameter_m", 8.5 * FT);
  close("launch_rail_length_m", 60 * FT);
  close("derrick_height_m", 20 * FT);
  close("canard_arm_m", 7.32 * FT);
  close("canard_area_m2", 48 * FT * FT);
  close("rudder_area_m2", 20 * FT * FT);
  close("flight1_distance_m", 120 * FT);
  close("flight2_distance_m", 175 * FT);
  close("flight3_distance_m", 200 * FT);
  close("flight4_distance_m", 852 * FT);
  close("engine_displacement_m3", 201 * 0.0254 ** 3);
  close("engine_power_w", 12 * 745.7);
  close("engine_mass_kg", 180 * LB);
  close("chain_ratio", 23 / 8);
  close("airspeed_dec17_mps", 31 * 0.44704);
  close("wind_dec17_mps", 24 * 0.44704);
});

test("derived consistency: AR pair and gross-mass build-up", () => {
  const b = num("wingspan_m");
  const s = num("wing_area_both_m2");
  assert.ok(Math.abs(num("aspect_ratio_system") - (b * b) / s) < 0.02, "AR_system != b²/S_both");
  assert.ok(Math.abs(num("aspect_ratio_plane") - (b * b) / (s / 2)) < 0.04, "AR_plane != b²/S_one");
  const pilotKg = 145 * 0.45359237;
  assert.ok(
    Math.abs(num("gross_weight_kg") - (num("empty_weight_kg") + pilotKg)) < 0.5,
    "gross != empty + 145 lb pilot",
  );
});

test("demotion doctrine: the flagged discrepancies are NOT marked verified", () => {
  const expected: Record<string, string> = {
    canard_travel_deg: "absent-by-verification",
    engine_mass_kg: "tunable-with-provenance",
    rudder_area_m2: "tunable-with-provenance",
    derrick_height_m: "tunable-with-provenance",
    flight3_distance_m: "tunable-with-provenance",
    wing_gap_m: "tunable-with-provenance",
    canard_arm_m: "tunable-with-provenance",
  };
  for (const [id, status] of Object.entries(expected)) {
    assert.equal(doc.values[id]!.status, status, `${id} must stay demoted`);
  }
  // Every non-verified numeric value must document WHY: variants or an
  // explicit convention caveat.
  for (const [id, v] of Object.entries(doc.values)) {
    if (v.status === "tunable-with-provenance") {
      assert.ok(
        v.variants.length > 0 ||
          /UNSTATED|unresolved|undocumented|NO source|NO distance|later account/i.test(
            v.convention + v.sources.join(" "),
          ),
        `${id}: demoted without documented disagreement`,
      );
    }
  }
});

test("specific verified corrections are on the record", () => {
  assert.match(doc.values.canard_travel_deg!.sources.join(" "), /have not been reported/);
  assert.match(doc.values.derrick_height_m!.convention, /DROP DISTANCE/);
  assert.match(doc.values.prop_efficiency!.convention, /Wilbur-calculated/);
  assert.match(doc.values.flight1_distance_m!.sources.join(" "), /about 100 feet/);
  assert.match(doc.values.anhedral_intent!.convention, /PROMOTION/);
  assert.match(doc.values.wind_dec17_mps!.convention, /WindReference/);
  assert.ok(doc.no_claims.length > 40);
});
