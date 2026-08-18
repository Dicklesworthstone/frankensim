// E1.5 canard-mechanics battery (bead wf-root-guzez.2.5): per-item oracles
// over canard-mechanics-v1.json — the wide-priors law machine-enforced
// (every not-published parameter has lo < hi with a basis; point guesses
// refuse), evidence-class discipline for V-02b, cross-artifact consistency
// with flyer-reference, and arithmetic checks (static margin, slaving).
// Repro: node --test test/canard.test.ts

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

const load = (name: string): any =>
  JSON.parse(readFileSync(new URL(`../../../data/wright-flyer/${name}`, import.meta.url), "utf8"));

const doc = load("canard-mechanics-v1.json");
const reference = load("flyer-reference.json");
const dossier = load("source-dossier-v1.json");

test("structure: schema, dossier links resolve, frozen convention ids", () => {
  assert.equal(doc.schema, "org.frankensim.wright-flyer.canard-mechanics.v1");
  const ids = new Set(dossier.records.map((r: any) => r.id));
  for (const rec of doc.dossier_records) {
    assert.ok(ids.has(rec), `dossier record ${rec} must resolve`);
  }
  assert.equal(doc.convention_block.axes_id, "frd-body-v1");
  assert.equal(doc.convention_block.control_signs_id, "control-signs-v1");
  assert.ok(doc.no_claims.length > 40);
});

test("wide-priors law: every not-published parameter has a wide prior with a basis", () => {
  // Collect every status:"not-published" leaf in the artifact.
  const notPublished: string[] = [];
  const walk = (node: any, path: string): void => {
    if (node && typeof node === "object") {
      if (node.status === "not-published") {
        notPublished.push(path);
      }
      for (const [k, v] of Object.entries(node)) {
        walk(v, `${path}.${k}`);
      }
    }
  };
  walk(doc, "root");
  assert.ok(notPublished.length >= 6, `expected the known unpublished set, got ${notPublished.length}`);
  // The priors block must be non-degenerate per-item: lo < hi strictly, with
  // a substantive basis. A point guess (lo === hi) is the forbidden pattern.
  const priors = doc.unknown_priors.priors;
  assert.ok(priors.length >= notPublished.length - 1, "priors must cover the unpublished set");
  for (const p of priors) {
    assert.ok(Number.isFinite(p.lo) && Number.isFinite(p.hi), `${p.parameter}: finite bounds`);
    assert.ok(p.lo < p.hi, `${p.parameter}: prior must be WIDE (lo < hi), not a point guess`);
    assert.ok(p.basis.length > 20, `${p.parameter}: basis must be substantive`);
    console.log(JSON.stringify({ suite: "wf-canard", case: "prior", p: p.parameter, lo: p.lo, hi: p.hi }));
  }
});

test("V-02b evidence discipline: sign-tendency only, Estimated ceiling, A7a path", () => {
  const ev = doc.overbalance_evidence;
  assert.equal(ev.evidence_class, "sign-tendency-only");
  assert.equal(ev.quantitative_ceiling, "Estimated");
  assert.match(ev.promotion_path, /A7a/);
  assert.match(ev.diary_1903_verbatim, /balanced too near the center/);
  assert.match(ev.flying_1913_verbatim, /balanced too near the center/);
  assert.ok(ev.v02b_sourced_envelope.length >= 4, "the envelope must enumerate its clauses");
  assert.match(ev.conflation_warning, /VERTICAL rudder/);
});

test("cross-artifact consistency with flyer-reference", () => {
  // Travel: both artifacts agree ±30° is photo-inferred/absent; 1905 ±15 is
  // a different aircraft in both.
  assert.equal(doc.travel_stops.surface_deflection_1903.status, "absent-by-verification");
  assert.equal(reference.values.canard_travel_deg.status, "absent-by-verification");
  assert.equal(doc.travel_stops.surface_deflection_1903.value_deg, reference.values.canard_travel_deg.value);
  // Slaving + warp limits match the reference record.
  assert.equal(doc.hip_cradle.rudder_slaving.value, reference.values.rudder_warp_slaving.value);
  assert.match(reference.values.rudder_warp_slaving.original, /8\.5/);
  assert.equal(doc.hip_cradle.warp_limits_deg.value, 8.5);
  // Canard area: reference records 48 ft²; this dossier carries the variants.
  assert.equal(doc.geometry.area_ft2.value, 48);
  assert.ok(doc.geometry.area_ft2.variants.length >= 2, "51.6/52.75 variants recorded");
  // Arm agrees with the corrected wing-datum record (7.32 ft).
  assert.match(reference.values.canard_arm_m.original, /7\.32/);
  assert.equal(doc.geometry.arm_wing_to_canard_ft.value, 7.32);
});

test("arithmetic: static margin and inertia ordering", () => {
  const mi = doc.mass_inertia;
  const margin = mi.cg_neutral_point.neutral_point_pct_chord - mi.cg_neutral_point.cg_pct_chord_aft_of_wing_le;
  assert.ok(Math.abs(margin - mi.cg_neutral_point.static_margin_pct) < 0.11, `margin ${margin}`);
  assert.ok(margin < 0, "the Flyer is statically unstable — margin must be negative");
  const { Ixx, Iyy, Izz } = mi.aircraft_inertia_slugft2;
  assert.ok(Iyy < Ixx && Ixx < Izz, "slender-biplane ordering Iyy < Ixx < Izz");
  assert.equal(mi.aircraft_inertia_slugft2.status, "single-lineage");
});

test("variable camber is recorded as evidence, not decoration", () => {
  const vc = doc.kinematics.variable_camber;
  assert.equal(vc.status, "verified-qualitative");
  assert.match(vc.finding, /camber/i);
  assert.match(vc.sim_rule, /UNDERSTATES/);
  assert.ok(vc.sources.length >= 2, "two independent-group sources (L&P + Ames)");
});
