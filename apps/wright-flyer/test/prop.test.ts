// E1.6 prop-package battery (bead wf-root-guzez.2.6): per-item oracles over
// prop-geometry-v1.json — 1911 station-table physical consistency, static
// anchor arithmetic, J-arithmetic re-derived against flyer-reference,
// calibration/holdout partition discipline, and the 1903-absence honesty
// structure. Repro: node --test test/prop.test.ts

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

const load = (name: string): any =>
  JSON.parse(readFileSync(new URL(`../../../data/wright-flyer/${name}`, import.meta.url), "utf8"));

const doc = load("prop-geometry-v1.json");
const reference = load("flyer-reference.json");
const dossier = load("source-dossier-v1.json");

test("artifact structure: schema, dossier link, convention block, no_claims", () => {
  assert.equal(doc.schema, "org.frankensim.wright-flyer.prop-geometry.v1");
  assert.ok(dossier.records.some((r: any) => r.id === doc.dossier_record), "dossier link resolves");
  assert.equal(doc.convention_block.axes_id, "frd-body-v1");
  assert.ok(doc.no_claims.length > 40);
  for (const key of ["CT", "CP", "J", "eta"]) {
    assert.ok(doc.coefficient_conventions[key].length > 5, `${key} convention missing`);
  }
});

test("1903 radial geometry is absent-by-verification with real reconstruction paths", () => {
  const r1903 = doc.radial_geometry_1903;
  assert.equal(r1903.status, "absent-by-verification");
  assert.ok(r1903.reconstruction_paths.length >= 3, "must name the reconstruction options");
  for (const p of r1903.reconstruction_paths) {
    assert.ok(p.evidence_ceiling.length > 0, "every path declares its evidence ceiling");
  }
  assert.match(r1903.sim_rule, /never be promoted/, "no promotion-by-1911-agreement");
  assert.match(doc.radial_geometry_1911_bentend.warning, /NOT the 1903/);
});

test("1911 station table is physically consistent per-row", () => {
  const st = doc.radial_geometry_1911_bentend.stations;
  assert.equal(st.length, 8, "eight published stations (0.96/1.0 reuse the 0.9 section)");
  const R_in = 51.0;
  let prevR = 0;
  let prevPhi = 90;
  for (const row of st) {
    assert.ok(row.r_over_R > prevR && row.r_over_R <= 1.0, "r/R strictly ascending in (0,1]");
    // Twist must decrease outboard (helix geometry) — per-row oracle.
    assert.ok(row.phi_deg < prevPhi, `phi must decrease outboard at r/R=${row.r_over_R}`);
    // Mean thickness from area/chord stays in a plausible wooden-blade band.
    const tMean = row.area_in2 / row.chord_in;
    assert.ok(tMean > 0.2 && tMean < 1.0, `mean thickness ${tMean} in implausible`);
    // Perimeter of a thin section is a bit over twice the chord.
    assert.ok(
      row.perimeter_in > 2 * row.chord_in && row.perimeter_in < 2.6 * row.chord_in,
      `perimeter/chord ratio off at r/R=${row.r_over_R}`,
    );
    prevR = row.r_over_R;
    prevPhi = row.phi_deg;
  }
  // Geometric pitch 2*pi*r*tan(phi): a real Wright blade is near-constant-
  // pitch over the working span — check spread over r/R 0.4..0.9 is < 25%.
  const pitches = st
    .filter((r: any) => r.r_over_R >= 0.4)
    .map((r: any) => 2 * Math.PI * r.r_over_R * R_in * Math.tan((r.phi_deg * Math.PI) / 180));
  const maxP = Math.max(...pitches);
  const minP = Math.min(...pitches);
  assert.ok((maxP - minP) / maxP < 0.25, `pitch spread ${(maxP - minP) / maxP} too wide`);
  console.log(
    JSON.stringify({ suite: "wf-prop", case: "stations", pitch_in_range: [minP, maxP] }),
  );
});

test("static anchors: pair-to-per-prop arithmetic and repro agreement band", () => {
  // 132-136 lb pair => 66-68 per prop; LFST repro 64.2 within 6% of the low end.
  const perPropLo = 132 / 2;
  const perPropHi = 136 / 2;
  assert.ok(perPropLo === 66 && perPropHi === 68);
  const repro = 64.2;
  assert.ok(Math.abs(repro - perPropLo) / perPropLo < 0.06, "repro vs Wrights band");
  assert.match(doc.performance_fixtures.static_1903.wrights_bench, /132-136/);
  assert.match(doc.performance_fixtures.static_1903.lfst_reproduction, /64\.2/);
});

test("1911 eta curve: monotone to the peak, peak 0.87 at J=1.15, then falls", () => {
  const eta = doc.performance_fixtures.curves_1911_bentend.eta_vs_J as [number, number][];
  const peak = eta.reduce((a, b) => (b[1] > a[1] ? b : a));
  assert.deepEqual(peak, [1.15, 0.87]);
  const peakIdx = eta.findIndex((p) => p[0] === 1.15);
  for (let i = 1; i <= peakIdx; i++) {
    assert.ok(eta[i]![1] > eta[i - 1]![1], `eta must rise to the peak at J=${eta[i]![0]}`);
  }
  for (let i = peakIdx + 1; i < eta.length; i++) {
    assert.ok(eta[i]![1] < eta[peakIdx]![1], "eta must fall past the peak");
  }
  // CP monotone decreasing over the published range (per-row).
  const cp = doc.performance_fixtures.curves_1911_bentend.cp_vs_J as [number, number][];
  for (let i = 1; i < cp.length; i++) {
    assert.ok(cp[i]![1] < cp[i - 1]![1], `CP must decrease at J=${cp[i]![0]}`);
  }
});

test("J arithmetic re-derived against flyer-reference values", () => {
  const D = doc.overall_geometry.diameter_m.value;
  assert.ok(Math.abs(D - reference.values.prop_diameter_m.value) / D < 0.005, "diameter consistent");
  const n350 = 350 / 60;
  const jRail = reference.values.wind_dec17_mps.value / (n350 * D);
  assert.ok(Math.abs(jRail - doc.operating_points.J_rail_dec17.value) < 0.005, `J_rail ${jRail}`);
  const jLift = reference.values.airspeed_dec17_mps.value / (n350 * D);
  assert.ok(Math.abs(jLift - doc.operating_points.J_liftoff_dec17.value) < 0.005, `J_lift ${jLift}`);
  const jDesign = reference.values.wind_dec17_mps.value / ((330 / 60) * D);
  assert.ok(Math.abs(jDesign - doc.operating_points.J_design_wilbur.value) < 0.005);
  // The plan's rail band 0.7-0.8 must contain J_rail.
  assert.ok(jRail > 0.7 && jRail < 0.8);
  console.log(JSON.stringify({ suite: "wf-prop", case: "J", jRail, jLift, jDesign }));
});

test("partition discipline: 1903 anchors are holdout, 1911 curves calibration", () => {
  const fx = doc.performance_fixtures;
  assert.equal(fx.static_1903.partition, "holdout");
  assert.equal(fx.static_1904.partition, "holdout");
  assert.equal(fx.eta_band_1903.partition, "holdout");
  assert.equal(fx.curves_1911_bentend.partition, "calibration");
  assert.match(fx.eta_band_1903.sim_rule, /forbidden/);
  assert.match(fx.eta_band_1903.lineage, /SECONDARY/);
  // Every fixture must carry a partition label (per-item, no silent additions).
  for (const [name, f] of Object.entries(fx)) {
    assert.ok(
      ["calibration", "holdout"].includes((f as any).partition),
      `${name} missing partition label`,
    );
  }
});
