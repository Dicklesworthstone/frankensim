// E1.8 air-state battery (bead wf-root-guzez.2.8): re-computes every DERIVED
// value from the recorded measurements (rho via ideal gas, mu via Sutherland,
// Re from rho*V*c/mu), enforces the WindReference schema on every record,
// verifies the ensemble pre-registration is wide distributions (never point
// guesses) and that NO gust trace can ever live in this artifact.
// Repro: node --test test/air.test.ts

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

const load = (name: string): any =>
  JSON.parse(readFileSync(new URL(`../../../data/wright-flyer/${name}`, import.meta.url), "utf8"));

const doc = load("air-state-v1.json");
const reference = load("flyer-reference.json");
const dossier = load("source-dossier-v1.json");

test("structure: schema, dossier links, wind-reference convention binding", () => {
  assert.equal(doc.schema, "org.frankensim.wright-flyer.air-state.v1");
  const ids = new Set(dossier.records.map((r: any) => r.id));
  for (const rec of doc.dossier_records) {
    assert.ok(ids.has(rec), `${rec} must resolve in the source dossier`);
  }
  assert.equal(doc.convention_block.wind_reference, "wind-reference-v1 records below");
  assert.ok(doc.no_claims.length > 40);
});

test("derived air state re-computed from the recorded measurements", () => {
  // rho = p/(R*T): 30.1 inHg -> Pa, 34 F -> K.
  const p_pa = 30.1 * 3386.39;
  const t_k = ((34 - 32) * 5) / 9 + 273.15;
  const rho = p_pa / (287.05 * t_k);
  assert.ok(Math.abs(rho - doc.derived_air_state.rho_kg_m3.value) / rho < 0.003, `rho ${rho}`);
  // slug/ft3 conversion (1 slug/ft3 = 515.3788 kg/m3).
  assert.ok(Math.abs(rho / 515.3788 - doc.derived_air_state.rho_kg_m3.slug_ft3) < 2e-5);
  // The cold-day density must exceed standard sea level by ~5-7%.
  assert.ok(rho / 1.225 > 1.04 && rho / 1.225 < 1.08, "density advantage band");
  // Sutherland viscosity at T.
  const mu = (1.458e-6 * Math.pow(t_k, 1.5)) / (t_k + 110.4);
  assert.ok(Math.abs(mu - doc.derived_air_state.mu_kg_m_s.value) / mu < 0.005, `mu ${mu}`);
  assert.ok(Math.abs(mu / rho - doc.derived_air_state.nu_m2_s.value) / (mu / rho) < 0.01);
  // Reynolds band from flyer-reference chord and the stated speed range.
  const c = reference.values.chord_m.value;
  const [reLo, reHi] = doc.derived_air_state.reynolds_flight.value_band;
  const reAt = (v: number): number => (rho * v * c) / mu;
  assert.ok(reAt(13.4) > reLo * 0.98 && reAt(15.2) < reHi * 1.02, `Re ${reAt(13.9)}`);
  console.log(JSON.stringify({ suite: "wf-air", case: "derived", rho, mu, re_139: reAt(13.9) }));
});

test("every WindReference record carries the frozen schema with priors for nulls", () => {
  assert.ok(doc.wind_reference_records.length >= 2, "both 1903 instruments");
  for (const rec of doc.wind_reference_records) {
    for (const field of ["instrument", "height_m_or_null", "averaging_interval_s_or_null", "provenance"]) {
      assert.ok(field in rec, `record missing ${field}`);
    }
    if (rec.height_m_or_null === null) {
      assert.ok(
        Array.isArray(rec.height_prior_m) && rec.height_prior_m[0] < rec.height_prior_m[1],
        `${rec.instrument}: null height must carry a wide prior`,
      );
    }
  }
  const instruments = doc.wind_reference_records.map((r: any) => r.instrument);
  assert.ok(instruments.includes("wrights-hand-richard"));
  assert.ok(instruments.includes("government-weather-bureau"));
});

test("negative findings and attributions are recorded honestly", () => {
  assert.equal(doc.dec17_1903_records.wind.bureau_ledger.status, "verified-negative");
  assert.match(doc.dec17_1903_records.wind.bureau_ledger.finding, /No Record/);
  assert.match(doc.dec17_1903_records.temperature.ice_on_puddles.note, /1913/);
  assert.equal(doc.dec17_1903_records.station_metadata.anemometer.status, "not-published");
  assert.match(doc.derived_air_state.reynolds_flight.literature_note, /order-of-magnitude convention/);
});

test("ensemble pre-registration: wide distributions, no gust trace anywhere", () => {
  const pre = doc.dec17_ensemble_preregistration;
  assert.ok(pre.parameters.length >= 6, "the parameter set must be complete");
  for (const p of pre.parameters) {
    assert.equal(p.wide, true, `${p.name} must be declared wide`);
    assert.ok(p.distribution.length > 10, `${p.name} needs a real distribution`);
  }
  assert.match(pre.doctrine, /never a claimed reconstruction/i);
  assert.match(pre.gust_realization_rule, /PhysicalRealizationAlgorithmId/);
  // THE HARD GATE: no key anywhere in the artifact may look like a gust
  // trace/time-series payload (arrays of numbers under gust-ish names).
  const offenders: string[] = [];
  const walk = (node: any, path: string): void => {
    if (node && typeof node === "object") {
      for (const [k, v] of Object.entries(node)) {
        if (/gust_?(trace|series|sequence)|wind_?(trace|series)/i.test(k) && Array.isArray(v)) {
          offenders.push(`${path}.${k}`);
        }
        walk(v, `${path}.${k}`);
      }
    }
  };
  walk(doc, "root");
  assert.deepEqual(offenders, [], `gust traces are forbidden here: ${offenders.join(", ")}`);
});

test("stability declaration and surface priors are declared axes, not points", () => {
  assert.equal(doc.neutral_stability_declaration.declared_class, "neutral");
  assert.match(doc.neutral_stability_declaration.caveat, /uncertainty axis/);
  const z0 = doc.surface.z0_m_prior;
  assert.ok(z0.lo < z0.hi && z0.lo > 0, "z0 prior wide and positive");
  // Consistent with the flyer-reference roughness band (overlapping ranges).
  assert.ok(z0.hi >= 1e-3 && z0.lo <= 1e-3, "must overlap the reference band");
  const d = doc.surface.displacement_height_m_prior;
  assert.ok(d.lo < d.hi, "displacement prior wide");
});
