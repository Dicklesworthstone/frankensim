// E1.1 dossier battery (bead wf-root-guzez.2.1): per-record oracles over the
// evidence registry — never totals-only. Enforces the lineage doctrine: full
// EvidenceLineageId preimage on every record, independence-group discipline
// (derived sources share their origin's group; same-group corroboration is
// forbidden), the fitting block, and the plan's named claim prohibitions.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

interface DossierRecord {
  id: string;
  role: string;
  citations: string[];
  access: string;
  independence_group: string;
  derivation_steps: string[];
  permitted_claims: string[];
  forbidden_claims: string[];
  calibration_or_holdout: string;
  convention_uncertainty: string;
  unverified_residue: string[];
}

const doc = JSON.parse(
  readFileSync(
    new URL("../../../data/wright-flyer/source-dossier-v1.json", import.meta.url),
    "utf8",
  ),
) as { schema: string; attribution_correction: string; records: DossierRecord[]; no_claims: string };

const byId = new Map(doc.records.map((r) => [r.id, r]));
const get = (id: string): DossierRecord => {
  const r = byId.get(id);
  assert.ok(r, `record ${id} missing`);
  return r;
};

test("every record carries the complete lineage preimage (per-item, not totals)", () => {
  assert.equal(doc.schema, "org.frankensim.wright-flyer.source-dossier.v1");
  assert.equal(new Set(doc.records.map((r) => r.id)).size, doc.records.length, "duplicate ids");
  for (const r of doc.records) {
    for (const key of [
      "role",
      "access",
      "independence_group",
      "calibration_or_holdout",
      "convention_uncertainty",
    ] as const) {
      assert.ok(r[key].length > 0, `${r.id}.${key} empty`);
    }
    assert.ok(r.citations.length > 0, `${r.id}: no citations`);
    assert.ok(r.derivation_steps.length > 0, `${r.id}: no derivation steps`);
    assert.ok(r.forbidden_claims.length > 0, `${r.id}: every source has a claim boundary`);
    console.log(
      JSON.stringify({ suite: "wf-dossier", case: "lineage", id: r.id, group: r.independence_group }),
    );
  }
});

test("A2 is partitioned by origin into >=3 independence groups", () => {
  const a2 = doc.records.filter((r) => r.id.startsWith("a2-"));
  assert.ok(a2.length >= 4, "A2 family too small");
  const groups = new Set(a2.map((r) => r.independence_group));
  assert.ok(groups.size >= 3, `A2 groups ${[...groups].join(",")} — need >=3`);
});

test("derived source shares its origin's group and may not corroborate it", () => {
  const deters = get("a2-simmodels-deters");
  const ames = get("a2-ames-fullscale-1999");
  assert.equal(deters.independence_group, ames.independence_group, "derived => same group");
  assert.ok(
    deters.forbidden_claims.some((c) => c.toLowerCase().includes("corroboration")),
    "derived source must forbid corroborating its origin",
  );
});

test("synthesized stall may not validate Wright-specific deep stall (plan law)", () => {
  const synth = get("a2-synthesized-stall");
  assert.equal(synth.independence_group, "synthesized-internal");
  assert.ok(
    synth.forbidden_claims.some((c) => c.includes("deep stall")),
    "the deep-stall prohibition must be spelled on the record",
  );
  assert.ok(
    synth.forbidden_claims.some((c) => c.includes("Estimated")),
    "evidence-color ceiling must be declared",
  );
});

test("fitting block: no measured record is partitioned yet; anchors are holdout-only", () => {
  for (const r of doc.records) {
    assert.ok(r.calibration_or_holdout.length > 0);
    assert.ok(
      !["calibration", "both"].includes(r.calibration_or_holdout),
      `${r.id}: nothing may be marked fittable before E1.2 partitions it`,
    );
  }
  assert.match(get("a3-dec17-records").calibration_or_holdout, /holdout/);
  assert.match(get("a4-culick-stability").calibration_or_holdout, /holdout/);
});

test("future campaigns (A7a/A7b) carry zero permitted claims", () => {
  for (const id of ["a7a-canard-rig", "a7b-replica-fixture"]) {
    const r = get(id);
    assert.equal(r.permitted_claims.length, 0, `${id} claims before execution`);
    assert.match(r.calibration_or_holdout, /pre-registered/);
  }
});

test("correlated-analyst axis: Culick group discounts Ames corroboration explicitly", () => {
  const culick = get("a4-culick-stability");
  assert.notEqual(culick.independence_group, get("a2-ames-fullscale-1999").independence_group);
  assert.ok(
    culick.forbidden_claims.some((c) => c.includes("correlated-analyst")),
    "shared-personnel discount must be on the record",
  );
});

test("the attribution correction is present and repoints prop data", () => {
  assert.match(doc.attribution_correction, /2004-0211/);
  assert.match(doc.attribution_correction, /a2-props-bentend/);
  const props = get("a2-props-bentend");
  assert.ok(props.citations.some((c) => c.includes("2.2944")), "JoA prop citation");
  assert.ok(doc.no_claims.length > 40);
});
