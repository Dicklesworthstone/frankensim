// E0.9a schema-closure battery (bead wf-root-guzez.1.9.1): machine version of
// the Round-5/6 preimage audit over the FROZEN artifact — every identity name
// referenced inside a formula must itself be defined (formula, leaf closure,
// basis field, or explicit primitive), and domain separators must be unique
// and versioned.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

const schemaPath = new URL(
  "../../../data/wright-flyer/replay-identity-schema-v1.json",
  import.meta.url,
);
const doc = JSON.parse(readFileSync(schemaPath, "utf8")) as {
  schema: string;
  status: string;
  domain_separators: Record<string, string>;
  run_identity_basis_v1: { fields: string[] };
  formulas: Record<string, string>;
  leaf_identity_closures: Record<string, string>;
  replay_envelope_v1_fields: string[];
  checkpoint_state_v1_fields: string[];
  no_claims: string;
};

test("artifact parses, is FROZEN, and carries no_claims", () => {
  assert.equal(doc.schema, "org.frankensim.wright-flyer.replay-identity-schema.v1");
  assert.equal(doc.status, "FROZEN");
  assert.ok(doc.no_claims.length > 40, "no_claims must be substantive");
  console.log(JSON.stringify({ suite: "wf-identity-schema", case: "frozen", schema: doc.schema }));
});

test("domain separators are unique, versioned, and namespaced", () => {
  const seps = Object.values(doc.domain_separators);
  assert.equal(new Set(seps).size, seps.length, "separator collision");
  for (const s of seps) {
    assert.match(s, /^fs-flyer\/[a-z0-9-]+\/v\d+$/, `malformed separator ${s}`);
  }
});

test("preimage closure: every *Id referenced in a formula is defined somewhere", () => {
  const defined = new Set<string>([
    ...Object.keys(doc.formulas),
    ...Object.keys(doc.leaf_identity_closures),
    // Basis-field ids (snake_case → PascalCase) defined by the identity table:
    "PhysicalScenarioId",
    "ModelId",
    "ArtifactId",
    // Schema-id primitives (versioned literals, not derived hashes):
    "InputTraceV1",
    "RunIdentityBasisV1",
    "CheckpointStateV1",
  ]);
  const referenced = new Set<string>();
  for (const [name, text] of Object.entries(doc.formulas)) {
    for (const m of text.matchAll(/\b([A-Z][A-Za-z]*Id|[A-Z][A-Za-z]*V1)\b/g)) {
      if (m[1] !== name) {
        referenced.add(m[1]!);
      }
    }
  }
  const undefinedRefs = [...referenced].filter((r) => !defined.has(r));
  assert.deepEqual(undefinedRefs, [], `undefined identity references: ${undefinedRefs.join(", ")}`);
  console.log(
    JSON.stringify({ suite: "wf-identity-schema", case: "closure", referenced: referenced.size }),
  );
});

test("basis has exactly the seven Round-4 fields; RunSpecId excludes the tick-0 digest", () => {
  assert.equal(doc.run_identity_basis_v1.fields.length, 7);
  assert.ok(doc.run_identity_basis_v1.fields.includes("accepted_tick0_state_digest"));
  assert.match(doc.formulas.RunSpecId!, /excludes accepted_tick0_state_digest/);
  assert.match(doc.formulas.RunId!, /RunIntentId, InputTraceId/, "two-input RunId");
});

test("checkpoint carries the eight Round-3 algorithmic-history items", () => {
  const hist = doc.checkpoint_state_v1_fields.filter((f) => f.startsWith("algorithmic_history:"));
  assert.equal(hist.length, 8, `expected 8 algorithmic-history items, got ${hist.length}`);
  assert.ok(
    doc.replay_envelope_v1_fields.some((f) => f.includes("ATTACHMENT")),
    "acquisition trace must be labeled an attachment",
  );
});
