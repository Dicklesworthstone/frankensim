// E10.2-ii battery (bead wf-root-guzez.11.6.2): badge laws. A pass
// state REQUIRES a receipt digest (the green-without-evidence
// falsifier refuses, never silently demotes); unknown verdicts are
// NO-DATA never blank; passthrough verbatim; malformed digests
// refuse; caps at cap AND cap+1.
// Repro: node --test test/evidenceBadges.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_BADGES,
  buildEvidenceBadges,
  type ReceiptRow,
} from "../src/evidenceBadges.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-evidencebadges","case":"${kase}",${payload}}`);
}

const DIGEST = "b17cb8e3620e8fc8c7134cec9a4176d5603ad3a3a8708edae987c1ed18bc8d46";

function row(overrides: Partial<ReceiptRow>): ReceiptRow {
  return {
    caseId: "V-08b1",
    verdict: "pass",
    receiptDigest: DIGEST,
    comparisonClass: "formulation-band-0.15",
    ...overrides,
  };
}

test("real receipts render evidenced badges verbatim", () => {
  const r = buildEvidenceBadges([
    row({}),
    row({ caseId: "V-10", verdict: "reported-only", receiptDigest: DIGEST }),
    row({ caseId: "H-02c", verdict: "compatibility-check-cannot-promote", receiptDigest: undefined }),
  ]);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.equal(r.value[0]?.state, "evidenced-pass");
  assert.equal(r.value[0]?.receiptDigest, DIGEST);
  assert.equal(r.value[1]?.state, "reported");
  // Unknown verdict: NO-DATA, verbatim verdict preserved — never
  // blank, never green.
  assert.equal(r.value[2]?.state, "no-data");
  assert.equal(r.value[2]?.verdict, "compatibility-check-cannot-promote");
  assert.equal(r.value[2]?.receiptDigest, null);
  jlog("verbatim", `"badges":3`);
});

test("FALSIFIER: pass without a receipt digest refuses (never demotes silently)", () => {
  const forged = buildEvidenceBadges([row({ receiptDigest: undefined })]);
  assert.ok(!forged.ok);
  assert.equal(!forged.ok && forged.refusal.code, "badge-pass-without-receipt");
  jlog("green-without-evidence", `"code":"badge-pass-without-receipt"`);
});

test("malformed digests and missing case ids refuse", () => {
  const short = buildEvidenceBadges([row({ receiptDigest: "deadbeef" })]);
  assert.ok(!short.ok && short.refusal.code === "badge-digest-malformed");
  const upper = buildEvidenceBadges([row({ receiptDigest: DIGEST.toUpperCase() })]);
  assert.ok(!upper.ok && upper.refusal.code === "badge-digest-malformed");
  const unnamed = buildEvidenceBadges([row({ caseId: "  " })]);
  assert.ok(!unnamed.ok && unnamed.refusal.code === "badge-case-missing");
  jlog("malformed", `"refused":3`);
});

test("caps at cap AND cap+1", () => {
  const mk = (n: number) => Array.from({ length: n }, (_, i) => row({ caseId: `V-${i}` }));
  assert.ok(buildEvidenceBadges(mk(MAX_BADGES)).ok, "AT cap");
  const over = buildEvidenceBadges(mk(MAX_BADGES + 1));
  assert.ok(!over.ok && over.refusal.code === "badge-count-exceeded");
  jlog("caps", `"max":${MAX_BADGES}`);
});
