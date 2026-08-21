// E8.3b battery (bead wf-root-guzez.9.4): the scorecard-to-badge
// bridge over the SIX REAL registered WF rows (digests exactly as
// the frozen corpus carries them). Every row with a receipt surfaces
// Evidenced with its link; the stripped-digest twin surfaces
// "Estimated" — never blank, never a link; interpolation-vs-
// extrapolation is LIVE at the context boundary (AT = interpolated,
// past = EXTRAPOLATED in the sentence); missing coordinates refuse;
// caps at cap AND cap+1.
// Repro: node --test test/scorecardBridge.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_ROWS,
  bridgeScorecard,
  type WfScorecardRow,
} from "../src/scorecardBridge.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-scorecardbridge","case":"${kase}",${payload}}`);
}

/** The six registered rows, digests verbatim from the frozen corpus. */
function registeredRows(): WfScorecardRow[] {
  return [
    {
      datasetId: "wf-a-image-onplane-residual",
      metric: "wall-normal-residual-rel",
      receiptDigest: "cebf414b1ba1b5086b71afb372ab0b3f8bebf39f056e066d34f334b3d827f503",
      contextName: "ground-clearance-m",
      contextLo: 0.5,
      contextHi: 3.0,
    },
    {
      datasetId: "wf-a-wakeref-wagner-start",
      metric: "wagner-start-ratio",
      receiptDigest: "289bbe393d8b79dfddcbc92becfe4e42b0eb1a501f9561601a165926a21ce1f1",
      contextName: "freestream-mps",
      contextLo: 13,
      contextHi: 13,
    },
    {
      datasetId: "wf-a-farfield-v10-shape",
      metric: "v10-shape-rms",
      receiptDigest: "bc6c175b04fe249f5747bd286aa11a994f02c696bbfd524a9730e6d6a073b448",
      contextName: "overlap-ticks",
      contextLo: 237,
      contextHi: 237,
    },
    {
      datasetId: "wf-a-blade-cover-capsules",
      metric: "capsules-per-blade",
      receiptDigest: "6a7422ddc89749610ebc58a0aad41fb75522dd58362e0be2efc73ea4d1d7c28d",
      contextName: "geometry-uncertainty-m",
      contextLo: 0.01,
      contextHi: 0.01,
    },
    {
      datasetId: "wf-a-harness-worst-rel",
      metric: "harness-worst-abs-rel",
      receiptDigest: "289bbe393d8b79dfddcbc92becfe4e42b0eb1a501f9561601a165926a21ce1f1",
      contextName: "pinned-alpha-rad",
      contextLo: 0.03,
      contextHi: 0.07,
    },
    {
      datasetId: "wf-a-h07-slope-recovery",
      metric: "posterior-slope-recovery",
      receiptDigest: "b17cb8e3620e8fc8c7134cec9a4176d5603ad3a3a8708edae987c1ed18bc8d46",
      contextName: "lofo-folds",
      contextLo: 4,
      contextHi: 4,
    },
  ];
}

const IN_DOMAIN_POINT = {
  "ground-clearance-m": 1.5,
  "freestream-mps": 13,
  "overlap-ticks": 237,
  "geometry-uncertainty-m": 0.01,
  "pinned-alpha-rad": 0.05,
  "lofo-folds": 4,
};

test("every registered row surfaces Evidenced with its actual receipt", () => {
  const r = bridgeScorecard(registeredRows(), IN_DOMAIN_POINT);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.equal(r.value.length, 6);
  for (const row of r.value) {
    assert.equal(row.standing, "evidenced", row.datasetId);
    assert.equal(row.query, "interpolated", row.datasetId);
    assert.match(row.sentence, /evidenced by frozen-registry receipt/);
    assert.ok(row.receiptLink?.startsWith("receipt:"), row.datasetId);
  }
  // The links carry the EXACT frozen digests (per-row oracle).
  assert.equal(
    r.value[0]?.receiptLink,
    "receipt:cebf414b1ba1b5086b71afb372ab0b3f8bebf39f056e066d34f334b3d827f503",
  );
  jlog("evidenced", `"rows":6`);
});

test("no receipt -> Estimated, never blank, never a link", () => {
  const stripped = registeredRows().map((row) => {
    const { receiptDigest: _omitted, ...rest } = row;
    return rest as WfScorecardRow;
  });
  const r = bridgeScorecard(stripped, IN_DOMAIN_POINT);
  assert.ok(r.ok);
  if (!r.ok) return;
  for (const row of r.value) {
    assert.equal(row.standing, "estimated", row.datasetId);
    assert.match(row.sentence, /Estimated — no receipt in the frozen registry/);
    assert.ok(row.sentence.length > 20, "never blank");
    assert.equal(row.receiptLink, null, "never a link without a receipt");
  }
  jlog("estimated", `"rows":6`);
});

test("interpolation-vs-extrapolation is LIVE at the context boundary", () => {
  const rows = registeredRows();
  // AT the boundary: interpolated (closed interval).
  const at = bridgeScorecard(rows, { ...IN_DOMAIN_POINT, "ground-clearance-m": 3.0 });
  assert.ok(at.ok && at.value[0]?.query === "interpolated");
  // Past the boundary: EXTRAPOLATED, and the sentence says so.
  const past = bridgeScorecard(rows, { ...IN_DOMAIN_POINT, "ground-clearance-m": 3.000001 });
  assert.ok(past.ok);
  if (past.ok) {
    assert.equal(past.value[0]?.query, "extrapolated");
    assert.match(past.value[0]?.sentence ?? "", /extrapolated/);
    // Still evidenced — the receipt exists; only the QUERY standing
    // changed (the two statuses are independent, both displayed).
    assert.equal(past.value[0]?.standing, "evidenced");
  }
  // Below the low bound too.
  const low = bridgeScorecard(rows, { ...IN_DOMAIN_POINT, "pinned-alpha-rad": 0.02 });
  assert.ok(low.ok && low.value[4]?.query === "extrapolated");
  jlog("live-standing", `"at":"interpolated","past":"extrapolated"`);
});

test("refusals and caps", () => {
  // Missing operating coordinate refuses (a badge computed at
  // nowhere is a lie).
  const { "lofo-folds": _dropped, ...partial } = IN_DOMAIN_POINT;
  const missing = bridgeScorecard(registeredRows(), partial);
  assert.ok(!missing.ok && missing.refusal.code === "bridge-point-invalid");
  // Malformed digest refuses.
  const bad = registeredRows();
  const withBad: WfScorecardRow[] = [{ ...bad[0]!, receiptDigest: "nope" }];
  const badR = bridgeScorecard(withBad, IN_DOMAIN_POINT);
  assert.ok(!badR.ok && badR.refusal.code === "bridge-digest-malformed");
  // Caps: 64 rows admits, 65 refuses; empty refuses; malformed row.
  const mk = (n: number): WfScorecardRow[] =>
    Array.from({ length: n }, (_, i) => ({
      datasetId: `d${i}`,
      metric: "m",
      contextName: "x",
      contextLo: 0,
      contextHi: 1,
    }));
  assert.ok(bridgeScorecard(mk(MAX_ROWS), { x: 0.5 }).ok, "AT cap");
  const over = bridgeScorecard(mk(MAX_ROWS + 1), { x: 0.5 });
  assert.ok(!over.ok && over.refusal.code === "bridge-rows-invalid");
  const empty = bridgeScorecard([], { x: 0.5 });
  assert.ok(!empty.ok && empty.refusal.code === "bridge-rows-invalid");
  const disordered = bridgeScorecard(
    [{ datasetId: "d", metric: "m", contextName: "x", contextLo: 2, contextHi: 1 }],
    { x: 0.5 },
  );
  assert.ok(!disordered.ok && disordered.refusal.code === "bridge-row-malformed");
  jlog("caps", `"max_rows":${MAX_ROWS}`);
});
