// E7.4b battery (bead wf-root-guzez.8.6): the lateral view's laws on
// the LINKAGE-DECOUPLED fixture. 1901-mode adverse yaw REPRODUCED
// (sign law over commanded ticks; the coupled twin is NOT adverse —
// the rudder compensation is visible); broken decompositions refuse
// (never patched); spiral indicator vs closed form; A/B linkage
// admission (realization + prefix held fixed, falsifiers executed);
// non-vacuity refusal; caps at cap AND cap+1.
// Repro: node --test test/lateralView.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_WINDOW,
  admitLinkagePair,
  adverseYawVerdict,
  spiralIndicator,
  validateDecomposition,
  type YawDecomposition,
} from "../src/lateralView.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-lateral","case":"${kase}",${payload}}`);
}

/** The 1901-mode fixture: warp right, induced-drag yaw LEFT
 * (opposing), rudder silent (decoupled). */
function decoupledRows(n: number): YawDecomposition[] {
  return Array.from({ length: n }, (_, t) => {
    const warp = 0.08 * Math.sin(0.05 * t + 0.3);
    const induced = -14.0 * warp; // adverse: opposes the command
    const profile = 0.4;
    return {
      tick: t,
      warpCommandRad: warp,
      loadedTwistRad: 0.62 * warp, // aeroelastic loss, reported
      inducedDragYawNm: induced,
      rudderYawNm: 0,
      profileYawNm: profile,
      netYawNm: induced + 0 + profile,
    };
  });
}

/** The 1902+ coupled twin: the linked rudder overcompensates. */
function coupledRows(n: number): YawDecomposition[] {
  return decoupledRows(n).map((r) => {
    const rudder = -1.6 * r.inducedDragYawNm; // linkage compensation
    return {
      ...r,
      rudderYawNm: rudder,
      netYawNm: r.inducedDragYawNm + rudder + r.profileYawNm,
    };
  });
}

test("1901-mode adverse yaw REPRODUCED on the decoupled fixture; coupled twin is not adverse", () => {
  const decoupled = adverseYawVerdict(decoupledRows(200));
  assert.ok(decoupled.ok);
  if (decoupled.ok) {
    assert.equal(decoupled.value.adverse, true, "decoupled 1901 mode IS adverse");
    assert.ok(decoupled.value.meanSignProduct < -0.99);
    assert.ok(decoupled.value.commandedTicks > 150);
  }
  // The coupled twin: the induced component STILL opposes (physics
  // unchanged) but the NET yaw now follows the command — the view's
  // attribution shows the rudder doing it. Verdict on the induced
  // component stays adverse; the net-sign check flips.
  const coupled = coupledRows(200);
  const v = adverseYawVerdict(coupled);
  assert.ok(v.ok && v.value.adverse, "the induced component is adverse in BOTH riggings");
  let netFollows = 0;
  let commanded = 0;
  for (const r of coupled) {
    if (Math.abs(r.warpCommandRad) > 1e-6) {
      commanded += 1;
      if (Math.sign(r.netYawNm - r.profileYawNm) === Math.sign(r.warpCommandRad)) netFollows += 1;
    }
  }
  assert.ok(netFollows / commanded > 0.99, "linked rudder turns net yaw proverse");
  // Loaded twist is reported distinct from the command (aeroelastic
  // loss visible).
  const rows = decoupledRows(4);
  assert.ok(Math.abs((rows[1]?.loadedTwistRad ?? 0) / (rows[1]?.warpCommandRad ?? 1) - 0.62) < 1e-12);
  jlog("adverse-yaw", `"decoupled_adverse":true,"coupled_net_proverse":true`);
});

test("broken decompositions refuse, never patched", () => {
  const rows = decoupledRows(8);
  const broken = rows.map((r, i) => (i === 5 ? { ...r, netYawNm: r.netYawNm + 0.5 } : r));
  const v = validateDecomposition(broken);
  assert.ok(!v.ok && v.refusal.code === "lateral-decomposition-broken");
  assert.ok(!v.ok && v.refusal.message.includes("tick 5"), "first divergence localized");
  // Non-vacuity: a window with no commanded ticks refuses.
  const idle = rows.map((r) => ({ ...r, warpCommandRad: 0, inducedDragYawNm: 0, netYawNm: r.profileYawNm }));
  const nv = adverseYawVerdict(idle);
  assert.ok(!nv.ok && nv.refusal.code === "lateral-no-commanded-ticks");
  jlog("broken-split", `"code":"lateral-decomposition-broken"`);
});

test("spiral indicator matches closed form; caps at cap AND cap+1", () => {
  const div = spiralIndicator(0.09);
  assert.ok(div.ok);
  if (div.ok) {
    assert.ok(div.value.divergent);
    assert.ok(Math.abs((div.value.timeToDoubleS ?? 0) - Math.LN2 / 0.09) < 1e-15);
  }
  const conv = spiralIndicator(-0.02);
  assert.ok(conv.ok && !conv.value.divergent && conv.value.timeToDoubleS === null);
  const bad = spiralIndicator(Number.NaN);
  assert.ok(!bad.ok && bad.refusal.code === "lateral-pole-invalid");
  assert.ok(validateDecomposition(decoupledRows(MAX_WINDOW)).ok, "AT cap");
  const over = validateDecomposition(decoupledRows(MAX_WINDOW + 1));
  assert.ok(!over.ok && over.refusal.code === "lateral-window-invalid");
  jlog("spiral", `"t2":${Math.LN2 / 0.09}`);
});

test("A/B linkage admission holds realization + prefix fixed", () => {
  const coupled = { linkageCoupled: true, realizationId: "real-1", inputPrefixDigest: "p1" };
  const decoupled = { linkageCoupled: false, realizationId: "real-1", inputPrefixDigest: "p1" };
  const ok = admitLinkagePair(coupled, decoupled);
  assert.ok(ok.ok);
  if (ok.ok) {
    assert.equal(ok.value.coupled.linkageCoupled, true);
    assert.equal(ok.value.decoupled.linkageCoupled, false);
  }
  // Order-independent normalization.
  const swapped = admitLinkagePair(decoupled, coupled);
  assert.ok(swapped.ok && swapped.value.coupled.linkageCoupled);
  // FALSIFIERS: varied realization, varied prefix, degenerate pair.
  const r2 = admitLinkagePair(coupled, { ...decoupled, realizationId: "real-2" });
  assert.ok(!r2.ok && r2.refusal.code === "linkage-ab-realization-mismatch");
  const p2 = admitLinkagePair(coupled, { ...decoupled, inputPrefixDigest: "p2" });
  assert.ok(!p2.ok && p2.refusal.code === "linkage-ab-prefix-mismatch");
  const same = admitLinkagePair(coupled, { ...coupled });
  assert.ok(!same.ok && same.refusal.code === "linkage-ab-degenerate");
  const empty = admitLinkagePair({ ...coupled, realizationId: "" }, decoupled);
  assert.ok(!empty.ok && empty.refusal.code === "linkage-ab-realization-mismatch");
  jlog("ab-linkage", `"falsifiers":"executed"`);
});
