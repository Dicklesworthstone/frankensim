// E8.3a battery (bead wf-root-guzez.9.3): applicability plumbing.
// Intersection with per-BOUND limiting-subsystem attribution
// verified on boundary fixtures; empty intersections NAME the
// conflicting pair; AT-bound is inside (closed intervals, cap AND
// one-ulp-past behavior); badge composition laws (green is EARNED:
// evidenced AND inside — the evidenced-but-outside falsifier goes
// amber naming the limiter; NO-DATA is the honest gray empty state);
// receipt link only with a digest; caps at cap AND cap+1.
// Repro: node --test test/applicability.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_AXES,
  MAX_SUBSYSTEMS,
  composeBadge,
  intersectDomains,
  standingAt,
  type ApplicabilityDomain,
} from "../src/applicability.ts";
import type { EvidenceBadge } from "../src/evidenceBadges.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-applicability","case":"${kase}",${payload}}`);
}

const DIGEST = "b17cb8e3620e8fc8c7134cec9a4176d5603ad3a3a8708edae987c1ed18bc8d46";

/** Boundary fixture: aero binds the speed lo, structure the hi;
 * pilot binds alpha hi, aero the lo. */
function domains(): ApplicabilityDomain[] {
  return [
    {
      subsystem: "aero",
      axes: [
        { name: "speed", lo: 8.0, hi: 20.0 },
        { name: "alpha", lo: -0.05, hi: 0.3 },
      ],
    },
    {
      subsystem: "structure",
      axes: [
        { name: "speed", lo: 0.0, hi: 16.0 },
        { name: "alpha", lo: -0.2, hi: 0.4 },
      ],
    },
    {
      subsystem: "pilot",
      axes: [
        { name: "speed", lo: 5.0, hi: 18.0 },
        { name: "alpha", lo: -0.1, hi: 0.12 },
      ],
    },
  ];
}

test("intersection attributes each bound to its limiting subsystem", () => {
  const r = intersectDomains(domains());
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.equal(r.value.kind, "domain");
  if (r.value.kind !== "domain") return;
  const speed = r.value.axes[0];
  const alpha = r.value.axes[1];
  // Boundary fixtures verified per bound.
  assert.equal(speed?.lo, 8.0);
  assert.equal(speed?.loLimitedBy, "aero");
  assert.equal(speed?.hi, 16.0);
  assert.equal(speed?.hiLimitedBy, "structure");
  assert.equal(alpha?.lo, -0.05);
  assert.equal(alpha?.loLimitedBy, "aero");
  assert.equal(alpha?.hi, 0.12);
  assert.equal(alpha?.hiLimitedBy, "pilot");
  jlog("limiting", `"speed":"aero/structure","alpha":"aero/pilot"`);
});

test("empty intersections NAME the conflicting pair", () => {
  const conflicting: ApplicabilityDomain[] = [
    { subsystem: "a", axes: [{ name: "x", lo: 5, hi: 10 }] },
    { subsystem: "b", axes: [{ name: "x", lo: 0, hi: 4 }] },
  ];
  const r = intersectDomains(conflicting);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.equal(r.value.kind, "empty");
  if (r.value.kind !== "empty") return;
  assert.equal(r.value.axis, "x");
  assert.deepEqual(r.value.conflict, ["a", "b"]);
  // Standing in an empty domain refuses with the conflict named.
  const s = standingAt(r.value, { x: 4.5 });
  assert.ok(!s.ok && s.refusal.code === "applicability-empty-domain");
  assert.ok(!s.ok && s.refusal.message.includes("a") && s.refusal.message.includes("b"));
  jlog("empty", `"conflict":["a","b"]`);
});

test("AT a bound is inside; one past names the limiting subsystem", () => {
  const r = intersectDomains(domains());
  assert.ok(r.ok);
  if (!r.ok || r.value.kind !== "domain") return;
  const inter = r.value;
  // AT both speed bounds: inside (closed intervals).
  for (const speed of [8.0, 16.0]) {
    const s = standingAt(inter, { speed, alpha: 0.05 });
    assert.ok(s.ok && s.value.inside, `speed ${speed} is AT a bound = inside`);
  }
  // One ulp-ish past each bound: outside, correctly attributed.
  const low = standingAt(inter, { speed: 7.999_999, alpha: 0.05 });
  assert.ok(low.ok);
  if (low.ok && !low.value.inside) {
    assert.equal(low.value.axis, "speed");
    assert.equal(low.value.bound, "lo");
    assert.equal(low.value.limitedBy, "aero");
  } else {
    assert.fail("must be outside-lo");
  }
  const high = standingAt(inter, { speed: 16.000_001, alpha: 0.05 });
  assert.ok(high.ok);
  if (high.ok && !high.value.inside) {
    assert.equal(high.value.limitedBy, "structure");
  } else {
    assert.fail("must be outside-hi");
  }
  const alphaHigh = standingAt(inter, { speed: 12, alpha: 0.121 });
  assert.ok(alphaHigh.ok);
  if (alphaHigh.ok && !alphaHigh.value.inside) {
    assert.equal(alphaHigh.value.limitedBy, "pilot");
  } else {
    assert.fail("alpha must be outside-hi");
  }
  // Missing coordinate refuses.
  const missing = standingAt(inter, { speed: 12 });
  assert.ok(!missing.ok && missing.refusal.code === "applicability-point-invalid");
  jlog("bounds", `"at_bound":"inside","past_bound":"attributed"`);
});

function badge(state: EvidenceBadge["state"], digest: string | null): EvidenceBadge {
  return {
    caseId: "V-08b1",
    state,
    receiptDigest: digest,
    comparisonClass: "formulation-band-0.15",
    verdict: state === "evidenced-pass" ? "pass" : state === "reported" ? "reported-only" : "none",
  };
}

test("badge composition: green is EARNED, amber names the limiter, gray is honest", () => {
  const subsystems = ["aero", "structure", "pilot"];
  // Green: evidenced AND inside.
  const green = composeBadge(badge("evidenced-pass", DIGEST), { inside: true }, subsystems);
  assert.equal(green.color, "green");
  assert.equal(green.receiptLink, `receipt:${DIGEST}`);
  assert.match(green.sentence, /within declared applicability/);
  // FALSIFIER: evidenced but OUTSIDE -> amber, sentence names the
  // limiting subsystem and bound (never green).
  const amber = composeBadge(
    badge("evidenced-pass", DIGEST),
    { inside: false, axis: "speed", bound: "hi", limitedBy: "structure" },
    subsystems,
  );
  assert.equal(amber.color, "amber");
  assert.match(amber.sentence, /limited by structure/);
  assert.match(amber.sentence, /speed hi bound/);
  // NO-DATA: honest gray empty state, never blank, no link.
  const gray = composeBadge(badge("no-data", null), { inside: true }, subsystems);
  assert.equal(gray.color, "gray");
  assert.match(gray.sentence, /no receipt yet .* nothing is claimed/);
  assert.equal(gray.receiptLink, null);
  // Reported: gray with the verdict verbatim, never a pass claim.
  const reported = composeBadge(badge("reported", DIGEST), { inside: true }, subsystems);
  assert.equal(reported.color, "gray");
  assert.match(reported.sentence, /reported, not a pass claim/);
  assert.deepEqual(reported.subsystems, subsystems, "breakdown carried");
  jlog("compose", `"green":"earned","amber":"attributed","gray":"honest"`);
});

test("caps at cap AND cap+1 and axis-set mismatch refuses", () => {
  const mkDomains = (n: number): ApplicabilityDomain[] =>
    Array.from({ length: n }, (_, i) => ({
      subsystem: `s${i}`,
      axes: [{ name: "x", lo: 0, hi: 1 }],
    }));
  assert.ok(intersectDomains(mkDomains(MAX_SUBSYSTEMS)).ok, "AT cap");
  const over = intersectDomains(mkDomains(MAX_SUBSYSTEMS + 1));
  assert.ok(!over.ok && over.refusal.code === "applicability-subsystems-invalid");
  const mkAxes = (n: number): ApplicabilityDomain => ({
    subsystem: "s",
    axes: Array.from({ length: n }, (_, i) => ({ name: `a${i}`, lo: 0, hi: 1 })),
  });
  assert.ok(intersectDomains([mkAxes(MAX_AXES)]).ok, "AT axis cap");
  const overA = intersectDomains([mkAxes(MAX_AXES + 1)]);
  assert.ok(!overA.ok && overA.refusal.code === "applicability-axes-invalid");
  const mismatched = intersectDomains([
    { subsystem: "a", axes: [{ name: "x", lo: 0, hi: 1 }] },
    { subsystem: "b", axes: [{ name: "y", lo: 0, hi: 1 }] },
  ]);
  assert.ok(!mismatched.ok && mismatched.refusal.code === "applicability-axes-mismatched");
  const malformed = intersectDomains([
    { subsystem: "a", axes: [{ name: "x", lo: 2, hi: 1 }] },
  ]);
  assert.ok(!malformed.ok && malformed.refusal.code === "applicability-axes-invalid");
  jlog("caps", `"subsystems":${MAX_SUBSYSTEMS},"axes":${MAX_AXES}`);
});
