// E8.2-ii battery (bead wf-root-guzez.9.2.2): family grouping is
// engine-label passthrough (per-pole oracle, fixed family order);
// four-state teaching projection with the REPORTED residual and the
// actuator-mode falsifier (a mode dressed as rigid is caught by the
// caption law); polar redraw verbatim with the no-interpolation
// grid law; design-diff cards refuse unnamed attribution; caps at
// cap AND cap+1 throughout.
// Repro: node --test test/eigenView.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  FOUR_STATE_LABELS,
  MAX_CARD_METRICS,
  MAX_POLAR_POINTS,
  MAX_POLES,
  RIGID_SHARE_FLOOR,
  groupModeFamilies,
  polarDiff,
  teachingProjection,
  validateCard,
  validatePolar,
  type PolarPoint,
  type PublishedLabeledPole,
} from "../src/eigenView.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-eigenview","case":"${kase}",${payload}}`);
}

const ENGINE_LABELS = ["u", "w", "q", "theta", "dc", "dc_rate", "omega_rotor"] as const;

function pole(re: number, im: number, family: PublishedLabeledPole["family"]): PublishedLabeledPole {
  return { re, im, family, attributionShift: 0.1 };
}

test("family grouping: engine labels verbatim, fixed order, per-pole membership", () => {
  const poles = [
    pole(-0.5, 3.1, "rigid"),
    pole(0.35, 3.1, "rigid"),
    pole(-8, 12, "actuator"),
    pole(-0.9, 0, "rotor"),
  ];
  const r = groupModeFamilies(poles);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.deepEqual(
    r.value.map((g) => g.family),
    ["rigid", "actuator", "rotor"],
    "fixed family order, empty families dropped",
  );
  assert.equal(r.value[0]?.poles.length, 2);
  assert.equal(r.value[0]?.poles[1]?.re, 0.35, "order within family preserved");
  // Caps.
  const many = Array.from({ length: MAX_POLES }, () => pole(-1, 0, "rigid"));
  assert.ok(groupModeFamilies(many).ok, "AT cap");
  const over = groupModeFamilies([...many, pole(-1, 0, "rigid")]);
  assert.ok(!over.ok && over.refusal.code === "pole-count-invalid");
  const bad = groupModeFamilies([pole(Number.NaN, 0, "rigid")]);
  assert.ok(!bad.ok && bad.refusal.code === "pole-invalid");
  jlog("families", `"order":"rigid,actuator,rotor"`);
});

test("teaching projection: rigid mode projects cleanly, residual reported", () => {
  // A pure rigid phugoid-like shape: all content in u/w/q/theta.
  const mags = Float64Array.from([0.8, 0.5, 0.2, 0.27, 0, 0, 0]);
  const r = teachingProjection(ENGINE_LABELS, mags);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.ok(Math.abs(r.value.rigidShare - 1) < 1e-15);
  assert.ok(r.value.beyondFourStateContent < 1e-15);
  assert.deepEqual(
    r.value.components.map((c) => c.label),
    [...FOUR_STATE_LABELS],
  );
  assert.match(r.value.caption, /four-state projection/);
  jlog("projection-rigid", `"share":${r.value.rigidShare}`);
});

test("FALSIFIER: an actuator-dominated mode is never dressed as a rigid mode", () => {
  // Content almost entirely in the canard actuator states.
  const mags = Float64Array.from([0.02, 0.01, 0.03, 0.02, 0.9, 0.43, 0]);
  const r = teachingProjection(ENGINE_LABELS, mags);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.ok(r.value.rigidShare < RIGID_SHARE_FLOOR, `share ${r.value.rigidShare}`);
  assert.match(r.value.caption, /mostly NOT a rigid mode/);
  // The residual is REPORTED and large — hiding it would teach the
  // wrong airplane.
  assert.ok(r.value.beyondFourStateContent > 0.99);
  // Shape refusals.
  const short = teachingProjection(["u"], Float64Array.from([1, 2]));
  assert.ok(!short.ok && short.refusal.code === "projection-shape-mismatched");
  const zero = teachingProjection(ENGINE_LABELS, new Float64Array(7));
  assert.ok(!zero.ok && zero.refusal.code === "projection-vector-degenerate");
  jlog("projection-falsifier", `"beyond":${r.value.beyondFourStateContent}`);
});

test("polar redraw + design diff: verbatim, no interpolation", () => {
  const mk = (bump: number): PolarPoint[] =>
    Array.from({ length: 12 }, (_, i) => ({
      alphaRad: -0.05 + i * 0.02,
      cl: 0.3 + i * 0.08 + bump,
      cd: 0.03 + i * 0.004,
    }));
  const before = mk(0);
  const after = mk(0.05);
  assert.ok(validatePolar(before).ok);
  const d = polarDiff(before, after);
  assert.ok(d.ok);
  if (d.ok) {
    for (const p of d.value) {
      assert.ok(Math.abs(p.dCl - 0.05) < 1e-15, "exact delta at matched alpha");
      assert.equal(p.dCd, 0);
    }
  }
  // Mismatched grids REFUSE (no interpolation ever).
  const shifted = mk(0).map((p) => ({ ...p, alphaRad: p.alphaRad + 1e-9 }));
  const r = polarDiff(before, shifted);
  assert.ok(!r.ok && r.refusal.code === "polar-grids-mismatched");
  const shorter = polarDiff(before, before.slice(0, 11));
  assert.ok(!shorter.ok && shorter.refusal.code === "polar-grids-mismatched");
  // Unordered refuses; caps at cap AND cap+1.
  const unordered = validatePolar([before[1]!, before[0]!]);
  assert.ok(!unordered.ok && unordered.refusal.code === "polar-unordered");
  const big = (n: number): PolarPoint[] =>
    Array.from({ length: n }, (_, i) => ({ alphaRad: i * 1e-3, cl: 0, cd: 0 }));
  assert.ok(validatePolar(big(MAX_POLAR_POINTS)).ok, "AT cap");
  const overP = validatePolar(big(MAX_POLAR_POINTS + 1));
  assert.ok(!overP.ok && overP.refusal.code === "polar-length-invalid");
  jlog("polar", `"max_points":${MAX_POLAR_POINTS}`);
});

test("design-diff cards: named attribution required, metrics verbatim, caps", () => {
  const card = {
    change: "canard area +10%",
    attribution: "static margin shift via canard volume coefficient",
    metrics: [
      { metric: "time-to-double [s]", before: 1.98, after: 3.4 },
      { metric: "trim dc [rad]", before: 0.021, after: 0.017 },
    ],
  };
  const r = validateCard(card);
  assert.ok(r.ok);
  if (r.ok) {
    assert.ok(Object.is(r.value.metrics[0]?.before, 1.98), "metrics verbatim");
  }
  const unnamed = validateCard({ ...card, attribution: "  " });
  assert.ok(!unnamed.ok && unnamed.refusal.code === "card-attribution-missing");
  const empty = validateCard({ ...card, metrics: [] });
  assert.ok(!empty.ok && empty.refusal.code === "card-metrics-invalid");
  const mk = (n: number) =>
    Array.from({ length: n }, (_, i) => ({ metric: `m${i}`, before: 0, after: 1 }));
  assert.ok(validateCard({ ...card, metrics: mk(MAX_CARD_METRICS) }).ok, "AT cap");
  const over = validateCard({ ...card, metrics: mk(MAX_CARD_METRICS + 1) });
  assert.ok(!over.ok && over.refusal.code === "card-metrics-invalid");
  const nonfinite = validateCard({
    ...card,
    metrics: [{ metric: "x", before: Number.NaN, after: 0 }],
  });
  assert.ok(!nonfinite.ok && nonfinite.refusal.code === "card-metrics-invalid");
  jlog("cards", `"max_metrics":${MAX_CARD_METRICS}`);
});
