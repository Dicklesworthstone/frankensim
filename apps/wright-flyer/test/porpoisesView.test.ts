// E7.4 battery (bead wf-root-guzez.8.5): the flagship view's laws on
// a PINNED porpoising fixture. Time-to-double vs closed form (the
// one allowed derivation); delay indicator recovers a known shift
// exactly; attribution ranking with the mislabeling falsifier and
// the displayed residual; saturation/reversal events verbatim with
// order refusal; A/B common-prefix semantics with the inconsistent-
// receipt falsifier EXECUTED both ways.
// Repro: node --test test/porpoisesView.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_ATTRIBUTION_RESIDUAL,
  MAX_EVENTS,
  MAX_LAG_TICKS,
  abAnnotation,
  attributionView,
  estimateDelayTicks,
  poleIndicator,
  validateEvents,
} from "../src/porpoisesView.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-porpoises","case":"${kase}",${payload}}`);
}

// ---- The pinned porpoising incident fixture -----------------------
// The historical-mode PIO class: an unstable phugoid-like pair at
// sigma = 0.35 1/s, omega = 3.1 rad/s (period ~2.03 s, doubling
// ~1.98 s), pilot lag 9 ticks, canard saturating on alternate swings.
const FIXTURE_SIGMA = 0.35;
const FIXTURE_OMEGA = 3.1;
const FIXTURE_LAG_TICKS = 9;
const N = 480;

function fixtureCommand(): Float64Array {
  const c = new Float64Array(N);
  for (let i = 0; i < N; i += 1) {
    const t = i / 120;
    c[i] = Math.exp(FIXTURE_SIGMA * t) * Math.sin(FIXTURE_OMEGA * t) * 0.05;
  }
  return c;
}

function fixtureActual(): Float64Array {
  const cmd = fixtureCommand();
  const a = new Float64Array(N);
  for (let i = 0; i < N; i += 1) {
    a[i] = cmd[Math.max(0, i - FIXTURE_LAG_TICKS)] ?? 0;
  }
  return a;
}

test("pole indicator: time-to-double matches closed form on the fixture pole", () => {
  const r = poleIndicator({ reSigmaPerS: FIXTURE_SIGMA, imOmegaRadPerS: FIXTURE_OMEGA });
  assert.ok(r.ok);
  if (!r.ok) return;
  const expected = Math.LN2 / FIXTURE_SIGMA;
  assert.ok(Math.abs((r.value.timeToDoubleS ?? 0) - expected) < 1e-15);
  assert.ok(Math.abs((r.value.periodS ?? 0) - (2 * Math.PI) / FIXTURE_OMEGA) < 1e-15);
  // Stable pole: NO time-to-double (never a negative doubling time).
  const stable = poleIndicator({ reSigmaPerS: -0.2, imOmegaRadPerS: 1 });
  assert.ok(stable.ok && stable.value.timeToDoubleS === null);
  // Neutral: also none.
  const neutral = poleIndicator({ reSigmaPerS: 0, imOmegaRadPerS: 1 });
  assert.ok(neutral.ok && neutral.value.timeToDoubleS === null);
  const bad = poleIndicator({ reSigmaPerS: Number.NaN, imOmegaRadPerS: 0 });
  assert.ok(!bad.ok && bad.refusal.code === "pole-invalid");
  jlog("pole", `"t2":${expected}`);
});

test("delay indicator recovers the fixture lag exactly", () => {
  const r = estimateDelayTicks(fixtureCommand(), fixtureActual());
  assert.ok(r.ok);
  // Positive lag: actual[i + lag] aligns with command[i], i.e. the
  // actual trails the command by 9 ticks.
  assert.equal(r.ok && r.value, FIXTURE_LAG_TICKS, "actual trails command by 9 ticks");
  // Window refusal.
  const short = estimateDelayTicks(new Float64Array(10), new Float64Array(10));
  assert.ok(!short.ok && short.refusal.code === "delay-window-invalid");
  jlog("delay", `"lag_ticks":${FIXTURE_LAG_TICKS},"max_lag":${MAX_LAG_TICKS}`);
});

test("attribution: ranked passthrough, residual displayed, mislabel falsifier", () => {
  const shares = [
    { component: "pilot-delay", share: 0.46 },
    { component: "canard-authority", share: 0.31 },
    { component: "airframe-pole", share: 0.22 },
  ];
  const r = attributionView(shares);
  assert.ok(r.ok);
  if (!r.ok) return;
  assert.equal(r.value.dominant, "pilot-delay");
  assert.deepEqual(
    r.value.ranked.map((s) => s.component),
    ["pilot-delay", "canard-authority", "airframe-pole"],
  );
  assert.ok(Math.abs(r.value.residual - 0.01) < 1e-12, "residual displayed, not hidden");
  // FALSIFIER: swapped shares flip the dominant label — the per-item
  // ranking oracle catches the mislabeling.
  const swapped = attributionView([
    { component: "pilot-delay", share: 0.22 },
    { component: "canard-authority", share: 0.31 },
    { component: "airframe-pole", share: 0.46 },
  ]);
  assert.ok(swapped.ok && swapped.value.dominant === "airframe-pole");
  assert.notEqual(swapped.ok && swapped.value.dominant, r.value.dominant);
  // Broken split refuses (residual beyond the bound).
  const broken = attributionView([{ component: "x", share: 0.5 }]);
  assert.ok(!broken.ok && broken.refusal.code === "attribution-residual-exceeded");
  const empty = attributionView([]);
  assert.ok(!empty.ok && empty.refusal.code === "attribution-empty");
  jlog("attribution", `"residual_bound":${MAX_ATTRIBUTION_RESIDUAL}`);
});

test("events verbatim: order refusal, cap AND cap+1", () => {
  const ok = validateEvents([
    { tick: 10, kind: "saturation-enter" },
    { tick: 22, kind: "saturation-exit" },
    { tick: 22, kind: "command-reversal" },
  ]);
  assert.ok(ok.ok, "equal ticks are fine (same-tick events)");
  const unordered = validateEvents([
    { tick: 22, kind: "saturation-enter" },
    { tick: 10, kind: "saturation-exit" },
  ]);
  assert.ok(!unordered.ok && unordered.refusal.code === "events-unordered");
  const mk = (n: number) =>
    Array.from({ length: n }, (_, i) => ({ tick: i, kind: "command-reversal" as const }));
  assert.ok(validateEvents(mk(MAX_EVENTS)).ok, "AT cap");
  const over = validateEvents(mk(MAX_EVENTS + 1));
  assert.ok(!over.ok && over.refusal.code === "event-count-exceeded");
  jlog("events", `"cap":${MAX_EVENTS}`);
});

test("A/B common-prefix semantics + inconsistent-receipt falsifier", () => {
  // Two traces identical for 200 ticks, then diverging (atmosphere
  // realization held fixed — the divergence is the counterfactual).
  const a = fixtureCommand();
  const b = Float64Array.from(a);
  for (let i = 200; i < N; i += 1) b[i] = (b[i] ?? 0) * 1.02;
  const good = abAnnotation({ commonPrefixTicks: 200, divergent: true }, a, b);
  assert.ok(good.ok);
  if (good.ok) {
    assert.equal(good.value.sharedUntilTick, 200);
    assert.equal(good.value.divergenceFromTick, 200);
  }
  // FALSIFIER 1: a receipt claiming a LONGER common prefix than the
  // traces support is refused (divergence annotations may never be
  // hidden behind a stale receipt).
  const stale = abAnnotation({ commonPrefixTicks: 300, divergent: true }, a, b);
  assert.ok(!stale.ok && stale.refusal.code === "ab-receipt-inconsistent");
  // FALSIFIER 2: a receipt claiming identity over divergent traces.
  const denial = abAnnotation({ commonPrefixTicks: N, divergent: false }, a, b);
  assert.ok(!denial.ok && denial.refusal.code === "ab-receipt-inconsistent");
  // Identical traces + identical receipt: no divergence annotation.
  const same = abAnnotation({ commonPrefixTicks: N, divergent: false }, a, Float64Array.from(a));
  assert.ok(same.ok && same.ok && same.value.divergenceFromTick === null);
  // A receipt may honestly claim a SHORTER prefix than the last
  // identical sample (conservative receipts allowed): annotations
  // simply begin earlier — never refused.
  const conservative = abAnnotation({ commonPrefixTicks: 150, divergent: true }, a, b);
  assert.ok(conservative.ok && conservative.value.sharedUntilTick === 150);
  jlog("ab", `"falsifiers":"executed-both"`);
});
