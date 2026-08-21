// E5.2a facade battery (bead wf-root-guzez.6.3.1): the pure seam
// between the wasm engine's JSON envelopes and the worker protocol.
// Shape tests use contract-exact synthetic envelopes PLUS hostile
// twins (missing field, non-finite, unknown phase, malformed refusal);
// the LIVE section drives the REAL nodejs wasm pkg when WF_PKG points
// at it (fs_flyer_wasm.js) and proves parser-vs-engine agreement and
// the payload round trip on actual output — loudly reported either way.
// Repro: WF_PKG=<pkg>/fs_flyer_wasm.js node --test test/simFacade.test.ts

import assert from "node:assert/strict";
import { createRequire } from "node:module";
import { test } from "node:test";

import {
  fillPayload,
  parseDigestEnvelope,
  parseInitEnvelope,
  parseStepEnvelope,
} from "../src/sim/engineFacade.ts";
import {
  MODE_FIXED,
  P_ASSIST,
  P_DC_RAD,
  P_GUST_W_MPS,
  P_H_M,
  P_OMEGA_RAD_S,
  P_PHASE,
  P_Q_RAD_S,
  P_THETA_RAD,
  P_U_MPS,
  P_W_MPS,
  P_WARP_RAD,
  P_X_M,
  PAYLOAD_F64S,
  PHASE_CODES,
  payloadLayoutHash,
} from "../src/sim/protocol.ts";

const jlog = (payload: Record<string, unknown>): void => {
  console.info(JSON.stringify({ suite: "wf-e52a-facade", ...payload }));
};

const STEP_OK =
  '{"ok":{"tick":42,"phase":"airborne","x_m":10.5,"h_m":3.25,"u_mps":13.5,"w_mps":0.5,' +
  '"q_rad_s":0.1,"theta_rad":0.08,"dc_rad":0.12,"warp_rad":0.01,"omega_prop_rad_s":50.5,' +
  '"gust_w_mps":0.02,"assist_active":false}}';

test("phase codes mirror the native payload codes exactly", () => {
  assert.deepEqual(PHASE_CODES, {
    "on-rail": 0,
    airborne: 1,
    "ended:ground-contact": 2,
    "ended:rail-end-without-lift": 3,
    "ended:max-ticks": 4,
    "ended:envelope-exceeded": 5,
  });
});

test("payload layout hash is pinned (layout identity, not crypto)", () => {
  const h = payloadLayoutHash();
  assert.equal(h, payloadLayoutHash(), "deterministic");
  assert.notEqual(h, payloadLayoutHash("wf-snapshot-v2:something-else"), "descriptor-sensitive");
  jlog({ case: "layout-hash", value: h });
});

test("step envelope parses per-field and fills the frozen payload order", () => {
  const step = parseStepEnvelope(STEP_OK);
  assert.equal(step.kind, "ok");
  if (step.kind !== "ok") return;
  assert.equal(step.tick, 42);
  assert.equal(step.phase, "airborne");
  assert.equal(step.ended, false);
  const out = new Float64Array(PAYLOAD_F64S);
  fillPayload(step, out);
  // Per-slot oracle — NEVER a totals/spread check (sum tests are blind
  // to permutation).
  assert.equal(out[P_X_M], 10.5);
  assert.equal(out[P_H_M], 3.25);
  assert.equal(out[P_U_MPS], 13.5);
  assert.equal(out[P_W_MPS], 0.5);
  assert.equal(out[P_Q_RAD_S], 0.1);
  assert.equal(out[P_THETA_RAD], 0.08);
  assert.equal(out[P_DC_RAD], 0.12);
  assert.equal(out[P_WARP_RAD], 0.01);
  assert.equal(out[P_OMEGA_RAD_S], 50.5);
  assert.equal(out[P_GUST_W_MPS], 0.02);
  assert.equal(out[P_ASSIST], 0);
  assert.equal(out[P_PHASE], 1);
});

test("terminal phases parse as ended with the envelope receipt carried", () => {
  const json = STEP_OK.replace('"phase":"airborne"', '"phase":"ended:envelope-exceeded"').replace(
    ',"assist_active":false',
    ',"assist_active":true,"envelope_refusal_code":"PropAirframeCouplingDidNotConverge"',
  );
  const step = parseStepEnvelope(json);
  assert.equal(step.kind, "ok");
  if (step.kind !== "ok") return;
  assert.equal(step.ended, true);
  assert.equal(step.envelopeRefusalCode, "PropAirframeCouplingDidNotConverge");
  const out = new Float64Array(PAYLOAD_F64S);
  fillPayload(step, out);
  assert.equal(out[P_PHASE], 5);
  assert.equal(out[P_ASSIST], 1);
});

test("hostile twins refuse fail-closed (malformed, never a guess)", () => {
  // Missing field.
  assert.equal(parseStepEnvelope(STEP_OK.replace('"h_m":3.25,', "")).kind, "malformed");
  // Non-finite number survives JSON as null → wrong type.
  assert.equal(parseStepEnvelope(STEP_OK.replace("3.25", "null")).kind, "malformed");
  // Unknown phase word.
  assert.equal(
    parseStepEnvelope(STEP_OK.replace('"airborne"', '"ended:new-thing"')).kind,
    "malformed",
  );
  // Not JSON at all.
  assert.equal(parseStepEnvelope("<html>proxy error</html>").kind, "malformed");
  // Refusal with a broken shape is malformed, NOT silently accepted.
  assert.equal(parseStepEnvelope('{"refusal":{"code":7}}').kind, "malformed");
  // Well-formed refusal parses as refusal.
  const refusal = parseStepEnvelope(
    '{"refusal":{"code":"run-ended","message":"m","ranked_repairs":["r"]}}',
  );
  assert.equal(refusal.kind, "refusal");
  if (refusal.kind === "refusal") {
    assert.equal(refusal.refusal.code, "run-ended");
  }
  // Init envelope twins.
  assert.equal(parseInitEnvelope('{"ok":{"run_intent_id":"x"}}').kind, "malformed");
  // Digest: 64-hex or refuse (63 chars = cap-1 twin).
  assert.equal(parseDigestEnvelope(`{"ok":{"digest":"${"a".repeat(64)}"}}`), "a".repeat(64));
  const short = parseDigestEnvelope(`{"ok":{"digest":"${"a".repeat(63)}"}}`);
  assert.equal(typeof short === "string" ? "string" : short.kind, "malformed");
});

// --------------------------------------------------------------------------
// LIVE section: the real wasm engine (nodejs target) when WF_PKG is set.
// --------------------------------------------------------------------------
const pkgPath = process.env.WF_PKG;
if (pkgPath === undefined) {
  // Loud, structured, and impossible to mistake for coverage.
  jlog({ case: "live-engine", skipped: true, reason: "WF_PKG unset — parser-vs-engine agreement NOT verified in this run" });
} else {
  const require = createRequire(import.meta.url);
  const wasm = require(pkgPath);

  test("LIVE: init envelope from the real engine parses with identity fields", () => {
    const init = parseInitEnvelope(
      wasm.flyer_engine_init(1903n, 1.294, 11.0, MODE_FIXED, 0, 18.3, 24n, false, false),
    );
    assert.equal(init.kind, "ok");
    if (init.kind !== "ok") return;
    assert.match(init.tick0Digest, /^[0-9a-f]{64}$/);
    assert.ok(init.trimVMps > 5.0);
    jlog({ case: "live-init", tick0: init.tick0Digest.slice(0, 16) });
  });

  test("LIVE: every step parses; payload round-trips; terminal in-band", () => {
    let last: ReturnType<typeof parseStepEnvelope> | null = null;
    const out = new Float64Array(PAYLOAD_F64S);
    for (let i = 0; i < 24; i += 1) {
      const step = parseStepEnvelope(wasm.flyer_engine_step(false, 0, 0));
      assert.equal(step.kind, "ok", JSON.stringify(step));
      if (step.kind !== "ok") return;
      fillPayload(step, out);
      assert.equal(out[P_PHASE], PHASE_CODES[step.phase]);
      assert.equal(out[P_U_MPS], step.state.uMps, "bitwise field carry");
      last = step;
    }
    assert.ok(last !== null && last.kind === "ok" && last.ended, "max-ticks terminal parsed");
    const digest = parseDigestEnvelope(wasm.flyer_engine_digest());
    assert.equal(typeof digest, "string", "real digest parses");
    jlog({ case: "live-lifecycle", digest: typeof digest === "string" ? digest.slice(0, 16) : "?" });
  });

  test("LIVE: real refusal envelopes parse as refusals with stable codes", () => {
    const past = parseStepEnvelope(wasm.flyer_engine_step(false, 0, 0));
    assert.equal(past.kind, "refusal");
    if (past.kind === "refusal") {
      assert.equal(past.refusal.code, "run-ended");
    }
    const bad = parseInitEnvelope(wasm.flyer_engine_init(1n, 1.294, 11.0, 3, 0, 18.3, 24n, false, false));
    assert.equal(bad.kind, "refusal");
    if (bad.kind === "refusal") {
      assert.equal(bad.refusal.code, "mode-invalid");
    }
  });
}
