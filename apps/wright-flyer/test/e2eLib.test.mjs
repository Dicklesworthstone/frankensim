// Unit battery for e2e/lib.mjs (bead frankensim-xsz8b).
// Run: node --test test/e2eLib.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import {
  JsonlCapture,
  LOG_CAP,
  RefusalCodes,
  redactRecord,
  extractReceipts,
  extractQosStates,
  countLatencySamples,
  countActuatedSamples,
  compareRuns,
  digestLooksReal,
} from "../e2e/lib.mjs";

function line(obj) {
  return JSON.stringify(obj);
}

const HEX64 = "ab".repeat(32);

test("capture keeps only brace-prefixed lines and ignores chatter", () => {
  const cap = new JsonlCapture();
  assert.equal(cap.push("plain console text").ignored, true);
  const out = cap.push(line({ stage: "capability-probe", sharedArrayBuffer: true }));
  assert.equal(out.ok, true);
  assert.deepEqual(cap.lines(), [{ stage: "capability-probe", sharedArrayBuffer: true }]);
});

test("capture enforces cap exactly at cap and refuses at cap+1 (typed)", () => {
  const cap = new JsonlCapture(3);
  for (let i = 0; i < 3; i += 1) {
    assert.equal(cap.push(line({ i })).ok, true);
  }
  // First push BEYOND the cap is refused and dropped.
  const refused = cap.push(line({ i: 3 }));
  assert.equal(refused.ok, false);
  assert.equal(refused.code, RefusalCodes.LOG_CAP_EXCEEDED);
  // And so is every subsequent one.
  assert.equal(cap.push(line({ i: 4 })).code, RefusalCodes.LOG_CAP_EXCEEDED);
  assert.equal(cap.lines().length, 3);
});

test("capture records malformed JSONL as typed refusal, not a throw", () => {
  const cap = new JsonlCapture();
  const out = cap.push("{not json");
  assert.equal(out.ok, false);
  assert.equal(out.code, RefusalCodes.MALFORMED_JSONL);
  assert.equal(cap.refusals().length, 1);
});

test("redaction strips secret-looking keys at any depth", () => {
  const redacted = redactRecord({
    stage: "x",
    token: "t",
    nested: { apiKey2: "k", fine: 1 },
    listNote: "authorization",
  });
  assert.equal(redacted.token, "[REDACTED]");
  assert.equal(redacted.nested.apiKey2, "[REDACTED]");
  assert.equal(redacted.listNote, "authorization"); // value text is not a key
});

const GOOD_READY = { suite: "wright-flyer-app", stage: "sim-ready", tick0Digest: HEX64, runIntentId: "r1" };
const GOOD_TERMINAL = { suite: "wright-flyer-app", stage: "sim-terminal", phase: "ended:envelope-exceeded", tick: 1390, digest: "cd".repeat(32) };

test("extractReceipts accepts an ordered healthy boot", () => {
  const receipts = extractReceipts([
    line({ stage: "capability-probe", crossOriginIsolated: true }),
    line(GOOD_READY),
    line(GOOD_TERMINAL),
  ].map((l) => JSON.parse(l)));
  assert.equal(receipts.ok, true);
  assert.equal(receipts.tick0Digest, HEX64);
  assert.equal(receipts.finalDigest, GOOD_TERMINAL.digest);
  assert.equal(receipts.refusals.length, 0);
});

test("extractReceipts refuses missing sim-ready with typed MISSING_STAGE + MISSING_DIGEST", () => {
  const receipts = extractReceipts([
    line({ stage: "capability-probe" }),
    line(GOOD_TERMINAL),
  ].map((l) => JSON.parse(l)));
  assert.equal(receipts.ok, false);
  assert.ok(receipts.errors.some((e) => e.code === RefusalCodes.MISSING_STAGE && e.stage === "sim-ready"));
  assert.ok(receipts.errors.some((e) => e.code === RefusalCodes.MISSING_DIGEST && e.field === "tick0Digest"));
});

test("extractQosStates returns ordered governor states only", () => {
  const states = extractQosStates([
    { stage: "capability-probe" },
    { stage: "qos", state: "normal", enterConstrainedMs: 22 },
    { stage: "qos", state: "constrained" },
    { stage: "sim-terminal", digest: "x" },
    { stage: "qos", state: "normal" },
  ]);
  assert.deepEqual(states, ["normal", "constrained", "normal"]);
});

test("countLatencySamples counts only samples with applied_tick", () => {
  const n = countLatencySamples([
    { suite: "wf-input-latency", seq: 1, applied_tick: 12 },
    { suite: "wf-input-latency", seq: 2, applied_tick: null },
    { suite: "other" },
    { suite: "wf-input-latency", seq: 3, applied_tick: 15 },
  ]);
  assert.equal(n, 2);
});

test("countActuatedSamples counts applied NON-neutral commands only", () => {
  const n = countActuatedSamples([
    { suite: "wf-input-latency", seq: 1, applied_tick: 12, lever_n: -220, warp_rad: 0 },
    { suite: "wf-input-latency", seq: 2, applied_tick: 13, lever_n: 0, warp_rad: 0 },
    { suite: "wf-input-latency", seq: 3, applied_tick: null, lever_n: -220, warp_rad: 0 },
    { suite: "wf-input-latency", seq: 4, applied_tick: 15, lever_n: 0, warp_rad: 8.5 },
  ]);
  assert.equal(n, 2, "neutral heartbeat and unapplied sample excluded");
});

test("extractReceipts refuses out-of-order stages", () => {
  const receipts = extractReceipts([GOOD_TERMINAL, GOOD_READY].map((r) => ({ ...r })));
  assert.equal(receipts.ok, false);
  assert.ok(receipts.errors.some((e) => e.code === RefusalCodes.STAGE_OUT_OF_ORDER));
});

test("extractReceipts surfaces sim-refusal lines as run refusals", () => {
  const receipts = extractReceipts([
    { stage: "capability-probe" },
    { ...GOOD_READY },
    { stage: "sim-refusal", at: "step", code: "E_TEST", message: "boom" },
    { ...GOOD_TERMINAL },
  ]);
  assert.equal(receipts.refusals.length, 1);
  assert.equal(receipts.refusals[0].code, "E_TEST");
});

test("compareRuns: equal digests => IDENTICAL; any digest mismatch => DIVERGENT with fields", () => {
  const a = { tick0Digest: "a".repeat(20), finalDigest: "b".repeat(20) };
  assert.deepEqual(compareRuns(a, { ...a }).verdict, "IDENTICAL");
  const divergent = compareRuns(a, { tick0Digest: "a".repeat(20), finalDigest: "f".repeat(20) });
  assert.equal(divergent.verdict, "DIVERGENT");
  assert.deepEqual(divergent.fields, ["finalDigest"]);
});

test("digestLooksReal rejects placeholders, zero-digests, and short strings", () => {
  assert.equal(digestLooksReal("0".repeat(64)), false);
  assert.equal(digestLooksReal("deadbeef"), false); // too short
  assert.equal(digestLooksReal("not a digest!!"), false);
  assert.equal(digestLooksReal(null), false);
  assert.equal(digestLooksReal("b3".repeat(32)), true);
  assert.equal(digestLooksReal("base64url_digest-99"), true);
});

test("default cap constant is the workspace-bounded size", () => {
  assert.ok(LOG_CAP >= 512 && LOG_CAP <= 100000);
});
