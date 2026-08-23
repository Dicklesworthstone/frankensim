// E6.4a pure harness helpers (bead frankensim-xsz8b, leaf of
// wf-root-guzez.7.4). No DOM, no browser: everything here is unit-testable
// via `node --test test/e2eLib.test.mjs`.
//
// Receipts are the ENGINE's own chained digests surfaced by the app's
// console JSONL (`sim-ready` carries tick0Digest; `sim-terminal` carries
// { phase, tick, digest }). Sampled KPIs (snapshot counts, path lengths)
// are frame-cadence artifacts and are deliberately NOT determinism
// evidence.

export const LOG_CAP = 2000;

export const RefusalCodes = Object.freeze({
  LOG_CAP_EXCEEDED: "LOG_CAP_EXCEEDED",
  MALFORMED_JSONL: "MALFORMED_JSONL",
  MISSING_STAGE: "MISSING_STAGE",
  STAGE_OUT_OF_ORDER: "STAGE_OUT_OF_ORDER",
  MISSING_DIGEST: "MISSING_DIGEST",
  BAD_DIGEST_SHAPE: "BAD_DIGEST_SHAPE",
  RUN_DIVERGENT: "RUN_DIVERGENT",
  FALSIFIER_FAILED: "FALSIFIER_FAILED",
});

const SECRET_KEY_RE = /token|secret|password|authorization|cookie|api[-_]?key/i;

/** Drop values whose key looks secret-bearing; keep everything else. */
export function redactRecord(record) {
  if (record === null || typeof record !== "object" || Array.isArray(record)) {
    return record;
  }
  const out = {};
  for (const [key, value] of Object.entries(record)) {
    out[key] = SECRET_KEY_RE.test(key) ? "[REDACTED]" : redactRecord(value);
  }
  return out;
}

/**
 * Bounded console-JSONL capture. Accepts raw console text, keeps only
 * parseable `{...}` objects, redacts sensitive keys, and enforces the cap:
 * pushes AT the cap are accepted; the first push BEYOND the cap records a
 * typed refusal and is dropped (the capture stays bounded, never throws).
 */
export class JsonlCapture {
  #lines = [];
  #refusals = [];
  #capExceeded = false;

  constructor(cap = LOG_CAP) {
    if (!Number.isInteger(cap) || cap < 1) {
      throw new RangeError(`cap must be a positive integer, got ${cap}`);
    }
    this.#cap = cap;
  }

  #cap;

  /** Raw console text in; typed outcome out. */
  push(rawText) {
    if (typeof rawText !== "string") {
      return { ok: false, code: RefusalCodes.MALFORMED_JSONL };
    }
    if (!rawText.startsWith("{")) {
      return { ok: true, ignored: true };
    }
    let record;
    try {
      record = JSON.parse(rawText);
    } catch {
      this.#refusals.push({ code: RefusalCodes.MALFORMED_JSONL, at: this.#lines.length });
      return { ok: false, code: RefusalCodes.MALFORMED_JSONL };
    }
    if (this.#capExceeded) {
      this.#refusals.push({ code: RefusalCodes.LOG_CAP_EXCEEDED });
      return { ok: false, code: RefusalCodes.LOG_CAP_EXCEEDED };
    }
    this.#lines.push(redactRecord(record));
    if (this.#lines.length >= this.#cap) {
      // Cap reached: the NEXT push is the refusal (cap AND cap+1 semantics).
      this.#capExceeded = true;
    }
    return { ok: true, ignored: false };
  }

  lines() {
    return [...this.#lines];
  }

  refusals() {
    return [...this.#refusals];
  }
}

function stageOf(record) {
  return typeof record?.stage === "string" ? record.stage : null;
}

export const EXPECTED_STAGE_ORDER = Object.freeze([
  "capability-probe",
  "sim-ready",
  "sim-terminal",
]);

/**
 * Pull the run receipts out of one boot's captured lines. Typed refusals on
 * missing stages/digests or wrong order; per-item checks, never totals.
 */
export function extractReceipts(lines) {
  const errors = [];
  const stagesSeen = lines.map(stageOf).filter((s) => s !== null);
  const ordered = EXPECTED_STAGE_ORDER.filter((want) => stagesSeen.includes(want));
  for (let i = 1; i < ordered.length; i += 1) {
    if (stagesSeen.indexOf(ordered[i - 1]) > stagesSeen.indexOf(ordered[i])) {
      errors.push({
        code: RefusalCodes.STAGE_OUT_OF_ORDER,
        expected: ordered,
        observed: stagesSeen,
      });
      break;
    }
  }
  for (const want of EXPECTED_STAGE_ORDER) {
    if (!stagesSeen.includes(want)) {
      errors.push({ code: RefusalCodes.MISSING_STAGE, stage: want });
    }
  }
  const ready = lines.find((r) => stageOf(r) === "sim-ready") ?? null;
  const terminal = lines.find((r) => stageOf(r) === "sim-terminal") ?? null;
  const capabilityProbe = lines.find((r) => stageOf(r) === "capability-probe") ?? null;
  const tick0Digest = ready?.tick0Digest;
  if (typeof tick0Digest !== "string" || tick0Digest.length === 0) {
    errors.push({ code: RefusalCodes.MISSING_DIGEST, field: "tick0Digest" });
  }
  const finalDigest = terminal?.digest;
  if (typeof finalDigest !== "string" || finalDigest.length === 0) {
    errors.push({ code: RefusalCodes.MISSING_DIGEST, field: "finalDigest" });
  }
  const refusals = lines
    .filter((r) => stageOf(r) === "sim-refusal")
    .map((r) => ({ stage: r.at ?? "unknown", code: r.code ?? "unknown", message: r.message }));
  return {
    ok: errors.length === 0,
    errors,
    capabilityProbe,
    tick0Digest: typeof tick0Digest === "string" ? tick0Digest : null,
    finalDigest: typeof finalDigest === "string" ? finalDigest : null,
    terminalPhase: terminal?.phase ?? null,
    terminalTick: typeof terminal?.tick === "number" ? terminal.tick : null,
    refusals,
  };
}

/** Engine digests are opaque strings; reject placeholders and zero-digests. */
export function digestLooksReal(digest) {
  return (
    typeof digest === "string" &&
    digest.length >= 16 &&
    /^[0-9a-zA-Z_-]+$/.test(digest) &&
    !/^0+$/.test(digest)
  );
}

/**
 * Compare two boots' receipts. Verdict IDENTICAL requires equal tick0 and
 * final engine digests; anything else is DIVERGENT with a typed field list.
 * The comparator never sees KPIs — digests only.
 */
export function compareRuns(a, b) {
  const fields = [];
  if (a.tick0Digest !== b.tick0Digest) fields.push("tick0Digest");
  if (a.finalDigest !== b.finalDigest) fields.push("finalDigest");
  return fields.length === 0
    ? { verdict: "IDENTICAL", fields }
    : {
        verdict: "DIVERGENT",
        fields,
        detail: {
          a: { tick0: a.tick0Digest, final: a.finalDigest },
          b: { tick0: b.tick0Digest, final: b.finalDigest },
        },
      };
}
