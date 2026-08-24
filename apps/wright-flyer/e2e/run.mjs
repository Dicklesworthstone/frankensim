// E6.4 determinism suite entry (beads frankensim-xsz8b, y473e, nty3a).
//
//   node e2e/run.mjs                # BOTH isolation rows, double-boot each
//   node e2e/run.mjs --row=sab      # SAB row only (vite preview, COOP/COEP)
//   node e2e/run.mjs --row=fallback # headerless static dist row
//   node e2e/run.mjs --row=human    # scripted keyboard input row
//   node e2e/run.mjs --row=qos      # CPU-throttled QoS escalation row
//   node e2e/run.mjs --all          # every row
//   node e2e/run.mjs --falsifier    # default vs ?flight=2 MUST diverge
//   node e2e/run.mjs --negative     # each row gate MUST refuse its negative control
//
// Receipts are ENGINE digests only for the isolation rows (sim-ready
// tick0Digest, sim-terminal digest). Human-row efficacy is proven by
// input-latency samples carrying applied_tick; cross-boot digest identity
// is NOT asserted there because control admission is device-time-derived.
// Exit: 0 pass · 3 RUN_DIVERGENT · 4 FALSIFIER/NEGATIVE failed · 5 setup refusal.
import path from "node:path";
import { writeFile } from "node:fs/promises";
import { startVitePreview, startStaticDist, makeArtifactDir } from "./serve.mjs";
import { bootOnce, resolveChromeBin, BootRefusal } from "./boot.mjs";
import {
  countLatencySamples,
  countActuatedSamples,
  digestLooksReal,
  extractQosStates,
  extractReceipts,
  RefusalCodes,
} from "./lib.mjs";

const argv = process.argv.slice(2);
const falsifierMode = argv.includes("--falsifier");
const rowArg = argv.find((a) => a.startsWith("--row="))?.slice(6) ?? "both";
const negFlag = argv.find((a) => a.startsWith("--negative"));
// Bare --negative runs BOTH gate negatives; --negative=human|qos narrows it.
const negativeGates =
  negFlag === undefined ? null : (negFlag.split("=")[1] ?? "both") === "both" ? ["human", "qos"] : [negFlag.split("=")[1]];
if (
  !["sab", "fallback", "human", "qos", "both", "all"].includes(rowArg) ||
  (negFlag !== undefined && ![undefined, "human", "qos", "both"].includes(negFlag.split("=")[1]))
) {
  console.error(JSON.stringify({ suite: "wf-e2e", verdict: "REFUSED", code: "BAD_ROW", row: rowArg, negative: negFlag ?? null }));
  process.exit(5);
}

function fail(code, payload, exitCode = 3) {
  console.error(JSON.stringify({ suite: "wf-e2e", verdict: "FAIL", code, ...payload }));
  console.error(`repro: cd apps/wright-flyer && node e2e/run.mjs${falsifierMode ? " --falsifier" : ` --row=${rowArg}`}`);
  process.exit(exitCode);
}

if (!resolveChromeBin()) {
  console.error(
    JSON.stringify({ suite: "wf-e2e", verdict: "REFUSED", code: "CHROME_NOT_FOUND", hint: "set WF_CHROME_BIN" }),
  );
  process.exit(5);
}

async function runRow(kind, artifactDir) {
  const serve = kind === "fallback" ? await startStaticDist() : await startVitePreview();
  try {
    let queryA = "sim=1";
    let queryB = falsifierMode ? "sim=1&flight=2" : "sim=1";
    let bootOpts = {};
    if (kind === "human") {
      queryA = queryB = "sim=1&mode=human";
      // Hold a fixed key set through the whole flight: a saturated plateau
      // is the strongest wall-clock-reproducible DOM script available.
      bootOpts = { holdKeys: ["ArrowUp", "KeyW"], timeoutMs: 90000 };
    }
    if (kind === "qos") {
      bootOpts = { cpuThrottleRate: 10, timeoutMs: 180000 };
    }

    const a = await bootOnce({ baseUrl: serve.baseUrl, query: queryA, ...bootOpts });
    const receiptsA = extractReceipts(a.lines);
    const b = await bootOnce({ baseUrl: serve.baseUrl, query: queryB, ...bootOpts });
    const receiptsB = extractReceipts(b.lines);

    const boots = [
      ["A", a, receiptsA],
      ["B", b, receiptsB],
    ];
    for (const [name, boot, receipts] of boots) {
      if (!receipts.ok) fail("BOOT_INCOMPLETE", { row: kind, boot: name, errors: receipts.errors });
      for (const field of ["tick0Digest", "finalDigest"]) {
        if (!digestLooksReal(receipts[field])) {
          fail("BAD_DIGEST_SHAPE", { row: kind, boot: name, field, observed: receipts[field] });
        }
      }
      if (receipts.refusals.length > 0) {
        fail("SIM_REFUSAL", { row: kind, boot: name, refusals: receipts.refusals });
      }
      await writeFile(
        path.join(artifactDir, `${kind}-boot-${name}.jsonl`),
        boot.lines.map((l) => JSON.stringify(l)).join("\n") + "\n",
      );
    }

    // Per-row gates beyond shared receipt health.
    if (kind === "sab" && receiptsA.capabilityProbe?.crossOriginIsolated !== true) {
      fail("TRANSPORT_MISDECLARED", { row: kind, expected: "crossOriginIsolated=true" });
    }
    if (kind === "fallback" && receiptsA.capabilityProbe?.crossOriginIsolated === true) {
      fail("TRANSPORT_MISDECLARED", { row: kind, expected: "crossOriginIsolated=false (headers leaked?)" });
    }
    if (kind === "human") {
      // Input-path EFFICACY: injected keydowns must drive a NON-neutral
      // command into an engine tick. The pump heartbeats neutrals even
      // with no keys held, so admitted-but-neutral samples prove nothing.
      const samples = Math.min(countActuatedSamples(a.lines), countActuatedSamples(b.lines));
      if (samples < 1) {
        fail("INPUT_PATH_SILENT", {
          row: kind,
          detail: "no non-neutral control admitted at an engine tick — injected keys never actuated",
        });
      }
      return {
        row: kind,
        verdict: "PASS",
        latencySamples: samples,
        tick0Digest: receiptsA.tick0Digest,
        finalDigestA: receiptsA.finalDigest,
        finalDigestB: receiptsB.finalDigest,
        terminalPhase: receiptsA.terminalPhase,
        note: "digest identity not asserted under wall-clock inputs (by design)",
      };
    }
    if (kind === "qos") {
      const states = extractQosStates(a.lines);
      if (!states.includes("constrained")) {
        fail("QOS_NOT_OBSERVED", { row: kind, statesSeen: states });
      }
      return {
        row: kind,
        verdict: "PASS",
        qosStates: states,
        tick0Digest: receiptsA.tick0Digest,
        finalDigest: receiptsA.finalDigest,
        terminalPhase: receiptsA.terminalPhase,
      };
    }

    const comparison = compareRuns(receiptsA, receiptsB);
    if (!falsifierMode && comparison.verdict !== "IDENTICAL") {
      fail(RefusalCodes.RUN_DIVERGENT, { row: kind, comparison });
    }
    if (falsifierMode && comparison.verdict !== "DIVERGENT") {
      fail(RefusalCodes.FALSIFIER_FAILED, { row: kind, comparison }, 4);
    }
    return {
      row: kind,
      verdict: falsifierMode ? "FALSIFIER-OK" : "PASS",
      tick0Digest: receiptsA.tick0Digest,
      finalDigest: receiptsA.finalDigest,
      terminalPhase: receiptsA.terminalPhase,
      divergentFields: comparison.fields,
      degradedProbe: receiptsA.capabilityProbe?.degraded ?? null,
    };
  } finally {
    serve.stop();
  }
}

/**
 * Per-gate falsifier negatives (bead frankensim-nty3a): every row gate must
 * REFUSE its negative control, or the gate cannot fire and is decoration.
 *   human: NO injected keys -> zero admitted-control latency samples must
 *          hold (INPUT_PATH_SILENT condition). Any sample means the input
 *          path admitted controls nobody sent.
 *   qos:   NO CDP throttling -> "constrained" must stay absent
 *          (QOS_NOT_OBSERVED condition). Escalation on an unthrottled boot
 *          means the gate cannot tell throttle from environment.
 */
async function runNegativeGate(gate, artifactDir) {
  const serve = await startVitePreview();
  try {
    const query = gate === "human" ? "sim=1&mode=human" : "sim=1";
    const bootOpts = gate === "human" ? { timeoutMs: 90000 } : { timeoutMs: 180000 };
    const a = await bootOnce({ baseUrl: serve.baseUrl, query, ...bootOpts });
    await writeFile(
      path.join(artifactDir, `negative-${gate}.jsonl`),
      a.lines.map((l) => JSON.stringify(l)).join("\n") + "\n",
    );
    if (gate === "human") {
      const samples = countActuatedSamples(a.lines);
      if (samples !== 0) {
        fail("NEGATIVE_FAILED", { gate, detail: `${samples} actuated samples without injected keys` }, 4);
      }
      return { gate, verdict: "NEGATIVE-OK", refuses: "INPUT_PATH_SILENT" };
    }
    const states = extractQosStates(a.lines);
    if (states.includes("constrained")) {
      fail("NEGATIVE_FAILED", { gate, statesSeen: states, detail: "constrained observed with no CDP throttling" }, 4);
    }
    return { gate, verdict: "NEGATIVE-OK", refuses: "QOS_NOT_OBSERVED", qosStates: states };
  } finally {
    serve.stop();
  }
}

const artifactDir = makeArtifactDir();
try {
  await writeFile(path.join(artifactDir, "mode.json"), JSON.stringify({ falsifierMode, rowArg, negativeGates }) + "\n");
} catch {}

try {
  const rows = negativeGates ?? (falsifierMode
    ? ["sab"]
    : rowArg === "both"
      ? ["sab", "fallback"]
      : rowArg === "all"
        ? ["sab", "fallback", "human", "qos"]
        : [rowArg]);
  const summaries = [];
  for (const row of rows) {
    summaries.push(negativeGates === null ? await runRow(row, artifactDir) : await runNegativeGate(row, artifactDir));
  }
  console.log(
    JSON.stringify({
      suite: "wf-e2e",
      verdict: negativeGates !== null ? "NEGATIVE-OK" : falsifierMode ? "FALSIFIER-OK" : "PASS",
      rows: summaries,
      artifacts: artifactDir,
    }),
  );
} catch (error) {
  if (error instanceof BootRefusal) fail(error.code, { message: error.message });
  fail("UNEXPECTED", { message: String(error?.stack ?? error) });
}
