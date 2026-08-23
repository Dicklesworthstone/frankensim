// E6.4 determinism suite entry (beads frankensim-xsz8b, frankensim-y473e).
//
//   node e2e/run.mjs                # BOTH isolation rows, double-boot each
//   node e2e/run.mjs --row=sab      # SAB row only (vite preview, COOP/COEP)
//   node e2e/run.mjs --row=fallback # headerless static dist row
//   node e2e/run.mjs --falsifier    # default vs ?flight=2 MUST diverge
//
// Receipts are ENGINE digests only (sim-ready tick0Digest, sim-terminal
// digest). Snapshot counts and KPIs are frame-cadence artifacts.
// Exit: 0 pass · 3 RUN_DIVERGENT · 4 FALSIFIER_FAILED · 5 setup refusal.
import path from "node:path";
import { writeFile } from "node:fs/promises";
import { startVitePreview, startStaticDist, makeArtifactDir } from "./serve.mjs";
import { bootOnce, resolveChromeBin, BootRefusal } from "./boot.mjs";
import { compareRuns, digestLooksReal, extractReceipts, RefusalCodes } from "./lib.mjs";

const argv = process.argv.slice(2);
const falsifierMode = argv.includes("--falsifier");
const rowArg = argv.find((a) => a.startsWith("--row="))?.slice(6) ?? "both";
if (!["sab", "fallback", "both"].includes(rowArg)) {
  console.error(JSON.stringify({ suite: "wf-e2e", verdict: "REFUSED", code: "BAD_ROW", row: rowArg }));
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

/**
 * One row = one server + two cold boots + comparator + per-row transport
 * assertions. Returns a summary object; throws typed refusals on failure.
 */
async function runRow(kind, artifactDir) {
  const serve =
    kind === "sab"
      ? await startVitePreview()
      : await startStaticDist();
  try {
    const queryA = "sim=1";
    const queryB = falsifierMode ? "sim=1&flight=2" : "sim=1";

    const a = await bootOnce({ baseUrl: serve.baseUrl, query: queryA });
    const receiptsA = extractReceipts(a.lines);
    const b = await bootOnce({ baseUrl: serve.baseUrl, query: queryB });
    const receiptsB = extractReceipts(b.lines);

    for (const [name, boot, receipts] of [["A", a, receiptsA], ["B", b, receiptsB]]) {
      if (!receipts.ok) fail("BOOT_INCOMPLETE", { row: kind, boot: name, errors: receipts.errors });
      for (const field of ["tick0Digest", "finalDigest"]) {
        if (!digestLooksReal(receipts[field])) {
          fail("BAD_DIGEST_SHAPE", { row: kind, boot: name, field, observed: receipts[field] });
        }
      }
      if (receipts.refusals.length > 0) {
        fail("SIM_REFUSAL", { row: kind, boot: name, refusals: receipts.refusals });
      }
      // Transport honesty per row: the probe line must declare the mode we
      // actually served (the app's own runtime declaration).
      const isolated = receipts.capabilityProbe?.crossOriginIsolated === true;
      if (kind === "sab" && !isolated) {
        fail("TRANSPORT_MISDECLARED", { row: kind, expected: "crossOriginIsolated=true" });
      }
      if (kind === "fallback" && isolated) {
        fail("TRANSPORT_MISDECLARED", { row: kind, expected: "crossOriginIsolated=false (headers leaked?)" });
      }
      await writeFile(
        path.join(artifactDir, `${kind}-boot-${name}.jsonl`),
        boot.lines.map((l) => JSON.stringify(l)).join("\n") + "\n",
      );
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

let artifactDir = makeArtifactDir();
try {
  await writeFile(path.join(artifactDir, "mode.json"), JSON.stringify({ falsifierMode, rowArg }) + "\n");
} catch {}

try {
  const rows = falsifierMode ? ["sab"] : rowArg === "both" ? ["sab", "fallback"] : [rowArg];
  const summaries = [];
  for (const row of rows) {
    summaries.push(await runRow(row, artifactDir));
  }
  console.log(
    JSON.stringify({
      suite: "wf-e2e",
      verdict: falsifierMode ? "FALSIFIER-OK" : "PASS",
      rows: summaries,
      artifacts: artifactDir,
    }),
  );
} catch (error) {
  if (error instanceof BootRefusal) fail(error.code, { message: error.message });
  fail("UNEXPECTED", { message: String(error?.stack ?? error) });
}
