// E6.4a determinism suite entry (bead frankensim-xsz8b).
//
//   node e2e/run.mjs              # two fresh boots, ?sim=1: digests MUST match
//   node e2e/run.mjs --falsifier  # default vs ?flight=2: digests MUST differ
//
// Exit codes: 0 pass · 3 RUN_DIVERGENT · 4 FALSIFIER_FAILED · 5 setup refusal.
import path from "node:path";
import { writeFile, mkdir } from "node:fs/promises";
import { startVitePreview, makeArtifactDir } from "./serve.mjs";
import { bootOnce, resolveChromeBin, BootRefusal } from "./boot.mjs";
import {
  compareRuns,
  digestLooksReal,
  extractReceipts,
  RefusalCodes,
} from "./lib.mjs";

const args = new Set(process.argv.slice(2));
const falsifierMode = args.has("--falsifier");

function fail(code, payload) {
  console.error(JSON.stringify({ suite: "wf-e2e", verdict: "FAIL", code, ...payload }));
  console.error(`repro: cd apps/wright-flyer && node e2e/run.mjs${falsifierMode ? " --falsifier" : ""}`);
  process.exit(code === RefusalCodes.FALSIFIER_FAILED ? 4 : 3);
}

if (!resolveChromeBin()) {
  console.error(
    JSON.stringify({
      suite: "wf-e2e",
      verdict: "REFUSED",
      code: "CHROME_NOT_FOUND",
      hint: "set WF_CHROME_BIN",
    }),
  );
  process.exit(5);
}

const serve = await startVitePreview();
let artifactDir = null;
try {
  const bootQueryA = "sim=1";
  const bootQueryB = falsifierMode ? "sim=1&flight=2" : "sim=1";

  const a = await bootOnce({ baseUrl: serve.baseUrl, query: bootQueryA });
  const receiptsA = extractReceipts(a.lines);
  const b = await bootOnce({ baseUrl: serve.baseUrl, query: bootQueryB });
  const receiptsB = extractReceipts(b.lines);

  for (const [name, receipts] of [["A", receiptsA], ["B", receiptsB]]) {
    if (!receipts.ok) {
      fail("BOOT_INCOMPLETE", { boot: name, errors: receipts.errors });
    }
    for (const field of ["tick0Digest", "finalDigest"]) {
      if (!digestLooksReal(receipts[field])) {
        fail("BAD_DIGEST_SHAPE", { boot: name, field, observed: receipts[field] });
      }
    }
    if (receipts.refusals.length > 0) {
      fail("SIM_REFUSAL", { boot: name, refusals: receipts.refusals });
    }
  }

  const comparison = compareRuns(receiptsA, receiptsB);
  artifactDir = makeArtifactDir();
  await mkdir(path.join(artifactDir), { recursive: true });
  await writeFile(path.join(artifactDir, "boot-A.jsonl"), a.lines.map((l) => JSON.stringify(l)).join("\n") + "\n");
  await writeFile(path.join(artifactDir, "boot-B.jsonl"), b.lines.map((l) => JSON.stringify(l)).join("\n") + "\n");
  await writeFile(
    path.join(artifactDir, "summary.json"),
    JSON.stringify({ mode: falsifierMode ? "falsifier" : "determinism", comparison, queries: [bootQueryA, bootQueryB] }, null, 2) + "\n",
  );

  if (!falsifierMode && comparison.verdict !== "IDENTICAL") {
    fail(RefusalCodes.RUN_DIVERGENT, { comparison, artifacts: artifactDir });
  }
  if (falsifierMode && comparison.verdict !== "DIVERGENT") {
    // The comparator saw two DIFFERENT scenarios as identical: the receipt
    // channel or comparator is broken. This mode exists to catch exactly that.
    fail(RefusalCodes.FALSIFIER_FAILED, { comparison, artifacts: artifactDir });
  }

  console.log(
    JSON.stringify({
      suite: "wf-e2e",
      verdict: falsifierMode ? "FALSIFIER-OK" : "PASS",
      mode: falsifierMode ? "falsifier" : "determinism",
      tick0Digest: receiptsA.tick0Digest,
      finalDigest: receiptsA.finalDigest,
      terminalPhase: receiptsA.terminalPhase,
      divergentFields: comparison.fields,
      artifacts: artifactDir,
    }),
  );
} catch (error) {
  if (error instanceof BootRefusal) {
    fail(error.code, { message: error.message });
  }
  fail("UNEXPECTED", { message: String(error?.stack ?? error) });
} finally {
  serve.stop();
}
