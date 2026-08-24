// E0.6b: browser performance matrix rows via the E6.4 automation harness
// (bead frankensim-wf-root-guzez.1.6.2).
//
//   node e2e/bench.mjs                       # full 2x2 variant matrix
//   node e2e/bench.mjs --out <path.json>     # artifact destination
//
// Variants: isolation {sab, fallback} x contention {uncontended,
// cpu-throttled-4x}. Each variant boots a FRESH headless Chrome on the
// matching transport and runs the shared §7.2 kernel suite in-page
// (?bench=1&repeat=3): rep 0 is COLD, reps 1..2 are WARM.
//
// Gates (typed):
//   BENCH_INCOMPLETE  a variant lacks its cold/warm rows or sane percentiles
//   KERNEL_MISSING    an expected kernel row is absent for the transport
//                     (seqlock only under SAB; pool-pack everywhere)
// Exit 0 writes the v2 artifact carrying a sha256 content hash over the
// canonical JSON of every other field. Exit 4 gate failure, 5 setup refusal.
//
// Honesty: this is ONE host (shared, loaded dev box) x ONE browser
// (system Chrome headless) — a first matrix ROW, not the supported device
// matrix. Contended = CDP CPU throttle as a declared proxy for background
// load, not a second device. Other devices stay NO-DATA (rows carry the
// standing no_data entries from the kernel module).

import path from "node:path";
import { createHash } from "node:crypto";
import { writeFile } from "node:fs/promises";
import { startVitePreview, startStaticDist } from "./serve.mjs";
import { bootOnce, resolveChromeBin, BootRefusal } from "./boot.mjs";

const argv = process.argv.slice(2);
const outIdx = argv.indexOf("--out");
const outPath =
  outIdx >= 0 && argv[outIdx + 1]
    ? argv[outIdx + 1]
    : path.join("data", "wright-flyer", "perf-baseline-browser-chrome-dev-host.json");

function fail(code, payload, exitCode = 5) {
  console.error(JSON.stringify({ suite: "wf-bench-matrix", verdict: "FAIL", code, ...payload }));
  console.error(
    `repro: cd apps/wright-flyer && node e2e/bench.mjs${outIdx >= 0 ? ` --out ${outPath}` : ""}`,
  );
  process.exit(exitCode);
}

if (!resolveChromeBin()) {
  fail("CHROME_NOT_FOUND", { hint: "set WF_CHROME_BIN" });
}

/** Per-transport expected kernel set: seqlock needs SAB (isolation only). */
const KERNEL_NAMES_EXPECTED = {
  sab: [
    "biot-savart-2000x60",
    "dense-lu-80",
    "transcendental-4096",
    "seqlock-publish-256f64",
    "pool-pack-512f64",
    "f32-downconvert-98k",
  ],
  fallback: [
    "biot-savart-2000x60",
    "dense-lu-80",
    "transcendental-4096",
    "pool-pack-512f64",
    "f32-downconvert-98k",
  ],
};

const REPEATS = 3;
const TIMEOUT_MS = 240000;

const VARIANTS = [
  { isolation: "sab", serveKind: "vite", contention: "uncontended", throttle: null },
  { isolation: "sab", serveKind: "vite", contention: "cpu-throttled-4x", throttle: 4 },
  { isolation: "fallback", serveKind: "static", contention: "uncontended", throttle: null },
  { isolation: "fallback", serveKind: "static", contention: "cpu-throttled-4x", throttle: 4 },
];

async function runVariant(variant) {
  const serve = variant.serveKind === "static" ? await startStaticDist() : await startVitePreview();
  try {
    const boot = await bootOnce({
      baseUrl: serve.baseUrl,
      query: `bench=1&repeat=${REPEATS}`,
      timeoutMs: TIMEOUT_MS,
      cpuThrottleRate: variant.throttle,
      bench: true,
    });
    const benchLines = boot.lines.filter((r) => r?.suite === "wf-bench");
    if (!benchLines.some((r) => r.event === "suite-complete")) {
      fail("BENCH_INCOMPLETE", { variant, detail: "suite-complete event absent" }, 4);
    }
    const rows = benchLines
      .filter((r) => typeof r.row?.name === "string" && r.row.p50_us !== undefined)
      .map((r) => ({
        ...r.row,
        temperature: r.temperature === "warm" ? "warm" : "cold",
        rep: r.rep,
        isolation: variant.isolation,
        contention: variant.contention,
      }));
    if (rows.length === 0) {
      fail("BENCH_INCOMPLETE", { variant, detail: "no wf-bench rows captured" }, 4);
    }
    const expected = new Set(KERNEL_NAMES_EXPECTED[variant.isolation]);
    const got = new Set(rows.map((r) => r.name));
    for (const name of expected) {
      if (!got.has(name)) {
        fail("KERNEL_MISSING", { variant, kernel: name }, 4);
      }
      const per = rows.filter((r) => r.name === name);
      const colds = per.filter((r) => r.temperature === "cold").length;
      const warms = per.filter((r) => r.temperature === "warm").length;
      if (colds !== 1 || warms !== REPEATS - 1) {
        fail("BENCH_INCOMPLETE", { variant, kernel: name, colds, warms, warmExpected: REPEATS - 1 }, 4);
      }
      for (const r of per) {
        if (!(r.p50_us > 0 && r.p50_us <= r.p95_us && r.p95_us <= r.p99_us && r.opsPerSec > 0)) {
          fail("BENCH_INCOMPLETE", { variant, kernel: name, row: r, detail: "non-monotone percentiles" }, 4);
        }
      }
    }
    const unexpected = [...got].filter((n) => !expected.has(n));
    if (unexpected.length > 0) {
      fail("KERNEL_MISSING", { variant, detail: "unexpected kernels present", unexpected }, 4);
    }
    const noData = benchLines
      .filter((r) => r.no_data?.name)
      .map((r) => ({ ...r.no_data, isolation: variant.isolation, contention: variant.contention }));
    console.error(
      JSON.stringify({
        suite: "wf-bench-matrix",
        event: "variant-done",
        variant,
        kernels: [...got],
        noDataNames: [...new Set(noData.map((n) => n.name))],
      }),
    );
    return { rows, noData };
  } finally {
    serve.stop();
  }
}

const allRows = [];
const allNoData = [];
for (const variant of VARIANTS) {
  try {
    const { rows, noData } = await runVariant(variant);
    allRows.push(...rows);
    allNoData.push(...noData);
  } catch (error) {
    if (error instanceof BootRefusal) {
      fail(error.code, { variant: variant.isolation + "/" + variant.contention, message: error.message }, 4);
    }
    throw error;
  }
}

const artifact = {
  schema: "org.frankensim.wright-flyer.perf-baseline.v2",
  bead: "frankensim-wf-root-guzez.1.6.2",
  host: {
    kind: "chrome-headless-dev-host",
    browser: "system Chrome (headless, via puppeteer-core)",
    note: "SHARED loaded development host — indicative matrix ROW only; qualified quiet devices remain NO-DATA",
  },
  methodology: {
    repeats_per_variant: REPEATS,
    temperature: "rep 0 cold (fresh boot), reps 1..N warm (same page)",
    contention_proxy:
      "CDP Emulation.setCPUThrottlingRate(4x) — declared proxy for background load, not a second device",
    kernel_source: "src/bench/kernels.ts (identical to the E0.6a node suite)",
  },
  rows: allRows,
  no_data: allNoData,
};
const canonical = JSON.stringify(artifact, null, 1);
const contentHash = createHash("sha256").update(canonical).digest("hex");
const finalArtifact = { ...artifact, content_hash: contentHash };

await writeFile(outPath, JSON.stringify(finalArtifact, null, 1) + "\n");
console.log(
  JSON.stringify({
    suite: "wf-bench-matrix",
    verdict: "PASS",
    variants: VARIANTS.length,
    rows: allRows.length,
    artifact: outPath,
    content_hash: contentHash,
  }),
);
