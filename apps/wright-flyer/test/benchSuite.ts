// E0.6a bench suite: synthetic kernels shaped like the §7.2 budget rows.
// Run: node test/benchSuite.ts [--json out.json]
// Every row is a HOST measurement; unmeasurable-headless rows emit NO-DATA.
//
// E0.6b note: the kernel bodies moved to src/bench/kernels.ts so the SAME
// suite runs inside a real browser (?bench=1 via e2e/bench.mjs). This CLI
import { runBenchSuite, KERNEL_NAMES } from "../src/bench/kernels.ts";
import type { NoDataRow } from "../src/bench/kernels.ts";
import { writeFileSync } from "node:fs";
import { cpus } from "node:os";
import type { BenchResult } from "../src/bench/harness.ts";

const jlog = (obj: object): void =>
  console.log(JSON.stringify({ suite: "wf-bench", ...obj }));

const { rows: results, noData } = runBenchSuite();

// Node-host artifact keeps its own environment NO-DATA row: the browser
// matrix is measured by the E0.6b driver, not here.
const nodeOnlyNoData: NoDataRow[] = [
  { name: "browser-device-matrix", reason: "measured separately by e2e/bench.mjs (E0.6b)" },
];
const allNoData: NoDataRow[] = [...noData, ...nodeOnlyNoData];

for (const r of results as BenchResult[]) {
  jlog({ row: r });
}
for (const nd of allNoData) {
  jlog({ no_data: nd });
}

const artifact = {
  schema: "org.frankensim.wright-flyer.perf-baseline.v1",
  bead: "frankensim-wf-root-guzez.1.6.1",
  host: {
    kind: "node-headless",
    node: process.version,
    cpu: cpus()[0]?.model ?? "unknown",
    cores: cpus().length,
    note: "SHARED loaded development host — indicative baseline only; acceptance rows come from qualified quiet devices (E0.6b)",
  },
  kernel_order: KERNEL_NAMES,
  rows: results,
  no_data: allNoData,
};
const outIdx = process.argv.indexOf("--json");
if (outIdx >= 0 && process.argv[outIdx + 1]) {
  writeFileSync(process.argv[outIdx + 1]!, JSON.stringify(artifact, null, 1));
  jlog({ event: "artifact-written", path: process.argv[outIdx + 1] });
}
