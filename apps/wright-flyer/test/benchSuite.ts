// E0.6a bench suite: synthetic kernels shaped like the §7.2 budget rows.
// Run: node test/benchSuite.ts [--json out.json]
// Every row is a HOST measurement; unmeasurable-headless rows emit NO-DATA.

import { writeFileSync } from "node:fs";
import { cpus } from "node:os";
import { bench, type BenchResult } from "../src/bench/harness.ts";
import {
  SeqlockWriter,
  seqlockBytes,
  type SeqlockLayout,
} from "../src/transport/seqlock.ts";
import { TransferablePool } from "../src/transport/pool.ts";

const jlog = (obj: object): void =>
  console.log(JSON.stringify({ suite: "wf-bench", ...obj }));

const results: BenchResult[] = [];
const noData: { name: string; reason: string }[] = [];

// 1. Biot–Savart-class kernel: 2,000 targets × 60 neighbors (the Tier-B wake
// feedback shape from §7.2/Appendix B): rsqrt-dominated 3-D kernel.
{
  const n = 2000;
  const k = 60;
  const px = new Float64Array(n).map(() => Math.random() * 10);
  const py = new Float64Array(n).map(() => Math.random() * 10);
  const pz = new Float64Array(n).map(() => Math.random() * 10);
  const g = new Float64Array(n).map(() => Math.random());
  const out = new Float64Array(3 * n);
  results.push(
    bench("biot-savart-2000x60", n * k, () => {
      for (let i = 0; i < n; i += 1) {
        let ux = 0;
        let uy = 0;
        let uz = 0;
        for (let j = 0; j < k; j += 1) {
          const s = (i + j * 31) % n;
          const dx = px[s]! - px[i]!;
          const dy = py[s]! - py[i]!;
          const dz = pz[s]! - pz[i]!;
          const r2 = dx * dx + dy * dy + dz * dz + 1e-4;
          const inv = g[s]! / (r2 * Math.sqrt(r2));
          ux += dy * dz * inv;
          uy += dz * dx * inv;
          uz += dx * dy * inv;
        }
        out[3 * i] = ux;
        out[3 * i + 1] = uy;
        out[3 * i + 2] = uz;
      }
    }, 30),
  );
}

// 2. Dense 80×80 LU factor+solve (the coupled multisurface solve shape).
{
  const n = 80;
  const a0 = new Float64Array(n * n).map(() => Math.random() - 0.5);
  for (let i = 0; i < n; i += 1) {
    a0[i * n + i] = 10 + Math.random(); // diagonally dominant
  }
  const b0 = new Float64Array(n).map(() => Math.random());
  results.push(
    bench("dense-lu-80", 1, () => {
      const a = Float64Array.from(a0);
      const b = Float64Array.from(b0);
      for (let kk = 0; kk < n; kk += 1) {
        const piv = a[kk * n + kk]!;
        for (let i = kk + 1; i < n; i += 1) {
          const f = a[i * n + kk]! / piv;
          a[i * n + kk] = f;
          for (let j = kk + 1; j < n; j += 1) {
            a[i * n + j] = a[i * n + j]! - f * a[kk * n + j]!;
          }
          b[i] = b[i]! - f * b[kk]!;
        }
      }
      for (let i = n - 1; i >= 0; i -= 1) {
        let s = b[i]!;
        for (let j = i + 1; j < n; j += 1) {
          s -= a[i * n + j]! * b[j]!;
        }
        b[i] = s / a[i * n + i]!;
      }
    }, 60),
  );
}

// 3. Transcendental batch (atmosphere modal sums; JS Math as the stand-in —
// fs-math det:: wasm rows join in E0.6c).
{
  const x = new Float64Array(4096).map((_, i) => (i % 631) / 100);
  const out = new Float64Array(4096);
  results.push(
    bench("transcendental-4096", 4096, () => {
      for (let i = 0; i < x.length; i += 1) {
        out[i] = Math.sin(x[i]!) * Math.exp(-0.1 * x[i]!) + Math.cos(0.5 * x[i]!);
      }
    }, 60),
  );
}

// 4. Seqlock publication (256-f64 state snapshot, the E0.7 ring).
{
  const layout: SeqlockLayout = { slots: 3, payloadF64s: 256 };
  const sab = new SharedArrayBuffer(seqlockBytes(layout));
  const writer = new SeqlockWriter(sab, layout, 1, 1, 1);
  let tick = 0;
  results.push(
    bench("seqlock-publish-256f64", 1, () => {
      tick += 1;
      writer.publish(tick, (p) => {
        for (let i = 0; i < p.length; i += 1) {
          p[i] = tick + i;
        }
      });
    }, 200),
  );
}

// 5. Transferable pack/copy (512-f64 field snapshot fallback path).
{
  const pool = new TransferablePool(64, 8 * 512);
  const src = new Float64Array(512).map(() => Math.random());
  results.push(
    bench("pool-pack-512f64", 1, () => {
      const b = pool.pack(src);
      if (b) {
        pool.acknowledge(b);
      }
    }, 200),
  );
}

// 6. f64→f32 downconvert (32³ field export shape).
{
  const src = new Float64Array(32 * 32 * 32 * 3).map(() => Math.random());
  const dst = new Float32Array(src.length);
  results.push(
    bench("f32-downconvert-98k", src.length, () => {
      dst.set(src as unknown as ArrayLike<number>);
    }, 60),
  );
}

// Unmeasurable headless: explicit NO-DATA (never silent omission).
noData.push(
  { name: "float32-gpu-upload", reason: "requires a real browser/GPU context (E0.6c)" },
  { name: "wasm-aero-kernels", reason: "real kernels land with E4.2/E4.5/E4.7 (E0.6c)" },
  { name: "browser-device-matrix", reason: "E0.6b via the E6.4 automation harness" },
);

for (const r of results) {
  jlog({ row: r });
}
for (const nd of noData) {
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
  rows: results,
  no_data: noData,
};
const outIdx = process.argv.indexOf("--json");
if (outIdx >= 0 && process.argv[outIdx + 1]) {
  writeFileSync(process.argv[outIdx + 1]!, JSON.stringify(artifact, null, 1));
  jlog({ event: "artifact-written", path: process.argv[outIdx + 1] });
}
