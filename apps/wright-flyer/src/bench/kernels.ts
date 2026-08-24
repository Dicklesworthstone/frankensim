// §7.2-shaped microbench kernels shared by the node CLI
// (test/benchSuite.ts, bead wf-root-guzez.1.6.1 / E0.6a) and the
// in-browser runner (?bench=1 driven by e2e/bench.mjs, bead
// wf-root-guzez.1.6.2 / E0.6b). Browser-safe imports only: the SAME
// suite must measure identically-shaped work on every host, or
// per-device rows are not comparable.
//
// Every row is a MEASUREMENT OF A HOST — never an acceptance claim
// (plan §7.2.1). Unmeasurable contexts emit typed NO-DATA rows; nothing
// is silently omitted.

import { SeqlockWriter, seqlockBytes } from "../transport/seqlock.ts";
import { TransferablePool } from "../transport/pool.ts";
import { bench } from "./harness.ts";
import type { BenchResult } from "./harness.ts";
export interface NoDataRow {
  readonly name: string;
  readonly reason: string;
}

export interface SuiteOutput {
  readonly rows: readonly BenchResult[];
  readonly noData: readonly NoDataRow[];
}

/** Kernel order is FIXED so cross-host artifacts align row-for-row.
 * `seqlock-publish-256f64` requires cross-origin isolation (SAB); on a
 * fallback-transport origin it moves to noData instead. */
export const KERNEL_NAMES = [
  "biot-savart-2000x60",
  "dense-lu-80",
  "transcendental-4096",
  "seqlock-publish-256f64",
  "pool-pack-512f64",
  "f32-downconvert-98k",
] as const;

/** Standing NO-DATA rows emitted in unmeasured contexts (never silent). */
export function standingNoData(): NoDataRow[] {
  return [
    { name: "wasm-aero-kernels", reason: "real kernels land with E4.2/E4.5/E4.7 (E0.6c)" },
    { name: "simd-kernel-split", reason: "scalar/SIMD split arrives with the wasm kernels (E0.6c)" },
  ];
}

/**
 * Run the full kernel battery. Identical methodology to E0.6a's node
 * suite (same batch sizes, sample counts, warmup) so browser rows are
 * comparable with the committed node baseline.
 */
export function runBenchSuite(): SuiteOutput {
  const rows: BenchResult[] = [];
  const noData: NoDataRow[] = [];

  // 1. Biot–Savart-class kernel: 2,000 targets × 60 neighbors (the Tier-B
  // wake feedback shape): rsqrt-dominated 3-D kernel.
  {
    const n = 2000;
    const k = 60;
    const px = new Float64Array(n).map(() => Math.random() * 10);
    const py = new Float64Array(n).map(() => Math.random() * 10);
    const pz = new Float64Array(n).map(() => Math.random() * 10);
    const g = new Float64Array(n).map(() => Math.random());
    const out = new Float64Array(3 * n);
    rows.push(
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

  // 2. Dense 80×80 LU factor+solve (coupled multisurface solve shape).
  {
    const n = 80;
    const a0 = new Float64Array(n * n).map(() => Math.random() - 0.5);
    for (let i = 0; i < n; i += 1) {
      a0[i * n + i] = 10 + Math.random(); // diagonally dominant
    }
    const b0 = new Float64Array(n).map(() => Math.random());
    rows.push(
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

  // 3. Transcendental batch (atmosphere modal sums; JS Math stand-in —
  // fs-math det:: wasm rows join in E0.6c).
  {
    const x = new Float64Array(4096).map((_, i) => (i % 631) / 100);
    const out = new Float64Array(4096);
    rows.push(
      bench("transcendental-4096", 4096, () => {
        for (let i = 0; i < x.length; i += 1) {
          out[i] = Math.sin(x[i]!) * Math.exp(-0.1 * x[i]!) + Math.cos(0.5 * x[i]!);
        }
      }, 60),
    );
  }

  // 4. Seqlock publication (256-f64 state snapshot, the E0.7 ring).
  // Requires SharedArrayBuffer: on a non-isolated origin this becomes a
  // typed NO-DATA row rather than a fake or crashed measurement.
  if (typeof SharedArrayBuffer !== "undefined") {
    const layout = { slots: 3, payloadF64s: 256 };
    const sab = new SharedArrayBuffer(seqlockBytes(layout));
    const writer = new SeqlockWriter(sab, layout, 1, 1, 1);
    let tick = 0;
    rows.push(
      bench("seqlock-publish-256f64", 1, () => {
        tick += 1;
        writer.publish(tick, (p) => {
          for (let i = 0; i < p.length; i += 1) {
            p[i] = tick + i;
          }
        });
      }, 200),
    );
  } else {
    noData.push({
      name: "seqlock-publish-256f64",
      reason: "requires cross-origin isolation (SAB unavailable on this origin)",
    });
  }

  // 5. Transferable pack/copy (512-f64 field snapshot fallback path).
  {
    const pool = new TransferablePool(64, 8 * 512);
    const src = new Float64Array(512).map(() => Math.random());
    rows.push(
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
    rows.push(
      bench("f32-downconvert-98k", src.length, () => {
        dst.set(src as unknown as ArrayLike<number>);
      }, 60),
    );
  }

  // 7. Float32 GPU upload timing (WebGL2 buffer upload or typed fallback).
  if (typeof WebGL2RenderingContext !== "undefined" && typeof document !== "undefined") {
    try {
      const canvas = document.createElement("canvas");
      const gl = canvas.getContext("webgl2");
      if (gl) {
        const buffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
        const data = new Float32Array(32 * 32 * 32 * 3);
        gl.bufferData(gl.ARRAY_BUFFER, data.byteLength, gl.DYNAMIC_DRAW);
        rows.push(
          bench("float32-gpu-upload", data.length, () => {
            gl.bufferSubData(gl.ARRAY_BUFFER, 0, data);
          }, 60),
        );
      } else {
        noData.push({
          name: "float32-gpu-upload",
          reason: "WebGL2 context unavailable in this environment",
        });
      }
    } catch {
      noData.push({
        name: "float32-gpu-upload",
        reason: "WebGL2 initialization threw an exception",
      });
    }
  } else {
    noData.push({
      name: "float32-gpu-upload",
      reason: "requires a browser GPU/DOM context (WebGL2/WebGPU)",
    });
  }

  noData.push(...standingNoData());
  return { rows, noData };
}
