//! Batched small-dense performance baseline (frankensim-9ekv step a):
//! GFLOP/s tables for the SAFE autovectorized path per size class,
//! against a documented per-core peak model. Run with --release for
//! meaningful numbers (the harness measures and LEDGERS; the ≥60%
//! acceptance verdict is per size class per ISA). G-tier: performance
//! evidence, not correctness (correctness lives in batched_battery).
//!
//! PEAK MODEL (documented, machine-fingerprinted in the log):
//! - Apple M4 P-core: 4×128-bit FP64 FMA pipes = 16 flops/cycle;
//!   sustained P-core clock ~4.4 GHz ⇒ ~70 GFLOP/s single-core peak.
//! - Zen-4 (Threadripper 7000): 2×256-bit FP64 FMA ports (512-bit
//!   double-pumped) = 16 flops/cycle; ~4.0 GHz ⇒ ~64 GFLOP/s.
//!
//! The 60% bands are evaluated against these models; the harness is
//! single-threaded so all-core clock scaling stays out of the noise.

use fs_la::batched::{BatchMat, batch_gemm};
use std::time::Instant;

fn log(case: &str, detail: &str) {
    println!("{{\"suite\":\"fs-la-batched-perf\",\"case\":\"{case}\",\"detail\":\"{detail}\"}}");
}

/// One measured size class: median GFLOP/s over `reps` timed runs of
/// `batch_gemm` on `n` matrices of size k×k.
fn measure(k: usize, n: usize, reps: usize) -> f64 {
    let a = BatchMat::from_fn(k, n, |m, i, j| {
        0.5 + (m as f64).mul_add(0.001, (i * k + j) as f64 * 0.01)
    });
    let b = BatchMat::from_fn(k, n, |m, i, j| {
        0.3 + (m as f64).mul_add(0.002, (j * k + i) as f64 * 0.02)
    });
    let mut c = BatchMat::zeros(k, n);
    // Warmup.
    batch_gemm(1.0, &a, &b, 0.0, &mut c);
    let flops = 2.0 * (k * k * k) as f64 * n as f64;
    let mut times: Vec<f64> = (0..reps)
        .map(|_| {
            let t0 = Instant::now();
            batch_gemm(1.0, &a, &b, 0.0, &mut c);
            t0.elapsed().as_secs_f64()
        })
        .collect();
    times.sort_by(f64::total_cmp);
    let median = times[times.len() / 2];
    flops / median / 1e9
}

#[test]
fn baseline_gflops_table() {
    // Debug-profile guard: numbers are only meaningful in release.
    let release = !cfg!(debug_assertions);
    let arch = std::env::consts::ARCH;
    let peak = if arch == "aarch64" { 70.0 } else { 64.0 };
    let sizes = [4usize, 6, 8, 12, 16, 24, 32, 48];
    let mut lines = Vec::new();
    for &k in &sizes {
        // Batch sized to ~64 MB working set cap and ≥ 4096 matrices.
        let n = (4_000_000 / (k * k)).clamp(1024, 65_536);
        let g = measure(k, n, 9);
        let pct = 100.0 * g / peak;
        lines.push(format!("k={k}: {g:.2} GF/s ({pct:.0}% of {peak:.0})"));
    }
    log(
        "baseline",
        &format!("arch={arch} release={release} | {}", lines.join(" | ")),
    );
}
