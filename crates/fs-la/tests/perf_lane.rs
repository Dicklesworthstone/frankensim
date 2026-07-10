//! The fs-la GEMM PERF LANE (bead xdgf): packed BLIS-style gemm_f64
//! throughput against the MEASURED machine peak (fs-roofline axes) at
//! large square sizes, with the ≥75%-of-peak attainment gate. Run
//! explicitly in release:
//! `cargo test -p fs-la --release --test perf_lane -- --ignored --nocapture`
//!
//! The microkernel is the fs-simd capsule (NEON on aarch64), bitwise-
//! identical to the scalar twin — the golden 0x1d7a_a3c6_b631_7ef0 is
//! tier-invariant, verified by the gemm test suite, not here.

use fs_la::{gemm_f64, gemm_f64_parallel};
use fs_roofline::MachineAxes;

/// Best-of-3 measured GFLOP/s (2·m·n·k flops per GEMM).
fn measure(n: usize, reps: usize) -> f64 {
    let a: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.13).sin()).collect();
    let b: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.31).cos()).collect();
    let mut c = vec![0.0f64; n * n];
    gemm_f64(n, n, n, 1.0, &a, &b, 0.0, &mut c); // warm
    let mut best = f64::INFINITY;
    for _ in 0..3 {
        let t0 = std::time::Instant::now();
        for _ in 0..reps {
            gemm_f64(n, n, n, 1.0, &a, &b, 0.0, &mut c);
        }
        best = best.min(t0.elapsed().as_secs_f64() / reps as f64);
    }
    2.0 * (n * n * n) as f64 / best / 1e9
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn gemm_attainment() {
    let axes = MachineAxes::probe();
    println!(
        "{{\"metric\":\"axes\",\"cpu\":\"{}\",\"peak_single_gflops\":{:.1}}}",
        axes.cpu_brand, axes.peak_single_gflops
    );
    // Size ladder, ledgered; the gate reads the large square sizes.
    let mut large_best = 0.0f64;
    for &n in &[128usize, 256, 512, 1024] {
        let reps = (256 / (n / 128)).max(1) / 8 + 1;
        let gflops = measure(n, reps);
        let att = gflops / axes.peak_single_gflops;
        println!(
            "{{\"metric\":\"gemm-f64\",\"n\":{n},\"gflops\":{gflops:.2},\
             \"attainment_single\":{att:.3}}}"
        );
        if n >= 512 {
            large_best = large_best.max(att);
        }
    }
    // THE GATE: >= 75% of measured single-thread peak at large square
    // sizes on THIS machine (blocking constants MR=8 NR=4 KC=256
    // MC=128 NC=512). The second-ISA (x86-64/AVX) row is ARMED
    // PENDING hardware, per the recorded fleet census.
    println!(
        "{{\"metric\":\"gemm-gate\",\"attainment\":{large_best:.3},\"floor\":0.75,\
         \"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        large_best >= 0.75,
        "large-square gemm_f64 clears 75% of measured peak: {large_best:.3}"
    );
}

/// Best-of-3 all-core GFLOP/s via row-band parallel GEMM.
fn measure_parallel(n: usize, reps: usize, threads: usize) -> f64 {
    let a: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.13).sin()).collect();
    let b: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.31).cos()).collect();
    let mut c = vec![0.0f64; n * n];
    gemm_f64_parallel(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads); // warm
    let mut best = f64::INFINITY;
    for _ in 0..3 {
        let t0 = std::time::Instant::now();
        for _ in 0..reps {
            gemm_f64_parallel(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads);
        }
        best = best.min(t0.elapsed().as_secs_f64() / reps as f64);
    }
    2.0 * (n * n * n) as f64 / best / 1e9
}

/// The ALL-CORE attainment row (bead xlvx item 3): row-band parallel
/// GEMM against the measured all-core FMA axis. REPORT row by default;
/// FS_LA_ROOFLINE_GATE=1 asserts >= 0.5 (parallel GEMM leaves more on
/// the table than single-thread — memory bandwidth and band tails —
/// so the all-core floor is honest, not aspirational).
#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn gemm_attainment_all_core() {
    let threads = std::env::var("FS_LA_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(8, std::num::NonZero::get));
    let axes = MachineAxes::probe();
    println!(
        "{{\"metric\":\"axes-all-core\",\"cpu\":\"{}\",\"threads\":{threads},\"peak_all_core_gflops\":{:.1}}}",
        axes.cpu_brand, axes.peak_all_core_gflops
    );
    for n in [512usize, 1024, 2048] {
        let reps = if n >= 2048 { 1 } else { 3 };
        let g = measure_parallel(n, reps, threads);
        let att = g / axes.peak_all_core_gflops;
        println!(
            "{{\"metric\":\"gemm-f64-parallel\",\"n\":{n},\"gflops\":{g:.2},\"attainment_all_core\":{att:.3}}}"
        );
        if n == 2048 && std::env::var("FS_LA_ROOFLINE_GATE").as_deref() == Ok("1") {
            assert!(
                att >= 0.5,
                "all-core GEMM attainment {att:.3} below the 50% floor"
            );
        }
    }
}
