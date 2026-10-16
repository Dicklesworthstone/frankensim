//! SIMD tier-selection AUDIT (bead fz2.2): the chosen dispatch tier
//! must be at least as fast as the always-available scalar reference,
//! PER KERNEL — "chosen tier = fastest tier, verified", not assumed
//! from the ISA name. Rows report the measured ratio for the ledger;
//! the gate flags any kernel where the chosen tier LOSES to scalar
//! (a wrong dispatch choice, e.g. an intrinsic path the autovectorizer
//! already beats).
//! Run: `cargo test -p fs-simd --release --test tier_audit -- --ignored --nocapture`

use std::time::Instant;

/// Best-of-5 wall seconds for `reps` calls of `f`.
fn best(reps: usize, mut f: impl FnMut()) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..5 {
        let t0 = Instant::now();
        for _ in 0..reps {
            f();
        }
        best = best.min(t0.elapsed().as_secs_f64());
    }
    best
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn chosen_tier_beats_scalar_per_kernel() {
    let ops = fs_simd::ops();
    let n = 1 << 14; // L1/L2-resident: measures the kernel, not DRAM
    let x: Vec<f64> = (0..n).map(|i| (i as f64) * 1e-4 + 0.5).collect();
    let y0: Vec<f64> = (0..n).map(|i| (i as f64).mul_add(-3e-5, 1.0)).collect();
    let z: Vec<f64> = (0..n).map(|i| (i as f64) * 7e-5 - 0.25).collect();
    let reps = 2000;

    let mut rows: Vec<(&str, f64, f64)> = Vec::new();

    // axpy
    let mut y = y0.clone();
    let chosen = best(reps, || (ops.axpy)(1.000_1, &x, &mut y));
    let mut y = y0.clone();
    let scalar = best(reps, || fs_simd::scalar::axpy(1.000_1, &x, &mut y));
    rows.push(("axpy", scalar, chosen));

    // scale
    let mut y = y0.clone();
    let chosen = best(reps, || (ops.scale)(1.000_000_1, &mut y));
    let mut y = y0.clone();
    let scalar = best(reps, || fs_simd::scalar::scale(1.000_000_1, &mut y));
    rows.push(("scale", scalar, chosen));

    // mul_elem
    let mut out = vec![0.0f64; n];
    let chosen = best(reps, || (ops.mul_elem)(&x, &y0, &mut out));
    let scalar = best(reps, || fs_simd::scalar::mul_elem(&x, &y0, &mut out));
    rows.push(("mul_elem", scalar, chosen));

    // fma3
    let chosen = best(reps, || (ops.fma3)(&x, &y0, &z, &mut out));
    let scalar = best(reps, || fs_simd::scalar::fma3(&x, &y0, &z, &mut out));
    rows.push(("fma3", scalar, chosen));

    // dot
    let chosen = best(reps, || {
        std::hint::black_box((ops.dot)(&x, &y0));
    });
    let scalar = best(reps, || {
        std::hint::black_box(fs_simd::scalar::dot(&x, &y0));
    });
    rows.push(("dot", scalar, chosen));

    // sum
    let chosen = best(reps, || {
        std::hint::black_box((ops.sum)(&x));
    });
    let scalar = best(reps, || {
        std::hint::black_box(fs_simd::scalar::sum(&x));
    });
    rows.push(("sum", scalar, chosen));

    // mk8x4 (the GEMM register microkernel — the load-bearing one)
    let kc = 256;
    let a_panel: Vec<f64> = (0..8 * kc).map(|i| (i as f64) * 1e-5).collect();
    let b_panel: Vec<f64> = (0..4 * kc).map(|i| (i as f64) * -2e-5 + 0.1).collect();
    let mut acc = [[0.0f64; 4]; 8];
    let chosen = best(reps, || (ops.mk8x4_f64)(&a_panel, &b_panel, kc, &mut acc));
    let mut acc2 = [[0.0f64; 4]; 8];
    let scalar = best(reps, || {
        fs_simd::scalar::mk8x4_f64(&a_panel, &b_panel, kc, &mut acc2);
    });
    std::hint::black_box((acc, acc2));
    rows.push(("mk8x4_f64", scalar, chosen));

    let mut losers = Vec::new();
    for (kernel, scalar_s, chosen_s) in &rows {
        let ratio = scalar_s / chosen_s.max(1e-12);
        println!(
            "{{\"metric\":\"tier-audit\",\"tier\":\"{:?}\",\"kernel\":\"{kernel}\",\
             \"chosen_over_scalar\":{ratio:.2}}}",
            ops.tier
        );
        // 0.9: the chosen tier may TIE scalar within noise (identical
        // codegen), but a real loss means the dispatch picked wrong.
        if ratio < 0.9 {
            losers.push((*kernel, ratio));
        }
    }
    assert!(
        losers.is_empty(),
        "chosen tier {:?} LOSES to scalar on {losers:?} — dispatch is not picking the \
         fastest tier for these kernels",
        ops.tier
    );
}
