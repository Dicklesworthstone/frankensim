//! The fs-feec HIGH-ORDER PERF LANE (bead cwjn): sum-factorized apply
//! throughput at p = 4 against the MEASURED machine peak (fs-roofline
//! axes), plus the apply-throughput-vs-p sweep. Run explicitly in
//! release:
//! `cargo test -p fs-feec --release --test perf_lane -- --ignored --nocapture`
//!
//! GOLDEN CONSTRAINT: this lane only MEASURES the existing apply — the
//! 0xaaf1_076a_196c_6902 output golden is untouched by construction.

use fs_feec::highorder::hex::TensorSpace;
use fs_math::det;
use fs_roofline::MachineAxes;

/// FLOPs per element per apply for degree r (p = r + 1): 9 axis
/// contractions of 2·p⁴ each, plus 3·p³ accumulate adds.
fn flops_per_element(r: usize) -> f64 {
    let p = (r + 1) as f64;
    18.0 * det::powi(p, 4) + 3.0 * det::powi(p, 3)
}

fn measure_apply(m: usize, r: usize, reps: usize) -> (f64, f64) {
    let space = TensorSpace::new(m, r);
    let n = space.ndof();
    let u: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
    // Warm.
    let mut sink = space.apply_stiffness(&u)[0];
    // Best of 3 trials: the attainment claim is about machine
    // capability, so scheduler/thermal noise must not deflate it.
    let mut best = f64::INFINITY;
    for _ in 0..3 {
        let t0 = std::time::Instant::now();
        for _ in 0..reps {
            sink += space.apply_stiffness(&u)[n / 2];
        }
        best = best.min(t0.elapsed().as_secs_f64());
    }
    let elements = (m * m * m * reps) as f64;
    let gflops = elements * flops_per_element(r) / best / 1e9;
    (gflops, sink)
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn sum_factorized_attainment() {
    let axes = MachineAxes::probe();
    println!(
        "{{\"metric\":\"axes\",\"cpu\":\"{}\",\"peak_single_gflops\":{:.1}}}",
        axes.cpu_brand, axes.peak_single_gflops
    );
    // The p-sweep table (r = 1..6), ledgered.
    for r in 1..=6usize {
        let m = (48 / (r + 1)).max(6);
        let (gflops, sink) = measure_apply(m, r, 3);
        println!(
            "{{\"metric\":\"feec-apply\",\"r\":{r},\"m\":{m},\"gflops\":{gflops:.2},\
             \"attainment_single\":{:.3},\"sink\":{sink:.3}}}",
            gflops / axes.peak_single_gflops
        );
    }
    // THE GATE at p = 4 (r = 3, per the bead's p-convention: degree-4
    // tensor basis = 4 points per axis): >= 30% of measured
    // single-thread peak on THIS machine; the second-ISA gate is armed
    // and waits for x86 hardware (the rand-lane precedent).
    let (gflops, _) = measure_apply(12, 3, 6);
    let attainment = gflops / axes.peak_single_gflops;
    println!(
        "{{\"metric\":\"feec-gate\",\"r\":3,\"gflops\":{gflops:.2},\
         \"attainment\":{attainment:.3},\"floor\":0.30,\"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        attainment >= 0.30,
        "the p=4 sum-factorized apply clears 30% of measured peak: {attainment:.3}"
    );
}
