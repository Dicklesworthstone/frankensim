//! The fs-feec HIGH-ORDER PERF LANE (bead cwjn): sum-factorized apply
//! throughput at p = 4 against the MEASURED machine peak (fs-roofline
//! axes), plus the apply-throughput-vs-p sweep. Run explicitly in
//! release:
//! `FRANKENSIM_BASELINE_STORE=<jsonl> FRANKENSIM_FIRMWARE_ID=<id> cargo test
//! -p fs-feec --release --test perf_lane -- --ignored --nocapture`
//!
//! GOLDEN CONSTRAINT: this lane only MEASURES the existing apply — the
//! 0xaaf1_076a_196c_6902 output golden is untouched by construction.

use fs_feec::highorder::hex::TensorSpace;
use fs_math::det;
use fs_roofline::{
    AxisBaselinePolicy, BaselineIdentity, BaselineStore, MachineAxes, days_since_epoch_now,
};

fn fail_invalid_environment(reason: &str, attainment: Option<f64>) -> ! {
    match attainment {
        Some(value) => println!(
            "{{\"metric\":\"feec-gate\",\"verdict\":\"environment_invalid\",\
             \"reason\":\"{reason}\",\"attainment\":{value:.3},\
             \"machine\":\"{}-{}\"}}",
            std::env::consts::OS,
            std::env::consts::ARCH
        ),
        None => println!(
            "{{\"metric\":\"feec-gate\",\"verdict\":\"environment_invalid\",\
             \"reason\":\"{reason}\",\"attainment\":null,\
             \"machine\":\"{}-{}\"}}",
            std::env::consts::OS,
            std::env::consts::ARCH
        ),
    }
    panic!("FEEC roofline evidence rejected: {reason}");
}

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
    println!("{{\"metric\":\"axes-pre\",\"axes\":{}}}", axes.to_jsonl());
    // Environment validity (bead 1n61): implausible axes mean the probe
    // itself was contaminated (contended/throttled machine), so BOTH the
    // numerator and denominator of attainment are garbage — refuse to
    // gate rather than emit a vacuous pass or a false failure.
    if let Some(reason) = axes.plausibility_error() {
        fail_invalid_environment(reason, None);
    }
    let baseline_path = std::env::var("FRANKENSIM_BASELINE_STORE")
        .unwrap_or_else(|_| panic!("FRANKENSIM_BASELINE_STORE is required for a citable gate"));
    let firmware = std::env::var("FRANKENSIM_FIRMWARE_ID")
        .unwrap_or_else(|_| panic!("FRANKENSIM_FIRMWARE_ID is required for a citable gate"));
    let baseline_text = std::fs::read_to_string(&baseline_path)
        .unwrap_or_else(|error| panic!("cannot read baseline store {baseline_path:?}: {error}"));
    let baseline_store = BaselineStore::from_jsonl(&baseline_text)
        .unwrap_or_else(|error| panic!("invalid baseline store: {error}"));
    let identity = BaselineIdentity::current(&axes, firmware)
        .unwrap_or_else(|error| panic!("invalid baseline identity: {error}"));
    let now_day = days_since_epoch_now()
        .unwrap_or_else(|error| panic!("cannot establish baseline age: {error}"));
    let baseline_policy = AxisBaselinePolicy::new(
        baseline_store.for_fingerprint(axes.fingerprint),
        &identity,
        now_day,
    );
    // The p-sweep table (r = 1..6), ledgered.
    for r in 1..=6usize {
        let m = (48 / (r + 1)).max(6);
        let (gflops, sink) = measure_apply(m, r, 3);
        if !gflops.is_finite() || gflops < 0.0 {
            fail_invalid_environment("non-finite or negative FEEC throughput", None);
        }
        println!(
            "{{\"metric\":\"feec-apply\",\"r\":{r},\"m\":{m},\"gflops\":{gflops:.2},\
             \"attainment_single\":{:.3},\"sink\":{sink:.3}}}",
            gflops / axes.peak_single_gflops
        );
    }
    // THE GATE at p = 4 (r = 3, per the bead's p-convention: degree-4
    // tensor basis = 4 points per axis): >= 30% of measured
    // single-thread peak on THIS machine. Bead cwjn requires separately
    // admitted rows on both reference ISAs before the cross-ISA claim lands.
    let (gflops, _) = measure_apply(12, 3, 6);
    let attainment = gflops / axes.peak_single_gflops;
    if !attainment.is_finite() || attainment < 0.0 {
        fail_invalid_environment("non-finite or negative FEEC attainment", None);
    }
    println!(
        "{{\"metric\":\"feec-gate\",\"r\":3,\"gflops\":{gflops:.2},\
         \"attainment\":{attainment:.3},\"floor\":0.30,\"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    let post_axes = MachineAxes::probe();
    println!(
        "{{\"metric\":\"axes-post\",\"axes\":{}}}",
        post_axes.to_jsonl()
    );
    let baseline_verdict = baseline_policy.verdict(&axes, &post_axes);
    println!("{}", baseline_policy.receipt_json(&axes, &post_axes));
    if !baseline_verdict.trusted() {
        fail_invalid_environment("historical baseline admission rejected", Some(attainment));
    }
    // Over-roof poisoning (bead 1n61): a kernel "beating" the probed
    // peak by >1.5x means the PEAK probe was contaminated, not that the
    // kernel is fast — the whole measurement is invalid, both directions.
    if attainment > 1.5 {
        fail_invalid_environment(
            "attainment exceeds the credible roofline band",
            Some(attainment),
        );
    }
    assert!(
        attainment >= 0.30,
        "the p=4 sum-factorized apply clears 30% of measured peak: {attainment:.3}"
    );
}
