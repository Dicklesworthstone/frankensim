//! The GEMM autotune-loop drills (bead yqug): measure → cache → model →
//! dispatch closed end-to-end, with the failure paths exercised one by
//! one — cold start, ledger warm start, stale row, invalid row, cancel,
//! pin failure, and pinned replay. The oracle-tolerance lane runs
//! explicitly (`--ignored`, release) like every wall-clock gate.

use fs_exec::{CancelGate, GemmBlockPlan, TuneSource, Tuner};
use fs_session::{GemmTuneError, gemm_f64_session, gemm_kernel_key, gemm_shape_class};

const FP_THIS: u64 = 0x00AA_11BB_22CC_33DD;
const THREADS: usize = 4;

/// A problem big enough to exercise the parallel path (m >= 256).
const M: usize = 320;
const N: usize = 288;
const K: usize = 300;

fn fill(buf: &mut [f64], salt: u64) {
    for (i, slot) in buf.iter_mut().enumerate() {
        let mut z = (i as u64)
            .wrapping_add(salt)
            .wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z ^= z >> 31;
        *slot = (z >> 11) as f64 / 9_007_199_254_740_992.0 - 0.5;
    }
}

fn problem() -> (Vec<f64>, Vec<f64>) {
    let mut a = vec![0.0f64; M * K];
    let mut b = vec![0.0f64; K * N];
    fill(&mut a, 1);
    fill(&mut b, 2);
    (a, b)
}

fn serial_reference(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut c = vec![0.0f64; M * N];
    fs_la::gemm_f64(M, N, K, 1.0, a, b, 0.0, &mut c);
    c
}

fn temp_ledger(tag: &str) -> (std::path::PathBuf, fs_ledger::Ledger) {
    let dir = std::env::temp_dir().join(format!("fs-gemm-tune-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let ledger =
        fs_ledger::Ledger::open(dir.join("tune.led").to_str().expect("utf8")).expect("ledger");
    (dir, ledger)
}

/// COLD DRILL: no pins, no rows, no ledger — the loop must sweep once,
/// record a tuned row, dispatch bitwise-identically to the serial
/// reference, and NOT re-sweep on the second call.
#[test]
fn cold_start_sweeps_once_and_matches_serial_bits() {
    let (a, b) = problem();
    let reference = serial_reference(&a, &b);
    let mut tuner = Tuner::cold(FP_THIS);
    let gate = CancelGate::new();
    let mut c = vec![f64::NAN; M * N]; // beta = 0 must overwrite garbage
    let first = gemm_f64_session(
        &mut tuner, None, &gate, THREADS, M, N, K, 1.0, &a, &b, 0.0, &mut c,
    )
    .expect("cold dispatch");
    assert!(first.swept, "cold start measures");
    assert_eq!(first.source, TuneSource::Tuned);
    assert_eq!(first.kernel, gemm_kernel_key());
    assert_eq!(first.shape_class, gemm_shape_class(M, N, K));
    assert_eq!(
        c.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        reference.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        "autotuned dispatch is bitwise the serial contract"
    );
    // Second call: row cached, no sweep, same plan.
    let mut c2 = vec![0.0f64; M * N];
    let second = gemm_f64_session(
        &mut tuner, None, &gate, THREADS, M, N, K, 1.0, &a, &b, 0.0, &mut c2,
    )
    .expect("warm dispatch");
    assert!(!second.swept, "cached row answers the second call");
    assert_eq!(second.plan, first.plan);
    assert_eq!(second.source, TuneSource::Tuned);
    // Decisions were recorded for study pinning.
    assert!(
        tuner
            .decisions()
            .iter()
            .any(|d| d.kernel.starts_with(&gemm_kernel_key()) && d.params == first.plan.canonical())
    );
}

/// LEDGER WARM-START DRILL: session 1 sweeps and writes through; a fresh
/// tuner on the SAME machine seeds from the ledger and never re-measures.
#[test]
fn ledger_cache_warm_starts_a_fresh_session() {
    let (dir, ledger) = temp_ledger("warm");
    let (a, b) = problem();
    let gate = CancelGate::new();
    let mut c = vec![0.0f64; M * N];
    let mut tuner1 = Tuner::cold(FP_THIS);
    let first = gemm_f64_session(
        &mut tuner1,
        Some(&ledger),
        &gate,
        THREADS,
        M,
        N,
        K,
        1.0,
        &a,
        &b,
        0.0,
        &mut c,
    )
    .expect("session 1");
    assert!(first.swept);

    let mut tuner2 = Tuner::cold(FP_THIS);
    let mut c2 = vec![0.0f64; M * N];
    let second = gemm_f64_session(
        &mut tuner2,
        Some(&ledger),
        &gate,
        THREADS,
        M,
        N,
        K,
        1.0,
        &a,
        &b,
        0.0,
        &mut c2,
    )
    .expect("session 2");
    assert!(!second.swept, "the ledger row seeds the fresh session");
    assert_eq!(second.source, TuneSource::Tuned);
    assert_eq!(second.plan, first.plan, "cache replays the same plan");
    let _ = std::fs::remove_dir_all(&dir);
}

/// STALE DRILL: a ledger row whose embedded machine fingerprint differs
/// from this tuner's is refused at adoption and the loop re-measures —
/// another machine's timings never dispatch here.
#[test]
fn stale_other_machine_row_is_refused_and_remeasured() {
    let (dir, ledger) = temp_ledger("stale");
    let (a, b) = problem();
    let gate = CancelGate::new();
    let kernel = gemm_kernel_key();
    let class = gemm_shape_class(M, N, K);
    // Session on machine A records the row under A's fingerprint...
    let fp_other = 0x0DEA_D0BE_EF00_0001_u64;
    let mut tuner_a = Tuner::cold(fp_other);
    let mut c = vec![0.0f64; M * N];
    gemm_f64_session(
        &mut tuner_a,
        None,
        &gate,
        THREADS,
        M,
        N,
        K,
        1.0,
        &a,
        &b,
        0.0,
        &mut c,
    )
    .expect("machine A session");
    let row_a = tuner_a.row_json(&kernel, &class).expect("row recorded");
    // ...but the cache entry lands under THIS machine's key (a corrupt
    // or copied cache — exactly what the fingerprint check exists for).
    ledger
        .tune_put(
            &kernel,
            &class,
            &FP_THIS.to_le_bytes(),
            "\"gemm-block-plan\"",
            &row_a,
        )
        .expect("plant stale row");
    let mut tuner = Tuner::cold(FP_THIS);
    let report = gemm_f64_session(
        &mut tuner,
        Some(&ledger),
        &gate,
        THREADS,
        M,
        N,
        K,
        1.0,
        &a,
        &b,
        0.0,
        &mut c,
    )
    .expect("stale-cache session");
    assert!(
        report.swept,
        "the stale row is refused; the loop re-measures"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// INVALID DRILL: a cache row that is not a canonical tune-row line is
/// refused (fail closed) and the loop re-measures.
#[test]
fn invalid_cache_row_is_refused_and_remeasured() {
    let (dir, ledger) = temp_ledger("invalid");
    let (a, b) = problem();
    let gate = CancelGate::new();
    let kernel = gemm_kernel_key();
    let class = gemm_shape_class(M, N, K);
    ledger
        .tune_put(
            &kernel,
            &class,
            &FP_THIS.to_le_bytes(),
            "\"gemm-block-plan\"",
            "{\"not\":\"a tune row\"}",
        )
        .expect("plant invalid row");
    let mut tuner = Tuner::cold(FP_THIS);
    let mut c = vec![0.0f64; M * N];
    let report = gemm_f64_session(
        &mut tuner,
        Some(&ledger),
        &gate,
        THREADS,
        M,
        N,
        K,
        1.0,
        &a,
        &b,
        0.0,
        &mut c,
    )
    .expect("invalid-cache session");
    assert!(
        report.swept,
        "garbage never dispatches; the loop re-measures"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// CANCEL DRILL: a requested gate aborts the sweep BEFORE any candidate
/// runs — structured error, no row recorded, output untouched.
#[test]
fn cancelled_sweep_records_nothing_and_leaves_c_untouched() {
    let (a, b) = problem();
    let gate = CancelGate::new();
    gate.request();
    let mut tuner = Tuner::cold(FP_THIS);
    let sentinel = 42.5f64;
    let mut c = vec![sentinel; M * N];
    let err = gemm_f64_session(
        &mut tuner, None, &gate, THREADS, M, N, K, 1.0, &a, &b, 0.0, &mut c,
    )
    .expect_err("cancelled");
    assert!(matches!(err, GemmTuneError::Cancelled));
    assert!(
        !tuner.has_row(&gemm_kernel_key(), &gemm_shape_class(M, N, K)),
        "no partial evidence survives a cancel"
    );
    assert!(
        c.iter().all(|&x| x == sentinel),
        "nothing dispatched under an unselected plan"
    );
}

/// PIN-FAILURE DRILL: non-canonical params and mis-keyed pins fail
/// closed with structured errors instead of resolving to a default that
/// would be falsely recorded as pinned.
#[test]
fn pin_failures_are_structured_refusals() {
    let mut tuner = Tuner::cold(FP_THIS);
    let kernel = gemm_kernel_key();
    for bad in [
        "mc=banana,nc-cap=2048",
        "mc=33,nc-cap=2048",   // off-lattice mc
        "mc=32,nc-cap=100000", // unbounded nc_cap
        "mc=32",               // missing member
        "mc=032,nc-cap=2048",  // non-canonical spelling
    ] {
        assert!(tuner.pin(&kernel, bad).is_err(), "{bad:?} must be refused");
    }
    // A gemm plan pinned under a non-gemm kernel key is refused too.
    assert!(
        tuner
            .pin_gemm_blocking("stencil7-f32", GemmBlockPlan::COLD_START)
            .is_err()
    );
    // And the canonical spelling round-trips through the replay path.
    tuner
        .pin(&kernel, "mc=64,nc-cap=512")
        .expect("canonical pin");
}

/// REPLAY DRILL: pinning the recorded decision reproduces the plan with
/// Pinned provenance, runs NO measurement, and produces the same bits.
#[test]
fn pinned_replay_skips_measurement_and_reproduces_bits() {
    let (a, b) = problem();
    let gate = CancelGate::new();
    let mut tuner1 = Tuner::cold(FP_THIS);
    let mut c1 = vec![0.0f64; M * N];
    let live = gemm_f64_session(
        &mut tuner1,
        None,
        &gate,
        THREADS,
        M,
        N,
        K,
        1.0,
        &a,
        &b,
        0.0,
        &mut c1,
    )
    .expect("live session");
    // Replay: a second session pins the recorded decision's params.
    let mut tuner2 = Tuner::cold(FP_THIS);
    tuner2
        .pin(gemm_kernel_key(), live.plan.canonical())
        .expect("replay pin");
    let mut c2 = vec![0.0f64; M * N];
    let replay = gemm_f64_session(
        &mut tuner2,
        None,
        &gate,
        THREADS,
        M,
        N,
        K,
        1.0,
        &a,
        &b,
        0.0,
        &mut c2,
    )
    .expect("replay session");
    assert!(!replay.swept, "a pinned path never re-measures");
    assert_eq!(replay.source, TuneSource::Pinned);
    assert_eq!(replay.plan, live.plan);
    assert_eq!(
        c1.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        c2.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        "replay reproduces the live bits exactly"
    );
}

/// Every plan on the sweep lattice dispatches bitwise-identically to the
/// serial contract (the bead's bit-neutral-loop clause, checked from the
/// consumer side on an awkward non-multiple shape).
#[test]
fn every_lattice_plan_matches_serial_bits() {
    let (a, b) = problem();
    let reference = serial_reference(&a, &b);
    for mc in [16usize, 32, 64, 128] {
        for nc_cap in [512usize, 2048] {
            let mut c = vec![0.0f64; M * N];
            fs_la::gemm_f64_parallel_with(
                M,
                N,
                K,
                1.0,
                &a,
                &b,
                0.0,
                &mut c,
                THREADS,
                mc,
                N.min(nc_cap),
            );
            assert_eq!(
                c.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                reference.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                "mc={mc} nc_cap={nc_cap}"
            );
        }
    }
}

/// ORACLE LANE (wall-clock; run explicitly in release):
/// `cargo test -p fs-session --release --test gemm_tune -- --ignored --nocapture`
///
/// The live loop's selected plan must be within the declared tolerance
/// of the EXHAUSTIVE oracle: every lattice candidate re-measured at the
/// real problem size with best-of-3. Tolerance is declared here, not
/// fitted: 25% — candidate timing gaps on the measured xlvx landscape
/// exceed 2x, so a quarter margin separates model error from noise
/// while staying honest on a shared host. Second-ISA status: armed;
/// this lane reports its gate on every machine that runs it and the
/// x86 counterpart runs when an x86 host picks it up (the rand-lane
/// precedent).
#[test]
#[ignore = "wall-clock oracle lane: run explicitly in release with --ignored"]
fn live_choice_is_within_declared_tolerance_of_exhaustive_oracle() {
    const TOLERANCE: f64 = 0.25;
    let (a, b) = problem();
    let gate = CancelGate::new();
    let mut tuner = Tuner::cold(FP_THIS);
    let mut c = vec![0.0f64; M * N];
    let live = gemm_f64_session(
        &mut tuner, None, &gate, THREADS, M, N, K, 1.0, &a, &b, 0.0, &mut c,
    )
    .expect("live session");
    // Exhaustive oracle at the REAL size, best-of-3 per candidate.
    let mut oracle: Vec<(u64, String)> = Vec::new();
    for mc in [16usize, 32, 64, 128] {
        for nc_cap in [512usize, 2048] {
            let mut best = u64::MAX;
            for _ in 0..3 {
                c.fill(0.0);
                let t0 = std::time::Instant::now();
                fs_la::gemm_f64_parallel_with(
                    M,
                    N,
                    K,
                    1.0,
                    &a,
                    &b,
                    0.0,
                    &mut c,
                    THREADS,
                    mc,
                    N.min(nc_cap),
                );
                best = best.min(t0.elapsed().as_nanos() as u64);
            }
            oracle.push((best, format!("mc={mc},nc-cap={nc_cap}")));
        }
    }
    oracle.sort();
    let oracle_best = &oracle[0];
    let live_entry = oracle
        .iter()
        .find(|(_, params)| *params == live.plan.canonical())
        .expect("the live plan is on the lattice");
    let ratio = live_entry.0 as f64 / oracle_best.0 as f64;
    println!(
        "{{\"metric\":\"gemm-autotune-oracle\",\"live\":\"{}\",\"oracle_best\":\"{}\",\
         \"ratio\":{ratio:.3},\"tolerance\":{TOLERANCE},\"machine\":\"{}-{}\"}}",
        live.plan.canonical(),
        oracle_best.1,
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        ratio <= 1.0 + TOLERANCE,
        "live plan {} is {ratio:.3}x the oracle best {} (tolerance {TOLERANCE})",
        live.plan.canonical(),
        oracle_best.1
    );
}
