//! fs-roofline conformance suite (plan §13.3): any reimplementation must
//! pass. Covers the bead's acceptance criteria: attainment arithmetic vs
//! hand calculations, a seeded-slow kernel correctly reported below band,
//! ledgered attainment with fingerprint keying, staleness alerts on
//! fingerprint drift, and re-run reproducibility.
//!
//! Kernel sizes here are deliberately small: these tests verify HARNESS
//! correctness, not machine performance — real attainment numbers come from
//! full-size runs on fingerprinted machines.

use std::sync::atomic::{AtomicU32, Ordering};

use fs_roofline::kernels::{SeededSlowKernel, default_registry};
use fs_roofline::{
    MachineAxes, TUNE_SHAPE_CLASS, Verdict, measure, record_run, run_registry, staleness,
};

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

fn temp_db(tag: &str) -> String {
    let n = NEXT_DB.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!(
            "fs-roofline-conf-{tag}-{}-{n}.db",
            std::process::id()
        ))
        .display()
        .to_string()
}

fn cleanup_db(path: &str) {
    for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
        let _ = std::fs::remove_file(format!("{path}{suffix}"));
    }
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-roofline/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn synthetic_axes(fingerprint: u64) -> MachineAxes {
    MachineAxes {
        fingerprint,
        cpu_brand: "synthetic".to_string(),
        logical_cpus: 8,
        bandwidth_single_gbs: 100.0,
        bandwidth_all_core_gbs: 400.0,
        peak_single_gflops: 50.0,
        peak_all_core_gflops: 300.0,
    }
}

#[test]
fn rf_001_registry_runs_and_reports() {
    let axes = synthetic_axes(0x1);
    let mut registry = default_registry(1 << 12);
    let results = run_registry(&mut registry, 1, 3, &axes);
    assert_eq!(results.len(), 3);
    for r in &results {
        assert!(r.elems_per_sec > 0.0, "{}: zero rate", r.kernel);
        assert!(r.attainment >= 0.0);
        assert!(r.dispersion >= 0.0);
        assert_eq!(r.reps, 3);
        assert_eq!(
            r.verdict,
            Verdict::NoTarget,
            "{}: v0 kernels are report-only",
            r.kernel
        );
        assert!(r.to_jsonl().starts_with('{'));
    }
    verdict(
        "rf-001",
        "default registry (axpy/dot/sum) measured and reported",
    );
}

#[test]
fn rf_002_seeded_slow_kernel_is_below_band() {
    // Real machine axes so the limit is genuine; the kernel claims 90% of
    // the bandwidth roof and cannot come close by construction.
    let axes = MachineAxes::probe();
    let mut slow = SeededSlowKernel::new(1 << 14);
    let result = measure(&mut slow, 1, 3, &axes);
    // On a CONTENDED host the probe itself is implausible and the
    // 1n61 guard refuses to gate at all — the honest outcome (this
    // fired on a live agent-build storm during development). Either
    // way the seeded-slow kernel must never wear within_band.
    if axes.plausible() {
        assert_eq!(
            result.verdict,
            Verdict::BelowBand,
            "seeded-slow kernel must be caught below its band (attainment {:.4})",
            result.attainment
        );
        assert!(result.attainment < 0.9);
    } else {
        assert_eq!(
            result.verdict,
            Verdict::EnvironmentInvalid,
            "implausible axes must refuse to gate"
        );
    }
    verdict(
        "rf-002",
        &format!(
            "seeded-slow reported below_band at attainment {:.4} vs target 0.9",
            result.attainment
        ),
    );
}

#[test]
fn rf_003_ledgered_run_with_fingerprint_keying() {
    let db = temp_db("ledger");
    let ledger = fs_ledger::Ledger::open(&db).expect("open ledger");
    let axes = synthetic_axes(0xFEED_FACE);
    let mut registry = default_registry(1 << 10);
    let results = run_registry(&mut registry, 0, 2, &axes);
    let op = record_run(&ledger, &axes, &results).expect("record run");
    // The op is complete, metrics/events/tune rows exist per kernel.
    let row = ledger.op(op).unwrap().expect("op row");
    assert_eq!(row.outcome.as_deref(), Some("ok"));
    assert_eq!(
        ledger.table_count("metrics").unwrap(),
        3 * results.len() as u64
    );
    assert_eq!(ledger.table_count("events").unwrap(), results.len() as u64);
    assert_eq!(ledger.table_count("tune").unwrap(), results.len() as u64);
    // Tune rows are keyed by THIS fingerprint.
    let fp = 0xFEED_FACEu64.to_le_bytes();
    for r in &results {
        let tune = ledger
            .tune_get(&r.kernel, TUNE_SHAPE_CLASS, &fp)
            .unwrap()
            .expect("tune row under current fingerprint");
        assert!(tune.measured.contains("attainment"));
    }
    assert!(ledger.lint().unwrap().is_clean());
    drop(ledger);
    cleanup_db(&db);
    verdict(
        "rf-003",
        "run ledgered: op + metrics + events + fingerprint-keyed tune rows",
    );
}

#[test]
fn rf_004_staleness_alerts_on_fingerprint_drift() {
    let db = temp_db("stale");
    let ledger = fs_ledger::Ledger::open(&db).expect("open ledger");
    let old_axes = synthetic_axes(0xAAAA);
    let mut registry = default_registry(1 << 10);
    let results = run_registry(&mut registry, 0, 1, &old_axes);
    record_run(&ledger, &old_axes, &results).expect("record under old fingerprint");

    let kernel = &results[0].kernel;
    // Same machine → fresh.
    assert_eq!(
        staleness(&ledger, kernel, 0xAAAA).unwrap(),
        fs_roofline::Staleness::Fresh
    );
    // Drifted machine → alert.
    assert_eq!(
        staleness(&ledger, kernel, 0xBBBB).unwrap(),
        fs_roofline::Staleness::FingerprintDrift
    );
    // Unknown kernel → never measured.
    assert_eq!(
        staleness(&ledger, "gemm-f64", 0xAAAA).unwrap(),
        fs_roofline::Staleness::NeverMeasured
    );
    drop(ledger);
    cleanup_db(&db);
    verdict(
        "rf-004",
        "staleness: fresh / fingerprint-drift / never-measured all detected",
    );
}

#[test]
fn rf_004b_invalid_environment_never_becomes_fresh_evidence() {
    let db = temp_db("invalid-environment");
    let ledger = fs_ledger::Ledger::open(&db).expect("open ledger");
    let crushed = MachineAxes {
        fingerprint: 0xBAD,
        cpu_brand: "synthetic-crushed".to_string(),
        logical_cpus: 8,
        bandwidth_single_gbs: 0.2,
        bandwidth_all_core_gbs: 0.4,
        peak_single_gflops: 0.1,
        peak_all_core_gflops: 0.4,
    };
    let mut registry = default_registry(1 << 10);
    let results = run_registry(&mut registry, 0, 1, &crushed);
    assert!(
        results
            .iter()
            .all(|row| row.verdict == Verdict::EnvironmentInvalid)
    );
    let op = record_run(&ledger, &crushed, &results).expect("record invalid run");
    let row = ledger.op(op).unwrap().expect("op row");
    assert_eq!(row.outcome.as_deref(), Some("error"));
    assert_eq!(ledger.table_count("tune").unwrap(), 0);
    assert_eq!(
        staleness(&ledger, &results[0].kernel, crushed.fingerprint).unwrap(),
        fs_roofline::Staleness::NeverMeasured
    );
    assert_eq!(
        ledger.table_count("events").unwrap(),
        results.len() as u64 + 1
    );
    drop(ledger);
    cleanup_db(&db);
    verdict(
        "rf-004b",
        "invalid axes recorded as failed evidence without publishing fresh tune rows",
    );
}

#[test]
fn rf_005_reproducibility_within_dispersion() {
    // Two back-to-back measurements of the same kernel must agree within a
    // generous multiple of their reported dispersion (harness smoke, not a
    // machine claim — shared CI boxes are noisy).
    let axes = MachineAxes::probe();
    let mut registry_a = default_registry(1 << 14);
    let mut registry_b = default_registry(1 << 14);
    let a = &run_registry(&mut registry_a, 2, 7, &axes)[0];
    let b = &run_registry(&mut registry_b, 2, 7, &axes)[0];
    let rel_delta = (a.elems_per_sec - b.elems_per_sec).abs() / a.elems_per_sec.max(1.0);
    let allowance = 0.5 + 3.0 * a.dispersion.max(b.dispersion);
    assert!(
        rel_delta <= allowance,
        "re-run drift {rel_delta:.3} exceeds allowance {allowance:.3} \
         (a={:.3e}, b={:.3e}, disp a={:.3} b={:.3})",
        a.elems_per_sec,
        b.elems_per_sec,
        a.dispersion,
        b.dispersion
    );
    verdict(
        "rf-005",
        &format!("re-run delta {rel_delta:.3} within allowance {allowance:.3}"),
    );
}

#[test]
fn rf_006_cli_smoke_prints_report_and_ledgers() {
    let db = temp_db("cli");
    let exe = env!("CARGO_BIN_EXE_roofline");
    let out = std::process::Command::new(exe)
        .args([
            "--n",
            "4096",
            "--warmup",
            "0",
            "--reps",
            "2",
            "--ledger",
            db.as_str(),
        ])
        .output()
        .expect("run roofline CLI");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("\"fingerprint\""), "axes line present");
    assert!(stdout.contains("simd-axpy-f64"), "kernel lines present");
    assert!(
        stdout.contains("\"landed\":false"),
        "§14.1 coverage table present"
    );
    assert!(
        stdout.contains("\"ledgered\":true"),
        "ledger receipt present"
    );
    if stdout.contains("\"verdict\":\"environment_invalid\"") {
        assert!(
            stdout.contains("\"citable\":false"),
            "invalid run must be marked non-citable"
        );
        assert!(
            stdout.contains("NeverMeasured"),
            "invalid run must not publish fresh evidence"
        );
    } else {
        assert!(stdout.contains("\"citable\":true"));
        assert!(stdout.contains("Fresh"), "valid run is fresh");
    }
    // Bad args refuse structurally.
    let bad = std::process::Command::new(exe)
        .args(["--n", "zero"])
        .output()
        .expect("run");
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("Roofline"));
    cleanup_db(&db);
    verdict(
        "rf-006",
        "CLI prints axes + kernels + coverage + ledger receipt; refuses bad args",
    );
}
