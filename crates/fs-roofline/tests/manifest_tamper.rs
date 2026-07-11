//! Adversarial battery for the op-bound ordered result manifest (bead
//! gp3.15): a writer holding ordinary Ledger APIs must not be able to alter
//! recorded roofline evidence while keeping the finalized run receipt.
//!
//! The named attack: replace one result payload plus its matching artifact
//! and params (all self-consistent, edge included) while retaining the old
//! run receipt. Before the manifest, `validate_roofline_row` accepted that
//! row and staleness classified the altered evidence as Fresh. Now staleness
//! recomputes the whole receipt from the manifest and the stored rows, so
//! every tamper shape below must classify as CorruptEvidence — and honest
//! rerun history must stay Fresh.

use std::sync::atomic::{AtomicU32, Ordering};

use fs_ledger::{EdgeRole, Ledger};
use fs_roofline::kernels::default_registry;
use fs_roofline::{
    AxisBaselinePolicy, BaselineAxes, BaselineCandidate, BaselineIdentity, MachineAxes, Staleness,
    finalize_registry_tuning, promote_baseline, record_run, run_registry, staleness_at,
};

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

const FINGERPRINT: u64 = 0xBEEF;

fn temp_db(tag: &str) -> String {
    let n = NEXT_DB.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!(
            "fs-roofline-tamper-{tag}-{}-{n}.db",
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

fn trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
    let identity =
        BaselineIdentity::current(axes, "test-firmware").expect("valid synthetic identity");
    let candidates: Vec<_> = (0_u64..3)
        .map(|ordinal| {
            BaselineCandidate::from_receipt(
                axes.clone(),
                identity.clone(),
                fs_blake3::hash_domain(
                    "fs-roofline.tamper-baseline-source.v1",
                    &ordinal.to_le_bytes(),
                ),
            )
            .expect("valid synthetic candidate")
        })
        .collect();
    let baseline = promote_baseline(
        &candidates,
        "test-operator",
        "deterministic tamper fixture",
        20_000,
        90,
    )
    .expect("valid synthetic baseline");
    (baseline, identity)
}

struct RecordedRun {
    ledger: Ledger,
    baseline: BaselineAxes,
    kernels: Vec<(String, String)>,
    recorded_at: i64,
}

/// Run the default registry once and record it; returns everything a tamper
/// test needs to read rows back and probe staleness.
fn recorded_run(db: &str) -> RecordedRun {
    let ledger = Ledger::open(db).expect("open ledger");
    let axes = synthetic_axes(FINGERPRINT);
    let (baseline, identity) = trusted_baseline(&axes);
    let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
    let mut registry = default_registry(1 << 10);
    let mut results = run_registry(&mut registry, 0, 1, &axes);
    let mut finalized = finalize_registry_tuning(&mut registry, &axes, &axes, policy, &results)
        .expect("finalize run");
    let op = record_run(&ledger, &axes, &axes, policy, &mut finalized, &mut results)
        .expect("record run");
    let recorded_at = ledger
        .op(op)
        .unwrap()
        .expect("recorded op")
        .t_end
        .expect("finished op");
    let kernels = results
        .iter()
        .map(|r| (r.kernel.clone(), r.version.clone()))
        .collect();
    RecordedRun {
        ledger,
        baseline,
        kernels,
        recorded_at,
    }
}

fn probe(run: &RecordedRun, kernel: &str, version: &str) -> Staleness {
    staleness_at(
        &run.ledger,
        kernel,
        version,
        FINGERPRINT,
        Some(run.baseline.content_hash()),
        run.recorded_at + 1,
    )
    .expect("staleness probe")
}

/// The stored roofline row for one kernel (there is exactly one per run).
fn roofline_row(ledger: &Ledger, kernel: &str) -> fs_ledger::TuneRow {
    let rows = ledger.tune_rows(kernel).expect("tune rows");
    let mut roofline: Vec<_> = rows
        .into_iter()
        .filter(|r| r.shape_class.contains(":run="))
        .collect();
    assert_eq!(roofline.len(), 1, "expected one roofline row for {kernel}");
    roofline.pop().expect("row")
}

/// Self-consistently replace `row`'s payload with `new_measured`: store the
/// artifact, add the op edge, and rewrite params to cite the new artifact —
/// exactly what a malicious Ledger writer can do without touching the op row.
fn splice_payload(ledger: &Ledger, row: &fs_ledger::TuneRow, new_measured: &str) {
    let old_hash = fs_ledger::hash_bytes(row.measured.as_bytes()).to_string();
    let new_hash = fs_ledger::hash_bytes(new_measured.as_bytes());
    let artifact = ledger
        .put_artifact(
            "roofline-benchmark-result",
            new_measured.as_bytes(),
            Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}"),
        )
        .expect("store forged artifact");
    assert_eq!(artifact.hash, new_hash);
    let op: i64 = row
        .params
        .split_once("\"op\":")
        .and_then(|(_, rest)| rest.split_once(','))
        .and_then(|(digits, _)| digits.parse().ok())
        .expect("op id in params");
    ledger
        .link(op, &new_hash, EdgeRole::Out)
        .expect("forged edge");
    let forged_params = row.params.replace(&old_hash, &new_hash.to_string());
    assert_ne!(forged_params, row.params, "artifact hash must be rewritten");
    ledger
        .tune_put(
            &row.kernel,
            &row.shape_class,
            &row.machine,
            &forged_params,
            new_measured,
        )
        .expect("overwrite row");
}

/// Alter one numeric field inside the measured payload, leaving the receipt
/// prefix, reps, and verdict fields (everything the per-row checks pin)
/// untouched.
fn altered_measured(measured: &str) -> String {
    let (before, after) = measured
        .split_once("\"dispersion\":")
        .expect("dispersion field");
    let end = after.find([',', '}']).expect("field end");
    let forged = format!("{before}\"dispersion\":9.5e-1{}", &after[end..]);
    assert_ne!(forged, measured);
    forged
}

#[test]
fn replacing_one_payload_while_keeping_the_receipt_poisons_the_whole_run() {
    let db = temp_db("splice");
    let run = recorded_run(&db);
    let (kernel_a, version_a) = run.kernels[0].clone();
    let (kernel_b, version_b) = run.kernels[1].clone();
    assert_eq!(probe(&run, &kernel_a, &version_a), Staleness::Fresh);
    assert_eq!(probe(&run, &kernel_b, &version_b), Staleness::Fresh);

    let row = roofline_row(&run.ledger, &kernel_a);
    splice_payload(&run.ledger, &row, &altered_measured(&row.measured));

    // The tampered row itself can no longer prove receipt membership...
    assert_eq!(
        probe(&run, &kernel_a, &version_a),
        Staleness::CorruptEvidence
    );
    // ...and every sibling row of the same finalized run is poisoned too:
    // their receipts recompute over the altered payload.
    assert_eq!(
        probe(&run, &kernel_b, &version_b),
        Staleness::CorruptEvidence
    );
    cleanup_db(&db);
}

#[test]
fn rows_added_beyond_the_manifest_are_corrupt_evidence() {
    let db = temp_db("added");
    let run = recorded_run(&db);
    let (kernel_a, version_a) = run.kernels[0].clone();

    // Forge a whole extra kernel row citing the recorded op and receipt.
    let row = roofline_row(&run.ledger, &kernel_a);
    let ghost = "ghost-kernel";
    let ghost_measured = row.measured.replace(
        &format!("\"kernel\":\"{kernel_a}\""),
        &format!("\"kernel\":\"{ghost}\""),
    );
    assert_ne!(ghost_measured, row.measured);
    let ghost_hash = fs_ledger::hash_bytes(ghost_measured.as_bytes());
    run.ledger
        .put_artifact(
            "roofline-benchmark-result",
            ghost_measured.as_bytes(),
            Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}"),
        )
        .expect("store ghost artifact");
    let op: i64 = row
        .params
        .split_once("\"op\":")
        .and_then(|(_, rest)| rest.split_once(','))
        .and_then(|(digits, _)| digits.parse().ok())
        .expect("op id in params");
    run.ledger
        .link(op, &ghost_hash, EdgeRole::Out)
        .expect("ghost edge");
    let ghost_params = row.params.replace(
        &fs_ledger::hash_bytes(row.measured.as_bytes()).to_string(),
        &ghost_hash.to_string(),
    );
    run.ledger
        .tune_put(
            ghost,
            &row.shape_class,
            &row.machine,
            &ghost_params,
            &ghost_measured,
        )
        .expect("insert ghost row");

    // The forged row is self-consistent under every per-row check, but the
    // op-bound manifest has no entry for it.
    assert_eq!(probe(&run, ghost, &version_a), Staleness::CorruptEvidence);
    // The legitimate rows of the run remain untouched and fresh.
    assert_eq!(probe(&run, &kernel_a, &version_a), Staleness::Fresh);
    cleanup_db(&db);
}

#[test]
fn identical_rerun_history_stays_fresh() {
    let db = temp_db("rerun");
    let first = recorded_run(&db);
    let (kernel_a, version_a) = first.kernels[0].clone();
    assert_eq!(probe(&first, &kernel_a, &version_a), Staleness::Fresh);

    // A second honest run against the same ledger: both runs' rows must keep
    // validating (two ops, two manifests, disjoint shape classes).
    let axes = synthetic_axes(FINGERPRINT);
    let (baseline2, identity2) = trusted_baseline(&axes);
    let policy = AxisBaselinePolicy::new(Some(&baseline2), &identity2, 20_010);
    let mut registry = default_registry(1 << 10);
    let mut results = run_registry(&mut registry, 0, 1, &axes);
    let mut finalized = finalize_registry_tuning(&mut registry, &axes, &axes, policy, &results)
        .expect("finalize second run");
    let op2 = record_run(
        &first.ledger,
        &axes,
        &axes,
        policy,
        &mut finalized,
        &mut results,
    )
    .expect("record second run");
    let rerecorded_at = first
        .ledger
        .op(op2)
        .unwrap()
        .expect("second op")
        .t_end
        .expect("finished second op");

    for (kernel, version) in &first.kernels {
        assert_eq!(
            staleness_at(
                &first.ledger,
                kernel,
                version,
                FINGERPRINT,
                Some(first.baseline.content_hash()),
                rerecorded_at + 1,
            )
            .expect("staleness probe"),
            Staleness::Fresh
        );
    }
    cleanup_db(&db);
}

#[test]
fn pre_manifest_v2_rows_are_retired_as_corrupt() {
    let db = temp_db("v2row");
    let run = recorded_run(&db);
    let (kernel_a, version_a) = run.kernels[0].clone();
    assert_eq!(probe(&run, &kernel_a, &version_a), Staleness::Fresh);

    // Rewrite the row's params to the retired v2 schema tag. A v2 row cannot
    // prove membership in an op-bound manifest, so the explicit migration
    // classifies it as corrupt instead of grandfathering it.
    let row = roofline_row(&run.ledger, &kernel_a);
    let v2_params = row.params.replace(
        "\"schema\":\"fs-roofline-ledger-row-v3\"",
        "\"schema\":\"fs-roofline-ledger-row-v2\"",
    );
    assert_ne!(v2_params, row.params);
    run.ledger
        .tune_put(
            &row.kernel,
            &row.shape_class,
            &row.machine,
            &v2_params,
            &row.measured,
        )
        .expect("downgrade row schema");

    assert_eq!(
        probe(&run, &kernel_a, &version_a),
        Staleness::CorruptEvidence
    );
    cleanup_db(&db);
}
