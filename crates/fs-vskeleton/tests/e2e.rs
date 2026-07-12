//! The PV exit criteria as executable tests (milestone-pv):
//! 1. Same study twice → identical artifact hashes (determinism).
//! 2. Replay from the ledger reproduces the run.
//! 3. Corrupting a ledgered artifact makes replay fail LOUDLY.
//! 4. Adjoint gradient matches central differences (checked inside the run —
//!    a failing check aborts the study).
//! 5. Structured, teaching errors on bad studies.

use std::sync::atomic::AtomicU32;

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

fn temp_db() -> String {
    let n = NEXT_DB.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!("fs-vskeleton-e2e-{}-{n}.db", std::process::id()))
        .display()
        .to_string()
}

const STUDY: &str = r#"(study "pv-plate-hole-v1"
  (seed 0x5EED0001)
  (grid 33)
  (budget (cg-iters 2000000))
  (hole-radius 0.25)
  (opt-steps 3)
  (step-size 0.15)
  (volume-weight 0.05))"#;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-vskeleton/e2e\",\"case\":\"{case}\",\"verdict\":\"pass\",\"detail\":\"{detail}\"}}"
    );
}

#[test]
fn pv_001_deterministic_rerun_hash_equality() {
    let (db_a, db_b) = (temp_db(), temp_db());
    let a = fs_vskeleton::run_study(STUDY, &db_a).expect("run a");
    let b = fs_vskeleton::run_study(STUDY, &db_b).expect("run b");
    assert_eq!(
        a.artifact_hashes, b.artifact_hashes,
        "artifact hashes must be identical"
    );
    assert_eq!(a.report, b.report, "reports must be byte-identical");
    assert!(!a.artifact_hashes.is_empty());
    verdict(
        "pv-001",
        &format!(
            "{} artifacts bit-identical across reruns",
            a.artifact_hashes.len()
        ),
    );
    let _ = std::fs::remove_file(&db_a);
    let _ = std::fs::remove_file(&db_b);
}

#[test]
fn pv_002_replay_reproduces_ledger() {
    let db = temp_db();
    let outcome = fs_vskeleton::run_study(STUDY, &db).expect("run");
    fs_vskeleton::replay(&db).expect("replay must reproduce");
    verdict(
        "pv-002",
        &format!("replay matched {} artifacts", outcome.artifact_hashes.len()),
    );
    let _ = std::fs::remove_file(&db);
}

#[test]
fn pv_003_corrupted_ledger_fails_loudly() {
    let db = temp_db();
    fs_vskeleton::run_study(STUDY, &db).expect("run");
    let led = fs_vskeleton::ledger::MiniLedger::open(&db).expect("open");
    led.corrupt_first_artifact_for_test().expect("corrupt");
    let err = fs_vskeleton::replay(&db).expect_err("tampered ledger must not replay");
    assert!(
        err.contains("LedgerCorruption"),
        "loud corruption verdict expected: {err}"
    );
    verdict("pv-003", "byte corruption detected and refused");
    let _ = std::fs::remove_file(&db);
}

#[test]
fn pv_004_objective_improves_and_gradient_checks_pass() {
    let db = temp_db();
    let o = fs_vskeleton::run_study(STUDY, &db).expect("run");
    assert_eq!(o.objective_trace.len(), 3);
    assert!(
        o.objective_trace.last().unwrap() < o.objective_trace.first().unwrap(),
        "projected GD must reduce the objective: {:?}",
        o.objective_trace
    );
    assert!(o.gradient_check_rel_err.iter().all(|&e| e < 1e-4));
    assert!(o.report.contains("gradient checks: 3 / 3 passed"));
    verdict(
        "pv-004",
        &format!(
            "J: {:.6e} -> {:.6e}; worst grad rel err {:.2e}",
            o.objective_trace[0],
            o.objective_trace[2],
            o.gradient_check_rel_err.iter().copied().fold(0.0, f64::max)
        ),
    );
    let _ = std::fs::remove_file(&db);
}

#[test]
fn pv_005_bad_studies_teach() {
    let db = temp_db();
    // Missing budget: the P4 message.
    let e = fs_vskeleton::run_study(
        r#"(study "x" (seed 1) (grid 33) (hole-radius 0.25)
            (opt-steps 1) (step-size 0.1) (volume-weight 0.05))"#,
        &db,
    )
    .expect_err("budgets are mandatory");
    assert!(e.contains("budgets are mandatory"), "{e}");
    // Tiny budget: enforcement, not advice.
    let e = fs_vskeleton::run_study(
        r#"(study "x" (seed 1) (grid 33) (budget (cg-iters 3)) (hole-radius 0.25)
            (opt-steps 1) (step-size 0.1) (volume-weight 0.05))"#,
        &db,
    )
    .expect_err("budget must be enforced");
    assert!(e.contains("BudgetExhausted"), "{e}");
    verdict(
        "pv-005",
        "missing and exhausted budgets both refused with guidance",
    );
    let _ = std::fs::remove_file(&db);
}

#[test]
fn blake3_content_addresses_are_64_hex_and_domain_separated() {
    // Bead frankensim-ynsl: the FNV placeholder (16 hex) is retired; the
    // v2 format uses domain-separated BLAKE3 (64 hex). The domain string
    // matters: the same bytes under a different domain must not collide
    // with artifact addresses.
    let h = fs_vskeleton::ledger::content_hash(b"payload");
    assert_eq!(h.len(), 64, "BLAKE3 hex width: {h}");
    assert!(h.chars().all(|c| c.is_ascii_hexdigit()));
    assert_ne!(
        h,
        fs_blake3::hash_bytes(b"payload").to_hex(),
        "artifact addresses must be domain-separated from plain hashing"
    );
    println!("{{\"suite\":\"fs-vskeleton\",\"case\":\"hash-shape\",\"verdict\":\"pass\",\"hash\":\"{h}\"}}");
}

#[test]
fn pre_v2_ledger_is_version_refused_with_teaching_error() {
    // Bead frankensim-ynsl: a ledger holding artifacts but no format
    // version is FNV-era data; opening it under v2 must refuse with the
    // migration named, never silently misread 16-hex addresses as v2.
    let db = temp_db();
    {
        // Forge a v1-shaped ledger: schema without the meta stamp, one
        // FNV-style artifact row.
        let raw = fsqlite::Connection::open(&db).expect("raw open");
        raw.execute("CREATE TABLE artifacts(hash TEXT PRIMARY KEY, kind TEXT, bytes BLOB)")
            .expect("v1 ddl");
        raw.prepare("INSERT INTO artifacts(hash, kind, bytes) VALUES (?1, ?2, ?3)")
            .expect("prepare")
            .execute_with_params(&[
                fsqlite::SqliteValue::Text("00000000cbf29ce4".into()),
                fsqlite::SqliteValue::Text("field".into()),
                fsqlite::SqliteValue::Blob(b"legacy".to_vec().into()),
            ])
            .expect("v1 row");
    }
    let err = fs_vskeleton::ledger::MiniLedger::open(&db)
        .err()
        .expect("pre-v2 ledger must refuse");
    assert!(
        err.contains("LedgerFormatMismatch") && err.contains("fresh ledger"),
        "teaching refusal expected, got: {err}"
    );
    println!("{{\"suite\":\"fs-vskeleton\",\"case\":\"v1-refusal\",\"verdict\":\"pass\",\"detail\":\"{}\"}}",
        err.split(':').next().unwrap_or(""));
    let _ = std::fs::remove_file(&db); // test temp file cleanup, same as temp_db siblings
}

#[test]
fn future_format_ledger_is_version_refused() {
    let db = temp_db();
    {
        let l = fs_vskeleton::ledger::MiniLedger::open(&db).expect("fresh v2");
        l.put_artifact("field", b"bytes").expect("put");
    }
    {
        let raw = fsqlite::Connection::open(&db).expect("raw");
        raw.execute("UPDATE vskeleton_meta SET value = '99' WHERE key = 'format_version'")
            .expect("forge future version");
    }
    let err = fs_vskeleton::ledger::MiniLedger::open(&db)
        .err()
        .expect("future format must refuse");
    assert!(err.contains("v99"), "names the found version: {err}");
    let _ = std::fs::remove_file(&db);
}
