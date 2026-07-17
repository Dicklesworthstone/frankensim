//! The PV exit criteria as executable tests (milestone-pv):
//! 1. Same study twice → identical artifact hashes (determinism).
//! 2. Replay from the ledger reproduces the run.
//! 3. Corrupting a ledgered artifact makes replay fail LOUDLY.
//! 4. Adjoint gradient matches central differences (checked inside the run —
//!    a failing check aborts the study).
//! 5. Structured, teaching errors on bad studies.

use std::sync::atomic::AtomicU32;

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

const E2E_SUITE: &str = "fs-vskeleton/e2e";
const INFRA_SUITE: &str = "fs-vskeleton";
const STUDY_INPUT_SEED: u64 = 0x5EED_0001;
const BAD_STUDY_INPUT_SEED: u64 = 1;
const FIXED_INPUT_SEED: u64 = 0;

fn temp_db() -> String {
    let n = NEXT_DB.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!("fs-vskeleton-e2e-{}-{n}.db", std::process::id()))
        .display()
        .to_string()
}

fn study() -> String {
    format!(
        r#"(study "pv-plate-hole-v1"
  (seed {STUDY_INPUT_SEED:#X})
  (grid 33)
  (budget (cg-iters 2000000))
  (hole-radius 0.25)
  (opt-steps 3)
  (step-size 0.15)
  (volume-weight 0.05))"#
    )
}

fn verdict(suite: &str, case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(suite, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: suite.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("vertical-skeleton verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line)
        .expect("vertical-skeleton verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

#[test]
fn pv_001_deterministic_rerun_hash_equality() {
    let (db_a, db_b) = (temp_db(), temp_db());
    let study = study();
    let a = fs_vskeleton::run_study(&study, &db_a).expect("run a");
    let b = fs_vskeleton::run_study(&study, &db_b).expect("run b");
    let pass = a.artifact_hashes == b.artifact_hashes
        && a.report == b.report
        && !a.artifact_hashes.is_empty();
    verdict(
        E2E_SUITE,
        "pv-001",
        pass,
        &format!(
            "{} artifacts bit-identical across reruns",
            a.artifact_hashes.len()
        ),
        STUDY_INPUT_SEED,
    );
    let _ = std::fs::remove_file(&db_a);
    let _ = std::fs::remove_file(&db_b);
}

#[test]
fn pv_002_replay_reproduces_ledger() {
    let db = temp_db();
    let study = study();
    let outcome = fs_vskeleton::run_study(&study, &db).expect("run");
    fs_vskeleton::replay(&db).expect("replay must reproduce");
    verdict(
        E2E_SUITE,
        "pv-002",
        true,
        &format!("replay matched {} artifacts", outcome.artifact_hashes.len()),
        STUDY_INPUT_SEED,
    );
    let _ = std::fs::remove_file(&db);
}

#[test]
fn pv_003_corrupted_ledger_fails_loudly() {
    let db = temp_db();
    let study = study();
    fs_vskeleton::run_study(&study, &db).expect("run");
    let led = fs_vskeleton::ledger::MiniLedger::open(&db).expect("open");
    led.corrupt_first_artifact_for_test().expect("corrupt");
    let err = fs_vskeleton::replay(&db).expect_err("tampered ledger must not replay");
    verdict(
        E2E_SUITE,
        "pv-003",
        err.contains("LedgerCorruption"),
        &format!("byte corruption detected and refused: {err}"),
        STUDY_INPUT_SEED,
    );
    let _ = std::fs::remove_file(&db);
}

#[test]
fn pv_004_objective_improves_and_gradient_checks_pass() {
    let db = temp_db();
    let study = study();
    let o = fs_vskeleton::run_study(&study, &db).expect("run");
    let objective_start = o.objective_trace.first().copied().unwrap_or(f64::NAN);
    let objective_end = o.objective_trace.last().copied().unwrap_or(f64::NAN);
    let worst_gradient_error = o.gradient_check_rel_err.iter().copied().fold(0.0, f64::max);
    let pass = o.objective_trace.len() == 3
        && objective_end < objective_start
        && o.gradient_check_rel_err.iter().all(|&e| e < 1e-4)
        && o.report.contains("gradient checks: 3 / 3 passed");
    verdict(
        E2E_SUITE,
        "pv-004",
        pass,
        &format!(
            "J: {:.6e} -> {:.6e}; worst grad rel err {:.2e}",
            objective_start, objective_end, worst_gradient_error
        ),
        STUDY_INPUT_SEED,
    );
    let _ = std::fs::remove_file(&db);
}

#[test]
fn pv_005_bad_studies_teach() {
    let db = temp_db();
    // Missing budget: the P4 message.
    let missing_budget = format!(
        r#"(study "x" (seed {BAD_STUDY_INPUT_SEED}) (grid 33) (hole-radius 0.25)
            (opt-steps 1) (step-size 0.1) (volume-weight 0.05))"#
    );
    let e = fs_vskeleton::run_study(&missing_budget, &db).expect_err("budgets are mandatory");
    let missing_budget_teaches = e.contains("budgets are mandatory");
    // Tiny budget: enforcement, not advice.
    let tiny_budget = format!(
        r#"(study "x" (seed {BAD_STUDY_INPUT_SEED}) (grid 33) (budget (cg-iters 3)) (hole-radius 0.25)
            (opt-steps 1) (step-size 0.1) (volume-weight 0.05))"#
    );
    let e = fs_vskeleton::run_study(&tiny_budget, &db).expect_err("budget must be enforced");
    verdict(
        E2E_SUITE,
        "pv-005",
        missing_budget_teaches && e.contains("BudgetExhausted"),
        "missing and exhausted budgets both refused with guidance",
        BAD_STUDY_INPUT_SEED,
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
    let pass = h.len() == 64
        && h.chars().all(|c| c.is_ascii_hexdigit())
        && h != fs_blake3::hash_bytes(b"payload").to_hex();
    verdict(
        INFRA_SUITE,
        "hash-shape",
        pass,
        &format!("64-hex domain-separated BLAKE3 artifact address: {h}"),
        FIXED_INPUT_SEED,
    );
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
    // Regression: the teaching string was line-wrapped INSIDE the literal, so it
    // rendered with long embedded space runs. It must read cleanly.
    let pass =
        err.contains("LedgerFormatMismatch") && err.contains("fresh ledger") && !err.contains("  ");
    verdict(
        INFRA_SUITE,
        "v1-refusal",
        pass,
        &format!(
            "pre-v2 ledger refused with teaching diagnostic: {}",
            err.split(':').next().unwrap_or("")
        ),
        FIXED_INPUT_SEED,
    );
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
    assert!(
        !err.contains("  "),
        "teaching message must not have garbled space runs: {err}"
    );
    let _ = std::fs::remove_file(&db);
}
