//! Constellation smoke test: fs-exec's foundational dependency is asupersync
//! (plan §4.1/§12 contract: structured-concurrency scopes, cancel-correctness,
//! budgets). This test exercises the REAL library — the adapter contract's
//! semantics get exercised as the two-lane executor lands; this proves the
//! dependency wiring, version, and the Budget vocabulary the Cx will carry.

use asupersync::types::Budget;

#[test]
fn asupersync_links_and_budget_vocabulary_holds() {
    // Budget is the P4 primitive fs-exec's Cx will thread through kernels.
    let infinite = Budget::INFINITE;
    let zero = Budget::ZERO;
    let minimal = Budget::MINIMAL;
    assert!(infinite.poll_quota > minimal.poll_quota);
    assert_eq!(zero.poll_quota, 0);
    assert_eq!(
        minimal.poll_quota, 100,
        "cleanup budget contract (bounded cancellation drain)"
    );
    let detail = format!(
        "poll quotas: inf>{} min={} zero={}",
        minimal.poll_quota, minimal.poll_quota, zero.poll_quota
    );
    let mut emitter = fs_obs::Emitter::new("fs-exec/constellation", "asupersync-budget");
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-exec/constellation".to_string(),
            case: "asupersync-budget".to_string(),
            pass: true,
            detail,
            seed: 0,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("constellation verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("constellation verdict must use the fs-obs wire schema");
    println!("{line}");
}
