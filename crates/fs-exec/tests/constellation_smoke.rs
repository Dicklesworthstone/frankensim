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
    println!(
        "{{\"suite\":\"fs-exec/constellation\",\"case\":\"asupersync-budget\",\"verdict\":\"pass\",\"detail\":\"poll quotas: inf>{} min={} zero={}\"}}",
        minimal.poll_quota, minimal.poll_quota, zero.poll_quota
    );
}
