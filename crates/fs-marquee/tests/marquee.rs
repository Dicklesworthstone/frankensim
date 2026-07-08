//! fs-marquee contract conformance.
//!
//! The current crate is intentionally an L6 admission/status shell. The
//! `marquee` feature may name the future frontier lane, but it must not
//! expose an unproven runner or make simulation, rendering, ledger, or
//! filesystem side effects part of this crate's contract.

use fs_marquee::{MarqueeStatus, VERSION, scope_summary, status};

fn expected_status() -> MarqueeStatus {
    if cfg!(feature = "marquee") {
        MarqueeStatus::SmokeRunnerAvailable
    } else {
        MarqueeStatus::Disabled
    }
}

#[test]
fn marquee_status_matches_feature_gate() {
    assert_eq!(status(), expected_status());
    assert!(!VERSION.is_empty());
}

#[test]
fn marquee_scope_keeps_no_runner_boundary_explicit() {
    let summary = scope_summary();
    assert!(summary.contains("raw SDF"));
    assert!(summary.contains("CutFEM"));
    // The smoke runner shipped; the nightly golden lane is the
    // remaining no-claim boundary.
    assert!(summary.contains("nightly golden pending"));
}
