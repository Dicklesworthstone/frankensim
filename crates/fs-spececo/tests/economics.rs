//! Battery for the certified-speculation economics (addendum Proposal 9).
//! Covers the fail-safe decision rule, telemetry (accept rate, no
//! divide-by-zero, negative warm-starts are not wins), and drift demotion with
//! a minimum sample count + hysteresis (no flapping).

use fs_spececo::{Decision, DriftDetector, ProposerTelemetry, SolveRecord, decide};

#[test]
fn decide_accepts_only_a_finite_nonneg_bound_within_tolerance() {
    // bound meets tolerance -> accept outright (boundary: equal is accepted).
    assert_eq!(decide(4.0, 5.0), Decision::AcceptOutright);
    assert_eq!(decide(5.0, 5.0), Decision::AcceptOutright);
    // bound exceeds tolerance -> warm-start.
    assert_eq!(decide(5.0001, 5.0), Decision::WarmStart);
}

#[test]
fn decide_fails_safe_on_garbage_bounds() {
    // a non-finite bound never accepts (a bad proposer can only waste a check).
    assert_eq!(decide(f64::NAN, 10.0), Decision::WarmStart);
    assert_eq!(decide(f64::INFINITY, f64::INFINITY), Decision::WarmStart);
    // a negative "bound" is nonsense -> warm-start, never accept.
    assert_eq!(decide(-1.0, 10.0), Decision::WarmStart);
    // a non-finite tolerance never accepts.
    assert_eq!(decide(1.0, f64::NAN), Decision::WarmStart);
}

#[test]
fn telemetry_accumulates_accepts_and_savings() {
    let mut t = ProposerTelemetry::new();
    t.record(&SolveRecord::new("neighbor", "Re-2e5", true, 0.1, 0)); // accepted outright
    t.record(&SolveRecord::new("neighbor", "Re-2e5", false, 0.9, 50)); // warm-start saved 50
    t.record(&SolveRecord::new("neighbor", "Re-2e5", false, 0.9, -10)); // warm-start HURT
    let s = t.stats("neighbor", "Re-2e5").unwrap();
    assert_eq!(s.attempts, 3);
    assert_eq!(s.accepts, 1);
    // a negative warm-start is NOT a positive save...
    assert_eq!(s.positive_saves, 1);
    // ...but it does lower the net total (0 + 50 - 10 = 40).
    assert_eq!(s.net_iterations_saved, 40);
    assert!((t.accept_rate("neighbor", "Re-2e5") - 1.0 / 3.0).abs() < 1e-12);
    assert!((t.mean_iterations_saved("neighbor", "Re-2e5") - 40.0 / 3.0).abs() < 1e-12);
}

#[test]
fn accept_rate_of_unknown_pair_is_zero_not_a_panic() {
    let t = ProposerTelemetry::new();
    assert_eq!(
        t.accept_rate("nobody", "nowhere").to_bits(),
        0.0f64.to_bits()
    );
    assert_eq!(
        t.mean_iterations_saved("nobody", "nowhere").to_bits(),
        0.0f64.to_bits()
    );
    assert!(t.stats("nobody", "nowhere").is_none());
}

fn record_n(t: &mut ProposerTelemetry, proposer: &str, regime: &str, n: u64, accepted: bool) {
    for _ in 0..n {
        t.record(&SolveRecord::new(proposer, regime, accepted, 0.5, 0));
    }
}

#[test]
fn drift_does_not_demote_below_the_minimum_sample_count() {
    // a handful of rejects must NOT demote — noise, not signal.
    let mut t = ProposerTelemetry::new();
    record_n(&mut t, "surrogate", "Re-2e5", 3, false); // 3 rejects, all
    let mut d = DriftDetector::new(10, 0.2, 0.5);
    assert!(
        !d.update(&t, "surrogate", "Re-2e5"),
        "3 samples is below min"
    );
    assert!(!d.is_demoted("surrogate", "Re-2e5"));
}

#[test]
fn drift_demotes_on_collapse_and_holds_via_hysteresis() {
    let mut t = ProposerTelemetry::new();
    let mut d = DriftDetector::new(10, 0.2, 0.5);
    // 10 rejects -> accept-rate 0.0 <= 0.2, past min samples -> demote.
    record_n(&mut t, "surrogate", "Re-2e5", 10, false);
    assert!(d.update(&t, "surrogate", "Re-2e5"));
    assert!(d.is_demoted("surrogate", "Re-2e5"));

    // partial recovery into the hysteresis band (rate ~0.286, in (0.2, 0.5]):
    // once demoted, this must NOT re-promote (no flapping).
    record_n(&mut t, "surrogate", "Re-2e5", 4, true); // 14 attempts, 4 accepts
    let rate = t.accept_rate("surrogate", "Re-2e5");
    assert!(
        rate > 0.2 && rate <= 0.5,
        "rate {rate} should be in the hysteresis band"
    );
    assert!(
        d.update(&t, "surrogate", "Re-2e5"),
        "must stay demoted in the band"
    );

    // full recovery above the upper threshold -> restore.
    record_n(&mut t, "surrogate", "Re-2e5", 12, true); // 26 attempts, 16 accepts -> ~0.615
    assert!(t.accept_rate("surrogate", "Re-2e5") > 0.5);
    assert!(
        !d.update(&t, "surrogate", "Re-2e5"),
        "must restore above the upper threshold"
    );
    assert!(!d.is_demoted("surrogate", "Re-2e5"));
}

#[test]
#[should_panic(expected = "hysteresis band inverted")]
fn drift_new_rejects_an_inverted_band() {
    let _ = DriftDetector::new(10, 0.5, 0.2);
}

#[test]
#[should_panic(expected = "hysteresis band inverted")]
fn drift_new_rejects_a_zero_width_band() {
    // Regression: `new` used `restore_above >= demote_below`, admitting an
    // equal-threshold band. That band has NO dead-zone, so a rate at the shared
    // threshold T flaps: undemoted & rate<=T -> demote, then demoted & rate>T ->
    // restore, toggling every update — exactly the flapping hysteresis is meant
    // to prevent. The band must be STRICT (`demote_below < restore_above`).
    let _ = DriftDetector::new(10, 0.5, 0.5);
}

#[test]
fn telemetry_and_drift_are_deterministic() {
    let build = || {
        let mut t = ProposerTelemetry::new();
        record_n(&mut t, "p", "r", 7, false);
        record_n(&mut t, "p", "r", 3, true);
        t
    };
    let a = build();
    let b = build();
    assert_eq!(a.stats("p", "r"), b.stats("p", "r"));
    let mut da = DriftDetector::new(5, 0.2, 0.5);
    let mut db = DriftDetector::new(5, 0.2, 0.5);
    assert_eq!(
        da.update(&a, "p", "r"),
        db.update(&b, "p", "r"),
        "drift decisions are deterministic"
    );
}
