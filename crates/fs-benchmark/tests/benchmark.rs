//! Battery for the wedge-vertical benchmark corpus (addendum Proposal 7).
//! Verifies the datasets exist and are colored, the measurement helpers
//! compute the addendum kill numbers, the Governance-Rule-2 discharge (every
//! instrumented proposal references a real dataset), the deterministic content
//! digest, and the completeness audit.

use fs_benchmark::{
    BENCHMARK_VERSION, ColorRank, accept_rate, audit, conflict_rate, corpus_digest, design_tasks,
    edit_traces, instrumented_proposals, merge_trials, mms_battery, query_set, rate, speedup,
    win_rate,
};

#[test]
fn every_dataset_is_populated() {
    assert_eq!(query_set().len(), 3);
    assert_eq!(design_tasks().len(), 3);
    assert_eq!(edit_traces().len(), 2);
    assert_eq!(mms_battery().len(), 2);
    assert_eq!(merge_trials().len(), 2);
}

#[test]
fn reference_answers_carry_their_color() {
    // the three color classes are all represented in the query set.
    let ranks: Vec<ColorRank> = query_set().iter().map(|q| q.reference_color).collect();
    assert!(ranks.contains(&ColorRank::Verified));
    assert!(ranks.contains(&ColorRank::Validated));
    assert!(ranks.contains(&ColorRank::Estimated));
    // every query has a finite reference answer + positive cost.
    for q in query_set() {
        assert!(q.reference_answer.is_finite() && q.reference_cost > 0.0);
    }
}

#[test]
fn the_measurement_helpers_compute_the_kill_numbers() {
    // Proposal 8 / 2 speedup: a 1000->400 cost is 2.5x (clears the >=2x bar).
    assert!((speedup(1000.0, 400.0) - 2.5).abs() < 1e-12);
    assert!(speedup(1000.0, 400.0) >= 2.0);
    assert_eq!(speedup(1000.0, 0.0).to_bits(), 0.0f64.to_bits()); // divide-by-zero guard
    // Proposal 1 win rate: 4 of 5 tasks won = 0.8 (clears >=0.70).
    assert!((win_rate(&[true, true, false, true, true]) - 0.8).abs() < 1e-12);
    assert_eq!(win_rate(&[]).to_bits(), 0.0f64.to_bits());
    // Proposal 9 accept-rate.
    assert!((accept_rate(35, 100) - 0.35).abs() < 1e-12);
    assert_eq!(rate(1, 0).to_bits(), 0.0f64.to_bits());
}

#[test]
fn the_harmonic_conflict_rate_is_computed_from_merge_trials() {
    // Proposal 10 kill metric: <25% of merges surface harmonic conflicts.
    let t = &merge_trials()[0]; // 6 of 40
    assert!((conflict_rate(t) - 0.15).abs() < 1e-12);
    assert!(conflict_rate(t) < 0.25);
}

#[test]
fn the_skippable_fraction_comes_from_edit_traces() {
    // Proposal 2: the known-correct skip set (raise-fin-count: 96 of 120).
    let e = &edit_traces()[0];
    assert!((rate(e.correct_skips, e.total_ops) - 0.8).abs() < 1e-12);
}

#[test]
fn governance_rule_2_is_discharged_for_every_instrumented_proposal() {
    let ips = instrumented_proposals();
    assert_eq!(ips.len(), 8);
    // the six+ proposals named in the bead are all present.
    let props: Vec<&str> = ips.iter().map(|p| p.proposal).collect();
    for p in ["8", "1", "2", "D", "F", "9", "10", "A"] {
        assert!(props.contains(&p), "proposal {p} not instrumented");
    }
    // every instrumented proposal declares a kill metric AND a dataset.
    for ip in ips {
        assert!(!ip.kill_metric.is_empty());
        assert!(!ip.dataset.is_empty());
    }
}

#[test]
fn the_corpus_digest_is_deterministic() {
    // bit-stable across calls -> measurements are replayable.
    assert_eq!(corpus_digest(), corpus_digest());
    assert_ne!(corpus_digest(), 0);
}

#[test]
fn the_audit_is_complete() {
    let a = audit();
    assert!(a.ok(), "gaps: {:?}", a.gaps);
    assert_eq!(a.version, BENCHMARK_VERSION);
    assert_eq!(a.instrumented, 8);
}
