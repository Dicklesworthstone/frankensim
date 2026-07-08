//! Battery for go-to-market wedge selection (addendum Proposal 7). Verifies
//! the chosen beachhead, the four-criteria scoring, the named second/third
//! verticals with their proposal-exercise mapping, the measurable cycle-time
//! kill criterion, the audit, and the negative doctrine.

use fs_wedge::{
    CHT_BASELINE, STRONG_THRESHOLD, WEDGE_DOCTRINE, WedgeCriterion, audit, chosen_wedge,
    four_criteria, to_json, verticals,
};

#[test]
fn the_beachhead_is_conjugate_heat_transfer() {
    let w = chosen_wedge();
    assert_eq!(w.name, "conjugate-heat-transfer");
    assert_eq!(w.rank, 1);
    // it exercises incremental re-solve (2), adjoints (1), the ladder (3),
    // and the evidence package (12).
    assert!(w.exercises.contains(&"2") && w.exercises.contains(&"3"));
}

#[test]
fn the_chosen_wedge_is_strong_on_all_four_criteria() {
    let w = chosen_wedge();
    for c in four_criteria() {
        assert!(
            w.score(c) >= STRONG_THRESHOLD,
            "{} weak on {}",
            w.name,
            c.label()
        );
    }
    assert!(w.weakest_criterion_score() >= STRONG_THRESHOLD);
}

#[test]
fn three_verticals_are_ranked_with_proposal_mappings() {
    let vs = verticals();
    assert_eq!(vs.len(), 3);
    let mut ranks: Vec<u8> = vs.iter().map(|v| v.rank).collect();
    ranks.sort_unstable();
    assert_eq!(ranks, vec![1, 2, 3]);
    // second vertical exercises Proposal 1; third exercises 11 and 4.
    let aero = vs
        .iter()
        .find(|v| v.name == "aeroelastic-screening")
        .unwrap();
    assert_eq!(aero.rank, 2);
    assert!(aero.exercises.contains(&"1"));
    let am = vs
        .iter()
        .find(|v| v.name == "additive-manufacturing-distortion")
        .unwrap();
    assert_eq!(am.rank, 3);
    assert!(am.exercises.contains(&"11") && am.exercises.contains(&"4"));
    // every vertical names at least one exercised proposal.
    assert!(vs.iter().all(|v| !v.exercises.is_empty()));
}

#[test]
fn the_cycle_time_kill_criterion_is_measurable() {
    assert!((CHT_BASELINE.baseline_days - 5.0).abs() < 1e-12);
    assert!((CHT_BASELINE.target_reduction - 3.0).abs() < 1e-12);
    assert_eq!(CHT_BASELINE.kill_within_quarters, 2);
    // a 1.5-day cycle is a 3.33x reduction -> meets the criterion.
    assert!(CHT_BASELINE.meets_kill_criterion(1.5));
    // a 2-day cycle is only 2.5x -> does not.
    assert!(!CHT_BASELINE.meets_kill_criterion(2.0));
    // guard against divide-by-zero.
    assert!(!CHT_BASELINE.meets_kill_criterion(0.0));
}

#[test]
fn the_audit_is_complete() {
    let a = audit();
    assert!(a.ok(), "gaps: {:?}", a.gaps);
    assert!(a.passed("chosen-strong-on-all"));
    assert!(a.passed("ranks-complete"));
    assert!(a.passed("all-exercise-proposals"));
    assert!(a.passed("kill-criterion-measurable"));
    assert_eq!(a.checks.len(), 4);
}

#[test]
fn the_negative_doctrine_is_stated() {
    // the load-bearing anti-pattern: don't sell against peak single-physics.
    assert!(
        WEDGE_DOCTRINE
            .to_lowercase()
            .contains("peak single-physics")
    );
    // criterion labels are unique.
    let labels: Vec<&str> = WedgeCriterion::ALL.iter().map(|c| c.label()).collect();
    let mut sorted = labels.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), labels.len());
}

#[test]
fn json_is_well_formed_and_deterministic() {
    let j = to_json();
    assert_eq!(j, to_json());
    assert!(j.starts_with('{') && j.ends_with('}'));
    assert!(j.contains("conjugate-heat-transfer"));
    assert!(j.contains("\"target_reduction\":3"));
    assert_eq!(j.matches("\"rank\":").count(), 3);
}
