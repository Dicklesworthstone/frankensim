//! Battery for the addendum doctrine (P1–P8 + 4 governance rules) and the
//! nineteen-proposal governance registry + audit.

use fs_govern::{governance_audit, principles, proposals, proposals_json, rules};
use std::collections::BTreeSet;

#[test]
fn eight_principles_p1_through_p8() {
    let ps = principles();
    assert_eq!(ps.len(), 8);
    let ids: Vec<&str> = ps.iter().map(|p| p.id).collect();
    assert_eq!(ids, vec!["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]);
    for p in ps {
        assert!(!p.name.is_empty() && !p.statement.is_empty(), "{}", p.id);
    }
}

#[test]
fn four_governance_rules_numbered_1_through_4() {
    let rs = rules();
    assert_eq!(rs.len(), 4);
    let nums: Vec<u8> = rs.iter().map(|r| r.number).collect();
    assert_eq!(nums, vec![1, 2, 3, 4]);
    // Rule 2 is the kill-criteria-enforcement rule.
    assert!(
        rs[1]
            .statement
            .contains("unmeasured survival is not survival")
    );
    for r in rs {
        assert!(!r.name.is_empty() && !r.statement.is_empty());
    }
}

#[test]
fn all_nineteen_proposals_present_with_unique_ids() {
    let ps = proposals();
    assert_eq!(ps.len(), 19);
    let ids: BTreeSet<&str> = ps.iter().map(|p| p.id).collect();
    assert_eq!(ids.len(), 19, "ids must be unique");
    let expected: BTreeSet<&str> = [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "A", "B", "C", "D",
        "E", "F",
    ]
    .into_iter()
    .collect();
    assert_eq!(ids, expected);
}

#[test]
fn proposals_are_in_descending_composite_order() {
    let ps = proposals();
    // top is Certified speculation (850), bottom is the spacetime complex (590).
    assert_eq!(ps.first().unwrap().id, "9");
    assert_eq!(ps.first().unwrap().mean, 850);
    assert_eq!(ps.last().unwrap().id, "4");
    // non-increasing means.
    for w in ps.windows(2) {
        assert!(w[0].mean >= w[1].mean, "{} !>= {}", w[0].mean, w[1].mean);
    }
    // means are in range.
    for p in ps {
        assert!(p.mean <= 1000);
    }
}

#[test]
fn every_proposal_declares_a_kill_metric_and_owner() {
    for p in proposals() {
        assert!(!p.kill_metric.is_empty(), "proposal {} kill_metric", p.id);
        assert!(
            p.owning_bead.starts_with("frankensim-"),
            "proposal {} owner should be a bead id, got {}",
            p.id,
            p.owning_bead
        );
        assert!(!p.phase.is_empty());
    }
}

#[test]
fn governance_audit_is_complete_with_honest_zero_instrumented() {
    let a = governance_audit();
    assert_eq!(a.total, 19);
    assert_eq!(
        a.with_kill_metric_and_owner, 19,
        "every proposal must declare a kill metric + owner"
    );
    assert!(a.ok(), "no gaps: {:?}", a.gaps);
    // honest baseline: no kill metric is instrumented yet.
    assert_eq!(a.instrumented, 0);
}

#[test]
fn owners_map_to_the_expected_beads() {
    let find = |id: &str| proposals().iter().find(|p| p.id == id).unwrap();
    // the Goodhart guard (D) and interface types (13) are the ones this
    // session's author implemented — check their owning beads.
    assert_eq!(find("D").owning_bead, "frankensim-epic-epistype-qmao.5");
    assert_eq!(find("13").owning_bead, "frankensim-epic-selfknow-knh1.1");
    assert_eq!(find("3").owning_bead, "frankensim-epic-epistype-qmao.1");
}

#[test]
fn proposals_json_is_well_formed_and_deterministic() {
    let j = proposals_json();
    assert!(j.starts_with('[') && j.ends_with(']'));
    assert_eq!(j.matches("\"id\":").count(), 19);
    assert!(j.contains("\"mean\":850"));
    assert!(j.contains("frankensim-epic-epistype-qmao.5")); // Goodhart guard owner
    assert!(!j.contains(",,"));
    assert_eq!(proposals_json(), proposals_json());
}
