//! Battery for the addendum crate registry (Proposal 7 contract discipline /
//! the crate-contracts governance bead): completeness, unique names, expected
//! inventory, owner mapping, the audit, and deterministic JSON.

use fs_govern::{addendum_crates, crate_audit, crates_json};
use std::collections::BTreeSet;

#[test]
fn registry_lists_the_net_new_addendum_crates() {
    let cs = addendum_crates();
    assert_eq!(cs.len(), 7);
    let names: BTreeSet<&str> = cs.iter().map(|c| c.name).collect();
    assert_eq!(names.len(), 7, "names must be unique");
    let expected: BTreeSet<&str> = [
        "fs-iface",
        "fs-ladder",
        "fs-probe",
        "fs-spececo",
        "fs-verify",
        "fs-recompute",
        "fs-govern",
    ]
    .into_iter()
    .collect();
    assert_eq!(names, expected);
}

#[test]
fn every_crate_declares_purpose_owner_layer_and_no_claim() {
    for c in addendum_crates() {
        assert!(!c.purpose.is_empty(), "{} purpose", c.name);
        assert!(!c.owning_proposal.is_empty(), "{} owning_proposal", c.name);
        assert!(!c.layer.is_empty(), "{} layer", c.name);
        assert!(!c.no_claim.is_empty(), "{} no_claim", c.name);
    }
}

#[test]
fn owners_map_to_the_expected_proposals() {
    let find = |name: &str| addendum_crates().iter().find(|c| c.name == name).unwrap();
    assert_eq!(find("fs-iface").owning_proposal, "13");
    assert_eq!(find("fs-ladder").owning_proposal, "3");
    assert_eq!(find("fs-probe").owning_proposal, "3");
    assert_eq!(find("fs-spececo").owning_proposal, "9");
    assert_eq!(find("fs-verify").owning_proposal, "9");
    assert_eq!(find("fs-recompute").owning_proposal, "2");
    assert_eq!(find("fs-govern").owning_proposal, "governance");
}

#[test]
fn crate_audit_is_complete() {
    let a = crate_audit();
    assert_eq!(a.total, 7);
    assert_eq!(
        a.complete, 7,
        "every crate must declare purpose + owner + no-claim"
    );
    assert!(a.ok(), "no gaps: {:?}", a.gaps);
}

#[test]
fn crates_json_is_well_formed_and_deterministic() {
    let j = crates_json();
    assert!(j.starts_with('[') && j.ends_with(']'));
    assert_eq!(j.matches("\"name\":").count(), 7);
    assert!(j.contains("fs-iface") && j.contains("fs-spececo"));
    assert!(j.contains("\"no_claim\":"));
    assert!(!j.contains(",,"));
    assert_eq!(crates_json(), crates_json());
}
