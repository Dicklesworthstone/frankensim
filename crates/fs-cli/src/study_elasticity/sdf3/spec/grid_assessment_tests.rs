//! Source admission for a bounded, physical two-grid affine experiment.
use super::*;

const ASSESSED: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-grid-assessment.fsim"
));
const FIXED: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-responses.fsim"
));

#[test]
fn assessment_admits_one_physical_enrichment_with_the_complete_original_family() {
    let assessed = parse(ASSESSED).unwrap();
    let fixed = parse(FIXED).unwrap();
    let policy = assessed.affine_grid_assessment.unwrap();
    assert_eq!(policy.max_transfer_terms, 500_000);
    assert_eq!(policy.compliance_tolerance, 0.02);
    assert_eq!(policy.response_tolerance, 0.02);
    assert_eq!(policy.stress_tolerance, 0.02);
    assert_eq!(policy.inventory_tolerance, 0.0001);
    assert_eq!(assessed.level, fixed.level);
    assert_eq!(assessed.max_level, assessed.level + 1);
    assert_eq!(assessed.leaves, 64);
    assert_eq!(assessed.affine_responses, fixed.affine_responses);
    assert_eq!(
        assessed.affine_load_space.unwrap().constraint_policy(),
        fixed.affine_load_space.unwrap().constraint_policy()
    );
    assert!(fixed.affine_grid_assessment.is_none());
    assert!(parse(&FIXED.replace(":maximum-level 1", ":maximum-level 2")).is_err());
}

#[test]
fn assessment_refuses_unfunded_grids_and_invalid_change_policies_before_work() {
    for (old, new) in [
        (":maximum-level 2", ":maximum-level 1"),
        (":maximum-leaves 64", ":maximum-leaves 8"),
        (":memory-bytes 1073741824", ":memory-bytes 134217728"),
        (":max-transfer-terms 500000", ":max-transfer-terms 0"),
        (":max-transfer-terms 500000", ":max-transfer-terms 4000001"),
        (
            ":compliance-change-tolerance 0.02",
            ":compliance-change-tolerance -0.01",
        ),
        (
            ":response-change-tolerance 0.02",
            ":response-change-tolerance 1000001.0",
        ),
        (
            ":stress-change-tolerance 0.02",
            ":stress-change-tolerance -0.01",
        ),
        (
            ":inventory-change-tolerance 0.0001",
            ":inventory-change-tolerance -0.01",
        ),
    ] {
        assert!(ASSESSED.contains(old));
        assert!(parse(&ASSESSED.replace(old, new)).is_err(), "{new}");
    }
    let zero = ASSESSED
        .replace("change-tolerance 0.02", "change-tolerance 0.0")
        .replace("change-tolerance 0.0001", "change-tolerance 0.0");
    assert!(parse(&zero).is_ok());
}

#[test]
fn assessment_section_is_unique_final_and_bound_to_affine_material_design() {
    let start = ASSESSED.find("  (affine-grid-assessment").unwrap();
    let section = ASSESSED[start..].trim_end().strip_suffix(')').unwrap();
    let prefix = &ASSESSED[..start];
    let duplicate = format!("{prefix}{section}\n{section})");
    assert!(parse(&duplicate).is_err());
    let extra = format!("{prefix}{section}\n  (unrecognized-assessment-field :value 1))");
    assert!(parse(&extra).is_err());
    let before_responses = prefix.find("  (affine-responses").unwrap();
    let reordered = format!(
        "{}{}\n{})",
        &prefix[..before_responses],
        section,
        &prefix[before_responses..]
    );
    assert!(parse(&reordered).is_err());
    let mut wrong_owner = parse(ASSESSED).unwrap();
    wrong_owner.affine_load_space = None;
    assert!(
        wrong_owner
            .affine_grid_assessment
            .unwrap()
            .validate(&wrong_owner)
            .is_err()
    );
}
