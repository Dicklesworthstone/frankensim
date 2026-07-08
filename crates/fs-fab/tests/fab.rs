//! Battery for fabrication & code constraints (fs-fab). Covers the constraint
//! families (detection + margins), AISC/ACI-class rules, repair suggestions,
//! differentiable gradient gates, the manufacturability report (localization +
//! fantasy-artifact prevention), and the Evidence-typed cost/carbon model.

use fs_fab::{
    Color, ConstraintKind, Differentiability, bolt_spacing_aisc, check_all, cnc_tool_radius,
    draft_angle, embodied_carbon, evaluate, member_length_transport, min_feature_size,
    overhang_angle, process_cost, rebar_spacing_aci,
};

#[test]
fn constraint_families_detect_violations_and_margins() {
    // additive overhang (<= 45deg): 30 ok, 60 violated.
    let oh = overhang_angle(45.0);
    assert!(oh.satisfied(30.0) && !oh.satisfied(60.0));
    assert!((oh.margin(30.0) - 15.0).abs() < 1e-12);
    // additive min feature (>= 0.5mm): 0.3 violated, 1.0 ok.
    assert!(!min_feature_size(0.5).satisfied(0.3));
    // subtractive tool reachability (concave radius >= tool 3mm).
    assert!(!cnc_tool_radius(3.0).satisfied(2.0));
    // casting draft (>= 3deg).
    assert!(!draft_angle(3.0).satisfied(1.0));
}

#[test]
fn aisc_and_aci_rules_encode_the_published_limits() {
    // AISC-360 §J3.3: min bolt spacing = 2⅔ · diameter. For a 19mm bolt -> 50.67mm.
    let bolt = bolt_spacing_aisc(19.0);
    assert!((bolt.limit - (8.0 / 3.0) * 19.0).abs() < 1e-12);
    assert_eq!(bolt.kind, ConstraintKind::Code("AISC-360"));
    assert!(!bolt.satisfied(45.0) && bolt.satisfied(60.0));
    // ACI-318 §25.2.1: min clear spacing = max(25mm, bar diameter).
    let rebar = rebar_spacing_aci(16.0);
    assert!((rebar.limit - 25.0).abs() < 1e-12);
    assert!((rebar_spacing_aci(32.0).limit - 32.0).abs() < 1e-12);
    assert_eq!(rebar.kind, ConstraintKind::Code("ACI-318"));
}

#[test]
fn a_violation_carries_a_repair_suggestion() {
    let oh = overhang_angle(45.0);
    // satisfied -> no repair.
    assert!(oh.repair(30.0).is_none());
    // violated -> move the feature to the limit.
    let r = oh.repair(60.0).unwrap();
    assert!((r.target_value - 45.0).abs() < 1e-12);
    assert!((r.delta + 15.0).abs() < 1e-12); // -15deg
    // an AtLeast violation repairs upward.
    let mf = min_feature_size(0.5).repair(0.3).unwrap();
    assert!((mf.delta - 0.2).abs() < 1e-12);
}

#[test]
fn differentiable_constraints_pass_a_gradient_gate() {
    let e = 1e-6;
    for c in [
        overhang_angle(45.0),
        cnc_tool_radius(3.0),
        draft_angle(3.0),
        bolt_spacing_aisc(19.0),
    ] {
        let g = c.margin_gradient().expect("differentiable");
        let fd = (c.margin(10.0 + e) - c.margin(10.0 - e)) / (2.0 * e);
        assert!((fd - g).abs() < 1e-6, "{} gradient", c.name);
    }
    // AtMost gradient is -1, AtLeast is +1.
    assert_eq!(overhang_angle(45.0).margin_gradient(), Some(-1.0));
    assert_eq!(min_feature_size(0.5).margin_gradient(), Some(1.0));
    // a DISCRETE constraint has no gradient.
    let ml = member_length_transport(12.0);
    assert_eq!(ml.differentiability, Differentiability::Discrete);
    assert_eq!(ml.margin_gradient(), None);
}

#[test]
fn the_report_localizes_violations_and_prevents_fantasy_artifacts() {
    let design = vec![
        (overhang_angle(45.0), 30.0),          // ok
        (min_feature_size(0.5), 0.3),          // VIOLATED
        (draft_angle(3.0), 5.0),               // ok
        (member_length_transport(12.0), 15.0), // VIOLATED
    ];
    let report = check_all(&design);
    // a design with any violation is NOT manufacturable (no fantasy artifact).
    assert!(!report.feasible);
    // the violations are localized by name.
    assert!(report.violations.contains(&"additive-min-feature"));
    assert!(report.violations.contains(&"member-transport-length"));
    assert_eq!(report.violations.len(), 2);
    // a fully-feasible design passes.
    let ok = check_all(&[(overhang_angle(45.0), 20.0), (draft_angle(3.0), 5.0)]);
    assert!(ok.feasible && ok.violations.is_empty());
    // one result carries its repair.
    assert!(evaluate(&min_feature_size(0.5), 0.3).repair.is_some());
}

#[test]
fn cost_and_carbon_are_evidence_typed_estimates() {
    // 100 units at rate 2 with 10% relative uncertainty.
    let cost = process_cost(100.0, 2.0, 0.1);
    assert!((cost.mean - 200.0).abs() < 1e-12);
    assert!((cost.std() - 20.0).abs() < 1e-12);
    // a modeled quantity is estimated-color.
    assert!(matches!(cost.color(), Color::Estimated { .. }));
    // embodied carbon likewise.
    let carbon = embodied_carbon(500.0, 1.9, 0.15);
    assert!((carbon.mean - 950.0).abs() < 1e-9);
    assert!(carbon.std() > 0.0);
}

#[test]
fn constraint_evaluation_is_deterministic() {
    let c = bolt_spacing_aisc(19.0);
    assert_eq!(evaluate(&c, 55.0), evaluate(&c, 55.0));
}
