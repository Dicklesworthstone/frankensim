//! Conformance + admission battery for the declarative query language v0
//! (addendum Proposal 8). Covers: the wedge-vertical QoI menu is
//! expressible; QoI metadata is correct; well-posed queries admit; the five
//! ill-posed classes each reject with a distinct teaching finding;
//! dimensional consistency is enforced; admission is deterministic; and the
//! `(query …)` IR surface round-trips under `same_shape`.

use fs_ir::Node;
use fs_ir::admission::Severity;
use fs_ir::query::{FieldRegistry, Qoi, Query, Target};
use fs_qty::Dims;

// von Mises stress: Pa = kg·m⁻¹·s⁻²  → [M,KG,S,K,A]
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
// temperature: K
const KELVIN: Dims = Dims([0, 0, 0, 1, 0, 0]);
// seconds (time) — a deliberately WRONG dimension for a stress tolerance.
const SECONDS: Dims = Dims([0, 0, 1, 0, 0, 0]);
// volume: m³
const VOLUME: Dims = Dims([3, 0, 0, 0, 0, 0]);

fn design() -> FieldRegistry {
    FieldRegistry::new()
        .with("vm", PA)
        .with("temperature", KELVIN)
}

fn ok_max_query() -> Query {
    // "is max von Mises stress under 180 MPa, to within 5 MPa, for $50 in 30s?"
    Query::new(
        Qoi::MaxOverRegion {
            field: "vm".to_string(),
            region: "bracket".to_string(),
        },
        Target::Tolerance {
            value: 5.0e6,
            dims: PA,
        },
        50.0,
        30.0,
    )
}

fn rejects_for<'a>(adm: &'a fs_ir::QueryAdmission, check: &str) -> Vec<&'a str> {
    adm.findings
        .iter()
        .filter(|f| f.severity == Severity::Reject && f.check == check)
        .map(|f| f.check)
        .collect()
}

#[test]
fn wedge_menu_is_expressible_with_correct_metadata() {
    // max-over-region: nonlinear, adjoint-estimable, ladder applies.
    let m = Qoi::MaxOverRegion {
        field: "vm".into(),
        region: "r".into(),
    };
    assert_eq!(m.kind_name(), "max-over-region");
    assert_eq!(m.field(), "vm");
    let meta = m.meta();
    assert!(!meta.linear && meta.adjoint_available && meta.ladder_applicable);
    assert!(!m.is_probabilistic());

    // integral: LINEAR in the field (the DWR sweet spot).
    let i = Qoi::Integral {
        field: "temperature".into(),
        region: "r".into(),
    };
    assert!(i.meta().linear, "a spatial integral must be linear");

    // exceedance: probabilistic, needs an environment.
    let e = Qoi::Exceedance {
        field: "vm".into(),
        region: "r".into(),
        threshold: 180.0e6,
        threshold_dims: PA,
        environment: "site-hazard".into(),
    };
    assert!(e.is_probabilistic() && !e.meta().adjoint_available);
}

#[test]
fn value_dims_follow_the_functional() {
    // max inherits the field's dims;
    let m = Qoi::MaxOverRegion {
        field: "vm".into(),
        region: "r".into(),
    };
    assert_eq!(m.value_dims(PA), PA);
    // integral multiplies by volume;
    let i = Qoi::Integral {
        field: "temperature".into(),
        region: "r".into(),
    };
    assert_eq!(i.value_dims(KELVIN), KELVIN.plus(VOLUME));
    // exceedance is a dimensionless probability.
    let e = Qoi::Exceedance {
        field: "vm".into(),
        region: "r".into(),
        threshold: 1.0,
        threshold_dims: PA,
        environment: "h".into(),
    };
    assert_eq!(e.value_dims(PA), Dims::NONE);
}

#[test]
fn well_posed_query_admits() {
    let adm = ok_max_query().admit(&design());
    assert!(
        adm.admitted,
        "well-posed query must admit: {}",
        adm.diagnosis()
    );
    assert_eq!(adm.value_dims, Some(PA));
    assert!(
        adm.findings.is_empty(),
        "no findings expected: {:?}",
        adm.findings
    );
}

#[test]
fn exceedance_with_matching_threshold_admits() {
    let q = Query::new(
        Qoi::Exceedance {
            field: "vm".into(),
            region: "bracket".into(),
            threshold: 180.0e6,
            threshold_dims: PA,
            environment: "site-hazard".into(),
        },
        Target::Confidence(0.95),
        50.0,
        30.0,
    );
    assert!(q.admit(&design()).admitted);
}

// ---- the five ill-posed classes (test-hardening round 3) ------------------

#[test]
fn reject_zero_budget() {
    let mut q = ok_max_query();
    q.budget_usd = 0.0;
    let adm = q.admit(&design());
    assert!(!adm.admitted);
    assert_eq!(rejects_for(&adm, "query.budget").len(), 1);
    // teaching: the reject carries at least one ranked fix.
    assert!(adm.findings.iter().all(|f| !f.fixes.is_empty()));
}

#[test]
fn reject_past_deadline() {
    let mut q = ok_max_query();
    q.deadline_s = 0.0;
    let adm = q.admit(&design());
    assert!(!adm.admitted);
    assert_eq!(rejects_for(&adm, "query.deadline").len(), 1);
    // negative is also rejected.
    q.deadline_s = -3.0;
    assert!(!q.admit(&design()).admitted);
}

#[test]
fn reject_hundred_percent_confidence() {
    let q = Query::new(
        Qoi::MaxOverRegion {
            field: "vm".into(),
            region: "bracket".into(),
        },
        Target::Confidence(1.0),
        50.0,
        30.0,
    );
    let adm = q.admit(&design());
    assert!(!adm.admitted, "100% confidence is uncertifiable");
    assert_eq!(rejects_for(&adm, "query.confidence").len(), 1);
    // 0% is rejected too; 0.999 is admitted.
    let zero = Query::new(q.qoi.clone(), Target::Confidence(0.0), 50.0, 30.0);
    assert!(!zero.admit(&design()).admitted);
    let near = Query::new(q.qoi.clone(), Target::Confidence(0.999), 50.0, 30.0);
    assert!(near.admit(&design()).admitted);
}

#[test]
fn reject_field_absent_from_design() {
    let q = Query::new(
        Qoi::MaxOverRegion {
            field: "displacement".into(), // not in the design
            region: "bracket".into(),
        },
        Target::Tolerance {
            value: 1.0,
            dims: PA,
        },
        50.0,
        30.0,
    );
    let adm = q.admit(&design());
    assert!(!adm.admitted);
    assert_eq!(rejects_for(&adm, "query.field").len(), 1);
    // the teaching fix lists the fields that DO exist.
    let fix = &adm.findings[0].fixes[0].action;
    assert!(
        fix.contains("vm") && fix.contains("temperature"),
        "fix teaches available fields: {fix}"
    );
}

#[test]
fn reject_self_contradictory_dimensions() {
    // asking for the max von Mises STRESS to within "5 seconds" is nonsense.
    let q = Query::new(
        Qoi::MaxOverRegion {
            field: "vm".into(),
            region: "bracket".into(),
        },
        Target::Tolerance {
            value: 5.0,
            dims: SECONDS,
        },
        50.0,
        30.0,
    );
    let adm = q.admit(&design());
    assert!(!adm.admitted);
    assert_eq!(rejects_for(&adm, "query.dimensions").len(), 1);
}

#[test]
fn reject_exceedance_threshold_off_dimension() {
    let q = Query::new(
        Qoi::Exceedance {
            field: "vm".into(),
            region: "bracket".into(),
            threshold: 10.0,
            threshold_dims: KELVIN, // threshold is not a stress
            environment: "h".into(),
        },
        Target::Confidence(0.95),
        50.0,
        30.0,
    );
    assert!(!q.admit(&design()).admitted);
}

#[test]
fn integral_tolerance_must_carry_volume_dims() {
    // ∫ temperature dV has dims K·m³; a tolerance in plain K is wrong.
    let bad = Query::new(
        Qoi::Integral {
            field: "temperature".into(),
            region: "chip".into(),
        },
        Target::Tolerance {
            value: 1.0,
            dims: KELVIN,
        },
        50.0,
        30.0,
    );
    assert!(
        !bad.admit(&design()).admitted,
        "K tolerance on a K·m³ integral is wrong"
    );
    // the correct tolerance carries volume dims.
    let good = Query::new(
        Qoi::Integral {
            field: "temperature".into(),
            region: "chip".into(),
        },
        Target::Tolerance {
            value: 1.0,
            dims: KELVIN.plus(VOLUME),
        },
        50.0,
        30.0,
    );
    assert!(good.admit(&design()).admitted);
}

#[test]
fn multiple_faults_all_reported_together() {
    // zero budget AND bad dims AND 100% confidence in one query.
    let q = Query::new(
        Qoi::MaxOverRegion {
            field: "vm".into(),
            region: "bracket".into(),
        },
        Target::ToleranceAndConfidence {
            value: 1.0,
            dims: SECONDS,
            confidence: 1.0,
        },
        0.0,
        30.0,
    );
    let adm = q.admit(&design());
    assert!(!adm.admitted);
    assert!(!rejects_for(&adm, "query.budget").is_empty());
    assert!(!rejects_for(&adm, "query.confidence").is_empty());
    assert!(!rejects_for(&adm, "query.dimensions").is_empty());
}

#[test]
fn admission_is_deterministic() {
    let q = ok_max_query();
    let d = design();
    // identical inputs → identical verdict (findings included), replayable.
    assert_eq!(q.admit(&d), q.admit(&d));
    let bad = {
        let mut b = q.clone();
        b.budget_usd = 0.0;
        b
    };
    assert_eq!(bad.admit(&d), bad.admit(&d));
}

// ---- IR surface round-trip ------------------------------------------------

#[test]
fn tolerance_query_round_trips_through_ir() {
    let q = ok_max_query();
    let node = q.to_node();
    let back = Query::from_node(&node).expect("query parses");
    assert_eq!(back.qoi, q.qoi);
    assert_eq!(back.target, q.target);
    // exact round-trip of the same literals: compare by bits (crate idiom).
    assert_eq!(back.budget_usd.to_bits(), q.budget_usd.to_bits());
    assert_eq!(back.deadline_s.to_bits(), q.deadline_s.to_bits());
    // and the emitted node equals a re-emission (same_shape isomorphism).
    assert!(node.same_shape(&back.to_node()));
}

#[test]
fn exceedance_and_confidence_round_trip() {
    let q = Query::new(
        Qoi::Exceedance {
            field: "vm".into(),
            region: "bracket".into(),
            threshold: 180.0e6,
            threshold_dims: PA,
            environment: "site-hazard".into(),
        },
        Target::ToleranceAndConfidence {
            value: 1.0e6,
            dims: Dims::NONE, // exceedance value is a probability; tolerance on it is dimensionless
            confidence: 0.95,
        },
        50.0,
        30.0,
    );
    let back = Query::from_node(&q.to_node()).expect("parses");
    assert_eq!(back.qoi, q.qoi);
    assert_eq!(back.target, q.target);
}

#[test]
fn non_query_form_is_a_teaching_error() {
    let node = Node::synthetic(fs_ir::NodeKind::List(vec![Node::synthetic(
        fs_ir::NodeKind::Symbol("study".into()),
    )]));
    let err = Query::from_node(&node).expect_err("not a query");
    // the error teaches how to shape a query.
    assert!(
        err.hint.contains("query"),
        "hint teaches the query form: {}",
        err.hint
    );
}
