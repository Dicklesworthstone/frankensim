//! Conformance + admission battery for the declarative query language v0
//! (addendum Proposal 8). Covers: the wedge-vertical QoI menu is
//! expressible; QoI metadata is correct; well-posed queries admit; the five
//! ill-posed classes each reject with a distinct teaching finding;
//! dimensional consistency is enforced; admission is deterministic; and the
//! `(query …)` IR surface round-trips under `same_shape`.

use fs_ir::admission::Severity;
use fs_ir::query::{FieldRegistry, Qoi, Query, Target};
use fs_ir::{Node, NodeKind, Span, VersionedProgram, json, sexpr};
use fs_qty::Dims;

// von Mises stress: Pa = kg·m⁻¹·s⁻²  → [M,KG,S,K,A,MOL]
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
    assert_eq!(m.value_dims(PA), Some(PA));
    // integral multiplies by volume;
    let i = Qoi::Integral {
        field: "temperature".into(),
        region: "r".into(),
    };
    assert_eq!(i.value_dims(KELVIN), KELVIN.checked_plus(VOLUME),);
    // An integral over a field already at the exponent ceiling refuses
    // typed instead of clamping (bead sj31i.11).
    let ceiling = Qoi::Integral {
        field: "hostile".into(),
        region: "r".into(),
    };
    assert_eq!(ceiling.value_dims(Dims([127, 0, 0, 0, 0, 0])), None);
    // exceedance is a dimensionless probability.
    let e = Qoi::Exceedance {
        field: "vm".into(),
        region: "r".into(),
        threshold: 1.0,
        threshold_dims: PA,
        environment: "h".into(),
    };
    assert_eq!(e.value_dims(PA), Some(Dims::NONE));
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
            dims: KELVIN.checked_plus(VOLUME).expect("in-range exponents"),
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
    let node = q.to_node().expect("query serializes");
    let back = Query::from_node(&node).expect("query parses");
    assert_eq!(back.qoi, q.qoi);
    assert_eq!(back.target, q.target);
    // exact round-trip of the same literals: compare by bits (crate idiom).
    assert_eq!(back.budget_usd.to_bits(), q.budget_usd.to_bits());
    assert_eq!(back.deadline_s.to_bits(), q.deadline_s.to_bits());
    // and the emitted node equals a re-emission (same_shape isomorphism).
    assert!(node.same_shape(&back.to_node().expect("query reserializes")));
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
    let back = Query::from_node(&q.to_node().expect("query serializes")).expect("parses");
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

fn assert_query_semantics(expected: &Query, actual: &Query) {
    match (&expected.qoi, &actual.qoi) {
        (
            Qoi::Exceedance {
                field: expected_field,
                region: expected_region,
                threshold: expected_threshold,
                threshold_dims: expected_dims,
                environment: expected_environment,
            },
            Qoi::Exceedance {
                field: actual_field,
                region: actual_region,
                threshold: actual_threshold,
                threshold_dims: actual_dims,
                environment: actual_environment,
            },
        ) => {
            assert_eq!(actual_field, expected_field);
            assert_eq!(actual_region, expected_region);
            assert_eq!(actual_threshold.to_bits(), expected_threshold.to_bits());
            assert_eq!(actual_dims, expected_dims);
            assert_eq!(actual_environment, expected_environment);
        }
        _ => assert_eq!(actual.qoi, expected.qoi),
    }
    match (&expected.target, &actual.target) {
        (
            Target::Tolerance {
                value: expected_value,
                dims: expected_dims,
            },
            Target::Tolerance {
                value: actual_value,
                dims: actual_dims,
            },
        ) => {
            assert_eq!(actual_value.to_bits(), expected_value.to_bits());
            assert_eq!(actual_dims, expected_dims);
        }
        (Target::Confidence(expected), Target::Confidence(actual)) => {
            assert_eq!(actual.to_bits(), expected.to_bits())
        }
        (
            Target::ToleranceAndConfidence {
                value: expected_value,
                dims: expected_dims,
                confidence: expected_confidence,
            },
            Target::ToleranceAndConfidence {
                value: actual_value,
                dims: actual_dims,
                confidence: actual_confidence,
            },
        ) => {
            assert_eq!(actual_value.to_bits(), expected_value.to_bits());
            assert_eq!(actual_dims, expected_dims);
            assert_eq!(actual_confidence.to_bits(), expected_confidence.to_bits());
        }
        _ => panic!("query target variant changed during round-trip"),
    }
    assert_eq!(actual.budget_usd.to_bits(), expected.budget_usd.to_bits());
    assert_eq!(actual.deadline_s.to_bits(), expected.deadline_s.to_bits());
}

#[test]
fn every_qoi_and_target_round_trips_through_both_syntaxes() {
    let qois = [
        Qoi::MaxOverRegion {
            field: "vm".into(),
            region: "bracket".into(),
        },
        Qoi::Integral {
            field: "temperature".into(),
            region: "volume".into(),
        },
        Qoi::Exceedance {
            field: "vm".into(),
            region: "bracket".into(),
            threshold: 180.0e6,
            threshold_dims: PA,
            environment: "site-hazard".into(),
        },
    ];
    let targets = [
        Target::Tolerance {
            value: 5.0e6,
            dims: PA,
        },
        Target::Confidence(0.95),
        Target::ToleranceAndConfidence {
            value: 0.01,
            dims: Dims::NONE,
            confidence: 0.99,
        },
    ];

    for qoi in qois {
        for target in &targets {
            let query = Query::new(qoi.clone(), target.clone(), 50.25, 30.5);
            let node = query.to_node().expect("finite query serializes");
            let artifact = VersionedProgram::try_current(node).expect("valid versioned query");
            let sexpr_bytes = artifact
                .print_sexpr_checked()
                .expect("valid versioned s-expression");
            let json_bytes = artifact
                .print_json_checked()
                .expect("valid versioned JSON mapping");
            let sexpr_back = VersionedProgram::parse_sexpr(&sexpr_bytes)
                .expect("canonical s-expression reparses");
            let json_back =
                VersionedProgram::parse_json(&json_bytes).expect("canonical JSON reparses");

            assert_eq!(sexpr_back.print_sexpr(), sexpr_bytes);
            assert_eq!(json_back.print_sexpr(), sexpr_bytes);
            assert_eq!(sexpr_back.print_json(), json_bytes);
            assert_eq!(json_back.print_json(), json_bytes);
            assert_query_semantics(
                &query,
                &Query::from_node(sexpr_back.program()).expect("s-expression query recognizes"),
            );
            assert_query_semantics(
                &query,
                &Query::from_node(json_back.program()).expect("JSON query recognizes"),
            );
        }
    }
}

#[test]
fn query_grammar_rejects_duplicate_unknown_dangling_and_trailing_fields() {
    let cases = [
        (
            "(query (max :field \"vm\" :field \"other\" :region \"r\") (confidence 0.95) (budget 50) (deadline 30s))",
            "duplicate query.qoi.max.field",
        ),
        (
            "(query (max :field \"vm\" :region \"r\" :gpu \"yes\") (confidence 0.95) (budget 50) (deadline 30s))",
            "unknown query.qoi.max.gpu",
        ),
        (
            "(query (max :field \"vm\" :region) (confidence 0.95) (budget 50) (deadline 30s))",
            "dangling query.qoi.max.region",
        ),
        (
            "(query (max :field \"vm\" :region \"r\" trailing) (confidence 0.95) (budget 50) (deadline 30s))",
            "unexpected positional argument",
        ),
        (
            "(query (max :field \"vm\" :region \"r\") (confidence 0.95 extra) (budget 50) (deadline 30s))",
            "confidence takes exactly one",
        ),
        (
            "(query (max :field \"vm\" :region \"r\") (tolerance 5MPa :confidence 0.95 :confidence 0.9) (budget 50) (deadline 30s))",
            "duplicate query.target.tolerance.confidence",
        ),
        (
            "(query (max :field \"vm\" :region \"r\") (confidence 0.95) (budget 50) (budget 60) (deadline 30s))",
            "duplicate query.budget",
        ),
        (
            "(query (max :field \"vm\" :region \"r\") (confidence 0.95) (budget 50 extra) (deadline 30s))",
            "budget takes exactly one",
        ),
        (
            "(query (max :field \"vm\" :region \"r\") (confidence 0.95) (budget 50))",
            "missing query.deadline",
        ),
        (
            "(query (max :field \"vm\" :region \"r\") (confidence 0.95) (deadline 30s) (budget 50))",
            "exact canonical schema",
        ),
        (
            "(query (max :field \"vm\" :region \"r\") (confidence 0.95) (budget 9007199254740993) (deadline 30s))",
            "cannot be represented exactly",
        ),
    ];
    for (source, expected) in cases {
        let node = sexpr::parse(source).expect("malformed query remains valid syntax");
        let error = Query::from_node(&node).expect_err("ambiguous query must refuse");
        assert!(
            error.detail.contains(expected),
            "{source}: expected {expected:?}, got {:?}",
            error.detail
        );
    }
}

#[test]
fn forged_ast_atoms_refuse_with_exact_tree_paths() {
    let invalid = Node {
        kind: NodeKind::List(vec![
            Node::synthetic(NodeKind::Symbol("query".into())),
            Node {
                kind: NodeKind::Float(f64::NAN),
                span: Span::new(10, 13),
            },
        ]),
        span: Span::new(0, 14),
    };
    let error = invalid.validate().expect_err("NaN must refuse");
    assert_eq!(error.span, Span::new(10, 13));
    assert!(error.detail.contains("$[1]"));
    assert!(sexpr::print_checked(&invalid).is_err());
    assert!(json::print_checked(&invalid).is_err());

    for kind in [
        NodeKind::Symbol("1not-a-symbol".into()),
        NodeKind::Keyword(String::new()),
        NodeKind::Qty {
            value: 2.0,
            dims: Dims([1, 0, 0, 0, 0, 0]),
            text: "1m".into(),
        },
    ] {
        assert!(Node::try_synthetic(kind).is_err());
    }

    let inverted = Node {
        kind: NodeKind::Int(1),
        span: Span::new(9, 4),
    };
    let error = inverted.validate().expect_err("inverted span must refuse");
    assert_eq!(error.span, Span::new(4, 9));

    let mut deep = Node::synthetic(NodeKind::Int(1));
    for _ in 0..258 {
        deep = Node::synthetic(NodeKind::List(vec![deep]));
    }
    assert_eq!(
        deep.validate()
            .expect_err("forged over-deep AST must refuse")
            .kind
            .code(),
        "IrTooDeep"
    );
}

#[test]
fn quantity_identity_uses_one_canonical_encoder() {
    let short = sexpr::parse("5MPa").expect("quantity parses");
    let base = sexpr::parse("5000000Pa").expect("equivalent quantity parses");
    assert!(short.same_shape(&base));
    assert_eq!(
        sexpr::print_checked(&short).expect("canonical print"),
        sexpr::print_checked(&base).expect("canonical print")
    );
    assert_eq!(
        json::print_checked(&short).expect("canonical print"),
        json::print_checked(&base).expect("canonical print")
    );

    let dimensionless = Node::quantity(0.95, Dims::NONE).expect("canonical quantity");
    let printed = sexpr::print_checked(&dimensionless).expect("canonical print");
    assert!(
        printed.ends_with("rad"),
        "dimensionless qty stays typed: {printed}"
    );
    assert!(
        sexpr::parse(&printed)
            .expect("reparse")
            .same_shape(&dimensionless)
    );
}
