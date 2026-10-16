//! G3 relation-engine battery (bead frankensim-2uce): canonical declarations,
//! finite tolerance semantics, relation reuse, joint input/transform shrinking,
//! and structured final-counterexample receipts.

use fs_propcheck::metamorphic::{
    CanonicalRelation, RelationCase, Tolerance, ToleranceError, adjoint_finite_difference,
    check_relation, conversion_path_independence, refinement_monotonicity,
    regime_scaling_coherence, rigid_motion, unit_rescaling,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-propcheck/metamorphic\",\"case\":\"{case}\",\
         \"verdict\":\"pass\",\"detail\":\"{detail}\"}}"
    );
}

#[test]
fn canonical_relation_labels_and_tolerances_are_fail_closed() {
    let transform = |input: &i64, _: &i64| *input;
    let compare = |base: &i64, transformed: &i64, _: &i64, tolerance: Tolerance| {
        tolerance.evaluate_scalar(*base as f64, *transformed as f64)
    };
    let rigid = rigid_motion("rigid", Tolerance::Exact, transform, compare);
    let units = unit_rescaling("units", Tolerance::Exact, transform, compare);
    let refine = refinement_monotonicity("refine", Tolerance::Exact, transform, compare);
    let adjoint = adjoint_finite_difference("adjoint-fd", Tolerance::Exact, transform, compare);
    let paths = conversion_path_independence("paths", Tolerance::Exact, transform, compare);
    let regimes = regime_scaling_coherence("regimes", Tolerance::Exact, transform, compare);
    assert_eq!(rigid.kind(), CanonicalRelation::RigidMotion);
    assert_eq!(units.kind(), CanonicalRelation::UnitRescaling);
    assert_eq!(refine.kind(), CanonicalRelation::RefinementMonotonicity);
    assert_eq!(adjoint.kind(), CanonicalRelation::AdjointFiniteDifference);
    assert_eq!(paths.kind(), CanonicalRelation::ConversionPathIndependence);
    assert_eq!(regimes.kind(), CanonicalRelation::RegimeScalingCoherence);

    assert!(Tolerance::Exact.evaluate_scalar(1.0, 1.0).admitted());
    assert!(
        !Tolerance::Exact.evaluate_scalar(0.0, -0.0).admitted(),
        "exact relation semantics preserve signed-zero bits"
    );
    assert!(
        Tolerance::Absolute { max_abs: 0.25 }
            .evaluate_scalar(2.0, 2.25)
            .admitted()
    );
    assert!(
        Tolerance::AbsoluteRelative {
            max_abs: 1e-12,
            max_relative: 0.01,
        }
        .evaluate_scalar(100.0, 101.0)
        .admitted()
    );
    assert!(
        Tolerance::NonIncreasing { max_increase: 0.1 }
            .evaluate_scalar(5.0, 5.100_000_1)
            .margin()
            < 0.0
    );
    assert!(
        !Tolerance::Absolute { max_abs: 1.0 }
            .evaluate_scalar(f64::NAN, 0.0)
            .admitted()
    );
    assert_eq!(
        Tolerance::Absolute { max_abs: -1.0 }.validate(),
        Err(ToleranceError::Negative("max_abs"))
    );
    assert_eq!(
        Tolerance::AbsoluteRelative {
            max_abs: 0.0,
            max_relative: f64::INFINITY,
        }
        .validate(),
        Err(ToleranceError::NonFinite("max_relative"))
    );
    verdict(
        "canonical-tolerances",
        "six stable relation kinds; exact/absolute/relative/monotone limits fail closed",
    );
}

#[test]
fn one_declared_path_relation_is_reused_across_operators() {
    let relation = conversion_path_independence(
        "sign-symmetric-route",
        Tolerance::Exact,
        |input: &i64, route: &i64| if *route < 0 { -*input } else { *input },
        |base: &i64, transformed: &i64, _: &i64, tolerance: Tolerance| {
            tolerance.evaluate_scalar(*base as f64, *transformed as f64)
        },
    );
    let generate = |stream: &mut fs_propcheck::Stream| {
        let input = stream.int_in(-1_000, 1_000);
        let route = if stream.next_u64().is_multiple_of(2) {
            -1
        } else {
            1
        };
        RelationCase::new(input, route)
    };
    let absolute = |input: &i64| input.abs();
    let square = |input: &i64| input * input;

    check_relation(
        "absolute-value",
        0x63_33_41,
        256,
        generate,
        &absolute,
        &relation,
    );
    check_relation(
        "integer-square",
        0x63_33_42,
        256,
        generate,
        &square,
        &relation,
    );
    verdict(
        "relation-reuse",
        "one path-independence declaration exercised two operators across 512 cases",
    );
}

#[test]
fn planted_relation_violation_shrinks_input_and_transform_and_receipts_context() {
    let relation = rigid_motion(
        "planted-translation",
        Tolerance::Exact,
        |input: &i64, delta: &i64| {
            *input + *delta + if *input >= 100 && *delta >= 7 { 1 } else { 0 }
        },
        |base: &i64, transformed: &i64, delta: &i64, tolerance: Tolerance| {
            tolerance.evaluate_scalar((*base + *delta) as f64, *transformed as f64)
        },
    );
    let identity = |input: &i64| *input;
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        check_relation(
            "identity",
            0x6333_0BAD,
            1,
            |_| RelationCase::new(100_000_i64, 5_000_i64),
            &identity,
            &relation,
        );
    }));
    let payload = result.expect_err("the planted relation break must fail");
    let message = payload
        .downcast_ref::<String>()
        .cloned()
        .unwrap_or_else(|| "non-string panic".to_string());
    assert!(
        message.contains("input: 100") && message.contains("transform: 7"),
        "joint shrink kernel missing: {message}"
    );
    assert!(
        message.contains("relation-violation"),
        "failure kind: {message}"
    );

    let artifact_marker = "replay artifact: ";
    let artifact_tail = message
        .split_once(artifact_marker)
        .expect("successful failure-artifact path")
        .1;
    let artifact_path = artifact_tail
        .rsplit_once(", ")
        .expect("artifact path before shrink count")
        .0;
    let artifact = std::fs::read_to_string(artifact_path).expect("read relation replay artifact");
    let row = artifact
        .lines()
        .rev()
        .find(|row| row.contains("\"property\":\"identity::planted-translation\""))
        .expect("relation failure row");
    for field in [
        "\"failure_kind\":\"relation-violation\"",
        "\"relation_id\":\"planted-translation\"",
        "\"relation_kind\":\"rigid-motion\"",
        "\"transform\":\"7\"",
        "\"base_output\":\"100\"",
        "\"transformed_output\":\"108\"",
        "\"tolerance\":\"Exact\"",
        "\"margin\":\"-1.00000000000000000e0\"",
    ] {
        assert!(row.contains(field), "missing {field} in {row}");
    }
    verdict(
        "joint-relation-shrink",
        "planted break minimized to input=100/transform=7 with structured final outputs",
    );
}
