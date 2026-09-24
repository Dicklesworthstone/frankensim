use super::*;

const FIXTURE: &str = include_str!("../../../../../examples/marquee/bracket-3d-adaptive.fsim");

#[test]
fn g0_sdf3_schema_preserves_all_inputs_and_rejects_unknown_or_unbounded_work() {
    let spec = spec::parse(FIXTURE).expect("complete explicit 3-D input");
    assert_eq!(spec.schedule.len(), 3);
    assert_eq!(spec.loads.len(), 2);
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    for bad in [
        FIXTURE.replace(":maximum-leaves 2048", ":maximum-leaves 1000000000"),
        FIXTURE.replace(":quadrature-points 1000000", ":quadrature-points 4000000"),
        FIXTURE.replace(":filter-radius-m 0.15", ":filter-radius-m 0.0"),
        FIXTURE.replace(":weight 0.7", ":weight -0.7"),
        FIXTURE.replace(":height-m 0.7", ":height-m 0.7 :ignored-geometry 4"),
        FIXTURE.replace(":poissons-ratio 0.3", ":poissons-ratio 0.49"),
        FIXTURE.replace(":initial-density 0.5", ":initial-density 1.1"),
        FIXTURE.replace(":volume-fraction 0.5", ":volume-fraction 0.0"),
        FIXTURE.replace(":schedule ((1.0 1.0) (2.0 2.0) (3.0 8.0))", ":schedule ()"),
    ] {
        assert!(
            spec::parse(&bad).is_err(),
            "must refuse invalid declared physics or work"
        );
    }
    let changed = spec::parse(&FIXTURE.replace(":weight 0.7", ":weight 0.8")).unwrap();
    assert_ne!(
        spec.id, changed.id,
        "independent load weights bind study identity"
    );
}

#[test]
fn g4_sdf3_precancellation_performs_no_geometry_or_equilibrium() {
    let mut spec = spec::parse(FIXTURE).unwrap();
    let gate = CancelGate::new();
    gate.request();
    let error = match compute(&spec, &gate) {
        Ok(_) => panic!("pre-cancelled study must refuse"),
        Err(error) => error,
    };
    assert_eq!(error.exit, exit::CANCELLED);
    assert_eq!(error.code, "cli-study-sdf3-cancelled");
    spec.wall_s = 1e-99;
    let error = match compute(&spec, &CancelGate::new()) {
        Ok(_) => panic!("expired wall allowance must refuse"),
        Err(error) => error,
    };
    assert_eq!(error.exit, exit::BUDGET);
    assert_eq!(error.code, "cli-study-sdf3-wall-budget");
}

#[test]
fn g4_sdf3_enrichment_setup_budget_is_not_a_numerical_failure() {
    let setup = |error| {
        AdaptiveContinuationError3::Goal(GoalRefinementError3::Preconditioner(
            AdaptivePreconditionError3::Coarse(error),
        ))
    };
    for reason in [
        "adaptive coarse space/setup",
        "transfer dimensions",
        "vector interpolation entries",
        "Galerkin operator applications",
    ] {
        assert!(refinement_budget(&setup(TwoLevelError::Budget(reason))));
    }
    for error in [
        TwoLevelError::NotPositiveDefinite,
        TwoLevelError::Nonsymmetric,
        TwoLevelError::Invalid("nonfinite fine operator action"),
        TwoLevelError::Cancelled,
    ] {
        assert!(!refinement_budget(&setup(error)));
    }
    assert!(refinement_budget(&AdaptiveContinuationError3::Background(
        OctreeError3::LeafBudget
    )));
    assert!(refinement_budget(&AdaptiveContinuationError3::Goal(
        GoalRefinementError3::Physics(ElasticityError3::Quadrature(QuadratureError3::PointBudget))
    )));
}
