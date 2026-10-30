//! G1/G3/G4: source density reaches actual equilibrium and catalog sizing.
//! Manufactured material values test the nodal truss model, not a building code.
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Pressure};
use fs_truss::material::{MaterialLayout, TrussMaterial, body_load::BodyLoadDesignError};
use fs_truss::{GroundLimits, GroundStructure, LayoutCase, LayoutLimits, PdhgSettings};

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate);
    });
}

fn material(density: f64) -> TrussMaterial {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, density),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 100.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic body-load fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "axial test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    TrussMaterial::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn settings() -> PdhgSettings {
    PdhgSettings {
        max_iters: 50000,
        gap_tol: 1e-9,
        check_every: 50,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8 * b.abs().max(1.0), "{a} != {b}");
}
fn vertical(cx: &Cx<'_>, hanging: bool) -> (GroundStructure, LayoutCase) {
    let nodes = if hanging {
        [[0.0, 1.0], [0.0, 0.0]]
    } else {
        [[0.0, 0.0], [0.0, 1.0]]
    };
    let ground =
        GroundStructure::try_from_parts(&nodes, &[(0, 1)], &[1.0], GroundLimits::default(), cx)
            .unwrap();
    let case = LayoutCase::try_new(
        vec![[true, true], [true, false]],
        vec![[0.0, 2.0], [0.0, -1.0]],
        2,
    )
    .unwrap();
    (ground, case)
}

#[test]
fn heterogeneous_member_mass_preserves_nodal_force_moment_and_supported_loads() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[1.0, 2.0], [4.0, 2.0], [4.0, 6.0]],
            &[(0, 1), (0, 2), (1, 2)],
            &[3.0, 5.0, 4.0],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        let case = LayoutCase::try_new(
            vec![[true, true], [false, false], [true, true]],
            vec![[12.0, -7.0], [1.0, 2.0], [-3.0, 4.0]],
            3,
        )
        .unwrap();
        let states = [material(2.0), material(8.0)];
        let layout = MaterialLayout::try_new(
            &ground,
            &case,
            &states,
            &[0, 0, 1],
            LayoutLimits::default(),
            cx,
        )
        .unwrap();
        let body = layout
            .body_loads(&[0.2, 0.0, 0.3], [3.0, -4.0], cx)
            .unwrap();
        for (actual, expected) in body.member_mass_kg.iter().zip([1.2, 0.0, 9.6]) {
            near(*actual, expected);
        }
        near(body.mass_kg, 10.8);
        near(body.total_force_n[0], 32.4);
        near(body.total_force_n[1], -43.2);
        for (actual, expected) in
            body.nodal_forces_n
                .iter()
                .zip([[1.8, -2.4], [16.2, -21.6], [14.4, -19.2]])
        {
            near(actual[0], expected[0]);
            near(actual[1], expected[1]);
        }
        let moment: f64 = ground
            .nodes()
            .iter()
            .zip(&body.nodal_forces_n)
            .map(|(x, f)| x[0] * f[1] - x[1] * f[0])
            .sum();
        // Uniform member centers (2.5,2) and (4,4), independently integrated.
        near(moment, 1.2 * (-16.0) + 9.6 * (-28.0));
        assert_eq!(layout.load_case(), &case);
        let absent = layout.body_loads(&[0.0; 3], [3.0, -4.0], cx).unwrap();
        assert_eq!(absent.mass_kg, 0.0);
        assert_eq!(absent.nodal_forces_n, vec![[0.0; 2]; 3]);
    });
}

#[test]
fn density_only_changes_reach_self_weight_sections_and_actual_tension_or_compression() {
    with_cx(|cx, _| {
        let catalog = [
            0.011, 0.014, 0.018, 0.024, 0.04, 0.05, 0.06, 0.07, 0.08, 0.16,
        ];
        for hanging in [true, false] {
            let (ground, case) = vertical(cx, hanging);
            let mut outcomes = Vec::new();
            for density in [2.0, 8.0] {
                let states = [material(density)];
                let layout = MaterialLayout::try_new(
                    &ground,
                    &case,
                    &states,
                    &[0],
                    LayoutLimits::default(),
                    cx,
                )
                .unwrap();
                let result = layout
                    .solve_with_body_loads([0.0, -9.81], &catalog, 1e-3, settings(), 8, cx)
                    .unwrap();
                let expected_area = match (hanging, density == 2.0) {
                    (true, true) => 0.014,
                    (true, false) => 0.018,
                    (false, true) => 0.05,
                    (false, false) => 0.07,
                };
                assert_eq!(result.areas_m2, vec![expected_area]);
                let mass = density * expected_area;
                let half_weight = 0.5 * mass * 9.81;
                near(result.body_loads.mass_kg, mass);
                near(result.design.mass_kg.unwrap(), mass);
                near(result.load_case.loads()[0][1], 2.0 - half_weight);
                near(result.load_case.loads()[1][1], -1.0 - half_weight);
                assert_eq!(result.load_case.supported(), case.supported());
                let expected_force = if hanging {
                    1.0 + half_weight
                } else {
                    -1.0 - half_weight
                };
                near(result.design.audit.members[0].force, expected_force);
                near(
                    result.solution.forces[0] - result.solution.forces[1],
                    expected_force,
                );
                if hanging {
                    assert_eq!(result.design.audit.members[0].area_buckling, 0.0);
                } else {
                    assert!(result.design.audit.members[0].area_buckling > 0.0);
                }
                assert!(result.iterations >= 2);
                assert!(result.design.audit.all_pass);
                // Re-evaluating the accepted sections gives the very same body
                // forces: no earlier iteration's gravity remains in the case.
                assert_eq!(
                    layout
                        .body_loads(&result.areas_m2, [0.0, -9.81], cx)
                        .unwrap(),
                    result.body_loads
                );
                assert_eq!(layout.load_case(), &case);
                outcomes.push((expected_area, result.load_case.loads()[1][1]));
            }
            assert!(outcomes[1].0 > outcomes[0].0);
            assert!(outcomes[1].1 < outcomes[0].1);
        }
    });
}

#[test]
fn numerical_shape_catalog_iteration_and_cancellation_fail_without_changing_the_case() {
    with_cx(|cx, gate| {
        let (ground, case) = vertical(cx, true);
        let states = [material(2.0)];
        let layout =
            MaterialLayout::try_new(&ground, &case, &states, &[0], LayoutLimits::default(), cx)
                .unwrap();
        for areas in [vec![], vec![-1.0], vec![f64::NAN], vec![f64::MAX]] {
            assert!(layout.body_loads(&areas, [0.0, -9.81], cx).is_err());
        }
        assert!(layout.body_loads(&[0.1], [f64::INFINITY, 0.0], cx).is_err());
        assert!(
            layout
                .body_loads(&[f64::from_bits(1)], [0.0, f64::MIN_POSITIVE], cx)
                .is_err()
        );
        assert!(matches!(
            layout.solve_with_body_loads([0.0, -9.81], &[0.011], 1e-3, settings(), 8, cx),
            Err(BodyLoadDesignError::InvalidSizing { iteration: 1 })
        ));
        assert!(matches!(
            layout.solve_with_body_loads([0.0, -9.81], &[0.011, 0.014], 1e-3, settings(), 1, cx),
            Err(BodyLoadDesignError::IterationLimit {
                maximum_iterations: 1
            })
        ));
        assert!(
            layout
                .solve_with_body_loads([0.0, -9.81], &[0.011, 0.014], 1e-3, settings(), 0, cx)
                .is_err()
        );
        assert!(matches!(
            layout.solve_with_body_loads(
                [0.0, -9.81],
                &[0.011, 0.014],
                1e-3,
                PdhgSettings {
                    max_iters: 1,
                    check_every: 1,
                    ..settings()
                },
                8,
                cx
            ),
            Err(BodyLoadDesignError::SolverNotConverged { iteration: 0, .. })
        ));
        gate.request();
        assert!(layout.body_loads(&[0.1], [0.0, -9.81], cx).is_err());
        assert!(matches!(
            layout.solve_with_body_loads([0.0, -9.81], &[0.011, 0.014], 1e-3, settings(), 8, cx),
            Err(BodyLoadDesignError::Construction(
                fs_truss::TrussConstructionError::Cancelled { .. }
            ))
        ));
        assert_eq!(layout.load_case(), &case);
        assert_eq!(case.loads(), &[[0.0, 2.0], [0.0, -1.0]]);
    });
}
