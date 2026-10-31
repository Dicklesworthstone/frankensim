//! G1/G3/G4: the public frame path uses one material assignment for yield,
//! stiffness, catalog mass and self-weight. Synthetic cards qualify the numerical
//! model, not measured alloys, real buildings or a building-code allowable.

use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_frame::layout::material::{BodyLoadControls, MaterialLayoutControls, layout_and_size};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Pressure};
use fs_truss::material::{TrussMaterial, body_load::BodyLoadDesignError};
use fs_truss::{GroundLimits, GroundStructure, LayoutCase, PdhgSettings, TrussConstructionError};

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 13,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate);
    });
}

fn material(density: f64, youngs: f64, yield_stress: f64) -> TrussMaterial {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, density),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, youngs),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, yield_stress),
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
                    source: "synthetic frame material".into(),
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
            process: "frame numerical reference".into(),
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

fn controls() -> MaterialLayoutControls {
    MaterialLayoutControls {
        solver: PdhgSettings {
            max_iters: 50_000,
            gap_tol: 1e-9,
            check_every: 50,
        },
        ..MaterialLayoutControls::default()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-15),
        "{actual} != {expected}",
    );
}

fn two_members(cx: &Cx<'_>, force: f64) -> (GroundStructure, LayoutCase) {
    let ground = GroundStructure::try_from_parts(
        &[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        &[(0, 1), (1, 2)],
        &[1.0, 1.0],
        GroundLimits::default(),
        cx,
    )
    .unwrap();
    let case = LayoutCase::try_new(
        vec![[true, true], [false, true], [false, true]],
        vec![[0.0; 2], [0.0; 2], [force, 0.0]],
        3,
    )
    .unwrap();
    (ground, case)
}

#[test]
fn heterogeneous_sources_drive_yield_buckling_sizing_and_mass() {
    with_cx(|cx, _| {
        let materials = [
            material(8000.0, 200e9, 300e6),
            material(7900.0, 50e9, 100e6),
        ];
        let assignments = [0, 1];
        assert_ne!(
            materials[0].resolved().identity(),
            materials[1].resolved().identity()
        );
        for (force, catalog, expected) in [
            (1.0, [4e-9, 1.1e-8, 2e-8], [4e-9, 1.1e-8]),
            (-1.0, [3e-6, 6e-6, 1e-5], [3e-6, 6e-6]),
        ] {
            let (ground, case) = two_members(cx, force);
            let report = layout_and_size(
                &ground,
                &case,
                &materials,
                &assignments,
                &catalog,
                controls(),
                cx,
            )
            .unwrap();
            assert_eq!(report.assignments, assignments);
            assert_eq!(report.materials, materials.as_slice());
            assert!(std::ptr::eq(report.ground, &ground));
            assert_eq!(report.prescribed_case, &case);
            assert_eq!(report.controls, controls());
            let design = &report.design;
            assert_eq!(design.load_case, case);
            assert_eq!(design.areas_m2, expected);
            assert_eq!(design.iterations, 0);
            assert_eq!(design.body_loads.total_force_n, [0.0; 2]);
            assert!(design.design.audit.all_pass);
            near(design.solution.volume_m3, 1.0 / 300e6 + 1.0 / 100e6, 1e-7);
            let expected_mass = expected[0] * 8000.0 + expected[1] * 7900.0;
            near(design.design.mass_kg.unwrap(), expected_mass, 1e-12);
            near(design.body_loads.mass_kg, expected_mass, 1e-12);
            for member in &design.design.audit.members {
                near(member.force, force, 1e-8);
                let state = &report.materials[report.assignments[member.index]];
                near(
                    design.body_loads.member_mass_kg[member.index],
                    member.area_catalog * member.length * state.density_kg_m3(),
                    1e-12,
                );
            }
            if force < 0.0 {
                near(
                    design.design.audit.members[1].area_buckling,
                    2.0 * design.design.audit.members[0].area_buckling,
                    1e-12,
                );
            }
        }
    });
}

#[test]
fn density_swap_changes_actual_self_weight_load_and_catalog_section_once() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[0.0, 1.0], [0.0, 0.0]],
            &[(0, 1)],
            &[1.0],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        let case = LayoutCase::try_new(
            vec![[true, true], [true, false]],
            vec![[0.0, 2.0], [0.0, -1.0]],
            2,
        )
        .unwrap();
        let original = case.clone();
        let catalog = [0.011, 0.014, 0.018, 0.024];
        let policy = MaterialLayoutControls {
            body_loads: Some(BodyLoadControls {
                acceleration_m_s2: [0.0, -9.81],
                maximum_iterations: 8,
            }),
            ..controls()
        };
        for (density, area) in [(2.0, 0.014), (8.0, 0.018)] {
            let materials = [material(density, 1000.0, 100.0)];
            let report =
                layout_and_size(&ground, &case, &materials, &[0], &catalog, policy, cx).unwrap();
            let design = &report.design;
            let mass = density * area;
            let half_weight = 0.5 * mass * 9.81;
            assert_eq!(design.areas_m2, [area]);
            assert!(design.iterations >= 2);
            near(design.design.mass_kg.unwrap(), mass, 1e-12);
            near(design.body_loads.total_force_n[1], -mass * 9.81, 1e-12);
            near(design.load_case.loads()[0][1], 2.0 - half_weight, 1e-12);
            near(design.load_case.loads()[1][1], -1.0 - half_weight, 1e-12);
            near(
                design.design.audit.members[0].force,
                1.0 + half_weight,
                1e-8,
            );
            assert_eq!(design.load_case.supported(), case.supported());
            assert_eq!(report.prescribed_case, &original);
            assert_eq!(case, original);
        }
    });
}

#[test]
fn invalid_assignments_unconverged_solve_and_missing_catalog_refuse() {
    with_cx(|cx, _| {
        let materials = [material(8000.0, 200e9, 300e6)];
        let (ground, case) = two_members(cx, 1.0);
        let original = case.clone();
        for assignments in [&[0][..], &[0, 1][..]] {
            assert!(matches!(
                layout_and_size(
                    &ground,
                    &case,
                    &materials,
                    assignments,
                    &[4e-9],
                    controls(),
                    cx
                ),
                Err(BodyLoadDesignError::Construction(_)),
            ));
        }
        let short = MaterialLayoutControls {
            solver: PdhgSettings {
                max_iters: 1,
                check_every: 1,
                ..controls().solver
            },
            ..controls()
        };
        assert!(matches!(
            layout_and_size(&ground, &case, &materials, &[0, 0], &[4e-9], short, cx),
            Err(BodyLoadDesignError::SolverNotConverged { iteration: 0, .. }),
        ));
        assert!(matches!(
            layout_and_size(
                &ground,
                &case,
                &materials,
                &[0, 0],
                &[1e-12],
                controls(),
                cx
            ),
            Err(BodyLoadDesignError::InvalidSizing { iteration: 0 }),
        ));
        assert_eq!(case, original);
    });
}

#[test]
fn cancellation_refuses_before_publishing_source_bound_design() {
    with_cx(|cx, gate| {
        let materials = [material(8000.0, 200e9, 300e6)];
        let (ground, case) = two_members(cx, 1.0);
        let original = case.clone();
        gate.request();
        assert!(matches!(
            layout_and_size(&ground, &case, &materials, &[0, 0], &[4e-9], controls(), cx),
            Err(BodyLoadDesignError::Construction(
                TrussConstructionError::Cancelled { .. }
            )),
        ));
        assert_eq!(case, original);
    });
}
