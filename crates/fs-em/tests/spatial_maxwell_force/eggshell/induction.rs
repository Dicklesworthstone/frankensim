//! A compact real current loop distinguishes retained induction from a static reset.
use super::*;
use fs_em::spatial::vector_potential::{
    induction::{
        SpatialInduction, SpatialInductionControls, SpatialInductionLimits,
        nonlinear::{SpatialNonlinearInduction, SpatialNonlinearInductionControls},
    },
    nonlinear::{NonlinearVectorControls, SpatialNonlinearVectorMagnetostatics},
};
use fs_material::{
    conductor::{
        ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY, ResolvedConductorState,
    },
    field_storage::MAGNETIC_PERMEABILITY_PROPERTY,
};

fn conducting_source() -> (SpatialMagneticMaterial, ResolvedConductorState) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1.0,
        ),
        (
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            5.0,
        ),
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
                    source: "manufactured conducting current loop".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "induction force fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let source = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor = ResolvedConductorState::try_from_material(source.clone()).unwrap();
    (SpatialMagneticMaterial::Linear(source), conductor)
}
fn preparation() -> SpatialInductionLimits {
    SpatialInductionLimits {
        maximum_unknowns: 20000,
        maximum_conducting_cells: 4000,
        maximum_assembly_terms: 1_000_000,
    }
}
fn evolution() -> SpatialInductionControls {
    SpatialInductionControls {
        relative_tolerance: 2e-10,
        absolute_residual_a: 1e-11,
        restart: 64,
        maximum_cycles: 100,
        maximum_krylov_scalars: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 3e-8,
        maximum_cell_magnetic_peclet: 1.0,
    }
}
fn interior(problem: &SpatialVectorMagnetostatics) -> Vec<usize> {
    problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .filter_map(|(cell, &tet)| {
            center(problem.positions_m(), tet)
                .iter()
                .all(|x| x.abs() < 1.0)
                .then_some(cell)
        })
        .collect()
}

#[test]
fn retained_current_loop_torque_survives_source_removal_then_decays_in_both_owners() {
    let (positions, tets) = mesh(8);
    let origin = positions.iter().position(|&p| p == [0.0; 3]).unwrap() as u32;
    let (metal, conductor) = conducting_source();
    let air = SpatialMagneticMaterial::Linear(material(1.0));
    let mut materials = Vec::new();
    let mut conductors = Vec::new();
    let mut currents = Vec::new();
    let mut moment = 0.0;
    for &tet in &tets {
        let active = tet.contains(&origin);
        materials.push(if active { metal.clone() } else { air.clone() });
        conductors.push(active.then(|| conductor.clone()));
        if active {
            // N is the original central-vertex P1 hat. J=curl(N e_z) has
            // continuous normal trace, zero divergence and compact support.
            // Compute grad N from its zero face, independently of FEEC data.
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != origin).collect();
            let point = |v: u32| {
                let p = positions[v as usize];
                Vec3::new(p[0], p[1], p[2])
            };
            let a = point(face[0]);
            let normal = (point(face[1]) - a).cross(point(face[2]) - a);
            let gradient = normal.scale(1.0 / normal.dot(point(origin) - a));
            currents.push([gradient.y, -gradient.x, 0.0]);
            // For compact J=curl(N e_z), m=(1/2) integral r cross J
            // equals e_z integral N, and integral_tet N=V/4 exactly.
            moment += signed_volume(&positions, tet).abs() / 4.0;
        } else {
            currents.push([0.0; 3]);
        }
    }
    assert!((moment - 0.125).abs() < 2e-15);
    let problem =
        SpatialVectorMagnetostatics::new(positions, tets, materials, limits(), &mut || false)
            .unwrap();
    let boundary = vector_trace(&problem, |p| [0.0, 0.0, 0.4 * p[1]]);
    let initial = problem
        .solve(&currents, &boundary, vector_controls(), &mut || false)
        .unwrap();
    let mut state = SpatialInduction::new(
        &initial,
        conductors.clone(),
        preparation(),
        2e-9,
        &mut || false,
    )
    .unwrap();
    let shell = problem
        .prepare_eggshell(
            &weights(problem.positions_m(), false),
            force_limits(),
            &mut || false,
        )
        .unwrap();
    let selected = interior(&problem);
    let shell_load = |s: &SpatialInduction<'_>| {
        s.eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false)
            .unwrap()
    };
    let initial_load = shell_load(&state);
    let expected = Wrench::new(Vec3::new(0.0, 0.4 * moment, 0.0), Vec3::default());
    eprintln!("current-loop moment={moment} initial={initial_load:?} reference={expected:?}");
    // Raw coarse-mesh Maxwell stress in a finite box is an estimate, not the
    // infinite-space exact moment law. This still detects sign/scale errors.
    assert_wrench(initial_load, expected, 0.02 * expected.torque.y);
    assert!(initial_load.force.norm() < 2e-8);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0; 3]; problem.mesh().vertex_count];
    let duration = 0.02;
    let first = state
        .step_candidate(&zero, &velocity, 2.0 * duration, evolution(), &mut || false)
        .unwrap();
    assert!(first.step().heat_j > 0.0);
    let retained = shell_load(first.state());
    assert!(retained.torque.y > 0.3 * initial_load.torque.y);
    assert!(retained.torque.y < initial_load.torque.y);
    let surface = first
        .state()
        .maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert!(surface.torque.y > 0.2 * expected.torque.y);
    assert!(surface.force.norm() < 2e-8);
    // The identical endpoint source statically re-solved has no loop current
    // and no torque. A load bridge that resets history cannot pass this test.
    let reset = problem
        .solve(&zero, &boundary, vector_controls(), &mut || false)
        .unwrap();
    assert_wrench(
        reset
            .eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false)
            .unwrap(),
        Wrench::default(),
        2e-8,
    );

    let nonlinear = SpatialNonlinearVectorMagnetostatics::new(
        problem.positions_m().to_vec(),
        problem.mesh().tets.clone(),
        problem
            .materials()
            .iter()
            .map(|m| match m {
                SpatialMagneticMaterial::Linear(m) => SpatialNonlinearMaterial::Linear(m.clone()),
                SpatialMagneticMaterial::Recoil(_) => unreachable!(),
            })
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let nonlinear_initial = nonlinear
        .solve(
            &currents,
            &boundary,
            initial.potential_wb(),
            NonlinearVectorControls {
                field: vector_controls(),
                relative_residual_tolerance: 2e-10,
                maximum_newton_iterations: 20,
                maximum_backtracks: 20,
            },
            &mut || false,
        )
        .unwrap();
    let nonlinear_state = SpatialNonlinearInduction::new(
        &nonlinear_initial,
        conductors,
        preparation(),
        2e-9,
        &mut || false,
    )
    .unwrap();
    // The linear owner's ramp J_old -> 0 over 2dt, starting at equilibrium,
    // has exactly the same endpoint as backward Euler with zero current over
    // dt: (M+dt K) delta_A_next=M delta_A_old. Match physical fields despite
    // their different time and source-work conventions; do not equate heat.
    let nonlinear_first = nonlinear_state
        .step_candidate(
            &zero,
            &velocity,
            duration,
            SpatialNonlinearInductionControls {
                induction: evolution(),
                maximum_newton_iterations: 20,
                maximum_backtracks: 20,
            },
            &mut || false,
        )
        .unwrap();
    let nonlinear_shell = nonlinear
        .prepare_eggshell(
            &weights(nonlinear.positions_m(), false),
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert_wrench(
        nonlinear_first
            .state()
            .eggshell_wrench(&nonlinear_shell, &Se3::identity(), [0.0; 3], &mut || false)
            .unwrap(),
        retained,
        3e-8,
    );
    assert_wrench(
        nonlinear_first
            .state()
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap(),
        surface,
        3e-8,
    );

    state = first.into_parts().0;
    let mut previous_torque = retained.torque.y;
    for _ in 0..4 {
        let next = state
            .step_candidate(&zero, &velocity, duration, evolution(), &mut || false)
            .unwrap();
        state = next.into_parts().0;
        let load = shell_load(&state);
        assert!(load.torque.y > 0.0 && load.torque.y < previous_torque);
        assert!(load.force.norm() < 2e-8);
        previous_torque = load.torque.y;
    }
    assert!(previous_torque < 0.8 * retained.torque.y);
    eprintln!(
        "retained torque={} later torque={previous_torque}",
        retained.torque.y
    );

    let original = shell_load(&state);
    let pose = Se3::from_parts(
        So3::exp(So3Tangent {
            angular: Vec3::new(0.3, -0.2, 0.4),
        })
        .unwrap(),
        Vec3::new(4.0, -3.0, 2.0),
    )
    .unwrap();
    let origin = Vec3::new(0.4, -0.7, 0.2);
    let transformed = original.transform_by(&pose);
    let expected = Wrench::new(
        transformed.torque - origin.cross(transformed.force),
        transformed.force,
    );
    let mut polls = 0;
    let moved = state
        .eggshell_wrench(&shell, &pose, [origin.x, origin.y, origin.z], &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert_wrench(moved, expected, 2e-11);
    assert!(polls > problem.mesh().tets.len());
    let before = state.potential_wb().to_vec();
    let time = state.time_s();
    let mut seen = 0;
    assert!(matches!(
        state.eggshell_wrench(&shell, &pose, [origin.x, origin.y, origin.z], &mut || {
            seen += 1;
            seen == polls - 1
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    assert_eq!(state.potential_wb(), before);
    assert_eq!(state.time_s(), time);
    assert_eq!(state.impressed_current_density_a_m2(), zero);
    assert_eq!(shell_load(&state), original);
}

#[test]
fn transient_collars_refuse_idle_conductors_impressed_current_and_different_original_sources() {
    let (positions, tets) = mesh(4);
    // This tetrahedron lies immediately outside x=1, with its x=1 face
    // shared by the closed interior region. Both load functionals encounter it.
    let collar = tets
        .iter()
        .position(|tet| {
            tet.iter()
                .filter(|&&v| positions[v as usize][0] == 1.0)
                .count()
                == 3
                && center(&positions, *tet)[0] > 1.0
                && tet.iter().all(|&v| {
                    positions[v as usize][1].abs() <= 1.0 && positions[v as usize][2].abs() <= 1.0
                })
        })
        .unwrap();
    // All scenarios need a real conductor for induction preparation. Select
    // an enclosed tetrahedron that shares no face with the surface collar.
    let core = tets
        .iter()
        .position(|&tet| {
            center(&positions, tet).iter().all(|x| x.abs() < 1.0)
                && (0..3).all(|axis| {
                    [-1.0, 1.0].into_iter().all(|boundary| {
                        tet.iter()
                            .filter(|&&v| positions[v as usize][axis] == boundary)
                            .count()
                            < 3
                    })
                })
        })
        .unwrap();
    let card = common::source(
        PropertyValue::Scalar {
            value: 1.0,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        ValidityDomain::unconstrained(),
    );
    let changed_query = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    let distinct = common::resolve(&card, &changed_query);
    let (metal, conductor) = conducting_source();
    for scenario in 0..4 {
        let mut materials = vec![SpatialMagneticMaterial::Linear(material(1.0)); tets.len()];
        let mut conductors = vec![None; tets.len()];
        materials[core] = metal.clone();
        conductors[core] = Some(conductor.clone());
        if scenario == 0 {
            materials[collar] = metal.clone();
            conductors[collar] = Some(conductor.clone());
        } else if scenario == 2 {
            materials[collar] = SpatialMagneticMaterial::Linear(distinct.clone());
        }
        let problem = SpatialVectorMagnetostatics::new(
            positions.clone(),
            tets.clone(),
            materials,
            limits(),
            &mut || false,
        )
        .unwrap();
        let current = scenario == 1;
        let boundary = vector_trace(&problem, |p| {
            [
                0.0,
                0.0,
                if current {
                    -0.25 * (p[0] * p[0] + p[1] * p[1])
                } else {
                    0.4 * p[1]
                },
            ]
        });
        let currents = vec![[0.0, 0.0, if current { 1.0 } else { 0.0 }]; tets.len()];
        let initial = problem
            .solve(&currents, &boundary, vector_controls(), &mut || false)
            .unwrap();
        let state = SpatialInduction::new(&initial, conductors, preparation(), 2e-9, &mut || false)
            .unwrap();
        // In scenario 0 all impressed currents are zero and the accepted field
        // is stationary. Conductivity itself must still prohibit the collar.
        let shell = problem
            .prepare_eggshell(
                &weights(problem.positions_m(), false),
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let distributed = state.eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false);
        let surface = state.maxwell_wrench(
            &interior(&problem),
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        );
        if scenario == 3 {
            // The identical fixture with a clean insulating collar is valid.
            assert_wrench(distributed.unwrap(), Wrench::default(), 2e-8);
            assert_wrench(surface.unwrap(), Wrench::default(), 2e-8);
        } else {
            assert!(
                matches!(distributed, Err(SpatialMagneticError::Invalid(_))),
                "eggshell admitted collar scenario {scenario}"
            );
            assert!(
                matches!(surface, Err(SpatialMagneticError::Invalid(_))),
                "closed surface admitted collar scenario {scenario}"
            );
        }
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(state.impressed_current_density_a_m2(), currents);
    }
}
