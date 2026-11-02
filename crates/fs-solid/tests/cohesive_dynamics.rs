//! Original P1 mass/stiffness, source-bound fracture and inertial continuation.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget as ExecBudget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::cohesive::{BilinearCohesive, CohesiveViscosity};
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::Pressure;
use fs_solid::cohesive::dynamics::{DynamicsConfig, DynamicsError};
use fs_solid::cohesive::p1::{MaterialFacePair, P1Interfaces, TraceBudget};
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, Error, JumpTerm};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticMaterial, TetLinearElasticProblem,
    TetMaterialField,
};
use fs_solver::{CsrOp, LinearOp, NewtonKrylovConfig};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 93,
                kernel_id: 93,
                tile: 0,
                iteration: 0,
            },
            ExecBudget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn config() -> DynamicsConfig {
    DynamicsConfig {
        solve: Config {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-12,
                linear_restart: 24,
                forcing_maximum: 0.05,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 60,
            max_operator_evaluations: 4000,
            force_tolerance_n: 1e-8,
            local_energy_tolerance_j_m2: 1e-10,
            total_cohesive_energy_tolerance_j: 1e-10,
            maximum_displacement_m: 0.02,
        },
        difference_step_m: 1e-8,
        maximum_step_s: 0.02,
        maximum_velocity_m_s: 2.0,
        maximum_steps: 1000,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}
fn budget() -> Budget {
    Budget {
        max_dofs: 24,
        max_patches: 3,
        max_jump_terms: 54,
        max_krylov_basis_entries: 50_000,
    }
}
fn law() -> BilinearCohesive {
    BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap()
}
fn patch(law: &BilinearCohesive) -> CohesivePatch<'_> {
    CohesivePatch::new(
        7,
        0.01,
        law,
        vec![JumpTerm {
            dof: 0,
            local: [1.0, 0.0, 0.0],
        }],
    )
    .unwrap()
}
fn scalar_tet(cx: &Cx<'_>) -> TetElasticAssembly {
    // Actual consistent P1 assembly: rho=6, E=60, nu=0, V=1/6.
    // Only full DOF 3 is free: M=.1 kg, K=10 N/m.
    let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [93; 32]).unwrap();
    let fixed: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
    TetLinearElasticProblem {
        nodes_m: &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        tetrahedra: &[[0, 1, 2, 3]],
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::default(),
    }
    .assemble(cx)
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn original_consistent_mass_and_cohesive_stiffness_have_second_order_motion() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law();
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let mut errors = Vec::new();
        let velocity = 0.01;
        let frequency = 200.0_f64.sqrt(); // (10 bulk + .01*1000 bond)/.1 mass.
        for count in [25, 50, 100] {
            let mut motion = system
                .bonded_state()
                .into_dynamics(&mass, vec![velocity], 0.0, config())
                .unwrap();
            for _ in 0..count {
                let step = motion
                    .advance(&[0.0], &[[0.0; 3]], 0.25 / count as f64, cx)
                    .unwrap();
                close(
                    step.kinetic_after_j + step.bulk_stored_after_j + step.cohesive_stored_after_j,
                    0.5 * 0.1 * velocity * velocity,
                    1e-12,
                );
                close(step.energy_residual_j, 0.0, 1e-12);
                assert_eq!(step.cohesive_energy_j.dissipated, 0.0);
                assert_eq!(step.viscous_dissipated_j, 0.0);
                motion = step.motion;
            }
            assert_eq!(motion.state().points()[0].damage(), 0.0);
            close(motion.time_s(), 0.25, 1e-14);
            errors.push(
                (motion.state().displacements_m()[0]
                    - velocity / frequency * (frequency * 0.25).sin())
                .abs(),
            );
        }
        for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
            assert!((3.8..4.2).contains(&ratio), "{errors:?}");
        }
    });
}

#[test]
fn prescribed_jump_external_force_and_physical_viscosity_have_independent_work() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law().with_viscosity(CohesiveViscosity::new(4.0, 7.0, 2.0).unwrap());
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let motion = system
            .bonded_state()
            .into_dynamics(&mass, vec![0.003], 0.0, config())
            .unwrap();
        let h = 0.01;
        let g = 0.0002;
        let force = 0.02;
        let m = 0.1;
        let bulk = 10.0;
        let bond = 10.0;
        let dashpot = 0.04;
        // Independent scalar balance with midpoint elastic forces and exact
        // surviving-bond dashpot integral. All jumps stay below initiation.
        let u = (force + 2.0 * m * 0.003 / h - (0.5 * bond + dashpot / h) * g)
            / (2.0 * m / (h * h) + 0.5 * (bulk + bond) + dashpot / h);
        let v = 2.0 * u / h - 0.003;
        let mean_bond_force = 0.5 * bond * (u + g) + dashpot * (u + g) / h;
        let step = motion.advance(&[force], &[[g, 0.0, 0.0]], h, cx).unwrap();
        close(step.motion.state().displacements_m()[0], u, 1e-13);
        close(step.motion.velocity_m_s()[0], v, 1e-11);
        close(step.external_work_j, force * u, 1e-14);
        close(step.prescribed_jump_work_j, mean_bond_force * g, 1e-14);
        close(
            step.viscous_dissipated_j,
            dashpot * (u + g) * (u + g) / h,
            1e-14,
        );
        close(step.mean_traction_pa[0][0] * 0.01, mean_bond_force, 1e-11);
        close(step.energy_residual_j, 0.0, 1e-13);
        assert_eq!(step.cohesive_energy_j.dissipated, 0.0);
        // Original accepted trajectory has not been modified by the candidate.
        assert_eq!(motion.state().displacements_m(), &[0.0]);
        assert_eq!(motion.velocity_m_s(), &[0.003]);
        assert_eq!(motion.time_s(), 0.0);
    });
}

#[test]
fn adopting_a_softened_solid_preserves_damage_and_releases_its_actual_energy() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law();
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let u0 = 0.002;
        let traction = (0.004 - u0) / 0.003;
        let accepted = system
            .bonded_state()
            .solve(
                &[10.0 * u0 + 0.01 * traction],
                &[[0.0; 3]],
                config().solve,
                cx,
            )
            .unwrap()
            .state;
        close(accepted.displacements_m()[0], u0, 1e-13);
        assert!((0.0..1.0).contains(&accepted.points()[0].damage()));
        let old = accepted.points()[0];
        let motion = accepted
            .into_dynamics(&mass, vec![-0.001], 3.0, config())
            .unwrap();
        assert_eq!(motion.state().points()[0], old);
        let h = 0.001;
        let k = 10.0 + (1.0 - old.damage()) * 10.0;
        let expected = ((2.0 * 0.1 / (h * h) - 0.5 * k) * u0 + 2.0 * 0.1 * (-0.001) / h)
            / (2.0 * 0.1 / (h * h) + 0.5 * k);
        let step = motion.advance(&[0.0], &[[0.0; 3]], h, cx).unwrap();
        close(step.motion.state().displacements_m()[0], expected, 1e-12);
        assert_eq!(step.motion.state().points()[0].damage(), old.damage());
        assert_eq!(
            step.motion.state().points()[0].dissipated_j_m2(),
            old.dissipated_j_m2()
        );
        assert_eq!(step.cohesive_energy_j.dissipated, 0.0);
        close(step.energy_residual_j, 0.0, 1e-12);
        assert!(step.kinetic_after_j > step.kinetic_before_j);
    });
}

fn source() -> ResolvedBilinearCohesive {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1200.0,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 0.12),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            0.00006,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2400.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured inertial debonding; no measured calibration".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.into(),
            phase: "solid".into(),
            process: "manufactured dynamics".into(),
            revision: 1,
        },
        texture_frame: format!("{chemistry} fixed reference frame"),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "declared cohesive bond".into(),
            third_body: None,
            environment: "fixed 300 K".into(),
            history: "initially bonded".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap()
}
fn body_momenta(mass: &dyn LinearOp, velocity: &[f64]) -> [[f64; 3]; 2] {
    let mut momentum = vec![0.0; velocity.len()];
    mass.apply(velocity, &mut momentum);
    let mut bodies = [[0.0; 3]; 2];
    for (node, values) in momentum.chunks_exact(3).enumerate() {
        for c in 0..3 {
            bodies[node / 4][c] += values[c];
        }
    }
    bodies
}

#[test]
fn source_bound_deformable_p1_bodies_fracture_then_keep_separate_momenta() {
    with_cx(|cx| {
        let nodes = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., -1.],
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
        ];
        let tets = [[0, 2, 1, 3], [4, 5, 6, 7]];
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [94; 32]).unwrap();
        let assembled = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::default(),
        }
        .assemble(cx)
        .unwrap();
        close(assembled.total_mass_kg, 2.0, 1e-14);
        let source = source();
        let faces = P1Interfaces::from_interfaces(
            &nodes,
            &tets,
            &assembled.free_dofs,
            &[MaterialFacePair {
                id: 41,
                side_a: [2, 0, 1],
                side_b: [5, 6, 4],
                source: &source,
                surface_materials: [
                    &source.surfaces()[0].material,
                    &source.surfaces()[1].material,
                ],
            }],
            TraceBudget {
                maximum_nodes: 8,
                maximum_tetrahedra: 2,
                maximum_interfaces: 1,
                minimum_face_quality: 1e-10,
            },
            cx,
        )
        .unwrap();
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let system = CohesiveSystem::new(&elastic, faces.patches().to_vec(), budget()).unwrap();
        let initial: Vec<_> = (0..8)
            .flat_map(|node| {
                if node < 4 {
                    [-0.1, 0.05, -0.3]
                } else {
                    [0.2, -0.03, 0.5]
                }
            })
            .collect();
        let mut motion = system
            .bonded_state()
            .into_dynamics(&mass, initial.clone(), 0.0, config())
            .unwrap();
        let initial_momenta = body_momenta(&mass, &initial);
        let initial_energy = 0.5
            * (0.1_f64.powi(2)
                + 0.05_f64.powi(2)
                + 0.3_f64.powi(2)
                + 0.2_f64.powi(2)
                + 0.03_f64.powi(2)
                + 0.5_f64.powi(2));
        let mut fracture = 0.0;
        for _ in 0..20 {
            let step = motion
                .advance(&[0.0; 24], &[[0.0; 3]; 3], 1e-4, cx)
                .unwrap();
            fracture += step.cohesive_energy_j.dissipated;
            close(
                step.kinetic_after_j
                    + step.bulk_stored_after_j
                    + step.cohesive_stored_after_j
                    + fracture,
                initial_energy,
                2e-9,
            );
            assert_eq!(step.viscous_dissipated_j, 0.0);
            close(step.prescribed_jump_work_j, 0.0, 1e-15);
            motion = step.motion;
            let momenta = body_momenta(&mass, motion.velocity_m_s());
            for c in 0..3 {
                close(
                    momenta[0][c] + momenta[1][c],
                    initial_momenta[0][c] + initial_momenta[1][c],
                    2e-9,
                );
            }
        }
        assert!(motion.state().points().iter().all(|point| point.failed()));
        close(fracture, 0.5 * 0.00006, 1e-11);
        let separated_momenta = body_momenta(&mass, motion.velocity_m_s());
        let separation = motion.state().points()[0].jump_m()[0];
        for _ in 0..10 {
            let step = motion
                .advance(&[0.0; 24], &[[0.0; 3]; 3], 1e-4, cx)
                .unwrap();
            assert_eq!(step.cohesive_energy_j.dissipated, 0.0);
            assert!(step.mean_traction_pa.iter().all(|&mean| mean == [0.0; 3]));
            let momenta = body_momenta(&mass, step.motion.velocity_m_s());
            for body in 0..2 {
                for c in 0..3 {
                    close(momenta[body][c], separated_momenta[body][c], 2e-9);
                }
            }
            motion = step.motion;
        }
        assert!(motion.state().points()[0].jump_m()[0] > separation);
        assert!(motion.velocity_m_s().iter().any(|v| v.abs() > 0.2));
        for patch in system.patches() {
            assert_eq!(patch.source().unwrap().resolved(), source.resolved());
        }

        // Physical body forces reverse the separated bodies. The same failed
        // quadrature histories recover normal compression without rebonding,
        // removing inertia, inventing shear friction or adding fracture heat.
        let acceleration: Vec<_> = (0..8)
            .flat_map(|node| {
                if node < 4 {
                    [0.0, 0.0, 50.0]
                } else {
                    [0.0, 0.0, -50.0]
                }
            })
            .collect();
        let mut force = vec![0.0; 24];
        mass.apply(&acceleration, &mut force);
        for _ in 0..20 {
            let step = motion.advance(&force, &[[0.0; 3]; 3], 0.001, cx).unwrap();
            assert_eq!(step.cohesive_energy_j.dissipated, 0.0);
            assert_eq!(step.viscous_dissipated_j, 0.0);
            close(step.energy_residual_j, 0.0, 2e-9);
            motion = step.motion;
        }
        for point in motion.state().points() {
            assert!(point.failed() && point.jump_m()[0] < 0.0);
            assert!(point.traction_pa()[0] < 0.0);
            assert_eq!(&point.traction_pa()[1..], &[0.0; 2]);
        }
    });
}

#[test]
fn late_cancellation_domain_and_work_budget_preserve_the_complete_motion() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law();
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let motion = system
            .bonded_state()
            .into_dynamics(&mass, vec![0.01], 2.0, config())
            .unwrap();
        let saved = motion.state().points().to_vec();
        let mut complete_polls = 0;
        let expected = motion
            .advance_controlled(&[0.02], &[[0.0; 3]], 0.01, &mut || {
                complete_polls += 1;
                false
            })
            .unwrap();
        let mut polls = 0;
        assert!(matches!(
            motion.advance_controlled(&[0.02], &[[0.0; 3]], 0.01, &mut || {
                polls += 1;
                polls == complete_polls
            }),
            Err(DynamicsError::Model(Error::Cancelled))
        ));
        assert_eq!(polls, complete_polls);
        assert!(matches!(
            motion.advance(&[0.02], &[[0.03, 0.0, 0.0]], 0.01, cx),
            Err(DynamicsError::Model(Error::Material { .. }))
        ));
        let mut limited = config();
        limited.solve.max_operator_evaluations = 1;
        let limited_motion = motion
            .state()
            .clone()
            .into_dynamics(&mass, vec![0.01], 2.0, limited)
            .unwrap();
        assert!(matches!(
            limited_motion.advance(&[0.02], &[[0.0; 3]], 0.01, cx),
            Err(DynamicsError::Model(Error::Budget(_)))
        ));
        assert_eq!(motion.state().points(), saved);
        assert_eq!(motion.state().displacements_m(), &[0.0]);
        assert_eq!(motion.velocity_m_s(), &[0.01]);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.steps(), 0);
        let retry = motion.advance(&[0.02], &[[0.0; 3]], 0.01, cx).unwrap();
        assert_eq!(
            retry.motion.state().displacements_m(),
            expected.motion.state().displacements_m()
        );
        assert_eq!(retry.motion.velocity_m_s(), expected.motion.velocity_m_s());
        assert_eq!(
            retry.motion.state().points(),
            expected.motion.state().points()
        );
    });
}
