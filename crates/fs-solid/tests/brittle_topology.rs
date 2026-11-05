//! Explicit original-mesh seams feed the native consistent-mass cohesive body
//! dynamics. Manufactured coefficients verify mechanics, not material validity.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget as ExecBudget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::cohesive::BilinearCohesive;
use fs_solid::brittle::topology::{FractureTopologyBudget, InteriorSeam, TetFractureTopology};
use fs_solid::cohesive::dynamics::DynamicsConfig;
use fs_solid::cohesive::p1::{P1Interfaces, TraceBudget};
use fs_solid::cohesive::{Budget, CohesiveSystem, Config};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
};
use fs_solid::loads3::TetLoadGeometry3;
use fs_solver::{CsrOp, LinearOp, NewtonKrylovConfig};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 73,
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
            maximum_displacement_m: 0.1,
        },
        difference_step_m: 1e-8,
        maximum_step_s: 0.02,
        maximum_velocity_m_s: 2.0,
        maximum_steps: 100,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn original_tetrahedra_split_fracture_separate_and_reclose_without_losing_inertia() {
    with_cx(|cx| {
        let nodes = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., -1.],
            [0., 0., 1.],
        ];
        let tets = [[0, 2, 1, 3], [0, 1, 2, 4]];
        let geometry =
            TetLoadGeometry3::new(&nodes, &tets, TetAssemblyBudget::default(), 8, cx).unwrap();
        let topology = TetFractureTopology::new(
            &geometry,
            &[InteriorSeam {
                id: 73,
                nodes: [2, 0, 1],
                side_a_element: 0,
            }],
            FractureTopologyBudget {
                maximum_output_nodes: 8,
                maximum_face_entries: 8,
                maximum_seams: 1,
            },
            cx,
        )
        .unwrap();
        assert_eq!(topology.element_fragments(), &[0, 1]);
        let velocity: Vec<_> = (0..5).flat_map(|_| [0.1, -0.05, 0.02]).collect();
        let transfer = topology
            .transfer_motion(&[0.0; 15], &velocity, &[6.0; 2], cx)
            .unwrap();
        assert_eq!(transfer.before, transfer.after);

        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [73; 32]).unwrap();
        let assembled = TetLinearElasticProblem {
            nodes_m: topology.nodes_m(),
            tetrahedra: topology.tetrahedra(),
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::default(),
        }
        .assemble(cx)
        .unwrap();
        close(assembled.total_mass_kg, transfer.after.mass_kg, 1e-14);
        let law = BilinearCohesive::new(1200.0, 0.12, 0.00006, 2400.0, 0.02).unwrap();
        let faces = P1Interfaces::new(
            topology.nodes_m(),
            topology.tetrahedra(),
            &assembled.free_dofs,
            &[topology.seams()[0].face_pair(&law)],
            TraceBudget {
                maximum_nodes: 8,
                maximum_tetrahedra: 2,
                maximum_interfaces: 1,
                minimum_face_quality: 1e-10,
            },
            cx,
        )
        .unwrap();
        close(
            faces.points().iter().map(|p| p.area_weight_m2).sum(),
            0.5,
            1e-15,
        );
        assert!(
            faces
                .jumps(&transfer.displacement_m, &mut || false)
                .unwrap()
                .iter()
                .all(|jump| *jump == [0.0; 3])
        );
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let system = CohesiveSystem::new(
            &elastic,
            faces.patches().to_vec(),
            Budget {
                max_dofs: 24,
                max_patches: 3,
                max_jump_terms: 54,
                max_krylov_basis_entries: 50_000,
            },
        )
        .unwrap();
        let mut mv = vec![0.0; 24];
        mass.apply(&transfer.velocity_m_s, &mut mv);
        close(
            0.5 * mv
                .iter()
                .zip(&transfer.velocity_m_s)
                .map(|(p, v)| p * v)
                .sum::<f64>(),
            transfer.after.kinetic_energy_j,
            1e-14,
        );
        let initial_energy = transfer.after.kinetic_energy_j;
        let initial_momentum = transfer.after.linear_momentum_kg_m_s;
        let mut motion = system
            .bonded_state()
            .into_dynamics(&mass, transfer.velocity_m_s, 0.0, config())
            .unwrap();

        // Actual equal/opposite body forces create separation from a continuous
        // original velocity field. No artificial fragment velocity is assigned.
        let mut acceleration = vec![0.0; 24];
        for (element, tet) in topology.tetrahedra().iter().enumerate() {
            for &node in tet {
                acceleration[3 * node + 2] = if element == 0 { -50.0 } else { 50.0 };
            }
        }
        let mut pull = vec![0.0; 24];
        mass.apply(&acceleration, &mut pull);
        let push: Vec<_> = pull.iter().map(|f| -f).collect();
        let mut work = 0.0;
        let mut fracture = 0.0;
        let mut compression_seen = false;
        for step_index in 0..56 {
            let force = if step_index < 12 {
                &pull[..]
            } else if step_index < 20 {
                &[0.0; 24][..]
            } else {
                &push[..]
            };
            let step = motion.advance(force, &[[0.0; 3]; 3], 0.001, cx).unwrap();
            work += step.external_work_j;
            fracture += step.cohesive_energy_j.dissipated;
            close(
                step.kinetic_after_j
                    + step.bulk_stored_after_j
                    + step.cohesive_stored_after_j
                    + fracture,
                initial_energy + work,
                2e-8,
            );
            close(step.energy_residual_j, 0.0, 2e-9);
            assert_eq!(step.viscous_dissipated_j, 0.0);
            assert_eq!(step.friction_dissipated_j, 0.0);
            mass.apply(step.motion.velocity_m_s(), &mut mv);
            for c in 0..3 {
                close(
                    mv.chunks_exact(3).map(|p| p[c]).sum(),
                    initial_momentum[c],
                    2e-8,
                );
            }
            if step_index >= 12 {
                assert!(step.motion.state().points().iter().all(|p| p.failed()));
                close(step.cohesive_energy_j.dissipated, 0.0, 1e-15);
            }
            if (12..20).contains(&step_index) {
                assert!(
                    step.mean_traction_pa
                        .iter()
                        .all(|traction| *traction == [0.0; 3])
                );
            }
            if step_index >= 20 {
                compression_seen |= step
                    .motion
                    .state()
                    .points()
                    .iter()
                    .any(|p| p.jump_m()[0] < 0.0 && p.traction_pa()[0] < 0.0);
            }
            motion = step.motion;
        }
        close(fracture, 0.5 * law.fracture_energy_j_m2(), 1e-11);
        assert!(
            compression_seen,
            "the same failed interface must resist reclosure"
        );
        assert!(motion.state().points().iter().all(|p| p.failed()));
        assert!(motion.velocity_m_s().iter().any(|v| v.abs() > 0.2));
    });
}
