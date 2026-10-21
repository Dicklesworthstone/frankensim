//! Anchor/actuator force objectives through the original full transient tape.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::TransientAdjointError3;
use fs_solid::loads3::transient::adjoint::support::NodalMotionCotangent3;
use fs_solid::loads3::transient::adjoint::support::trajectory::{
    NodalMotionObservation3, PrescribedRecordingBudget3, ReactionObservation3,
};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];
const SUPPORT: [f64; 4] = [0.0001, 0.00012, 0.00020, 0.00022];
const FORCE: [f64; 4] = [0.0005, 0.002, -0.001, 0.0015];
const FIXED_FORCE: [f64; 4] = [0.0003, -0.0002, 0.0001, -0.0004];
const STEPS: [f64; 3] = [0.025, 0.04, 0.02];

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 107,
                kernel_id: 107,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn zero_seed() -> NodalMotionCotangent3 {
    NodalMotionCotangent3 {
        q: ZERO.to_vec(),
        v: ZERO.to_vec(),
        a: ZERO.to_vec(),
    }
}
fn config(eta: f64, zeta: f64) -> TransientConfig3 {
    TransientConfig3 {
        mass_damping_per_s: eta,
        stiffness_damping_s: zeta,
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-12,
            maximum_iterations: 10_000,
        },
        ..Default::default()
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-6 * (1. + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

// Stiffness, density, eta, zeta, support amplitude, free force, fixed force.
fn evaluate(p: [f64; 7], cx: &Cx<'_>) -> (f64, [f64; 7]) {
    let material = TetElasticMaterial::try_new(6. * p[1], 60. * p[0], 0.25, [107; 32]).unwrap();
    let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    let thermal = TetThermalStrainState::try_new(
        301.,
        300.,
        [0.0007, 0., 0., 0., 0., 0.],
        material.material_state_identity,
        fs_blake3::ContentHash([107; 32]),
        0.01,
    )
    .unwrap();
    let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
    let initial = geometry
        .assemble(
            &[],
            VolumeLoad3::None,
            &[(0, FIXED_FORCE[0] * p[6]), (3, FORCE[0] * p[5])],
            cx,
        )
        .unwrap();
    let mut q = ZERO;
    q[0][0] = SUPPORT[0] * p[4];
    q[1][0] = 0.00015;
    let mut v = ZERO;
    v[0][0] = -0.0005 * p[4];
    v[1][0] = 0.001;
    let mut a = ZERO;
    a[0][0] = 0.002 * p[4];
    let run = initial
        .begin_transient_with_supports(
            &problem,
            Some(TetThermalStrainField::Uniform(&thermal)),
            0.,
            &q,
            &v,
            &a,
            config(p[2], p[3]),
            cx,
        )
        .unwrap();
    let mut recorded = run
        .record_prescribed_controls(
            Default::default(),
            Default::default(),
            Default::default(),
            cx,
        )
        .unwrap();
    recorded
        .enable_material_recording(Default::default(), cx)
        .unwrap();
    let mut states = vec![recorded.run().nodal_motion(cx).unwrap()];
    let mut reactions = vec![
        recorded
            .run()
            .recover_response(fs_blake3::ContentHash([107; 32]), 0.05, cx)
            .unwrap()
            .support_reactions_n,
    ];
    for (i, h) in STEPS.into_iter().enumerate() {
        let load = geometry
            .assemble(
                &[],
                VolumeLoad3::None,
                &[(0, FIXED_FORCE[i + 1] * p[6]), (3, FORCE[i + 1] * p[5])],
                cx,
            )
            .unwrap();
        recorded
            .advance_with_supports(&load, &[(0, SUPPORT[i + 1] * p[4])], h, cx)
            .unwrap();
        states.push(recorded.run().nodal_motion(cx).unwrap());
        reactions.push(
            recorded
                .run()
                .recover_response(fs_blake3::ContentHash([107; 32]), 0.05, cx)
                .unwrap()
                .support_reactions_n,
        );
    }
    let mut value = 0.;
    let mut motion_rows = Vec::new();
    let mut reaction_rows = Vec::new();
    for (state, motion) in states.iter().enumerate() {
        let residual = (motion[0][1][0] - 0.0004) / 0.001;
        value += 0.5 * residual * residual;
        motion_rows.push(NodalMotionObservation3 {
            state,
            node: 1,
            q: [residual / 0.001, 0., 0.],
            v: [0.; 3],
            a: [0.; 3],
        });
        for (dof, target, scale) in [(0, 0.004, 0.02), (6, -0.003, 0.03)] {
            let residual = (reactions[state][dof / 3][dof % 3] - target) / scale;
            value += 0.5 * residual * residual;
            // Duplicate sparse rows exercise deterministic accumulation.
            reaction_rows.push(ReactionObservation3 {
                state,
                dof,
                value: residual / scale * 0.25,
            });
            reaction_rows.push(ReactionObservation3 {
                state,
                dof,
                value: residual / scale * 0.75,
            });
        }
    }
    let result = recorded
        .reverse_design_with_reactions(&zero_seed(), &motion_rows, &reaction_rows, cx)
        .unwrap();
    let support_index = recorded.adjoint().fixed_dofs().binary_search(&0).unwrap();
    let support_gradient = result
        .design
        .controls
        .displacements
        .iter()
        .zip(SUPPORT)
        .map(|(row, s)| row.values[support_index] * s)
        .sum::<f64>()
        - 0.0005 * result.design.controls.initial.v[0][0]
        + 0.002 * result.design.controls.initial.a[0][0];
    let force_gradient = result
        .design
        .controls
        .forces
        .iter()
        .zip(FORCE)
        .map(|(r, f)| r.values[0] * f)
        .sum();
    let fixed_gradient = result
        .fixed_forces
        .iter()
        .filter(|r| r.dof == 0)
        .map(|r| r.value * FIXED_FORCE[r.state])
        .sum();
    assert_eq!(result.fixed_forces.len(), reaction_rows.len());
    for (actual, input) in result.fixed_forces.iter().zip(&reaction_rows) {
        assert_eq!(
            (actual.state, actual.dof, actual.value),
            (input.state, input.dof, -input.value)
        );
    }
    close(
        result.thermal_fixed_force[support_index],
        result
            .fixed_forces
            .iter()
            .filter(|r| r.dof == 0)
            .map(|r| r.value)
            .sum(),
    );
    (
        value,
        [
            result.design.material.element_scales[0][0] / p[0],
            result.design.material.element_scales[0][1] / p[1],
            result.design.material.rayleigh[0],
            result.design.material.rayleigh[1],
            support_gradient,
            force_gradient,
            fixed_gradient,
        ],
    )
}

#[test]
fn g3_reaction_and_motion_objectives_match_seven_physical_design_perturbations() {
    with_cx(false, |cx| {
        for p in [
            [1., 1., 0.4, 0.006, 1., 1., 1.],
            [1.6, 0.7, 0.8, 0.009, 0.8, 1.2, 0.6],
        ] {
            let (_, gradient) = evaluate(p, cx);
            for axis in 0..7 {
                let h = 1e-5 * p[axis];
                let mut plus = p;
                let mut minus = p;
                plus[axis] += h;
                minus[axis] -= h;
                close(
                    gradient[axis],
                    (evaluate(plus, cx).0 - evaluate(minus, cx).0) / (2. * h),
                );
            }
        }
    });
}

#[test]
fn g3_fully_prescribed_reactions_retain_direct_material_and_initial_motion_derivatives() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0.25, [108; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry
            .assemble(&[], VolumeLoad3::None, &[(3, 0.001)], cx)
            .unwrap();
        let mut q = ZERO;
        let mut v = ZERO;
        let mut a = ZERO;
        q[1][0] = 0.001;
        v[1][0] = 0.02;
        a[1][0] = 0.03;
        let run = load
            .begin_transient_with_supports(&problem, None, 0., &q, &v, &a, config(0.4, 0.006), cx)
            .unwrap();
        let direct = run
            .reaction_pullback(&[(3, 2.)], Default::default(), cx)
            .unwrap();
        let mut recorded = run
            .record_prescribed_controls(
                Default::default(),
                Default::default(),
                Default::default(),
                cx,
            )
            .unwrap();
        recorded
            .enable_material_recording(Default::default(), cx)
            .unwrap();
        let row = ReactionObservation3 {
            state: 0,
            dof: 3,
            value: 2.,
        };
        let result = recorded
            .reverse_design_with_reactions(&zero_seed(), &[], &[row], cx)
            .unwrap();
        assert_eq!(result.design.material, direct.material);
        assert!(
            result.design.material.element_scales[0]
                .iter()
                .all(|v| *v != 0.)
        );
        assert!(result.design.material.rayleigh.iter().all(|v| *v != 0.));
        assert!(result.design.controls.forces[0].values.is_empty());
        assert_eq!(result.design.controls.initial.v, direct.motion.v);
        assert_eq!(result.design.controls.initial.a, direct.motion.a);
        assert_eq!(result.design.controls.initial.q, ZERO);
        for dof in 0..12 {
            close(
                result.design.controls.displacements[0].values[dof],
                direct.motion.q[dof / 3][dof % 3],
            );
        }
        assert_eq!(result.fixed_forces[0].value, -2.);
        assert_eq!(result.thermal_fixed_force[3], -2.);
        let empty = recorded
            .reverse_design_with_reactions(&zero_seed(), &[], &[], cx)
            .unwrap();
        assert_eq!(
            empty.design,
            recorded.reverse_design(&zero_seed(), cx).unwrap()
        );
    });
}

#[test]
fn g4_reaction_observations_are_bounded_validated_and_cancellation_leaves_tape_reusable() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0.25, [109; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: Default::default(),
        };
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let run = load
            .begin_transient_with_supports(
                &problem,
                None,
                0.,
                &ZERO,
                &ZERO,
                &ZERO,
                config(0.4, 0.006),
                cx,
            )
            .unwrap();
        let row = ReactionObservation3 {
            state: 0,
            dof: 0,
            value: 1.,
        };
        let mut recorded = run
            .clone()
            .record_prescribed_controls(
                PrescribedRecordingBudget3 {
                    maximum_observations: 1,
                    ..Default::default()
                },
                Default::default(),
                Default::default(),
                cx,
            )
            .unwrap();
        assert!(
            recorded
                .reverse_design_with_reactions(&zero_seed(), &[], &[row], cx)
                .is_err()
        );
        recorded
            .enable_material_recording(Default::default(), cx)
            .unwrap();
        let before = recorded.run().state().clone();
        for bad in [
            ReactionObservation3 { state: 1, ..row },
            ReactionObservation3 { dof: 3, ..row },
            ReactionObservation3 { dof: 12, ..row },
            ReactionObservation3 {
                value: f64::NAN,
                ..row
            },
        ] {
            assert!(
                recorded
                    .reverse_design_with_reactions(&zero_seed(), &[], &[bad], cx)
                    .is_err()
            );
        }
        assert!(matches!(
            recorded.reverse_design_with_reactions(&zero_seed(), &[], &[row, row], cx),
            Err(TransientAdjointError3::Budget)
        ));
        let motion = NodalMotionObservation3 {
            state: 0,
            node: 1,
            q: [1., 0., 0.],
            v: [0.; 3],
            a: [0.; 3],
        };
        assert!(matches!(
            recorded.reverse_design_with_reactions(&zero_seed(), &[motion], &[row], cx),
            Err(TransientAdjointError3::Budget)
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                recorded.reverse_design_with_reactions(&zero_seed(), &[], &[row], cancelled),
                Err(TransientAdjointError3::Cancelled)
            ))
        });
        let good = recorded
            .reverse_design_with_reactions(&zero_seed(), &[], &[row], cx)
            .unwrap();
        assert_eq!(
            good,
            recorded
                .clone()
                .reverse_design_with_reactions(&zero_seed(), &[], &[row], cx)
                .unwrap()
        );
        assert_eq!(recorded.run().state(), &before);
        // Base controls use 49 scalars; 11 fixed thermal values plus one force
        // observation need 61. Recording still fits the smaller allowance.
        let mut tight = run
            .record_prescribed_controls(
                PrescribedRecordingBudget3 {
                    maximum_gradient_values: 60,
                    ..Default::default()
                },
                Default::default(),
                Default::default(),
                cx,
            )
            .unwrap();
        tight
            .enable_material_recording(Default::default(), cx)
            .unwrap();
        assert!(matches!(
            tight.reverse_design_with_reactions(&zero_seed(), &[], &[row], cx),
            Err(TransientAdjointError3::Budget)
        ));
    });
}
