//! Stress and strength observations through accepted driven thermoelastic motion.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::TransientAdjointError3;
use fs_solid::loads3::transient::adjoint::support::NodalMotionCotangent3;
use fs_solid::loads3::transient::adjoint::support::stress::von_mises_limit;
use fs_solid::loads3::transient::adjoint::support::trajectory::{
    ElementStressObservation3, NodalMotionObservation3, PrescribedRecordingBudget3,
};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];
const FRAME: fs_blake3::ContentHash = fs_blake3::ContentHash([117; 32]);
const SUPPORT: [f64; 4] = [0.0001, 0.00012, 0.00020, 0.00022];
const FORCE: [f64; 4] = [0.0005, 0.002, -0.001, 0.0015];
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
                seed: 117,
                kernel_id: 117,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn seed() -> NodalMotionCotangent3 {
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

// Independent stiffness, density, Rayleigh eta/zeta, support and force amplitudes.
fn evaluate(p: [f64; 6], cx: &Cx<'_>) -> (f64, [f64; 6]) {
    let material = TetElasticMaterial::try_new(6. * p[1], 60. * p[0], 0.25, [117; 32]).unwrap();
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
        [0.0007, 0.0002, -0.0001, 0.0003, -0.0002, 0.0001],
        material.material_state_identity,
        FRAME,
        0.01,
    )
    .unwrap();
    let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
    let initial = geometry
        .assemble(&[], VolumeLoad3::None, &[(3, FORCE[0] * p[5])], cx)
        .unwrap();
    let mut q = ZERO;
    let mut v = ZERO;
    let mut a = ZERO;
    q[0][0] = SUPPORT[0] * p[4];
    q[1][0] = 0.00015;
    v[0][0] = -0.0005 * p[4];
    v[1][0] = 0.001;
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
    let mut stresses = vec![
        *recorded
            .run()
            .recover_response(FRAME, 0.05, cx)
            .unwrap()
            .stress
            .elements()[0]
            .stress_mandel_pa(),
    ];
    let mut motion = vec![recorded.run().nodal_motion(cx).unwrap()];
    for (i, h) in STEPS.into_iter().enumerate() {
        let load = geometry
            .assemble(&[], VolumeLoad3::None, &[(3, FORCE[i + 1] * p[5])], cx)
            .unwrap();
        recorded
            .advance_with_supports(&load, &[(0, SUPPORT[i + 1] * p[4])], h, cx)
            .unwrap();
        stresses.push(
            *recorded
                .run()
                .recover_response(FRAME, 0.05, cx)
                .unwrap()
                .stress
                .elements()[0]
                .stress_mandel_pa(),
        );
        motion.push(recorded.run().nodal_motion(cx).unwrap());
    }
    let target = [0.002, -0.003, 0.001, -0.001, 0.004, -0.002];
    let mut value = 0.;
    let mut stress_rows = Vec::new();
    let mut motion_rows = Vec::new();
    for (state, &sigma) in stresses.iter().enumerate() {
        let (strength, strength_seed) = von_mises_limit(sigma, 0.08).unwrap();
        value += 0.1 * strength;
        let values: [f64; 6] = core::array::from_fn(|c| {
            let scale = 0.03 + 0.005 * c as f64;
            let residual = (sigma[c] - target[c]) / scale;
            value += 0.5 * residual * residual;
            residual / scale + 0.1 * strength_seed[c]
        });
        stress_rows.push(ElementStressObservation3 {
            state,
            element: 0,
            values: values.map(|v| 0.25 * v),
        });
        stress_rows.push(ElementStressObservation3 {
            state,
            element: 0,
            values: values.map(|v| 0.75 * v),
        });
        let residual = (motion[state][0][1][0] - 0.0004) / 0.001;
        value += 0.5 * residual * residual;
        motion_rows.push(NodalMotionObservation3 {
            state,
            node: 1,
            q: [residual / 0.001, 0., 0.],
            v: [0.; 3],
            a: [0.; 3],
        });
    }
    stress_rows.reverse(); // Endpoint observations need not be chronological.
    let result = recorded
        .reverse_design_with_stress(&seed(), &motion_rows, &stress_rows, 0.05, cx)
        .unwrap();
    let slot = recorded.adjoint().fixed_dofs().binary_search(&0).unwrap();
    let support = result
        .controls
        .displacements
        .iter()
        .zip(SUPPORT)
        .map(|(r, s)| r.values[slot] * s)
        .sum::<f64>()
        - 0.0005 * result.controls.initial.v[0][0]
        + 0.002 * result.controls.initial.a[0][0];
    let force = result
        .controls
        .forces
        .iter()
        .zip(FORCE)
        .map(|(r, f)| r.values[0] * f)
        .sum();
    (
        value,
        [
            result.material.element_scales[0][0] / p[0],
            result.material.element_scales[0][1] / p[1],
            result.material.rayleigh[0],
            result.material.rayleigh[1],
            support,
            force,
        ],
    )
}

#[test]
fn g3_mixed_stress_strength_and_motion_objectives_match_six_design_directions() {
    with_cx(false, |cx| {
        for p in [
            [1., 1., 0.4, 0.006, 1., 1.],
            [1.6, 0.7, 0.8, 0.009, 0.8, 1.2],
        ] {
            let (_, gradient) = evaluate(p, cx);
            for axis in 0..6 {
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
fn g1_fully_prescribed_intermediate_stress_retains_direct_stiffness_and_endpoint_control() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0.25, [118; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
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
        recorded
            .advance_with_supports(&load, &[(3, 0.001)], 0.02, cx)
            .unwrap();
        let weights = [2., -1., 0.5, 0.3, -0.2, 0.1];
        let direct = recorded
            .run()
            .stress_pullback(&[(0, weights)], 0.05, Default::default(), cx)
            .unwrap();
        recorded
            .advance_with_supports(&load, &[(3, 0.002)], 0.03, cx)
            .unwrap();
        let row = ElementStressObservation3 {
            state: 1,
            element: 0,
            values: weights,
        };
        let result = recorded
            .reverse_design_with_stress(&seed(), &[], &[row], 0.05, cx)
            .unwrap();
        assert_eq!(result.material, direct.material);
        assert!(result.material.element_scales[0][0].abs() > 0.);
        assert_eq!(result.material.element_scales[0][1], 0.);
        assert_eq!(result.material.rayleigh, [0.; 2]);
        assert!(result.controls.forces.iter().all(|r| r.values.is_empty()));
        for dof in 0..12 {
            close(
                result.controls.displacements[1].values[dof],
                direct.motion.q[dof / 3][dof % 3],
            );
            close(result.controls.displacements[0].values[dof], 0.);
            close(result.controls.displacements[2].values[dof], 0.);
        }
    });
}

#[test]
fn g4_stress_observations_refuse_invalid_inputs_and_leave_the_tape_reusable() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0.25, [119; 32]).unwrap();
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
        let mut q = ZERO;
        q[1][0] = 0.001;
        let run = load
            .begin_transient_with_supports(
                &problem,
                None,
                0.,
                &q,
                &ZERO,
                &ZERO,
                config(0.4, 0.006),
                cx,
            )
            .unwrap();
        let mut recorded = run
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
        let row = ElementStressObservation3 {
            state: 0,
            element: 0,
            values: [1., 0., 0., 0., 0., 0.],
        };
        assert!(
            recorded
                .reverse_design_with_stress(&seed(), &[], &[row], 0.05, cx)
                .is_err()
        );
        recorded
            .enable_material_recording(Default::default(), cx)
            .unwrap();
        let before = recorded.run().state().clone();
        for ceiling in [0., -1., f64::NAN, f64::INFINITY, 1e-5] {
            assert!(
                recorded
                    .reverse_design_with_stress(&seed(), &[], &[row], ceiling, cx)
                    .is_err()
            );
        }
        let zero = ElementStressObservation3 {
            values: [0.; 6],
            ..row
        };
        assert!(
            recorded
                .reverse_design_with_stress(&seed(), &[], &[zero], 1e-5, cx)
                .is_err()
        );
        for bad in [
            ElementStressObservation3 { state: 1, ..row },
            ElementStressObservation3 { element: 1, ..row },
            ElementStressObservation3 {
                values: [f64::NAN; 6],
                ..row
            },
        ] {
            assert!(
                recorded
                    .reverse_design_with_stress(&seed(), &[], &[bad], 0.05, cx)
                    .is_err()
            );
        }
        assert!(matches!(
            recorded.reverse_design_with_stress(&seed(), &[], &[row, row], 0.05, cx),
            Err(TransientAdjointError3::Budget)
        ));
        let motion = NodalMotionObservation3 {
            state: 0,
            node: 0,
            q: [0.; 3],
            v: [0.; 3],
            a: [0.; 3],
        };
        assert!(matches!(
            recorded.reverse_design_with_stress(&seed(), &[motion], &[row], 0.05, cx),
            Err(TransientAdjointError3::Budget)
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                recorded.reverse_design_with_stress(&seed(), &[], &[row], 0.05, cancelled),
                Err(TransientAdjointError3::Cancelled)
            ))
        });
        let result = recorded
            .reverse_design_with_stress(&seed(), &[], &[row], 0.05, cx)
            .unwrap();
        assert_eq!(
            result,
            recorded
                .clone()
                .reverse_design_with_stress(&seed(), &[], &[row], 0.05, cx)
                .unwrap()
        );
        assert_eq!(recorded.run().state(), &before);
        assert_eq!(
            recorded
                .reverse_design_with_stress(&seed(), &[motion], &[], 0.05, cx)
                .unwrap(),
            recorded
                .reverse_design_with_observations(&seed(), &[motion], cx)
                .unwrap()
        );
    });
}
