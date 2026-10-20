//! Full driven-support design gradients: material, load and support controls.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::support::trajectory::{
    NodalMotionObservation3, PrescribedMaterialRecordingBudget3, PrescribedRecordingBudget3,
};
use fs_solid::loads3::transient::adjoint::support::{NodalMotionCotangent3, SupportAdjointBudget3};
use fs_solid::loads3::transient::adjoint::{TransientAdjointConfig3, TransientAdjointError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];
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
                seed: 97,
                kernel_id: 97,
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
        (a - b).abs() <= 3e-6 * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

// Independent stiffness, density, eta, zeta, support amplitude, force amplitude.
fn evaluate(p: [f64; 6], cx: &Cx<'_>) -> (f64, [f64; 6]) {
    let material = TetElasticMaterial::try_new(6.0 * p[1], 60.0 * p[0], 0.25, [97; 32]).unwrap();
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
        fs_blake3::ContentHash([97; 32]),
        0.01,
    )
    .unwrap();
    let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
    let initial = geometry
        .assemble(&[], VolumeLoad3::None, &[(3, FORCE[0] * p[5])], cx)
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
            PrescribedRecordingBudget3::default(),
            TransientAdjointConfig3::default(),
            SupportAdjointBudget3::default(),
            cx,
        )
        .unwrap();
    recorded
        .enable_material_recording(Default::default(), cx)
        .unwrap();
    let mut states = vec![recorded.run().nodal_motion(cx).unwrap()];
    for (i, h) in STEPS.into_iter().enumerate() {
        let load = geometry
            .assemble(&[], VolumeLoad3::None, &[(3, FORCE[i + 1] * p[5])], cx)
            .unwrap();
        recorded
            .advance_with_supports(&load, &[(0, SUPPORT[i + 1] * p[4])], h, cx)
            .unwrap();
        states.push(recorded.run().nodal_motion(cx).unwrap());
    }
    let mut value = 0.;
    let mut observations = Vec::new();
    for (state, motion) in states.iter().enumerate() {
        let mut row = NodalMotionObservation3 {
            state,
            node: 1,
            q: [0.; 3],
            v: [0.; 3],
            a: [0.; 3],
        };
        for (field, target, scale) in [(0, 0.0004, 0.001), (1, -0.002, 0.02), (2, 0.004, 0.3)] {
            let residual = (motion[field][1][0] - target) / scale;
            value += 0.5 * residual * residual;
            match field {
                0 => row.q[0] = residual / scale,
                1 => row.v[0] = residual / scale,
                _ => row.a[0] = residual / scale,
            }
        }
        observations.push(row);
        let residual = (motion[0][0][0] - 0.0003) / 0.001;
        value += 0.5 * residual * residual;
        observations.push(NodalMotionObservation3 {
            state,
            node: 0,
            q: [residual / 0.001, 0., 0.],
            v: [0.; 3],
            a: [0.; 3],
        });
    }
    let result = recorded
        .reverse_design_with_observations(&zero_seed(), &observations, cx)
        .unwrap();
    assert_eq!(
        result.controls,
        recorded
            .reverse_with_observations(&zero_seed(), &observations, cx)
            .unwrap()
    );
    let support_index = recorded.adjoint().fixed_dofs().binary_search(&0).unwrap();
    let support_gradient = result
        .controls
        .displacements
        .iter()
        .zip(SUPPORT)
        .map(|(r, s)| r.values[support_index] * s)
        .sum::<f64>()
        - 0.0005 * result.controls.initial.v[0][0]
        + 0.002 * result.controls.initial.a[0][0];
    let force_gradient = result
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
            support_gradient,
            force_gradient,
        ],
    )
}

#[test]
fn g3_one_reverse_matches_material_and_authored_control_changes_across_a_full_trajectory() {
    with_cx(false, |cx| {
        for p in [
            [1., 1., 0.4, 0.006, 1., 1.],
            [1.6, 0.7, 0.8, 0.009, 0.8, 1.2],
        ] {
            let (_, gradient) = evaluate(p, cx);
            for axis in 0..6 {
                let step = 1e-5 * p[axis];
                let mut plus = p;
                let mut minus = p;
                plus[axis] += step;
                minus[axis] -= step;
                let difference = (evaluate(plus, cx).0 - evaluate(minus, cx).0) / (2. * step);
                close(gradient[axis], difference);
            }
        }
    });
}

#[test]
fn g4_material_recording_is_opt_in_bounded_atomic_and_forkable() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [98; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
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
        let mut unrecorded = recorded.clone();
        assert!(unrecorded.reverse_design(&zero_seed(), cx).is_err());
        unrecorded
            .advance_with_supports(&load, &[], 0.01, cx)
            .unwrap();
        assert!(
            unrecorded
                .enable_material_recording(Default::default(), cx)
                .is_err()
        );
        let before = recorded.run().state().clone();
        assert!(matches!(
            recorded.enable_material_recording(
                PrescribedMaterialRecordingBudget3 {
                    maximum_motion_values: 35,
                    ..Default::default()
                },
                cx
            ),
            Err(TransientAdjointError3::Budget)
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                recorded.enable_material_recording(Default::default(), cancelled),
                Err(TransientAdjointError3::Cancelled)
            ))
        });
        assert_eq!(recorded.run().state(), &before);
        recorded
            .enable_material_recording(
                PrescribedMaterialRecordingBudget3 {
                    maximum_motion_values: 72,
                    ..Default::default()
                },
                cx,
            )
            .unwrap();
        assert!(
            recorded
                .enable_material_recording(Default::default(), cx)
                .is_err()
        );
        assert!(
            recorded
                .advance_with_supports(&load, &[(0, 0.001), (0, 0.002)], 0.01, cx)
                .is_err()
        );
        assert_eq!(recorded.sample_count(), 1);
        let mut fork = recorded.clone();
        recorded
            .advance_with_supports(&load, &[(0, 0.0001)], 0.01, cx)
            .unwrap();
        fork.advance_with_supports(&load, &[(0, 0.0001)], 0.01, cx)
            .unwrap();
        let mut seed = zero_seed();
        seed.a[1][0] = 1.;
        let expected = fork.reverse_design(&seed, cx).unwrap();
        assert_eq!(recorded.reverse_design(&seed, cx).unwrap(), expected);
        let accepted = recorded.run().nodal_motion(cx).unwrap();
        assert!(matches!(
            recorded.advance_with_supports(&load, &[], 0.01, cx),
            Err(TransientAdjointError3::Budget)
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                recorded.reverse_design(&seed, cancelled),
                Err(TransientAdjointError3::Cancelled)
            ))
        });
        assert_eq!(recorded.run().nodal_motion(cx).unwrap(), accepted);
        assert_eq!(recorded.sample_count(), 2);
        assert_eq!(recorded.reverse_design(&seed, cx).unwrap(), expected);
    });
}

#[test]
fn g1_fully_prescribed_motion_has_zero_material_derivative() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0.25, [99; 32]).unwrap();
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
        let mut seed = zero_seed();
        seed.q[0][0] = 1.;
        seed.v[1][1] = -0.2;
        seed.a[2][2] = 0.03;
        for _ in 0..2 {
            recorded
                .advance_with_supports(&load, &[(0, 0.0001)], 0.01, cx)
                .unwrap();
        }
        let result = recorded.reverse_design(&seed, cx).unwrap();
        assert_eq!(result.material.element_scales, vec![[0.; 2]]);
        assert_eq!(result.material.rayleigh, [0.; 2]);
        assert_eq!(result.controls.total_adjoint_iterations, 0);
    });
}
