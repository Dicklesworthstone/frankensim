//! G1/G3/G4: physical acceleration histories through the original trajectory adjoint.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::{FrameTree, WORLD};
use fs_solid::linear3::*;
use fs_solid::loads3::history::dynamics::events::EventLoadHistory3;
use fs_solid::loads3::history::supports::SupportAcceleration3;
use fs_solid::loads3::history::{ACCELERATION_DIMS, TetLoadHistory3, VolumeHistory3};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::support::trajectory::{
    ElementStressObservation3, NodalMotionObservation3, PrescribedMaterialRecordingBudget3,
    PrescribedRecordingBudget3, ReactionObservation3, RecordedPrescribedTransient3,
};
use fs_solid::loads3::transient::adjoint::support::{NodalMotionCotangent3, SupportAdjointBudget3};
use fs_solid::loads3::transient::adjoint::trajectory::ControlForceKind3;
use fs_solid::loads3::transient::adjoint::{TransientAdjointConfig3, TransientAdjointError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];
const ALL: [usize; 4] = [0, 1, 2, 3];
const FIXED: [usize; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
const FRAME: ContentHash = ContentHash([121; 32]);
type FullMotion = [Vec<[f64; 3]>; 3];

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
                seed: 121,
                kernel_id: 121,
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
fn problem<'a>(
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::default(),
    }
}
fn config(eta: f64, zeta: f64) -> TransientConfig3 {
    TransientConfig3 {
        mass_damping_per_s: eta,
        stiffness_damping_s: zeta,
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-12,
            maximum_iterations: 10000,
        },
        ..Default::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1. + actual.abs().max(expected.abs())),
        "{actual:e} != {expected:e}"
    );
}
fn capture(
    run: &RecordedPrescribedTransient3<'_, '_, '_, '_>,
    cx: &Cx<'_>,
) -> (FullMotion, f64, [f64; 6]) {
    let response = run.run().recover_response(FRAME, 0.05, cx).unwrap();
    (
        run.run().nodal_motion(cx).unwrap(),
        response.support_reactions_n[0][0],
        *response.stress.elements()[0].stress_mandel_pa(),
    )
}

// Five load samples; four interval/jump support controls; initial free q/v and
// support q/v/a; physical stiffness, density and both Rayleigh coefficients.
fn mixed(p: [f64; 18], stress_objective: bool, cx: &Cx<'_>) -> (f64, [f64; 18]) {
    let material = TetElasticMaterial::try_new(6. * p[15], 60. * p[14], 0.25, [121; 32]).unwrap();
    let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
    let problem = problem(&material, &fixed);
    let thermal = TetThermalStrainState::try_new(
        301.,
        300.,
        [0.0007, 0., 0., 0., 0., 0.],
        material.material_state_identity,
        FRAME,
        0.01,
    )
    .unwrap();
    let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
    let load = |value| {
        geometry
            .assemble(&[], VolumeLoad3::None, &[(0, 0.4 * value), (3, value)], cx)
            .unwrap()
    };
    let (mut q, mut v, mut a) = (ZERO, ZERO, ZERO);
    q[1][0] = p[9];
    v[1][0] = p[10];
    q[0][0] = p[11];
    v[0][0] = p[12];
    a[0][0] = p[13];
    let initial = load(p[0]);
    let run = initial
        .begin_transient_with_supports(
            &problem,
            Some(TetThermalStrainField::Uniform(&thermal)),
            0.,
            &q,
            &v,
            &a,
            config(p[16], p[17]),
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
    let mut states = vec![capture(&recorded, cx)];
    recorded
        .advance_with_support_accelerations(&load(p[1]), &[(0, p[5])], 0.025, cx)
        .unwrap();
    states.push(capture(&recorded, cx));
    recorded
        .apply_support_acceleration_jump(&load(p[2]), &[(0, p[6])], cx)
        .unwrap();
    states.push(capture(&recorded, cx));
    assert_eq!(states[1].0[0], states[2].0[0]);
    assert_eq!(states[1].0[1], states[2].0[1]);
    recorded
        .advance_with_supports(&load(p[3]), &[(0, p[7])], 0.04, cx)
        .unwrap();
    states.push(capture(&recorded, cx));
    recorded
        .advance_with_support_accelerations(&load(p[4]), &[(0, p[8])], 0.02, cx)
        .unwrap();
    states.push(capture(&recorded, cx));
    assert_eq!(recorded.sample_count(), 5);
    assert_eq!(recorded.run().state().step, 3);
    let times: Vec<_> = recorded.accepted_times().collect();
    assert_eq!(times[1], times[2]);
    let mut value = 0.;
    let mut observations = Vec::new();
    for (state, (motion, _, _)) in states.iter().enumerate() {
        for (node, channel, target, scale) in [
            (1, 0, 0.0003, 0.001),
            (1, 1, 0.002, 0.02),
            (1, 2, 0.005, 0.3),
            (0, 0, 0.0002, 0.001),
            (0, 2, 0.01, 0.3),
        ] {
            let r = (motion[channel][node][0] - target) / scale;
            value += 0.5 * r * r;
            let mut row = NodalMotionObservation3 {
                state,
                node,
                q: [0.; 3],
                v: [0.; 3],
                a: [0.; 3],
            };
            match channel {
                0 => row.q[0] = r / scale,
                1 => row.v[0] = r / scale,
                _ => row.a[0] = r / scale,
            }
            observations.push(row);
        }
    }
    observations.reverse();
    let (design, fixed_forces) = if stress_objective {
        let mut rows = Vec::new();
        for state in [0, 1, 2, 4] {
            let mut weights = [0.; 6];
            for (channel, target) in [(0, 0.002), (3, -0.001)] {
                let r = (states[state].2[channel] - target) / 0.05;
                value += 0.5 * r * r;
                weights[channel] = r / 0.05;
            }
            rows.push(ElementStressObservation3 {
                state,
                element: 0,
                values: weights,
            });
        }
        (
            recorded
                .reverse_design_with_stress(&zero_seed(), &observations, &rows, 0.05, cx)
                .unwrap(),
            Vec::new(),
        )
    } else {
        let mut rows = Vec::new();
        for state in [0, 1, 2, 4] {
            let r = (states[state].1 - 0.004) / 0.02;
            value += 0.5 * r * r;
            rows.push(ReactionObservation3 {
                state,
                dof: 0,
                value: r / 0.02,
            });
        }
        let result = recorded
            .reverse_design_with_reactions(&zero_seed(), &observations, &rows, cx)
            .unwrap();
        (result.design, result.fixed_forces)
    };
    let c = &design.controls;
    assert_eq!(c.forces.len(), 5);
    assert_eq!(c.displacements.len(), 2);
    assert_eq!(c.accelerations.len(), 3);
    assert_eq!(c.accelerations[1].kind, ControlForceKind3::RightLimit);
    assert_eq!(c.forces[2].kind, ControlForceKind3::RightLimit);
    let mut gradient = [0.; 18];
    for (i, force) in c.forces.iter().enumerate() {
        gradient[i] = force.values[0];
    }
    for row in fixed_forces {
        gradient[row.state] += 0.4 * row.value;
    }
    gradient[5] = c.accelerations[0].values[0];
    gradient[6] = c.accelerations[1].values[0];
    gradient[7] = c.displacements[1].values[0];
    gradient[8] = c.accelerations[2].values[0];
    gradient[9] = c.initial.q[1][0];
    gradient[10] = c.initial.v[1][0];
    gradient[11] = c.displacements[0].values[0];
    gradient[12] = c.initial.v[0][0];
    gradient[13] = c.initial.a[0][0];
    gradient[14] = design.material.element_scales[0][0] / p[14];
    gradient[15] = design.material.element_scales[0][1] / p[15];
    gradient[16] = design.material.rayleigh[0];
    gradient[17] = design.material.rayleigh[1];
    assert_eq!(c.initial.a[1][0], 0.);
    (value, gradient)
}

fn check_mixed(stress_objective: bool, cx: &Cx<'_>) {
    let p = [
        0.0005, 0.002, -0.001, 0.0015, 0.003, 0.03, -0.012, 0.0002, 0.02, 0.00015, 0.001, 0.0001,
        -0.0005, 0.002, 1.1, 0.9, 0.4, 0.006,
    ];
    let (_, gradient) = mixed(p, stress_objective, cx);
    for i in 0..p.len() {
        let epsilon = 1e-6 * p[i].abs().max(0.01);
        let (mut plus, mut minus) = (p, p);
        plus[i] += epsilon;
        minus[i] -= epsilon;
        let numerical = (mixed(plus, stress_objective, cx).0
            - mixed(minus, stress_objective, cx).0)
            / (2. * epsilon);
        close(gradient[i], numerical, 4e-6);
    }
}

#[test]
fn g3_mixed_acceleration_jumps_and_displacement_retain_all_control_material_and_reaction_gradients()
{
    with_cx(false, |cx| check_mixed(false, cx));
}

#[test]
fn g3_acceleration_and_jump_trajectory_stress_gradients_match_original_material_motion() {
    with_cx(false, |cx| check_mixed(true, cx));
}

fn payload(times: &[f64], values: &[f64], interpolation: TableInterpolation) -> VectorPayload {
    VectorPayload::new(
        PayloadMeta::new(
            QuantityContract::Dimensions(ACCELERATION_DIMS),
            PayloadId::new("basis/world-cartesian").unwrap(),
            WORLD,
            OrientationParity::Even,
            ReferenceSemantics::Continuous,
        )
        .unwrap(),
        SampleSource::table(
            times
                .iter()
                .map(|&t| QtyAny::new(t, Dims([0, 0, 1, 0, 0, 0])))
                .collect(),
            values
                .iter()
                .map(|&a| {
                    [a, 0., 0.]
                        .map(|v| QtyAny::new(v, ACCELERATION_DIMS))
                        .to_vec()
                })
                .collect(),
            interpolation,
            OutsideDomainPolicy::Clamp,
        )
        .unwrap(),
    )
    .unwrap()
}
fn history(frames: &FrameTree) -> TetLoadHistory3<'_> {
    TetLoadHistory3 {
        surfaces: &[],
        volume: VolumeHistory3::Fixed(VolumeLoad3::None),
        nodal_forces_n: &[],
        frames,
    }
}

#[test]
fn g1_recorded_typed_pulse_has_exact_left_right_motion_and_analytic_amplitude_gradient() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [121; 32]).unwrap();
        let problem = problem(&material, &FIXED);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let amplitude = 1. / 64.;
        let source = payload(
            &[0., 0.125, 0.25, 0.5],
            &[0., amplitude, 0., 0.],
            TableInterpolation::StepLeft,
        );
        let supports = [SupportAcceleration3 {
            nodes: &ALL,
            components: [true, false, false],
            acceleration: &source,
        }];
        let schedule =
            EventLoadHistory3::new_with_support_accelerations(history(&frames), &supports, 16, cx)
                .unwrap();
        let run = schedule
            .begin_transient_with_support_accelerations(
                &geometry,
                &problem,
                None,
                0.,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap()
            .into_run();
        let mut recorded = run
            .record_prescribed_controls(
                Default::default(),
                Default::default(),
                Default::default(),
                cx,
            )
            .unwrap();
        let mut jumps = Vec::new();
        while recorded.run().state().time_s < 0.5 {
            let report = schedule
                .advance_recorded_accelerations_to(&mut recorded, 0.5, 0.2, cx)
                .unwrap();
            if let Some(jump) = report.jump {
                jumps.push(jump.time_s);
            }
        }
        assert_eq!(jumps, vec![0.125, 0.25]);
        assert_eq!(recorded.sample_count(), 7);
        assert_eq!(recorded.run().state().step, 4);
        let full = recorded.run().nodal_motion(cx).unwrap();
        let q_per_a = 0.5 * (0.375_f64.powi(2) - 0.25_f64.powi(2));
        close(full[0][0][0], amplitude * q_per_a, 1e-13);
        close(full[1][0][0], amplitude * 0.125, 1e-13);
        let mut seed = zero_seed();
        seed.q[0][0] = 3.;
        seed.v[0][0] = full[1][0][0];
        let gradient = recorded.reverse(&seed, cx).unwrap();
        let derivative: f64 = gradient
            .accelerations
            .iter()
            .filter(|row| {
                (row.kind == ControlForceKind3::RightLimit && row.time_s == 0.125)
                    || (row.kind == ControlForceKind3::IntervalEndpoint
                        && row.time_s > 0.125
                        && row.time_s <= 0.25)
            })
            .map(|row| [0, 3, 6, 9].into_iter().map(|d| row.values[d]).sum::<f64>())
            .sum();
        close(derivative, 3. * q_per_a + full[1][0][0] * 0.125, 1e-12);
        assert_eq!(gradient.displacements.len(), 1);
        assert_eq!(gradient.accelerations.len(), 6);
    });
}

#[test]
fn g4_late_jump_refusal_and_recording_limits_preserve_primal_motion_and_material_tape() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [121; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let problem = problem(&material, &fixed);
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
                Default::default(),
                cx,
            )
            .unwrap();
        let mut recorded = run
            .clone()
            .record_prescribed_controls(
                PrescribedRecordingBudget3 {
                    maximum_steps: 2,
                    ..Default::default()
                },
                Default::default(),
                Default::default(),
                cx,
            )
            .unwrap();
        recorded
            .enable_material_recording(Default::default(), cx)
            .unwrap();
        let before = recorded.run().state().clone();
        let motion = recorded.run().nodal_motion(cx).unwrap();
        let gradient = recorded.reverse_design(&zero_seed(), cx).unwrap();
        assert!(
            recorded
                .advance_to_support_acceleration_jump(
                    &load,
                    &[(0, 0.02)],
                    &load,
                    &[(3, 1.)],
                    0.125,
                    cx
                )
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert_eq!(
                recorded.advance_to_support_acceleration_jump(
                    &load,
                    &[(0, 0.02)],
                    &load,
                    &[(0, 0.)],
                    0.125,
                    cancelled
                ),
                Err(TransientAdjointError3::Cancelled)
            );
        });
        assert_eq!(recorded.run().state(), &before);
        assert_eq!(recorded.run().nodal_motion(cx).unwrap(), motion);
        assert_eq!(recorded.sample_count(), 1);
        assert_eq!(recorded.reverse_design(&zero_seed(), cx).unwrap(), gradient);

        let mut limited = run
            .record_prescribed_controls(
                Default::default(),
                Default::default(),
                Default::default(),
                cx,
            )
            .unwrap();
        limited
            .enable_material_recording(
                PrescribedMaterialRecordingBudget3 {
                    maximum_motion_values: 72,
                    ..Default::default()
                },
                cx,
            )
            .unwrap();
        assert_eq!(
            limited.advance_to_support_acceleration_jump(
                &load,
                &[(0, 0.02)],
                &load,
                &[(0, 0.)],
                0.125,
                cx
            ),
            Err(TransientAdjointError3::Budget)
        );
        assert_eq!(limited.sample_count(), 1);
        recorded
            .advance_to_support_acceleration_jump(&load, &[(0, 0.02)], &load, &[(0, 0.)], 0.125, cx)
            .unwrap();
        let accepted = recorded.run().nodal_motion(cx).unwrap();
        assert_eq!(
            recorded.apply_support_acceleration_jump(&load, &[(0, 0.01)], cx),
            Err(TransientAdjointError3::Budget)
        );
        assert_eq!(recorded.sample_count(), 3);
        assert_eq!(recorded.run().nodal_motion(cx).unwrap(), accepted);
    });
}

#[test]
fn g3_typed_smooth_acceleration_matches_direct_recording_and_jump_discards_old_acceleration() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6., 60., 0., [121; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let history = history(&frames);
        let source = payload(&[0., 0.5], &[0.02, 0.03], TableInterpolation::Linear);
        let supports = [SupportAcceleration3 {
            nodes: &[0],
            components: [true, false, false],
            acceleration: &source,
        }];
        let run = history
            .begin_transient_with_support_acceleration_history(
                &geometry,
                &problem,
                None,
                &supports,
                0.,
                &ZERO,
                &ZERO,
                config(0.4, 0.006),
                cx,
            )
            .unwrap();
        let mut recorded = run
            .record_prescribed_controls(
                Default::default(),
                TransientAdjointConfig3::default(),
                SupportAdjointBudget3::default(),
                cx,
            )
            .unwrap();
        let mut direct = recorded.clone();
        recorded
            .advance_history_with_support_accelerations(&history, &supports, 0.125, cx)
            .unwrap();
        let endpoint = history
            .at_with_support_accelerations(&geometry, &problem, &supports, 0.125, cx)
            .unwrap();
        let load = geometry
            .assemble(&[], VolumeLoad3::None, endpoint.nodal_forces_n(), cx)
            .unwrap();
        direct
            .advance_with_support_accelerations(&load, endpoint.prescribed_dofs_m_s2(), 0.125, cx)
            .unwrap();
        assert_eq!(
            recorded.run().nodal_motion(cx).unwrap(),
            direct.run().nodal_motion(cx).unwrap()
        );
        let mut seed = zero_seed();
        seed.q[1][0] = 0.3;
        seed.v[0][0] = -0.2;
        seed.a[1][0] = 0.7;
        seed.a[0][0] = -0.3;
        assert_eq!(
            recorded.reverse(&seed, cx).unwrap(),
            direct.reverse(&seed, cx).unwrap()
        );
        let jump = recorded
            .adjoint()
            .pullback_acceleration_jump(&seed, cx)
            .unwrap();
        assert_eq!(jump.previous.a, ZERO);
        // Independent consistent tet mass M_ff=.1, M_fc=.05: lambda=.7/.1,
        // da_free/da_fixed=-.5; direct constrained acceleration seed is -.3.
        close(jump.force[0], 7., 1e-12);
        close(jump.prescribed_acceleration[0], -0.3 - 0.5 * 0.7, 1e-12);
    });
}
