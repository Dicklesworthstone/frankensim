//! G1/G3/G4: authored acceleration records drive the original inertial supports.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Vec3;
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, TimeSignal, WORLD};
use fs_solid::linear3::*;
use fs_solid::loads3::history::dynamics::events::EventLoadHistory3;
use fs_solid::loads3::history::dynamics::events::adaptive::{
    AdaptiveHistoryError3, TransientTimeControl3,
};
use fs_solid::loads3::history::supports::{DISPLACEMENT_DIMS, SupportAcceleration3};
use fs_solid::loads3::history::{ACCELERATION_DIMS, TetLoadHistory3, VolumeHistory3};
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.; 3]; 4];
const ALL: [usize; 4] = [0, 1, 2, 3];
const FIXED: [usize; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
const FRAME: ContentHash = ContentHash([105; 32]);

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
                seed: 105,
                kernel_id: 105,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(6., 60., 0., [105; 32]).unwrap()
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
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap()
}
fn history(frames: &FrameTree) -> TetLoadHistory3<'_> {
    TetLoadHistory3 {
        surfaces: &[],
        volume: VolumeHistory3::Fixed(VolumeLoad3::None),
        nodal_forces_n: &[],
        frames,
    }
}
fn table(
    dims: Dims,
    frame: FrameId,
    times: &[f64],
    values: &[f64],
    interpolation: TableInterpolation,
) -> VectorPayload {
    VectorPayload::new(
        PayloadMeta::new(
            QuantityContract::Dimensions(dims),
            PayloadId::new(if frame == WORLD {
                "basis/world-cartesian"
            } else {
                "basis/frame-cartesian"
            })
            .unwrap(),
            frame,
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
                .map(|&x| [x, 0., 0.].map(|v| QtyAny::new(v, dims)).to_vec())
                .collect(),
            interpolation,
            OutsideDomainPolicy::Clamp,
        )
        .unwrap(),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-11,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_native_acceleration_switch_preserves_motion_energy_and_full_force_ownership() {
    with_cx(false, |cx| {
        let material = material();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let before_load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let after_load = geometry
            .assemble(&[], VolumeLoad3::None, &[(0, 0.034), (3, 0.012)], cx)
            .unwrap();
        let mut q = ZERO;
        q[0][0] = 0.001;
        q[1][0] = 0.0003;
        let mut v = ZERO;
        v[0][0] = 0.002;
        v[1][0] = -0.001;
        let mut a = ZERO;
        a[0][0] = 0.004;
        let config = TransientConfig3 {
            mass_damping_per_s: 0.4,
            stiffness_damping_s: 0.006,
            ..Default::default()
        };
        let mut run = before_load
            .begin_transient_with_supports(&problem, None, 0., &q, &v, &a, config, cx)
            .unwrap();
        let before = run.state().clone();
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        let next_acceleration = 0.02;
        let jump = run
            .apply_support_acceleration_jump(&after_load, &[(0, next_acceleration)], cx)
            .unwrap();
        // Original partition: M_ff=.1, M_fc=.05. q/v and both Rayleigh
        // contributions are unchanged; only f_f and M_fc*a_c change.
        let delta_free = (0.012 - 0.05 * (next_acceleration - a[0][0])) / 0.1;
        close(
            run.state().acceleration[0],
            before.acceleration[0] + delta_free,
        );
        close(jump.acceleration_change_l2_m_s2, delta_free.abs());
        assert_eq!(run.state().displacement, before.displacement);
        assert_eq!(run.state().velocity, before.velocity);
        assert_eq!(run.state().time_s, before.time_s);
        assert_eq!(run.state().step, before.step);
        assert_eq!(run.state().kinetic_energy_j, before.kinetic_energy_j);
        assert_eq!(run.state().strain_energy_j, before.strain_energy_j);
        assert_eq!(run.state().total_energy_j, before.total_energy_j);
        let full = run.nodal_motion(cx).unwrap();
        assert_eq!(full[0], q);
        assert_eq!(full[1], v);
        close(full[2][0][0], next_acceleration);
        let accepted = run.recover_response(FRAME, 0.01, cx).unwrap();
        close(
            accepted.support_reactions_n[0][0] - response.support_reactions_n[0][0],
            0.1 * (next_acceleration - a[0][0]) + 0.05 * delta_free - 0.034,
        );
        assert_eq!(accepted.external_force_n[0][0], 0.034);
        assert_eq!(accepted.external_force_n[1][0], 0.012);
        assert!(accepted.free_force_residual_n < 1e-11);
        for invalid in [
            vec![(3, 1.)],
            vec![(12, 1.)],
            vec![(0, 1.), (0, 1.)],
            vec![(1, 1.), (0, 1.)],
            vec![(0, f64::NAN)],
            vec![(0, f64::INFINITY)],
        ] {
            assert!(
                run.apply_support_acceleration_jump(&before_load, &invalid, cx)
                    .is_err()
            );
            assert_eq!(run.nodal_motion(cx).unwrap(), full);
            assert_eq!(run.recover_response(FRAME, 0.01, cx).unwrap(), accepted);
        }
        with_cx(true, |cancelled| {
            assert_eq!(
                run.apply_support_acceleration_jump(&before_load, &[], cancelled),
                Err(TransientError3::Cancelled)
            );
        });
        assert_eq!(run.recover_response(FRAME, 0.01, cx).unwrap(), accepted);
        // Omitting an acceleration stops neither position nor velocity.
        run.apply_support_acceleration_jump(&after_load, &[], cx)
            .unwrap();
        let reset = run.nodal_motion(cx).unwrap();
        assert_eq!(reset[0], q);
        assert_eq!(reset[1], v);
        assert_eq!(reset[2][0][0], 0.);
    });
}

#[test]
fn g1_held_acceleration_pulse_has_no_precursor_motion_and_retains_support_work() {
    with_cx(false, |cx| {
        let material = material();
        let problem = problem(&material, &FIXED);
        let geometry = geometry(cx);
        let frames = FrameTree::new();
        let amplitude = 1.0 / 64.0;
        let payload = table(
            ACCELERATION_DIMS,
            WORLD,
            &[0., 0.125, 0.25, 0.5],
            &[0., amplitude, 0., 0.],
            TableInterpolation::StepLeft,
        );
        let rows = [SupportAcceleration3 {
            nodes: &ALL,
            components: [true, false, false],
            acceleration: &payload,
        }];
        let schedule =
            EventLoadHistory3::new_with_support_accelerations(history(&frames), &rows, 16, cx)
                .unwrap();
        assert_eq!(schedule.knots_s(), &[0., 0.125, 0.25, 0.5]);
        assert!(
            EventLoadHistory3::new_with_support_accelerations(history(&frames), &rows, 15, cx)
                .is_err()
        );
        let q0 = [[1.0 / 1024.0, 0., 0.]; 4];
        let mut run = schedule
            .begin_transient_with_support_accelerations(
                &geometry,
                &problem,
                None,
                0.,
                &q0,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        let mut q = q0[0][0];
        let mut v = 0.;
        let mut work = 0.;
        let mut jumps = Vec::new();
        while run.run().state().time_s < 0.5 {
            let report = run.advance_to(0.5, 1., cx).unwrap();
            let a = if report.start_time_s == 0.125 {
                amplitude
            } else {
                0.
            };
            let h = report.end_time_s - report.start_time_s;
            q += h * v + 0.5 * h * h * a;
            v += h * a;
            let motion = run.run().nodal_motion(cx).unwrap();
            for node in 0..4 {
                close(motion[0][node][0], q);
                close(motion[1][node][0], v);
            }
            if report.end_time_s == 0.125 {
                assert_eq!(motion[0], q0);
                assert_eq!(motion[1], ZERO);
                assert_eq!(report.interval.support_work_j, 0.);
            }
            if let Some(jump) = report.jump {
                jumps.push(jump.time_s);
            }
            work += report.interval.support_work_j;
            close(report.interval.balance_error_j, 0.);
            let current_acceleration = if report.end_time_s == 0.125 {
                amplitude
            } else {
                0.
            };
            close(motion[2][0][0], current_acceleration);
            let response = run.run().recover_response(FRAME, 0.01, cx).unwrap();
            close(
                response.support_reactions_n.iter().map(|r| r[0]).sum(),
                current_acceleration,
            );
        }
        assert_eq!(jumps, vec![0.125, 0.25]);
        close(work, 0.5 * v * v); // Original tetra has mass rho*V=1 kg.
        close(run.run().state().total_energy_j, work);
    });
}

#[test]
fn g1_combined_base_and_body_acceleration_events_converge_to_the_independent_oscillator() {
    with_cx(false, |cx| {
        let material = material();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let frames = FrameTree::new();
        let a = 1.0 / 64.0;
        let b = 1.0 / 128.0;
        let base = table(
            ACCELERATION_DIMS,
            WORLD,
            &[0., 0.125, 0.25, 0.5],
            &[0., a, 0., 0.],
            TableInterpolation::StepLeft,
        );
        let body = table(
            ACCELERATION_DIMS,
            WORLD,
            &[0., 0.0625, 0.375, 0.5],
            &[0., b, 0., 0.],
            TableInterpolation::StepLeft,
        );
        let rows = [SupportAcceleration3 {
            nodes: &[0],
            components: [true, false, false],
            acceleration: &base,
        }];
        let history = TetLoadHistory3 {
            volume: VolumeHistory3::Acceleration(&body),
            ..history(&frames)
        };
        let schedule =
            EventLoadHistory3::new_with_support_accelerations(history, &rows, 32, cx).unwrap();
        assert_eq!(schedule.knots_s(), &[0., 0.0625, 0.125, 0.25, 0.375, 0.5]);
        let mut errors = Vec::new();
        for cap in [1.0 / 32.0, 1.0 / 64.0, 1.0 / 128.0] {
            let mut run = schedule
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
                .unwrap();
            let mut work = 0.;
            let mut jumps = Vec::new();
            while run.run().state().time_s < 0.5 {
                let report = run.advance_to(0.5, cap, cx).unwrap();
                work += report.interval.mechanical_work_j + report.interval.support_work_j;
                close(report.interval.balance_error_j, 0.);
                if let Some(jump) = report.jump {
                    let expected = if jump.time_s == 0.125 || jump.time_s == 0.25 {
                        0.5 * a
                    } else {
                        2.5 * b
                    };
                    close(jump.acceleration_change_l2_m_s2, expected);
                    jumps.push(jump.time_s);
                }
            }
            assert_eq!(jumps, vec![0.0625, 0.125, 0.25, 0.375]);
            close(run.run().state().total_energy_j, work);
            // M_ff=.1, M_fc=.05, K_ff=10, K_fc=-10. The support
            // acceleration step response solves q''+100q=100s-.5s''.
            let base_step = |t: f64| 0.5 * a * t * t - 0.015 * a * (1. - (10. * t).cos());
            let body_step = |t: f64| 0.025 * b * (1. - (10. * t).cos());
            let expected =
                base_step(0.375) - base_step(0.25) + body_step(0.4375) - body_step(0.125);
            errors.push((run.run().state().displacement[0] - expected).abs());
            let full = run.run().nodal_motion(cx).unwrap();
            close(
                full[0][0][0],
                0.5 * a * (0.375_f64.powi(2) - 0.25_f64.powi(2)),
            );
            close(full[1][0][0], a * 0.125);
        }
        assert!(
            errors[1] < errors[0] / 3.5 && errors[2] < errors[1] / 3.5,
            "{errors:?}"
        );
        assert!(errors[2] < 5e-7, "{errors:?}");
    });
}

#[test]
fn g3_clamped_acceleration_values_keep_the_physical_frame_clock_and_initial_state() {
    with_cx(false, |cx| {
        let material = material();
        let problem = problem(&material, &FIXED);
        let geometry = geometry(cx);
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(1),
            name: "acceleration axes".into(),
            parent: WORLD,
            motion: FrameMotion::Rotating {
                axis: [0., 0., 1.],
                center: Vec3::new(0., 0., 0.),
                rate: QtyAny::new(2., Dims([0, 0, -1, 0, 0, 0])),
            },
        });
        let amplitude = 1.0 / 64.0;
        let payload = table(
            ACCELERATION_DIMS,
            FrameId(1),
            &[0., 0.1],
            &[amplitude, amplitude],
            TableInterpolation::Linear,
        );
        let rows = [SupportAcceleration3 {
            nodes: &ALL,
            components: [true; 3],
            acceleration: &payload,
        }];
        let history = history(&frames);
        let q = [[0.001, 0.0005, 0.]; 4];
        let v = [[0.002, -0.001, 0.]; 4];
        let mut run = history
            .begin_transient_with_support_acceleration_history(
                &geometry,
                &problem,
                None,
                &rows,
                1.,
                &q,
                &v,
                Default::default(),
                cx,
            )
            .unwrap();
        let initial = run.nodal_motion(cx).unwrap();
        assert_eq!(initial[0], q);
        assert_eq!(initial[1], v);
        close(initial[2][0][0], amplitude * 2.0_f64.cos());
        close(initial[2][0][1], amplitude * 2.0_f64.sin());
        let h = 0.125;
        run.advance_history_with_support_accelerations(&history, &rows, h, cx)
            .unwrap();
        let next = run.nodal_motion(cx).unwrap();
        let acceleration = [amplitude * 2.25_f64.cos(), amplitude * 2.25_f64.sin(), 0.];
        for c in 0..3 {
            close(next[2][0][c], acceleration[c]);
            close(
                next[1][0][c],
                v[0][c] + 0.5 * h * (initial[2][0][c] + acceleration[c]),
            );
            close(next[0][0][c], q[0][c] + 0.5 * h * (v[0][c] + next[1][0][c]));
        }
        assert_eq!(run.state().time_s, 1.125);
    });
}

#[test]
fn g1_frame_tilt_switches_acceleration_after_the_left_interval_without_a_velocity_kick() {
    with_cx(false, |cx| {
        let material = material();
        let problem = problem(&material, &FIXED);
        let geometry = geometry(cx);
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(1),
            name: "held acceleration axes".into(),
            parent: WORLD,
            motion: FrameMotion::Tilt {
                axis: [0., 0., 1.],
                center: Vec3::new(0., 0., 0.),
                angle: TimeSignal::Table {
                    times: vec![0., 0.25],
                    values: vec![0., core::f64::consts::FRAC_PI_2],
                    dims: Dims::NONE,
                    interp: fs_scenario::signal::Interp::Hold,
                },
            },
        });
        let amplitude = 1.0 / 64.0;
        let payload = table(
            ACCELERATION_DIMS,
            FrameId(1),
            &[0., 0.5],
            &[amplitude, amplitude],
            TableInterpolation::Linear,
        );
        let rows = [SupportAcceleration3 {
            nodes: &ALL,
            components: [true; 3],
            acceleration: &payload,
        }];
        let schedule =
            EventLoadHistory3::new_with_support_accelerations(history(&frames), &rows, 12, cx)
                .unwrap();
        let mut run = schedule
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
            .unwrap();
        let first = run.advance_to(0.5, 1., cx).unwrap();
        assert_eq!(first.end_time_s, 0.25);
        assert!(first.jump.is_some());
        let full = run.run().nodal_motion(cx).unwrap();
        close(full[0][0][0], 0.5 * amplitude * 0.25 * 0.25);
        close(full[1][0][0], amplitude * 0.25);
        assert_eq!(full[0][0][1], 0.);
        assert_eq!(full[1][0][1], 0.);
        close(full[2][0][0], 0.);
        close(full[2][0][1], amplitude);
        let second = run.advance_to(0.5, 1., cx).unwrap();
        let full = run.run().nodal_motion(cx).unwrap();
        close(full[1][0][0], amplitude * 0.25);
        close(full[1][0][1], amplitude * 0.25);
        close(
            run.run().state().total_energy_j,
            first.interval.support_work_j + second.interval.support_work_j,
        );
    });
}

#[test]
fn g3_adaptive_comparison_includes_integrated_acceleration_driven_supports() {
    with_cx(false, |cx| {
        let material = material();
        let problem = problem(&material, &FIXED);
        let geometry = geometry(cx);
        let frames = FrameTree::new();
        let jerk = 0.02;
        let payload = table(
            ACCELERATION_DIMS,
            WORLD,
            &[0., 1.],
            &[0., jerk],
            TableInterpolation::Linear,
        );
        let rows = [SupportAcceleration3 {
            nodes: &ALL,
            components: [true, false, false],
            acceleration: &payload,
        }];
        let schedule =
            EventLoadHistory3::new_with_support_accelerations(history(&frames), &rows, 8, cx)
                .unwrap();
        let mut run = schedule
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
            .unwrap();
        assert!(run.run().state().displacement.is_empty());
        let before = run.run().state().clone();
        let strict = TransientTimeControl3 {
            absolute_displacement_m: 1e-10,
            absolute_velocity_m_s: 1e-10,
            absolute_acceleration_m_s2: 1e-10,
            relative: 0.,
            maximum_step_s: 0.125,
            maximum_attempts: 1,
            ..Default::default()
        };
        assert!(matches!(
            run.advance_adaptive_to(0.125, 0.125, strict, cx),
            Err(AdaptiveHistoryError3::Rejected { .. })
        ));
        assert_eq!(run.run().state(), &before);
        assert_eq!(
            run.run().nodal_motion(cx).unwrap(),
            [ZERO.to_vec(), ZERO.to_vec(), ZERO.to_vec()]
        );
        let mut manual = run.clone();
        let first = manual.advance_to(0.0625, 0.0625, cx).unwrap();
        let second = manual.advance_to(0.125, 0.0625, cx).unwrap();
        let loose = TransientTimeControl3 {
            absolute_displacement_m: 0.01,
            absolute_velocity_m_s: 0.01,
            absolute_acceleration_m_s2: 0.01,
            ..strict
        };
        let accepted = run.advance_adaptive_to(0.125, 0.125, loose, cx).unwrap();
        assert_eq!(accepted.first, first);
        assert_eq!(accepted.second, second);
        assert_eq!(
            run.run().nodal_motion(cx).unwrap(),
            manual.run().nodal_motion(cx).unwrap()
        );
        close(
            run.run().nodal_motion(cx).unwrap()[0][0][0],
            3.0 / 16.0 * jerk * 0.125_f64.powi(3),
        );
    });
}

#[test]
fn g4_bad_units_bindings_and_a_late_acceleration_switch_leave_the_complete_run_unchanged() {
    with_cx(false, |cx| {
        let material = material();
        let problem = problem(&material, &FIXED);
        let geometry = geometry(cx);
        let frames = FrameTree::new();
        let history = history(&frames);
        let good = table(
            ACCELERATION_DIMS,
            WORLD,
            &[0., 0.25],
            &[0., 0.01],
            TableInterpolation::Linear,
        );
        let wrong = table(
            DISPLACEMENT_DIMS,
            WORLD,
            &[0., 0.25],
            &[0., 0.01],
            TableInterpolation::Linear,
        );
        let rows = [SupportAcceleration3 {
            nodes: &ALL,
            components: [true, false, false],
            acceleration: &good,
        }];
        let mut run = history
            .begin_transient_with_support_acceleration_history(
                &geometry,
                &problem,
                None,
                &rows,
                0.,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        let before = run.state().clone();
        let full = run.nodal_motion(cx).unwrap();
        for bad in [
            vec![SupportAcceleration3 {
                acceleration: &wrong,
                ..rows[0]
            }],
            vec![rows[0], rows[0]],
            vec![SupportAcceleration3 {
                nodes: &[4],
                ..rows[0]
            }],
        ] {
            assert!(
                run.advance_history_with_support_accelerations(&history, &bad, 0.125, cx)
                    .is_err()
            );
            assert_eq!(run.state(), &before);
            assert_eq!(run.nodal_motion(cx).unwrap(), full);
        }
        with_cx(true, |cancelled| {
            assert!(
                run.advance_history_with_support_accelerations(&history, &rows, 0.125, cancelled)
                    .is_err()
            );
        });
        assert_eq!(run.state(), &before);

        // The LEFT interval is finite and admissible. Only the RIGHT support
        // reaction overflows (nodal mass 2.5 kg times MAX acceleration), so
        // refusal must also roll back the already-computed private interval.
        let heavy = TetElasticMaterial::try_new(60., 60., 0., [106; 32]).unwrap();
        let heavy_problem = self::problem(&heavy, &FIXED);
        let pulse = table(
            ACCELERATION_DIMS,
            WORLD,
            &[0., 0.125],
            &[0., f64::MAX],
            TableInterpolation::StepLeft,
        );
        let rows = [SupportAcceleration3 {
            nodes: &ALL,
            components: [true, false, false],
            acceleration: &pulse,
        }];
        let schedule =
            EventLoadHistory3::new_with_support_accelerations(history, &rows, 8, cx).unwrap();
        let velocity = [[0.001, 0., 0.]; 4];
        let mut event = schedule
            .begin_transient_with_support_accelerations(
                &geometry,
                &heavy_problem,
                None,
                0.,
                &ZERO,
                &velocity,
                Default::default(),
                cx,
            )
            .unwrap();
        let before = event.run().state().clone();
        let full = event.run().nodal_motion(cx).unwrap();
        assert!(event.advance_to(0.25, 1., cx).is_err());
        assert_eq!(event.run().state(), &before);
        assert_eq!(event.run().nodal_motion(cx).unwrap(), full);
        let mut fork = event.clone();
        assert_eq!(
            event.advance_to(0.0625, 1., cx).unwrap(),
            fork.advance_to(0.0625, 1., cx).unwrap()
        );
    });
}
