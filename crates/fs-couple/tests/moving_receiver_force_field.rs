//! G1/G3 complete past-force fields on the original moving microphone clock.
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource, ReflectedTrajectoryBudget, TimedSourceKinematics,
    TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    TimedCompactSource,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::observation::moving::field::CompactFieldSample;
use fs_couple::pcm_wav::observation::moving::force::{ForceHistoryQuadrature, MovingForceRenderer};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}
fn path(times: &[f64], motion: impl Fn(f64) -> SourceKinematics) -> PreparedKinematicTrajectory {
    let knots: Vec<_> = times
        .iter()
        .map(|&source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect();
    PreparedKinematicTrajectory::new(
        medium(),
        &knots,
        0.8,
        TrajectoryBudget {
            maximum_segments: 16,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 8192,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}
fn stationary(position_m: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m,
        velocity_m_s: [0.0; 3],
        acceleration_m_s2: [0.0; 3],
    }
}
fn receiver(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            12.0 - 8.0 * t + 3.0 * t * t,
            7.0 + 2.0 * t,
            5.0 + 0.5 * t * t,
        ],
        velocity_m_s: [-8.0 + 6.0 * t, 2.0, t],
        acceleration_m_s2: [6.0, 0.0, 1.0],
    }
}
fn pair_config() -> MovingReceiverConfig {
    MovingReceiverConfig {
        receiver_time_interval_s: [0.08, 0.18],
        minimum_distance_m: 0.001,
        maximum_time_residual_s: 1e-12,
    }
}
fn pair_budget() -> MovingReceiverBudget {
    MovingReceiverBudget {
        maximum_segment_pairs: 256,
        maximum_bound_nodes: 32_768,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}
fn pair<'a>(
    source: &'a PreparedKinematicTrajectory,
    receiver: &'a PreparedKinematicTrajectory,
) -> PreparedMovingReceiverSource<'a> {
    PreparedMovingReceiverSource::new(source, receiver, pair_config(), pair_budget(), &mut || {
        false
    })
    .unwrap()
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.10,
        sample_rate_hz: 1000,
        samples: 17,
        max_block: 17,
        receiver_bandwidth_hz: 20.0,
        maximum_time_residual_s: 1e-12,
    }
}
fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 8,
        maximum_source_samples: 256,
        maximum_output_samples: 4096,
        maximum_payload_bytes: 1_000_000,
        maximum_source_receiver_pairs: 65_536,
    }
}
fn quadrature() -> ForceHistoryQuadrature {
    ForceHistoryQuadrature {
        absolute_tolerance_n_s_m3: 1e-13,
        maximum_subdivision_depth: 24,
        maximum_evaluations: 1_000_000,
    }
}
fn strength(t: f64, direction: [f64; 3]) -> CompactSourceStrength {
    CompactSourceStrength {
        force_on_fluid_n: direction.map(|d| d * (t * t + 0.3 * t * t * t)),
        force_rate_on_fluid_n_s: direction.map(|d| d * (2.0 * t + 0.9 * t * t)),
        ..Default::default()
    }
}
fn history(direction: [f64; 3]) -> [TimedCompactSource; 3] {
    [0.0, 0.1, 0.3].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: strength(source_time_s, direction),
    })
}
fn norm(x: [f64; 3]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

// Closed-form stationary Green field, evaluated at the CURRENT receiver point.
// The microphone's velocity affects emission time only, never field amplitude.
fn stationary_field(
    source: [f64; 3],
    x: [f64; 3],
    time: f64,
    direction: [f64; 3],
) -> (f64, [f64; 3]) {
    let r = core::array::from_fn(|i| x[i] - source[i]);
    let distance = norm(r);
    let n = r.map(|v| v / distance);
    let tau = time - distance / 300.0;
    assert!(tau > 0.0);
    let value = strength(tau, direction);
    let f = value.force_on_fluid_n;
    let fdot = value.force_rate_on_fluid_n_s;
    let impulse = direction.map(|d| d * (tau.powi(3) / 3.0 + 0.075 * tau.powi(4)));
    let four_pi = 4.0 * core::f64::consts::PI;
    let pressure = (dot(n, f) / distance.powi(2) + dot(n, fdot) / (300.0 * distance)) / four_pi;
    let velocity = core::array::from_fn(|i| {
        ((3.0 * n[i] * dot(n, impulse) - impulse[i]) / distance.powi(3)
            + (3.0 * n[i] * dot(n, f) - f[i]) / (300.0 * distance.powi(2))
            + n[i] * dot(n, fdot) / (300.0 * 300.0 * distance))
            / (four_pi * 1.2)
    });
    (pressure, velocity)
}

#[test]
fn moving_reflected_array_matches_analytic_current_position_history_and_cross_flux() {
    let positions = [[2.0, 1.0, 0.5], [1.0, 3.0, 1.0]];
    let directions = [[1.0, -0.4, 0.7], [-0.3, 0.8, 0.2]];
    let sources = positions.map(|x| path(&[0.0, 0.04, 0.2, 0.3], |_| stationary(x)));
    let receiver = path(&[0.0, 0.09, 0.18, 0.3], receiver);
    let walls = MovingReflectionObserver::new(
        &receiver,
        &[
            IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
            IdealAcousticPlane::new(
                [0.0; 3],
                [0.0, 1.0, 0.0],
                IdealPlaneBoundary::PressureRelease,
            )
            .unwrap(),
        ],
        4,
        ReflectedTrajectoryBudget {
            maximum_payload_bytes: 1_000_000,
            maximum_bound_nodes: 8192,
            maximum_subdivision_depth: 20,
        },
        &mut || false,
    )
    .unwrap();
    let motions: Vec<_> = sources
        .iter()
        .map(|source| {
            PreparedReflectedMovingReceiverSource::new(
                source,
                &walls,
                pair_config(),
                ReflectedTrajectoryBudget {
                    maximum_payload_bytes: 1_000_000,
                    maximum_bound_nodes: 8192,
                    maximum_subdivision_depth: 20,
                },
                pair_budget(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let histories = directions.map(history);
    let tracks: Vec<_> = motions
        .iter()
        .zip(&histories)
        .map(|(motion, samples)| CompactSourceHistory {
            motion,
            samples: &samples[..],
            source_bandwidth_hz: 1.0,
        })
        .collect();
    let microphone = FirstOrder::new(0.25, [0.0, 1.0, 0.0]).unwrap();
    let mut renderer = MovingForceRenderer::with_moving_receiver_reflections(
        &tracks,
        microphone,
        config(),
        budget(),
        quadrature(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(renderer.source_count(), 2);
    assert_eq!(renderer.source_path_count(), 8);
    assert_eq!(renderer.quadrature_evaluations(), 0);
    let mut fields = [CompactFieldSample::default(); 17];
    renderer.field_block(&mut fields).unwrap();
    let mut largest_cross_term = 0.0_f64;
    for field in fields {
        let x = receiver
            .kinematics_at(field.receiver_time_s)
            .unwrap()
            .position_m;
        let mut pressure = 0.0;
        let mut velocity = [0.0; 3];
        let mut diagonal_flux = [0.0; 3];
        for source in 0..2 {
            for mask in 0..4 {
                let reflection = [
                    if mask & 1 == 0 { 1.0 } else { -1.0 },
                    if mask & 2 == 0 { 1.0 } else { -1.0 },
                    1.0,
                ];
                let parity = reflection[1];
                let virtual_x = core::array::from_fn(|i| reflection[i] * x[i]);
                let (p, u) = stationary_field(
                    positions[source],
                    virtual_x,
                    field.receiver_time_s,
                    directions[source],
                );
                let p = parity * p;
                let u: [f64; 3] = core::array::from_fn(|i| parity * reflection[i] * u[i]);
                pressure += p;
                for i in 0..3 {
                    velocity[i] += u[i];
                    diagonal_flux[i] += p * u[i];
                }
            }
        }
        close(field.pressure_pa, pressure, 1e-17, 3e-9);
        for i in 0..3 {
            close(field.particle_velocity_m_s[i], velocity[i], 1e-20, 3e-9);
            close(
                field.instantaneous_intensity_w_m2[i],
                pressure * velocity[i],
                1e-25,
                6e-9,
            );
            largest_cross_term =
                largest_cross_term.max((pressure * velocity[i] - diagonal_flux[i]).abs());
        }
        close(
            field.microphone_pressure_pa,
            0.25 * pressure - 0.75 * 360.0 * velocity[1],
            1e-17,
            3e-9,
        );
    }
    assert!(largest_cross_term > 1e-15);
}

fn curved(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            2.0 + 20.0 * t + 60.0 * t * t,
            1.0 + 4.0 * t - 10.0 * t * t,
            0.2,
        ],
        velocity_m_s: [20.0 + 120.0 * t, 4.0 - 20.0 * t, 0.0],
        acceleration_m_s2: [120.0, -20.0, 0.0],
    }
}
fn close_receiver(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [6.0 - 2.0 * t + 7.0 * t * t, 2.0 + t, 1.0],
        velocity_m_s: [-2.0 + 14.0 * t, 1.0, 0.0],
        acceleration_m_s2: [14.0, 0.0, 0.0],
    }
}
// Independent bisection and dense fixed Simpson integration on analytic inputs.
fn reference_history(x: [f64; 3], receiver_time: f64) -> (f64, [f64; 3]) {
    let (mut lo, mut hi) = (0.0_f64, 0.3_f64);
    for _ in 0..80 {
        let t = lo.midpoint(hi);
        let r = norm(core::array::from_fn(|i| x[i] - curved(t).position_m[i]));
        if t + r / 300.0 < receiver_time {
            lo = t;
        } else {
            hi = t;
        }
    }
    let tau = lo.midpoint(hi);
    let steps = 8192;
    let h = tau / steps as f64;
    let mut total = [0.0; 3];
    for j in 0..=steps {
        let t = j as f64 * h;
        let r = core::array::from_fn(|i| x[i] - curved(t).position_m[i]);
        let distance = norm(r);
        let n = r.map(|v| v / distance);
        let f = strength(t, [1.0, -0.4, 0.7]).force_on_fluid_n;
        let weight = if j == 0 || j == steps {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for i in 0..3 {
            total[i] += weight * (3.0 * n[i] * dot(n, f) - f[i]) / distance.powi(3);
        }
    }
    (tau, total.map(|v| v * h / 3.0))
}

#[test]
fn accelerating_source_and_microphone_integrate_the_same_past_at_each_arrival() {
    let source = path(&[0.0, 0.03, 0.07, 0.11, 0.18, 0.3], curved);
    let receiver = path(&[0.0, 0.09, 0.18, 0.3], close_receiver);
    let motion = pair(&source, &receiver);
    // Many geometry knots inside ONE original cubic force interval.
    let samples = [0.0, 0.3].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: strength(source_time_s, [1.0, -0.4, 0.7]),
    });
    let track = CompactSourceHistory {
        motion: &motion,
        samples: &samples,
        source_bandwidth_hz: 1.0,
    };
    let clock = MovingPressureConfig {
        first_receiver_time_s: 0.12,
        samples: 9,
        max_block: 9,
        ..config()
    };
    let mut pressure_owner =
        MovingCompactRenderer::with_moving_receiver(&[track], clock, budget(), &mut || false)
            .unwrap();
    let mut pressures = [0.0; 9];
    pressure_owner.block(&mut pressures).unwrap();
    let mut errors = [0.0_f64; 3];
    for (rung, tolerance) in [1e-6, 1e-10, 1e-14].into_iter().enumerate() {
        let mut renderer = MovingForceRenderer::with_moving_receiver(
            &[track],
            FirstOrder::default(),
            clock,
            budget(),
            ForceHistoryQuadrature {
                absolute_tolerance_n_s_m3: tolerance,
                ..quadrature()
            },
            &mut || false,
        )
        .unwrap();
        let mut fields = [CompactFieldSample::default(); 9];
        renderer.field_block(&mut fields).unwrap();
        for (index, field) in fields.iter().enumerate() {
            close(field.pressure_pa, pressures[index], 1e-18, 3e-12);
            let x = close_receiver(field.receiver_time_s).position_m;
            let (tau, j) = reference_history(x, field.receiver_time_s);
            let observer = MovingCompactObserver::new(
                medium(),
                x,
                MovingCompactLimits {
                    maximum_mach: 0.8,
                    minimum_distance_m: 0.001,
                },
            )
            .unwrap();
            let expected = observer
                .prepare_force_emission(tau, curved(tau))
                .unwrap()
                .evaluate(strength(tau, [1.0, -0.4, 0.7]), j)
                .unwrap();
            for i in 0..3 {
                errors[rung] = errors[rung].max(
                    (field.particle_velocity_m_s[i] - expected.particle_velocity_m_s[i]).abs(),
                );
            }
        }
    }
    assert!(errors[2] < 2e-12, "{errors:?}");
    assert!(errors[2] <= errors[0] + 1e-14, "{errors:?}");
}

#[test]
fn mixed_volume_force_sources_keep_weak_fields_after_coherent_cancellation() {
    let source = path(&[0.0, 0.3], |_| stationary([2.0, 1.0, 0.5]));
    let receiver = path(&[0.0, 0.3], receiver);
    let motion = pair(&source, &receiver);
    let weak = history([2.0_f64.powi(-60), -2.0_f64.powi(-61), 0.0]);
    let mut mixed = weak;
    let mut opposite = history([0.0; 3]);
    for sample in &mut mixed {
        sample.strength.volume_flow_m3_s = 1.0 + sample.source_time_s;
        sample.strength.volume_flow_rate_m3_s2 = 1.0;
    }
    for sample in &mut opposite {
        sample.strength.volume_flow_m3_s = -(1.0 + sample.source_time_s);
        sample.strength.volume_flow_rate_m3_s2 = -1.0;
    }
    let track = |samples| CompactSourceHistory {
        motion: &motion,
        samples,
        source_bandwidth_hz: 1.0,
    };
    let render = |tracks: &[_]| {
        let mut renderer = MovingForceRenderer::with_moving_receiver(
            tracks,
            FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap(),
            config(),
            budget(),
            quadrature(),
            &mut || false,
        )
        .unwrap();
        let mut output = [CompactFieldSample::default(); 17];
        renderer.field_block(&mut output).unwrap();
        output
    };
    let expected = render(&[track(&weak)]);
    let actual = render(&[track(&mixed), track(&opposite)]);
    for (a, e) in actual.iter().zip(&expected) {
        assert!(e.pressure_pa.abs() > 1e-26);
        close(a.pressure_pa, e.pressure_pa, 0.0, 2e-8);
        for i in 0..3 {
            close(
                a.particle_velocity_m_s[i],
                e.particle_velocity_m_s[i],
                1e-35,
                2e-8,
            );
        }
        close(
            a.microphone_pressure_pa,
            e.microphone_pressure_pa,
            1e-34,
            3e-8,
        );
    }
}

#[test]
fn moving_history_blocks_replay_exactly_and_failed_work_never_publishes() {
    let source = path(&[0.0, 0.03, 0.11, 0.3], curved);
    let receiver = path(&[0.0, 0.3], close_receiver);
    let motion = pair(&source, &receiver);
    let samples = history([1.0, -0.4, 0.7]);
    let track = CompactSourceHistory {
        motion: &motion,
        samples: &samples,
        source_bandwidth_hz: 1.0,
    };
    let make = |quad| {
        MovingForceRenderer::with_moving_receiver(
            &[track],
            FirstOrder::new(0.5, [0.0, 1.0, 0.0]).unwrap(),
            config(),
            budget(),
            quad,
            &mut || false,
        )
        .unwrap()
    };
    let mut whole = make(quadrature());
    let mut expected = [CompactFieldSample::default(); 17];
    whole.field_block(&mut expected).unwrap();
    let mut split = make(quadrature());
    let mut actual = [CompactFieldSample::default(); 17];
    for chunk in actual.chunks_mut(4) {
        split.field_block(chunk).unwrap();
    }
    assert_eq!(actual, expected);
    assert_eq!(
        split.quadrature_evaluations(),
        whole.quadrature_evaluations()
    );

    let sentinel = CompactFieldSample {
        pressure_pa: 123.0,
        ..Default::default()
    };
    let mut output = [sentinel; 17];
    let mut retry = make(quadrature());
    let mut polls = 0;
    assert!(
        retry
            .field_block_controlled(&mut output, &mut || {
                polls += 1;
                polls == 40
            })
            .is_err()
    );
    assert_eq!(output, [sentinel; 17]);
    assert_eq!(retry.samples_rendered(), 0);
    let failed_work = retry.quadrature_evaluations();
    assert!(failed_work > 0);
    retry.field_block(&mut output).unwrap();
    assert_eq!(output, expected);
    assert_eq!(
        retry.quadrature_evaluations(),
        whole.quadrature_evaluations() + failed_work
    );

    let mut limited = make(ForceHistoryQuadrature {
        maximum_evaluations: 4,
        ..quadrature()
    });
    let mut audio = [123.0; 1];
    assert!(limited.block(&mut audio).is_err());
    assert_eq!(audio, [123.0]);
    assert_eq!(limited.samples_rendered(), 0);
    assert_eq!(limited.quadrature_evaluations(), 4);
}

#[test]
fn causal_clearance_does_not_hide_a_collision_with_the_complete_force_past() {
    let source = path(&[0.0, 0.1, 0.4, 1.0], |t| SourceKinematics {
        position_m: [20.0 * t, 0.0, 0.0],
        velocity_m_s: [20.0, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    });
    let receiver = path(&[0.69, 0.7, 0.71], |t| SourceKinematics {
        position_m: [1.0, t - 0.7, 0.0],
        velocity_m_s: [0.0, 1.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    });
    let motion = PreparedMovingReceiverSource::new(
        &source,
        &receiver,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.69, 0.71],
            ..pair_config()
        },
        pair_budget(),
        &mut || false,
    )
    .unwrap();
    let samples = [0.0, 1.0].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: strength(source_time_s, [0.0, 0.0, 1.0]),
    });
    let track = CompactSourceHistory {
        motion: &motion,
        samples: &samples,
        source_bandwidth_hz: 0.1,
    };
    let clock = MovingPressureConfig {
        first_receiver_time_s: 0.69,
        sample_rate_hz: 100,
        samples: 3,
        max_block: 3,
        ..config()
    };
    let mut pressure =
        MovingCompactRenderer::with_moving_receiver(&[track], clock, budget(), &mut || false)
            .unwrap();
    pressure.block(&mut [0.0; 3]).unwrap();
    let mut renderer = MovingForceRenderer::with_moving_receiver(
        &[track],
        FirstOrder::default(),
        clock,
        budget(),
        quadrature(),
        &mut || false,
    )
    .unwrap();
    let sentinel = CompactFieldSample {
        pressure_pa: 123.0,
        ..Default::default()
    };
    let mut output = [sentinel; 3];
    assert!(renderer.field_block(&mut output).is_err());
    assert_eq!(output, [sentinel; 3]);
    assert_eq!(renderer.samples_rendered(), 0);
    assert!(renderer.quadrature_evaluations() > 0);
}
