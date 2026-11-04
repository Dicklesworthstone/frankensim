//! G1/G3 moving microphones: retarded clocks and Eulerian pressure stay distinct.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, PreparedSourceTrajectory, TimedSourceKinematics,
    TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

#[path = "moving_receiver/uniform_flow.rs"]
mod uniform_flow;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}

fn path_budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 16,
        maximum_payload_bytes: 1_000_000,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}

fn observation_budget() -> MovingReceiverBudget {
    MovingReceiverBudget {
        maximum_segment_pairs: 256,
        maximum_bound_nodes: 65_536,
        maximum_subdivision_depth: 24,
        maximum_root_iterations: 96,
    }
}

fn config(times: [f64; 2], clearance: f64) -> MovingReceiverConfig {
    MovingReceiverConfig {
        receiver_time_interval_s: times,
        minimum_distance_m: clearance,
        maximum_time_residual_s: 1e-13,
    }
}

fn samples(times: &[f64], motion: impl Fn(f64) -> SourceKinematics) -> Vec<TimedSourceKinematics> {
    times
        .iter()
        .map(|&source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect()
}

fn path(times: &[f64], motion: impl Fn(f64) -> SourceKinematics) -> PreparedKinematicTrajectory {
    PreparedKinematicTrajectory::new(
        medium(),
        &samples(times, motion),
        0.8,
        path_budget(),
        &mut || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, relative: f64, floor: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(floor),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn fixed(position_m: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m,
        ..SourceKinematics::default()
    }
}

fn source(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            0.2 + t * (12.0 + t * (8.0 - 4.0 * t)),
            -0.3 + t * (-6.0 + t * (5.0 + 3.0 * t)),
            0.1 + t * (4.0 - 3.0 * t),
        ],
        velocity_m_s: [
            12.0 + t * (16.0 - 12.0 * t),
            -6.0 + t * (10.0 + 9.0 * t),
            4.0 - 6.0 * t,
        ],
        acceleration_m_s2: [16.0 - 24.0 * t, 10.0 + 18.0 * t, -6.0],
    }
}

fn microphone(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            8.0 + t * (-9.0 + t * (4.0 + 2.0 * t)),
            5.0 + t * (4.0 - 2.0 * t),
            -3.0 + t * (3.0 + 2.0 * t * t),
        ],
        velocity_m_s: [-9.0 + t * (8.0 + 6.0 * t), 4.0 - 4.0 * t, 3.0 + 6.0 * t * t],
        acceleration_m_s2: [8.0 + 12.0 * t, -4.0, 12.0 * t],
    }
}

fn strength(t: f64) -> CompactSourceStrength {
    let (s, c) = (37.0 * t).sin_cos();
    let (s2, c2) = (23.0 * t).sin_cos();
    CompactSourceStrength {
        volume_flow_m3_s: 0.003 * s + 0.001 * c2,
        volume_flow_rate_m3_s2: 0.111 * c - 0.023 * s2,
        force_on_fluid_n: [2.0 * s, -3.0 * c, 0.5 * s2],
        force_rate_on_fluid_n_s: [74.0 * c, 111.0 * s, 11.5 * c2],
    }
}

fn distance(x: [f64; 3], y: [f64; 3]) -> f64 {
    (0..3).map(|i| (x[i] - y[i]).powi(2)).sum::<f64>().sqrt()
}

// Independent bisection over the explicit polynomial, not the prepared owner.
fn emission_at(x: [f64; 3], receiver_time: f64) -> f64 {
    let mut lo = -0.3;
    let mut hi = 0.4;
    assert!(lo + distance(x, source(lo).position_m) / 300.0 < receiver_time);
    assert!(hi + distance(x, source(hi).position_m) / 300.0 > receiver_time);
    for _ in 0..90 {
        let mid = lo.midpoint(hi);
        if mid + distance(x, source(mid).position_m) / 300.0 < receiver_time {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo.midpoint(hi)
}

fn potentials(x: [f64; 3], receiver_time: f64) -> [f64; 4] {
    let tau = emission_at(x, receiver_time);
    let state = source(tau);
    let r = distance(x, state.position_m);
    let radial_velocity = (0..3)
        .map(|i| (x[i] - state.position_m[i]) * state.velocity_m_s[i] / r)
        .sum::<f64>();
    let denominator = 4.0 * PI * r * (1.0 - radial_velocity / 300.0);
    let qf = strength(tau);
    [
        qf.volume_flow_m3_s / denominator,
        qf.force_on_fluid_n[0] / denominator,
        qf.force_on_fluid_n[1] / denominator,
        qf.force_on_fluid_n[2] / denominator,
    ]
}

fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
}

#[test]
fn curved_microphone_matches_independently_differentiated_retarded_potentials() {
    let source_path = path(&[-0.3, -0.04, 0.11, 0.4], source);
    let receiver_path = path(&[0.0, 0.07, 0.2], microphone);
    let observed = PreparedMovingReceiverSource::new(
        &source_path,
        &receiver_path,
        config([0.0, 0.2], 0.1),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    for t in [0.01, 0.06, 0.12, 0.19] {
        let x = microphone(t).position_m;
        let tau = emission_at(x, t);
        let actual = observed.prepare_receiver(t, &mut || false).unwrap();
        let pressure = actual
            .pressure
            .evaluate(strength(actual.pressure.source_time_s()))
            .unwrap();
        // The time derivative HOLDS x FIXED at this receiver event. Taking a
        // material derivative along the microphone would invent extra pressure.
        let monopole = 1.2 * derivative(|dt| potentials(x, t + dt)[0], 1e-5);
        let dipole = -(0..3)
            .map(|axis| {
                derivative(
                    |dx| {
                        let mut displaced = x;
                        displaced[axis] += dx;
                        potentials(displaced, t)[axis + 1]
                    },
                    1e-3,
                )
            })
            .sum::<f64>();
        near(pressure.source_time_s, tau, 2e-11, 0.1);
        near(pressure.monopole_pressure_pa, monopole, 2e-7, 1e-6);
        near(pressure.dipole_pressure_pa, dipole, 2e-7, 1e-6);
        near(pressure.pressure_pa, monopole + dipole, 2e-7, 1e-6);
        let clock_rate = derivative(
            |dt| emission_at(microphone(t + dt).position_m, t + dt),
            1e-5,
        );
        near(actual.observed_source_time_rate, clock_rate, 2e-9, 1.0);
        assert!(actual.observed_source_time_rate <= observed.maximum_source_time_rate());
        assert_eq!(pressure.receiver_time_s, t);
        assert!(pressure.retarded_time_residual_s <= 1e-13);
        for axis in 0..3 {
            near(
                actual.receiver_kinematics.position_m[axis],
                x[axis],
                2e-13,
                1.0,
            );
        }
    }
}

#[test]
fn stationary_microphone_recovers_the_existing_curved_source_owner() {
    let x = [8.0, 5.0, -3.0];
    let endpoints = samples(&[-0.3, -0.04, 0.11, 0.4], source);
    let stationary = MovingCompactObserver::new(
        medium(),
        x,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.1,
        },
    )
    .unwrap();
    let old = PreparedSourceTrajectory::new(stationary, &endpoints, path_budget(), &mut || false)
        .unwrap();
    let source_path =
        PreparedKinematicTrajectory::new(medium(), &endpoints, 0.8, path_budget(), &mut || false)
            .unwrap();
    let receiver_path = PreparedKinematicTrajectory::new(
        medium(),
        &samples(&[0.0, 0.2], |_| fixed(x)),
        0.0,
        path_budget(),
        &mut || false,
    )
    .unwrap();
    let observed = PreparedMovingReceiverSource::new(
        &source_path,
        &receiver_path,
        config([0.0, 0.2], 0.1),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    for t in [0.0, 0.03, 0.11, 0.2] {
        let old_sample = old.prepare_receiver(t, 1e-13, &mut || false).unwrap();
        let new_sample = observed.prepare_receiver(t, &mut || false).unwrap();
        let a = old_sample
            .evaluate(strength(old_sample.source_time_s()))
            .unwrap();
        let b = new_sample
            .pressure
            .evaluate(strength(new_sample.pressure.source_time_s()))
            .unwrap();
        near(a.source_time_s, b.source_time_s, 3e-11, 0.1);
        near(a.pressure_pa, b.pressure_pa, 3e-10, 1e-6);
        near(
            new_sample.observed_source_time_rate,
            old_sample.source_time_rate(),
            3e-12,
            1.0,
        );
    }
}

#[test]
fn moving_microphone_changes_phase_rate_without_a_spurious_pressure_gain() {
    let source_path = path(&[-0.1, 0.2], |_| SourceKinematics::default());
    for speed in [-60.0, 60.0] {
        let receiver_path = path(&[0.0, 0.1], |t| SourceKinematics {
            position_m: [12.0 + speed * t, 0.0, 0.0],
            velocity_m_s: [speed, 0.0, 0.0],
            acceleration_m_s2: [0.0; 3],
        });
        let observed = PreparedMovingReceiverSource::new(
            &source_path,
            &receiver_path,
            config([0.0, 0.1], 0.1),
            observation_budget(),
            &mut || false,
        )
        .unwrap();
        for t in [0.0, 0.03, 0.1] {
            let r = 12.0 + speed * t;
            let tau = (1.0 - speed / 300.0) * t - 12.0 / 300.0;
            let prepared = observed.prepare_receiver(t, &mut || false).unwrap();
            let qdot = 0.1 * (60.0 * tau).cos();
            let actual = prepared
                .pressure
                .evaluate(CompactSourceStrength {
                    volume_flow_m3_s: (0.1 / 60.0) * (60.0 * tau).sin(),
                    volume_flow_rate_m3_s2: qdot,
                    ..CompactSourceStrength::default()
                })
                .unwrap();
            near(actual.source_time_s, tau, 2e-12, 0.1);
            near(actual.pressure_pa, 1.2 * qdot / (4.0 * PI * r), 2e-12, 1e-8);
            near(
                prepared.observed_source_time_rate,
                1.0 - speed / 300.0,
                2e-13,
                1.0,
            );
            assert_eq!(actual.source_time_rate, 1.0);
        }
    }
}

#[test]
fn comoving_paths_preserve_source_convection_and_retarded_clearance() {
    let source_path = path(&[-0.1, 0.3], |t| SourceKinematics {
        position_m: [60.0 * t, 0.0, 0.0],
        velocity_m_s: [60.0, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    });
    let receiver_path = path(&[0.0, 0.2], |t| SourceKinematics {
        position_m: [12.0 + 60.0 * t, 0.0, 0.0],
        velocity_m_s: [60.0, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    });
    // Simultaneous separation is 12 m; every retarded separation is 15 m.
    // Full spatial path ranges overlap, so admission must prune pairs that
    // cannot lie on the acoustic light cone before judging their clearance.
    let observed = PreparedMovingReceiverSource::new(
        &source_path,
        &receiver_path,
        config([0.0, 0.2], 13.0),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    for t in [0.0, 0.07, 0.2] {
        let prepared = observed.prepare_receiver(t, &mut || false).unwrap();
        let actual = prepared
            .pressure
            .evaluate(CompactSourceStrength {
                volume_flow_rate_m3_s2: 2.0,
                ..CompactSourceStrength::default()
            })
            .unwrap();
        near(actual.source_time_s, t - 0.05, 3e-12, 0.1);
        near(actual.distance_m, 15.0, 3e-12, 1.0);
        near(prepared.observed_source_time_rate, 1.0, 3e-13, 1.0);
        near(actual.source_time_rate, 1.25, 3e-13, 1.0);
        near(
            actual.pressure_pa,
            1.2 * 2.0 / (4.0 * PI * 15.0 * 0.8_f64.powi(2)),
            3e-12,
            1e-6,
        );
    }
}

#[test]
fn endpoint_clearance_cannot_hide_a_receiver_crossing_the_source_between_samples() {
    let source_path = path(&[-0.1, 1.0], |_| SourceKinematics::default());
    let receiver_path = path(&[0.0, 1.0], |t| SourceKinematics {
        position_m: [(2.0 * t - 1.0).powi(2), 0.0, 0.0],
        velocity_m_s: [8.0 * t - 4.0, 0.0, 0.0],
        acceleration_m_s2: [8.0, 0.0, 0.0],
    });
    assert!(
        PreparedMovingReceiverSource::new(
            &source_path,
            &receiver_path,
            config([0.0, 1.0], 0.05),
            observation_budget(),
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn finite_source_and_receiver_domains_and_original_medium_are_required() {
    let source_path = path(&[-0.3, 0.4], source);
    let receiver_path = path(&[0.0, 0.2], microphone);
    let observed = PreparedMovingReceiverSource::new(
        &source_path,
        &receiver_path,
        config([0.0, 0.2], 0.1),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    for t in [0.0_f64.next_down(), 0.2_f64.next_up(), f64::NAN] {
        assert!(observed.prepare_receiver(t, &mut || false).is_err());
    }
    assert!(
        PreparedMovingReceiverSource::new(
            &source_path,
            &receiver_path,
            config([-0.01, 0.2], 0.1),
            observation_budget(),
            &mut || false,
        )
        .is_err()
    );
    let truncated = path(&[0.0, 0.4], source);
    assert!(
        PreparedMovingReceiverSource::new(
            &truncated,
            &receiver_path,
            config([0.0, 0.2], 0.1),
            observation_budget(),
            &mut || false,
        )
        .is_err()
    );
    let other_medium = PreparedKinematicTrajectory::new(
        FluidAcousticState::try_new(1.3, 300.0).unwrap(),
        &samples(&[0.0, 0.2], microphone),
        0.8,
        path_budget(),
        &mut || false,
    )
    .unwrap();
    assert!(
        PreparedMovingReceiverSource::new(
            &source_path,
            &other_medium,
            config([0.0, 0.2], 0.1),
            observation_budget(),
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn explicit_work_limits_and_cancellation_do_not_mutate_borrowed_paths() {
    let source_path = path(&[-0.3, -0.04, 0.11, 0.4], source);
    let receiver_path = path(&[0.0, 0.07, 0.2], microphone);
    for budget in [
        MovingReceiverBudget {
            maximum_segment_pairs: 0,
            ..observation_budget()
        },
        MovingReceiverBudget {
            maximum_bound_nodes: 0,
            ..observation_budget()
        },
        MovingReceiverBudget {
            maximum_root_iterations: 0,
            ..observation_budget()
        },
    ] {
        assert!(
            PreparedMovingReceiverSource::new(
                &source_path,
                &receiver_path,
                config([0.0, 0.2], 0.1),
                budget,
                &mut || false,
            )
            .is_err()
        );
    }
    assert!(matches!(
        PreparedMovingReceiverSource::new(
            &source_path,
            &receiver_path,
            config([0.0, 0.2], 0.1),
            observation_budget(),
            &mut || true,
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let observed = PreparedMovingReceiverSource::new(
        &source_path,
        &receiver_path,
        config([0.0, 0.2], 0.1),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    let t = 0.137;
    let source_before = source_path.kinematics_at(0.071).unwrap();
    let receiver_before = receiver_path.kinematics_at(t).unwrap();
    let mut polls = 0;
    let expected = observed
        .prepare_receiver(t, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 1);
    let mut calls = 0;
    let stop_at = (polls / 2).max(1);
    assert!(matches!(
        observed.prepare_receiver(t, &mut || {
            calls += 1;
            calls == stop_at
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let retry = observed.prepare_receiver(t, &mut || false).unwrap();
    assert_eq!(source_path.kinematics_at(0.071).unwrap(), source_before);
    assert_eq!(receiver_path.kinematics_at(t).unwrap(), receiver_before);
    assert_eq!(
        retry.pressure.source_time_s().to_bits(),
        expected.pressure.source_time_s().to_bits()
    );
    assert_eq!(
        retry.observed_source_time_rate.to_bits(),
        expected.observed_source_time_rate.to_bits()
    );
    assert_eq!(
        retry
            .pressure
            .evaluate(strength(retry.pressure.source_time_s()))
            .unwrap(),
        expected
            .pressure
            .evaluate(strength(expected.pressure.source_time_s()))
            .unwrap(),
    );
}
