//! G1/G3 finite accelerating trajectories, continuous admission and causal clocks.

#[path = "moving_trajectory/reflection.rs"]
mod reflection;

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    PreparedSourceTrajectory, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    UniformSourceMotion,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

fn observer(
    position: [f64; 3],
    c: f64,
    maximum_mach: f64,
    clearance: f64,
) -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, c).unwrap(),
        position,
        MovingCompactLimits {
            maximum_mach,
            minimum_distance_m: clearance,
        },
    )
    .unwrap()
}

fn budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 64,
        maximum_payload_bytes: 1_000_000,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}

fn near(actual: f64, expected: f64, relative: f64, floor: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(floor),
        "actual={actual:e}, expected={expected:e}"
    );
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

fn prepare(
    observer: MovingCompactObserver,
    samples: &[TimedSourceKinematics],
) -> PreparedSourceTrajectory {
    PreparedSourceTrajectory::new(observer, samples, budget(), &mut || false).unwrap()
}

fn strength(t: f64) -> CompactSourceStrength {
    let (s, c) = (25.0 * t).sin_cos();
    let (s2, c2) = (17.0 * t).sin_cos();
    CompactSourceStrength {
        volume_flow_m3_s: 0.003 * s + 0.001 * c2,
        volume_flow_rate_m3_s2: 0.075 * c - 0.017 * s2,
        force_on_fluid_n: [2.0 * s, -3.0 * c, 0.5 * s2],
        force_rate_on_fluid_n_s: [50.0 * c, 75.0 * s, 8.5 * c2],
    }
}

#[test]
fn stationary_and_uniform_paths_recover_the_existing_analytic_observer() {
    for velocity in [[0.0; 3], [45.0, -10.0, 4.0], [-30.0, 8.0, -2.0]] {
        let stationary = velocity == [0.0; 3];
        let microphone = observer(
            [8.0, 3.0, 2.0],
            300.0,
            if stationary { 0.0 } else { 0.8 },
            0.05,
        );
        let uniform = UniformSourceMotion::new(0.0, [0.2, -0.1, 0.4], velocity).unwrap();
        let knots = samples(&[-0.1, 0.1, 0.3], |t| uniform.kinematics_at(t).unwrap());
        let trajectory = prepare(microphone, &knots);
        assert_eq!(trajectory.segment_count(), 2);
        assert_eq!(trajectory.source_time_interval_s(), [-0.1, 0.3]);
        assert_eq!(trajectory.observer().medium(), microphone.medium());
        if stationary {
            near(trajectory.maximum_source_time_rate(), 1.0, 2e-14, 1.0);
        }
        for receiver_time in [0.01, 0.08, 0.17] {
            let analytic = microphone
                .prepare_uniform(receiver_time, &uniform, 1e-12)
                .unwrap();
            let solved = trajectory
                .prepare_receiver(receiver_time, 1e-12, &mut || false)
                .unwrap();
            near(solved.source_time_s(), analytic.source_time_s(), 1e-10, 0.1);
            near(solved.distance_m(), analytic.distance_m(), 1e-10, 1.0);
            assert!(solved.source_time_rate() <= trajectory.maximum_source_time_rate());
            let actual = solved.evaluate(strength(solved.source_time_s())).unwrap();
            let expected = analytic
                .evaluate(strength(analytic.source_time_s()))
                .unwrap();
            near(actual.pressure_pa, expected.pressure_pa, 1e-9, 1e-6);
            assert_eq!(actual.receiver_time_s, receiver_time);
            assert!(actual.retarded_time_residual_s <= 1e-12);
        }
    }
}

// A general degree-five world-frame trajectory with independently written
// first and second derivatives. Any knot partition represents this same curve.
fn polynomial(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            0.2 + t * (10.0 + t * (20.0 + t * (-15.0 + t * (7.0 + 9.0 * t)))),
            -0.3 + t * (-6.0 + t * (11.0 + t * (8.0 + t * (-3.0 + 2.0 * t)))),
            0.1 + t * (3.0 + t * (-4.0 + t * (3.0 + t * (5.0 - 6.0 * t)))),
        ],
        velocity_m_s: [
            10.0 + t * (40.0 + t * (-45.0 + t * (28.0 + 45.0 * t))),
            -6.0 + t * (22.0 + t * (24.0 + t * (-12.0 + 10.0 * t))),
            3.0 + t * (-8.0 + t * (9.0 + t * (20.0 - 30.0 * t))),
        ],
        acceleration_m_s2: [
            40.0 + t * (-90.0 + t * (84.0 + 180.0 * t)),
            22.0 + t * (48.0 + t * (-36.0 + 40.0 * t)),
            -8.0 + t * (18.0 + t * (60.0 - 120.0 * t)),
        ],
    }
}

#[test]
fn one_quintic_supplies_position_velocity_and_acceleration_on_irregular_knots() {
    let knots = samples(&[-0.2, -0.03, 0.11, 0.4], polynomial);
    let trajectory = prepare(observer([8.0, 5.0, 3.0], 300.0, 0.8, 0.05), &knots);
    for endpoint in &knots {
        assert_eq!(
            trajectory.kinematics_at(endpoint.source_time_s).unwrap(),
            endpoint.kinematics
        );
    }
    for i in 0..=120 {
        let t = -0.2 + f64::from(i) * 0.005;
        let actual = trajectory.kinematics_at(t.min(0.4)).unwrap();
        let expected = polynomial(t.min(0.4));
        for axis in 0..3 {
            near(
                actual.position_m[axis],
                expected.position_m[axis],
                2e-11,
                1.0,
            );
            near(
                actual.velocity_m_s[axis],
                expected.velocity_m_s[axis],
                2e-11,
                1.0,
            );
            near(
                actual.acceleration_m_s2[axis],
                expected.acceleration_m_s2[axis],
                2e-10,
                1.0,
            );
        }
    }
}

#[test]
fn known_emission_events_recover_causal_times_and_accelerating_pressure() {
    let microphone = observer([8.0, 5.0, 3.0], 300.0, 0.8, 0.05);
    let knots = samples(&[-0.2, -0.03, 0.11, 0.4], polynomial);
    let trajectory = prepare(microphone, &knots);
    for tau in [-0.17, 0.0, 0.07, 0.299] {
        let emitted = microphone.prepare_emission(tau, polynomial(tau)).unwrap();
        let solved = trajectory
            .prepare_receiver(emitted.receiver_time_s(), 1e-13, &mut || false)
            .unwrap();
        near(solved.source_time_s(), tau, 2e-11, 0.1);
        let actual = solved.evaluate(strength(solved.source_time_s())).unwrap();
        let expected = emitted.evaluate(strength(tau)).unwrap();
        near(
            actual.monopole_pressure_pa,
            expected.monopole_pressure_pa,
            2e-9,
            1e-7,
        );
        near(
            actual.dipole_pressure_pa,
            expected.dipole_pressure_pa,
            2e-9,
            1e-6,
        );
        near(actual.pressure_pa, expected.pressure_pa, 2e-9, 1e-6);
        assert!(actual.source_time_s < actual.receiver_time_s);
        assert!(actual.retarded_time_residual_s <= 1e-13);
    }
}

#[test]
fn zero_endpoint_speeds_cannot_hide_an_interior_sonic_excursion() {
    let microphone = observer([0.0, 10.0, 0.0], 4.0, 0.8, 0.1);
    let knots = |distance| {
        [
            TimedSourceKinematics {
                source_time_s: 0.0,
                kinematics: SourceKinematics::default(),
            },
            TimedSourceKinematics {
                source_time_s: 1.0,
                kinematics: SourceKinematics {
                    position_m: [distance, 0.0, 0.0],
                    ..SourceKinematics::default()
                },
            },
        ]
    };
    // Quintic smoothstep has maximum velocity 1.875*distance/duration.
    // Both endpoint states have exactly zero velocity and acceleration.
    let sonic = knots(3.0);
    assert!(sonic.iter().all(|s| {
        microphone
            .prepare_emission(s.source_time_s, s.kinematics)
            .is_ok()
    }));
    assert!(PreparedSourceTrajectory::new(microphone, &sonic, budget(), &mut || false).is_err());
    let admitted = prepare(microphone, &knots(1.0));
    near(
        admitted.kinematics_at(0.5).unwrap().velocity_m_s[0],
        1.875,
        2e-12,
        1.0,
    );
    assert!(admitted.maximum_source_time_rate() >= 1.0 / (1.0 - 1.875 / 4.0));
}

fn parabola(t: f64, impact: f64) -> SourceKinematics {
    let x = 2.0 * t - 1.0;
    SourceKinematics {
        position_m: [x, impact + x * x, 0.0],
        velocity_m_s: [2.0, 4.0 * x, 0.0],
        acceleration_m_s2: [0.0, 8.0, 0.0],
    }
}

#[test]
fn continuous_curved_clearance_is_resolved_by_bounded_hull_subdivision() {
    let microphone = observer([0.0; 3], 100.0, 0.8, 0.1);
    let near_pass = samples(&[0.0, 1.0], |t| parabola(t, 0.05));
    assert!(near_pass.iter().all(|s| {
        microphone
            .prepare_emission(s.source_time_s, s.kinematics)
            .is_ok()
    }));
    assert!(
        PreparedSourceTrajectory::new(microphone, &near_pass, budget(), &mut || false).is_err()
    );
    let clear_pass = samples(&[0.0, 1.0], |t| parabola(t, 0.25));
    // The unsplit degree-five position hull extends to y=0.05, whereas
    // the true path has y>=0.25. Refusal with a one-node budget is honest;
    // subdivision must then recover this valid curved path.
    assert!(
        PreparedSourceTrajectory::new(
            microphone,
            &clear_pass,
            TrajectoryBudget {
                maximum_bound_nodes: 1,
                ..budget()
            },
            &mut || false
        )
        .is_err()
    );
    let trajectory = prepare(microphone, &clear_pass);
    let closest = trajectory.kinematics_at(0.5).unwrap();
    near(closest.position_m[0], 0.0, 1e-12, 1.0);
    near(closest.position_m[1], 0.25, 1e-12, 1.0);
}

fn circle(t: f64) -> SourceKinematics {
    let (s, c) = (20.0 * t).sin_cos();
    SourceKinematics {
        position_m: [2.0 * c, 2.0 * s, 0.0],
        velocity_m_s: [-40.0 * s, 40.0 * c, 0.0],
        acceleration_m_s2: [-800.0 * c, -800.0 * s, 0.0],
    }
}

#[test]
fn circular_motion_refinement_recovers_center_pressure_cancellation() {
    // At the center of an exact constant-speed circle, R and kappa are
    // constant. Constant Q therefore has zero pressure: acceleration and
    // convective terms cancel. The sampled quintic approximation converges
    // toward that independent physical limit as its angular step decreases.
    let microphone = observer([0.0; 3], 100.0, 0.8, 0.5);
    let strength = CompactSourceStrength {
        volume_flow_m3_s: 0.01,
        ..CompactSourceStrength::default()
    };
    let mut errors = Vec::new();
    for count in [2_u32, 4, 8] {
        let times: Vec<_> = (0..=count)
            .map(|i| f64::from(i) * PI / (20.0 * f64::from(count)))
            .collect();
        let trajectory = prepare(microphone, &samples(&times, circle));
        let mut maximum = 0.0_f64;
        for interval in times.windows(2) {
            for u in [0.23, 0.37, 0.61, 0.79] {
                let tau = interval[0] + u * (interval[1] - interval[0]);
                let emitted = microphone
                    .prepare_emission(tau, trajectory.kinematics_at(tau).unwrap())
                    .unwrap();
                let sample = trajectory
                    .prepare_receiver(emitted.receiver_time_s(), 1e-13, &mut || false)
                    .unwrap()
                    .evaluate(strength)
                    .unwrap();
                maximum = maximum.max(sample.pressure_pa.abs());
            }
        }
        errors.push(maximum);
    }
    assert!(
        errors[0] > 1e-6,
        "coarse sampled-circle error must be visible: {errors:?}"
    );
    assert!(
        errors[1] < errors[0] / 8.0 && errors[2] < errors[1] / 8.0,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-7, "{errors:?}");
}

fn distance(x: [f64; 3], y: [f64; 3]) -> f64 {
    ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt()
}

fn independent_time(x: [f64; 3], receiver: f64) -> f64 {
    let mut lo = -0.2;
    let mut hi = 0.4;
    assert!(lo + distance(x, polynomial(lo).position_m) / 300.0 < receiver);
    assert!(hi + distance(x, polynomial(hi).position_m) / 300.0 > receiver);
    for _ in 0..100 {
        let middle = 0.5 * (lo + hi);
        if middle + distance(x, polynomial(middle).position_m) / 300.0 < receiver {
            lo = middle;
        } else {
            hi = middle;
        }
    }
    0.5 * (lo + hi)
}

fn potential(x: [f64; 3], receiver: f64) -> [f64; 4] {
    let tau = independent_time(x, receiver);
    let state = polynomial(tau);
    let r = distance(x, state.position_m);
    let radial: f64 = (0..3)
        .map(|i| (x[i] - state.position_m[i]) * state.velocity_m_s[i] / r)
        .sum();
    let denominator = 4.0 * PI * r * (1.0 - radial / 300.0);
    let source = strength(tau);
    [
        source.volume_flow_m3_s / denominator,
        source.force_on_fluid_n[0] / denominator,
        source.force_on_fluid_n[1] / denominator,
        source.force_on_fluid_n[2] / denominator,
    ]
}

fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
}

#[test]
fn curved_receiver_pressure_matches_independently_differentiated_retarded_potentials() {
    let x = [8.0, 5.0, 3.0];
    let trajectory = prepare(
        observer(x, 300.0, 0.8, 0.05),
        &samples(&[-0.2, -0.03, 0.11, 0.4], polynomial),
    );
    for receiver in [-0.05, 0.05, 0.2, 0.35] {
        let prepared = trajectory
            .prepare_receiver(receiver, 1e-13, &mut || false)
            .unwrap();
        near(
            prepared.source_time_s(),
            independent_time(x, receiver),
            2e-11,
            0.1,
        );
        let sample = prepared
            .evaluate(strength(prepared.source_time_s()))
            .unwrap();
        let monopole = 1.2 * derivative(|dt| potential(x, receiver + dt)[0], 2e-5);
        let mut dipole = 0.0;
        for axis in 0..3 {
            dipole -= derivative(
                |dx| {
                    let mut shifted = x;
                    shifted[axis] += dx;
                    potential(shifted, receiver)[axis + 1]
                },
                2e-3,
            );
        }
        near(sample.monopole_pressure_pa, monopole, 3e-7, 1e-8);
        near(sample.dipole_pressure_pa, dipole, 3e-7, 1e-8);
        near(sample.pressure_pa, monopole + dipole, 3e-7, 1e-8);
    }
}

#[test]
fn finite_history_work_refusals_and_cancellation_leave_the_trajectory_reusable() {
    let microphone = observer([8.0, 5.0, 3.0], 300.0, 0.8, 0.05);
    let knots = samples(&[-0.2, -0.03, 0.11, 0.4], polynomial);
    for limits in [
        TrajectoryBudget {
            maximum_segments: 2,
            ..budget()
        },
        TrajectoryBudget {
            maximum_payload_bytes: 1,
            ..budget()
        },
        TrajectoryBudget {
            maximum_bound_nodes: 0,
            ..budget()
        },
        TrajectoryBudget {
            maximum_root_iterations: 0,
            ..budget()
        },
    ] {
        assert!(PreparedSourceTrajectory::new(microphone, &knots, limits, &mut || false).is_err());
    }
    assert!(matches!(
        PreparedSourceTrajectory::new(microphone, &knots, budget(), &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let trajectory = prepare(microphone, &knots);
    let [first, last] = trajectory.receiver_time_interval_s();
    for t in [first.next_down(), last.next_up(), f64::NAN] {
        assert!(
            trajectory
                .prepare_receiver(t, 1e-12, &mut || false)
                .is_err()
        );
    }
    assert!(trajectory.kinematics_at((-0.2_f64).next_down()).is_err());
    assert!(trajectory.kinematics_at(0.4_f64.next_up()).is_err());
    let receiver = microphone
        .prepare_emission(0.137, polynomial(0.137))
        .unwrap()
        .receiver_time_s();
    let short = PreparedSourceTrajectory::new(
        microphone,
        &knots,
        TrajectoryBudget {
            maximum_root_iterations: 1,
            ..budget()
        },
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        short.prepare_receiver(receiver, 1e-14, &mut || false),
        Err(AeroacError::NotConverged { .. })
    ));
    let mut polls = 0;
    let expected = trajectory
        .prepare_receiver(receiver, 1e-13, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 1);
    let stop_at = (polls / 2).max(1);
    let mut calls = 0;
    assert!(matches!(
        trajectory.prepare_receiver(receiver, 1e-13, &mut || {
            calls += 1;
            calls == stop_at
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let retried = trajectory
        .prepare_receiver(receiver, 1e-13, &mut || false)
        .unwrap();
    assert_eq!(
        retried.evaluate(strength(retried.source_time_s())).unwrap(),
        expected
            .evaluate(strength(expected.source_time_s()))
            .unwrap()
    );
}

#[test]
fn large_epoch_inversion_uses_the_same_public_source_time_for_all_kinematics() {
    let epoch = 1e12;
    let native = |time: f64| {
        let s = time - epoch;
        SourceKinematics {
            position_m: [1.0 + 0.1 * s + 0.2 * s * s, 0.2 - 0.1 * s, 0.0],
            velocity_m_s: [0.1 + 0.4 * s, -0.1, 0.0],
            acceleration_m_s2: [0.4, 0.0, 0.0],
        }
    };
    let microphone = observer([10.0, 3.0, 0.0], 10.0, 0.8, 0.1);
    let trajectory = prepare(microphone, &samples(&[epoch, epoch + 1.0], native));
    let tau = epoch + 0.375;
    let arrival = microphone
        .prepare_emission(tau, native(tau))
        .unwrap()
        .receiver_time_s();
    let prepared = trajectory
        .prepare_receiver(arrival, 1e-4, &mut || false)
        .unwrap();
    assert!((prepared.source_time_s() - tau).abs() <= 2.0 * (tau.next_up() - tau));
    let same_time_state = trajectory.kinematics_at(prepared.source_time_s()).unwrap();
    let same_time_source = strength(prepared.source_time_s() - epoch);
    let actual = prepared.evaluate(same_time_source).unwrap();
    let direct = microphone
        .prepare_emission(prepared.source_time_s(), same_time_state)
        .unwrap()
        .evaluate(same_time_source)
        .unwrap();
    assert_eq!(actual.pressure_pa.to_bits(), direct.pressure_pa.to_bits());
    assert_eq!(actual.distance_m.to_bits(), direct.distance_m.to_bits());
    assert_eq!(actual.receiver_time_s, arrival);
    assert!(actual.retarded_time_residual_s <= 1e-4);
}
