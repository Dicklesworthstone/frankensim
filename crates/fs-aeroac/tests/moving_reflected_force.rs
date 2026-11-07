//! Full reflected Q/F fields preserve near-field velocity and coherent flux.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedReflectedMovingReceiverSource,
    PreparedReflectedSourceTrajectory, PreparedSourceTrajectory, ReflectedCompactField,
    ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}

fn observer(position: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium(),
        position,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.05,
        },
    )
    .unwrap()
}

fn budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 8,
        maximum_payload_bytes: 131_072,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}

fn bounds() -> ReflectedTrajectoryBudget {
    ReflectedTrajectoryBudget {
        maximum_payload_bytes: 1_000_000,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
    }
}

fn pairs() -> MovingReceiverBudget {
    MovingReceiverBudget {
        maximum_segment_pairs: 1024,
        maximum_bound_nodes: 65_536,
        maximum_subdivision_depth: 24,
        maximum_root_iterations: 96,
    }
}

fn config() -> MovingReceiverConfig {
    MovingReceiverConfig {
        receiver_time_interval_s: [0.0, 0.3],
        minimum_distance_m: 0.05,
        maximum_time_residual_s: 1e-13,
    }
}

fn knots(motion: impl Fn(f64) -> SourceKinematics) -> Vec<TimedSourceKinematics> {
    [-1.0, 0.0, 0.5]
        .into_iter()
        .map(|source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect()
}

fn fixed(position_m: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m,
        ..SourceKinematics::default()
    }
}

fn path(motion: impl Fn(f64) -> SourceKinematics) -> PreparedKinematicTrajectory {
    PreparedKinematicTrajectory::new(medium(), &knots(motion), 0.8, budget(), &mut || false)
        .unwrap()
}

fn close(actual: f64, expected: f64, floor: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * expected.abs().max(floor),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn drive(time: f64) -> CompactSourceStrength {
    // Smooth onset at t=-1; the force impulse below has zero initial value.
    let s = time + 1.0;
    CompactSourceStrength {
        volume_flow_m3_s: 0.002 * s * s,
        volume_flow_rate_m3_s2: 0.004 * s,
        force_on_fluid_n: [2.0 * s * s, -s * s * s, 0.7 * s * s],
        force_rate_on_fluid_n_s: [4.0 * s, -3.0 * s * s, 1.4 * s],
    }
}

fn impulse(time: f64) -> [f64; 3] {
    let s = time + 1.0;
    [
        2.0 * s.powi(3) / 3.0,
        -s.powi(4) / 4.0,
        0.7 * s.powi(3) / 3.0,
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

// Independent stationary Green-function derivatives: 1/r radiation, 1/r²
// induction and the 1/r³ force-impulse field, on each SOURCE image's clock.
fn stationary_field(
    x: [f64; 3],
    y: [f64; 3],
    time: f64,
    transform: impl Fn([f64; 3]) -> [f64; 3],
    sign: f64,
) -> (f64, [f64; 3]) {
    let d: [f64; 3] = core::array::from_fn(|i| x[i] - y[i]);
    let r = dot(d, d).sqrt();
    let n = d.map(|v| v / r);
    let tau = time - r / 300.0;
    let source = drive(tau);
    let force = transform(source.force_on_fluid_n);
    let rate = transform(source.force_rate_on_fluid_n_s);
    let history = transform(impulse(tau));
    let p = sign
        * (1.2 * source.volume_flow_rate_m3_s2 / r
            + dot(n, force) / r.powi(2)
            + dot(n, rate) / (300.0 * r))
        / (4.0 * PI);
    let u = core::array::from_fn(|i| {
        sign * (n[i]
            * (source.volume_flow_m3_s / r.powi(2) + source.volume_flow_rate_m3_s2 / (300.0 * r))
            + ((3.0 * n[i] * dot(n, history) - history[i]) / r.powi(3)
                + (3.0 * n[i] * dot(n, force) - force[i]) / (300.0 * r.powi(2))
                + n[i] * dot(n, rate) / (300.0_f64.powi(2) * r))
                / 1.2)
            / (4.0 * PI)
    });
    (p, u)
}

fn check_field(actual: ReflectedCompactField, pressure: f64, velocity: [f64; 3]) {
    close(actual.pressure_pa, pressure, 1e-5);
    for axis in 0..3 {
        close(actual.particle_velocity_m_s[axis], velocity[axis], 1e-7);
        close(
            actual.instantaneous_intensity_w_m2[axis],
            pressure * velocity[axis],
            1e-9,
        );
    }
    assert!(actual.maximum_time_residual_s <= 1e-13);
}

fn corner_vector(v: [f64; 3], mask: usize) -> [f64; 3] {
    core::array::from_fn(|i| if mask & (1 << i) == 0 { v[i] } else { -v[i] })
}

fn walls(parity: usize, count: usize) -> Vec<IdealAcousticPlane> {
    (0..count)
        .map(|i| {
            let mut normal = [0.0; 3];
            normal[i] = 1.0;
            IdealAcousticPlane::new(
                [0.0; 3],
                normal,
                if parity & (1 << i) == 0 {
                    IdealPlaneBoundary::Rigid
                } else {
                    IdealPlaneBoundary::PressureRelease
                },
            )
            .unwrap()
        })
        .collect()
}

#[test]
fn full_corner_fields_match_independent_source_images_and_keep_cross_path_flux() {
    let x = [7.0, 5.0, 4.0];
    let y = [1.0, 2.0, 3.0];
    let base =
        PreparedSourceTrajectory::new(observer(x), &knots(|_| fixed(y)), budget(), &mut || false)
            .unwrap();
    let mut saw_cross_path_flux = false;
    for parity in 0usize..8 {
        let images = IdealReflectionObserver::new(observer(x), &walls(parity, 3), 8).unwrap();
        let reflected =
            PreparedReflectedSourceTrajectory::new(&base, &images, bounds(), &mut || false)
                .unwrap();
        for time in [0.05, 0.137, 0.3] {
            let mut history_events = Vec::new();
            let actual = reflected
                .force_field_at(
                    time,
                    1e-13,
                    &mut |tau| Ok(drive(tau)),
                    &mut |index, virtual_observer, tau| {
                        history_events.push((index, virtual_observer.observer_position_m(), tau));
                        virtual_observer.force_history_integrand(y, impulse(tau))
                    },
                    &mut || false,
                )
                .unwrap();
            assert_eq!(actual.receiver_time_s, time);
            assert_eq!(history_events.len(), 8);
            let mut p = 0.0;
            let mut u = [0.0; 3];
            let mut separate_flux = [0.0; 3];
            for mask in 0usize..8 {
                let sign = if (parity & mask).count_ones() % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                let (path_p, path_u) = stationary_field(
                    x,
                    corner_vector(y, mask),
                    time,
                    |v| corner_vector(v, mask),
                    sign,
                );
                p += path_p;
                for axis in 0..3 {
                    u[axis] += path_u[axis];
                    separate_flux[axis] += path_p * path_u[axis];
                }
                assert_eq!(history_events[mask].0, mask);
                assert_eq!(history_events[mask].1, corner_vector(x, mask));
            }
            assert_ne!(history_events[0].2, history_events[7].2);
            check_field(actual, p, u);
            for axis in 0..3 {
                saw_cross_path_flux |=
                    (p * u[axis] - separate_flux[axis]).abs() > 0.1 * (p * u[axis]).abs();
            }
            close(
                actual.pressure_pa,
                reflected
                    .pressure_at(time, 1e-13, &mut |tau| Ok(drive(tau)), &mut || false)
                    .unwrap(),
                1e-5,
            );
        }
    }
    assert!(
        saw_cross_path_flux,
        "the full field must retain cross-path energy transport"
    );
}

fn wedge_vector(v: [f64; 3], index: usize, order: usize) -> [f64; 3] {
    let (s, c) = (2.0 * PI * (index / 2) as f64 / order as f64).sin_cos();
    let y = if index % 2 == 0 { v[1] } else { -v[1] };
    [c * v[0] - s * y, s * v[0] + c * y, v[2]]
}

#[test]
fn moving_microphones_retain_complete_wedge_force_vectors_and_original_clock() {
    let y = [3.0, 0.5, 1.0];
    let receiver = |t: f64| SourceKinematics {
        position_m: [8.0 - 2.0 * t + t * t, 1.1 + 0.3 * t, 3.0 + 0.2 * t],
        velocity_m_s: [-2.0 + 2.0 * t, 0.3, 0.2],
        acceleration_m_s2: [2.0, 0.0, 0.0],
    };
    let source_path = path(|_| fixed(y));
    let receiver_path = path(receiver);
    for (order, parity) in [(3usize, 0usize), (3, 3), (4, 1), (4, 2)] {
        let (s, c) = (PI / order as f64).sin_cos();
        let planes = core::array::from_fn(|i| {
            IdealAcousticPlane::new(
                [0.0; 3],
                if i == 0 {
                    [0.0, 1.0, 0.0]
                } else {
                    [s, -c, 0.0]
                },
                if parity & (1 << i) == 0 {
                    IdealPlaneBoundary::Rigid
                } else {
                    IdealPlaneBoundary::PressureRelease
                },
            )
            .unwrap()
        });
        let images = MovingReflectionObserver::new_wedge(
            &receiver_path,
            planes,
            order,
            2 * order,
            bounds(),
            &mut || false,
        )
        .unwrap();
        let reflected = PreparedReflectedMovingReceiverSource::new(
            &source_path,
            &images,
            config(),
            bounds(),
            pairs(),
            &mut || false,
        )
        .unwrap();
        for time in [0.05, 0.137, 0.3] {
            let actual = reflected
                .force_field_at(
                    time,
                    &mut |tau| Ok(drive(tau)),
                    &mut |_, virtual_observer, tau| {
                        virtual_observer.force_history_integrand(y, impulse(tau))
                    },
                    &mut || false,
                )
                .unwrap();
            let x = receiver(time).position_m;
            let first_sign = if parity & 1 == 0 { 1.0_f64 } else { -1.0 };
            let second_sign = if parity & 2 == 0 { 1.0_f64 } else { -1.0 };
            let mut p = 0.0;
            let mut u = [0.0; 3];
            for index in 0..2 * order {
                let sign = (first_sign * second_sign).powi((index / 2) as i32)
                    * if index % 2 == 0 { 1.0 } else { first_sign };
                let (path_p, path_u) = stationary_field(
                    x,
                    wedge_vector(y, index, order),
                    time,
                    |v| wedge_vector(v, index, order),
                    sign,
                );
                p += path_p;
                for axis in 0..3 {
                    u[axis] += path_u[axis];
                }
            }
            check_field(actual, p, u);
        }
    }
}

fn accelerating_source(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            2.0 + 0.2 * t + 0.1 * t * t,
            0.5 + 0.3 * t,
            1.0 + 0.05 * t * t,
        ],
        velocity_m_s: [0.2 + 0.2 * t, 0.3, 0.1 * t],
        acceleration_m_s2: [0.2, 0.0, 0.1],
    }
}

fn past_force(observer: &MovingCompactObserver, tau: f64) -> [f64; 3] {
    // Finite original history from t=-1, at this exact virtual microphone.
    // Boundary symmetry is independent of quadrature error; the same explicit
    // history is also supplied to the independent free-space field owner.
    const PANELS: usize = 128;
    let h = (tau + 1.0) / PANELS as f64;
    let mut sum = [0.0; 3];
    for i in 0..=PANELS {
        let t = -1.0 + h * i as f64;
        let value = observer
            .force_history_integrand(accelerating_source(t).position_m, drive(t).force_on_fluid_n)
            .unwrap();
        let weight = if i == 0 || i == PANELS {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for axis in 0..3 {
            sum[axis] += weight * value[axis];
        }
    }
    sum.map(|v| v * h / 3.0)
}

#[test]
fn accelerated_force_fields_enforce_walls_and_keep_motion_at_pressure_nodes() {
    let source_path = path(accelerating_source);
    let receiver = |t: f64| SourceKinematics {
        position_m: [0.0, 4.0 + 0.5 * t, 3.0 - 0.2 * t],
        velocity_m_s: [0.0, 0.5, -0.2],
        acceleration_m_s2: [0.0; 3],
    };
    let receiver_path = path(receiver);
    for parity in 0..2 {
        let planes = walls(parity, 1);
        let images =
            MovingReflectionObserver::new(&receiver_path, &planes, 2, bounds(), &mut || false)
                .unwrap();
        let reflected = PreparedReflectedMovingReceiverSource::new(
            &source_path,
            &images,
            config(),
            bounds(),
            pairs(),
            &mut || false,
        )
        .unwrap();
        for time in [0.05, 0.137, 0.3] {
            let actual = reflected
                .force_field_at(
                    time,
                    &mut |tau| Ok(drive(tau)),
                    &mut |_, microphone, tau| Ok(past_force(microphone, tau)),
                    &mut || false,
                )
                .unwrap();
            let microphone = observer(receiver(time).position_m);
            let tau = reflected
                .prepare_path(0, time, &mut || false)
                .unwrap()
                .pressure
                .source_time_s();
            let direct = microphone
                .prepare_force_emission(tau, accelerating_source(tau))
                .unwrap()
                .evaluate(drive(tau), past_force(&microphone, tau))
                .unwrap();
            let expected_p = if parity == 0 {
                2.0 * direct.sample.pressure_pa
            } else {
                0.0
            };
            let expected_u = core::array::from_fn(|i| {
                if (parity == 0 && i != 0) || (parity == 1 && i == 0) {
                    2.0 * direct.particle_velocity_m_s[i]
                } else {
                    0.0
                }
            });
            check_field(actual, expected_p, expected_u);
            assert_eq!(actual.instantaneous_intensity_w_m2[0], 0.0);
            if parity == 1 {
                assert_eq!(actual.pressure_pa, 0.0);
                assert!(
                    actual.particle_velocity_m_s[0].abs() > 1e-6,
                    "a pressure node retains the normal particle motion"
                );
            }
            // A moving microphone samples the same Eulerian field as a fixed
            // one at its instantaneous location, without an amplitude Doppler gain.
            let base = PreparedSourceTrajectory::new(
                microphone,
                &knots(accelerating_source),
                budget(),
                &mut || false,
            )
            .unwrap();
            let fixed_images = IdealReflectionObserver::new(microphone, &planes, 2).unwrap();
            let fixed_paths =
                PreparedReflectedSourceTrajectory::new(&base, &fixed_images, bounds(), &mut || {
                    false
                })
                .unwrap();
            let fixed_field = fixed_paths
                .force_field_at(
                    time,
                    1e-13,
                    &mut |tau| Ok(drive(tau)),
                    &mut |_, microphone, tau| Ok(past_force(microphone, tau)),
                    &mut || false,
                )
                .unwrap();
            check_field(
                fixed_field,
                actual.pressure_pa,
                actual.particle_velocity_m_s,
            );
        }
    }
}

#[test]
fn incomplete_or_nonfinite_history_and_cancellation_refuse_the_whole_field() {
    let x = [3.0, 0.0, 0.0];
    let y = [1.0, 0.0, 0.0];
    let base =
        PreparedSourceTrajectory::new(observer(x), &knots(|_| fixed(y)), budget(), &mut || false)
            .unwrap();
    let images = IdealReflectionObserver::new(observer(x), &walls(0, 1), 2).unwrap();
    let reflected =
        PreparedReflectedSourceTrajectory::new(&base, &images, bounds(), &mut || false).unwrap();
    let mut callbacks = 0;
    assert!(matches!(
        reflected.force_field_at(
            0.1,
            1e-13,
            &mut |t| Ok(drive(t)),
            &mut |_, _, _| {
                callbacks += 1;
                Ok([0.0; 3])
            },
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(callbacks, 0);
    assert!(
        reflected
            .force_field_at(
                2.0,
                1e-13,
                &mut |t| Ok(drive(t)),
                &mut |_, _, _| {
                    callbacks += 1;
                    Ok([0.0; 3])
                },
                &mut || false
            )
            .is_err()
    );
    assert_eq!(callbacks, 0);
    assert!(matches!(
        reflected.force_field_at(
            0.1,
            1e-13,
            &mut |t| Ok(drive(t)),
            &mut |index, _, _| {
                callbacks += 1;
                Ok(if index == 0 {
                    [0.0; 3]
                } else {
                    [f64::NAN, 0.0, 0.0]
                })
            },
            &mut || false
        ),
        Err(AeroacError::NonFinite { .. })
    ));
    assert_eq!(callbacks, 2);
    let missing = AeroacError::InvalidParameter {
        what: "reflected force history is unavailable",
    };
    assert_eq!(
        reflected
            .force_field_at(
                0.1,
                1e-13,
                &mut |t| Ok(drive(t)),
                &mut |index, _, _| if index == 0 {
                    Ok([0.0; 3])
                } else {
                    Err(missing.clone())
                },
                &mut || false
            )
            .unwrap_err(),
        missing
    );
}
