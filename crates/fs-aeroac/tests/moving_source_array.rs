//! Independently retarded moving source arrays retain full coherent fields.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver, MovingSourceArrayBudget,
    MovingSourceArrayField, PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedMovingSourceArray, PreparedReflectedMovingReceiverSource, ReflectedTrajectoryBudget,
    TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.4, 20.0).unwrap()
}

fn path(motion: impl Fn(f64) -> SourceKinematics) -> PreparedKinematicTrajectory {
    let knots: Vec<_> = [-1.0, 0.1, 0.8]
        .into_iter()
        .map(|source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect();
    PreparedKinematicTrajectory::new(
        medium(),
        &knots,
        0.8,
        TrajectoryBudget {
            maximum_segments: 8,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 16_384,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

fn pairs() -> MovingReceiverBudget {
    MovingReceiverBudget {
        maximum_segment_pairs: 1024,
        maximum_bound_nodes: 65_536,
        maximum_subdivision_depth: 24,
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

fn budget() -> MovingSourceArrayBudget {
    MovingSourceArrayBudget {
        maximum_sources: 8,
        maximum_source_paths: 64,
        maximum_retarded_iterations: 16_384,
    }
}

fn config(window: [f64; 2]) -> MovingReceiverConfig {
    MovingReceiverConfig {
        receiver_time_interval_s: window,
        minimum_distance_m: 0.05,
        maximum_time_residual_s: 1e-12,
    }
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

fn fixed(position_m: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m,
        ..SourceKinematics::default()
    }
}

fn motion(index: usize, t: f64) -> SourceKinematics {
    let position = [[1.0, 0.8, 0.5], [2.0, 1.4, 1.0], [1.5, 0.6, 1.7]][index];
    let velocity = [[0.3, 0.2, 0.0], [-0.25, 0.0, 0.2], [0.2, -0.1, -0.15]][index];
    let acceleration = [[0.24, 0.0, 0.2], [0.2, 0.3, 0.0], [0.0, 0.1, 0.0]][index];
    SourceKinematics {
        position_m: core::array::from_fn(|i| {
            position[i] + velocity[i] * t + 0.5 * acceleration[i] * t * t
        }),
        velocity_m_s: core::array::from_fn(|i| velocity[i] + acceleration[i] * t),
        acceleration_m_s2: acceleration,
    }
}

fn receiver(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [5.5 + 0.8 * t - 0.2 * t * t, 4.0 + 0.3 * t, 2.5 - 0.15 * t],
        velocity_m_s: [0.8 - 0.4 * t, 0.3, -0.15],
        acceleration_m_s2: [-0.4, 0.0, 0.0],
    }
}

fn drive(index: usize, t: f64) -> CompactSourceStrength {
    let s = t + 1.0;
    let q = [0.025, -0.017, 0.011][index];
    let a = [[2.0, -1.0, 0.3], [-0.6, 1.3, -0.8], [0.2, 0.5, 1.1]][index];
    let b = [[0.2, 0.1, -0.4], [0.3, -0.2, 0.1], [-0.1, 0.4, 0.2]][index];
    CompactSourceStrength {
        volume_flow_m3_s: q * s * s,
        volume_flow_rate_m3_s2: 2.0 * q * s,
        force_on_fluid_n: core::array::from_fn(|i| s * s * (a[i] + b[i] * s)),
        force_rate_on_fluid_n_s: core::array::from_fn(|i| s * (2.0 * a[i] + 3.0 * b[i] * s)),
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn image(v: [f64; 3], mask: usize) -> [f64; 3] {
    core::array::from_fn(|i| if mask & (1 << i) == 0 { v[i] } else { -v[i] })
}
fn sign(mask: usize, parity: usize) -> f64 {
    if (mask & parity).count_ones() % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}
fn walls(count: usize, parity: usize) -> Vec<IdealAcousticPlane> {
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

fn integrate(end: f64, f: impl Fn(f64) -> [f64; 3]) -> [f64; 3] {
    const PANELS: usize = 256;
    let h = (end + 1.0) / PANELS as f64;
    let mut sum = [0.0; 3];
    for i in 0..=PANELS {
        let value = f(-1.0 + h * i as f64);
        let weight = if i == 0 || i == PANELS {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for j in 0..3 {
            sum[j] += weight * value[j];
        }
    }
    sum.map(|v| v * h / 3.0)
}

fn history(index: usize, microphone: &MovingCompactObserver, tau: f64) -> [f64; 3] {
    integrate(tau, |t| {
        microphone
            .force_history_integrand(
                motion(index, t).position_m,
                drive(index, t).force_on_fluid_n,
            )
            .unwrap()
    })
}

// Independent SOURCE images, bisection roots and Green potentials. No native
// pressure coefficients, velocity formula or reflected-owner field is used.
// phi=-Q/(4*pi*R*kappa), A=F/(R*kappa), B=integral F(s)/R(s) ds.
fn potentials(x: [f64; 3], time: f64, count: usize, parity: usize) -> (f64, [f64; 3], [f64; 3]) {
    let mut phi = 0.0;
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    for source in 0..3 {
        for mask in 0..1 << count {
            let distance = |t: f64| {
                let y = image(motion(source, t).position_m, mask);
                norm(core::array::from_fn(|i| x[i] - y[i]))
            };
            let mut lo = -1.0;
            let mut hi = 0.8;
            assert!(lo + distance(lo) / 20.0 < time);
            assert!(hi + distance(hi) / 20.0 > time);
            for _ in 0..64 {
                let mid = 0.5 * (lo + hi);
                if mid + distance(mid) / 20.0 < time {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let tau = 0.5 * (lo + hi);
            let kinematics = motion(source, tau);
            let y = image(kinematics.position_m, mask);
            let velocity = image(kinematics.velocity_m_s, mask);
            let r = distance(tau);
            let n = core::array::from_fn(|i| (x[i] - y[i]) / r);
            let kappa = 1.0 - dot(n, velocity) / 20.0;
            let strength = drive(source, tau);
            let force = image(strength.force_on_fluid_n, mask);
            let impulse = integrate(tau, |t| {
                let f = image(drive(source, t).force_on_fluid_n, mask);
                f.map(|value| value / distance(t))
            });
            let parity = sign(mask, parity);
            phi -= parity * strength.volume_flow_m3_s / (4.0 * PI * r * kappa);
            for i in 0..3 {
                a[i] += parity * force[i] / (r * kappa);
                b[i] += parity * impulse[i];
            }
        }
    }
    (phi, a, b)
}

fn reference(x: [f64; 3], time: f64, count: usize, parity: usize, h: f64) -> (f64, [f64; 3]) {
    // Eulerian derivatives keep x fixed even when the physical receiver moves.
    let dt = h / 20.0;
    let mut p = -1.4
        * (potentials(x, time + dt, count, parity).0 - potentials(x, time - dt, count, parity).0)
        / (2.0 * dt);
    let center = potentials(x, time, count, parity).2;
    let mut u = [0.0; 3];
    for i in 0..3 {
        let mut plus = x;
        plus[i] += h;
        let mut minus = x;
        minus[i] -= h;
        let fp = potentials(plus, time, count, parity);
        let fm = potentials(minus, time, count, parity);
        p -= (fp.1[i] - fm.1[i]) / (8.0 * PI * h);
        u[i] = (fp.0 - fm.0) / (2.0 * h);
        for j in 0..3 {
            let hessian = if i == j {
                (fp.2[j] - 2.0 * center[j] + fm.2[j]) / (h * h)
            } else {
                let shifted = |di: f64, dj: f64| {
                    let mut point = x;
                    point[i] += di * h;
                    point[j] += dj * h;
                    potentials(point, time, count, parity).2[j]
                };
                (shifted(1.0, 1.0) - shifted(1.0, -1.0) - shifted(-1.0, 1.0) + shifted(-1.0, -1.0))
                    / (4.0 * h * h)
            };
            u[i] += hessian / (4.0 * PI * 1.4);
        }
    }
    (p, u)
}

fn near(actual: f64, expected: f64, relative: f64, floor: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(floor),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn field_error(actual: MovingSourceArrayField, expected: (f64, [f64; 3])) -> f64 {
    let mut error = (actual.pressure_pa - expected.0).abs();
    for i in 0..3 {
        error = error.max((actual.particle_velocity_m_s[i] - expected.1[i]).abs());
    }
    error
}

#[test]
fn independent_moving_sources_and_images_match_refined_green_potentials() {
    let sources: Vec<_> = (0..3).map(|i| path(|t| motion(i, t))).collect();
    let microphone = path(receiver);
    let windows = [[0.1, 0.6], [0.15, 0.55], [0.2, 0.5]];
    let direct: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(i, source)| {
            PreparedMovingReceiverSource::new(
                source,
                &microphone,
                config(windows[i]),
                pairs(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let planes = walls(2, 1);
    let images =
        MovingReflectionObserver::new(&microphone, &planes, 4, bounds(), &mut || false).unwrap();
    let reflected: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(i, source)| {
            PreparedReflectedMovingReceiverSource::new(
                source,
                &images,
                config(windows[i]),
                bounds(),
                pairs(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let free = PreparedMovingSourceArray::new(&direct, budget(), &mut || false).unwrap();
    let corner =
        PreparedMovingSourceArray::new_reflected(&reflected, budget(), &mut || false).unwrap();
    for (array, count, parity) in [(&free, 0, 0), (&corner, 2, 1)] {
        assert_eq!(array.source_count(), 3);
        assert_eq!(array.source_path_count(), 3 * (1 << count));
        assert_eq!(array.maximum_retarded_iterations(), 3 * (1 << count) * 96);
        assert_eq!(array.receiver_time_interval_s(), [0.2, 0.5]);
        assert!(core::ptr::eq(array.receiver(), &microphone));
        for time in [0.25, 0.45] {
            let x = receiver(time).position_m;
            let mut events = Vec::new();
            let actual = array
                .force_field_at(
                    time,
                    &mut |i, tau| Ok(drive(i, tau)),
                    &mut |i, image_index, observer, tau| {
                        events.push((i, image_index, tau));
                        for axis in 0..3 {
                            near(
                                observer.observer_position_m()[axis],
                                image(x, image_index)[axis],
                                1e-13,
                                1.0,
                            );
                        }
                        Ok(history(i, observer, tau))
                    },
                    &mut || false,
                )
                .unwrap();
            assert_eq!(events.len(), array.source_path_count());
            for (index, event) in events.iter().enumerate() {
                assert_eq!((event.0, event.1), (index >> count, index % (1 << count)));
            }
            assert_ne!(events[0].2, events[events.len() - 1].2);
            assert_eq!(actual.receiver_time_s, time);
            assert_eq!(actual.source_count, 3);
            assert_eq!(actual.source_path_count, array.source_path_count());
            assert!(actual.maximum_time_residual_s <= 1e-12);
            let coarse = reference(x, time, count, parity, 0.03);
            let fine = reference(x, time, count, parity, 0.015);
            assert!(
                field_error(actual, fine) < 0.4 * field_error(actual, coarse) + 1e-11,
                "independent Green-potential derivatives must converge toward the full native field"
            );
            let expected_p = (4.0 * fine.0 - coarse.0) / 3.0;
            near(actual.pressure_pa, expected_p, 2e-6, 1e-7);
            for i in 0..3 {
                let expected_u = (4.0 * fine.1[i] - coarse.1[i]) / 3.0;
                near(actual.particle_velocity_m_s[i], expected_u, 2e-6, 1e-7);
                near(
                    actual.instantaneous_intensity_w_m2[i],
                    expected_p * expected_u,
                    4e-6,
                    1e-10,
                );
            }
        }
    }
}

fn weak(t: f64) -> CompactSourceStrength {
    let s = t + 1.0;
    let force = 2.0_f64.powi(-60);
    CompactSourceStrength {
        force_on_fluid_n: [force * s * s, 0.0, 0.0],
        force_rate_on_fluid_n_s: [2.0 * force * s, 0.0, 0.0],
        ..CompactSourceStrength::default()
    }
}
fn cancellation_drive(index: usize, t: f64) -> CompactSourceStrength {
    let mut strength = if index == 0 {
        weak(t)
    } else {
        CompactSourceStrength::default()
    };
    let sign = if index == 0 { 1.0 } else { -1.0 };
    strength.volume_flow_m3_s = sign * (t + 1.0).powi(2);
    strength.volume_flow_rate_m3_s2 = sign * 2.0 * (t + 1.0);
    strength
}
fn weak_history(microphone: &MovingCompactObserver, tau: f64) -> [f64; 3] {
    microphone
        .force_history_integrand(
            [1.0, 0.0, 0.0],
            [2.0_f64.powi(-60) * (tau + 1.0).powi(3) / 3.0, 0.0, 0.0],
        )
        .unwrap()
}

#[test]
fn force_signal_survives_cross_source_cancellation_before_any_f64_rounding() {
    let source = path(|_| fixed([1.0, 0.0, 0.0]));
    let microphone = path(|_| fixed([5.0, 0.0, 0.0]));
    let direct: Vec<_> = (0..2)
        .map(|_| {
            PreparedMovingReceiverSource::new(
                &source,
                &microphone,
                config([0.0, 0.5]),
                pairs(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let array = PreparedMovingSourceArray::new(&direct, budget(), &mut || false).unwrap();
    let time = 0.3;
    let original = direct[0]
        .prepare_receiver(time, &mut || false)
        .unwrap()
        .pressure;
    let tau = original.source_time_s();
    let native = observer([5.0, 0.0, 0.0])
        .prepare_force_emission(tau, fixed([1.0, 0.0, 0.0]))
        .unwrap();
    let history = weak_history(&observer([5.0, 0.0, 0.0]), tau);
    let expected = native.evaluate(weak(tau), history).unwrap();
    let rounded_first = native
        .evaluate(cancellation_drive(0, tau), history)
        .unwrap();
    let rounded_second = native
        .evaluate(cancellation_drive(1, tau), [0.0; 3])
        .unwrap();
    assert_eq!(
        rounded_first.sample.pressure_pa + rounded_second.sample.pressure_pa,
        0.0
    );
    assert_eq!(
        rounded_first.particle_velocity_m_s[0] + rounded_second.particle_velocity_m_s[0],
        0.0
    );
    let actual = array
        .force_field_at(
            time,
            &mut |i, t| Ok(cancellation_drive(i, t)),
            &mut |i, _, receiver, t| {
                Ok(if i == 0 {
                    weak_history(receiver, t)
                } else {
                    [0.0; 3]
                })
            },
            &mut || false,
        )
        .unwrap();
    near(
        actual.pressure_pa,
        expected.sample.pressure_pa,
        2e-11,
        1e-100,
    );
    near(
        actual.particle_velocity_m_s[0],
        expected.particle_velocity_m_s[0],
        2e-11,
        1e-100,
    );
    near(
        actual.instantaneous_intensity_w_m2[0],
        expected.instantaneous_intensity_w_m2[0],
        4e-11,
        1e-100,
    );
    for parity in 0..2 {
        let planes = walls(1, parity);
        let images =
            MovingReflectionObserver::new(&microphone, &planes, 2, bounds(), &mut || false)
                .unwrap();
        let reflected: Vec<_> = (0..2)
            .map(|_| {
                PreparedReflectedMovingReceiverSource::new(
                    &source,
                    &images,
                    config([0.0, 0.5]),
                    bounds(),
                    pairs(),
                    &mut || false,
                )
                .unwrap()
            })
            .collect();
        let array =
            PreparedMovingSourceArray::new_reflected(&reflected, budget(), &mut || false).unwrap();
        let expected = reflected[0]
            .force_field_at(
                time,
                &mut |t| Ok(weak(t)),
                &mut |_, receiver, t| Ok(weak_history(receiver, t)),
                &mut || false,
            )
            .unwrap();
        let first = reflected[0]
            .force_field_at(
                time,
                &mut |t| Ok(cancellation_drive(0, t)),
                &mut |_, receiver, t| Ok(weak_history(receiver, t)),
                &mut || false,
            )
            .unwrap();
        let second = reflected[1]
            .force_field_at(
                time,
                &mut |t| Ok(cancellation_drive(1, t)),
                &mut |_, _, _| Ok([0.0; 3]),
                &mut || false,
            )
            .unwrap();
        assert_eq!(first.pressure_pa + second.pressure_pa, 0.0);
        let actual = array
            .force_field_at(
                time,
                &mut |i, t| Ok(cancellation_drive(i, t)),
                &mut |i, _, receiver, t| {
                    Ok(if i == 0 {
                        weak_history(receiver, t)
                    } else {
                        [0.0; 3]
                    })
                },
                &mut || false,
            )
            .unwrap();
        near(actual.pressure_pa, expected.pressure_pa, 2e-11, 1e-100);
        near(
            actual.particle_velocity_m_s[0],
            expected.particle_velocity_m_s[0],
            2e-11,
            1e-100,
        );
        near(
            actual.instantaneous_intensity_w_m2[0],
            expected.instantaneous_intensity_w_m2[0],
            4e-11,
            1e-100,
        );
    }
}

#[test]
fn admission_binds_original_receiver_boundary_and_complete_work_window() {
    let source = path(|t| motion(0, t));
    let microphone = path(receiver);
    let another = path(receiver);
    let make = |receiver, window| {
        PreparedMovingReceiverSource::new(&source, receiver, config(window), pairs(), &mut || false)
            .unwrap()
    };
    assert!(PreparedMovingSourceArray::new(&[], budget(), &mut || false).is_err());
    let mismatched = [make(&microphone, [0.1, 0.5]), make(&another, [0.1, 0.5])];
    assert!(PreparedMovingSourceArray::new(&mismatched, budget(), &mut || false).is_err());
    let uncovered = [make(&microphone, [0.1, 0.2]), make(&microphone, [0.3, 0.5])];
    assert!(PreparedMovingSourceArray::new(&uncovered, budget(), &mut || false).is_err());
    let admitted = [make(&microphone, [0.1, 0.4]), make(&microphone, [0.2, 0.5])];
    let array = PreparedMovingSourceArray::new(&admitted, budget(), &mut || false).unwrap();
    assert_eq!(array.receiver_time_interval_s(), [0.2, 0.4]);
    let exact = MovingSourceArrayBudget {
        maximum_sources: 2,
        maximum_source_paths: 2,
        maximum_retarded_iterations: array.maximum_retarded_iterations(),
    };
    assert!(PreparedMovingSourceArray::new(&admitted, exact, &mut || false).is_ok());
    for limited in [
        MovingSourceArrayBudget {
            maximum_sources: 1,
            ..exact
        },
        MovingSourceArrayBudget {
            maximum_source_paths: 1,
            ..exact
        },
        MovingSourceArrayBudget {
            maximum_retarded_iterations: exact.maximum_retarded_iterations - 1,
            ..exact
        },
    ] {
        assert!(PreparedMovingSourceArray::new(&admitted, limited, &mut || false).is_err());
    }
    let rigid =
        MovingReflectionObserver::new(&microphone, &walls(1, 0), 2, bounds(), &mut || false)
            .unwrap();
    let release =
        MovingReflectionObserver::new(&microphone, &walls(1, 1), 2, bounds(), &mut || false)
            .unwrap();
    let mixed = [
        PreparedReflectedMovingReceiverSource::new(
            &source,
            &rigid,
            config([0.1, 0.5]),
            bounds(),
            pairs(),
            &mut || false,
        )
        .unwrap(),
        PreparedReflectedMovingReceiverSource::new(
            &source,
            &release,
            config([0.1, 0.5]),
            bounds(),
            pairs(),
            &mut || false,
        )
        .unwrap(),
    ];
    assert!(PreparedMovingSourceArray::new_reflected(&mixed, budget(), &mut || false).is_err());
    assert!(matches!(
        PreparedMovingSourceArray::new(&admitted, exact, &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn last_source_failure_or_cancellation_never_publishes_a_partial_field() {
    let sources: Vec<_> = (0..2).map(|i| path(|t| motion(i, t))).collect();
    let microphone = path(receiver);
    let original: Vec<_> = sources
        .iter()
        .map(|source| {
            PreparedMovingReceiverSource::new(
                source,
                &microphone,
                config([0.2, 0.5]),
                pairs(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let array = PreparedMovingSourceArray::new(&original, budget(), &mut || false).unwrap();
    let mut callbacks = 0;
    assert!(
        array
            .force_field_at(
                0.1,
                &mut |i, t| {
                    callbacks += 1;
                    Ok(drive(i, t))
                },
                &mut |_, _, _, _| Ok([0.0; 3]),
                &mut || false
            )
            .is_err()
    );
    assert_eq!(callbacks, 0);
    let missing = AeroacError::InvalidParameter {
        what: "second source history is unavailable",
    };
    assert_eq!(
        array
            .force_field_at(
                0.3,
                &mut |i, t| Ok(drive(i, t)),
                &mut |i, _, receiver, t| if i == 0 {
                    Ok(history(i, receiver, t))
                } else {
                    Err(missing.clone())
                },
                &mut || false
            )
            .unwrap_err(),
        missing
    );
    assert!(matches!(
        array.force_field_at(
            0.3,
            &mut |i, t| Ok(drive(i, t)),
            &mut |i, _, receiver, t| Ok(if i == 0 {
                history(i, receiver, t)
            } else {
                [f64::NAN, 0.0, 0.0]
            }),
            &mut || false
        ),
        Err(AeroacError::NonFinite { .. })
    ));
    let mut polls = 0;
    let expected = array
        .force_field_at(
            0.3,
            &mut |i, t| Ok(drive(i, t)),
            &mut |i, _, receiver, t| Ok(history(i, receiver, t)),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        array.force_field_at(
            0.3,
            &mut |i, t| Ok(drive(i, t)),
            &mut |i, _, receiver, t| Ok(history(i, receiver, t)),
            &mut || {
                observed += 1;
                observed == polls
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let replay = array
        .force_field_at(
            0.3,
            &mut |i, t| Ok(drive(i, t)),
            &mut |i, _, receiver, t| Ok(history(i, receiver, t)),
            &mut || false,
        )
        .unwrap();
    assert_eq!(replay, expected);
}
