//! G1/G3 moving microphones through original, independently retarded images.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource, PreparedReflectedSourceTrajectory,
    PreparedSourceTrajectory, ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
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
fn path_budget() -> TrajectoryBudget {
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
fn config(window: [f64; 2]) -> MovingReceiverConfig {
    MovingReceiverConfig {
        receiver_time_interval_s: window,
        minimum_distance_m: 0.1,
        maximum_time_residual_s: 1e-13,
    }
}
fn knots(times: &[f64], motion: impl Fn(f64) -> SourceKinematics) -> Vec<TimedSourceKinematics> {
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
        &knots(times, motion),
        0.8,
        path_budget(),
        &mut || false,
    )
    .unwrap()
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
            3.0 + 6.0 * t + 2.0 * t * t,
            2.0 + 3.0 * t + t * t * t,
            1.0 + t * t,
        ],
        velocity_m_s: [6.0 + 4.0 * t, 3.0 + 3.0 * t * t, 2.0 * t],
        acceleration_m_s2: [4.0, 6.0 * t, 2.0],
    }
}
fn receiver(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            12.0 - 20.0 * t + 4.0 * t * t,
            8.0 + 5.0 * t - 0.5 * t * t,
            5.0 + t * t,
        ],
        velocity_m_s: [-20.0 + 8.0 * t, 5.0 - t, 2.0 * t],
        acceleration_m_s2: [8.0, -1.0, 2.0],
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
fn walls(parity: usize, count: usize) -> Vec<IdealAcousticPlane> {
    (0..count)
        .map(|axis| {
            let mut normal = [0.0; 3];
            normal[axis] = 1.0;
            IdealAcousticPlane::new(
                [0.0; 3],
                normal,
                if parity & (1 << axis) == 0 {
                    IdealPlaneBoundary::Rigid
                } else {
                    IdealPlaneBoundary::PressureRelease
                },
            )
            .unwrap()
        })
        .collect()
}
fn image(vector: [f64; 3], mask: usize) -> [f64; 3] {
    core::array::from_fn(|i| {
        if mask & (1 << i) == 0 {
            vector[i]
        } else {
            -vector[i]
        }
    })
}
fn distance(x: [f64; 3], y: [f64; 3]) -> f64 {
    (0..3).map(|i| (x[i] - y[i]).powi(2)).sum::<f64>().sqrt()
}
fn near(actual: f64, expected: f64, tolerance: f64, floor: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(floor),
        "{actual:e} != {expected:e}"
    );
}
fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
}
// Independent root of the explicitly mirrored analytic SOURCE polynomial.
fn emission(x: [f64; 3], t: f64, mask: usize) -> f64 {
    let mut lo: f64 = -0.3;
    let mut hi: f64 = 0.4;
    assert!(lo + distance(x, image(source(lo).position_m, mask)) / 300.0 < t);
    assert!(hi + distance(x, image(source(hi).position_m, mask)) / 300.0 > t);
    for _ in 0..90 {
        let tau = lo.midpoint(hi);
        if tau + distance(x, image(source(tau).position_m, mask)) / 300.0 < t {
            lo = tau;
        } else {
            hi = tau;
        }
    }
    lo.midpoint(hi)
}
fn potential(x: [f64; 3], t: f64, mask: usize) -> [f64; 4] {
    let tau = emission(x, t, mask);
    let state = source(tau);
    let y = image(state.position_m, mask);
    let v = image(state.velocity_m_s, mask);
    let r = distance(x, y);
    let radial: f64 = (0..3).map(|i| (x[i] - y[i]) * v[i] / r).sum();
    let denominator = 4.0 * PI * r * (1.0 - radial / 300.0);
    let strengths = strength(tau);
    let force = image(strengths.force_on_fluid_n, mask);
    [
        strengths.volume_flow_m3_s / denominator,
        force[0] / denominator,
        force[1] / denominator,
        force[2] / denominator,
    ]
}

#[test]
fn accelerating_images_match_independent_green_derivatives_and_both_moving_clocks() {
    let source = path(&[-0.3, -0.04, 0.11, 0.4], source);
    let receiver_path = path(&[0.0, 0.07, 0.2], receiver);
    let observer =
        MovingReflectionObserver::new(&receiver_path, &walls(5, 3), 8, bounds(), &mut || false)
            .unwrap();
    let reflected = PreparedReflectedMovingReceiverSource::new(
        &source,
        &observer,
        config([0.0, 0.2]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    let direct = PreparedMovingReceiverSource::new(
        &source,
        &receiver_path,
        config([0.0, 0.2]),
        pairs(),
        &mut || false,
    )
    .unwrap();
    assert!(core::ptr::eq(reflected.source(), &source));
    assert!(core::ptr::eq(reflected.receiver(), &receiver_path));
    for time in [0.01, 0.08, 0.19] {
        let x = receiver(time).position_m;
        let mut expected = 0.0;
        for mask in 0..8 {
            let sign = if (mask & 5_usize).count_ones() % 2 == 0 {
                1.0
            } else {
                -1.0
            };
            // Eulerian time differentiation keeps x fixed at this event.
            // Differentiating along the microphone would invent pressure gain.
            let q = 1.2 * derivative(|dt| potential(x, time + dt, mask)[0], 2e-5);
            let f: f64 = -(0..3)
                .map(|axis| {
                    derivative(
                        |dx| {
                            let mut shifted = x;
                            shifted[axis] += dx;
                            potential(shifted, time, mask)[axis + 1]
                        },
                        2e-3,
                    )
                })
                .sum::<f64>();
            let sample = reflected.prepare_path(mask, time, &mut || false).unwrap();
            let actual = sample
                .pressure
                .evaluate(strength(sample.pressure.source_time_s()))
                .unwrap();
            near(actual.monopole_pressure_pa, sign * q, 5e-7, 1e-9);
            near(actual.dipole_pressure_pa, sign * f, 5e-7, 1e-9);
            let history = strength(sample.pressure.source_time_s());
            let field = sample
                .pressure
                .evaluate_monopole_field(history.volume_flow_m3_s, history.volume_flow_rate_m3_s2)
                .unwrap();
            assert_eq!(field.sample.pressure_pa, actual.monopole_pressure_pa);
            for axis in 0..3 {
                let velocity = -sign
                    * derivative(
                        |dx| {
                            let mut shifted = x;
                            shifted[axis] += dx;
                            potential(shifted, time, mask)[0]
                        },
                        2e-3,
                    );
                near(field.particle_velocity_m_s[axis], velocity, 5e-7, 1e-10);
                near(
                    field.instantaneous_intensity_w_m2[axis],
                    sign * q * velocity,
                    1e-6,
                    1e-12,
                );
            }
            let endpoint = reflected
                .prepare_path_emission(mask, sample.pressure.source_time_s(), time)
                .unwrap()
                .evaluate_monopole_field(history.volume_flow_m3_s, history.volume_flow_rate_m3_s2)
                .unwrap();
            assert_eq!(field.particle_velocity_m_s, endpoint.particle_velocity_m_s);
            near(actual.source_time_s, emission(x, time, mask), 3e-11, 0.1);
            let rate = derivative(
                |dt| emission(receiver(time + dt).position_m, time + dt, mask),
                2e-5,
            );
            near(sample.observed_source_time_rate, rate, 3e-9, 1.0);
            assert!(
                sample.observed_source_time_rate
                    <= reflected.path_maximum_source_time_rate(mask).unwrap()
            );
            assert_eq!(
                sample.receiver_kinematics,
                receiver_path.kinematics_at(time).unwrap()
            );
            if mask == 0 {
                let original = direct.prepare_receiver(time, &mut || false).unwrap();
                assert_eq!(
                    actual,
                    original
                        .pressure
                        .evaluate(strength(original.pressure.source_time_s()))
                        .unwrap()
                );
                assert_eq!(
                    sample.observed_source_time_rate,
                    original.observed_source_time_rate
                );
            }
            expected += sign * (q + f);
        }
        near(
            reflected
                .pressure_at(time, &mut |tau| Ok(strength(tau)), &mut || false)
                .unwrap(),
            expected,
            2e-6,
            1e-8,
        );
    }
}

#[test]
fn reflected_doppler_uses_the_image_receiver_velocity() {
    let motion = |origin, t| SourceKinematics {
        position_m: [origin - 90.0 * t, 1.0, 1.0],
        velocity_m_s: [-90.0, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    };
    let source = path(&[-0.3, 0.03], |t| motion(3.0, t));
    let receiver = path(&[0.0, 0.02], |t| motion(20.0, t));
    let observer =
        MovingReflectionObserver::new(&receiver, &walls(0, 1), 2, bounds(), &mut || false).unwrap();
    let reflected = PreparedReflectedMovingReceiverSource::new(
        &source,
        &observer,
        config([0.0, 0.02]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    assert!(reflected.path_maximum_source_time_rate(0).unwrap() < 1.000_001);
    assert!(reflected.path_maximum_source_time_rate(1).unwrap() >= 1.3 / 0.7);
    for t in [0.0, 0.01, 0.02] {
        let direct = reflected.prepare_path(0, t, &mut || false).unwrap();
        let image = reflected.prepare_path(1, t, &mut || false).unwrap();
        near(
            direct.pressure.source_time_s(),
            t - 17.0 / 390.0,
            1e-11,
            0.1,
        );
        near(
            image.pressure.source_time_s(),
            (1.3 * t - 23.0 / 300.0) / 0.7,
            1e-11,
            0.1,
        );
        near(direct.observed_source_time_rate, 1.0, 2e-13, 1.0);
        near(image.observed_source_time_rate, 1.3 / 0.7, 2e-13, 1.0);
    }
}

#[test]
fn moving_microphone_retains_weak_force_after_image_path_cancellation() {
    let source = path(&[-0.2, 0.3], |_| fixed([1.0, 0.0, 0.0]));
    let receiver = path(&[0.0, 0.1, 0.2], |time| SourceKinematics {
        position_m: [3.0 + 10.0 * (time - 0.1), 0.0, 0.0],
        velocity_m_s: [10.0, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    });
    let weak_force = 2.0_f64.powi(-60);
    for parity in 0..2 {
        let observer =
            MovingReflectionObserver::new(&receiver, &walls(parity, 1), 2, bounds(), &mut || false)
                .unwrap();
        let reflected = PreparedReflectedMovingReceiverSource::new(
            &source,
            &observer,
            config([0.05, 0.15]),
            bounds(),
            pairs(),
            &mut || false,
        )
        .unwrap();
        let direct = reflected.prepare_path(0, 0.1, &mut || false).unwrap();
        let image = reflected.prepare_path(1, 0.1, &mut || false).unwrap();
        let sign = if parity == 0 { 1.0 } else { -1.0 };
        let early_time = image.pressure.source_time_s();
        let duration = direct.pressure.source_time_s() - early_time;
        assert!(duration > 0.0);
        assert_ne!(direct.observed_source_time_rate, 1.0);

        // The receiver crosses x=3 at this event. A physical quadratic
        // volume history cancels the strong monopoles at distances 2 and 4;
        // a weak steady normal force retains its two near-field pressures.
        let history = |time| {
            let u = (time - early_time) / duration;
            let early_rate = -2.0 * sign;
            let slope = 1.0 - early_rate;
            CompactSourceStrength {
                volume_flow_m3_s: duration * u * (early_rate + 0.5 * slope * u),
                volume_flow_rate_m3_s2: early_rate + slope * u,
                force_on_fluid_n: [weak_force, 0.0, 0.0],
                ..CompactSourceStrength::default()
            }
        };
        let direct_sample = direct
            .pressure
            .evaluate(history(direct.pressure.source_time_s()))
            .unwrap();
        let image_sample = image
            .pressure
            .evaluate(history(image.pressure.source_time_s()))
            .unwrap();
        assert_eq!(direct_sample.pressure_pa + image_sample.pressure_pa, 0.0);
        let expected = weak_force * (0.25 - sign * 0.0625) / (4.0 * PI);
        near(
            reflected
                .pressure_at(0.1, &mut |time| Ok(history(time)), &mut || false)
                .unwrap(),
            expected,
            1e-13,
            0.0,
        );
    }
}

#[test]
fn a_microphone_moving_on_a_release_wall_retains_exact_pressure_cancellation() {
    let source = path(&[-0.3, -0.04, 0.4], source);
    let receiver = path(&[0.0, 0.2], |t| SourceKinematics {
        position_m: [0.0, 8.0 + 5.0 * t + t * t, 5.0 - 2.0 * t],
        velocity_m_s: [0.0, 5.0 + 2.0 * t, -2.0],
        acceleration_m_s2: [0.0, 2.0, 0.0],
    });
    let observer =
        MovingReflectionObserver::new(&receiver, &walls(1, 1), 2, bounds(), &mut || false).unwrap();
    let reflected = PreparedReflectedMovingReceiverSource::new(
        &source,
        &observer,
        config([0.0, 0.2]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    for t in [0.0, 0.07, 0.2] {
        assert_eq!(
            reflected
                .pressure_at(t, &mut |tau| Ok(strength(tau)), &mut || false)
                .unwrap(),
            0.0
        );
    }
}

#[test]
fn a_stationary_microphone_recovers_the_original_fixed_reflection_owner() {
    let position = [12.0, 8.0, 5.0];
    let source_knots = knots(&[-0.3, -0.04, 0.11, 0.4], source);
    let source =
        PreparedKinematicTrajectory::new(medium(), &source_knots, 0.8, path_budget(), &mut || {
            false
        })
        .unwrap();
    let receiver = path(&[0.0, 0.2], |_| fixed(position));
    let observer =
        MovingReflectionObserver::new(&receiver, &walls(2, 2), 4, bounds(), &mut || false).unwrap();
    let reflected = PreparedReflectedMovingReceiverSource::new(
        &source,
        &observer,
        config([0.0, 0.2]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    let old_observer = MovingCompactObserver::new(
        medium(),
        position,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.1,
        },
    )
    .unwrap();
    let old_source =
        PreparedSourceTrajectory::new(old_observer, &source_knots, path_budget(), &mut || false)
            .unwrap();
    let old_observer = IdealReflectionObserver::new(old_observer, &walls(2, 2), 4).unwrap();
    let old =
        PreparedReflectedSourceTrajectory::new(&old_source, &old_observer, bounds(), &mut || false)
            .unwrap();
    for t in [0.01, 0.08, 0.19] {
        near(
            reflected
                .pressure_at(t, &mut |tau| Ok(strength(tau)), &mut || false)
                .unwrap(),
            old.pressure_at(t, 1e-13, &mut |tau| Ok(strength(tau)), &mut || false)
                .unwrap(),
            2e-9,
            1e-8,
        );
    }
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1], v[2]]
}
fn rotated(mut state: SourceKinematics) -> SourceKinematics {
    let shift = [4.0, -3.0, 2.0];
    let p = rotate(state.position_m);
    state.position_m = core::array::from_fn(|i| p[i] + shift[i]);
    state.velocity_m_s = rotate(state.velocity_m_s);
    state.acceleration_m_s2 = rotate(state.acceleration_m_s2);
    state
}
#[test]
fn moving_images_keep_source_vector_orientation_under_a_common_rigid_frame_change() {
    let a = path(&[-0.3, -0.04, 0.11, 0.4], source);
    let b = path(&[0.0, 0.07, 0.2], receiver);
    let ar = path(&[-0.3, -0.04, 0.11, 0.4], |t| rotated(source(t)));
    let br = path(&[0.0, 0.07, 0.2], |t| rotated(receiver(t)));
    let planes = walls(5, 3);
    let rotated_planes: Vec<_> = planes
        .iter()
        .map(|p| {
            IdealAcousticPlane::new([4.0, -3.0, 2.0], rotate(p.inward_normal()), p.boundary())
                .unwrap()
        })
        .collect();
    let observer = MovingReflectionObserver::new(&b, &planes, 8, bounds(), &mut || false).unwrap();
    let observer_r =
        MovingReflectionObserver::new(&br, &rotated_planes, 8, bounds(), &mut || false).unwrap();
    let p = PreparedReflectedMovingReceiverSource::new(
        &a,
        &observer,
        config([0.0, 0.2]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    let pr = PreparedReflectedMovingReceiverSource::new(
        &ar,
        &observer_r,
        config([0.0, 0.2]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    for t in [0.01, 0.08, 0.19] {
        let expected = p
            .pressure_at(t, &mut |tau| Ok(strength(tau)), &mut || false)
            .unwrap();
        let actual = pr
            .pressure_at(
                t,
                &mut |tau| {
                    let mut s = strength(tau);
                    s.force_on_fluid_n = rotate(s.force_on_fluid_n);
                    s.force_rate_on_fluid_n_s = rotate(s.force_rate_on_fluid_n_s);
                    Ok(s)
                },
                &mut || false,
            )
            .unwrap();
        near(actual, expected, 2e-9, 1e-8);
    }
}

#[test]
fn whole_path_confinement_complete_history_aggregate_work_and_cancellation_are_required() {
    let planes = walls(0, 1);
    let crossing = path(&[0.0, 1.0], |t| SourceKinematics {
        position_m: [-0.1 + (t - 0.5).powi(2), 8.0, 5.0],
        velocity_m_s: [2.0 * (t - 0.5), 0.0, 0.0],
        acceleration_m_s2: [2.0, 0.0, 0.0],
    });
    assert!(MovingReflectionObserver::new(&crossing, &planes, 2, bounds(), &mut || false).is_err());
    let receiver = path(&[0.0, 0.2], |_| fixed([12.0, 8.0, 5.0]));
    assert!(MovingReflectionObserver::new(&receiver, &planes, 1, bounds(), &mut || false).is_err());
    assert!(
        MovingReflectionObserver::new(&receiver, &planes.repeat(2), 4, bounds(), &mut || false)
            .is_err()
    );
    assert!(matches!(
        MovingReflectionObserver::new(&receiver, &planes, 2, bounds(), &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let observer =
        MovingReflectionObserver::new(&receiver, &planes, 2, bounds(), &mut || false).unwrap();
    let source = path(&[-0.3, 0.4], |_| fixed([3.0, 2.0, 1.0]));
    let budget = MovingReceiverBudget {
        maximum_segment_pairs: 2,
        maximum_bound_nodes: 2,
        ..pairs()
    };
    let admitted = PreparedReflectedMovingReceiverSource::new(
        &source,
        &observer,
        config([0.04, 0.1]),
        bounds(),
        budget,
        &mut || false,
    )
    .unwrap();
    assert!(
        PreparedReflectedMovingReceiverSource::new(
            &source,
            &observer,
            config([0.04, 0.1]),
            bounds(),
            MovingReceiverBudget {
                maximum_segment_pairs: 1,
                ..budget
            },
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedReflectedMovingReceiverSource::new(
            &source,
            &observer,
            config([0.04, 0.1]),
            bounds(),
            MovingReceiverBudget {
                maximum_bound_nodes: 1,
                ..budget
            },
            &mut || false
        )
        .is_err()
    );
    let late = path(&[-0.005, 0.2], |_| fixed([3.0, 2.0, 1.0]));
    assert!(
        PreparedMovingReceiverSource::new(
            &late,
            &receiver,
            config([0.04, 0.1]),
            pairs(),
            &mut || false
        )
        .is_ok()
    );
    assert!(
        PreparedReflectedMovingReceiverSource::new(
            &late,
            &observer,
            config([0.04, 0.1]),
            bounds(),
            pairs(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedReflectedMovingReceiverSource::new(
            &crossing,
            &observer,
            config([0.04, 0.1]),
            bounds(),
            pairs(),
            &mut || false
        )
        .is_err()
    );
    let before = admitted
        .pressure_at(0.07, &mut |tau| Ok(strength(tau)), &mut || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        admitted.pressure_at(0.07, &mut |tau| Ok(strength(tau)), &mut || {
            polls += 1;
            polls == 6
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        before,
        admitted
            .pressure_at(0.07, &mut |tau| Ok(strength(tau)), &mut || false)
            .unwrap()
    );
    assert!(admitted.prepare_path(2, 0.07, &mut || false).is_err());
    assert!(admitted.prepare_path(0, 0.03, &mut || false).is_err());
}

fn wedge_planes(order: usize, parity: usize) -> [IdealAcousticPlane; 2] {
    let (s, c) = (PI / order as f64).sin_cos();
    let normals = [
        [0.0, 1.0, 0.0],
        if order == 2 {
            [1.0, 0.0, 0.0]
        } else {
            [s, -c, 0.0]
        },
    ];
    core::array::from_fn(|i| {
        IdealAcousticPlane::new(
            [0.0; 3],
            normals[i],
            if parity & (1 << i) == 0 {
                IdealPlaneBoundary::Rigid
            } else {
                IdealPlaneBoundary::PressureRelease
            },
        )
        .unwrap()
    })
}
fn wedge_motion(mut value: SourceKinematics) -> SourceKinematics {
    value.position_m[1] *= 0.2;
    value.velocity_m_s[1] *= 0.2;
    value.acceleration_m_s2[1] *= 0.2;
    value
}
// Independent polar SOURCE image, rather than the production affine RECEIVER.
fn wedge_vector(v: [f64; 3], index: usize, order: usize) -> [f64; 3] {
    let (s, c) = (2.0 * PI * (index / 2) as f64 / order as f64).sin_cos();
    let y = if index % 2 == 0 { v[1] } else { -v[1] };
    [c * v[0] - s * y, s * v[0] + c * y, v[2]]
}
fn wedge_image(t: f64, index: usize, order: usize) -> SourceKinematics {
    let state = wedge_motion(source(t));
    SourceKinematics {
        position_m: wedge_vector(state.position_m, index, order),
        velocity_m_s: wedge_vector(state.velocity_m_s, index, order),
        acceleration_m_s2: wedge_vector(state.acceleration_m_s2, index, order),
    }
}

fn near_wedge_monopole(
    actual: &fs_aeroac::curle3d::moving::MovingMonopoleField,
    expected: &fs_aeroac::curle3d::moving::MovingMonopoleField,
) {
    near(actual.sample.pressure_pa, expected.sample.pressure_pa, 5e-10, 1e-8);
    for axis in 0..3 {
        near(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            5e-10,
            1e-8,
        );
        near(
            actual.instantaneous_intensity_w_m2[axis],
            expected.instantaneous_intensity_w_m2[axis],
            1e-9,
            1e-10,
        );
    }
}

#[test]
fn dihedral_moving_images_preserve_each_retarded_clock_dipole_and_coherent_sum() {
    let original = path(&[-0.3, -0.04, 0.11, 0.4], |t| wedge_motion(source(t)));
    let receiver_path = path(&[0.0, 0.07, 0.2], |t| wedge_motion(receiver(t)));
    for (order, parity) in [(3, 0), (3, 3), (4, 1), (6, 2)] {
        let observer = MovingReflectionObserver::new_wedge(
            &receiver_path,
            wedge_planes(order, parity),
            order,
            2 * order,
            bounds(),
            &mut || false,
        )
        .unwrap();
        let reflected = PreparedReflectedMovingReceiverSource::new(
            &original,
            &observer,
            config([0.0, 0.2]),
            bounds(),
            pairs(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(observer.wedge_order(), Some(order));
        assert_eq!(reflected.path_count(), 2 * order);
        assert!(core::ptr::eq(reflected.source(), &original));
        assert!(core::ptr::eq(reflected.receiver(), &receiver_path));
        let mut expected = [0.0; 3];
        for index in 0..2 * order {
            let image = path(&[-0.3, -0.04, 0.11, 0.4], |t| wedge_image(t, index, order));
            let reference = PreparedMovingReceiverSource::new(
                &image,
                &receiver_path,
                config([0.0, 0.2]),
                pairs(),
                &mut || false,
            )
            .unwrap();
            let first_sign = if parity & 1 == 0 { 1.0 } else { -1.0 };
            let step_sign = if parity == 1 || parity == 2 {
                -1.0_f64
            } else {
                1.0
            };
            let sign =
                step_sign.powi((index / 2) as i32) * if index % 2 == 0 { 1.0 } else { first_sign };
            for (sample_index, time) in [0.01, 0.08, 0.19].into_iter().enumerate() {
                let sample = reflected.prepare_path(index, time, &mut || false).unwrap();
                let oracle = reference.prepare_receiver(time, &mut || false).unwrap();
                let mut strengths = strength(oracle.pressure.source_time_s());
                strengths.force_on_fluid_n = wedge_vector(strengths.force_on_fluid_n, index, order);
                strengths.force_rate_on_fluid_n_s =
                    wedge_vector(strengths.force_rate_on_fluid_n_s, index, order);
                let oracle_pressure = oracle.pressure.evaluate(strengths).unwrap();
                let actual = sample
                    .pressure
                    .evaluate(strength(sample.pressure.source_time_s()))
                    .unwrap();
                let drive = strength(sample.pressure.source_time_s());
                let oracle_drive = strength(oracle.pressure.source_time_s());
                let field = sample
                    .pressure
                    .evaluate_monopole_field(
                        drive.volume_flow_m3_s,
                        drive.volume_flow_rate_m3_s2,
                    )
                    .unwrap();
                let oracle_field = oracle
                    .pressure
                    .evaluate_monopole_field(
                        sign * oracle_drive.volume_flow_m3_s,
                        sign * oracle_drive.volume_flow_rate_m3_s2,
                    )
                    .unwrap();
                near_wedge_monopole(&field, &oracle_field);
                let tau = sample.pressure.source_time_s();
                let emitted = reflected
                    .prepare_path_emission(index, tau, time)
                    .unwrap()
                    .evaluate_monopole_field(
                        drive.volume_flow_m3_s,
                        drive.volume_flow_rate_m3_s2,
                    )
                    .unwrap();
                let oracle_emitted = MovingCompactObserver::new(
                    medium(),
                    receiver_path.kinematics_at(time).unwrap().position_m,
                    MovingCompactLimits {
                        maximum_mach: 0.8,
                        minimum_distance_m: 0.1,
                    },
                )
                    .unwrap()
                    .prepare_emission(tau, image.kinematics_at(tau).unwrap())
                    .unwrap()
                    .evaluate_monopole_field(
                        sign * drive.volume_flow_m3_s,
                        sign * drive.volume_flow_rate_m3_s2,
                    )
                    .unwrap();
                near_wedge_monopole(&emitted, &oracle_emitted);
                if index == 0 {
                    assert_eq!(field.particle_velocity_m_s, oracle_field.particle_velocity_m_s);
                    assert_eq!(
                        field.instantaneous_intensity_w_m2,
                        oracle_field.instantaneous_intensity_w_m2
                    );
                }
                near(
                    actual.source_time_s,
                    oracle.pressure.source_time_s(),
                    2e-12,
                    0.1,
                );
                near(
                    actual.monopole_pressure_pa,
                    sign * oracle_pressure.monopole_pressure_pa,
                    3e-11,
                    1e-8,
                );
                near(
                    actual.dipole_pressure_pa,
                    sign * oracle_pressure.dipole_pressure_pa,
                    3e-11,
                    1e-8,
                );
                near(
                    sample.observed_source_time_rate,
                    oracle.observed_source_time_rate,
                    3e-12,
                    1.0,
                );
                assert!(
                    sample.observed_source_time_rate
                        <= reflected.path_maximum_source_time_rate(index).unwrap()
                );
                assert_eq!(
                    sample.receiver_kinematics,
                    receiver_path.kinematics_at(time).unwrap()
                );
                if index == 0 {
                    assert_eq!(actual, oracle_pressure);
                    assert_eq!(
                        sample.observed_source_time_rate,
                        oracle.observed_source_time_rate
                    );
                }
                expected[sample_index] += sign * oracle_pressure.pressure_pa;
            }
        }
        for (time, expected) in [0.01, 0.08, 0.19].into_iter().zip(expected) {
            near(
                reflected
                    .pressure_at(time, &mut |t| Ok(strength(t)), &mut || false)
                    .unwrap(),
                expected,
                3e-10,
                1e-8,
            );
        }
    }
}

#[test]
fn stationary_wedge_paths_recover_original_harmonic_phase_and_wall_signs() {
    use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
    use fs_math::c64::C64;
    let x = [8.0, 1.5, 3.0];
    let y = [3.0, 0.5, 1.0];
    let observer = MovingCompactObserver::new(
        medium(),
        x,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.1,
        },
    )
    .unwrap();
    let base = PreparedSourceTrajectory::new(
        observer,
        &knots(&[-0.3, 0.4], |_| fixed(y)),
        path_budget(),
        &mut || false,
    )
    .unwrap();
    let physical_source = path(&[-0.3, 0.4], |_| fixed(y));
    let microphone = path(&[0.0, 0.2], |_| fixed(x));
    let omega = 137.0;
    let phasor = CompactSourcePhasor {
        volume_flow_m3_s: C64::new(0.003, -0.001),
        force_on_fluid_n: [C64::new(1.3, -0.4), C64::new(-0.7, 0.9), C64::new(0.2, 0.8)],
    };
    let drive = |t: f64| {
        let (s, c) = (omega * t).sin_cos();
        let value = |z: C64| z.re * c + z.im * s;
        let rate = |z: C64| omega * (-z.re * s + z.im * c);
        CompactSourceStrength {
            volume_flow_m3_s: value(phasor.volume_flow_m3_s),
            volume_flow_rate_m3_s2: rate(phasor.volume_flow_m3_s),
            force_on_fluid_n: phasor.force_on_fluid_n.map(value),
            force_rate_on_fluid_n_s: phasor.force_on_fluid_n.map(rate),
        }
    };
    for (order, parity) in [(2, 1), (3, 0), (3, 3), (4, 2), (8, 0)] {
        let planes = wedge_planes(order, parity);
        let fixed_observer =
            IdealReflectionObserver::new_wedge(observer, planes, order, 2 * order).unwrap();
        let fixed_paths =
            PreparedReflectedSourceTrajectory::new(&base, &fixed_observer, bounds(), &mut || false)
                .unwrap();
        let moving_observer = MovingReflectionObserver::new_wedge(
            &microphone,
            planes,
            order,
            2 * order,
            bounds(),
            &mut || false,
        )
        .unwrap();
        let moving_paths = PreparedReflectedMovingReceiverSource::new(
            &physical_source,
            &moving_observer,
            config([0.0, 0.2]),
            bounds(),
            pairs(),
            &mut || false,
        )
        .unwrap();
        let harmonic = PreparedCompactSourceObserver::new(medium(), omega / (2.0 * PI), x, &[y], 1)
            .unwrap()
            .with_wedge(planes, order, 2 * order)
            .unwrap()
            .evaluate(&[phasor])
            .unwrap();
        assert_eq!(fixed_observer.wedge_order(), Some(order));
        assert_eq!(fixed_paths.path_count(), 2 * order);
        for index in 0..2 * order {
            let (_, sign, velocity_map) = fixed_observer.stationary_field_path(index, y).unwrap();
            let first_sign = if parity & 1 == 0 { 1.0 } else { -1.0 };
            let step_sign = if parity == 1 || parity == 2 { -1.0_f64 } else { 1.0 };
            let expected_sign = step_sign.powi((index / 2) as i32)
                * if index % 2 == 0 { 1.0 } else { first_sign };
            assert_eq!(sign, expected_sign);
            for column in 0..3 {
                let mut basis = [0.0; 3];
                basis[column] = 1.0;
                let expected = wedge_vector(basis, index, order);
                for row in 0..3 {
                    near(velocity_map[row][column], sign * expected[row], 3e-12, 1.0);
                }
            }
        }
        for time in [0.01, 0.08, 0.19] {
            for index in 0..2 * order {
                let fixed_sample = fixed_paths
                    .prepare_path(index, time, 1e-13, &mut || false)
                    .unwrap();
                let moving_sample = moving_paths
                    .prepare_path(index, time, &mut || false)
                    .unwrap();
                let fixed_drive = drive(fixed_sample.source_time_s());
                let moving_drive = drive(moving_sample.pressure.source_time_s());
                let field = fixed_sample
                    .evaluate_monopole_field(
                        fixed_drive.volume_flow_m3_s,
                        fixed_drive.volume_flow_rate_m3_s2,
                    )
                    .unwrap();
                let moving_field = moving_sample
                    .pressure
                    .evaluate_monopole_field(
                        moving_drive.volume_flow_m3_s,
                        moving_drive.volume_flow_rate_m3_s2,
                    )
                    .unwrap();
                let emission_field = fixed_paths
                    .prepare_path_emission(index, fixed_sample.source_time_s())
                    .unwrap()
                    .evaluate_monopole_field(
                        fixed_drive.volume_flow_m3_s,
                        fixed_drive.volume_flow_rate_m3_s2,
                    )
                    .unwrap();
                near_wedge_monopole(&field, &moving_field);
                near_wedge_monopole(&emission_field, &field);
            }
            let (s, c) = (omega * time).sin_cos();
            let expected = harmonic.pressure_pa.re * c + harmonic.pressure_pa.im * s;
            near(
                fixed_paths
                    .pressure_at(time, 1e-13, &mut |t| Ok(drive(t)), &mut || false)
                    .unwrap(),
                expected,
                5e-10,
                1e-8,
            );
            near(
                moving_paths
                    .pressure_at(time, &mut |t| Ok(drive(t)), &mut || false)
                    .unwrap(),
                expected,
                5e-10,
                1e-8,
            );
        }
    }
}

#[test]
fn wedge_preparation_enforces_geometry_full_group_work_and_cancellation() {
    let source = path(&[-0.3, 0.4], |_| fixed([3.0, 0.5, 1.0]));
    let receiver = path(&[0.0, 0.2], |_| fixed([8.0, 1.5, 3.0]));
    let prepare = |planes, order, cap| {
        MovingReflectionObserver::new_wedge(&receiver, planes, order, cap, bounds(), &mut || false)
    };
    assert!(prepare(wedge_planes(3, 1), 3, 6).is_err());
    assert!(prepare(wedge_planes(4, 0), 3, 8).is_err());
    assert!(prepare(wedge_planes(4, 0), 4, 7).is_err());
    assert!(prepare(wedge_planes(4, 0), 129, 258).is_err());
    let observer = prepare(wedge_planes(4, 0), 4, 8).unwrap();
    assert!(
        PreparedReflectedMovingReceiverSource::new(
            &source,
            &observer,
            config([0.0, 0.2]),
            bounds(),
            MovingReceiverBudget {
                maximum_segment_pairs: 7,
                ..pairs()
            },
            &mut || false
        )
        .is_err()
    );
    let complete = PreparedReflectedMovingReceiverSource::new(
        &source,
        &observer,
        config([0.0, 0.2]),
        bounds(),
        MovingReceiverBudget {
            maximum_segment_pairs: 8,
            ..pairs()
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(complete.path_count(), 8);
    let mut polls = 0;
    assert!(matches!(
        MovingReflectionObserver::new_wedge(
            &receiver,
            wedge_planes(8, 0),
            8,
            16,
            bounds(),
            &mut || {
                polls += 1;
                polls == 10
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        complete.pressure_at(0.1, &mut |t| Ok(strength(t)), &mut || {
            polls += 1;
            polls == 10
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        complete
            .pressure_at(0.1, &mut |t| Ok(strength(t)), &mut || false)
            .unwrap()
            .is_finite()
    );
}

#[test]
fn moving_wedge_affine_images_keep_a_translated_rotated_physical_frame() {
    let source_path = path(&[-0.3, 0.0, 0.4], |t| wedge_motion(source(t)));
    let receiver_path = path(&[0.0, 0.2], |t| wedge_motion(receiver(t)));
    let rotated_source = path(&[-0.3, 0.0, 0.4], |t| rotated(wedge_motion(source(t))));
    let rotated_receiver = path(&[0.0, 0.2], |t| rotated(wedge_motion(receiver(t))));
    let planes = wedge_planes(4, 1);
    let rotated_planes = planes.map(|plane| {
        IdealAcousticPlane::new(
            [4.0, -3.0, 2.0],
            rotate(plane.inward_normal()),
            plane.boundary(),
        )
        .unwrap()
    });
    let observer =
        MovingReflectionObserver::new_wedge(&receiver_path, planes, 4, 8, bounds(), &mut || false)
            .unwrap();
    let other = MovingReflectionObserver::new_wedge(
        &rotated_receiver,
        rotated_planes,
        4,
        8,
        bounds(),
        &mut || false,
    )
    .unwrap();
    let original = PreparedReflectedMovingReceiverSource::new(
        &source_path,
        &observer,
        config([0.0, 0.2]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    let rotated = PreparedReflectedMovingReceiverSource::new(
        &rotated_source,
        &other,
        config([0.0, 0.2]),
        bounds(),
        pairs(),
        &mut || false,
    )
    .unwrap();
    for time in [0.01, 0.08, 0.19] {
        let expected = original
            .pressure_at(time, &mut |t| Ok(strength(t)), &mut || false)
            .unwrap();
        let actual = rotated
            .pressure_at(
                time,
                &mut |t| {
                    let mut drive = strength(t);
                    drive.force_on_fluid_n = rotate(drive.force_on_fluid_n);
                    drive.force_rate_on_fluid_n_s = rotate(drive.force_rate_on_fluid_n_s);
                    Ok(drive)
                },
                &mut || false,
            )
            .unwrap();
        near(actual, expected, 3e-10, 1e-8);
    }
}
