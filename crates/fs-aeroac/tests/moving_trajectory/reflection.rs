//! G1/G3 independently delayed and coherently summed moving-source images.

use super::*;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, PreparedReflectedSourceTrajectory, ReflectedTrajectoryBudget,
};
use fs_aeroac::curle3d::{
    CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary, PreparedCompactSourceObserver,
};
use fs_math::c64::C64;

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

fn reflection_budget() -> ReflectedTrajectoryBudget {
    ReflectedTrajectoryBudget {
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
    }
}

fn motion(t: f64) -> SourceKinematics {
    let mut state = polynomial(t);
    state.position_m = state.position_m.map(|x| x + 3.0);
    state
}

fn source_path(x: [f64; 3]) -> PreparedSourceTrajectory {
    prepare(
        observer(x, 300.0, 0.8, 0.05),
        &samples(&[-0.2, -0.03, 0.11, 0.4], motion),
    )
}

fn harmonic_strength(phasor: CompactSourcePhasor, omega: f64, time: f64) -> CompactSourceStrength {
    let (s, c) = (omega * time).sin_cos();
    let value = |z: C64| z.re * c + z.im * s;
    let rate = |z: C64| omega * (-z.re * s + z.im * c);
    CompactSourceStrength {
        volume_flow_m3_s: value(phasor.volume_flow_m3_s),
        volume_flow_rate_m3_s2: rate(phasor.volume_flow_m3_s),
        force_on_fluid_n: phasor.force_on_fluid_n.map(value),
        force_rate_on_fluid_n_s: phasor.force_on_fluid_n.map(rate),
    }
}

#[test]
fn stationary_transient_images_recover_harmonic_corner_interference_and_each_delay() {
    let x = [7.0, 5.0, 4.0];
    let y = [1.0, 2.0, 3.0];
    let medium = FluidAcousticState::try_new(1.2, 300.0).unwrap();
    let base = prepare(
        observer(x, 300.0, 0.0, 0.05),
        &samples(&[-1.0, 1.0], |_| SourceKinematics {
            position_m: y,
            ..SourceKinematics::default()
        }),
    );
    let omega = 137.0;
    let phasor = CompactSourcePhasor {
        volume_flow_m3_s: C64::new(0.003, -0.001),
        force_on_fluid_n: [C64::new(1.3, -0.4), C64::new(-0.7, 0.9), C64::new(0.2, 0.8)],
    };
    for parity in 0..8 {
        let planes = walls(parity, 3);
        let images = IdealReflectionObserver::new(*base.observer(), &planes, 8).unwrap();
        let transient = PreparedReflectedSourceTrajectory::new(
            &base,
            &images,
            reflection_budget(),
            &mut || false,
        )
        .unwrap();
        let harmonic = PreparedCompactSourceObserver::new(medium, omega / (2.0 * PI), x, &[y], 1)
            .unwrap()
            .with_orthogonal_planes(&planes, 8)
            .unwrap()
            .evaluate(&[phasor])
            .unwrap();
        assert_eq!(transient.path_count(), 8);
        for time in [0.05, 0.137, 0.3] {
            let actual = transient
                .pressure_at(
                    time,
                    1e-14,
                    &mut |tau| Ok(harmonic_strength(phasor, omega, tau)),
                    &mut || false,
                )
                .unwrap();
            let (s, c) = (omega * time).sin_cos();
            near(
                actual,
                harmonic.pressure_pa.re * c + harmonic.pressure_pa.im * s,
                5e-11,
                1e-8,
            );
            let direct = base.prepare_receiver(time, 1e-14, &mut || false).unwrap();
            let unchanged = transient
                .prepare_path(0, time, 1e-14, &mut || false)
                .unwrap();
            assert_eq!(
                direct
                    .evaluate(harmonic_strength(phasor, omega, direct.source_time_s()))
                    .unwrap(),
                unchanged
                    .evaluate(harmonic_strength(phasor, omega, unchanged.source_time_s()))
                    .unwrap()
            );
            for mask in 0..8 {
                let image =
                    core::array::from_fn(|i| if mask & (1 << i) == 0 { y[i] } else { -y[i] });
                let expected_delay = distance(x, image) / 300.0;
                let path = transient
                    .prepare_path(mask, time, 1e-14, &mut || false)
                    .unwrap();
                near(path.source_time_s(), time - expected_delay, 2e-12, 0.1);
                near(path.propagation_delay_s(), expected_delay, 2e-14, 0.1);
                assert_eq!(path.receiver_time_s(), time);
            }
        }
    }
}

fn image_vector(v: [f64; 3], mask: usize) -> [f64; 3] {
    core::array::from_fn(|i| if mask & (1 << i) == 0 { v[i] } else { -v[i] })
}

fn image_time(x: [f64; 3], time: f64, mask: usize) -> f64 {
    let mut lo = -0.2_f64;
    let mut hi = 0.4;
    for _ in 0..100 {
        let tau = lo.midpoint(hi);
        if tau + distance(x, image_vector(motion(tau).position_m, mask)) / 300.0 < time {
            lo = tau;
        } else {
            hi = tau;
        }
    }
    lo.midpoint(hi)
}

// An independent SOURCE-image Green potential, with explicit force AND motion
// reflection. The implementation instead uses unchanged sources and receivers
// mirrored at preparation. Numerical differentiation supplies the pressure law.
fn image_potential(x: [f64; 3], time: f64, mask: usize) -> [f64; 4] {
    let tau = image_time(x, time, mask);
    let m = motion(tau);
    let y = image_vector(m.position_m, mask);
    let velocity = image_vector(m.velocity_m_s, mask);
    let r = distance(x, y);
    let radial: f64 = (0..3).map(|i| (x[i] - y[i]) * velocity[i] / r).sum();
    let d = 4.0 * PI * r * (1.0 - radial / 300.0);
    let s = strength(tau);
    let f = image_vector(s.force_on_fluid_n, mask);
    [s.volume_flow_m3_s / d, f[0] / d, f[1] / d, f[2] / d]
}

#[test]
fn accelerating_images_match_independently_differentiated_retarded_potentials() {
    let x = [15.0, 9.0, 8.0];
    let base = source_path(x);
    let images = IdealReflectionObserver::new(*base.observer(), &walls(2, 2), 4).unwrap();
    let reflected =
        PreparedReflectedSourceTrajectory::new(&base, &images, reflection_budget(), &mut || false)
            .unwrap();
    for time in [0.03, 0.12, 0.25] {
        let mut expected = 0.0;
        for mask in 0..4 {
            let sign = if mask & 2 == 0 { 1.0 } else { -1.0 };
            let monopole = 1.2 * derivative(|dt| image_potential(x, time + dt, mask)[0], 2e-5);
            let mut dipole = 0.0;
            for axis in 0..3 {
                dipole -= derivative(
                    |dx| {
                        let mut shifted = x;
                        shifted[axis] += dx;
                        image_potential(shifted, time, mask)[axis + 1]
                    },
                    2e-3,
                );
            }
            let path = reflected
                .prepare_path(mask, time, 1e-13, &mut || false)
                .unwrap();
            let observed = path.evaluate(strength(path.source_time_s())).unwrap();
            near(observed.monopole_pressure_pa, sign * monopole, 4e-7, 1e-9);
            near(observed.dipole_pressure_pa, sign * dipole, 4e-7, 1e-9);
            near(path.source_time_s(), image_time(x, time, mask), 2e-11, 0.1);
            let rate = derivative(|dt| image_time(x, time + dt, mask), 1e-5);
            near(path.source_time_rate(), rate, 2e-8, 1.0);
            expected += sign * (monopole + dipole);
        }
        let actual = reflected
            .pressure_at(time, 1e-13, &mut |tau| Ok(strength(tau)), &mut || false)
            .unwrap();
        near(actual, expected, 1e-6, 1e-8);
    }
}

#[test]
fn moving_sources_satisfy_release_and_rigid_boundaries_for_all_source_components() {
    let time = 0.1;
    let query = |x: [f64; 3], parity: usize, component: usize| {
        let base = source_path(x);
        let images = IdealReflectionObserver::new(*base.observer(), &walls(parity, 1), 2).unwrap();
        let reflected = PreparedReflectedSourceTrajectory::new(
            &base,
            &images,
            reflection_budget(),
            &mut || false,
        )
        .unwrap();
        reflected
            .pressure_at(
                time,
                1e-14,
                &mut |t| {
                    let mut s = CompactSourceStrength::default();
                    if component == 0 {
                        s.volume_flow_m3_s = (17.0 * t).sin() * 0.01;
                        s.volume_flow_rate_m3_s2 = (17.0 * t).cos() * 0.17;
                    } else {
                        s.force_on_fluid_n[component - 1] = (17.0 * t).sin();
                        s.force_rate_on_fluid_n_s[component - 1] = 17.0 * (17.0 * t).cos();
                    }
                    Ok(s)
                },
                &mut || false,
            )
            .unwrap()
    };
    for component in 0..4 {
        assert_eq!(query([0.0, 9.0, 8.0], 1, component), 0.0);
        let h = 0.002;
        let p = (0..4)
            .map(|i| query([f64::from(i) * h, 9.0, 8.0], 0, component))
            .collect::<Vec<_>>();
        let normal_derivative = (-11.0 * p[0] + 18.0 * p[1] - 9.0 * p[2] + 2.0 * p[3]) / (6.0 * h);
        assert!(
            normal_derivative.abs() <= 2e-7 * p[0].abs().max(1e-7),
            "{component}: {normal_derivative:e} {p:?}"
        );
    }
}

#[test]
fn admission_checks_the_whole_source_path_not_just_knots() {
    let x = [8.0, 8.0, 8.0];
    let microphone = observer(x, 300.0, 0.8, 0.05);
    let images = IdealReflectionObserver::new(microphone, &walls(0, 1), 2).unwrap();
    let polynomial = |minimum: f64, t: f64| SourceKinematics {
        position_m: [minimum + (t - 0.5).powi(2), 2.0, 3.0],
        velocity_m_s: [2.0 * (t - 0.5), 0.0, 0.0],
        acceleration_m_s2: [2.0, 0.0, 0.0],
    };
    let crossing = prepare(microphone, &samples(&[0.0, 1.0], |t| polynomial(-0.1, t)));
    assert!(
        PreparedReflectedSourceTrajectory::new(
            &crossing,
            &images,
            reflection_budget(),
            &mut || false
        )
        .is_err()
    );
    let safe = prepare(microphone, &samples(&[0.0, 1.0], |t| polynomial(0.01, t)));
    assert!(
        PreparedReflectedSourceTrajectory::new(
            &safe,
            &images,
            ReflectedTrajectoryBudget {
                maximum_subdivision_depth: 0,
                ..reflection_budget()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedReflectedSourceTrajectory::new(&safe, &images, reflection_budget(), &mut || false)
            .is_ok()
    );
    assert!(
        PreparedReflectedSourceTrajectory::new(
            &safe,
            &images,
            ReflectedTrajectoryBudget {
                maximum_bound_nodes: 0,
                ..reflection_budget()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedReflectedSourceTrajectory::new(&safe, &images, reflection_budget(), &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(IdealReflectionObserver::new(microphone, &walls(0, 3), 7).is_err());
    assert!(IdealReflectionObserver::new(microphone, &[], 8).is_err());
    assert!(IdealReflectionObserver::new(microphone, &[walls(0, 1)[0]; 2], 4).is_err());
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1], v[2]]
}

#[test]
fn reflected_accelerating_sources_keep_coordinate_covariance_and_original_history() {
    let x = [15.0, 9.0, 8.0];
    let base = source_path(x);
    let planes = walls(1, 3);
    let images = IdealReflectionObserver::new(*base.observer(), &planes, 8).unwrap();
    let reflected =
        PreparedReflectedSourceTrajectory::new(&base, &images, reflection_budget(), &mut || false)
            .unwrap();
    let shift = [21.0, -13.0, 7.0];
    let point = |v: [f64; 3]| {
        let v = rotate(v);
        core::array::from_fn(|i| v[i] + shift[i])
    };
    let transformed = prepare(
        observer(point(x), 300.0, 0.8, 0.05),
        &samples(&[-0.2, -0.03, 0.11, 0.4], |t| {
            let state = motion(t);
            SourceKinematics {
                position_m: point(state.position_m),
                velocity_m_s: rotate(state.velocity_m_s),
                acceleration_m_s2: rotate(state.acceleration_m_s2),
            }
        }),
    );
    let planes = planes
        .iter()
        .map(|p| {
            IdealAcousticPlane::new(point(p.point_m()), rotate(p.inward_normal()), p.boundary())
                .unwrap()
        })
        .collect::<Vec<_>>();
    let transformed_images =
        IdealReflectionObserver::new(*transformed.observer(), &planes, 8).unwrap();
    let transformed_reflection = PreparedReflectedSourceTrajectory::new(
        &transformed,
        &transformed_images,
        reflection_budget(),
        &mut || false,
    )
    .unwrap();
    for time in [0.03, 0.1, 0.25] {
        let expected = reflected
            .pressure_at(time, 1e-14, &mut |t| Ok(strength(t)), &mut || false)
            .unwrap();
        let actual = transformed_reflection
            .pressure_at(
                time,
                1e-14,
                &mut |t| {
                    let mut s = strength(t);
                    s.force_on_fluid_n = rotate(s.force_on_fluid_n);
                    s.force_rate_on_fluid_n_s = rotate(s.force_rate_on_fluid_n_s);
                    Ok(s)
                },
                &mut || false,
            )
            .unwrap();
        near(actual, expected, 2e-10, 1e-8);
    }
}

#[test]
fn reflection_refuses_missing_history_and_cancelled_queries_without_changing_the_source() {
    let base = source_path([15.0, 9.0, 8.0]);
    let images = IdealReflectionObserver::new(*base.observer(), &walls(1, 3), 8).unwrap();
    let reflected =
        PreparedReflectedSourceTrajectory::new(&base, &images, reflection_budget(), &mut || false)
            .unwrap();
    let [first, last] = reflected.receiver_time_interval_s();
    assert!(first > base.receiver_time_interval_s()[0]);
    for time in [first.next_down(), last.next_up(), f64::NAN] {
        assert!(
            reflected
                .pressure_at(time, 1e-12, &mut |t| Ok(strength(t)), &mut || false)
                .is_err()
        );
    }
    assert!(
        reflected
            .prepare_path(8, 0.1, 1e-12, &mut || false)
            .is_err()
    );
    let mut count = 0;
    let expected = reflected
        .pressure_at(0.1, 1e-13, &mut |t| Ok(strength(t)), &mut || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        reflected.pressure_at(0.1, 1e-13, &mut |t| Ok(strength(t)), &mut || {
            calls += 1;
            calls == count / 2
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        reflected
            .pressure_at(0.1, 1e-13, &mut |t| Ok(strength(t)), &mut || false)
            .unwrap()
            .to_bits(),
        expected.to_bits()
    );
    let mut images_read = 0;
    assert!(
        reflected
            .pressure_at(
                0.1,
                1e-13,
                &mut |t| {
                    images_read += 1;
                    if images_read == 2 {
                        Err(AeroacError::InvalidParameter {
                            what: "unsupported original history",
                        })
                    } else {
                        Ok(strength(t))
                    }
                },
                &mut || false
            )
            .is_err()
    );
}
