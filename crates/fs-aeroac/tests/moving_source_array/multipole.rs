//! G1/G3/G4: original mixed mechanisms share causal events and coherent fields.

use super::*;
use fs_aeroac::curle3d::moving::trajectory::{MovingMultipoleHistory, MovingMultipoleSource};
use fs_aeroac::curle3d::quadrupole::transient::TransientQuadrupoleSource;
use fs_math::dd::Dd;

fn stress(index: usize, time: f64) -> TransientQuadrupoleSource {
    let a = [0.7, -0.2, 0.4, 0.3, -0.1, 0.23];
    let b = [-0.11, 0.13, 0.09, -0.05, 0.07, 0.02];
    let scale = [1.0, -0.7, 0.3][index];
    let s = time + 1.0;
    TransientQuadrupoleSource {
        stress_integral_n_m_s: core::array::from_fn(|i| {
            scale * (a[i] * s.powi(4) / 4.0 + b[i] * s.powi(5) / 5.0)
        }),
        stress_moment_n_m: core::array::from_fn(|i| scale * (a[i] * s.powi(3) + b[i] * s.powi(4))),
        stress_rate_n_m_s: core::array::from_fn(|i| {
            scale * (3.0 * a[i] * s * s + 4.0 * b[i] * s.powi(3))
        }),
        stress_acceleration_n_m_s2: core::array::from_fn(|i| {
            scale * (6.0 * a[i] * s + 12.0 * b[i] * s * s)
        }),
    }
}

fn strength(index: usize, time: f64) -> MovingMultipoleSource {
    MovingMultipoleSource {
        compact: drive(index, time),
        stress: stress(index, time),
    }
}

fn qn(q: [f64; 6], n: [f64; 3]) -> [f64; 3] {
    let q = [[q[0], q[3], q[4]], [q[3], q[1], q[5]], [q[4], q[5], q[2]]];
    q.map(|row| dot(row, n))
}
fn stress_shape(q: [f64; 6], n: [f64; 3], a: f64, b: f64, c: f64) -> [f64; 3] {
    let qn = qn(q, n);
    let radial = a * dot(n, qn) - c * (q[0] + q[1] + q[2]);
    core::array::from_fn(|i| radial * n[i] - b * qn[i])
}

// Independent closed-form spatial derivatives of 1/(4*pi*r), including every
// induction order. These do not read native moving coefficients or image maps.
fn stationary_field(
    position: [f64; 3],
    microphone: [f64; 3],
    source: MovingMultipoleSource,
    force_impulse: [f64; 3],
) -> (f64, [f64; 3]) {
    let r_vec = core::array::from_fn(|i| microphone[i] - position[i]);
    let r = norm(r_vec);
    let n = r_vec.map(|v| v / r);
    let c = medium().sound_speed_m_s();
    let rho = medium().density_kg_m3();
    let f = source.compact;
    let q = source.stress;
    let radial = |q: [f64; 6]| 3.0 * dot(n, qn(q, n)) - q[0] - q[1] - q[2];
    let p = (rho * f.volume_flow_rate_m3_s2 / r
        + dot(n, f.force_on_fluid_n) / r.powi(2)
        + dot(n, f.force_rate_on_fluid_n_s) / (c * r)
        + radial(q.stress_moment_n_m) / r.powi(3)
        + radial(q.stress_rate_n_m_s) / (c * r.powi(2))
        + dot(n, qn(q.stress_acceleration_n_m_s2, n)) / (c * c * r))
        / (4.0 * PI);
    let h = stress_shape(q.stress_integral_n_m_s, n, 15.0, 6.0, 3.0);
    let t = stress_shape(q.stress_moment_n_m, n, 15.0, 6.0, 3.0);
    let dt = stress_shape(q.stress_rate_n_m_s, n, 6.0, 2.0, 1.0);
    let ddt = stress_shape(q.stress_acceleration_n_m_s2, n, 1.0, 0.0, 0.0);
    let u = core::array::from_fn(|i| {
        (n[i] * (f.volume_flow_m3_s / r.powi(2) + f.volume_flow_rate_m3_s2 / (c * r))
            + ((3.0 * n[i] * dot(n, force_impulse) - force_impulse[i]) / r.powi(3)
                + (3.0 * n[i] * dot(n, f.force_on_fluid_n) - f.force_on_fluid_n[i])
                    / (c * r.powi(2))
                + n[i] * dot(n, f.force_rate_on_fluid_n_s) / (c * c * r)
                + h[i] / r.powi(4)
                + t[i] / (c * r.powi(3))
                + dt[i] / (c * c * r.powi(2))
                + ddt[i] / (c.powi(3) * r))
                / rho)
            / (4.0 * PI)
    });
    (p, u)
}

fn stationary_history(
    position: [f64; 3],
    microphone: [f64; 3],
    force_impulse: [f64; 3],
    stress_integral: [f64; 6],
) -> MovingMultipoleHistory {
    let displacement = core::array::from_fn(|i| microphone[i] - position[i]);
    let r = norm(displacement);
    let n = displacement.map(|v| v / r);
    MovingMultipoleHistory {
        force_n_s_m3: core::array::from_fn(|i| {
            Dd::from_f64((3.0 * n[i] * dot(n, force_impulse) - force_impulse[i]) / r.powi(3))
        }),
        stress_n_s_m3: stress_shape(stress_integral, n, 15.0, 6.0, 3.0)
            .map(|v| Dd::from_f64(v / r.powi(4))),
    }
}

fn moving_history(
    index: usize,
    observer: &MovingCompactObserver,
    time: f64,
) -> MovingMultipoleHistory {
    MovingMultipoleHistory {
        force_n_s_m3: history(index, observer, time).map(Dd::from_f64),
        stress_n_s_m3: integrate(time, |s| {
            observer
                .quadrupole_history_integrand(
                    motion(index, s).position_m,
                    stress(index, s).stress_moment_n_m,
                )
                .unwrap()
        })
        .map(Dd::from_f64),
    }
}

fn tensor_image(q: [f64; 6], mask: usize) -> [f64; 6] {
    let axes = image([1.0; 3], mask);
    [
        q[0],
        q[1],
        q[2],
        q[3] * axes[0] * axes[1],
        q[4] * axes[0] * axes[2],
        q[5] * axes[1] * axes[2],
    ]
}

fn check_field(actual: MovingSourceArrayField, expected: (f64, [f64; 3]), tolerance: f64) {
    near(actual.pressure_pa, expected.0, tolerance, 1e-100);
    let scale = norm(expected.1).max(1e-100);
    for i in 0..3 {
        assert!((actual.particle_velocity_m_s[i] - expected.1[i]).abs() <= tolerance * scale);
        assert!(
            (actual.instantaneous_intensity_w_m2[i] - expected.0 * expected.1[i]).abs()
                <= 3.0 * tolerance * expected.0.abs() * scale
        );
    }
}

#[test]
fn moving_microphone_and_complete_images_match_independent_stationary_multipoles() {
    let positions: Vec<_> = (0..3).map(|i| motion(i, 0.0).position_m).collect();
    let sources: Vec<_> = positions.iter().map(|&p| path(|_| fixed(p))).collect();
    let microphone = path(receiver);
    let direct: Vec<_> = sources
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
    let images =
        MovingReflectionObserver::new(&microphone, &walls(2, 1), 4, bounds(), &mut || false)
            .unwrap();
    let reflected: Vec<_> = sources
        .iter()
        .map(|source| {
            PreparedReflectedMovingReceiverSource::new(
                source,
                &images,
                config([0.2, 0.5]),
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
    for (array, image_count, parity) in [(&free, 1, 0), (&corner, 4, 1)] {
        for time in [0.25, 0.45] {
            let actual = array
                .multipole_field_at(
                    time,
                    &mut |i, t| Ok(strength(i, t)),
                    &mut |i, _, observer, t| {
                        Ok(stationary_history(
                            positions[i],
                            observer.observer_position_m(),
                            integrate(t, |s| drive(i, s).force_on_fluid_n),
                            stress(i, t).stress_integral_n_m_s,
                        ))
                    },
                    &mut || false,
                )
                .unwrap();
            let x = receiver(time).position_m;
            let mut expected = (0.0, [0.0; 3]);
            for (i, &y) in positions.iter().enumerate() {
                for mask in 0..image_count {
                    let y = image(y, mask);
                    let tau = time - norm(core::array::from_fn(|axis| x[axis] - y[axis])) / 20.0;
                    let mut value = strength(i, tau);
                    value.compact.force_on_fluid_n = image(value.compact.force_on_fluid_n, mask);
                    value.compact.force_rate_on_fluid_n_s =
                        image(value.compact.force_rate_on_fluid_n_s, mask);
                    value.stress.stress_integral_n_m_s =
                        tensor_image(value.stress.stress_integral_n_m_s, mask);
                    value.stress.stress_moment_n_m =
                        tensor_image(value.stress.stress_moment_n_m, mask);
                    value.stress.stress_rate_n_m_s =
                        tensor_image(value.stress.stress_rate_n_m_s, mask);
                    value.stress.stress_acceleration_n_m_s2 =
                        tensor_image(value.stress.stress_acceleration_n_m_s2, mask);
                    let impulse = image(integrate(tau, |s| drive(i, s).force_on_fluid_n), mask);
                    let field = stationary_field(y, x, value, impulse);
                    expected.0 += sign(mask, parity) * field.0;
                    for axis in 0..3 {
                        expected.1[axis] += sign(mask, parity) * field.1[axis];
                    }
                }
            }
            assert_eq!(actual.source_path_count, 3 * image_count);
            check_field(actual, expected, 2e-10);
        }
    }
}

#[test]
fn accelerating_sources_in_a_wedge_share_native_clocks_and_full_coherent_intensity() {
    let sources: Vec<_> = (0..3).map(|i| path(|t| motion(i, t))).collect();
    let microphone = path(receiver);
    let angle = PI / 3.0;
    let planes = [
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
        IdealAcousticPlane::new(
            [0.0; 3],
            [angle.sin(), -angle.cos(), 0.0],
            IdealPlaneBoundary::Rigid,
        )
        .unwrap(),
    ];
    let images =
        MovingReflectionObserver::new_wedge(&microphone, planes, 3, 6, bounds(), &mut || false)
            .unwrap();
    let reflected: Vec<_> = sources
        .iter()
        .map(|source| {
            PreparedReflectedMovingReceiverSource::new(
                source,
                &images,
                config([0.2, 0.5]),
                bounds(),
                pairs(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let array =
        PreparedMovingSourceArray::new_reflected(&reflected, budget(), &mut || false).unwrap();
    assert_eq!(array.maximum_retarded_iterations(), 18 * 96);
    for time in [0.25, 0.45] {
        let mut source_events = Vec::new();
        let mut history_events = Vec::new();
        let actual = array
            .multipole_field_at(
                time,
                &mut |i, t| {
                    source_events.push((i, t));
                    Ok(strength(i, t))
                },
                &mut |i, path, observer, t| {
                    history_events.push((i, path, t));
                    assert_eq!(
                        observer.observer_position_m(),
                        reflected[i]
                            .force_history_path_at(path, time)
                            .unwrap()
                            .0
                            .observer_position_m()
                    );
                    Ok(moving_history(i, observer, t))
                },
                &mut || false,
            )
            .unwrap();
        assert_eq!(source_events.len(), 18);
        assert_eq!(history_events.len(), 18);
        for (slot, (&(i, t), &(j, path, s))) in
            source_events.iter().zip(&history_events).enumerate()
        {
            assert_eq!((i, j, path), (slot / 6, slot / 6, slot % 6));
            assert_eq!(t, s);
        }
        assert_ne!(source_events[0].1, source_events[17].1);
        let compact = array
            .force_field_at(
                time,
                &mut |i, t| Ok(drive(i, t)),
                &mut |i, _, observer, t| Ok(history(i, observer, t)),
                &mut || false,
            )
            .unwrap();
        let mut p = Dd::from_f64(compact.pressure_pa);
        let mut u = compact.particle_velocity_m_s.map(Dd::from_f64);
        let mut separate_intensity = compact.instantaneous_intensity_w_m2;
        for (i, source) in reflected.iter().enumerate() {
            for path in 0..source.path_count() {
                let field = source
                    .prepare_quadrupole_field_path(path, time, &mut || false)
                    .unwrap();
                let observer = source.force_history_path_at(path, time).unwrap().0;
                let history = moving_history(i, &observer, field.source_time_s());
                let (pp, uu) = field
                    .evaluate_dd(stress(i, field.source_time_s()), history.stress_n_s_m3)
                    .unwrap();
                p = p + pp;
                for axis in 0..3 {
                    u[axis] = u[axis] + uu[axis];
                    separate_intensity[axis] += (pp * uu[axis]).to_f64();
                }
            }
        }
        check_field(actual, (p.to_f64(), u.map(Dd::to_f64)), 5e-11);
        assert!((0..3).any(|i| {
            (actual.instantaneous_intensity_w_m2[i] - separate_intensity[i]).abs() > 1e-10
        }));
        assert!(actual.maximum_time_residual_s <= 1e-12);
    }
}

#[test]
fn stress_survives_volume_cancellation_before_source_rounding() {
    let y = [1.0, 0.0, 0.0];
    let x = [5.0, 0.0, 0.0];
    let source = path(|_| fixed(y));
    let microphone = path(|_| fixed(x));
    let pairs: Vec<_> = (0..2)
        .map(|_| {
            PreparedMovingReceiverSource::new(
                &source,
                &microphone,
                config([0.2, 0.5]),
                pairs(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let array = PreparedMovingSourceArray::new(&pairs, budget(), &mut || false).unwrap();
    let weak = |time| {
        let value = stress(0, time);
        let scale = 2.0_f64.powi(-60);
        TransientQuadrupoleSource {
            stress_integral_n_m_s: value.stress_integral_n_m_s.map(|v| v * scale),
            stress_moment_n_m: value.stress_moment_n_m.map(|v| v * scale),
            stress_rate_n_m_s: value.stress_rate_n_m_s.map(|v| v * scale),
            stress_acceleration_n_m_s2: value.stress_acceleration_n_m_s2.map(|v| v * scale),
        }
    };
    let value = |i, t: f64| MovingMultipoleSource {
        compact: CompactSourceStrength {
            volume_flow_m3_s: if i == 0 { 1.0 } else { -1.0 } * (t + 1.0).powi(2),
            volume_flow_rate_m3_s2: if i == 0 { 1.0 } else { -1.0 } * 2.0 * (t + 1.0),
            ..CompactSourceStrength::default()
        },
        stress: if i == 0 {
            weak(t)
        } else {
            TransientQuadrupoleSource::default()
        },
    };
    let histories = |i, observer: &MovingCompactObserver, t| {
        stationary_history(
            y,
            observer.observer_position_m(),
            [0.0; 3],
            if i == 0 {
                weak(t).stress_integral_n_m_s
            } else {
                [0.0; 6]
            },
        )
    };
    let time = 0.3;
    let actual = array
        .multipole_field_at(
            time,
            &mut |i, t| Ok(value(i, t)),
            &mut |i, _, observer, t| Ok(histories(i, observer, t)),
            &mut || false,
        )
        .unwrap();
    let first = PreparedMovingSourceArray::new(&pairs[..1], budget(), &mut || false)
        .unwrap()
        .multipole_field_at(
            time,
            &mut |_, t| Ok(value(0, t)),
            &mut |_, _, observer, t| Ok(histories(0, observer, t)),
            &mut || false,
        )
        .unwrap();
    let second = PreparedMovingSourceArray::new(&pairs[1..], budget(), &mut || false)
        .unwrap()
        .multipole_field_at(
            time,
            &mut |_, t| Ok(value(1, t)),
            &mut |_, _, observer, t| Ok(histories(1, observer, t)),
            &mut || false,
        )
        .unwrap();
    assert_eq!(first.pressure_pa + second.pressure_pa, 0.0);
    let expected = stationary_field(
        y,
        x,
        MovingMultipoleSource {
            compact: CompactSourceStrength::default(),
            stress: weak(time - 4.0 / 20.0),
        },
        [0.0; 3],
    );
    assert_ne!(expected.0, 0.0);
    check_field(actual, expected, 3e-11);
}

#[test]
fn mixed_fields_obey_rigid_and_release_boundaries_without_erasing_velocity() {
    let sources: Vec<_> = (0..2).map(|i| path(|t| motion(i, t))).collect();
    let microphone = path(|t| SourceKinematics {
        position_m: [0.0, 3.0 + 0.2 * t, 2.0 + 0.1 * t],
        velocity_m_s: [0.0, 0.2, 0.1],
        acceleration_m_s2: [0.0; 3],
    });
    for parity in 0..2 {
        let images =
            MovingReflectionObserver::new(&microphone, &walls(1, parity), 2, bounds(), &mut || {
                false
            })
            .unwrap();
        let reflected: Vec<_> = sources
            .iter()
            .map(|source| {
                PreparedReflectedMovingReceiverSource::new(
                    source,
                    &images,
                    config([0.2, 0.5]),
                    bounds(),
                    pairs(),
                    &mut || false,
                )
                .unwrap()
            })
            .collect();
        let array =
            PreparedMovingSourceArray::new_reflected(&reflected, budget(), &mut || false).unwrap();
        let field = array
            .multipole_field_at(
                0.3,
                &mut |i, t| Ok(strength(i, t)),
                &mut |i, _, observer, t| Ok(moving_history(i, observer, t)),
                &mut || false,
            )
            .unwrap();
        assert_eq!(field.instantaneous_intensity_w_m2[0], 0.0);
        if parity == 0 {
            assert_eq!(field.particle_velocity_m_s[0], 0.0);
            assert!(field.pressure_pa.abs() > 1e-8);
        } else {
            assert_eq!(field.pressure_pa, 0.0);
            assert!(field.particle_velocity_m_s[0].abs() > 1e-8);
            assert_eq!(&field.particle_velocity_m_s[1..], &[0.0; 2]);
        }
    }
}

#[test]
fn zero_stress_recovers_native_compact_fields_and_failures_publish_no_partial_result() {
    let sources: Vec<_> = (0..2).map(|i| path(|t| motion(i, t))).collect();
    let microphone = path(receiver);
    let pairs: Vec<_> = sources
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
    let array = PreparedMovingSourceArray::new(&pairs, budget(), &mut || false).unwrap();
    let original = array
        .force_field_at(
            0.3,
            &mut |i, t| Ok(drive(i, t)),
            &mut |i, _, observer, t| Ok(history(i, observer, t)),
            &mut || false,
        )
        .unwrap();
    let zero_stress = array
        .multipole_field_at(
            0.3,
            &mut |i, t| {
                Ok(MovingMultipoleSource {
                    compact: drive(i, t),
                    ..MovingMultipoleSource::default()
                })
            },
            &mut |i, _, observer, t| {
                Ok(MovingMultipoleHistory {
                    force_n_s_m3: history(i, observer, t).map(Dd::from_f64),
                    ..MovingMultipoleHistory::default()
                })
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(zero_stress, original);
    let mut calls = 0;
    assert!(
        array
            .multipole_field_at(
                0.1,
                &mut |i, t| {
                    calls += 1;
                    Ok(strength(i, t))
                },
                &mut |_, _, _, _| Ok(MovingMultipoleHistory::default()),
                &mut || false
            )
            .is_err()
    );
    assert_eq!(calls, 0);
    let missing = AeroacError::InvalidParameter {
        what: "last stress history is missing",
    };
    assert_eq!(
        array
            .multipole_field_at(
                0.3,
                &mut |i, t| Ok(strength(i, t)),
                &mut |i, _, observer, t| if i == 1 {
                    Err(missing.clone())
                } else {
                    Ok(moving_history(i, observer, t))
                },
                &mut || false
            )
            .unwrap_err(),
        missing
    );
    for bad_history in [false, true] {
        assert!(matches!(
            array.multipole_field_at(
                0.3,
                &mut |i, t| {
                    let mut value = strength(i, t);
                    if i == 1 && !bad_history {
                        value.stress.stress_rate_n_m_s[4] = f64::NAN;
                    }
                    Ok(value)
                },
                &mut |i, _, observer, t| {
                    let mut value = moving_history(i, observer, t);
                    if i == 1 && bad_history {
                        value.stress_n_s_m3[2] = Dd::from_f64(f64::NAN);
                    }
                    Ok(value)
                },
                &mut || false
            ),
            Err(AeroacError::NonFinite { .. })
        ));
    }
    let cancelled = std::cell::Cell::new(false);
    assert!(matches!(
        array.multipole_field_at(
            0.3,
            &mut |i, t| Ok(strength(i, t)),
            &mut |i, _, observer, t| {
                let value = moving_history(i, observer, t);
                if i == 1 {
                    cancelled.set(true);
                }
                Ok(value)
            },
            &mut || cancelled.get()
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let evaluate = || {
        array
            .multipole_field_at(
                0.3,
                &mut |i, t| Ok(strength(i, t)),
                &mut |i, _, observer, t| Ok(moving_history(i, observer, t)),
                &mut || false,
            )
            .unwrap()
    };
    assert_eq!(evaluate(), evaluate());
}
