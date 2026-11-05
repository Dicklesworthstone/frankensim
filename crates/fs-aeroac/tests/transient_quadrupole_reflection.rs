//! G1/G3/G4: causal tensor images, full velocity, physical histories and bounds.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::array::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, RetardedQuadrupoleField, TransientQuadrupoleSource,
};
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use std::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.23, 320.0).unwrap()
}

fn observer(x: [f64; 3], sources: &[[f64; 3]]) -> PreparedTransientQuadrupoleObserver {
    PreparedTransientQuadrupoleObserver::new(medium(), x, sources, sources.len()).unwrap()
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.abs().max(1e-25),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

// A causal C2 polynomial history, with its original integral and derivatives.
fn ramp(time: f64, tensor: [f64; 6]) -> TransientQuadrupoleSource {
    if time <= 0.0 {
        return TransientQuadrupoleSource::default();
    }
    let duration = 0.01;
    let s = time / duration;
    TransientQuadrupoleSource {
        stress_integral_n_m_s: tensor.map(|v| v * duration * s.powi(4) / 4.0),
        stress_moment_n_m: tensor.map(|v| v * s.powi(3)),
        stress_rate_n_m_s: tensor.map(|v| v * 3.0 * s.powi(2) / duration),
        stress_acceleration_n_m_s2: tensor.map(|v| v * 6.0 * s / (duration * duration)),
    }
}

fn tensor() -> [f64; 6] {
    [0.7, -0.3, 0.9, 0.23, -0.11, 0.4]
}

// Independent physical source-image construction: Q_image = s H Q H^T.
fn tensor_image(q: [f64; 6], matrix: [[f64; 3]; 3], sign: f64) -> [f64; 6] {
    let [xx, yy, zz, xy, xz, yz] = q;
    let q = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
    [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)].map(|(i, j)| {
        let mut value = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                value += matrix[i][a] * q[a][b] * matrix[j][b];
            }
        }
        sign * value
    })
}

fn check_sum(
    actual: RetardedQuadrupoleField,
    a: RetardedQuadrupoleField,
    b: RetardedQuadrupoleField,
) {
    close(
        actual.pressure_pa,
        a.pressure_pa + b.pressure_pa,
        a.pressure_pa.abs() + b.pressure_pa.abs(),
        4e-13,
    );
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            a.particle_velocity_m_s[axis] + b.particle_velocity_m_s[axis],
            a.particle_velocity_m_s[axis].abs() + b.particle_velocity_m_s[axis].abs(),
            4e-13,
        );
        assert_eq!(
            actual.instantaneous_intensity_w_m2[axis],
            actual.pressure_pa * actual.particle_velocity_m_s[axis]
        );
    }
}

#[test]
fn g1_transient_oblique_wall_preserves_tensor_images_and_distinct_arrivals() {
    let point = [0.03, -0.02, 0.04];
    let normal = [0.6, 0.0, 0.8];
    let source = [0.2, 0.3, 0.4];
    let receiver = [0.9, 0.7, 0.8];
    let distance: f64 = (0..3).map(|i| (source[i] - point[i]) * normal[i]).sum();
    let image_source = core::array::from_fn(|i| source[i] - 2.0 * distance * normal[i]);
    let matrix = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i == j {
                1.0 - 2.0 * normal[i] * normal[j]
            } else {
                -2.0 * normal[i] * normal[j]
            }
        })
    });
    let direct = observer(receiver, &[source]);
    let image = observer(receiver, &[image_source]);
    let d0 = direct.source_delay_s(0).unwrap();
    let d1 = image.source_delay_s(0).unwrap();
    assert!(d1 > d0);
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new(point, normal, boundary).unwrap();
        let reflected = direct.clone().with_orthogonal_planes(&[plane], 2).unwrap();
        let sign = if boundary == IdealPlaneBoundary::Rigid {
            1.0
        } else {
            -1.0
        };
        let image_tensor = tensor_image(tensor(), matrix, sign);
        assert_eq!(reflected.source_count(), 1);
        assert_eq!(reflected.path_count(), 2);
        assert_eq!(reflected.source_position_m(0), Some(source));
        assert_eq!(reflected.source_position_m(1), None);
        assert_eq!(reflected.source_delay_s(1), None);
        assert_eq!(reflected.path_source_index(0), Some(0));
        assert_eq!(reflected.path_source_index(1), Some(0));
        assert_eq!(reflected.reflection_planes(), &[plane]);
        assert_eq!(reflected.reflection_wedge_order(), None);
        let bounds = reflected.source_delay_bounds_s(0).unwrap();
        close(bounds[0], d0, d0, 3e-15);
        close(bounds[1], d1, d1, 3e-15);
        for time in [0.0, 0.5 * d0, 0.5 * (d0 + d1), d1 + 0.0007, d1 + 0.004] {
            let mut calls = Vec::new();
            let actual = reflected
                .evaluate_at_time(time, |index, emission| {
                    calls.push((index, emission));
                    Ok(ramp(emission, tensor()))
                })
                .unwrap();
            let a = direct
                .evaluate_at_time(time, |_, t| Ok(ramp(t, tensor())))
                .unwrap();
            let b = image
                .evaluate_at_time(time, |_, t| Ok(ramp(t, image_tensor)))
                .unwrap();
            check_sum(actual, a, b);
            assert_eq!(calls.len(), 2);
            for path in 0..2 {
                assert_eq!(calls[path].0, 0);
                assert_eq!(calls[path].1, time - reflected.path_delay_s(path).unwrap());
            }
            if time < d1 {
                assert_eq!(actual.pressure_pa, a.pressure_pa);
                assert_eq!(actual.particle_velocity_m_s, a.particle_velocity_m_s);
            }
            let explicit: Vec<_> = calls.iter().map(|(_, t)| ramp(*t, tensor())).collect();
            assert_eq!(actual, reflected.evaluate_at_retarded(&explicit).unwrap());
        }
        assert!(
            reflected
                .evaluate_at_retarded(&[ramp(0.004, tensor())])
                .is_err()
        );
    }
}

#[test]
fn g1_mixed_corner_enforces_every_tensor_component_at_the_physical_walls() {
    let boundaries = [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
        IdealPlaneBoundary::Rigid,
    ];
    let planes: [_; 3] = core::array::from_fn(|axis| {
        let mut normal = [0.0; 3];
        normal[axis] = 1.0;
        IdealAcousticPlane::new([0.0; 3], normal, boundaries[axis]).unwrap()
    });
    for face in 0..3 {
        let mut receiver = [0.6, 0.7, 0.8];
        receiver[face] = 0.0;
        let direct = observer(receiver, &[[0.2, 0.3, 0.4]]);
        let reflected = direct.clone().with_orthogonal_planes(&planes, 8).unwrap();
        assert_eq!(reflected.path_count(), 8);
        for component in 0..6 {
            let mut q = [0.0; 6];
            q[component] = 1.0;
            for time in [0.004, 0.008, 0.012] {
                let expected = direct
                    .evaluate_at_time(time, |_, t| Ok(ramp(t, q)))
                    .unwrap();
                let field = reflected
                    .evaluate_at_time(time, |_, t| Ok(ramp(t, q)))
                    .unwrap();
                if boundaries[face] == IdealPlaneBoundary::Rigid {
                    let scale = expected
                        .particle_velocity_m_s
                        .iter()
                        .fold(0.0_f64, |s, v| s.max(v.abs()));
                    close(field.particle_velocity_m_s[face], 0.0, scale, 5e-12);
                } else {
                    close(field.pressure_pa, 0.0, expected.pressure_pa, 5e-12);
                }
            }
        }
    }
}

fn real_at(value: C64, omega: f64, time: f64) -> f64 {
    let (s, c) = (omega * time).sin_cos();
    value.re * c + value.im * s
}

fn harmonic(q: [C64; 6], omega: f64, time: f64) -> TransientQuadrupoleSource {
    TransientQuadrupoleSource {
        stress_integral_n_m_s: q
            .map(|v| real_at(C64::new(-v.im / omega, v.re / omega), omega, time)),
        stress_moment_n_m: q.map(|v| real_at(v, omega, time)),
        stress_rate_n_m_s: q.map(|v| real_at(C64::new(omega * v.im, -omega * v.re), omega, time)),
        stress_acceleration_n_m_s2: q.map(|v| real_at(v.scale(-omega * omega), omega, time)),
    }
}

fn wedge() -> [IdealAcousticPlane; 2] {
    let (s, c) = (PI / 4.0).sin_cos();
    [
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
        IdealAcousticPlane::new([0.0; 3], [s, -c, 0.0], IdealPlaneBoundary::PressureRelease)
            .unwrap(),
    ]
}

#[test]
fn g1_wedge_transients_match_harmonics_and_receiver_momentum() {
    let positions = [[0.2, 0.06, -0.1], [0.3, 0.09, 0.2]];
    let receiver = [0.7, 0.15, 0.3];
    let phasors = [
        tensor().map(|q| C64::new(q, 0.2 * q)),
        tensor().map(|q| C64::new(-0.3 * q, 0.7 * q)),
    ];
    let transient = observer(receiver, &positions)
        .with_wedge(wedge(), 4, 16)
        .unwrap();
    assert_eq!(transient.source_count(), 2);
    assert_eq!(transient.path_count(), 16);
    assert_eq!(transient.reflection_wedge_order(), Some(4));
    for frequency in [43.0, 537.0] {
        let omega = TAU * frequency;
        let field = PreparedQuadrupoleObserver::new(medium(), frequency, receiver, &positions, 2)
            .unwrap()
            .with_wedge(wedge(), 4, 16)
            .unwrap()
            .evaluate(&phasors.map(|stress_moment_n_m| QuadrupolePhasor { stress_moment_n_m }))
            .unwrap();
        for sample in 0..12 {
            let time = (0.13 + sample as f64 / 12.0) / frequency;
            let observed = transient
                .evaluate_at_time(time, |i, t| Ok(harmonic(phasors[i], omega, t)))
                .unwrap();
            close(
                observed.pressure_pa,
                real_at(field.pressure_pa, omega, time),
                field.pressure_pa.abs(),
                4e-11,
            );
            for axis in 0..3 {
                close(
                    observed.particle_velocity_m_s[axis],
                    real_at(field.particle_velocity_m_s[axis], omega, time),
                    field.particle_velocity_m_s[axis].abs(),
                    4e-11,
                );
            }
        }
    }
    // Independent receiver-coordinate pressure derivatives require the physical
    // velocity pullback H, not H^-1; the wedge contains non-symmetric rotations.
    let at = |x, t| {
        observer(x, &positions)
            .with_wedge(wedge(), 4, 16)
            .unwrap()
            .evaluate_at_time(t, |i, s| Ok(ramp(s, tensor().map(|q| q * (i + 1) as f64))))
            .unwrap()
    };
    let time = 0.01;
    let dt = 1e-8;
    let h = 2e-5;
    let earlier = at(receiver, time - dt);
    let later = at(receiver, time + dt);
    for axis in 0..3 {
        let mut plus = receiver;
        plus[axis] += h;
        let mut minus = receiver;
        minus[axis] -= h;
        let gradient = (at(plus, time).pressure_pa - at(minus, time).pressure_pa) / (2.0 * h);
        let inertia = medium().density_kg_m3()
            * (later.particle_velocity_m_s[axis] - earlier.particle_velocity_m_s[axis])
            / (2.0 * dt);
        close(inertia, -gradient, gradient, 3e-7);
    }
}

#[test]
fn g4_history_bounds_budgets_and_refusals_preserve_original_source_identity() {
    let positions = [[0.2, 0.3, 0.4], [0.3, 0.2, 0.1]];
    let plane =
        IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], IdealPlaneBoundary::Rigid).unwrap();
    let direct = observer([0.5, 0.6, 0.7], &positions);
    assert!(direct.clone().with_orthogonal_planes(&[plane], 3).is_err());
    assert!(direct.clone().with_orthogonal_planes(&[], 4).is_err());
    assert!(
        direct
            .clone()
            .with_orthogonal_planes(&[plane, plane], 8)
            .is_err()
    );
    assert!(
        observer([0.5, -0.6, 0.7], &positions)
            .with_orthogonal_planes(&[plane], 4)
            .is_err()
    );
    assert!(
        observer([0.5, 0.6, 0.7], &[[0.2, 0.0, 0.4]])
            .with_orthogonal_planes(&[plane], 2)
            .is_err()
    );
    for stop in [1, 3, 6] {
        let mut polls = 0;
        assert!(matches!(
            direct
                .clone()
                .with_orthogonal_planes_cancellable(&[plane], 4, || {
                    polls += 1;
                    polls == stop
                }),
            Err(AeroacError::Cancelled { .. })
        ));
    }
    let reflected = direct.clone().with_orthogonal_planes(&[plane], 4).unwrap();
    assert!(
        reflected
            .clone()
            .with_orthogonal_planes(&[plane], 8)
            .is_err()
    );
    assert_eq!(
        (0..4)
            .map(|p| reflected.path_source_index(p).unwrap())
            .collect::<Vec<_>>(),
        [0, 1, 0, 1]
    );
    for source in 0..2 {
        let bounds = reflected.source_delay_bounds_s(source).unwrap();
        let direct_delay = reflected.path_delay_s(source).unwrap();
        let image_delay = reflected.path_delay_s(source + 2).unwrap();
        assert_eq!(
            bounds,
            [direct_delay.min(image_delay), direct_delay.max(image_delay)]
        );
        assert_eq!(reflected.source_position_m(source), Some(positions[source]));
    }
    assert_eq!(reflected.path_source_index(4), None);
    assert_eq!(reflected.path_delay_s(4), None);
    assert_eq!(reflected.source_delay_bounds_s(2), None);
    let before = reflected
        .evaluate_at_time(0.008, |_, t| Ok(ramp(t, tensor())))
        .unwrap();
    let mut calls = 0;
    let refusal = AeroacError::InvalidParameter {
        what: "missing reflected source history",
    };
    assert_eq!(
        reflected.evaluate_at_time(0.008, |_, t| {
            calls += 1;
            if calls == 3 {
                Err(refusal.clone())
            } else {
                Ok(ramp(t, tensor()))
            }
        }),
        Err(refusal)
    );
    assert!(matches!(
        reflected.evaluate_at_time_cancellable(0.008, |_, _| panic!("cancelled"), || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        reflected
            .evaluate_at_time(0.008, |_, t| Ok(ramp(t, tensor())))
            .unwrap(),
        before
    );
    let empty = observer([0.5, 0.6, 0.7], &[])
        .with_orthogonal_planes(&[plane], 0)
        .unwrap();
    assert_eq!(empty.source_count(), 0);
    assert_eq!(empty.path_count(), 0);
    assert_eq!(empty.path_source_index(0), None);
    assert_eq!(empty.source_delay_bounds_s(0), None);
    let silent = empty
        .evaluate_at_time(0.0, |_, _| panic!("no physical sources"))
        .unwrap();
    assert_eq!(silent.pressure_pa, 0.0);
    assert_eq!(silent.particle_velocity_m_s, [0.0; 3]);
}

#[test]
fn g3_reflected_dd_clock_preserves_arrivals_and_source_identity_at_large_epochs() {
    let positions = [[0.2, 0.06, -0.1], [0.3, 0.09, 0.2]];
    let reflected = observer([0.7, 0.15, 0.3], &positions)
        .with_wedge(wedge(), 4, 16)
        .unwrap();
    // At 2^48 seconds, an f64 timestamp is spaced by 0.0625 s. Every path
    // delay and the 48 kHz sample interval are smaller, but remain live in DD.
    let epoch = Dd::from_f64(2.0_f64.powi(48));
    let source_tensor = |i| tensor().map(|q| q * (i + 1) as f64);
    for relative_time in [0.001, 0.003, 0.008, 0.008 + 1.0 / 48_000.0] {
        let receiver_time = epoch + Dd::from_f64(relative_time);
        let explicit: Vec<_> = (0..reflected.path_count())
            .map(|path| {
                ramp(
                    relative_time - reflected.path_delay_s(path).unwrap(),
                    source_tensor(reflected.path_source_index(path).unwrap()),
                )
            })
            .collect();
        let expected = reflected.evaluate_at_retarded(&explicit).unwrap();
        let mut callbacks = Vec::new();
        let actual = reflected
            .evaluate_at_clock_dd(receiver_time, |physical_source, emission_time| {
                callbacks.push((physical_source, emission_time));
                Ok(ramp(
                    (emission_time - epoch).to_f64(),
                    source_tensor(physical_source),
                ))
            })
            .unwrap();
        assert_eq!(callbacks.len(), reflected.path_count());
        for (path, (physical_source, emission_time)) in callbacks.into_iter().enumerate() {
            assert_eq!(Some(physical_source), reflected.path_source_index(path));
            assert!(emission_time.lt(receiver_time));
            let delay = reflected.path_delay_s(path).unwrap();
            close(
                (receiver_time - emission_time).to_f64(),
                delay,
                delay,
                3e-15,
            );
        }
        close(
            actual.pressure_pa,
            expected.pressure_pa,
            expected.pressure_pa,
            4e-13,
        );
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
                4e-13,
            );
            assert_eq!(
                actual.instantaneous_intensity_w_m2[axis],
                actual.pressure_pa * actual.particle_velocity_m_s[axis]
            );
        }
    }
    assert!(
        reflected
            .evaluate_at_time(epoch.to_f64(), |_, _| panic!("f64 delay is unresolved"))
            .is_err()
    );
}
