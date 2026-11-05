//! G1/G3/G4 checks of causal stress radiation, near-field momentum and coherence.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, TransientQuadrupoleSource,
};
use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn medium(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}

fn close(actual: f64, expected: f64, scale: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn real_at(value: C64, omega: f64, time: f64) -> f64 {
    let (s, c) = (omega * time).sin_cos();
    value.re * c + value.im * s
}

fn harmonic_source(q: [C64; 6], omega: f64, time: f64) -> TransientQuadrupoleSource {
    TransientQuadrupoleSource {
        stress_integral_n_m_s: q
            .map(|v| real_at(C64::new(-v.im / omega, v.re / omega), omega, time)),
        stress_moment_n_m: q.map(|v| real_at(v, omega, time)),
        stress_rate_n_m_s: q.map(|v| real_at(C64::new(omega * v.im, -omega * v.re), omega, time)),
        stress_acceleration_n_m_s2: q.map(|v| real_at(v.scale(-omega * omega), omega, time)),
    }
}

fn tensor() -> [f64; 6] {
    [0.7, -0.3, 0.9, 0.23, -0.11, 0.4]
}

// A C2 stress pulse, with its actual causal integral retained after the pulse.
fn pulse(time: f64) -> TransientQuadrupoleSource {
    let duration = 0.012;
    let s = time / duration;
    let [integral, q, rate, acceleration] = if s <= 0.0 {
        [0.0; 4]
    } else if s >= 1.0 {
        [duration / 140.0, 0.0, 0.0, 0.0]
    } else {
        [
            duration
                * (s.powi(4) / 4.0 - 3.0 * s.powi(5) / 5.0 + s.powi(6) / 2.0 - s.powi(7) / 7.0),
            s.powi(3) * (1.0 - s).powi(3),
            (3.0 * s.powi(2) - 12.0 * s.powi(3) + 15.0 * s.powi(4) - 6.0 * s.powi(5)) / duration,
            (6.0 * s - 36.0 * s.powi(2) + 60.0 * s.powi(3) - 30.0 * s.powi(4)) / duration.powi(2),
        ]
    };
    TransientQuadrupoleSource {
        stress_integral_n_m_s: tensor().map(|v| integral * v),
        stress_moment_n_m: tensor().map(|v| q * v),
        stress_rate_n_m_s: tensor().map(|v| rate * v),
        stress_acceleration_n_m_s2: tensor().map(|v| acceleration * v),
    }
}

#[test]
fn g1_transient_array_matches_harmonic_pressure_velocity_and_mean_flux() {
    let source_positions = [[0.0; 3], [-0.013, 0.021, -0.008]];
    let phasors = [
        [
            C64::new(0.7, 0.2),
            C64::new(-0.3, 0.4),
            C64::new(0.9, -0.1),
            C64::new(0.23, -0.3),
            C64::new(-0.11, 0.3),
            C64::new(0.4, 0.2),
        ],
        [
            C64::new(-0.1, 0.5),
            C64::new(0.4, -0.3),
            C64::new(-0.2, 0.7),
            C64::new(-0.13, 0.1),
            C64::new(0.2, -0.4),
            C64::new(0.1, 0.3),
        ],
    ];
    let frequency = 173.0;
    let omega = TAU * frequency;
    for acoustic in [medium(1.23, 317.0), medium(997.0, 1490.0)] {
        for radius in [0.003, 0.2, 12.0] {
            let position = [0.6 * radius, 0.8 * radius, 0.0];
            let transient =
                PreparedTransientQuadrupoleObserver::new(acoustic, position, &source_positions, 2)
                    .unwrap();
            let harmonic = PreparedQuadrupoleObserver::new(
                acoustic,
                frequency,
                position,
                &source_positions,
                2,
            )
            .unwrap()
            .evaluate(&phasors.map(|stress_moment_n_m| QuadrupolePhasor { stress_moment_n_m }))
            .unwrap();
            assert_eq!(transient.acoustic(), acoustic);
            assert_eq!(transient.observer_position_m(), position);
            assert_eq!(transient.source_count(), 2);
            assert_eq!(transient.source_position_m(1), Some(source_positions[1]));
            assert_eq!(transient.source_delay_s(2), None);
            let mut mean_flux = [0.0; 3];
            let mut mean_magnitude = [0.0; 3];
            for sample in 0..64 {
                let time = (0.137 + sample as f64 / 64.0) / frequency;
                let sources = core::array::from_fn::<_, 2, _>(|i| {
                    harmonic_source(
                        phasors[i],
                        omega,
                        time - transient.source_delay_s(i).unwrap(),
                    )
                });
                let field = transient.evaluate_at_retarded(&sources).unwrap();
                close(
                    field.pressure_pa,
                    real_at(harmonic.pressure_pa, omega, time),
                    harmonic.pressure_pa.abs(),
                    2e-12,
                );
                for axis in 0..3 {
                    close(
                        field.particle_velocity_m_s[axis],
                        real_at(harmonic.particle_velocity_m_s[axis], omega, time),
                        harmonic.particle_velocity_m_s[axis].abs(),
                        4e-12,
                    );
                    assert_eq!(
                        field.instantaneous_intensity_w_m2[axis],
                        field.pressure_pa * field.particle_velocity_m_s[axis]
                    );
                    mean_flux[axis] += field.instantaneous_intensity_w_m2[axis] / 64.0;
                    mean_magnitude[axis] += field.instantaneous_intensity_w_m2[axis].abs() / 64.0;
                }
            }
            for axis in 0..3 {
                close(
                    mean_flux[axis],
                    harmonic.active_intensity_w_m2[axis],
                    mean_magnitude[axis],
                    4e-12,
                );
            }
        }
    }
}

#[test]
fn g1_transient_pulse_is_the_retarded_green_hessian() {
    let position = [1.3, -0.8, 0.6];
    let acoustic = medium(1.2, 320.0);
    let observer =
        PreparedTransientQuadrupoleObserver::new(acoustic, position, &[[0.0; 3]], 1).unwrap();
    let time = observer.source_delay_s(0).unwrap() + 0.004;
    let expected = observer
        .evaluate_at_retarded(&[pulse(0.004)])
        .unwrap()
        .pressure_pa;
    // Independently differentiate Qij(t-|x|/c)/(4 pi |x|) in receiver space.
    let potential = |x: [f64; 3], entry: usize| {
        let r = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        pulse(time - r / 320.0).stress_moment_n_m[entry] / (4.0 * PI * r)
    };
    let axes = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
    let errors = [0.0016, 0.0008, 0.0004].map(|h| {
        let mut hessian = 0.0;
        for (entry, (i, j)) in axes.into_iter().enumerate() {
            if i == j {
                let mut plus = position;
                let mut minus = position;
                plus[i] += h;
                minus[i] -= h;
                hessian += (potential(plus, entry) - 2.0 * potential(position, entry)
                    + potential(minus, entry))
                    / (h * h);
            } else {
                for a in [-1.0, 1.0] {
                    for b in [-1.0, 1.0] {
                        let mut x = position;
                        x[i] += a * h;
                        x[j] += b * h;
                        hessian += a * b * potential(x, entry) / (2.0 * h * h);
                    }
                }
            }
        }
        (hessian - expected).abs()
    });
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.4 * errors[1], "{errors:?}");
    assert!(
        errors[2] < 2e-7 * expected.abs(),
        "{errors:?}, expected={expected:e}"
    );
    let undelayed = observer.evaluate_at_retarded(&[pulse(time)]).unwrap();
    assert!((undelayed.pressure_pa - expected).abs() > 0.1 * expected.abs());
}

#[test]
fn g1_transient_momentum_holds_in_all_three_cartesian_directions() {
    let rho = 1.2;
    let acoustic = medium(rho, 320.0);
    let position = [1.3, -0.8, 0.6];
    let at = |x, time| {
        let observer =
            PreparedTransientQuadrupoleObserver::new(acoustic, x, &[[0.0; 3]], 1).unwrap();
        observer
            .evaluate_at_retarded(&[pulse(time - observer.source_delay_s(0).unwrap())])
            .unwrap()
    };
    let time = 0.009;
    let dt = 1e-8;
    let earlier = at(position, time - dt);
    let later = at(position, time + dt);
    let h = 2e-5;
    for axis in 0..3 {
        let mut plus = position;
        let mut minus = position;
        plus[axis] += h;
        minus[axis] -= h;
        let gradient = (at(plus, time).pressure_pa - at(minus, time).pressure_pa) / (2.0 * h);
        let inertia = rho
            * (later.particle_velocity_m_s[axis] - earlier.particle_velocity_m_s[axis])
            / (2.0 * dt);
        close(inertia, -gradient, gradient, 3e-8);
    }
}

#[test]
fn g1_causal_onset_and_original_integral_determine_the_initial_near_field() {
    let observer = PreparedTransientQuadrupoleObserver::new(
        medium(1.2, 320.0),
        [1.3, -0.8, 0.6],
        &[[0.0; 3]],
        1,
    )
    .unwrap();
    let delay = observer.source_delay_s(0).unwrap();
    for time in [0.0, 0.5 * delay, delay] {
        let field = observer
            .evaluate_at_retarded(&[pulse(time - delay)])
            .unwrap();
        assert_eq!(field.pressure_pa, 0.0);
        assert_eq!(field.particle_velocity_m_s, [0.0; 3]);
    }
    let accepted = pulse(0.02);
    let after = observer.evaluate_at_retarded(&[accepted]).unwrap();
    assert_eq!(after.pressure_pa, 0.0);
    assert!(after.particle_velocity_m_s.iter().any(|v| v.abs() > 1e-8));
    assert_eq!(after.instantaneous_intensity_w_m2, [0.0; 3]);
    // Dropping earlier source history changes velocity even after p has passed.
    let reset = TransientQuadrupoleSource {
        stress_integral_n_m_s: [0.0; 6],
        ..accepted
    };
    assert_eq!(
        observer
            .evaluate_at_retarded(&[reset])
            .unwrap()
            .particle_velocity_m_s,
        [0.0; 3]
    );
    assert_eq!(observer.evaluate_at_retarded(&[accepted]).unwrap(), after);
}

#[test]
fn g1_shear_pressure_null_retains_transverse_induction_velocity() {
    let rho = 1.2;
    let c = 320.0;
    let r = 0.7_f64;
    let observer =
        PreparedTransientQuadrupoleObserver::new(medium(rho, c), [r, 0.0, 0.0], &[[0.0; 3]], 1)
            .unwrap();
    let mut source = TransientQuadrupoleSource::default();
    source.stress_integral_n_m_s[3] = 4.0;
    source.stress_moment_n_m[3] = 2.0;
    source.stress_rate_n_m_s[3] = 3.0;
    source.stress_acceleration_n_m_s2[3] = 5.0;
    let field = observer.evaluate_at_retarded(&[source]).unwrap();
    assert_eq!(field.pressure_pa, 0.0);
    assert_eq!(field.particle_velocity_m_s[0], 0.0);
    assert_eq!(field.particle_velocity_m_s[2], 0.0);
    let expected =
        -(24.0 / r.powi(4) + 12.0 / (c * r.powi(3)) + 6.0 / (c * c * r * r)) / (4.0 * PI * rho);
    close(field.particle_velocity_m_s[1], expected, expected, 2e-15);
}

#[test]
fn g3_radiative_directivity_spreading_and_material_scaling_are_physical() {
    let source = TransientQuadrupoleSource {
        stress_acceleration_n_m_s2: [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ..Default::default()
    };
    for (rho, c) in [(1.2, 320.0), (997.0, 1490.0)] {
        for r in [2.0, 4.0] {
            let n = [0.6, 0.8, 0.0];
            let observer = PreparedTransientQuadrupoleObserver::new(
                medium(rho, c),
                n.map(|v| v * r),
                &[[0.0; 3]],
                1,
            )
            .unwrap();
            let field = observer.evaluate_at_retarded(&[source]).unwrap();
            let expected = 2.0 * 0.36 / (4.0 * PI * c * c * r);
            close(field.pressure_pa, expected, expected, 2e-15);
            for axis in 0..3 {
                close(
                    field.particle_velocity_m_s[axis],
                    n[axis] * expected / (rho * c),
                    expected / (rho * c),
                    2e-15,
                );
            }
            close(observer.source_delay_s(0).unwrap(), r / c, r / c, 2e-15);
        }
    }
    let sample = pulse(0.004);
    let at_density = |rho| {
        PreparedTransientQuadrupoleObserver::new(
            medium(rho, 320.0),
            [1.3, -0.8, 0.6],
            &[[0.0; 3]],
            1,
        )
        .unwrap()
        .evaluate_at_retarded(&[sample])
        .unwrap()
    };
    let a = at_density(1.0);
    let b = at_density(1000.0);
    assert_eq!(a.pressure_pa, b.pressure_pa);
    for axis in 0..3 {
        close(
            1000.0 * b.particle_velocity_m_s[axis],
            a.particle_velocity_m_s[axis],
            a.particle_velocity_m_s[axis],
            2e-15,
        );
    }
}

#[test]
fn g3_coherent_sums_keep_cross_terms_and_small_cancellation_residuals() {
    let acoustic = medium(1.2, 320.0);
    let single =
        PreparedTransientQuadrupoleObserver::new(acoustic, [1.0, 0.0, 0.0], &[[0.0; 3]], 1)
            .unwrap();
    let source = |value| TransientQuadrupoleSource {
        stress_moment_n_m: [value, 0.0, 0.0, 0.0, 0.0, 0.0],
        ..Default::default()
    };
    let expected = single.evaluate_at_retarded(&[source(1.0)]).unwrap();
    let pair =
        PreparedTransientQuadrupoleObserver::new(acoustic, [1.0, 0.0, 0.0], &[[0.0; 3]; 2], 2)
            .unwrap();
    let doubled = pair.evaluate_at_retarded(&[source(1.0); 2]).unwrap();
    assert_eq!(doubled.pressure_pa, 2.0 * expected.pressure_pa);
    assert_eq!(
        doubled.instantaneous_intensity_w_m2[0],
        4.0 * expected.instantaneous_intensity_w_m2[0]
    );
    let null = pair
        .evaluate_at_retarded(&[source(1.0), source(-1.0)])
        .unwrap();
    assert_eq!(null.pressure_pa, 0.0);
    assert_eq!(null.particle_velocity_m_s, [0.0; 3]);
    let triple =
        PreparedTransientQuadrupoleObserver::new(acoustic, [1.0, 0.0, 0.0], &[[0.0; 3]; 3], 3)
            .unwrap();
    let residual = triple
        .evaluate_at_retarded(&[source(1e16), source(1.0), source(-1e16)])
        .unwrap();
    close(
        residual.pressure_pa,
        expected.pressure_pa,
        expected.pressure_pa,
        2e-15,
    );
    close(
        residual.particle_velocity_m_s[0],
        expected.particle_velocity_m_s[0],
        expected.particle_velocity_m_s[0],
        2e-15,
    );
}

#[test]
fn g4_admission_and_cancellation_never_publish_partial_fields() {
    let acoustic = medium(1.2, 320.0);
    for position in [[0.0; 3], [f64::NAN, 0.0, 0.0], [1e-200, 0.0, 0.0]] {
        assert!(
            PreparedTransientQuadrupoleObserver::new(acoustic, position, &[[0.0; 3]], 1).is_err()
        );
    }
    assert!(PreparedTransientQuadrupoleObserver::new(acoustic, [1.0; 3], &[[0.0; 3]], 0).is_err());
    assert!(matches!(
        PreparedTransientQuadrupoleObserver::new_cancellable(
            acoustic,
            [1.0; 3],
            &[[0.0; 3]],
            1,
            || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let empty = PreparedTransientQuadrupoleObserver::new(acoustic, [0.0; 3], &[], 0).unwrap();
    assert_eq!(empty.evaluate_at_retarded(&[]).unwrap().pressure_pa, 0.0);
    let observer = PreparedTransientQuadrupoleObserver::new(
        acoustic,
        [0.001, 0.0, 0.0],
        &[[0.0; 3]; 257],
        257,
    )
    .unwrap();
    let accepted = vec![pulse(0.004); 257];
    let baseline = observer.evaluate_at_retarded(&accepted).unwrap();
    assert!(observer.evaluate_at_retarded(&accepted[..256]).is_err());
    let mut polls = 0;
    assert!(matches!(
        observer.evaluate_at_retarded_cancellable(&accepted, || {
            polls += 1;
            polls == 3
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    for field in 0..4 {
        let mut invalid = accepted.clone();
        match field {
            0 => invalid[256].stress_integral_n_m_s[5] = f64::NAN,
            1 => invalid[256].stress_moment_n_m[5] = f64::NAN,
            2 => invalid[256].stress_rate_n_m_s[5] = f64::NAN,
            _ => invalid[256].stress_acceleration_n_m_s2[5] = f64::NAN,
        }
        assert!(observer.evaluate_at_retarded(&invalid).is_err());
    }
    let mut overflowing = accepted.clone();
    overflowing[256].stress_moment_n_m[0] = f64::MAX;
    assert!(observer.evaluate_at_retarded(&overflowing).is_err());
    assert_eq!(observer.evaluate_at_retarded(&accepted).unwrap(), baseline);
}

#[test]
fn g1_history_callback_preserves_distinct_causal_arrivals_and_refusals() {
    let acoustic = medium(1.2, 320.0);
    let positions = [[0.1, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let observer =
        PreparedTransientQuadrupoleObserver::new(acoustic, [0.0; 3], &positions, 2).unwrap();
    let first =
        PreparedTransientQuadrupoleObserver::new(acoustic, [0.0; 3], &positions[..1], 1).unwrap();
    let time = 0.002;
    let mut callbacks = Vec::new();
    let actual = observer
        .evaluate_at_time(time, |index, emission_time| {
            callbacks.push((index, emission_time));
            Ok(pulse(emission_time))
        })
        .unwrap();
    assert_eq!(callbacks.len(), 2);
    assert_eq!(callbacks[0].0, 0);
    assert_eq!(callbacks[1].0, 1);
    assert!(callbacks[0].1 > 0.0 && callbacks[1].1 < 0.0);
    let expected = first
        .evaluate_at_retarded(&[pulse(time - 0.1 / 320.0)])
        .unwrap();
    assert_eq!(actual, expected);
    let later = 0.011;
    let explicit: Vec<_> = (0..2)
        .map(|i| pulse(later - observer.source_delay_s(i).unwrap()))
        .collect();
    assert_eq!(
        observer
            .evaluate_at_time(later, |_, t| Ok(pulse(t)))
            .unwrap(),
        observer.evaluate_at_retarded(&explicit).unwrap()
    );
    for bad_time in [f64::NAN, 1e30] {
        let mut invoked = false;
        assert!(
            observer
                .evaluate_at_time(bad_time, |_, _| {
                    invoked = true;
                    Ok(TransientQuadrupoleSource::default())
                })
                .is_err()
        );
        assert!(!invoked);
    }
    let source_error = AeroacError::InvalidParameter {
        what: "physical history has no sample",
    };
    assert_eq!(
        observer.evaluate_at_time(time, |_, _| Err(source_error.clone())),
        Err(source_error)
    );
    let mut invoked = false;
    assert!(matches!(
        observer.evaluate_at_time_cancellable(
            time,
            |_, _| {
                invoked = true;
                Ok(TransientQuadrupoleSource::default())
            },
            || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(!invoked);
}
