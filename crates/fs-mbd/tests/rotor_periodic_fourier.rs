//! G1/G3: original-coordinate Fourier forcing, resonant growth and switched physics.
use fs_mbd::rotordynamics::periodic::{
    MAX_ROTOR_FORCE_TONES, RotorFloquetBudget, RotorFourierForce, RotorHarmonicForce,
    RotorPeriodicSegment, periodic_floquet, prepare_forced_cycle, prepare_fourier_cycle,
};
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::f64::consts::TAU;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn oscillator(m: f64, c: f64, k: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[m], &[c], &[k], &[0.0], || false).unwrap()
}

// Independent scalar dynamic-stiffness algebra, with F_hat=Fc-i*Fs.
fn particular(
    m: f64,
    c: f64,
    k: f64,
    dc: f64,
    tones: &[(usize, f64, f64)],
    fundamental: f64,
    time: f64,
) -> (f64, f64) {
    let (mut q, mut v) = (dc / k, 0.0);
    for &(order, cosine, sine) in tones {
        let w = order as f64 * fundamental;
        let (a, b) = (k - m * w * w, c * w);
        let real = (a * cosine - b * sine) / (a * a + b * b);
        let imag = -(b * cosine + a * sine) / (a * a + b * b);
        let (sin, cos) = (w * time).sin_cos();
        q += real * cos - imag * sin;
        v -= w * (real * sin + imag * cos);
    }
    (q, v)
}

#[test]
fn steady_and_multiple_harmonics_recover_the_complete_resonant_orbit() {
    let model = oscillator(1.0, 0.04, 1.0);
    let fractions = [0.11, 0.29, 0.60];
    let segments: Vec<_> = fractions
        .iter()
        .map(|fraction| RotorPeriodicSegment {
            model: &model,
            duration_seconds: fraction * TAU,
        })
        .collect();
    let unscaled = [(1, 0.4, -0.2), (3, -0.7, 0.3), (7, 0.2, 0.1)];
    for amplitude in [1e-80, 1.0, 1e80] {
        let steady = [0.9 * amplitude];
        let cosine = [[0.4 * amplitude], [-0.7 * amplitude], [0.2 * amplitude]];
        let sine = [[-0.2 * amplitude], [0.3 * amplitude], [0.1 * amplitude]];
        let tones: Vec<_> = (0..3)
            .map(|i| RotorHarmonicForce {
                cosine: &cosine[i],
                sine: &sine[i],
            })
            .collect();
        let forces = vec![
            RotorFourierForce {
                steady: &steady,
                harmonics: &tones
            };
            3
        ];
        let plan = prepare_fourier_cycle(
            1.0,
            &segments,
            &forces,
            &[1, 3, 7],
            1e-12,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        let orbit = plan.periodic_orbit(1e-9, 1e8, || false).unwrap();
        for (&t, state) in orbit
            .trajectory
            .boundary_times_seconds
            .iter()
            .zip(&orbit.trajectory.boundary_states)
        {
            let (q, v) = particular(1.0, 0.04, 1.0, 0.9, &unscaled, 1.0, t);
            near(state.displacement[0] / amplitude, q, 3e-9);
            near(state.velocity[0] / amplitude, v, 3e-9);
        }
        assert_eq!(plan.floquet_report().multipliers.len(), 2);
        near(
            plan.floquet_report().spectral_radius,
            (-0.02 * TAU).exp(),
            2e-11,
        );
        assert!(orbit.relative_fixed_point_residual < 1e-9);
        assert!(orbit.relative_closure_residual < 1e-9);
    }
}

#[test]
fn full_mutual_mass_integrates_steady_and_rotating_loads_in_original_coordinates() {
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.5, 0.5, 1.0],
        &[0.0; 4],
        &[0.0; 4],
        &[0.0; 4],
        || false,
    )
    .unwrap();
    let fundamental = 1.7;
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: TAU / (4.0 * fundamental),
    }; 4];
    let steady = [0.3, -0.2];
    let cosine = [[1.3, 0.4], [-0.5, 0.7]];
    let sine = [[0.2, -0.7], [0.8, -0.1]];
    let tones: Vec<_> = (0..2)
        .map(|i| RotorHarmonicForce {
            cosine: &cosine[i],
            sine: &sine[i],
        })
        .collect();
    let forces = [RotorFourierForce {
        steady: &steady,
        harmonics: &tones,
    }; 4];
    let orders = [1, 4];
    let plan = prepare_fourier_cycle(
        0.0,
        &segments,
        &forces,
        &orders,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let initial = RotorTransientState {
        displacement: vec![0.2, -0.4],
        velocity: vec![-0.1, 0.3],
    };
    let motion = plan.advance(&initial, 1e-10, || false).unwrap();
    let acceleration = |f: &[f64]| {
        [
            (f[0] - 0.5 * f[1]) / 1.75,
            (-0.5 * f[0] + 2.0 * f[1]) / 1.75,
        ]
    };
    let dc = acceleration(&steady);
    for (&t, state) in motion
        .boundary_times_seconds
        .iter()
        .zip(&motion.boundary_states)
    {
        for i in 0..2 {
            let mut q = initial.displacement[i] + t * initial.velocity[i] + 0.5 * t * t * dc[i];
            let mut v = initial.velocity[i] + t * dc[i];
            for j in 0..2 {
                let w = orders[j] as f64 * fundamental;
                let (sin, cos) = (w * t).sin_cos();
                let ac = acceleration(&cosine[j])[i];
                let az = acceleration(&sine[j])[i];
                q += ac * (1.0 - cos) / (w * w) + az * (t / w - sin / (w * w));
                v += ac * sin / w + az * (1.0 - cos) / w;
            }
            near(state.displacement[i], q, 3e-11);
            near(state.velocity[i], v, 3e-11);
        }
    }
    assert_eq!(initial.displacement, vec![0.2, -0.4]);
    assert_eq!(plan.floquet_report().multipliers.len(), 4);
    assert!(plan.periodic_orbit(1e-9, 1e8, || false).is_err());
}

#[test]
fn coefficient_and_amplitude_switches_keep_one_common_phase_for_every_tone() {
    let first = oscillator(1.0, 0.2, 3.0);
    let second = oscillator(1.0, 0.5, 7.0);
    let segments = [
        RotorPeriodicSegment {
            model: &first,
            duration_seconds: 0.17 * TAU,
        },
        RotorPeriodicSegment {
            model: &second,
            duration_seconds: 0.28 * TAU,
        },
        RotorPeriodicSegment {
            model: &first,
            duration_seconds: 0.55 * TAU,
        },
    ];
    let dc = [[0.2], [-0.4], [0.7]];
    let cosine = [[[0.4], [-0.8]], [[-0.3], [0.6]], [[0.1], [0.9]]];
    let sine = [[[0.3], [0.1]], [[0.2], [-0.4]], [[-0.2], [0.5]]];
    let tones: Vec<Vec<_>> = (0..3)
        .map(|i| {
            (0..2)
                .map(|j| RotorHarmonicForce {
                    cosine: &cosine[i][j],
                    sine: &sine[i][j],
                })
                .collect()
        })
        .collect();
    let forces: Vec<_> = (0..3)
        .map(|i| RotorFourierForce {
            steady: &dc[i],
            harmonics: &tones[i],
        })
        .collect();
    let plan = prepare_fourier_cycle(
        0.0,
        &segments,
        &forces,
        &[1, 3],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let initial = RotorTransientState {
        displacement: vec![0.2],
        velocity: vec![-0.3],
    };
    let motion = plan.advance(&initial, 1e-10, || false).unwrap();
    let (mut time, mut q, mut v) = (0.0, 0.2, -0.3);
    for i in 0..3 {
        let (c, k) = if i == 1 { (0.5, 7.0) } else { (0.2, 3.0) };
        let h = segments[i].duration_seconds;
        let input = [
            (1, cosine[i][0][0], sine[i][0][0]),
            (3, cosine[i][1][0], sine[i][1][0]),
        ];
        let begin = particular(1.0, c, k, dc[i][0], &input, 1.0, time);
        let end = particular(1.0, c, k, dc[i][0], &input, 1.0, time + h);
        let rate = 0.5 * c;
        let w = f64::sqrt(k - rate * rate);
        let a = q - begin.0;
        let b = (v - begin.1 + rate * a) / w;
        let (sin, cos) = (w * h).sin_cos();
        let decay = (-rate * h).exp();
        q = end.0 + decay * (a * cos + b * sin);
        v = end.1 + decay * ((w * b - rate * a) * cos - (w * a + rate * b) * sin);
        time += h;
        near(motion.boundary_states[i + 1].displacement[0], q, 3e-10);
        near(motion.boundary_states[i + 1].velocity[0], v, 3e-10);
    }
    let homogeneous =
        periodic_floquet(0.0, &segments, 1e-12, RotorFloquetBudget::default(), || {
            false
        })
        .unwrap();
    for (a, b) in homogeneous
        .monodromy
        .iter()
        .zip(&plan.floquet_report().monodromy)
    {
        near(*a, *b, 2e-12);
    }
}

#[test]
fn lossless_resonant_tone_grows_with_the_steady_and_off_resonant_response() {
    let model = oscillator(1.0, 0.0, 1.0);
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: TAU / 4.0,
    }; 4];
    let tones = [
        RotorHarmonicForce {
            cosine: &[1.0],
            sine: &[0.0],
        },
        RotorHarmonicForce {
            cosine: &[0.6],
            sine: &[0.0],
        },
    ];
    let forces = [RotorFourierForce {
        steady: &[0.3],
        harmonics: &tones,
    }; 4];
    let plan = prepare_fourier_cycle(
        0.0,
        &segments,
        &forces,
        &[1, 3],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let mut initial = RotorTransientState {
        displacement: vec![0.0],
        velocity: vec![0.0],
    };
    for cycle in 0..2 {
        let motion = plan.advance(&initial, 1e-9, || false).unwrap();
        for (&local, state) in motion
            .boundary_times_seconds
            .iter()
            .zip(&motion.boundary_states)
        {
            let t = local + cycle as f64 * TAU;
            let q = 0.3 * (1.0 - t.cos())
                + 0.5 * t * t.sin()
                + 0.6 / (1.0 - 9.0) * ((3.0 * t).cos() - t.cos());
            let v = 0.3 * t.sin()
                + 0.5 * (t.sin() + t * t.cos())
                + 0.6 / (1.0 - 9.0) * (-3.0 * (3.0 * t).sin() + t.sin());
            near(state.displacement[0], q, 2e-10);
            near(state.velocity[0], v, 2e-10);
        }
        initial = motion.boundary_states.last().unwrap().clone();
    }
    assert!(plan.periodic_orbit(1e-9, 1e8, || false).is_err());
}

#[test]
fn empty_and_silent_spectra_reuse_the_held_source_and_do_not_add_modes() {
    let model = oscillator(1.0, 0.2, 2.0);
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.4,
    }; 2];
    let held = prepare_forced_cycle(
        0.0,
        &segments,
        &[&[0.3], &[-0.2]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let orders: Vec<_> = (1..MAX_ROTOR_FORCE_TONES).chain([1024]).collect();
    let silent = vec![
        RotorHarmonicForce {
            cosine: &[0.0],
            sine: &[0.0]
        };
        MAX_ROTOR_FORCE_TONES
    ];
    for tones in [&[][..], &silent[..]] {
        let rows = [
            RotorFourierForce {
                steady: &[0.3],
                harmonics: tones,
            },
            RotorFourierForce {
                steady: &[-0.2],
                harmonics: tones,
            },
        ];
        let plan = prepare_fourier_cycle(
            0.0,
            &segments,
            &rows,
            &orders[..tones.len()],
            1e-12,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        assert_eq!(plan.floquet_report(), held.floquet_report());
        assert_eq!(plan.forcing_offset(), held.forcing_offset());
        let initial = RotorTransientState {
            displacement: vec![0.1],
            velocity: vec![0.2],
        };
        assert_eq!(
            plan.advance(&initial, 1e-10, || false).unwrap(),
            held.advance(&initial, 1e-10, || false).unwrap()
        );
    }
}

#[test]
fn complete_spectrum_budget_and_shape_admission_precede_publication() {
    let model = oscillator(1.0, 0.1, 2.0);
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.7,
    }];
    let tone = RotorHarmonicForce {
        cosine: &[0.4],
        sine: &[0.2],
    };
    let tones = [tone; 3];
    let rows = [RotorFourierForce {
        steady: &[0.1],
        harmonics: &tones,
    }];
    let default = RotorFloquetBudget::default();
    let limited = RotorFloquetBudget {
        max_multiply_terms: 70_000,
        ..default
    };
    // n=1: one held exponential and QR fit, four full source flows do not.
    prepare_forced_cycle(0.0, &segments, &[&[0.1]], 1e-12, limited, || false).unwrap();
    assert!(matches!(
        prepare_fourier_cycle(0.0, &segments, &rows, &[1, 2, 3], 1e-12, limited, || false),
        Err(RotorError::PeriodicBudget(_))
    ));
    assert!(matches!(
        prepare_fourier_cycle(
            0.0,
            &segments,
            &rows,
            &[1, 2, 3],
            1e-12,
            RotorFloquetBudget {
                max_dense_bytes: 1,
                ..default
            },
            || false
        ),
        Err(RotorError::PeriodicBudget(_))
    ));
    for orders in [
        &[0, 1, 2][..],
        &[1, 1, 2],
        &[2, 1, 3],
        &[1, 2, 1025],
        &[1, 2],
    ] {
        assert!(
            prepare_fourier_cycle(0.0, &segments, &rows, orders, 1e-12, default, || false).is_err()
        );
    }
    let excessive = vec![tone; MAX_ROTOR_FORCE_TONES + 1];
    let orders: Vec<_> = (1..=excessive.len()).collect();
    assert!(
        prepare_fourier_cycle(
            0.0,
            &segments,
            &[RotorFourierForce {
                steady: &[0.0],
                harmonics: &excessive
            }],
            &orders,
            1e-12,
            default,
            || false
        )
        .is_err()
    );
    let malformed = [RotorFourierForce {
        steady: &[f64::NAN],
        harmonics: &tones,
    }];
    assert!(
        prepare_fourier_cycle(
            0.0,
            &segments,
            &malformed,
            &[1, 2, 3],
            1e-12,
            default,
            || false
        )
        .is_err()
    );
    let wrong = [RotorFourierForce {
        steady: &[0.0, 0.0],
        harmonics: &tones,
    }];
    assert!(
        prepare_fourier_cycle(0.0, &segments, &wrong, &[1, 2, 3], 1e-12, default, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        prepare_fourier_cycle(0.0, &segments, &rows, &[1, 2, 3], 1e-12, default, || {
            polls += 1;
            polls == 28
        }),
        Err(RotorError::Cancelled)
    ));
    assert_eq!(model.mass(), &[1.0]);
    assert_eq!(model.stiffness(), &[2.0]);
}
