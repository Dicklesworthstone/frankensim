//! Original rotor period maps: parametric resonance, gyroscopics and coordinates.
use fs_math::c64::C64;
use fs_mbd::rotordynamics::periodic::{
    MAX_ROTOR_PERIODIC_SEGMENTS, RotorFloquetBudget, RotorPeriodicSegment, periodic_floquet,
};
use fs_mbd::rotordynamics::{
    MAX_ROTOR_DOFS, ROTOR_QR_SWEEPS_PER_WINDOW, RotorBearingModel, RotorError,
};
use std::cell::Cell;
use std::f64::consts::PI;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn spectrum_near(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let mut used = vec![false; actual.len()];
    for target in expected {
        let (index, error) = actual
            .iter()
            .enumerate()
            .filter(|(index, _)| !used[*index])
            .map(|(index, value)| (index, (value.re - target.re).hypot(value.im - target.im)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        assert!(error <= tolerance * (1.0 + target.re.hypot(target.im)));
        used[index] = true;
    }
}

fn oscillator(stiffness: f64, damping: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[1.0], &[damping], &[stiffness], &[0.0], || false).unwrap()
}

#[test]
fn quarter_period_stiffness_cycle_is_unstable_despite_stable_frozen_poles() {
    let soft = oscillator(1.0, 0.0);
    let stiff = oscillator(4.0, 0.0);
    for (model, frequency) in [(&soft, 1.0), (&stiff, 2.0)] {
        let poles = model
            .at_speed(7.0, || false)
            .unwrap()
            .poles(|| false)
            .unwrap();
        spectrum_near(
            &poles.poles_per_s,
            &[C64::new(0.0, frequency), C64::new(0.0, -frequency)],
            1e-14,
        );
    }
    let cycle = [
        RotorPeriodicSegment {
            model: &soft,
            duration_seconds: PI / 2.0,
        },
        RotorPeriodicSegment {
            model: &stiff,
            duration_seconds: PI / 4.0,
        },
    ];
    let result =
        periodic_floquet(7.0, &cycle, 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    // A quarter-period oscillator map is [[0,1/w],[-w,0]]. Chronological
    // multiplication E_stiff E_soft gives diag(-1/2,-2), not its reversal.
    for (actual, expected) in result.monodromy.iter().zip([-0.5, 0.0, 0.0, -2.0]) {
        near(*actual, expected, 4e-13);
    }
    spectrum_near(
        &result.multipliers,
        &[C64::new(-0.5, 0.0), C64::new(-2.0, 0.0)],
        4e-13,
    );
    near(result.period_seconds, 3.0 * PI / 4.0, 1e-15);
    near(result.spectral_radius, 2.0, 4e-13);
    near(
        result.growth_rate_per_s,
        2.0_f64.ln() / result.period_seconds,
        4e-13,
    );
    assert!(result.maximum_mass_residual <= 1e-12);
    assert_eq!(
        result,
        periodic_floquet(7.0, &cycle, 1e-12, RotorFloquetBudget::default(), || false).unwrap()
    );
}

#[test]
fn constant_full_inertia_gyro_cycle_matches_poles_and_segment_refinement() {
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.5, 0.5, 1.5],
        &[0.25, 0.0625, 0.0625, 0.1875],
        &[8.0, 2.0, 2.0, 6.0],
        &[0.0, 0.75, -0.75, 0.0],
        || false,
    )
    .unwrap();
    let duration = 0.75;
    for spin in [-3.0, 3.0] {
        let single = [RotorPeriodicSegment {
            model: &model,
            duration_seconds: duration,
        }];
        let result = periodic_floquet(spin, &single, 1e-12, RotorFloquetBudget::default(), || {
            false
        })
        .unwrap();
        let poles = model
            .at_speed(spin, || false)
            .unwrap()
            .poles(|| false)
            .unwrap();
        let expected: Vec<_> = poles
            .poles_per_s
            .iter()
            .map(|pole| {
                let amplitude = (pole.re * duration).exp();
                let (sin, cos) = (pole.im * duration).sin_cos();
                C64::new(amplitude * cos, amplitude * sin)
            })
            .collect();
        spectrum_near(&result.multipliers, &expected, 2e-11);
        let refined = [RotorPeriodicSegment {
            model: &model,
            duration_seconds: duration / 4.0,
        }; 4];
        let split = periodic_floquet(spin, &refined, 1e-12, RotorFloquetBudget::default(), || {
            false
        })
        .unwrap();
        for (actual, expected) in split.monodromy.iter().zip(&result.monodromy) {
            near(*actual, *expected, 3e-12);
        }
        // Liouville: det(P)=exp(-T tr(M^-1 C)); skew gyro contributes zero.
        let determinant = result
            .multipliers
            .iter()
            .fold(C64::new(1.0, 0.0), |a, b| a * *b);
        near(determinant.re, (-0.25 * duration).exp(), 3e-12);
        near(determinant.im, 0.0, 3e-12);
        assert!(result.spectral_radius < 1.0);
    }
}

#[test]
fn signed_gyroscopic_cycle_retains_whirl_handedness_and_free_displacements() {
    let model = RotorBearingModel::new(
        2,
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &[0.0; 4],
        &[0.0, 0.5, -0.5, 0.0],
        || false,
    )
    .unwrap();
    let duration = 0.5;
    for spin in [-2.0, 2.0] {
        let result = periodic_floquet(
            spin,
            &[RotorPeriodicSegment {
                model: &model,
                duration_seconds: duration,
            }],
            1e-12,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        // v'=[-a*v_y, a*v_x]. Starting at q=0 and v=[1,0] has a
        // signed circular velocity and the integrated original displacement.
        let rate = 0.5 * spin;
        let (sin, cos) = (rate * duration).sin_cos();
        for (actual, expected) in [
            result.monodromy[2],
            result.monodromy[6],
            result.monodromy[10],
            result.monodromy[14],
        ]
        .into_iter()
        .zip([sin / rate, (1.0 - cos) / rate, cos, sin])
        {
            near(actual, expected, 4e-13);
        }
        near(
            result.monodromy[10].powi(2) + result.monodromy[14].powi(2),
            1.0,
            4e-13,
        );
        spectrum_near(
            &result.multipliers,
            &[
                C64::new(1.0, 0.0),
                C64::new(1.0, 0.0),
                C64::new(cos, sin),
                C64::new(cos, -sin),
            ],
            4e-13,
        );
        near(result.spectral_radius, 1.0, 4e-13);
    }
}

fn congruence(matrix: &[f64], transform: &[f64; 4]) -> Vec<f64> {
    let mut result = vec![0.0; 4];
    for i in 0..2 {
        for j in 0..2 {
            for p in 0..2 {
                for q in 0..2 {
                    result[2 * i + j] +=
                        transform[2 * p + i] * matrix[2 * p + q] * transform[2 * q + j];
                }
            }
        }
    }
    result
}

fn state_map(matrix: &[f64], state: &[f64; 4]) -> [f64; 4] {
    std::array::from_fn(|i| (0..4).map(|j| matrix[4 * i + j] * state[j]).sum())
}

fn change_coordinates(transform: &[f64; 4], state: &[f64; 4]) -> [f64; 4] {
    std::array::from_fn(|i| {
        let offset = 2 * (i / 2);
        transform[2 * (i % 2)] * state[offset] + transform[2 * (i % 2) + 1] * state[offset + 1]
    })
}

#[test]
fn active_cross_coupled_cycle_retains_full_coordinate_maps_and_phase_origin() {
    let mass = [1.0, 0.0, 0.0, 1.0];
    let gyro = [0.0, 0.5, -0.5, 0.0];
    let damping = [[0.125, 0.25, -0.125, -0.0625], [0.25, 0.125, 0.0625, 0.125]];
    let stiffness = [[4.0, 0.5, -0.25, 2.0], [3.0, -0.25, 0.75, 5.0]];
    let original: Vec<_> = (0..2)
        .map(|i| {
            RotorBearingModel::new(2, &mass, &damping[i], &stiffness[i], &gyro, || false).unwrap()
        })
        .collect();
    let transform = [2.0, 0.25, -0.5, 1.5];
    let transformed: Vec<_> = (0..2)
        .map(|i| {
            RotorBearingModel::new(
                2,
                &congruence(&mass, &transform),
                &congruence(&damping[i], &transform),
                &congruence(&stiffness[i], &transform),
                &congruence(&gyro, &transform),
                || false,
            )
            .unwrap()
        })
        .collect();
    let analyze = |models: &[RotorBearingModel]| {
        periodic_floquet(
            -1.25,
            &[
                RotorPeriodicSegment {
                    model: &models[0],
                    duration_seconds: 0.3,
                },
                RotorPeriodicSegment {
                    model: &models[1],
                    duration_seconds: 0.2,
                },
            ],
            1e-12,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap()
    };
    let a = analyze(&original);
    let b = analyze(&transformed);
    spectrum_near(&a.multipliers, &b.multipliers, 2e-11);
    let initial = [0.3, -0.5, 0.1, 0.7];
    let actual = state_map(&a.monodromy, &change_coordinates(&transform, &initial));
    let expected = change_coordinates(&transform, &state_map(&b.monodromy, &initial));
    for (actual, expected) in actual.into_iter().zip(expected) {
        near(actual, expected, 3e-12);
    }
    let shifted = periodic_floquet(
        -1.25,
        &[
            RotorPeriodicSegment {
                model: &original[1],
                duration_seconds: 0.2,
            },
            RotorPeriodicSegment {
                model: &original[0],
                duration_seconds: 0.3,
            },
        ],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    spectrum_near(&a.multipliers, &shifted.multipliers, 2e-11);
    assert!(
        a.monodromy
            .iter()
            .zip(&shifted.monodromy)
            .any(|(x, y)| (x - y).abs() > 1e-3)
    );
}

#[test]
fn free_coordinates_survive_at_the_existing_rotor_dimension_and_cycle_caps() {
    let n = MAX_ROTOR_DOFS;
    let mut mass = vec![0.0; n * n];
    for i in 0..n {
        mass[i * n + i] = 1.0;
    }
    let zero = vec![0.0; n * n];
    let model = RotorBearingModel::new(n, &mass, &zero, &zero, &zero, || false).unwrap();
    let cycle = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.125,
    }];
    let response =
        periodic_floquet(-2.0, &cycle, 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    let d = 2 * n;
    assert_eq!(response.monodromy.len(), d * d);
    assert_eq!(response.multipliers, vec![C64::new(1.0, 0.0); d]);
    for i in 0..d {
        for j in 0..d {
            let expected = if i == j {
                1.0
            } else if i < n && j == n + i {
                0.125
            } else {
                0.0
            };
            near(response.monodromy[i * d + j], expected, 1e-15);
        }
    }
    assert_eq!(response.spectral_radius, 1.0);
    near(response.growth_rate_per_s, 0.0, 1e-15);

    let scalar = oscillator(0.0, 0.0);
    let segment = RotorPeriodicSegment {
        model: &scalar,
        duration_seconds: 1.0,
    };
    let segments = vec![segment; MAX_ROTOR_PERIODIC_SEGMENTS];
    let budget = RotorFloquetBudget {
        max_segments: MAX_ROTOR_PERIODIC_SEGMENTS,
        ..RotorFloquetBudget::default()
    };
    let longest = periodic_floquet(0.0, &segments, 1e-12, budget, || false).unwrap();
    near(
        longest.monodromy[1],
        MAX_ROTOR_PERIODIC_SEGMENTS as f64,
        1e-15,
    );
    assert!(
        periodic_floquet(
            0.0,
            &vec![segment; MAX_ROTOR_PERIODIC_SEGMENTS + 1],
            1e-12,
            budget,
            || false,
        )
        .is_err()
    );
    assert!(RotorBearingModel::new(n + 1, &[], &[], &[], &[], || false,).is_err());
}

#[test]
fn mismatched_physics_malformed_clock_and_exhausted_budgets_refuse() {
    let base = oscillator(1.0, 0.0);
    let segment = RotorPeriodicSegment {
        model: &base,
        duration_seconds: 0.5,
    };
    let defaults = RotorFloquetBudget::default();
    assert!(periodic_floquet(0.0, &[], 1e-12, defaults, || false).is_err());
    for spin in [f64::NAN, f64::INFINITY] {
        assert!(periodic_floquet(spin, &[segment], 1e-12, defaults, || false).is_err());
    }
    for duration in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(
            periodic_floquet(
                0.0,
                &[RotorPeriodicSegment {
                    model: &base,
                    duration_seconds: duration
                },],
                1e-12,
                defaults,
                || false
            )
            .is_err()
        );
    }
    for tolerance in [0.0, 1.0, f64::NAN] {
        assert!(periodic_floquet(0.0, &[segment], tolerance, defaults, || false).is_err());
    }
    let different_mass = RotorBearingModel::new(
        1,
        &[f64::from_bits(1.0_f64.to_bits() + 1)],
        &[0.0],
        &[1.0],
        &[0.0],
        || false,
    )
    .unwrap();
    assert!(
        periodic_floquet(
            0.0,
            &[
                segment,
                RotorPeriodicSegment {
                    model: &different_mass,
                    duration_seconds: 0.5
                },
            ],
            1e-12,
            defaults,
            || false
        )
        .is_err()
    );
    let identity = [1.0, 0.0, 0.0, 1.0];
    let a =
        RotorBearingModel::new(2, &identity, &[0.0; 4], &identity, &[0.0; 4], || false).unwrap();
    let b = RotorBearingModel::new(
        2,
        &identity,
        &[0.0; 4],
        &identity,
        &[0.0, 0.5, -0.5, 0.0],
        || false,
    )
    .unwrap();
    assert!(
        periodic_floquet(
            0.0,
            &[
                RotorPeriodicSegment {
                    model: &a,
                    duration_seconds: 0.5
                },
                RotorPeriodicSegment {
                    model: &b,
                    duration_seconds: 0.5
                },
            ],
            1e-12,
            defaults,
            || false
        )
        .is_err()
    );
    assert!(
        periodic_floquet(
            0.0,
            &[
                segment,
                RotorPeriodicSegment {
                    model: &a,
                    duration_seconds: 0.5
                },
            ],
            1e-12,
            defaults,
            || false
        )
        .is_err()
    );
    for budget in [
        RotorFloquetBudget {
            max_segments: 0,
            ..defaults
        },
        RotorFloquetBudget {
            max_segments: MAX_ROTOR_PERIODIC_SEGMENTS + 1,
            ..defaults
        },
        RotorFloquetBudget {
            max_qr_sweeps_per_window: 0,
            ..defaults
        },
        RotorFloquetBudget {
            max_qr_sweeps_per_window: ROTOR_QR_SWEEPS_PER_WINDOW + 1,
            ..defaults
        },
        RotorFloquetBudget {
            max_dense_bytes: 0,
            ..defaults
        },
        RotorFloquetBudget {
            max_multiply_terms: 0,
            ..defaults
        },
    ] {
        assert!(matches!(
            periodic_floquet(0.0, &[segment], 1e-12, budget, || false),
            Err(RotorError::PeriodicBudget(_))
        ));
    }
    assert!(matches!(
        periodic_floquet(
            0.0,
            &[segment; 2],
            1e-12,
            RotorFloquetBudget {
                max_segments: 1,
                ..defaults
            },
            || false
        ),
        Err(RotorError::PeriodicBudget(_))
    ));
}

#[test]
fn original_equation_residual_and_unrepresentable_growth_refuse() {
    let coupled = RotorBearingModel::new(
        2,
        &[2.0, 0.5, 0.5, 1.0],
        &[0.1, 0.3, -0.2, 0.4],
        &[3.0, 0.2, -0.1, 2.0],
        &[0.0, 0.5, -0.5, 0.0],
        || false,
    )
    .unwrap();
    assert!(matches!(
        periodic_floquet(
            0.75,
            &[RotorPeriodicSegment {
                model: &coupled,
                duration_seconds: 0.25
            },],
            1e-30,
            RotorFloquetBudget::default(),
            || false
        ),
        Err(RotorError::Residual { .. })
    ));
    // Both exact multipliers are exp(-1000), below representable binary64.
    // The returned radius must never be zero with a manufactured log(0).
    let vanished = oscillator(1_000_000.0, 2_000.0);
    assert!(matches!(
        periodic_floquet(
            0.0,
            &[RotorPeriodicSegment {
                model: &vanished,
                duration_seconds: 1.0
            },],
            1e-12,
            RotorFloquetBudget::default(),
            || false
        ),
        Err(RotorError::Invalid(_))
    ));
    // The slow oscillator retains unit-radius multipliers; the fast pair
    // still cannot be reported as zeros hidden by that finite spectral radius.
    let mixed = RotorBearingModel::new(
        2,
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0, 0.0, 0.0, 2_000.0],
        &[1.0, 0.0, 0.0, 1_000_000.0],
        &[0.0; 4],
        || false,
    )
    .unwrap();
    assert!(matches!(
        periodic_floquet(
            0.0,
            &[RotorPeriodicSegment {
                model: &mixed,
                duration_seconds: 1.0
            },],
            1e-12,
            RotorFloquetBudget::default(),
            || false
        ),
        Err(RotorError::Invalid(_))
    ));
    let explosive = oscillator(-1_000_000.0, 0.0);
    assert!(matches!(
        periodic_floquet(
            0.0,
            &[RotorPeriodicSegment {
                model: &explosive,
                duration_seconds: 1.0
            },],
            1e-12,
            RotorFloquetBudget::default(),
            || false
        ),
        Err(RotorError::PeriodicExponential(_))
    ));
    let finite = oscillator(1.0, 0.0);
    assert!(
        periodic_floquet(
            0.0,
            &[RotorPeriodicSegment {
                model: &finite,
                duration_seconds: f64::MIN_POSITIVE
            },],
            1e-12,
            RotorFloquetBudget::default(),
            || false
        )
        .is_err()
    );
}

#[test]
fn cancellation_during_assembly_propagation_and_final_publication_is_atomic() {
    let a = oscillator(1.0, 0.125);
    let b = oscillator(4.0, -0.0625);
    let originals = [a.clone(), b.clone()];
    let cycle = [
        RotorPeriodicSegment {
            model: &a,
            duration_seconds: 0.25,
        },
        RotorPeriodicSegment {
            model: &b,
            duration_seconds: 0.125,
        },
    ];
    let calls = Cell::new(0);
    periodic_floquet(1.0, &cycle, 1e-12, RotorFloquetBudget::default(), || {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let completed = calls.get();
    assert!(completed > 20);
    for stop in [1, 3, completed / 2, completed - 1, completed] {
        let polls = Cell::new(0);
        assert_eq!(
            periodic_floquet(1.0, &cycle, 1e-12, RotorFloquetBudget::default(), || {
                polls.set(polls.get() + 1);
                polls.get() >= stop
            }),
            Err(RotorError::Cancelled)
        );
    }
    assert_eq!([a, b], originals);
}
