//! Quadratic waveforms retain physical p/u interference and authored statistics.
use fs_aeroac::AeroacError;
use fs_aeroac::observation::*;
use fs_fft::harmonic_window::{HarmonicProductTone, HarmonicWindow, HarmonicWindowProfile};
use fs_math::c64::C64;
use std::f64::consts::TAU;

type Line = HarmonicFieldLine<2>;

const PRODUCTS: [(AcousticFieldComponent, AcousticFieldComponent); 6] = [
    (
        AcousticFieldComponent::Pressure,
        AcousticFieldComponent::Pressure,
    ),
    (
        AcousticFieldComponent::Pressure,
        AcousticFieldComponent::Velocity { axis: 0 },
    ),
    (
        AcousticFieldComponent::Pressure,
        AcousticFieldComponent::Velocity { axis: 1 },
    ),
    (
        AcousticFieldComponent::Velocity { axis: 0 },
        AcousticFieldComponent::Velocity { axis: 0 },
    ),
    (
        AcousticFieldComponent::Velocity { axis: 0 },
        AcousticFieldComponent::Velocity { axis: 1 },
    ),
    (
        AcousticFieldComponent::Velocity { axis: 1 },
        AcousticFieldComponent::Velocity { axis: 1 },
    ),
];

fn limits() -> HarmonicProductSpectrumLimits {
    HarmonicProductSpectrumLimits {
        maximum_products: 32,
        maximum_generated_terms: 64,
        maximum_workspace_bytes: 64 * core::mem::size_of::<HarmonicProductTone>(),
    }
}

fn prepare(lines: &[Line]) -> PreparedHarmonicField<2> {
    PreparedHarmonicField::new(
        lines,
        HarmonicObservationBudget {
            max_lines: lines.len(),
            max_line_pairs: lines.len() * lines.len(),
        },
    )
    .unwrap()
}

fn scale(line: Line, value: C64) -> Line {
    Line {
        pressure_pa: line.pressure_pa * value,
        particle_velocity_m_s: line.particle_velocity_m_s.map(|a| a * value),
        ..line
    }
}

fn fixture() -> ([Line; 2], [Line; 2], [Line; 2]) {
    let mean = [
        Line {
            frequency_hz: 1.25,
            pressure_pa: C64::new(0.7, -0.2),
            particle_velocity_m_s: [C64::new(0.1, 0.4), C64::new(-0.3, 0.2)],
        },
        Line {
            frequency_hz: 2.5,
            pressure_pa: C64::new(-0.3, 0.5),
            particle_velocity_m_s: [C64::new(0.2, -0.1), C64::new(0.5, 0.3)],
        },
    ];
    let real = [
        scale(mean[0], C64::new(0.4, 0.2)),
        scale(mean[1], C64::new(-0.1, 0.6)),
    ];
    let proper = [
        scale(mean[0], C64::new(-0.7, 0.1)),
        scale(mean[1], C64::new(0.2, -0.8)),
    ];
    (mean, real, proper)
}

fn ensemble(mean: &[Line], real: &[Line], proper: &[Line]) -> PreparedHarmonicEnsemble<2> {
    PreparedHarmonicEnsemble::new(
        mean,
        &[
            HarmonicFieldFactor {
                lines: real,
                statistics: HarmonicFactorStatistics::Real,
            },
            HarmonicFieldFactor {
                lines: proper,
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
        ],
        HarmonicEnsembleBudget {
            max_factors: 2,
            max_total_lines: mean.len() + real.len() + proper.len(),
            max_line_pairs: mean.len() * mean.len()
                + real.len() * real.len()
                + proper.len() * proper.len(),
        },
    )
    .unwrap()
}

fn real_field(lines: &[Line], coefficient: C64, time: f64) -> [f64; 3] {
    let mut result = [0.0; 3];
    for line in lines {
        let (s, c) = (TAU * line.frequency_hz * time).sin_cos();
        for (value, amplitude) in result.iter_mut().zip([
            line.pressure_pa,
            line.particle_velocity_m_s[0],
            line.particle_velocity_m_s[1],
        ]) {
            let a = amplitude * coefficient;
            *value += a.re * c + a.im * s;
        }
    }
    result
}

fn products([p, u, v]: [f64; 3]) -> [f64; 6] {
    [p * p, p * u, p * v, u * u, u * v, v * v]
}

fn moments(value: &WindowedAcousticEnsemble<2>) -> [f64; 6] {
    [
        value.mean_square_pressure_pa2,
        value.mean_intensity_w_m2[0],
        value.mean_intensity_w_m2[1],
        value.velocity_second_moment_m2_s2[0][0],
        value.velocity_second_moment_m2_s2[0][1],
        value.velocity_second_moment_m2_s2[1][1],
    ]
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn signed_standing_wave_flux_and_broadband_beats_match_actual_real_fields() {
    // Equal counterpropagating waves can have pressure and velocity in
    // quadrature. Their instantaneous flux reverses while its DC part is zero.
    let standing = [Line {
        frequency_hz: 2.0,
        pressure_pa: C64::from_re(2.0),
        particle_velocity_m_s: [C64::new(0.0, 0.5), C64::ZERO],
    }];
    let owner = prepare(&standing);
    let flux = owner
        .product_spectrum(PRODUCTS[1].0, PRODUCTS[1].1, limits())
        .unwrap();
    assert_eq!(flux.spectrum().dc_dd().to_f64(), 0.0);
    close(flux.evaluate(0.0625).unwrap(), 0.5);
    close(flux.evaluate(0.1875).unwrap(), -0.5);

    let broadband = [
        standing[0],
        Line {
            frequency_hz: 3.0,
            pressure_pa: C64::new(0.4, -0.3),
            particle_velocity_m_s: [C64::new(-0.2, 0.1), C64::new(0.3, 0.2)],
        },
    ];
    let owner = prepare(&broadband);
    let spectra = PRODUCTS.map(|(a, b)| owner.product_spectrum(a, b, limits()).unwrap());
    assert!(
        spectra[1]
            .spectrum()
            .tones()
            .iter()
            .any(|tone| tone.frequency_hz_dd().to_f64() == 1.0)
    );
    assert!(
        spectra[1]
            .spectrum()
            .tones()
            .iter()
            .any(|tone| tone.frequency_hz_dd().to_f64() == 5.0)
    );
    let reversed = owner
        .product_spectrum(PRODUCTS[1].1, PRODUCTS[1].0, limits())
        .unwrap();
    for index in 0..97 {
        let time = -0.23 + index as f64 * 0.013;
        let expected = products(real_field(&broadband, C64::ONE, time));
        for (spectrum, expected) in spectra.iter().zip(expected) {
            close(spectrum.evaluate(time).unwrap(), expected);
        }
        close(reversed.evaluate(time).unwrap(), expected[1]);
    }
    assert_eq!(spectra[1].components(), [PRODUCTS[1].0, PRODUCTS[1].1]);
}

#[test]
fn mixed_statistics_spectrum_matches_the_exact_eight_state_ensemble_at_each_time() {
    let (mean, real, proper) = fixture();
    let owner = ensemble(&mean, &real, &proper);
    let spectra = PRODUCTS.map(|(a, b)| owner.product_spectrum(a, b, limits()).unwrap());
    // These eight equally likely states have exactly the authored first,
    // covariance, mixed and pseudocovariance moments; no sampling assumption.
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    for time in [-0.31, -0.13, 0.0, 0.17, 0.43, 0.91] {
        let mut expected = [0.0; 6];
        for sign in [-1.0, 1.0] {
            for phase in phases {
                let m = real_field(&mean, C64::ONE, time);
                let r = real_field(&real, C64::from_re(sign), time);
                let c = real_field(&proper, phase, time);
                let field = core::array::from_fn(|i| m[i] + r[i] + c[i]);
                for (sum, value) in expected.iter_mut().zip(products(field)) {
                    *sum += value / 8.0;
                }
            }
        }
        for (spectrum, value) in spectra.iter().zip(expected) {
            close(spectrum.evaluate(time).unwrap(), value);
        }
    }

    // A phase-locked real factor has the deterministic second harmonic.
    // A proper factor's expected square is constant, even at subcycle times.
    let locked = [Line {
        frequency_hz: 2.0,
        pressure_pa: C64::ONE,
        particle_velocity_m_s: [C64::ZERO; 2],
    }];
    let real_owner = ensemble(&[], &locked, &[]);
    let proper_owner = ensemble(&[], &[], &locked);
    let real_square = real_owner
        .product_spectrum(PRODUCTS[0].0, PRODUCTS[0].1, limits())
        .unwrap();
    let proper_square = proper_owner
        .product_spectrum(PRODUCTS[0].0, PRODUCTS[0].1, limits())
        .unwrap();
    close(real_square.evaluate(0.0).unwrap(), 1.0);
    close(real_square.evaluate(0.125).unwrap(), 0.0);
    close(proper_square.evaluate(0.0).unwrap(), 0.5);
    close(proper_square.evaluate(0.125).unwrap(), 0.5);
    assert!(proper_square.spectrum().tones().is_empty());
}

#[test]
fn spectrum_averages_match_original_owners_for_both_physical_shutters() {
    let (mean, real, proper) = fixture();
    let deterministic = prepare(&mean);
    let mixed = ensemble(&mean, &real, &proper);
    let proper_owner = PreparedProperHarmonicField::new(
        &[&proper],
        ProperHarmonicObservationBudget {
            max_modes: 1,
            max_total_lines: 2,
            max_line_pairs: 4,
        },
    )
    .unwrap();
    let deterministic_spectra =
        PRODUCTS.map(|(a, b)| deterministic.product_spectrum(a, b, limits()).unwrap());
    let mixed_spectra = PRODUCTS.map(|(a, b)| mixed.product_spectrum(a, b, limits()).unwrap());
    let proper_spectra =
        PRODUCTS.map(|(a, b)| proper_owner.product_spectrum(a, b, limits()).unwrap());
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        for (start, duration) in [(0.17, 0.43), (-0.11, 0.003)] {
            let window = HarmonicWindow::with_profile(start, duration, profile).unwrap();
            let d = deterministic.observe_window_moments(window).unwrap();
            let m = mixed.observe_window(window).unwrap();
            let p = proper_owner.observe_window(window).unwrap();
            let proper_values = [
                p.mean_square_pressure_pa2,
                p.mean_intensity_w_m2[0],
                p.mean_intensity_w_m2[1],
                p.velocity_covariance_m2_s2[0][0],
                p.velocity_covariance_m2_s2[0][1],
                p.velocity_covariance_m2_s2[1][1],
            ];
            for (spectra, values) in [
                (&deterministic_spectra, moments(&d)),
                (&mixed_spectra, moments(&m)),
                (&proper_spectra, proper_values),
            ] {
                for (spectrum, expected) in spectra.iter().zip(values) {
                    close(spectrum.observe_window(window).unwrap(), expected);
                }
            }
        }
    }
}

#[test]
fn original_grouped_dd_fields_survive_budget_and_cancellation_refusals() {
    let lines = [
        Line {
            frequency_hz: 2.0,
            pressure_pa: C64::from_re(1e16),
            particle_velocity_m_s: [C64::from_re(1e16); 2],
        },
        Line {
            frequency_hz: 2.0,
            pressure_pa: C64::ONE,
            particle_velocity_m_s: [C64::from_re(-2.0), C64::from_re(0.125)],
        },
        Line {
            frequency_hz: 2.0,
            pressure_pa: C64::from_re(-1e16),
            particle_velocity_m_s: [C64::from_re(-1e16); 2],
        },
    ];
    let owner = prepare(&lines);
    assert_eq!(owner.frequency_count(), 1);
    let exact = HarmonicProductSpectrumLimits {
        maximum_products: 1,
        maximum_generated_terms: 2,
        maximum_workspace_bytes: 2 * core::mem::size_of::<HarmonicProductTone>(),
    };
    let spectrum = owner
        .product_spectrum(PRODUCTS[1].0, PRODUCTS[1].1, exact)
        .unwrap();
    assert_eq!(spectrum.evaluate(0.0).unwrap(), -2.0);
    let window = HarmonicWindow::new(0.02, 0.07).unwrap();
    let before = owner.observe_window_moments(window).unwrap();
    close(
        spectrum.observe_window(window).unwrap(),
        before.mean_intensity_w_m2[0],
    );
    for budget in [
        HarmonicProductSpectrumLimits {
            maximum_products: 0,
            ..exact
        },
        HarmonicProductSpectrumLimits {
            maximum_generated_terms: 1,
            ..exact
        },
        HarmonicProductSpectrumLimits {
            maximum_workspace_bytes: 0,
            ..exact
        },
    ] {
        assert!(
            owner
                .product_spectrum(PRODUCTS[1].0, PRODUCTS[1].1, budget)
                .is_err()
        );
    }
    assert!(
        owner
            .product_spectrum(
                AcousticFieldComponent::Pressure,
                AcousticFieldComponent::Velocity { axis: 2 },
                exact
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.product_spectrum_cancellable(PRODUCTS[1].0, PRODUCTS[1].1, exact, || {
            polls += 1;
            polls == 5
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        spectrum.evaluate_cancellable(0.0, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        spectrum.observe_window_cancellable(window, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(spectrum.evaluate(f64::NAN).is_err());
    assert_eq!(
        moments(&before),
        moments(&owner.observe_window_moments(window).unwrap())
    );
    assert_eq!(spectrum.evaluate(0.0).unwrap(), -2.0);
    let silence = prepare(&[])
        .product_spectrum(
            PRODUCTS[0].0,
            PRODUCTS[0].1,
            HarmonicProductSpectrumLimits {
                maximum_products: 0,
                maximum_generated_terms: 0,
                maximum_workspace_bytes: 0,
            },
        )
        .unwrap();
    assert_eq!(silence.evaluate(0.37).unwrap(), 0.0);
}
