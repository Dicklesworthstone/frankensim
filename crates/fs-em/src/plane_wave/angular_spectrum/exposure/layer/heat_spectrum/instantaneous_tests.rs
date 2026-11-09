//! Instantaneous heat keeps coherent reservoir fields until the final square.
use super::*;

fn limits() -> HarmonicProductSpectrumLimits {
    HarmonicProductSpectrumLimits {
        maximum_products: 256,
        maximum_generated_terms: 512,
        maximum_workspace_bytes: 512 * core::mem::size_of::<HarmonicProductTone>(),
    }
}

fn append(
    banks: &mut Vec<HeatBank>,
    fields: &mut Vec<HeatLine>,
    law: HarmonicProductLaw,
    left_scale: Dd,
    lines: impl IntoIterator<Item = HeatLine>,
) {
    let start = fields.len();
    fields.extend(lines);
    banks.push(HeatBank {
        start,
        end: fields.len(),
        law,
        left_scale,
    });
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && actual > 0.0, "actual={actual:e}");
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn instantaneous_real_and_proper_heat_resolve_three_tone_dark_fields() {
    for law in [HarmonicProductLaw::Real, HarmonicProductLaw::ProperComplex] {
        let proper = law == HarmonicProductLaw::ProperComplex;
        for epoch in [0.0, 0.375] {
            let carrier = if proper { 562_949_953_421_312.0 } else { 1.0 };
            let mut banks = Vec::new();
            let mut fields = Vec::new();
            append(
                &mut banks,
                &mut fields,
                law,
                Dd::ONE,
                [1.0, -2.0, 1.0]
                    .into_iter()
                    .enumerate()
                    .map(|(i, amplitude)| {
                        // Author a null at epoch. A proper factor needs only its
                        // relative source phase; its common carrier cancels.
                        let authored_omega = if proper { i as f64 } else { carrier + i as f64 };
                        let (sine, cosine) =
                            det::sin_cos_dd(Dd::from_f64(authored_omega) * Dd::from_f64(epoch));
                        HeatLine {
                            omega: carrier + i as f64,
                            field: [
                                [
                                    cosine * Dd::from_f64(amplitude),
                                    sine * Dd::from_f64(amplitude),
                                ],
                                [Dd::ZERO; 2],
                                [Dd::ZERO; 2],
                            ],
                        }
                    }),
            );
            // The original signed-spectrum preparation still admits the owner.
            prepare_products(&banks, &fields, limits(), &mut || false).unwrap();
            for offset in if epoch == 0.0 {
                [1e-9, 1e-12]
            } else {
                [1e-8, 1e-10]
            } {
                let time = epoch + offset;
                let delta = (Dd::from_f64(time) - Dd::from_f64(epoch)).to_f64();
                // (1 - exp(-i delta))² times the common carrier. This oracle
                // never subtracts bright cosine powers near the dark field.
                let beat = 4.0 * (0.5 * delta).sin().powi(2);
                let expected = if proper {
                    0.5 * beat * beat
                } else {
                    (beat * (2.0 * delta).cos()).powi(2)
                };
                relative(
                    positive_heat_at_time(&banks, &fields, time, &mut || false).unwrap(),
                    expected,
                    3e-6,
                );
            }
            if epoch == 0.0 {
                assert_eq!(
                    positive_heat_at_time(&banks, &fields, 0.0, &mut || false).unwrap(),
                    0.0
                );
            }
        }
    }
}

#[test]
fn instantaneous_heat_keeps_native_sign_and_independent_physical_reservoirs() {
    let omega = [2.0_f64, 2.3];
    let sigma = 0.3;
    let delta_debye = 0.7_f64;
    let tau = 0.4;
    let delta_lorentz = 0.2_f64;
    let omega0 = 3.0;
    let gamma = 0.5;
    let electric = [
        [
            [C64::new(0.8, 0.2), C64::new(-0.1, 0.3), C64::new(0.2, -0.4)],
            [C64::new(-0.3, 0.6), C64::new(0.4, 0.1), C64::new(-0.2, 0.1)],
        ],
        [
            [
                C64::new(-0.8, -0.2),
                C64::new(0.2, -0.1),
                C64::new(-0.1, 0.2),
            ],
            [C64::new(0.3, -0.6), C64::new(0.1, 0.4), C64::new(0.2, -0.1)],
        ],
        [
            [C64::new(0.3, 0.1), C64::new(-0.2, 0.4), C64::new(0.1, -0.2)],
            [
                C64::new(-0.4, -0.2),
                C64::new(0.3, -0.1),
                C64::new(-0.1, 0.3),
            ],
        ],
    ];
    let mut banks = Vec::new();
    let mut fields = Vec::new();
    let root = Dd::from_f64(det::sqrt(sigma));
    for (source, electric) in electric.iter().enumerate() {
        for reservoir in 0..3 {
            append(
                &mut banks,
                &mut fields,
                if source == 2 {
                    HarmonicProductLaw::ProperComplex
                } else {
                    HarmonicProductLaw::Real
                },
                if reservoir == 0 {
                    Dd::from_f64(sigma) / (root * root)
                } else {
                    Dd::ONE
                },
                omega.into_iter().enumerate().map(|(i, omega)| {
                    let transfer = match reservoir {
                        0 => C64::from_re(root.to_f64()),
                        1 => {
                            C64::new(0.0, omega * (tau * delta_debye).sqrt())
                                / C64::new(1.0, omega * tau)
                        }
                        _ => {
                            C64::new(0.0, omega * omega0 * (gamma * delta_lorentz).sqrt())
                                / C64::new(omega0 * omega0 - omega * omega, gamma * omega)
                        }
                    };
                    HeatLine {
                        omega,
                        field: electric[i].map(|e| filtered_negative_time_parts(e, transfer)),
                    }
                }),
            );
        }
    }
    prepare_products(
        &banks,
        &fields,
        HarmonicProductSpectrumLimits {
            maximum_products: 512,
            maximum_generated_terms: 1024,
            maximum_workspace_bytes: 1024 * core::mem::size_of::<HarmonicProductTone>(),
        },
        &mut || false,
    )
    .unwrap();
    for time in [-0.3, 0.0, 0.17, 0.49, 1.7] {
        let mut expected = 0.0;
        // Enumerate the exact real +/- and proper four-quadrature states.
        // Reservoir currents use the original ODE, independently of the
        // production filtered-field helper and source-bank square.
        for sign in [-1.0, 1.0] {
            for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
                let mut e = [0.0; 3];
                let mut debye_current = [0.0; 3];
                let mut lorentz_current = [0.0; 3];
                for (source, coefficient) in electric.iter().zip([C64::ONE, C64::from_re(sign), z])
                {
                    for (i, omega) in omega.into_iter().enumerate() {
                        let (sine, cosine) = (omega * time).sin_cos();
                        let real = |a: C64| a.re * cosine - a.im * sine;
                        for axis in 0..3 {
                            let amplitude = source[i][axis] * coefficient;
                            let pd = amplitude.scale(delta_debye) / C64::new(1.0, omega * tau);
                            let pl = amplitude.scale(delta_lorentz * omega0 * omega0)
                                / C64::new(omega0 * omega0 - omega * omega, gamma * omega);
                            e[axis] += real(amplitude);
                            debye_current[axis] += real(C64::new(-pd.im, pd.re).scale(omega));
                            lorentz_current[axis] += real(C64::new(-pl.im, pl.re).scale(omega));
                        }
                    }
                }
                expected += (sigma * e.iter().map(|v| v * v).sum::<f64>()
                    + (tau / delta_debye) * debye_current.iter().map(|v| v * v).sum::<f64>()
                    + (gamma / (delta_lorentz * omega0 * omega0))
                        * lorentz_current.iter().map(|v| v * v).sum::<f64>())
                    / 8.0;
            }
        }
        relative(
            positive_heat_at_time(&banks, &fields, time, &mut || false).unwrap(),
            expected,
            5e-14,
        );
    }
}

#[test]
fn instantaneous_heat_applies_the_physical_scale_before_large_and_tiny_squares() {
    for law in [HarmonicProductLaw::Real, HarmonicProductLaw::ProperComplex] {
        let proper = law == HarmonicProductLaw::ProperComplex;
        for amplitude in [if proper { 1.5e154 } else { 1e154 }, 1e-161] {
            let mut banks = Vec::new();
            let mut fields = Vec::new();
            append(
                &mut banks,
                &mut fields,
                law,
                Dd::ONE,
                [HeatLine {
                    omega: 2.0,
                    field: [
                        [Dd::from_f64(amplitude), Dd::ZERO],
                        [Dd::ZERO; 2],
                        [Dd::ZERO; 2],
                    ],
                }],
            );
            prepare_products(&banks, &fields, limits(), &mut || false).unwrap();
            for time in [0.0, 0.17] {
                let expected = if proper {
                    (0.5 * amplitude) * amplitude
                } else {
                    (amplitude * (2.0_f64 * time).cos()).powi(2)
                };
                let actual = positive_heat_at_time(&banks, &fields, time, &mut || false).unwrap();
                assert!(actual.is_finite() && actual > 0.0);
                assert!(
                    (actual - expected).abs() <= 4.0 * f64::from_bits(1) + 2e-14 * expected,
                    "actual={actual:e}, expected={expected:e}"
                );
            }
        }
    }
}

#[test]
fn instantaneous_heat_refuses_invalid_time_and_cancellation_without_mutating_fields() {
    let fields = [HeatLine {
        omega: 2.0,
        field: [[Dd::ONE, Dd::from_f64(0.5)], [Dd::ZERO; 2], [Dd::ZERO; 2]],
    }];
    let banks = [HeatBank {
        start: 0,
        end: 1,
        law: HarmonicProductLaw::Real,
        left_scale: Dd::from_f64(0.3),
    }];
    let mut polls = 0;
    let before = positive_heat_at_time(&banks, &fields, 0.17, &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut observed = 0;
    assert!(matches!(
        positive_heat_at_time(&banks, &fields, 0.17, &mut || {
            observed += 1;
            observed == polls
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        before,
        positive_heat_at_time(&banks, &fields, 0.17, &mut || false).unwrap()
    );
    for time in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(positive_heat_at_time(&banks, &fields, time, &mut || false).is_err());
        assert!(positive_heat_at_time(&[], &[], time, &mut || false).is_err());
    }
    assert_eq!(
        positive_heat_at_time(&[], &[], 0.17, &mut || false).unwrap(),
        0.0
    );
}
