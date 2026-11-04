use super::*;

// Synthetic coefficients test the equation, not measured glass qualification.
const COEFFICIENTS: ThermoOpticCoefficients = ThermoOpticCoefficients {
    d: [1.2e-5, -2e-8, 3e-11],
    e: [4e-6, -5e-9],
    resonance_um: 0.18,
};
const SOURCE: ContentHash = ContentHash([17; 32]);

fn law() -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        COEFFICIENTS,
        300.0,
        [200.0, 500.0],
        [800.0, 2_000.0],
        SOURCE,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e} expected={expected:.17e}"
    );
}

#[test]
fn original_equation_matches_independent_high_precision_reference() {
    // 80-decimal-digit evaluation of unfactored TIE-19 eq. (3), followed by
    // independent symbolic temperature and wavelength differentiation.
    let sample = law()
        .prepare(1_550.0, 1.53, Some(1.57))
        .unwrap()
        .sample(337.5)
        .unwrap();
    close(sample.phase_index, 1.530_211_993_158_215_9, 3e-16);
    close(sample.group_index.unwrap(), 1.570_279_391_067_132_8, 3e-16);
    close(
        sample.derivative_temperature_per_k,
        5.326_805_804_928_410e-6,
        4e-21,
    );
    close(
        sample.second_derivative_temperature_per_k2,
        -1.641_911_509_518_227_2e-8,
        2e-23,
    );
}

#[test]
fn temperature_derivatives_match_independent_polynomial_and_finite_differences() {
    let prepared = law().prepare(1_550.0, 1.53, Some(1.57)).unwrap();
    let [d0, d1, d2] = COEFFICIENTS.d;
    let [e0, e1] = COEFFICIENTS.e;
    let q = 1.55_f64.powi(2) - COEFFICIENTS.resonance_um.powi(2);
    let amplitude = (1.53_f64.powi(2) - 1.0) / (2.0 * 1.53);
    let cubic = [
        amplitude * (d0 + e0 / q),
        amplitude * (d1 + e1 / q),
        amplitude * d2,
    ];
    for temperature in [210.0, 275.0, 300.0, 337.5, 490.0] {
        let dt = temperature - 300.0;
        let sample = prepared.sample(temperature).unwrap();
        let [a, b, c] = cubic;
        close(
            sample.phase_index,
            1.53 + a * dt + b * dt * dt + c * dt * dt * dt,
            5e-16,
        );
        close(
            sample.derivative_temperature_per_k,
            a + 2.0 * b * dt + 3.0 * c * dt * dt,
            8e-21,
        );
        close(
            sample.second_derivative_temperature_per_k2,
            2.0 * b + 6.0 * c * dt,
            3e-23,
        );
        let phase = |t| prepared.sample(t).unwrap().phase_index;
        let h = 0.5;
        let first = (phase(temperature - 2.0 * h) - 8.0 * phase(temperature - h)
            + 8.0 * phase(temperature + h)
            - phase(temperature + 2.0 * h))
            / (12.0 * h);
        let second =
            (phase(temperature + h) - 2.0 * phase(temperature) + phase(temperature - h)) / (h * h);
        close(sample.derivative_temperature_per_k, first, 8e-16);
        close(sample.second_derivative_temperature_per_k2, second, 2e-15);
    }
}

#[test]
fn group_index_differentiates_reference_dispersion_and_temperature_increment() {
    let law = law();
    for temperature in [210.0, 300.0, 480.0] {
        for wavelength_nm in [850.0, 1_550.0, 1_950.0] {
            let sample = |wavelength_nm: f64| {
                let wavelength_um = wavelength_nm * 1e-3;
                // Analytic Cauchy reference: n0 = a + b/lambda^2,
                // ng0 = a + 3b/lambda^2. It is intentionally dispersive.
                let dispersion = 0.04805 / (wavelength_um * wavelength_um);
                law.prepare(
                    wavelength_nm,
                    1.51 + dispersion,
                    Some(1.51 + 3.0 * dispersion),
                )
                .unwrap()
                .sample(temperature)
                .unwrap()
            };
            let h = 0.05;
            let derivative = (sample(wavelength_nm - 2.0 * h).phase_index
                - 8.0 * sample(wavelength_nm - h).phase_index
                + 8.0 * sample(wavelength_nm + h).phase_index
                - sample(wavelength_nm + 2.0 * h).phase_index)
                / (12.0 * h);
            let state = sample(wavelength_nm);
            close(
                state.group_index.unwrap(),
                state.phase_index - wavelength_nm * derivative,
                1e-10,
            );
        }
    }
}

#[test]
fn reference_temperature_preserves_bits_and_repeated_sampling_has_no_drift() {
    let n0 = f64::from_bits(1.53_f64.to_bits() + 1);
    let ng0 = f64::from_bits(1.57_f64.to_bits() + 1);
    let prepared = law().prepare(1_550.0, n0, Some(ng0)).unwrap();
    let reference = prepared.sample(300.0).unwrap();
    assert_eq!(reference.phase_index.to_bits(), n0.to_bits());
    assert_eq!(reference.group_index.unwrap().to_bits(), ng0.to_bits());
    let hot = prepared.sample(420.0).unwrap();
    for _ in 0..100 {
        assert_eq!(prepared.sample(420.0).unwrap(), hot);
        assert_eq!(prepared.sample(300.0).unwrap(), reference);
    }
    assert_eq!(prepared.wavelength_nm(), 1_550.0);
}

#[test]
fn absent_reference_dispersion_remains_unknown_at_every_temperature() {
    let prepared = law().prepare(1_550.0, 1.53, None).unwrap();
    for temperature in [200.0, 300.0, 500.0] {
        assert_eq!(prepared.sample(temperature).unwrap().group_index, None);
    }
}

#[test]
fn vacuum_phase_can_have_nonzero_group_temperature_response() {
    // At n0=1, A=0 but dA/dlambda need not vanish: an implementation that
    // discards group corrections with a zero phase prefactor is incorrect.
    let sample = law()
        .prepare(1_550.0, 1.0, Some(1.1))
        .unwrap()
        .sample(400.0)
        .unwrap();
    assert_eq!(sample.phase_index, 1.0);
    assert_eq!(sample.derivative_temperature_per_k, 0.0);
    assert_eq!(sample.second_derivative_temperature_per_k2, 0.0);
    assert!(sample.group_index.unwrap() > 1.1);
}

#[test]
fn active_poles_refuse_the_whole_inclusive_band() {
    for resonance_um in [0.8, 1.2, 2.0] {
        for e in [[1e-6, 0.0], [0.0, -1e-9]] {
            let coefficients = ThermoOpticCoefficients {
                resonance_um,
                e,
                ..COEFFICIENTS
            };
            assert_eq!(
                SchottThermoOpticLaw::try_new(
                    coefficients,
                    300.0,
                    [200.0, 500.0],
                    [800.0, 2_000.0],
                    SOURCE
                ),
                Err(ThermoOpticError::ResonanceInBand)
            );
        }
    }
}

#[test]
fn inactive_pole_is_irrelevant_and_signed_thermal_response_is_supported() {
    let coefficients = ThermoOpticCoefficients {
        d: [-2e-5, 0.0, 0.0],
        e: [0.0; 2],
        resonance_um: 1.55,
    };
    let law = SchottThermoOpticLaw::try_new(
        coefficients,
        300.0,
        [200.0, 500.0],
        [800.0, 2_000.0],
        SOURCE,
    )
    .unwrap();
    let sample = law
        .prepare(1_550.0, 1.53, Some(1.57))
        .unwrap()
        .sample(400.0)
        .unwrap();
    assert!(sample.phase_index < 1.53);
    assert!(sample.group_index.unwrap() < 1.57);
    assert!(sample.derivative_temperature_per_k < 0.0);
    assert_eq!(sample.second_derivative_temperature_per_k2, 0.0);
}

#[test]
fn explicit_infrared_domain_and_below_resonance_domain_are_admitted() {
    let infrared = SchottThermoOpticLaw::try_new(
        COEFFICIENTS,
        600.0,
        [550.0, 650.0],
        [3_000.0, 5_000.0],
        SOURCE,
    )
    .unwrap();
    assert!(
        infrared
            .prepare(4_000.0, 1.5, Some(1.6))
            .unwrap()
            .sample(625.0)
            .is_ok()
    );
    let below = SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            resonance_um: 3.0,
            ..COEFFICIENTS
        },
        300.0,
        [200.0, 500.0],
        [800.0, 2_000.0],
        SOURCE,
    )
    .unwrap();
    assert!(
        below
            .prepare(1_550.0, 1.5, Some(1.6))
            .unwrap()
            .sample(400.0)
            .is_ok()
    );
}

#[test]
fn unsupported_queries_refuse_and_source_endpoints_are_inclusive() {
    let law = law();
    for wavelength in [800.0, 2_000.0] {
        let prepared = law.prepare(wavelength, 1.53, Some(1.57)).unwrap();
        assert!(prepared.sample(200.0).is_ok());
        assert!(prepared.sample(500.0).is_ok());
        for temperature in [199.99, 500.01, f64::NAN, f64::INFINITY, 0.0] {
            assert_eq!(
                prepared.sample(temperature),
                Err(ThermoOpticError::TemperatureOutsideRange)
            );
        }
    }
    for wavelength in [799.99, 2_000.01, f64::NAN, f64::INFINITY, 0.0] {
        assert_eq!(
            law.prepare(wavelength, 1.53, None),
            Err(ThermoOpticError::WavelengthOutsideRange)
        );
    }
    for index in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            law.prepare(1_550.0, index, None),
            Err(ThermoOpticError::InvalidReferenceIndex)
        );
        assert_eq!(
            law.prepare(1_550.0, 1.53, Some(index)),
            Err(ThermoOpticError::InvalidReferenceIndex)
        );
    }
}

#[test]
fn malformed_source_coefficients_domains_and_identity_refuse() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for position in 0..6 {
            let mut coefficients = COEFFICIENTS;
            match position {
                0..3 => coefficients.d[position] = bad,
                3..5 => coefficients.e[position - 3] = bad,
                _ => coefficients.resonance_um = bad,
            }
            assert_eq!(
                SchottThermoOpticLaw::try_new(
                    coefficients,
                    300.0,
                    [200.0, 500.0],
                    [800.0, 2_000.0],
                    SOURCE
                ),
                Err(ThermoOpticError::InvalidCoefficients)
            );
        }
    }
    assert_eq!(
        SchottThermoOpticLaw::try_new(
            ThermoOpticCoefficients {
                resonance_um: -0.18,
                ..COEFFICIENTS
            },
            300.0,
            [200.0, 500.0],
            [800.0, 2_000.0],
            SOURCE
        ),
        Err(ThermoOpticError::InvalidCoefficients)
    );
    for bounds in [
        [0.0, 500.0],
        [300.0, 300.0],
        [500.0, 200.0],
        [f64::NAN, 500.0],
        [200.0, f64::INFINITY],
    ] {
        assert_eq!(
            SchottThermoOpticLaw::try_new(COEFFICIENTS, 300.0, bounds, [800.0, 2_000.0], SOURCE),
            Err(ThermoOpticError::InvalidTemperatureRange)
        );
        assert_eq!(
            SchottThermoOpticLaw::try_new(COEFFICIENTS, 300.0, [200.0, 500.0], bounds, SOURCE),
            Err(ThermoOpticError::InvalidWavelengthRange)
        );
    }
    for reference in [0.0, 100.0, 600.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            SchottThermoOpticLaw::try_new(
                COEFFICIENTS,
                reference,
                [200.0, 500.0],
                [800.0, 2_000.0],
                SOURCE
            ),
            Err(ThermoOpticError::InvalidReferenceTemperature)
        );
    }
    assert_eq!(
        SchottThermoOpticLaw::try_new(
            COEFFICIENTS,
            300.0,
            [200.0, 500.0],
            [800.0, 2_000.0],
            ContentHash([0; 32])
        ),
        Err(ThermoOpticError::InvalidSourceIdentity)
    );
}

#[test]
fn identity_binds_every_parameter_and_canonicalizes_signed_zero() {
    let original = law();
    assert_eq!(original.reference_temperature_k(), 300.0);
    assert_eq!(original.temperature_range_k(), [200.0, 500.0]);
    assert_eq!(original.wavelength_range_nm(), [800.0, 2_000.0]);
    assert_eq!(original.source_identity(), SOURCE);
    for position in 0..12 {
        let mut coefficients = COEFFICIENTS;
        let mut reference = 300.0;
        let mut temperature = [200.0, 500.0];
        let mut wavelength = [800.0, 2_000.0];
        let mut source = SOURCE;
        match position {
            0..3 => coefficients.d[position] *= 2.0,
            3..5 => coefficients.e[position - 3] *= 2.0,
            5 => coefficients.resonance_um *= 2.0,
            6 => reference += 1.0,
            7..9 => temperature[position - 7] += 1.0,
            9..11 => wavelength[position - 9] += 1.0,
            _ => source.0[0] += 1,
        }
        let changed =
            SchottThermoOpticLaw::try_new(coefficients, reference, temperature, wavelength, source)
                .unwrap();
        assert_ne!(original.identity(), changed.identity());
    }
    let zero = |value| {
        SchottThermoOpticLaw::try_new(
            ThermoOpticCoefficients {
                d: [value; 3],
                e: [value; 2],
                resonance_um: value,
            },
            300.0,
            [200.0, 500.0],
            [800.0, 2_000.0],
            SOURCE,
        )
        .unwrap()
    };
    assert_eq!(zero(0.0), zero(-0.0));
    assert_eq!(zero(0.0).identity(), zero(-0.0).identity());
    assert_ne!(zero(0.0).identity(), ContentHash([0; 32]));
}

#[test]
fn nonphysical_or_unrepresentable_response_refuses_without_clamping() {
    let with_d = |d| {
        SchottThermoOpticLaw::try_new(
            ThermoOpticCoefficients {
                d,
                e: [0.0; 2],
                resonance_um: 0.0,
            },
            300.0,
            [200.0, 500.0],
            [800.0, 2_000.0],
            SOURCE,
        )
        .unwrap()
    };
    assert_eq!(
        with_d([-1.0, 0.0, 0.0])
            .prepare(1_550.0, 1.5, None)
            .unwrap()
            .sample(400.0),
        Err(ThermoOpticError::NonPositiveIndex)
    );
    // Phase remains exactly one, but source group dispersion becomes negative.
    assert_eq!(
        with_d([-1.0, 0.0, 0.0])
            .prepare(1_550.0, 1.0, Some(2.0))
            .unwrap()
            .sample(400.0),
        Err(ThermoOpticError::NonPositiveIndex)
    );
    assert_eq!(
        with_d([0.0, 0.0, 1e304])
            .prepare(1_550.0, 1.5, None)
            .unwrap()
            .sample(500.0),
        Err(ThermoOpticError::NumericalRange)
    );
    assert_eq!(
        with_d([f64::MAX, 0.0, 0.0]).prepare(1_550.0, 4.0, None),
        Err(ThermoOpticError::NumericalRange)
    );
}
