use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e} expected={expected:.17e} tolerance={tolerance:.3e}"
    );
}

#[test]
fn g0_published_absolute_indices_cover_visible_infrared_and_temperature_endpoints() {
    // NASA NTRS 20070018851, Table 4. The source prints five decimal places;
    // this checks original absolute indices rather than a relative-air fit.
    let rows = [
        (400.0, [1.46899, 1.47054, 1.47059]),
        (500.0, [1.46129, 1.46275, 1.46279]),
        (600.0, [1.45704, 1.45846, 1.45850]),
        (800.0, [1.45236, 1.45373, 1.45377]),
        (1000.0, [1.44947, 1.45082, 1.45086]),
        (1500.0, [1.44369, 1.44502, 1.44507]),
        (2600.0, [1.42688, 1.42823, 1.42828]),
    ];
    for (column, temperature) in [30.0, 295.0, 300.0].into_iter().enumerate() {
        let law = Corning7980Silica::at_temperature(temperature).unwrap();
        for (wavelength, indices) in rows {
            close(
                law.evaluate(wavelength).unwrap().phase_index,
                indices[column],
                5.1e-6,
            );
        }
    }
}

#[test]
fn g1_analytic_dispersion_and_thermal_jet_match_independent_high_precision_controls() {
    // 70-digit Decimal evaluation of the printed Table-3 polynomials using
    // explicit powers; five-point numerical derivatives with a 1e-8 K/um
    // step. This oracle does not reuse the production Horner or derivative
    // formulas. The values are controls of the fit, not new measurements.
    let controls = [
        (
            30.0,
            400.0,
            [
                1.468_993_543_706_87,
                1.512_181_464_956_940_7,
                -1.079_698_031_251_771_2e-4,
                3.381_490_402_126_881e-6,
                -1.281_938_661_746_341e-7,
            ],
        ),
        (
            300.0,
            2600.0,
            [
                1.428_277_318_914_038_4,
                1.479_727_471_218_672_5,
                -1.978_852_011_716_696_8e-5,
                9.457_631_451_406_03e-6,
                6.333_263_814_503_085e-8,
            ],
        ),
        (
            222.25,
            1550.0,
            [
                1.443_885_066_753_925_3,
                1.462_442_027_459_719_8,
                -1.197_223_271_341_578_6e-5,
                6.762_005_789_946_712e-6,
                2.299_673_316_144_563_3e-8,
            ],
        ),
    ];
    for (temperature, wavelength, expected) in controls {
        let sample = Corning7980Silica::at_temperature(temperature)
            .unwrap()
            .evaluate(wavelength)
            .unwrap();
        for ((actual, expected), tolerance) in [
            sample.phase_index,
            sample.group_index,
            sample.wavelength_derivative_per_nm,
            sample.derivative_temperature_per_k,
            sample.second_derivative_temperature_per_k2,
        ]
        .into_iter()
        .zip(expected)
        .zip([7e-15, 2e-14, 2e-17, 1e-18, 1e-20])
        {
            close(actual, expected, tolerance);
        }
    }
}

#[test]
fn g3_derivatives_are_of_the_original_temperature_and_wavelength_law() {
    for temperature in [35.0, 160.0, 222.25, 295.0] {
        let law = Corning7980Silica::at_temperature(temperature).unwrap();
        for wavelength in [450.0, 1000.0, 1550.0, 2550.0] {
            let sample = law.evaluate(wavelength).unwrap();
            let phase = |t: f64| {
                Corning7980Silica::at_temperature(t)
                    .unwrap()
                    .evaluate(wavelength)
                    .unwrap()
                    .phase_index
            };
            let h = 0.25;
            let first = (phase(temperature - 2.0 * h) - 8.0 * phase(temperature - h)
                + 8.0 * phase(temperature + h)
                - phase(temperature + 2.0 * h))
                / (12.0 * h);
            let second = (-phase(temperature + 2.0 * h) + 16.0 * phase(temperature + h)
                - 30.0 * phase(temperature)
                + 16.0 * phase(temperature - h)
                - phase(temperature - 2.0 * h))
                / (12.0 * h * h);
            close(sample.derivative_temperature_per_k, first, 2e-13);
            close(sample.second_derivative_temperature_per_k2, second, 2e-13);

            let phase = |w: f64| law.evaluate(w).unwrap().phase_index;
            let h = 0.1;
            let dispersion = (phase(wavelength - 2.0 * h) - 8.0 * phase(wavelength - h)
                + 8.0 * phase(wavelength + h)
                - phase(wavelength + 2.0 * h))
                / (12.0 * h);
            close(sample.wavelength_derivative_per_nm, dispersion, 2e-14);
            close(
                sample.group_index,
                sample.phase_index - wavelength * dispersion,
                6e-11,
            );
        }
    }
}

#[test]
fn g0_source_domain_frozen_coefficients_and_identity_remain_explicit() {
    for temperature in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        29.999,
        300.001,
        310.0,
    ] {
        assert_eq!(
            Corning7980Silica::at_temperature(temperature),
            Err(SilicaRefractionError::TemperatureOutsideRange)
        );
    }
    let source = Corning7980Silica::source_identity();
    assert_ne!(source, ContentHash([0; 32]));
    assert_eq!(source, Corning7980Silica::source_identity());
    for temperature in [30.0, 160.0, 300.0] {
        let law = Corning7980Silica::at_temperature(temperature).unwrap();
        assert_eq!(law.temperature_k(), temperature);
        assert_eq!(law.wavelength_bounds_nm(), (400.0, 2600.0));
        assert_eq!(
            law.identity(),
            Corning7980Silica::at_temperature(temperature)
                .unwrap()
                .identity()
        );
        assert_ne!(
            law.identity(),
            Corning7980Silica::at_temperature(
                temperature - 1.0 + if temperature == 30.0 { 2.0 } else { 0.0 }
            )
            .unwrap()
            .identity()
        );
        let (strengths, resonances_um2) = law.sellmeier_coefficients();
        assert!(strengths.into_iter().all(|s| s.is_finite() && s > 0.0));
        assert!(
            resonances_um2
                .into_iter()
                .all(|c| c > 0.0 && !(0.16..=6.76).contains(&c))
        );
        for wavelength in [400.0, 2600.0] {
            let sample = law.evaluate(wavelength).unwrap();
            assert!(sample.phase_index > 1.0 && sample.group_index > sample.phase_index);
        }
        for wavelength in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            0.0,
            399.999,
            2600.001,
        ] {
            assert_eq!(
                law.evaluate(wavelength),
                Err(SilicaRefractionError::WavelengthOutsideRange)
            );
        }
    }
}
