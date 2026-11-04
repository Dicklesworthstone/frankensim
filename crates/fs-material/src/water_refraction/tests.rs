use super::*;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97State};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:.3e}"
    );
}

#[test]
fn g1_original_release_verification_with_independent_iapws95_density_tables() {
    // Refractive indices: original IAPWS R9-97, Table 3, PDF page 6.
    // Densities: A. H. Harvey, NISTIR 5078 (1998), Table 3, printed pages
    // 22, 23, 54 and 55 (PDF pages 28, 29, 60 and 61), IAPWS-95/SRD 10.
    // https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir5078.pdf
    // These INDEPENDENT densities are printed to five significant digits.
    // The tolerance includes their half-last-digit rounding, bounded by
    // 6e-4 m3/kg for dn/drho here, plus the original index's printed rounding.
    // No fitted density is recovered from the expected refractive index.
    let cases = [
        (273.15, 999.84, 0.005, [1.394527, 1.334344, 1.326135]),
        (373.15, 0.58967, 0.000005, [1.0002168, 1.0001876, 1.0001837]),
        (473.15, 0.46031, 0.000005, [1.0001683, 1.0001456, 1.0001427]),
        (773.15, 0.28039, 0.000005, [1.0001008, 1.0000871, 1.0000856]),
        (273.15, 1004.8, 0.05, [1.396526, 1.335969, 1.327710]),
        (373.15, 962.93, 0.005, [1.377286, 1.320084, 1.312577]),
        (473.15, 870.94, 0.005, [1.338299, 1.287891, 1.281529]),
        (773.15, 30.477, 0.0005, [1.0109906, 1.0094939, 1.0093267]),
    ];
    for (temperature, density, density_rounding, indices) in cases {
        for (wavelength, expected) in [226.5, 589.0, 1013.98].into_iter().zip(indices) {
            let sample = IapwsWaterRefraction::prepare(wavelength)
                .unwrap()
                .sample(temperature, density)
                .unwrap();
            let index_rounding = if expected < 1.1 { 5.0e-8 } else { 5.0e-7 };
            near(
                sample.phase_index,
                expected,
                6.0e-4 * density_rounding + index_rounding,
            );
            assert!(sample.group_index > sample.phase_index);
        }
    }
}

#[test]
fn g1_phase_group_and_state_derivatives_match_independent_differences() {
    // The independent numerical differences vary ORIGINAL source inputs;
    // there is no derivative implementation or coefficient copy in the oracle.
    for wavelength in [210.0, 400.0, 589.0, 1000.0, 1090.0] {
        let prepared = IapwsWaterRefraction::prepare(wavelength).unwrap();
        for (temperature, density) in [
            (275.0, 1030.0),
            (300.0, 996.0),
            (500.0, 12.0),
            (760.0, 450.0),
        ] {
            let sample = prepared.sample(temperature, density).unwrap();
            let difference = |offset: f64, f: &dyn Fn(f64) -> f64| {
                (f(-2.0 * offset) - 8.0 * f(-offset) + 8.0 * f(offset) - f(2.0 * offset))
                    / (12.0 * offset)
            };
            let wavelength_difference = difference(0.05, &|d| {
                IapwsWaterRefraction::prepare(wavelength + d)
                    .unwrap()
                    .sample(temperature, density)
                    .unwrap()
                    .phase_index
            });
            near(
                sample.wavelength_derivative_per_nm,
                wavelength_difference,
                7.0e-15,
            );
            near(
                sample.group_index,
                sample.phase_index - wavelength * wavelength_difference,
                8.0e-12,
            );
            let by_temperature = |d| prepared.sample(temperature + d, density).unwrap();
            let by_density = |d| prepared.sample(temperature, density + d).unwrap();
            near(
                sample.temperature_derivative_per_k,
                difference(0.1, &|d| by_temperature(d).phase_index),
                3.0e-15,
            );
            near(
                sample.density_derivative_m3_per_kg,
                difference(0.1, &|d| by_density(d).phase_index),
                3.0e-15,
            );
            near(
                sample.temperature_second_derivative_per_k2,
                difference(0.1, &|d| by_temperature(d).temperature_derivative_per_k),
                4.0e-19,
            );
            near(
                sample.density_second_derivative_m6_per_kg2,
                difference(0.1, &|d| by_density(d).density_derivative_m3_per_kg),
                4.0e-18,
            );
            near(
                sample.temperature_density_derivative_m3_per_kg_k,
                difference(0.1, &|d| by_temperature(d).density_derivative_m3_per_kg),
                4.0e-18,
            );
            near(
                sample.temperature_density_derivative_m3_per_kg_k,
                difference(0.1, &|d| by_density(d).temperature_derivative_per_k),
                4.0e-19,
            );
        }
    }
}

#[test]
fn g1_original_if97_states_retain_eos_density_and_isobaric_thermal_response() {
    let prepared = IapwsWaterRefraction::prepare(589.0).unwrap();
    for (pressure, temperature, region) in [
        (0.1e6, 300.0, If97Region::Region1),
        (10.0e6, 373.15, If97Region::Region1),
        (0.1e6, 500.0, If97Region::Region2),
        (10.0e6, 773.05, If97Region::Region2),
    ] {
        let state = If97State::at_pressure_temperature(
            Pressure::new(pressure),
            Temperature::new(temperature),
            region,
        )
        .unwrap();
        let actual = prepared.sample_if97_state(&state).unwrap();
        assert_eq!(
            actual,
            prepared
                .sample(temperature, state.density().value())
                .unwrap()
        );
        // d(n)/dT at fixed P includes the SAME EOS's density change. A constant
        // density optical partial alone must not masquerade as this derivative.
        let density_slope = -state.density().value() * state.thermal_expansion().value();
        let total = actual
            .density_derivative_m3_per_kg
            .mul_add(density_slope, actual.temperature_derivative_per_k);
        let shifted = |d| {
            let state = If97State::at_pressure_temperature(
                Pressure::new(pressure),
                Temperature::new(temperature + d),
                region,
            )
            .unwrap();
            prepared.sample_if97_state(&state).unwrap().phase_index
        };
        near(
            total,
            (shifted(0.05) - shifted(-0.05)) / 0.1,
            total.abs() * 3.0e-6 + 5.0e-13,
        );
        assert!((total - actual.temperature_derivative_per_k).abs() > total.abs() * 0.05);
    }
    let too_hot = If97State::at_pressure_temperature(
        Pressure::new(1.0e6),
        Temperature::new(800.0),
        If97Region::Region2,
    )
    .unwrap();
    assert_eq!(
        prepared.sample_if97_state(&too_hot),
        Err(WaterRefractionError::TemperatureOutsideRange)
    );
}

#[test]
fn g0_source_domain_and_vacuum_limit_remain_explicit() {
    assert_ne!(
        IapwsWaterRefraction::source_identity(),
        ContentHash([0; 32])
    );
    assert_eq!(
        IapwsWaterRefraction::source_identity(),
        IapwsWaterRefraction::source_identity()
    );
    for wavelength in [199.99, 1100.01, 0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            IapwsWaterRefraction::prepare(wavelength),
            Err(WaterRefractionError::WavelengthOutsideRange)
        );
    }
    for wavelength in WAVELENGTH_RANGE_NM {
        let prepared = IapwsWaterRefraction::prepare(wavelength).unwrap();
        assert_eq!(prepared.wavelength_nm(), wavelength);
        for temperature in TEMPERATURE_RANGE_K {
            for density in DENSITY_RANGE_KG_M3 {
                prepared.sample(temperature, density).unwrap();
            }
        }
    }
    let prepared = IapwsWaterRefraction::prepare(589.0).unwrap();
    for temperature in [261.14, 773.16, 0.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            prepared.sample(temperature, 1000.0),
            Err(WaterRefractionError::TemperatureOutsideRange)
        );
    }
    for density in [-1.0, 1060.01, f64::NAN, f64::INFINITY] {
        assert_eq!(
            prepared.sample(300.0, density),
            Err(WaterRefractionError::DensityOutsideRange)
        );
    }
    let vacuum = prepared.sample(300.0, 0.0).unwrap();
    assert_eq!(vacuum.phase_index, 1.0);
    assert_eq!(vacuum.phase_index_minus_one, 0.0);
    assert_eq!(vacuum.group_index, 1.0);
    assert_eq!(vacuum.wavelength_derivative_per_nm, 0.0);
    assert_eq!(vacuum.temperature_derivative_per_k, 0.0);
    assert_eq!(vacuum.temperature_second_derivative_per_k2, 0.0);
    let dilute = prepared.sample(300.0, 1.0e-15).unwrap();
    assert_eq!(dilute.phase_index, 1.0);
    assert!(dilute.phase_index_minus_one > 0.0);
    near(
        dilute.phase_index_minus_one / 1.0e-15,
        vacuum.density_derivative_m3_per_kg,
        1.0e-18,
    );
}
