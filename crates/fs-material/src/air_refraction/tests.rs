use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g1_original_nist_dry_and_moist_air_table() {
    // NIST documentation Tables 1/2 print phase index to nine decimal places.
    let dry = CiddorAir::try_new(0.0, 450.0).unwrap();
    for (temperature_k, p, wavelength, n) in [
        (293.15, 101_325.0, 633.0, 1.000271800),
        (293.15, 60_000.0, 633.0, 1.000160924),
        (293.15, 120_000.0, 633.0, 1.000321916),
        (323.15, 100_000.0, 633.0, 1.000243285),
        (278.15, 100_000.0, 633.0, 1.000282756),
        (233.15, 100_000.0, 633.0, 1.000337580),
        (293.15, 101_325.0, 1700.0, 1.000268479),
        (293.15, 101_325.0, 300.0, 1.000286581),
        (233.15, 120_000.0, 300.0, 1.000427233),
    ] {
        close(
            dry.prepare(wavelength)
                .unwrap()
                .sample(temperature_k, p)
                .unwrap()
                .phase_index,
            n,
            5.1e-10,
        );
    }
    // Independent A15--A18 water fraction from RH and table saturation pressure.
    // Printed saturation pressures are rounded to 1 Pa; allow that rounding too.
    for (t, p, wavelength, rh, saturation_pa, n) in [
        (50.0, 120_000.0, 633.0, 1.0, 12_351.0, 1.000287924),
        (40.0, 120_000.0, 633.0, 0.75, 7384.0, 1.000299418),
        (20.0, 100_000.0, 633.0, 1.0, 2339.0, 1.000267394),
        (40.0, 110_000.0, 1700.0, 1.0, 7384.0, 1.000270247),
        (40.0, 110_000.0, 300.0, 1.0, 7384.0, 1.000289000),
    ] {
        let x = (1.00062 + 3.14e-8 * p + 5.6e-7 * t * t) * rh * saturation_pa / p;
        let air = CiddorAir::try_new(x, 450.0).unwrap();
        close(
            air.prepare(wavelength)
                .unwrap()
                .sample(t + 273.15, p)
                .unwrap()
                .phase_index,
            n,
            7.5e-10,
        );
    }
}

#[test]
fn g1_pressure_temperature_and_spectral_jets_match_independent_differences() {
    for (x, co2, t, p, wavelength) in [
        (0.0, 450.0, 293.15, 101_325.0, 633.0),
        (0.025, 700.0, 313.15, 85_000.0, 410.0),
        (0.09, 2000.0, 353.15, 110_000.0, 1550.0),
    ] {
        let air = CiddorAir::try_new(x, co2).unwrap();
        let prepared = air.prepare(wavelength).unwrap();
        let m = prepared.sample(t, p).unwrap();
        let ht = 0.01;
        let tm = prepared.sample(t - ht, p).unwrap();
        let tp = prepared.sample(t + ht, p).unwrap();
        close(
            m.temperature_derivative_per_k,
            (tp.phase_index_minus_one - tm.phase_index_minus_one) / (2.0 * ht),
            3e-15,
        );
        close(
            m.temperature_second_derivative_per_k2,
            (tp.temperature_derivative_per_k - tm.temperature_derivative_per_k) / (2.0 * ht),
            3e-17,
        );
        let hp = 1.0;
        let pm = prepared.sample(t, p - hp).unwrap();
        let pp = prepared.sample(t, p + hp).unwrap();
        close(
            m.pressure_derivative_per_pa,
            (pp.phase_index_minus_one - pm.phase_index_minus_one) / (2.0 * hp),
            2e-19,
        );
        close(
            m.pressure_second_derivative_per_pa2,
            (pp.pressure_derivative_per_pa - pm.pressure_derivative_per_pa) / (2.0 * hp),
            2e-24,
        );
        close(
            m.temperature_pressure_derivative_per_k_pa,
            (pp.temperature_derivative_per_k - pm.temperature_derivative_per_k) / (2.0 * hp),
            2e-22,
        );
        let hl = 0.001;
        let lm = air.prepare(wavelength - hl).unwrap().sample(t, p).unwrap();
        let lp = air.prepare(wavelength + hl).unwrap().sample(t, p).unwrap();
        close(
            m.wavelength_derivative_per_nm,
            (lp.phase_index_minus_one - lm.phase_index_minus_one) / (2.0 * hl),
            1e-16,
        );
        close(
            m.wavelength_second_derivative_per_nm2,
            (lp.wavelength_derivative_per_nm - lm.wavelength_derivative_per_nm) / (2.0 * hl),
            1e-18,
        );
        close(
            m.group_index,
            m.phase_index - wavelength * m.wavelength_derivative_per_nm,
            3e-16,
        );
        assert!(m.group_index > m.phase_index);
    }
}

#[test]
fn g3_air_relative_conversion_keeps_the_group_product_rule() {
    let air = CiddorAir::try_new(0.012, 450.0).unwrap();
    let wavelength = 550.0;
    let relative = |lambda: f64| 1.5 + 0.004 / (lambda * 0.001).powi(2);
    let relative_group = 1.5 + 0.012 / (wavelength * 0.001_f64).powi(2);
    let sample = air
        .prepare(wavelength)
        .unwrap()
        .sample(293.15, 101_325.0)
        .unwrap();
    let (phase, group) = sample
        .absolute_from_relative(relative(wavelength), Some(relative_group))
        .unwrap();
    close(phase, sample.phase_index * relative(wavelength), 0.0);
    let absolute = |lambda| {
        relative(lambda)
            * air
                .prepare(lambda)
                .unwrap()
                .sample(293.15, 101_325.0)
                .unwrap()
                .phase_index
    };
    let h = 0.01;
    close(
        group.unwrap(),
        phase - wavelength * (absolute(wavelength + h) - absolute(wavelength - h)) / (2.0 * h),
        5e-11,
    );
    assert!(phase > relative(wavelength) + 4e-4);
    assert_eq!(
        sample
            .absolute_from_relative(relative(wavelength), None)
            .unwrap()
            .1,
        None
    );
    for value in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        assert!(sample.absolute_from_relative(value, None).is_err());
        assert!(sample.absolute_from_relative(1.5, Some(value)).is_err());
    }
}

#[test]
fn g0_original_domains_and_composition_identity_remain_explicit() {
    let dry = CiddorAir::try_new(0.0, 450.0).unwrap();
    let negative_zero = CiddorAir::try_new(-0.0, 450.0).unwrap();
    assert_eq!(dry, negative_zero);
    assert_eq!(dry.identity(), negative_zero.identity());
    assert_ne!(
        dry.identity(),
        CiddorAir::try_new(0.0, 451.0).unwrap().identity()
    );
    assert_ne!(
        dry.identity(),
        CiddorAir::try_new(0.001, 450.0).unwrap().identity()
    );
    for x in [f64::NAN, f64::INFINITY, -0.001, 0.201] {
        assert!(CiddorAir::try_new(x, 450.0).is_err());
    }
    for co2 in [f64::NAN, f64::INFINITY, -1.0, 2001.0] {
        assert!(CiddorAir::try_new(0.0, co2).is_err());
    }
    for wavelength in [299.0, 1701.0, f64::NAN] {
        assert!(dry.prepare(wavelength).is_err());
    }
    let prepared = dry.prepare(550.0).unwrap();
    for t in [233.14, 373.16, f64::NAN] {
        assert!(prepared.sample(t, 101_325.0).is_err());
    }
    for p in [59_999.0, 120_001.0, f64::NAN] {
        assert!(prepared.sample(293.15, p).is_err());
    }
    // These are numerical-domain corners, not assertions of vapor-only phase stability.
    for x in [0.0, 0.2] {
        for co2 in [0.0, 2000.0] {
            for wavelength in WAVELENGTH_RANGE_NM {
                let model = CiddorAir::try_new(x, co2)
                    .unwrap()
                    .prepare(wavelength)
                    .unwrap();
                for t in TEMPERATURE_RANGE_K {
                    for p in PRESSURE_RANGE_PA {
                        let m = model.sample(t, p).unwrap();
                        assert!(
                            [
                                m.phase_index,
                                m.phase_index_minus_one,
                                m.group_index,
                                m.wavelength_derivative_per_nm,
                                m.wavelength_second_derivative_per_nm2,
                                m.temperature_derivative_per_k,
                                m.pressure_derivative_per_pa,
                                m.temperature_second_derivative_per_k2,
                                m.pressure_second_derivative_per_pa2,
                                m.temperature_pressure_derivative_per_k_pa
                            ]
                            .into_iter()
                            .all(f64::is_finite)
                        );
                    }
                }
            }
        }
    }
}
