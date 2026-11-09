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

// Independent NIST A32--A41 density-ratio evaluation, retaining the original
// molar-mass/density equations instead of the production rational reduction.
fn original_refractivity([temperature, pressure, water, co2]: [f64; 4], wavelength: f64) -> f64 {
    let t = temperature - 273.15;
    let sigma2 = (1000.0 / wavelength).powi(2);
    let dry = 1e-8
        * (5_792_105.0 / (238.0185 - sigma2) + 167_917.0 / (57.362 - sigma2))
        * (1.0 + 5.34e-7 * (co2 - 450.0));
    let vapor = 1.022e-8
        * (295.235 + 2.6422 * sigma2 - 0.03238 * sigma2.powi(2) + 0.004028 * sigma2.powi(3));
    let z = 1.0
        - pressure / temperature
            * (1.58123e-6 - 2.9331e-8 * t
                + 1.1043e-10 * t * t
                + (5.707e-6 - 2.051e-8 * t) * water
                + (1.9898e-4 - 2.376e-6 * t) * water * water)
        + (pressure / temperature).powi(2) * (1.83e-11 - 0.765e-8 * water * water);
    let molar_mass = 0.0289635 + 1.2011e-8 * (co2 - 400.0);
    let dry_density = pressure * molar_mass * (1.0 - water) / (z * 8.314472 * temperature);
    let dry_reference = 101_325.0 * molar_mass / (0.9995922115 * 8.314472 * 288.15);
    let vapor_density = pressure * 0.018015 * water / (z * 8.314472 * temperature);
    dry * dry_density / dry_reference + vapor * vapor_density / 0.00985938
}

#[test]
fn g1_full_state_gradient_and_all_mixed_partials_match_original_density_law() {
    for state in [
        [293.15, 101_325.0, 0.01, 450.0],
        [313.15, 85_000.0, 0.025, 700.0],
        [353.15, 110_000.0, 0.09, 1750.0],
    ] {
        for wavelength in [300.0, 633.0, 1700.0] {
            let composition = CiddorAir::try_new(state[2], state[3]).unwrap();
            let spectrum = PreparedCiddorAirSpectrum::try_new(wavelength).unwrap();
            let jet = spectrum
                .sample_jet(state[0], state[1], composition)
                .unwrap();
            let value = |s| original_refractivity(s, wavelength);
            let center = value(state);
            close(jet.sample.phase_index_minus_one, center, 3e-19);
            for (i, (h, tolerance)) in [(0.003, 5e-15), (2.0, 3e-19), (2e-5, 3e-13), (5.0, 2e-20)]
                .into_iter()
                .enumerate()
            {
                let (mut plus, mut minus) = (state, state);
                plus[i] += h;
                minus[i] -= h;
                close(
                    jet.phase_gradient[i],
                    (value(plus) - value(minus)) / (2.0 * h),
                    tolerance,
                );
            }
            let steps = [0.02, 100.0, 0.0005, 10.0];
            for (i, tolerance) in [2e-15, 5e-23, 3e-11, 3e-21].into_iter().enumerate() {
                let (mut plus, mut minus) = (state, state);
                plus[i] += steps[i];
                minus[i] -= steps[i];
                close(
                    jet.phase_hessian_diagonal[i],
                    (value(plus) - 2.0 * center + value(minus)) / steps[i].powi(2),
                    tolerance,
                );
            }
            let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
            let tolerances = [3e-19, 3e-13, 3e-18, 3e-17, 3e-22, 3e-16];
            for (k, (i, j)) in pairs.into_iter().enumerate() {
                let (mut pp, mut pm, mut mp, mut mm) = (state, state, state, state);
                pp[i] += steps[i];
                pp[j] += steps[j];
                pm[i] += steps[i];
                pm[j] -= steps[j];
                mp[i] -= steps[i];
                mp[j] += steps[j];
                mm[i] -= steps[i];
                mm[j] -= steps[j];
                let reference =
                    (value(pp) - value(pm) - value(mp) + value(mm)) / (4.0 * steps[i] * steps[j]);
                close(jet.phase_hessian_off_diagonal[k], reference, tolerances[k]);
            }
            assert!(jet.phase_gradient[2] < 0.0);
            assert!(jet.phase_gradient[3] > 0.0);
            assert_eq!(jet.phase_hessian_diagonal[3], 0.0); // Exactly affine in dry CO2.
        }
    }
}

#[test]
fn g3_full_state_preparation_preserves_frozen_owner_and_admitted_boundaries() {
    for wavelength in WAVELENGTH_RANGE_NM {
        let spectrum = PreparedCiddorAirSpectrum::try_new(wavelength).unwrap();
        assert_eq!(spectrum.wavelength_nm(), wavelength);
        for water in [0.0, 0.2] {
            for co2 in [0.0, 2000.0] {
                let composition = CiddorAir::try_new(water, co2).unwrap();
                let frozen = composition.prepare(wavelength).unwrap();
                assert_eq!(spectrum.bind(composition), frozen);
                // Numerical-domain corners do not establish vapor phase stability.
                for t in TEMPERATURE_RANGE_K {
                    for p in PRESSURE_RANGE_PA {
                        let jet = spectrum.sample_jet(t, p, composition).unwrap();
                        assert_eq!(jet.sample, frozen.sample(t, p).unwrap());
                        assert!(
                            jet.phase_gradient
                                .into_iter()
                                .chain(jet.phase_hessian_diagonal)
                                .chain(jet.phase_hessian_off_diagonal)
                                .all(f64::is_finite)
                        );
                    }
                }
            }
        }
    }
    for wavelength in [299.0, 1701.0, f64::NAN] {
        assert!(PreparedCiddorAirSpectrum::try_new(wavelength).is_err());
    }
    let spectrum = PreparedCiddorAirSpectrum::try_new(633.0).unwrap();
    let composition = CiddorAir::try_new(0.01, 450.0).unwrap();
    for (t, p) in [
        (233.14, 101_325.0),
        (293.15, 120_001.0),
        (f64::NAN, 101_325.0),
    ] {
        assert!(spectrum.sample_jet(t, p, composition).is_err());
    }
}

#[test]
fn g5_frozen_ciddor_sample_bits_survive_full_state_extension() {
    // Captured from the original compiled owner before the spectral/state split.
    // The fixed-state lane must not silently change its scalar operation order.
    let states = [
        (0.0, 450.0, 293.15, 101_325.0, 633.0),
        (0.025, 700.0, 313.15, 85_000.0, 410.0),
        (0.09, 1750.0, 353.15, 110_000.0, 1550.0),
    ];
    let expected = [
        [
            0x3ff0011d00b634ca,
            0x3f31d00b634ca7c0,
            0x3ff001253e6e5458,
            0xbe4aa9c21e190737,
            0x3dd0aa5226693ace,
            0xbeaf33441ded63b3,
            0x3e270ce022b5152a,
            0x3e3b69f5ad06d37a,
            0x3c75886443bb4daa,
            0xbda43ebc5ec49fb3,
        ],
        [
            0x3ff000e3b10cd392,
            0x3f2c76219a723ba8,
            0x3ff000f44381004c,
            0xbe64b1e4beb86748,
            0x3df4ee9ad18ac4bb,
            0xbea751096a567f9e,
            0x3e25f2dbddf95f8e,
            0x3e33285750a7fbe2,
            0x3c6ca25cac2caf69,
            0xbda204e1a8a410a2,
        ],
        [
            0x3ff000faa11b907b,
            0x3f2f5423720f688a,
            0x3ff000fbdc217a04,
            0xbe0a03ca9d448244,
            0x3d79df4318f61449,
            0xbea6cb851575e9f7,
            0x3e22aadbae998283,
            0x3e30a8128e57f27b,
            0x3c57dc6f6e226c28,
            0xbd9b442d2f48b29e,
        ],
    ];
    for ((x, co2, t, p, wavelength), bits) in states.into_iter().zip(expected) {
        let m = CiddorAir::try_new(x, co2)
            .unwrap()
            .prepare(wavelength)
            .unwrap()
            .sample(t, p)
            .unwrap();
        let actual = [
            m.phase_index,
            m.phase_index_minus_one,
            m.group_index,
            m.wavelength_derivative_per_nm,
            m.wavelength_second_derivative_per_nm2,
            m.temperature_derivative_per_k,
            m.pressure_derivative_per_pa,
            m.temperature_second_derivative_per_k2,
            m.pressure_second_derivative_per_pa2,
            m.temperature_pressure_derivative_per_k_pa,
        ];
        assert_eq!(actual.map(f64::to_bits), bits);
    }
}
