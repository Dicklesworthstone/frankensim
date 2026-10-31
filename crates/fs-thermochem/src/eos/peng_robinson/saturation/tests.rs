use super::*;
use fs_qty::MolarMass;

fn methane() -> PengRobinsonComponentV1 {
    PengRobinsonComponentV1 {
        critical_temperature: Temperature::new(190.564),
        critical_pressure: Pressure::new(4.5992e6),
        acentric_factor: 0.01142,
        molar_mass: MolarMass::new(0.016_042_46),
    }
}
fn pure(component: PengRobinsonComponentV1) -> PengRobinsonMixtureV1 {
    PengRobinsonMixtureV1::new(vec![component], vec![1.0], vec![]).unwrap()
}
fn solve(component: PengRobinsonComponentV1, temperature: f64) -> PengRobinsonSaturationV1 {
    pure(component)
        .pure_saturation_at_temperature(
            Temperature::new(temperature),
            [
                Pressure::new(1.0e5),
                Pressure::new(0.98 * component.critical_pressure.value()),
            ],
            PengRobinsonSaturationOptionsV1::default(),
        )
        .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{actual} != {expected} (relative tolerance {relative})"
    );
}

#[test]
fn g1_saturation_matches_independent_maxwell_equal_area_solutions() {
    // Independent NumPy cubic roots + SciPy brentq on the Maxwell-area
    // integral, not on the production fugacity equation or its Newton step.
    for (t, p, vl, vv) in [
        (
            150.0,
            1_047_062.532_316_561_1,
            4.128_388_060_116_956e-5,
            9.710_793_482_353_756e-4,
        ),
        (
            180.0,
            3_309_034.443_368_976_5,
            5.962_818_439_365_452e-5,
            2.505_884_024_879_326e-4,
        ),
    ] {
        let component = methane();
        let saturation = solve(component, t);
        close(saturation.pressure().value(), p, 2.0e-9);
        close(saturation.liquid().molar_volume_m3_per_mol(), vl, 2.0e-9);
        close(saturation.vapor().molar_volume_m3_per_mol(), vv, 2.0e-9);
        assert_eq!(saturation.component(), &component);
        assert!(saturation.log_fugacity_residual().abs() <= 1e-10);
        assert!(saturation.relative_pressure_correction() <= 1e-10);
        assert!(saturation.liquid().pressure_derivatives().1 < 0.0);
        assert!(saturation.vapor().pressure_derivatives().1 < 0.0);
        assert!((2..=128).contains(&saturation.evaluations()));

        let tc = component.critical_temperature.value();
        let pc = component.critical_pressure.value();
        let omega = component.acentric_factor;
        let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        let alpha_base = 1.0 + kappa * (1.0 - (t / tc).sqrt());
        let a = 0.45724 * R * R * tc * tc / pc * alpha_base * alpha_base;
        let b = 0.07780 * R * tc / pc;
        let l = saturation.liquid().molar_volume_m3_per_mol();
        let v = saturation.vapor().molar_volume_m3_per_mol();
        let sqrt2 = 2.0_f64.sqrt();
        let ratio = |volume: f64| (volume + (1.0 - sqrt2) * b) / (volume + (1.0 + sqrt2) * b);
        let area = R * t * ((v - b) / (l - b)).ln()
            - a / (2.0 * sqrt2 * b) * (ratio(v) / ratio(l)).ln()
            - saturation.pressure().value() * (v - l);
        assert!(area.abs() < 2.0e-6, "Maxwell area = {area} J/mol");
    }
}

#[test]
fn g1_bracket_may_start_on_two_different_single_phase_branches() {
    let model = pure(methane());
    let temperature = Temperature::new(180.0);
    let bracket = [Pressure::new(1.0e5), Pressure::new(4.4e6)];
    for pressure in bracket {
        assert_eq!(
            model
                .evaluate_pt(temperature, pressure)
                .unwrap()
                .phase_roots()
                .len(),
            1
        );
    }
    let saturation = model
        .pure_saturation_at_temperature(
            temperature,
            bracket,
            PengRobinsonSaturationOptionsV1::default(),
        )
        .unwrap();
    close(
        saturation.pressure().value(),
        3_309_034.443_368_976_5,
        2.0e-9,
    );
    assert!(
        saturation.liquid().molar_volume_m3_per_mol()
            < saturation.vapor().molar_volume_m3_per_mol()
    );
}

#[test]
fn g1_clapeyron_and_latent_energy_use_the_original_departure_functions() {
    for (t, latent_h, latent_u, slope) in [
        (
            150.0,
            6621.875_622_767_658,
            5648.321_625_890_125,
            47_479.084_402_770_27,
        ),
        (
            180.0,
            3717.546_730_836_839,
            3085.652_791_849_621,
            108_153.612_305_477_6,
        ),
    ] {
        let saturation = solve(methane(), t);
        close(saturation.latent_enthalpy_j_mol(), latent_h, 2e-9);
        close(saturation.latent_internal_energy_j_mol(), latent_u, 2e-9);
        close(saturation.saturation_pressure_slope_pa_k(), slope, 2e-9);
        let step = 0.001;
        let differentiated = (solve(methane(), t + step).pressure().value()
            - solve(methane(), t - step).pressure().value())
            / (2.0 * step);
        close(
            saturation.saturation_pressure_slope_pa_k(),
            differentiated,
            2e-7,
        );
        let l = saturation.at_vapor_fraction(0.0).unwrap();
        let v = saturation.at_vapor_fraction(1.0).unwrap();
        close(
            t * (v.residual_entropy_j_mol_k() - l.residual_entropy_j_mol_k()),
            latent_h,
            2e-9,
        );
        close(
            saturation.latent_enthalpy_j_mol() - saturation.latent_internal_energy_j_mol(),
            saturation.pressure().value() * (v.molar_volume().value() - l.molar_volume().value()),
            2e-12,
        );
    }
}

#[test]
fn g0_volume_constraint_conserves_phase_amount_volume_and_energy() {
    let saturation = solve(methane(), 150.0);
    let liquid = saturation.liquid();
    let vapor = saturation.vapor();
    let vl = liquid.molar_volume_m3_per_mol();
    let vv = vapor.molar_volume_m3_per_mol();
    for quality in [0.0, 0.1, 0.37, 0.5, 0.94, 1.0] {
        let volume = (1.0 - quality) * vl + quality * vv;
        let mixture = saturation
            .at_molar_volume(MolarVolume::new(volume))
            .unwrap();
        assert!((mixture.vapor_fraction() - quality).abs() < 4e-16);
        close(mixture.molar_volume().value(), volume, 4e-16);
        close(
            mixture.mass_density_kg_m3() * volume,
            methane().molar_mass.value(),
            4e-16,
        );
        let weighted = |l: f64, v: f64| (1.0 - quality) * l + quality * v;
        close(
            mixture.residual_enthalpy_j_mol(),
            weighted(
                liquid.residual_enthalpy_j_mol(),
                vapor.residual_enthalpy_j_mol(),
            ),
            2e-15,
        );
        close(
            mixture.residual_internal_energy_j_mol(),
            weighted(
                liquid.residual_internal_energy_j_mol(),
                vapor.residual_internal_energy_j_mol(),
            ),
            2e-15,
        );
        // Adding the ideal reference once makes h-u=pv; residuals alone
        // necessarily differ by pv-RT, with no separate latent source.
        close(
            mixture.residual_enthalpy_j_mol() - mixture.residual_internal_energy_j_mol()
                + R * 150.0,
            saturation.pressure().value() * volume,
            2e-12,
        );
    }
    for invalid in [0.0, 0.99 * vl, 1.01 * vv, f64::NAN, f64::INFINITY] {
        assert!(
            saturation
                .at_molar_volume(MolarVolume::new(invalid))
                .is_err()
        );
    }
    for invalid in [-0.1, 1.1, f64::NAN] {
        assert!(saturation.at_vapor_fraction(invalid).is_err());
    }
}

#[test]
fn g3_another_original_component_changes_saturation_and_latent_heat() {
    let ethane = PengRobinsonComponentV1 {
        critical_temperature: Temperature::new(305.322),
        critical_pressure: Pressure::new(4.8722e6),
        acentric_factor: 0.0995,
        molar_mass: MolarMass::new(0.030_069_04),
    };
    let saturation = solve(ethane, 240.0);
    close(saturation.pressure().value(), 967_970.628_563_914_6, 2e-9);
    close(
        saturation.liquid().molar_volume_m3_per_mol(),
        6.093_232_803_685_145e-5,
        2e-9,
    );
    close(
        saturation.vapor().molar_volume_m3_per_mol(),
        0.001_721_935_656_465_691,
        2e-9,
    );
    close(
        saturation.latent_enthalpy_j_mol(),
        11_967.935_697_895_206,
        2e-9,
    );
    close(
        saturation.latent_internal_energy_j_mol(),
        10_360.133_262_029_187,
        2e-9,
    );
    close(
        saturation.saturation_pressure_slope_pa_k(),
        30_021.853_591_587_427,
        2e-9,
    );
    assert_eq!(saturation.component(), &ethane);
}

#[test]
fn g4_critical_domain_cancellation_and_whole_call_work_limits_refuse() {
    let model = pure(methane());
    let bracket = [Pressure::new(1e5), Pressure::new(4.4e6)];
    let options = PengRobinsonSaturationOptionsV1::default();
    let evaluate = |t, bracket, options| {
        model.pure_saturation_at_temperature(Temperature::new(t), bracket, options)
    };
    for t in [0.0, -10.0, 190.564, 200.0, f64::NAN] {
        assert!(evaluate(t, bracket, options).is_err());
    }
    // Rounded original PR constants put its numerical critical point slightly
    // below the supplied Tc. Do not fabricate two indistinguishable phases.
    assert!(matches!(
        evaluate(190.564 * (1.0 - 1e-12), bracket, options),
        Err(E::CriticalRegionUnresolved)
    ));
    for invalid in [
        [Pressure::new(0.0), bracket[1]],
        [bracket[1], bracket[0]],
        [bracket[0], Pressure::new(f64::NAN)],
    ] {
        assert!(evaluate(150.0, invalid, options).is_err());
    }
    assert!(matches!(
        evaluate(150.0, [Pressure::new(1.2e6), Pressure::new(2e6)], options),
        Err(E::NotBracketed)
    ));
    assert!(matches!(
        evaluate(
            150.0,
            bracket,
            PengRobinsonSaturationOptionsV1 {
                maximum_evaluations: 2,
                ..options
            }
        ),
        Err(E::NoConvergence { evaluations: 2, .. })
    ));
    assert!(
        evaluate(
            150.0,
            bracket,
            PengRobinsonSaturationOptionsV1 {
                maximum_evaluations: 513,
                ..options
            }
        )
        .is_err()
    );
    assert!(
        evaluate(
            150.0,
            bracket,
            PengRobinsonSaturationOptionsV1 {
                log_fugacity_tolerance: f64::NAN,
                ..options
            }
        )
        .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        model.pure_saturation_at_temperature_cancellable(
            Temperature::new(150.0),
            bracket,
            options,
            || {
                polls += 1;
                polls >= 4
            },
        ),
        Err(E::Cancelled)
    );
    assert_eq!(polls, 4);
    let mixture =
        PengRobinsonMixtureV1::new(vec![methane(), methane()], vec![0.5, 0.5], vec![]).unwrap();
    assert!(matches!(
        mixture.pure_saturation_at_temperature(Temperature::new(150.0), bracket, options),
        Err(E::Invalid(_))
    ));
}
