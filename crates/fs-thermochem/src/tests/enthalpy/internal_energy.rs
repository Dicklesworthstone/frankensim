use super::*;

#[test]
fn g0_internal_energy_change_retains_small_cv_and_formation_invariance() {
    for cp in [1.0_f64.next_up(), 2.5, 3.5] {
        for reference in [-1e24, 0.0, 1e24] {
            let mut coefficients = [0.0; 9];
            coefficients[2] = cp;
            coefficients[7] = reference;
            let source = model(coefficients);
            for (a, b) in [(500.0_f64, 500.0_f64.next_up()), (500.0, 1800.0)] {
                let actual = source.molar_internal_energy_change(
                    Temperature::new(a), Temperature::new(b),
                ).unwrap().value();
                let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * ((cp - 1.0) * (b - a));
                assert!(actual > 0.0);
                assert!((actual / expected - 1.0).abs() < 5e-15);
                assert_eq!(source.molar_internal_energy_change(
                    Temperature::new(b), Temperature::new(a),
                ).unwrap().value(), -actual);
                assert_eq!(source.molar_internal_energy_change(
                    Temperature::new(a), Temperature::new(a),
                ).unwrap().value(), 0.0);
            }
        }
    }
}

#[test]
fn g1_mixture_internal_energy_change_matches_variable_cv_and_mass_basis() {
    let ordinary = mixture(0.0);
    let shifted = mixture(1e24);
    for (a, b) in [(500.0_f64, 500.0_f64.next_up()), (500.0, 1800.0)] {
        let evaluate = |source: &FrozenIdealGasMixtureModelV1| source
            .specific_internal_energy_change(Temperature::new(a), Temperature::new(b))
            .unwrap().value();
        let actual = evaluate(&ordinary);
        assert_eq!(actual, evaluate(&shifted));
        let molar = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (
            (0.75 * 2.5 + 0.25 * 1.5) * (b - a) + 0.75 * 0.0005 * (b - a) * (b + a));
        let mass = 0.75 * 0.028_013_4 + 0.25 * 0.004;
        assert!((actual / (molar / mass) - 1.0).abs() < 5e-14);
        assert!((ordinary.molar_internal_energy_change(Temperature::new(a), Temperature::new(b))
            .unwrap().value() / molar - 1.0).abs() < 5e-14);
        let dh = ordinary.molar_enthalpy_change(Temperature::new(a), Temperature::new(b))
            .unwrap().value();
        assert!(((dh - molar) / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (b - a)) - 1.0).abs() < 1e-13);
    }
}

#[test]
fn g3_internal_energy_change_respects_raw_jumps_and_continuous_caloric_seams() {
    let mut low = [0.0; 9];
    low[2] = 3.5;
    low[7] = 1e12;
    let mut high = low;
    high[2] = 4.0;
    high[3] = 0.001;
    high[7] += 0.25;
    let raw = Nasa9StandardStateModelV1::new(
        SpeciesId::new("N2").unwrap(), MolarMass::new(0.028), convention(100_000.0),
        vec![
            Nasa9RegionV1::from_card(card(200.0, 600.0, 100_000.0, low)).unwrap(),
            Nasa9RegionV1::from_card(card(600.0, 1500.0, 100_000.0, high)).unwrap(),
        ],
    ).unwrap();
    let du = raw.molar_internal_energy_change(Temperature::new(500.0), Temperature::new(600.0))
        .unwrap().value();
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (3.0 * 600.0 + 0.0005 * 600.0 * 600.0 - 2.5 * 500.0 + 0.25);
    assert!((du / expected - 1.0).abs() < 5e-14);
    let continuous = raw.with_caloric_continuity(
        Temperature::new(500.0), MolarEnergyQuantityV1::new(1e5),
        MolarThermalQuantityV1::new(1e3),
    ).unwrap();
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (2.5 * 100.0 + 3.0 * 300.0 + 0.0005 * (900.0 - 600.0) * (900.0 + 600.0));
    let actual = continuous.molar_internal_energy_change(
        Temperature::new(500.0), Temperature::new(900.0),
    ).unwrap().value();
    assert!((actual / expected - 1.0).abs() < 5e-14);
    assert_eq!(continuous.molar_internal_energy_change(
        Temperature::new(900.0), Temperature::new(500.0),
    ).unwrap().value(), -actual);
}

#[test]
fn g0_internal_energy_change_refuses_invalid_endpoints_and_overflow() {
    let source = mixture(0.0);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 199.0, 6001.0] {
        for (a, b) in [(bad, 500.0), (500.0, bad), (bad, bad)] {
            assert!(source.specific_internal_energy_change(
                Temperature::new(a), Temperature::new(b),
            ).is_err());
        }
    }
    let mut coefficients = [0.0; 9];
    coefficients[6] = f64::MAX;
    assert!(matches!(model(coefficients).molar_internal_energy_change(
        Temperature::new(500.0), Temperature::new(900.0),
    ), Err(ThermochemErrorV1::NonFiniteEvaluation {
        property: Nasa9PropertyV1::InternalEnergy, ..
    })));
}
