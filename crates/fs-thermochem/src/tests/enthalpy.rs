use super::*;
use crate::mixture::FrozenIdealGasMixtureModelV1;

mod internal_energy;

#[test]
fn g0_enthalpy_change_retains_small_sensible_energy_under_large_formation_reference() {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[7] = 1e24;
    let source = model(coefficients);
    let (a, b) = (500.0, 500.0 + 1e-6);
    let absolute = |t| {
        source
            .evaluate(Temperature::new(t))
            .unwrap()
            .properties()
            .h()
            .value()
    };
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 3.5 * (b - a);
    assert!(
        ((absolute(b) - absolute(a)) - expected).abs() > 0.9 * expected,
        "absolute subtraction loses the sensible energy in reference-scale rounding"
    );
    let change = source
        .molar_enthalpy_change(Temperature::new(a), Temperature::new(b))
        .unwrap()
        .value();
    assert!((change / expected - 1.0).abs() < 5e-15);
    assert_eq!(
        source
            .molar_enthalpy_change(Temperature::new(b), Temperature::new(a))
            .unwrap()
            .value(),
        -change
    );
    assert_eq!(
        source
            .molar_enthalpy_change(Temperature::new(a), Temperature::new(a))
            .unwrap()
            .value(),
        0.0
    );
}

#[test]
fn g0_enthalpy_change_matches_each_integrated_nasa9_heat_capacity_term() {
    let (a, b) = (400.0_f64, 900.0_f64);
    let expected = [
        1.0 / a - 1.0 / b,
        (b / a).ln(),
        b - a,
        (b.powi(2) - a.powi(2)) / 2.0,
        (b.powi(3) - a.powi(3)) / 3.0,
        (b.powi(4) - a.powi(4)) / 4.0,
        (b.powi(5) - a.powi(5)) / 5.0,
        0.0,
        0.0,
    ];
    for channel in 0..9 {
        let mut coefficients = [0.0; 9];
        coefficients[channel] = 1.0;
        let actual = model(coefficients)
            .molar_enthalpy_change(Temperature::new(a), Temperature::new(b))
            .unwrap()
            .value()
            / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        assert!(
            (actual - expected[channel]).abs() <= 3e-14 * expected[channel].abs().max(1e-30),
            "channel {channel}: {actual} != {}",
            expected[channel]
        );
    }
}

#[test]
fn g0_enthalpy_change_preserves_distinct_region_coefficients_and_reference_jump() {
    let mut low = [0.0; 9];
    low[2] = 3.5;
    low[7] = 1e12;
    let mut high = low;
    high[2] = 4.0;
    high[3] = 0.001;
    high[7] += 0.25;
    let source = Nasa9StandardStateModelV1::new(
        SpeciesId::new("N2").unwrap(),
        MolarMass::new(0.028),
        convention(100_000.0),
        vec![
            Nasa9RegionV1::from_card(card(200.0, 600.0, 100_000.0, low)).unwrap(),
            Nasa9RegionV1::from_card(card(600.0, 1500.0, 100_000.0, high)).unwrap(),
        ],
    )
    .unwrap();
    // The selected upper region contributes both its new cp and its reference
    // discontinuity. Neither is replaced by integrating a smoothed source.
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (4.0 * 600.0 + 0.0005 * 600.0 * 600.0 - 3.5 * 500.0 + 0.25);
    let actual = source
        .molar_enthalpy_change(Temperature::new(500.0), Temperature::new(600.0))
        .unwrap()
        .value();
    assert!((actual / expected - 1.0).abs() < 5e-14);
    assert_eq!(
        source
            .molar_enthalpy_change(Temperature::new(600.0), Temperature::new(500.0))
            .unwrap()
            .value(),
        -actual
    );
}

fn mixture(reference: f64) -> FrozenIdealGasMixtureModelV1 {
    let mut a = [0.0; 9];
    a[2] = 3.5;
    a[3] = 0.001;
    a[7] = reference;
    let mut b = [0.0; 9];
    b[2] = 2.5;
    b[7] = -reference;
    let second = Nasa9StandardStateModelV1::new(
        SpeciesId::new("He").unwrap(),
        MolarMass::new(0.004),
        convention(100_000.0),
        vec![Nasa9RegionV1::from_card(card(200.0, 6000.0, 100_000.0, b)).unwrap()],
    )
    .unwrap();
    FrozenIdealGasMixtureModelV1::new(
        vec![model(a), second],
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    )
    .unwrap()
}

#[test]
fn g1_mixture_enthalpy_change_uses_composition_and_molar_mass_without_formation_offsets() {
    let ordinary = mixture(0.0);
    let shifted = mixture(1e24);
    for (a, b) in [(500.0, 500.0 + 1e-6), (500.0, 1800.0)] {
        let mass_specific = |source: &FrozenIdealGasMixtureModelV1| {
            source
                .specific_enthalpy_change(Temperature::new(a), Temperature::new(b))
                .unwrap()
                .value()
        };
        assert_eq!(mass_specific(&ordinary), mass_specific(&shifted));
        let enthalpy = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            * ((0.75 * 3.5 + 0.25 * 2.5) * (b - a) + 0.75 * 0.0005 * (b - a) * (b + a));
        let expected = enthalpy / (0.75 * 0.028_013_4 + 0.25 * 0.004);
        assert!((mass_specific(&ordinary) / expected - 1.0).abs() < 5e-14);
        let molar = ordinary
            .molar_enthalpy_change(Temperature::new(a), Temperature::new(b))
            .unwrap()
            .value();
        assert!((molar / enthalpy - 1.0).abs() < 5e-14);
        assert_eq!(
            ordinary
                .specific_enthalpy_change(Temperature::new(b), Temperature::new(a))
                .unwrap()
                .value(),
            -mass_specific(&ordinary)
        );
    }
}

#[test]
fn g0_enthalpy_difference_refuses_unsupported_endpoints_and_overflow() {
    let source = mixture(0.0);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 199.0, 6001.0] {
        for (a, b) in [(bad, 500.0), (500.0, bad), (bad, bad)] {
            assert!(
                source
                    .specific_enthalpy_change(Temperature::new(a), Temperature::new(b))
                    .is_err()
            );
        }
    }
    let mut coefficients = [0.0; 9];
    coefficients[6] = f64::MAX;
    assert!(matches!(
        model(coefficients).molar_enthalpy_change(Temperature::new(500.0), Temperature::new(900.0)),
        Err(ThermochemErrorV1::NonFiniteEvaluation {
            property: Nasa9PropertyV1::Enthalpy,
            ..
        })
    ));
}
