use super::*;
use crate::mixture::FrozenIdealGasMixtureModelV1;

#[test]
fn g0_entropy_change_retains_tiny_changes_despite_large_reference_entropy() {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[8] = 1e24;
    let source = model(coefficients);
    let (a, b) = (500.0, 500.0 + 1e-6);
    let old = source
        .evaluate(Temperature::new(a))
        .unwrap()
        .properties()
        .s()
        .value();
    let new = source
        .evaluate(Temperature::new(b))
        .unwrap()
        .properties()
        .s()
        .value();
    assert_eq!(old, new, "absolute entropy subtraction loses this change");
    let change = source
        .molar_entropy_change(Temperature::new(a), Temperature::new(b))
        .unwrap()
        .value();
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 3.5 * ((b - a) / a).ln_1p();
    assert!((change / expected - 1.0).abs() < 5e-14);
    let reversed = source
        .molar_entropy_change(Temperature::new(b), Temperature::new(a))
        .unwrap()
        .value();
    assert_eq!(reversed, -change);
    assert_eq!(
        source
            .molar_entropy_change(Temperature::new(a), Temperature::new(a))
            .unwrap()
            .value(),
        0.0
    );
}

#[test]
fn g0_entropy_change_matches_each_nasa9_channel() {
    let (a, b) = (400.0_f64, 900.0_f64);
    let expected = [
        0.5 * (1.0 / a.powi(2) - 1.0 / b.powi(2)),
        1.0 / a - 1.0 / b,
        (b / a).ln(),
        b - a,
        (b.powi(2) - a.powi(2)) / 2.0,
        (b.powi(3) - a.powi(3)) / 3.0,
        (b.powi(4) - a.powi(4)) / 4.0,
        0.0,
        0.0,
    ];
    for channel in 0..9 {
        let mut coefficients = [0.0; 9];
        coefficients[channel] = 1.0;
        let source = model(coefficients);
        let actual = source
            .molar_entropy_change(Temperature::new(a), Temperature::new(b))
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
fn g0_entropy_change_preserves_cross_region_reference_offsets() {
    let mut low = [0.0; 9];
    low[2] = 3.5;
    low[8] = 1e12;
    let mut high = low;
    high[8] += 0.25;
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
    // The source explicitly has a 0.25 R jump. Never smooth it away.
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (3.5 * (600.0_f64 / 500.0).ln() + 0.25);
    let actual = source
        .molar_entropy_change(Temperature::new(500.0), Temperature::new(600.0))
        .unwrap()
        .value();
    assert!((actual / expected - 1.0).abs() < 5e-14);
}

fn mixture_with_reference(reference: f64) -> FrozenIdealGasMixtureModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[8] = reference;
    let source = model(coefficients);
    FrozenIdealGasMixtureModelV1::new(
        vec![source],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap()
}

#[test]
fn g1_frozen_mixture_entropy_is_reference_invariant_and_contains_pressure_work() {
    let a = mixture_with_reference(0.0);
    let b = mixture_with_reference(1e24);
    let (t0, t1) = (500.0, 500.0 + 1e-6);
    let (p0, p1) = (1e5, 1e5 * (t1 / t0)); // same density
    let delta = |model: &FrozenIdealGasMixtureModelV1| {
        model
            .specific_entropy_change(
                Temperature::new(t0),
                Pressure::new(p0),
                Temperature::new(t1),
                Pressure::new(p1),
            )
            .unwrap()
            .value()
    };
    assert_eq!(delta(&a), delta(&b));
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028_013_4
        * (3.5 * ((t1 - t0) / t0).ln_1p() - ((p1 - p0) / p0).ln_1p());
    assert!((delta(&a) / expected - 1.0).abs() < 5e-14);
    let isothermal = a
        .molar_entropy_change(
            Temperature::new(500.0),
            Pressure::new(1e5),
            Temperature::new(500.0),
            Pressure::new(2e5),
        )
        .unwrap()
        .value();
    assert_close(
        isothermal,
        -UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 2.0_f64.ln(),
    );
}

#[test]
fn g0_entropy_difference_refuses_unsupported_endpoints_and_handles_extreme_pressure_ratio() {
    let source = mixture_with_reference(0.0);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 199.0, 6001.0] {
        assert!(
            source
                .specific_entropy_change(
                    Temperature::new(bad),
                    Pressure::new(1e5),
                    Temperature::new(bad),
                    Pressure::new(1e5)
                )
                .is_err()
        );
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            source
                .specific_entropy_change(
                    Temperature::new(500.0),
                    Pressure::new(bad),
                    Temperature::new(500.0),
                    Pressure::new(1e5)
                )
                .is_err()
        );
    }
    let value = source
        .molar_entropy_change(
            Temperature::new(500.0),
            Pressure::new(f64::MIN_POSITIVE),
            Temperature::new(500.0),
            Pressure::new(f64::MAX),
        )
        .unwrap()
        .value();
    let expected = -UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (f64::MAX.ln() - f64::MIN_POSITIVE.ln());
    assert!((value / expected - 1.0).abs() < 5e-14);
}

#[test]
fn g0_isochoric_entropy_retains_small_cv_and_cancels_reference_offsets() {
    for cp in [1.0_f64.next_up(), 2.5, 3.5] {
        for reference in [-1e24, 0.0, 1e24] {
            let mut coefficients = [0.0; 9];
            coefficients[2] = cp;
            coefficients[8] = reference;
            let source = model(coefficients);
            for (a, b) in [(500.0_f64, 500.0_f64.next_up()), (500.0, 1800.0)] {
                let actual = source
                    .molar_isochoric_entropy_change(Temperature::new(a), Temperature::new(b))
                    .unwrap()
                    .value();
                let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                    * ((cp - 1.0) * ((b - a) / a).ln_1p());
                assert!(actual > 0.0);
                assert!((actual / expected - 1.0).abs() < 5e-14);
                assert_eq!(
                    source
                        .molar_isochoric_entropy_change(Temperature::new(b), Temperature::new(a))
                        .unwrap()
                        .value(),
                    -actual
                );
                assert_eq!(
                    source
                        .molar_isochoric_entropy_change(Temperature::new(a), Temperature::new(a))
                        .unwrap()
                        .value(),
                    0.0
                );
            }
        }
    }
}

#[test]
fn g1_isochoric_mixture_entropy_integrates_variable_cv_on_its_mass_basis() {
    let mut nitrogen = [0.0; 9];
    nitrogen[2] = 3.5;
    nitrogen[3] = 0.001;
    nitrogen[8] = 1e24;
    let mut helium = [0.0; 9];
    helium[2] = 2.5;
    helium[8] = -1e24;
    let source = FrozenIdealGasMixtureModelV1::new(
        vec![
            model(nitrogen),
            Nasa9StandardStateModelV1::new(
                SpeciesId::new("He").unwrap(),
                MolarMass::new(0.004),
                convention(100_000.0),
                vec![Nasa9RegionV1::from_card(card(200.0, 6000.0, 100_000.0, helium)).unwrap()],
            )
            .unwrap(),
        ],
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    )
    .unwrap();
    for (a, b) in [(500.0_f64, 500.0_f64.next_up()), (500.0, 1800.0)] {
        let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            * ((0.75 * 2.5 + 0.25 * 1.5) * ((b - a) / a).ln_1p()
                + 0.75 * 0.001 * (b - a));
        let mass = 0.75 * 0.028_013_4 + 0.25 * 0.004;
        let molar = source
            .molar_isochoric_entropy_change(Temperature::new(a), Temperature::new(b))
            .unwrap()
            .value();
        let specific = source
            .specific_isochoric_entropy_change(Temperature::new(a), Temperature::new(b))
            .unwrap()
            .value();
        assert!((molar / expected - 1.0).abs() < 5e-14);
        assert!((specific / (expected / mass) - 1.0).abs() < 5e-14);
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 199.0, 6001.0] {
        assert!(source
            .specific_isochoric_entropy_change(Temperature::new(bad), Temperature::new(500.0))
            .is_err());
    }
}

#[test]
fn g3_isochoric_entropy_preserves_raw_jumps_and_continuous_source_seams() {
    let mut low = [0.0; 9];
    low[2] = 3.5;
    low[8] = 1e12;
    let mut high = low;
    high[2] = 4.0;
    high[3] = 0.001;
    high[8] += 0.25;
    let raw = Nasa9StandardStateModelV1::new(
        SpeciesId::new("N2").unwrap(),
        MolarMass::new(0.028),
        convention(100_000.0),
        vec![
            Nasa9RegionV1::from_card(card(200.0, 600.0, 100_000.0, low)).unwrap(),
            Nasa9RegionV1::from_card(card(600.0, 1500.0, 100_000.0, high)).unwrap(),
        ],
    )
    .unwrap();
    let actual = raw
        .molar_isochoric_entropy_change(Temperature::new(500.0), Temperature::new(600.0))
        .unwrap()
        .value();
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (3.0 * 600.0_f64.ln() + 0.001 * 600.0 - 2.5 * 500.0_f64.ln() + 0.25);
    assert!((actual / expected - 1.0).abs() < 5e-14);
    let continuous = raw
        .with_caloric_continuity(
            Temperature::new(500.0),
            MolarEnergyQuantityV1::new(1e5),
            MolarThermalQuantityV1::new(1e3),
        )
        .unwrap();
    let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (2.5 * (600.0_f64 / 500.0).ln()
            + 3.0 * (900.0_f64 / 600.0).ln()
            + 0.001 * (900.0 - 600.0));
    let actual = continuous
        .molar_isochoric_entropy_change(Temperature::new(500.0), Temperature::new(900.0))
        .unwrap()
        .value();
    assert!((actual / expected - 1.0).abs() < 5e-14);
}
