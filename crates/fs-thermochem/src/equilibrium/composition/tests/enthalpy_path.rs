//! G1/G3 conserved enthalpy and source-path regressions.
use super::*;

fn split_source(gap: f64) -> Nasa9StandardStateModelV1 {
    split_source_with_shift(gap, 10000.0)
}

fn split_source_with_shift(gap: f64, shift_j_mol: f64) -> Nasa9StandardStateModelV1 {
    let original = model("A", 0.01, 4.0, 0.0, 0.0);
    let mut lower = original.regions()[0].card().clone();
    lower.validity = ValidityDomain::unconstrained().with("T", 200.0, 1000.0);
    let mut upper = lower.clone();
    upper.validity = ValidityDomain::unconstrained().with("T", 1000.0 + gap, 6000.0);
    upper.parameters.get_mut("a7").unwrap().value += shift_j_mol / R;
    Nasa9StandardStateModelV1::new(
        id("A"),
        MolarMass::new(0.01),
        original.convention().clone(),
        vec![
            Nasa9RegionV1::from_card(lower).unwrap(),
            Nasa9RegionV1::from_card(upper).unwrap(),
        ],
    )
    .unwrap()
}

#[test]
fn g3_hp_refuses_an_unsupported_inlet_path_outside_the_product_bracket() {
    for gap in [0.0, 100.0] {
        let chemistry = owner(
            vec![split_source_with_shift(gap, -4000.0 * R)],
            &["H"],
            vec![vec![1]],
        );
        for (inlet, minimum, maximum, false_root) in [
            (1500.0, 200.0, 800.0, 500.0),
            (500.0, 1200.0, 1800.0, 1500.0),
        ] {
            let controls = EquilibriumEnthalpySolveV1 {
                minimum_temperature_k: minimum,
                maximum_temperature_k: maximum,
                ..hp_controls()
            };
            let feed = [(id("A"), Amount::new(1.0))];
            let target = chemistry
                .feed_enthalpy(Temperature::new(inlet), &feed, || false)
                .unwrap();
            // The bracket alone has a perfectly converged but physically false
            // root: the unsupported inlet-to-product chart jump supplies heat.
            let explicit = chemistry
                .equilibrate_hp_enthalpy(Pressure::new(100000.0), &feed, target, controls, || false)
                .unwrap();
            close(explicit.temperature().value(), false_root, 1e-12);
            assert!(matches!(
                chemistry.equilibrate_hp(
                    Pressure::new(100000.0),
                    &feed,
                    Temperature::new(inlet),
                    controls,
                    || false,
                ),
                Err(CompositionEquilibriumErrorV1::Invalid(_))
            ));
        }
    }
}

#[test]
fn g1_hp_inlet_path_within_one_original_chart_needs_no_reconciliation() {
    let chemistry = owner(vec![split_source(100.0)], &["H"], vec![vec![1]]);
    let controls = EquilibriumEnthalpySolveV1 {
        maximum_temperature_k: 800.0,
        ..hp_controls()
    };
    let result = chemistry
        .equilibrate_hp(
            Pressure::new(100000.0),
            &[(id("A"), Amount::new(1.0))],
            Temperature::new(500.0),
            controls,
            || false,
        )
        .unwrap();
    close(result.temperature().value(), 500.0, 1e-12);
}

#[test]
fn g1_hp_reconciles_inlet_path_without_demanding_absent_product_inlet_data() {
    let reactant = split_source(0.0)
        .with_caloric_continuity(
            Temperature::new(500.0),
            MolarEnergyQuantityV1::new(10001.0),
            crate::MolarThermalQuantityV1::new(1.0),
        )
        .unwrap();
    let original = model("B", 0.01, 4.0, -20000.0, -20000.0 / (R * REFERENCE_T));
    let mut card = original.regions()[0].card().clone();
    card.validity = ValidityDomain::unconstrained().with("T", 1000.0, 6000.0);
    let product = Nasa9StandardStateModelV1::new(
        id("B"),
        original.molar_mass(),
        original.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let chemistry = owner(vec![reactant, product], &["H"], vec![vec![1, 1]]);
    let controls = EquilibriumEnthalpySolveV1 {
        minimum_temperature_k: 1050.0,
        maximum_temperature_k: 1500.0,
        ..hp_controls()
    };
    // Explicit zero feed is admitted without evaluating B at the 600 K inlet,
    // where its original source has no support. Both products support the root.
    let result = chemistry
        .equilibrate_hp(
            Pressure::new(100000.0),
            &[(id("A"), Amount::new(1.0)), (id("B"), Amount::new(0.0))],
            Temperature::new(600.0),
            controls,
            || false,
        )
        .unwrap();
    let (mut low, mut high) = (1050.0, 1500.0);
    for _ in 0..80 {
        let temperature = 0.5 * (low + high);
        let product_fraction = 1.0 / (1.0 + (-20000.0 / (R * temperature)).exp());
        let residual = 4.0 * R * (temperature - 600.0) - 20000.0 * product_fraction;
        if residual < 0.0 {
            low = temperature;
        } else {
            high = temperature;
        }
    }
    close(result.temperature().value(), 0.5 * (low + high), 2e-9);
    assert!(result.enthalpy_residual().unwrap().value().abs() <= controls.enthalpy_tolerance_j);
}

#[test]
fn g3_hp_refuses_disconnected_source_charts_before_endpoint_acceptance() {
    for gap in [0.0, 100.0] {
        let chemistry = owner(vec![split_source(gap)], &["H"], vec![vec![1]]);
        let feed = [(id("A"), Amount::new(1.0))];
        // A target exactly on the lower bracket endpoint used to bypass the
        // unsupported/discontinuous rest of the caller's requested path.
        for inlet in [200.0, 1500.0] {
            assert!(matches!(
                chemistry.equilibrate_hp(
                    Pressure::new(100000.0),
                    &feed,
                    Temperature::new(inlet),
                    hp_controls(),
                    || false,
                ),
                Err(CompositionEquilibriumErrorV1::Invalid(_))
            ));
        }
    }
}

#[test]
fn g1_hp_crosses_original_charts_only_with_explicit_caloric_continuity() {
    let corrected = split_source(0.0)
        .with_caloric_continuity(
            Temperature::new(500.0),
            MolarEnergyQuantityV1::new(10001.0),
            crate::MolarThermalQuantityV1::new(1.0),
        )
        .unwrap();
    let chemistry = owner(vec![corrected], &["H"], vec![vec![1]]);
    let result = chemistry
        .equilibrate_hp_enthalpy(
            Pressure::new(100000.0),
            &[(id("A"), Amount::new(1.0))],
            Energy::new(4.0 * R * 1500.0),
            hp_controls(),
            || false,
        )
        .unwrap();
    close(result.temperature().value(), 1500.0, 1e-12);
    assert!(
        result.enthalpy_residual().unwrap().value().abs() <= hp_controls().enthalpy_tolerance_j
    );
}
