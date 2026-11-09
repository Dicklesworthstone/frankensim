//! G1/G3 conserved enthalpy and source-path regressions.
use super::*;

fn split_source(gap: f64) -> Nasa9StandardStateModelV1 {
    let original = model("A", 0.01, 4.0, 0.0, 0.0);
    let mut lower = original.regions()[0].card().clone();
    lower.validity = ValidityDomain::unconstrained().with("T", 200.0, 1000.0);
    let mut upper = lower.clone();
    upper.validity = ValidityDomain::unconstrained().with("T", 1000.0 + gap, 6000.0);
    upper.parameters.get_mut("a7").unwrap().value += 10000.0 / R;
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
