use super::*;

#[test]
fn molecular_background_keeps_original_heavy_calorics_outside_electron_temperature_range() {
    let (atomic, source) = initial(0.05, 30000.0, 4000.0, 0.0, 0.0);
    let original = atomic.thermodynamics();
    let translational = species("N2", 0.028 / NA, 1.0, 0.0);
    let mut card = translational.regions()[0].card().clone();
    card.parameters.get_mut("a2").unwrap().value = 3.5;
    card.validity = ValidityDomain::unconstrained().with("T", 1000.0, 8000.0);
    card.provenance.source =
        "manufactured frozen diatomic Cv=5R/2 background; only heavy-T support".into();
    let nitrogen = Nasa9StandardStateModelV1::new(
        id("N2"),
        MolarMass::new(0.028),
        original.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let ids = vec![id("H"), id("H+"), id("N2"), id("e")];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("N").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0, 0], vec![0, 0, 2, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, 0, -1]).unwrap(),
        vec![
            original.models()[0].clone(),
            original.models()[1].clone(),
            nitrogen,
            original.models()[2].clone(),
        ],
    )
    .unwrap();
    let state = FrozenIonizedGasStateV1::new(
        &bank,
        &[
            Amount::new(0.95),
            Amount::new(0.05),
            Amount::new(2.0),
            Amount::new(0.05),
        ],
        Volume::new(1.0),
        Temperature::new(4000.0),
        Temperature::new(30000.0),
        1e-12,
        || false,
    )
    .unwrap();
    let channel = AtomicIonizationV1::new(
        bank,
        ReactionId::new("H-impact-in-N2").unwrap(),
        id("H+"),
        source.forward_rate().clone(),
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap();
    let mut controls = controls();
    controls.maximum_heavy_temperature_k = 8000.0;
    let duration = 2e-8;
    let heavy_power = 1e7;
    let result = advance_atomic_ionization(
        &state,
        &channel,
        Power::new(0.0),
        Power::new(heavy_power),
        Time::new(duration),
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(result.state().amounts()[2], Amount::new(2.0));
    assert!(result.reaction_extent().value() > 0.0);
    assert!(result.state().electron_temperature().value() > 8000.0);
    assert!(
        (result.state().heavy_temperature().value() - 4000.0 - heavy_power * duration / (6.5 * R))
            .abs()
            < 1e-9
    );
    assert!(result.balance().energy_residual_j.abs() < controls.energy_tolerance_j);
}
