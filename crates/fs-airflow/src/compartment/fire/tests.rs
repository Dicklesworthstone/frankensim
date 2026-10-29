use super::*;
use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::compartment::ZoneSolve;
use fs_gas::multicomponent::{CaloricDomain, MixtureGas};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Area, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

fn model(name: &str, mass: f64, cp: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    let values = [0.0, 0.0, cp, slope, 0.0, 0.0, 0.0, 1000.0, 0.0];
    let mut parameters: std::collections::BTreeMap<_, _> = required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(values)
        .map(|((key, dims), value)| (key.into(), LawParameter { value, dims }))
        .collect();
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        sources: vec![],
        provenance: Provenance {
            source: "analytic prescribed-fire caloric fixture".into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100_000.0),
            ElementalReferenceIdV1::new("plume-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn room(interface: f64, lower_temperature: f64, lower_fractions: [f64; 2]) -> TwoZone {
    let models = vec![model("A", 0.028, 3.5, 0.002), model("B", 0.032, 4.0, 0.0)];
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 0]).unwrap(),
        models,
    )
    .unwrap();
    let gas = MixtureGas::new(
        bank,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(2000.0)),
        || false,
    )
    .unwrap();
    TwoZone::new(
        gas,
        Area::new(10.0),
        Length::new(3.0),
        Length::new(interface),
        Pressure::new(100_000.0),
        [
            Temperature::new(lower_temperature),
            Temperature::new(lower_temperature + 50.0),
        ],
        [lower_fractions.to_vec(), vec![0.0, 1.0]],
        || false,
    )
    .unwrap()
}
fn fire() -> PrescribedFire {
    PrescribedFire {
        total_power: Power::new(60_000.0),
        convective_power: Power::new(48_000.0),
        diameter: Length::new(0.3),
        base_height: Length::new(0.1),
    }
}
fn evaluate(fire: PrescribedFire, room: &TwoZone) -> PrescribedFireReport {
    fire.evaluate(room, Acceleration::new(9.81), || false)
        .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_prescribed_heat_and_identical_plume_transfer_conserve_extensives() {
    let mut room = room(1.2, 300.0, [0.75, 0.25]);
    let source = evaluate(fire(), &room);
    assert!(source.plume_mass_rate.value() > 0.0);
    assert_eq!(source.radiant_power.value(), 12_000.0);
    for i in 0..room.gas().species_count() {
        assert_eq!(
            source.rates[0].species_kg_s[i],
            -source.rates[1].species_kg_s[i]
        );
    }
    close(
        source.rates[0].referenced_energy_w + source.rates[1].referenced_energy_w,
        48_000.0,
    );
    let original_energy = room.physical_internal_energy().value();
    let original_mass = room.layer(Layer::Lower).mass_kg() + room.layer(Layer::Upper).mass_kg();
    let result = room
        .advance(Time::new(0.01), source.rates, ZoneSolve::default(), || {
            false
        })
        .unwrap();
    close(
        room.layer(Layer::Lower).mass_kg() + room.layer(Layer::Upper).mass_kg(),
        original_mass,
    );
    assert!((room.physical_internal_energy().value() - original_energy - 480.0).abs() < 1e-6);
    assert!(result.physical_energy_residual_j.abs() < 1e-6);
    assert_eq!(
        source.plume_estimate.unwrap().card().correlation,
        fs_convection::fire::FireCorrelationIdV1::HeskestadEntrainment
    );
}

#[test]
fn g1_plume_ambient_uses_actual_original_variable_cp_and_density() {
    let mut mass_rates = [0.0; 2];
    for (index, t) in [300.0, 600.0].into_iter().enumerate() {
        let zone = room(1.2, t, [0.75, 0.25]);
        let source = evaluate(fire(), &zone);
        let ambient = source.plume_ambient.unwrap();
        // Independent analytic properties of the authored original NASA bank.
        let cp = R * (0.75 * (3.5 + 0.002 * t) / 0.028 + 0.25 * 4.0 / 0.032);
        let rho = 100_000.0 / (R * t * (0.75 / 0.028 + 0.25 / 0.032));
        close(ambient.temperature().value(), t);
        close(ambient.density().value(), rho);
        close(ambient.specific_heat().value(), cp);
        mass_rates[index] = source.plume_mass_rate.value();
    }
    assert!((mass_rates[0] - mass_rates[1]).abs() > 0.01 * mass_rates[0]);
}

#[test]
fn g1_absent_upper_is_created_only_by_real_entrained_species() {
    let mut zone = room(3.0, 300.0, [1.0, 0.0]);
    assert_eq!(zone.layer(Layer::Upper).mass_kg(), 0.0);
    let source = evaluate(fire(), &zone);
    assert!(matches!(
        source.regime,
        PrescribedFireRegime::Entrainment(_)
    ));
    let duration = 0.001;
    zone.advance(
        Time::new(duration),
        source.rates,
        ZoneSolve::default(),
        || false,
    )
    .unwrap();
    close(
        zone.layer(Layer::Upper).mass_kg(),
        duration * source.plume_mass_rate.value(),
    );
    assert_eq!(zone.layer(Layer::Upper).species_masses()[1], 0.0);
    assert!(zone.layer(Layer::Upper).volume().value() > 0.0);
    assert!(zone.interface_height().value() < 3.0);
}

#[test]
fn g4_dormant_upper_source_invalid_regime_and_cancellation() {
    for (interface, base) in [(0.0, 0.1), (1.2, 1.2), (1.2, 1.5)] {
        let zone = room(interface, 300.0, [0.75, 0.25]);
        let mut fire = fire();
        fire.base_height = Length::new(base);
        let report = evaluate(fire, &zone);
        assert_eq!(report.regime, PrescribedFireRegime::InUpperLayer);
        assert_eq!(report.plume_mass_rate.value(), 0.0);
        assert_eq!(report.plume_estimate, None);
        assert_eq!(report.plume_ambient, None);
        assert_eq!(report.rates[0], LayerRates::default());
        assert!(report.rates[1].species_kg_s.iter().all(|m| *m == 0.0));
        assert_eq!(report.rates[1].referenced_energy_w, 48_000.0);
    }
    let zone = room(1.2, 300.0, [0.75, 0.25]);
    let mut source = fire();
    source.total_power = Power::new(0.0);
    source.convective_power = Power::new(0.0);
    let dormant = evaluate(source, &zone);
    assert_eq!(dormant.regime, PrescribedFireRegime::Dormant);
    assert_eq!(dormant.rates, [LayerRates::default(); 2]);
    source.total_power = Power::new(1.0);
    source.convective_power = Power::new(0.8);
    assert!(matches!(
        source.evaluate(&zone, Acceleration::new(9.81), || false),
        Err(PrescribedFireError::Correlation(_))
    ));
    source = fire();
    source.base_height = Length::new(3.0);
    assert!(
        source
            .evaluate(&zone, Acceleration::new(9.81), || false)
            .is_err()
    );
    source = fire();
    source.convective_power = Power::new(60_001.0);
    assert!(
        source
            .evaluate(&zone, Acceleration::new(9.81), || false)
            .is_err()
    );
    for stop in [1, 2, 4] {
        let mut calls = 0;
        assert_eq!(
            fire().evaluate(&zone, Acceleration::new(9.81), || {
                calls += 1;
                calls == stop
            }),
            Err(PrescribedFireError::Cancelled)
        );
    }
}
