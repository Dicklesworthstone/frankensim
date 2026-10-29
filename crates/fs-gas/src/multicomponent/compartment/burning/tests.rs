use super::*;
use crate::multicomponent::CaloricDomain;
use crate::multicomponent::compartment::ZoneSolve;
use crate::multicomponent::tests::source;
use fs_qty::{Area, Length, MassFlowRate, Pressure};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

const NAMES: [&str; 5] = ["CH4", "CO2", "H2O", "N2", "O2"];
const MASSES: [f64; 5] = [16.0, 44.0, 18.0, 28.0, 32.0];
const NU: [i128; 5] = [-1, 1, 2, 0, -2];
const FEED: [f64; 5] = [0.125, 0.0, 0.0, 0.625, 0.25];
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn reaction_id() -> ReactionId {
    ReactionId::new("manufactured-oxidation").unwrap()
}
fn fixture_mechanism(heat: f64, cp_shift: f64, catalyst: bool) -> IdealGasKineticsV1 {
    // Integer kg/mol masses and manufactured NASA polynomials: this verifies
    // chemistry algebra, not empirical methane thermochemistry or fire rates.
    let ids = NAMES.map(id).to_vec();
    let elements = ElementalMatrix::new(
        ["C", "H", "N", "O"]
            .map(|s| ElementId::new(s).unwrap())
            .to_vec(),
        ids.clone(),
        vec![
            vec![1, 1, 0, 0, 0],
            vec![4, 0, 2, 0, 0],
            vec![0, 0, 0, 2, 0],
            vec![0, 2, 1, 0, 2],
        ],
    )
    .unwrap();
    let formations = [0.0, -0.5 * heat, -0.25 * heat, 0.0, 0.0];
    let models = (0..5)
        .map(|i| {
            source(
                NAMES[i],
                MASSES[i],
                3.5 + if i == 1 { cp_shift } else { 0.0 },
                0.02,
                formations[i],
            )
        })
        .collect();
    let source = IdealGasEquilibriumV1::new(
        elements,
        ChargeVector::new(ids, vec![0; 5]).unwrap(),
        models,
    )
    .unwrap();
    let mut reactants = vec![(id("CH4"), 1), (id("O2"), 2)];
    let mut products = vec![(id("CO2"), 1), (id("H2O"), 2)];
    if catalyst {
        reactants.push((id("N2"), 1));
        products.push((id("N2"), 1));
    }
    let law = ArrheniusRateV1::new(
        Temperature::new(2.0),
        ReactionProgressRateQuantityV1::new(7.0),
        0.0,
        MolarEnergyQuantityV1::new(0.0),
    )
    .unwrap();
    let reaction = ElementaryReactionV1::new(reaction_id(), reactants, products, law).unwrap();
    IdealGasKineticsV1::new(
        source,
        "manufactured-global-closure",
        1,
        vec![reaction],
        0.0,
    )
    .unwrap()
}
fn room(m: &IdealGasKineticsV1, y: [f64; 5], whole_upper: bool) -> TwoZone {
    let gas = MixtureGas::new(
        m.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    TwoZone::new(
        gas,
        Area::new(0.5),
        Length::new(2.0),
        Length::new(if whole_upper { 0.0 } else { 1.0 }),
        Pressure::new(10.0 * R),
        [Temperature::new(2.0); 2],
        [y.to_vec(), y.to_vec()],
        || false,
    )
    .unwrap()
}
fn adapter(room: &TwoZone, m: &IdealGasKineticsV1) -> InventoryLimitedReaction {
    InventoryLimitedReaction::new(room.gas(), m.clone(), &reaction_id(), || false).unwrap()
}
fn evaluate(
    owner: &InventoryLimitedReaction,
    room: &TwoZone,
    extent: f64,
    dt: f64,
) -> BurningReport {
    owner
        .evaluate(
            room,
            Layer::Upper,
            Amount::new(extent),
            Time::new(dt),
            None,
            Temperature::new(2.0),
            || false,
        )
        .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:.17e}, expected={b:.17e}"
    );
}

#[test]
fn g0_oxygen_limits_one_conserved_extent_and_leaves_excess_fuel() {
    let mechanism = fixture_mechanism(40.0, 0.0, false);
    let mut zone = room(&mechanism, FEED, false);
    let owner = adapter(&zone, &mechanism);
    let original = zone.clone();
    let report = evaluate(&owner, &zone, 10.0, 0.25);
    assert_eq!(owner.mechanism(), &mechanism);
    assert_eq!(report.limiting_species, Some(4));
    assert_eq!(owner.species()[4], id("O2"));
    close(
        report.actual_extent.value(),
        zone.layer(Layer::Upper).species_masses()[4] / 64.0,
        1e-15,
    );
    assert!(report.actual_extent.value() > 0.0 && report.actual_extent.value() < 10.0);
    let elements = mechanism.thermodynamics().elemental_matrix();
    for row in 0..elements.row_count() {
        let exact: i128 = (0..5).map(|i| elements.get(row, i).unwrap() * NU[i]).sum();
        assert_eq!(exact, 0);
        close(
            sum((0..5).map(|i| {
                report.chemical_rates.species_kg_s[i] / MASSES[i]
                    * elements.get(row, i).unwrap() as f64
            })),
            0.0,
            2e-15,
        );
    }
    close(report.mass_balance_rate_kg_s, 0.0, 2e-14);
    zone.advance(
        Time::new(0.25),
        [LayerRates::default(), report.chemical_rates],
        ZoneSolve::default(),
        || false,
    )
    .unwrap();
    let before = original.layer(Layer::Upper).species_masses();
    let after = zone.layer(Layer::Upper).species_masses();
    assert!(after.iter().all(|m| *m >= 0.0));
    assert!(after[0] > 0.0 && after[4] <= 1e-13);
    assert_eq!(after[3].to_bits(), before[3].to_bits());
    for i in 0..5 {
        close(
            after[i] - before[i],
            MASSES[i] * NU[i] as f64 * report.actual_extent.value(),
            2e-14,
        );
    }
    close(
        zone.physical_internal_energy().value(),
        original.physical_internal_energy().value(),
        1e-10,
    );
}

#[test]
fn g1_original_formation_energy_heats_variable_cp_gas_exactly_once() {
    let mechanism = fixture_mechanism(40.0, 0.0, false);
    let mut zone = room(&mechanism, FEED, true);
    let owner = adapter(&zone, &mechanism);
    let before = zone.physical_internal_energy().value();
    let amount = sum(zone
        .layer(Layer::Upper)
        .species_masses()
        .iter()
        .zip(MASSES)
        .map(|(m, mass)| m / mass));
    let report = evaluate(&owner, &zone, 0.1, 0.3);
    zone.advance(
        Time::new(0.3),
        [LayerRates::default(), report.chemical_rates],
        ZoneSolve::default(),
        || false,
    )
    .unwrap();
    // Net molecular count is zero. Original Cv/R=2.5+0.02T, and the
    // reaction's original formation-energy change is exactly -40 R per mol.
    let c = 40.0 * report.actual_extent.value() / amount;
    let expected = 2.0 + 2.0 * c / (2.54 + (2.54_f64.powi(2) + 0.04 * c).sqrt());
    close(
        zone.layer(Layer::Upper).temperature().value(),
        expected,
        2e-9,
    );
    close(zone.physical_internal_energy().value(), before, 2e-10);
    close(report.reaction_enthalpy.value(), -40.0 * R, 2e-13);
    close(
        report.released_enthalpy.value(),
        40.0 * R * report.actual_extent.value(),
        1e-13,
    );
    assert_eq!(zone.layer(Layer::Lower).mass_kg(), 0.0);
    assert_eq!(zone.interface_height().value(), 0.0);
}

#[test]
fn g0_injected_fuel_is_available_immediately_without_repeating_inflow() {
    let mechanism = fixture_mechanism(40.0, 0.0, false);
    let mut zone = room(&mechanism, [0.0, 0.0, 0.0, 0.75, 0.25], false);
    let owner = adapter(&zone, &mechanism);
    let original = zone.clone();
    assert_eq!(evaluate(&owner, &zone, 1.0, 0.3).actual_extent.value(), 0.0);
    let incoming = zone
        .inflow(
            MassFlowRate::new(1.6),
            &[1.0, 0.0, 0.0, 0.0, 0.0],
            Temperature::new(2.0),
        )
        .unwrap();
    let report = owner
        .evaluate(
            &zone,
            Layer::Upper,
            Amount::new(1.0),
            Time::new(0.3),
            Some(&incoming),
            Temperature::new(2.0),
            || false,
        )
        .unwrap();
    assert_eq!(report.limiting_species, Some(0));
    close(report.actual_extent.value(), 0.03, 1e-16);
    assert!(report.chemical_rates.species_kg_s[0] < 0.0);
    let mut combined = incoming;
    combined.add(&report.chemical_rates).unwrap();
    zone.advance(
        Time::new(0.3),
        [LayerRates::default(), combined],
        ZoneSolve::default(),
        || false,
    )
    .unwrap();
    let total_mass =
        |z: &TwoZone| z.layer(Layer::Lower).mass_kg() + z.layer(Layer::Upper).mass_kg();
    close(total_mass(&zone) - total_mass(&original), 0.48, 2e-14);
    let incoming_physical_power = incoming.referenced_energy_w
        + sum((0..5).map(|i| incoming.species_kg_s[i] * zone.gas().energy_references()[i]));
    close(
        zone.physical_internal_energy().value() - original.physical_internal_energy().value(),
        0.3 * incoming_physical_power,
        2e-10,
    );
    assert!(zone.layer(Layer::Upper).species_masses()[0] <= 1e-15);
    assert!(zone.layer(Layer::Upper).species_masses()[1] > 0.0);
}

#[test]
fn g0_common_extent_rounds_inward_instead_of_clipping_oxygen() {
    let mechanism = fixture_mechanism(40.0, 0.0, false);
    let mut zone = room(&mechanism, FEED, false);
    let owner = adapter(&zone, &mechanism);
    let dt = 0.06405707420669589;
    let oxygen = zone.layer(Layer::Upper).species_masses()[4];
    let naive_rate = (oxygen / 64.0) / dt;
    assert!(oxygen + dt * (-64.0 * naive_rate) < 0.0);
    let report = evaluate(&owner, &zone, 10.0, dt);
    assert!(report.actual_extent.value() <= oxygen / 64.0);
    assert!(report.actual_extent.value() >= (oxygen / 64.0) * (1.0 - 4.0 * f64::EPSILON));
    let rate = report.chemical_rates.species_kg_s[0] / -16.0;
    for i in 0..5 {
        assert_eq!(
            report.chemical_rates.species_kg_s[i],
            MASSES[i] * NU[i] as f64 * rate
        );
    }
    zone.advance(
        Time::new(dt),
        [LayerRates::default(), report.chemical_rates],
        ZoneSolve::default(),
        || false,
    )
    .unwrap();
    assert!(zone.layer(Layer::Upper).species_masses()[4] >= 0.0);
}

#[test]
fn g0_power_uses_original_temperature_dependent_reaction_enthalpy() {
    let mechanism = fixture_mechanism(40.0, 1.5, false);
    let zone = room(&mechanism, FEED, false);
    let owner = adapter(&zone, &mechanism);
    for t in [2.0, 8.0] {
        let extent = owner
            .extent_for_power(
                Power::new(20.0),
                Time::new(0.3),
                Temperature::new(t),
                || false,
            )
            .unwrap();
        close(extent.value(), 6.0 / ((40.0 - 1.5 * t) * R), 1e-16);
        let report = owner
            .evaluate(
                &zone,
                Layer::Upper,
                extent,
                Time::new(0.3),
                None,
                Temperature::new(t),
                || false,
            )
            .unwrap();
        close(report.released_power.value(), 20.0, 1e-13);
        close(report.released_enthalpy.value(), 6.0, 1e-13);
    }
    assert!(
        owner
            .extent_for_power(
                Power::new(20.0),
                Time::new(0.3),
                Temperature::new(40.0),
                || false
            )
            .is_err()
    );
}

#[test]
fn g4_cancellation_zero_oxygen_and_refusals_preserve_owners() {
    let mechanism = fixture_mechanism(40.0, 0.0, false);
    let zone = room(&mechanism, FEED, false);
    let owner = adapter(&zone, &mechanism);
    let original = (zone.clone(), owner.clone());
    let mut calls = 0;
    let baseline = owner
        .evaluate(
            &zone,
            Layer::Upper,
            Amount::new(0.1),
            Time::new(0.3),
            None,
            Temperature::new(2.0),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut seen = 0;
        assert_eq!(
            owner.evaluate(
                &zone,
                Layer::Upper,
                Amount::new(0.1),
                Time::new(0.3),
                None,
                Temperature::new(2.0),
                || {
                    seen += 1;
                    seen == stop
                }
            ),
            Err(BurningError::Cancelled)
        );
    }
    let other = room(&fixture_mechanism(20.0, 0.0, false), FEED, false);
    assert_eq!(
        owner.evaluate(
            &other,
            Layer::Upper,
            Amount::new(0.1),
            Time::new(0.3),
            None,
            Temperature::new(2.0),
            || false
        ),
        Err(BurningError::SourceMismatch)
    );
    assert_eq!(
        InventoryLimitedReaction::new(other.gas(), mechanism.clone(), &reaction_id(), || false)
            .err(),
        Some(BurningError::SourceMismatch)
    );
    let overlap = fixture_mechanism(40.0, 0.0, true);
    let overlap_room = room(&overlap, FEED, false);
    assert!(matches!(
        InventoryLimitedReaction::new(overlap_room.gas(), overlap, &reaction_id(), || false),
        Err(BurningError::UnsupportedReaction(_))
    ));
    let no_oxygen = room(&mechanism, [0.125, 0.0, 0.0, 0.875, 0.0], false);
    let empty = evaluate(&owner, &no_oxygen, 1.0, 0.3);
    assert_eq!(empty.actual_extent.value(), 0.0);
    assert_eq!(empty.chemical_rates, LayerRates::default());
    let zero = evaluate(&owner, &zone, 0.0, 0.3);
    assert_eq!(zero.actual_extent.value(), 0.0);
    assert_eq!(zero.limiting_species, None);
    let mut invalid_incoming = LayerRates::default();
    invalid_incoming.species_kg_s[4] = -1.0;
    assert!(
        owner
            .evaluate(
                &zone,
                Layer::Upper,
                Amount::new(0.1),
                Time::new(0.3),
                Some(&invalid_incoming),
                Temperature::new(2.0),
                || false
            )
            .is_err()
    );
    for (extent, dt, temperature) in [
        (-1.0, 0.3, 2.0),
        (1.0, 0.0, 2.0),
        (1.0, 0.3, 100.0),
        (f64::NAN, 0.3, 2.0),
    ] {
        assert!(
            owner
                .evaluate(
                    &zone,
                    Layer::Upper,
                    Amount::new(extent),
                    Time::new(dt),
                    None,
                    Temperature::new(temperature),
                    || false
                )
                .is_err()
        );
    }
    assert_eq!((zone.clone(), owner.clone()), original);
    assert_eq!(evaluate(&owner, &zone, 0.1, 0.3), baseline);
}
