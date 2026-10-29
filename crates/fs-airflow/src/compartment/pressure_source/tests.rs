use super::super::*;
use fs_gas::multicomponent::MixtureGas;
use fs_gas::multicomponent::compartment::burning::{BurningReport, InventoryLimitedReaction};
use fs_qty::{Amount, Area};
use fs_thermochem::ReactionId;
use std::cell::Cell;

mod kinetics;

// Reuse the exact original NASA bank exercised by the public reacting example.
#[path = "../../../examples/ventilated_compartment_fire/nasa9.rs"]
mod nasa9;

fn room(gas: MixtureGas, fractions: Vec<f64>, temperature: f64) -> Room {
    Room {
        gas: TwoZone::new(
            gas,
            Area::new(1.0),
            Length::new(1.0),
            Length::new(0.0),
            Pressure::new(101325.0),
            [Temperature::new(temperature); 2],
            [fractions.clone(), fractions],
            || false,
        )
        .unwrap(),
        floor_elevation: Length::new(0.0),
    }
}

fn opening(width: f64) -> WallOpening {
    WallOpening::new(
        Length::new(0.0),
        Length::new(1.0),
        Length::new(width),
        stratified::DischargeCoefficient::sharp_edged().unwrap(),
    )
    .unwrap()
}

fn fixture(
    width: f64,
    fuel_fraction: f64,
    exhaust_rate: f64,
) -> (ZoneNetwork, InventoryLimitedReaction) {
    let (gas, mechanism, air) = nasa9::source().unwrap();
    let reaction = InventoryLimitedReaction::new(
        &gas,
        mechanism,
        &ReactionId::new(nasa9::REACTION).unwrap(),
        || false,
    )
    .unwrap();
    let mut initial: Vec<_> = air.iter().map(|y| (1.0 - fuel_fraction) * y).collect();
    initial[0] = fuel_fraction;
    let network = ZoneNetwork::new(
        vec![room(gas, initial, 300.0)],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: air,
        },
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: opening(width),
        }],
        vec![MechanicalVent {
            from: Endpoint::Room(0),
            from_layer: Layer::Upper,
            to: Endpoint::Ambient,
            to_layer: Layer::Upper,
            mass_rate: MassFlowRate::new(exhaust_rate),
            deposited_power: Power::new(0.0),
        }],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    (network, reaction)
}

fn burn(
    reaction: &InventoryLimitedReaction,
    rooms: &mut [Room],
    dt: Time,
    additional: &mut LayerRates,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<BurningReport, NetworkError> {
    let result = reaction
        .apply_split(
            &mut rooms[0].gas,
            Layer::Upper,
            Amount::new(1000.0),
            dt,
            Temperature::new(300.0),
            ZoneSolve::default(),
            cancelled,
        )?
        .burning;
    additional.add(&result.chemical_rates)?;
    Ok(result)
}

fn feed(rooms: &[Room], rates: &mut [[LayerRates; 2]]) -> Result<(), NetworkError> {
    rates[0][1].add(&rooms[0].gas.inflow(
        MassFlowRate::new(0.002),
        &[1.0, 0.0, 0.0, 0.0, 0.0],
        Temperature::new(300.0),
    )?)?;
    Ok(())
}

fn balance(report: &NetworkStepReport) {
    assert!(report.maximum_pressure_residual_pa <= 1e-5, "{report:?}");
    assert_eq!(report.maximum_split_pressure_change_pa, 0.0);
    assert!(report.mass_residual_kg.abs() < 1e-12, "{report:?}");
    assert!(report.species_residual_kg.iter().all(|v| v.abs() < 1e-12));
    assert!(report.physical_energy_residual_j.abs() < 5e-6, "{report:?}");
}

#[test]
fn g1_sealed_oxygen_limit_is_applied_once_while_pressure_includes_chemistry() {
    let (mut network, reaction) = fixture(0.0, 0.2, 0.0);
    let before = network.rooms()[0].gas.clone();
    let available = before.layer(Layer::Upper).species_masses()[4] / (2.0 * nasa9::MASSES[4]);
    let calls = Cell::new(0);
    let (report, (), chemistry) = network
        .step_with_pressure_sources(
            Time::new(0.02),
            NetworkControl::default(),
            |_, _, _, _| Ok(()),
            |rooms, dt, _, additional, cancelled| {
                calls.set(calls.get() + 1);
                burn(&reaction, rooms, dt, additional, cancelled)
            },
            || false,
        )
        .unwrap();
    balance(&report);
    assert!(
        calls.get() > 1,
        "the source must participate in pressure iteration"
    );
    assert!((chemistry.actual_extent.value() - available).abs() < 1e-12);
    let after = &network.rooms()[0].gas;
    for s in 0..5 {
        let actual = after.layer(Layer::Upper).species_masses()[s]
            - before.layer(Layer::Upper).species_masses()[s];
        let expected = nasa9::STOICHIOMETRY[s] * nasa9::MASSES[s] * available;
        assert!((actual - expected).abs() < 1e-12);
    }
    assert!(after.layer(Layer::Upper).species_masses()[4] < 1e-14);
    assert!(after.layer(Layer::Upper).species_masses()[0] > 0.0);
    assert!(after.pressure().value() > before.pressure().value());
    assert!(
        (after.physical_internal_energy().value() - before.physical_internal_energy().value())
            .abs()
            < 5e-6
    );
    assert_eq!(after.time(), report.elapsed);
}

#[test]
fn g0_oxygen_removed_by_outflow_is_unavailable_to_every_chemical_trial() {
    let (mut network, reaction) = fixture(0.0, 0.2, 0.5);
    let before = network.rooms()[0].gas.clone();
    let upper = before.layer(Layer::Upper);
    let oxygen_fraction = upper.species_masses()[4] / upper.mass_kg();
    let (report, (), chemistry) = network
        .step_with_pressure_sources(
            Time::new(0.1),
            NetworkControl::default(),
            |_, _, _, _| Ok(()),
            |rooms, dt, _, additional, cancelled| burn(&reaction, rooms, dt, additional, cancelled),
            || false,
        )
        .unwrap();
    balance(&report);
    let remaining_oxygen =
        upper.species_masses()[4] - report.elapsed.value() * 0.5 * oxygen_fraction;
    assert!(
        (chemistry.actual_extent.value() - remaining_oxygen / (2.0 * nasa9::MASSES[4])).abs()
            < 1e-12
    );
    assert!(network.rooms()[0].gas.layer(Layer::Upper).species_masses()[4] < 1e-14);
}

#[test]
fn g1_accepted_chemical_pressure_reproduces_the_actual_opening_mass_transfer() {
    let (mut network, reaction) = fixture(0.2, 0.0, 0.0);
    let before = network.rooms()[0].gas.clone();
    let (report, (), chemistry) = network
        .step_with_pressure_sources(
            Time::new(0.02),
            NetworkControl::default(),
            |rooms, _, rates, _| feed(rooms, rates),
            |rooms, dt, _, additional, cancelled| burn(&reaction, rooms, dt, additional, cancelled),
            || false,
        )
        .unwrap();
    balance(&report);
    let after = &network.rooms()[0].gas;
    let density = before.layer(Layer::Upper).density_kg_m3().unwrap();
    let p = after.pressure().value();
    let inside = TwoLayerColumn::new(
        Length::new(0.0),
        Pressure::new(p),
        Length::new(0.0),
        Density::new(density * p / before.pressure().value()),
        Density::new(density * p / before.pressure().value()),
    )
    .unwrap();
    let outside = TwoLayerColumn::new(
        Length::new(0.0),
        Pressure::new(101325.0),
        Length::new(0.0),
        Density::new(density),
        Density::new(density),
    )
    .unwrap();
    let slabs = stratified::wall_vent_flow(
        &opening(0.2),
        inside,
        outside,
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    let net_outflow: f64 = slabs.slabs().iter().map(|s| s.mass_rate.value()).sum();
    let actual_mass_change =
        after.layer(Layer::Upper).mass_kg() - before.layer(Layer::Upper).mass_kg();
    let expected_mass_change =
        report.elapsed.value() * (0.002 - net_outflow + chemistry.mass_balance_rate_kg_s);
    assert!((actual_mass_change - expected_mass_change).abs() < 2e-8);
    assert!(p > 101325.0 && p - 101325.0 < 10.0, "pressure {p}");
}

#[test]
fn g0_two_room_chemical_pressure_redistributes_gas_with_one_global_balance() {
    let (gas, mechanism, air) = nasa9::source().unwrap();
    let reaction = InventoryLimitedReaction::new(
        &gas,
        mechanism,
        &ReactionId::new(nasa9::REACTION).unwrap(),
        || false,
    )
    .unwrap();
    let rooms = vec![
        room(gas.clone(), air.clone(), 300.0),
        room(gas, air.clone(), 300.0),
    ];
    let initial_mass: f64 = rooms
        .iter()
        .map(|r| r.gas.layer(Layer::Upper).mass_kg())
        .sum();
    let initial_receiver_mass = rooms[1].gas.layer(Layer::Upper).mass_kg();
    let mut network = ZoneNetwork::new(
        rooms,
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: air,
        },
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Room(1),
            opening: opening(0.2),
        }],
        vec![],
        Acceleration::new(0.0),
        || false,
    )
    .unwrap();
    let (report, (), chemistry) = network
        .step_with_pressure_sources(
            Time::new(0.02),
            NetworkControl::default(),
            |rooms, _, rates, _| feed(rooms, rates),
            |rooms, dt, _, additional, cancelled| burn(&reaction, rooms, dt, additional, cancelled),
            || false,
        )
        .unwrap();
    balance(&report);
    let final_mass: f64 = network
        .rooms()
        .iter()
        .map(|r| r.gas.layer(Layer::Upper).mass_kg())
        .sum();
    assert!((final_mass - initial_mass - 0.002 * report.elapsed.value()).abs() < 1e-12);
    assert!(network.rooms()[1].gas.layer(Layer::Upper).mass_kg() > initial_receiver_mass);
    assert!(network.rooms()[0].gas.pressure().value() > network.rooms()[1].gas.pressure().value());
    assert!(network.rooms()[1].gas.pressure().value() > 101325.0);
    assert!(
        (chemistry.actual_extent.value() * nasa9::MASSES[0] - 0.002 * report.elapsed.value()).abs()
            < 1e-16
    );
    assert!(
        network
            .rooms()
            .iter()
            .all(|r| r.gas.time() == report.elapsed)
    );
}

#[test]
fn g4_cancellation_after_chemical_trial_preserves_all_accepted_rooms() {
    let (mut network, reaction) = fixture(0.2, 0.0, 0.0);
    let before = network.rooms().to_vec();
    let cancel = Cell::new(false);
    let result = network.step_with_pressure_sources(
        Time::new(0.02),
        NetworkControl::default(),
        |rooms, _, rates, _| feed(rooms, rates),
        |rooms, dt, _, additional, cancelled| {
            let chemistry = burn(&reaction, rooms, dt, additional, cancelled)?;
            cancel.set(true);
            Ok(chemistry)
        },
        || cancel.get(),
    );
    assert_eq!(result, Err(NetworkError::Cancelled));
    assert_eq!(network.rooms(), before);
}

#[test]
fn g4_retried_duration_rebuilds_preparation_and_publishes_only_accepted_source() {
    let (mut network, reaction) = fixture(0.0, 0.0, 0.0);
    let preparations = Cell::new(0);
    let (report, prepared_duration, chemistry) = network
        .step_with_pressure_sources(
            Time::new(0.08),
            NetworkControl::default(),
            |rooms, dt, rates, _| {
                preparations.set(preparations.get() + 1);
                feed(rooms, rates)?;
                Ok(dt)
            },
            |rooms, dt, prepared, additional, cancelled| {
                assert_eq!(dt, *prepared);
                let chemistry = burn(&reaction, rooms, dt, additional, cancelled)?;
                if dt.value() > 0.02 {
                    return Err(NetworkError::Invalid("authored source duration budget"));
                }
                Ok(chemistry)
            },
            || false,
        )
        .unwrap();
    balance(&report);
    assert_eq!(report.elapsed, Time::new(0.02));
    assert_eq!(prepared_duration, report.elapsed);
    assert_eq!(report.halvings, 2);
    assert_eq!(preparations.get(), 3);
    assert!((chemistry.actual_extent.value() * nasa9::MASSES[0] - 0.002 * 0.02).abs() < 1e-16);
    assert_eq!(network.rooms()[0].gas.time(), report.elapsed);
}

#[test]
fn g0_source_clock_change_and_source_axis_mismatch_refuse_without_publication() {
    for wrong_clock in [false, true] {
        let (mut network, _) = fixture(0.0, 0.0, 0.0);
        let before = network.rooms().to_vec();
        let result = network.step_with_pressure_sources(
            Time::new(0.02),
            NetworkControl {
                maximum_halvings: 0,
                ..NetworkControl::default()
            },
            |_, _, _, _| Ok(()),
            |rooms, dt, _, additional, cancelled| {
                if wrong_clock {
                    rooms[0].gas.advance(
                        dt,
                        [LayerRates::default(); 2],
                        ZoneSolve::default(),
                        cancelled,
                    )?;
                } else {
                    additional.species_kg_s[5] = 1.0;
                }
                Ok(())
            },
            || false,
        );
        assert!(matches!(result, Err(NetworkError::RetryLimit { .. })));
        assert_eq!(network.rooms(), before);
    }
}

#[test]
fn g3_zero_source_keeps_original_pressure_transport_and_accounting_identical() {
    let (mut original, _) = fixture(0.2, 0.0, 0.0);
    let (mut coupled, _) = fixture(0.2, 0.0, 0.0);
    let sources = [[LayerRates::default(); 2]];
    let expected = original
        .step(Time::new(0.02), &sources, NetworkControl::default(), || {
            false
        })
        .unwrap();
    let (actual, (), ()) = coupled
        .step_with_pressure_sources(
            Time::new(0.02),
            NetworkControl::default(),
            |_, _, _, _| Ok(()),
            |_, _, _, _, _| Ok(()),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(coupled.rooms(), original.rooms());
}
