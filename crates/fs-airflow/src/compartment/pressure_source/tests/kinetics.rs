use super::*;
use fs_gas::chemistry::IsochoricKineticsSolve;
use fs_gas::multicomponent::compartment::kinetics::{SplitKineticsReport, ZoneKinetics};
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{MolarEnergyQuantityV1, SpeciesId};

fn fixture() -> (ZoneNetwork, ZoneKinetics) {
    let (gas, _unused_inventory_closure, air) = nasa9::source().unwrap();
    // The NASA calorics are original. This separately constructed forward law
    // is an explicitly synthetic numerical coupling fixture, not methane
    // ignition/combustion data or the inventory adapter's unused placeholder.
    // The original kinetics owner derives the thermodynamic reverse rate.
    let synthetic = ElementaryReactionV1::new(
        ReactionId::new("synthetic-network-finite-rate-test").unwrap(),
        vec![
            (SpeciesId::new("CH4").unwrap(), 1),
            (SpeciesId::new("O2").unwrap(), 2),
        ],
        vec![
            (SpeciesId::new("CO2").unwrap(), 1),
            (SpeciesId::new("H2O").unwrap(), 2),
        ],
        ArrheniusRateV1::new(
            Temperature::new(300.0),
            ReactionProgressRateQuantityV1::new(1e-3),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        gas.source().clone(),
        "synthetic-pressure-coupling-regression-not-physical-kinetics",
        1,
        vec![synthetic],
        32.0 * f64::EPSILON,
    )
    .unwrap();
    let mut reactive: Vec<_> = air.iter().map(|y| 0.98 * y).collect();
    reactive[0] = 0.02;
    let rooms = vec![
        room(gas.clone(), reactive, 300.0),
        room(gas, air.clone(), 300.0),
    ];
    let kinetics = ZoneKinetics::new(
        &rooms[0].gas,
        mechanism,
        IsochoricKineticsSolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 3500.0,
            energy_tolerance_j: 1e-7,
            relative_tolerance: 1e-6,
            initial_step_s: 1e-3,
            maximum_step_s: 0.01,
            ..IsochoricKineticsSolve::default()
        },
        || false,
    )
    .unwrap();
    let network = ZoneNetwork::new(
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
    (network, kinetics)
}

fn inventory(rooms: &[Room]) -> ([f64; 5], f64) {
    let species = std::array::from_fn(|s| {
        rooms
            .iter()
            .map(|room| {
                room.gas.layer(Layer::Lower).species_masses()[s]
                    + room.gas.layer(Layer::Upper).species_masses()[s]
            })
            .sum()
    });
    (
        species,
        rooms
            .iter()
            .map(|room| room.gas.physical_internal_energy().value())
            .sum(),
    )
}

fn run(maximum_dt: f64, retry: bool) -> (Vec<Room>, NetworkStepReport, SplitKineticsReport) {
    let (mut network, mut kinetics) = fixture();
    let before = inventory(network.rooms());
    let receiver_mass = network.rooms()[1].gas.layer(Layer::Upper).mass_kg();
    let calls = Cell::new(0);
    let preparations = Cell::new(0);
    let (report, accepted_duration, chemistry) = network
        .step_with_pressure_sources(
            Time::new(maximum_dt),
            NetworkControl::default(),
            |_, duration, _, _| {
                preparations.set(preparations.get() + 1);
                Ok(duration)
            },
            |rooms, duration, prepared, additional, cancel| {
                calls.set(calls.get() + 1);
                assert_eq!(duration, *prepared);
                let chemistry = kinetics.apply_split(
                    &mut rooms[0].gas,
                    duration,
                    [true; 2],
                    ZoneSolve::default(),
                    &mut *cancel,
                )?;
                for rate in &chemistry.chemistry.chemical_rates {
                    additional.add(rate)?;
                }
                if retry && duration.value() > 0.01 {
                    return Err(NetworkError::Invalid(
                        "authored rejection after complete chemical trial",
                    ));
                }
                Ok(chemistry)
            },
            || false,
        )
        .unwrap();
    balance(&report);
    assert!(
        calls.get() > 2,
        "kinetics must participate in actual pressure iterations"
    );
    assert_eq!(preparations.get(), 1 + report.halvings);
    assert_eq!(accepted_duration, report.elapsed);
    assert!(accepted_duration.value() <= 0.01);
    assert_eq!(chemistry.chemistry.duration, accepted_duration);
    assert!(
        chemistry.chemistry.layers[1]
            .unwrap()
            .reactor
            .accepted_steps
            > 0
    );
    assert_eq!(chemistry.chemistry.layers[0], None);
    let after = inventory(network.rooms());
    let reacted_fuel = before.0[0] - after.0[0];
    assert!(
        reacted_fuel > 0.0 && reacted_fuel < 0.01 * before.0[0],
        "finite progress {reacted_fuel:e}"
    );
    for s in 0..5 {
        let source = report.elapsed.value() * chemistry.chemistry.chemical_rates[1].species_kg_s[s];
        assert!((after.0[s] - before.0[s] - source).abs() < 2e-13);
        assert!((report.external_species_change_kg[s] - source).abs() < 2e-13);
    }
    // Independently reconstruct C, H, N and O atoms from actual global masses.
    for atoms in [
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [4.0, 0.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 2.0, 0.0],
        [0.0, 2.0, 1.0, 0.0, 2.0],
    ] {
        let change: f64 = (0..5)
            .map(|s| atoms[s] * (after.0[s] - before.0[s]) / nasa9::MASSES[s])
            .sum();
        assert!(change.abs() < 2e-11, "element change {change:e}");
    }
    assert!((after.1 - before.1).abs() < 5e-6);
    assert!(report.external_physical_energy_j.abs() < 1e-7);
    assert!(network.rooms()[1].gas.layer(Layer::Upper).mass_kg() > receiver_mass);
    assert!(network.rooms()[0].gas.pressure().value() > network.rooms()[1].gas.pressure().value());
    assert!(network.rooms()[1].gas.pressure().value() > 101325.0);
    assert!(
        network
            .rooms()
            .iter()
            .all(|r| r.gas.time() == report.elapsed)
    );
    (network.rooms().to_vec(), report, chemistry)
}

#[test]
fn finite_rate_pressure_trials_close_actual_elements_energy_and_replay_after_retry() {
    let (retried, report, chemistry) = run(0.02, true);
    assert!(report.halvings > 0);
    // Use the actual accepted duration, including any further pressure/kinetic
    // retries selected by the original owners, rather than assuming one half.
    let (direct, direct_report, direct_chemistry) = run(report.elapsed.value(), false);
    assert_eq!(direct_report.halvings, 0);
    assert_eq!(retried, direct);
    assert_eq!(chemistry, direct_chemistry);
    assert_eq!(
        report.external_species_change_kg,
        direct_report.external_species_change_kg
    );
    assert_eq!(
        report.external_physical_energy_j,
        direct_report.external_physical_energy_j
    );
    assert_eq!(
        report.maximum_pressure_residual_pa,
        direct_report.maximum_pressure_residual_pa
    );
}
