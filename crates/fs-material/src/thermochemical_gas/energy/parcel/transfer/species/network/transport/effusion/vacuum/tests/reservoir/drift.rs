//! Same original species/caloric owner with qualified moving incidence.
use super::*;
use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::{
    drift::EffusiveDriftInterface, momentum::EffusiveMomentumInterface,
};

fn qualified(model: &ThermochemicalGasModel, maximum: f64) -> EffusiveDriftInterface {
    EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            interface(model, support(), 1.0),
            "synthetic uncollimated nondeflecting angular source",
        )
        .unwrap(),
        maximum,
        "synthetic collisionless drift qualification over original T/p support",
    )
    .unwrap()
}
fn path<'a>(
    source: &'a ClosedGasParcel,
    law: &'a EffusiveDriftInterface,
    velocity: [f64; 3],
) -> open::EffusionDriftingReservoirLink<'a> {
    open::EffusionDriftingReservoirLink {
        id: 9,
        parcel: 0,
        reservoir: source,
        interface: law,
        area_m2: AREA,
        relative_bulk_velocity_m_s: velocity,
        outward_unit_normal: [1.0, 0.0, 0.0],
    }
}
fn advance(
    gases: &BTreeMap<u64, ClosedGasParcel>,
    links: &[open::EffusionDriftingReservoirLink<'_>],
    dt: f64,
) -> open::ReservoirEffusionCandidate {
    open::step_candidate_with_drift(
        gases,
        GasTransportPaths::default(),
        &[],
        &[],
        links,
        dt,
        &endpoints(gases),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn zero_drift_reuses_the_original_comoving_caloric_and_species_transaction() {
    let gas = pure(2.5, 0.0, 400.0, 100_000.0);
    let source = pure(2.5, 0.0, 800.0, 240_000.0);
    let law = qualified(gas.model(), 2000.0);
    let gases = BTreeMap::from([(0, gas)]);
    let expected = open::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[],
        &[port(&source, law.momentum_interface().interface())],
        0.1,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let actual = advance(&gases, &[path(&source, &law, [0.0; 3])], 0.1);
    assert_eq!(actual.parcels(), expected.parcels());
    assert_eq!(actual.step().exchange, expected.step().exchange);
    assert!(actual.step().reservoirs.is_empty());
    assert_eq!(actual.step().drifting_reservoirs.len(), 1);
}

#[test]
fn tangential_drift_heats_once_without_changing_mass_or_losing_large_datum_precision() {
    let mut temperatures = Vec::new();
    for datum in [0.0, 1e12] {
        let gas = pure(2.5, datum, 400.0, 100_000.0);
        let source = gas.clone();
        let law = qualified(gas.model(), 2000.0);
        let gases = BTreeMap::from([(0, gas)]);
        let actual = advance(&gases, &[path(&source, &law, [0.0, 200.0, 0.0])], 0.1);
        let row = &actual.step().exchange.boundaries[0];
        assert_eq!(row.transfer.net_mass_into_gas_kg, 0.0);
        assert_eq!(row.transfer.incoming_mass_kg, row.transfer.outgoing_mass_kg);
        let heat = row.transfer.incoming_mass_kg * 0.5 * 200.0_f64.powi(2);
        near(row.transfer.relative_energy_into_gas_j, heat, 1e-11);
        let expected = 400.0 + heat / (gases[&0].mass_kg() * 1.5 * R / M);
        let temperature = actual.parcels()[&0].state().state().temperature;
        near(temperature, expected, 2e-8);
        temperatures.push(temperature);
        assert_eq!(actual.parcels()[&0].mass_kg(), gases[&0].mass_kg());
        assert_eq!(source, gases[&0]);
    }
    near(temperatures[0], temperatures[1], 2e-8);
}

#[test]
fn vacuum_comoving_and_drifting_ports_spend_one_original_gross_allowance() {
    let gas = pure(2.5, 0.0, 400.0, 100_000.0);
    let source = pure(2.5, 0.0, 400.0, 300_000.0);
    let law = qualified(gas.model(), 2000.0);
    let gases = BTreeMap::from([(0, gas)]);
    let saved = gases.clone();
    let base = law.momentum_interface().interface();
    let vents = [outlet(base, AREA)];
    let comoving = [port(&source, base)];
    let drifting = [open::EffusionDriftingReservoirLink {
        id: 11,
        ..path(&source, &law, [0.0, 100.0, 0.0])
    }];
    let q = gases[&0]
        .effusive_flux_to_vacuum(base, AREA, || false)
        .unwrap()
        .species[0]
        .outward
        .mass_rate_kg_s;
    let mut bounded = controls();
    bounded.exchange.maximum_species_fraction = 2.5 * q / gases[&0].mass_kg();
    assert!(
        open::step_candidate_with_drift(
            &gases,
            GasTransportPaths::default(),
            &vents,
            &[],
            &drifting,
            1.0,
            &endpoints(&gases),
            bounded,
            || false
        )
        .is_ok()
    );
    assert!(matches!(
        open::step_candidate_with_drift(
            &gases,
            GasTransportPaths::default(),
            &vents,
            &comoving,
            &drifting,
            1.0,
            &endpoints(&gases),
            bounded,
            || false
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        ))
    ));
    assert_eq!(gases, saved);
}

#[test]
fn tiny_normal_drift_retains_signed_mass_and_a_far_tail_never_becomes_closed() {
    let gas = pure(2.5, 0.0, 400.0, 100_000.0);
    let source = gas.clone();
    let law = qualified(gas.model(), 20_000.0);
    let gases = BTreeMap::from([(0, gas)]);
    let result = advance(&gases, &[path(&source, &law, [-1e-14, 0.0, 0.0])], 0.1);
    let net = result.step().exchange.boundaries[0]
        .transfer
        .net_mass_into_gas_kg;
    let expected = 0.5 * (source.mass_kg() / source.volume_m3()) * 1e-14 * AREA * 0.1;
    assert!(net > 0.0);
    near(net, expected, expected * 1e-12);

    let catalog = model(&[("A", M, 2.5), ("B", 0.004, 2.5)], vec![0.5, 0.5], 0.0);
    let source = ClosedGasParcel::new(
        catalog
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap(),
            )
            .unwrap(),
        point(400.0, 100_000.0),
        VOLUME,
    )
    .unwrap();
    let original = interface(&catalog, support(), 1.0);
    let transmission = catalog
        .thermodynamics()
        .components()
        .iter()
        .map(|s| {
            EffusiveSpeciesTransmission::try_new(
                s.species().clone(),
                if s.species().as_str() == "B" {
                    1.0
                } else {
                    0.0
                },
                "synthetic B-selective nondeflecting aperture",
            )
            .unwrap()
        })
        .collect();
    let law = EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            EffusiveInterface::try_new(catalog, transmission, support(), original.regime().clone())
                .unwrap(),
            "synthetic uncollimated angular law",
        )
        .unwrap(),
        20_000.0,
        "synthetic positive-tail drift qualification",
    )
    .unwrap();
    let speed = 12.0 * (R * 400.0 / 0.004).sqrt();
    let result = advance(&gases, &[path(&source, &law, [speed, 0.0, 0.0])], 0.1);
    let row = &result.step().exchange.boundaries[0].transfer;
    assert_eq!(row.species.as_str(), "B");
    assert_eq!(row.outgoing_mass_kg, 0.0);
    assert!(row.incoming_mass_kg > 0.0);
    assert_eq!(row.net_mass_into_gas_kg, row.incoming_mass_kg);
    assert!(result.parcels()[&0].species_mass_kg(&row.species).unwrap() > 0.0);
}

#[test]
fn drift_endpoint_source_and_last_cancellation_fail_atomically() {
    let gases = BTreeMap::from([(0, pure(2.5, 0.0, 400.0, 100_000.0))]);
    let saved = gases.clone();
    let source = gases[&0].clone();
    let narrow = EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            interface(
                source.model(),
                EffusiveSupport {
                    temperature_k: (200.0, 400.0001),
                    ..support()
                },
                1.0,
            ),
            "synthetic narrow original angular source",
        )
        .unwrap(),
        1000.0,
        "synthetic drift support",
    )
    .unwrap();
    assert!(
        open::step_candidate_with_drift(
            &gases,
            GasTransportPaths::default(),
            &[],
            &[],
            &[path(&source, &narrow, [0.0, 200.0, 0.0])],
            0.1,
            &endpoints(&gases),
            controls(),
            || false
        )
        .is_err()
    );
    let law = qualified(source.model(), 1000.0);
    let foreign = pure(3.5, 0.0, 400.0, 100_000.0);
    assert!(
        open::step_candidate_with_drift(
            &gases,
            GasTransportPaths::default(),
            &[],
            &[],
            &[path(&foreign, &law, [0.0; 3])],
            0.1,
            &endpoints(&gases),
            controls(),
            || false
        )
        .is_err()
    );
    let closed = open::EffusionDriftingReservoirLink {
        area_m2: 0.0,
        ..path(&foreign, &law, [2000.0, 0.0, 0.0])
    };
    assert_eq!(advance(&gases, &[closed], 0.1).parcels(), &gases);
    let links = [path(&source, &law, [-100.0, 200.0, 0.0])];
    let mut polls = 0;
    open::step_candidate_with_drift(
        &gases,
        GasTransportPaths::default(),
        &[],
        &[],
        &links,
        0.1,
        &endpoints(&gases),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(
        open::step_candidate_with_drift(
            &gases,
            GasTransportPaths::default(),
            &[],
            &[],
            &links,
            0.1,
            &endpoints(&gases),
            controls(),
            || {
                count += 1;
                count == polls
            }
        )
        .is_err()
    );
    assert_eq!(gases, saved);
}
