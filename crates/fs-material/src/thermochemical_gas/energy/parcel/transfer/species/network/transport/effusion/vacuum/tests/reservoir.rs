//! Original molecular sources, gross counterflow and one finite caloric solve.
use super::super::super::reservoir as open;
use super::*;

fn port<'a>(
    source: &'a ClosedGasParcel,
    law: &'a EffusiveInterface,
) -> open::EffusionReservoirLink<'a> {
    open::EffusionReservoirLink {
        id: 9,
        parcel: 0,
        reservoir: source,
        interface: law,
        area_m2: AREA,
    }
}

#[test]
fn equal_sources_keep_both_gross_beams_at_zero_net_mass_and_heat() {
    let gas = pure(2.5, 0.0, 400.0, 100_000.0);
    let source = gas.clone();
    let law = interface(gas.model(), support(), 1.0);
    let gases = BTreeMap::from([(0, gas)]);
    let saved = gases.clone();
    let result = open::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[],
        &[port(&source, &law)],
        0.1,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let beam = &result.step().reservoirs[0].flux.species[0];
    assert!(beam.left_to_right.mass_rate_kg_s > 0.0);
    assert_eq!(
        beam.left_to_right.mass_rate_kg_s,
        beam.right_to_left.mass_rate_kg_s
    );
    let row = &result.step().exchange.boundaries[0];
    assert!(row.transfer.incoming_mass_kg > 0.0);
    assert_eq!(row.transfer.incoming_mass_kg, row.transfer.outgoing_mass_kg);
    assert_eq!(row.transfer.net_mass_into_gas_kg, 0.0);
    assert_eq!(row.transfer.relative_energy_into_gas_j, 0.0);
    near(result.parcels()[&0].mass_kg(), gases[&0].mass_kg(), 1e-18);
    near(
        result.parcels()[&0].state().state().temperature,
        400.0,
        1e-8,
    );
    assert_eq!(gases, saved);
    assert_eq!(source, gases[&0]);
}

#[test]
fn hot_counterflow_transfers_heat_without_losing_it_to_formation_datums() {
    let mut temperatures = Vec::new();
    for datum in [0.0, 1e12] {
        let gas = pure(2.5, datum, 400.0, 100_000.0);
        let source = pure(2.5, datum, 1600.0, 200_000.0);
        let saved = source.clone();
        let law = interface(gas.model(), support(), 1.0);
        let gases = BTreeMap::from([(0, gas)]);
        let result = open::step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[],
            &[port(&source, &law)],
            0.1,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap();
        let row = &result.step().exchange.boundaries[0].transfer;
        let dm = row.net_mass_into_gas_kg;
        let energy = 0.5 * dm * R * 400.0 / M + 2.0 * row.incoming_mass_kg * R * 1200.0 / M;
        let expected = 400.0 + energy / ((gases[&0].mass_kg() + dm) * 1.5 * R / M);
        let actual = result.parcels()[&0].state().state().temperature;
        near(actual, expected, 2e-8);
        assert!(row.relative_energy_into_gas_j > 0.0);
        assert!(actual > 400.0);
        assert!(
            result.step().exchange.total_energy_residual_j.abs()
                <= controls().exchange.energy_tolerance_j
        );
        temperatures.push(actual);
        assert_eq!(source, saved);
    }
    near(temperatures[0], temperatures[1], 2e-8);
}

#[test]
fn supply_and_vacuum_share_gross_withdrawal_and_original_caloric_update() {
    let gas = pure(2.5, 0.0, 400.0, 100_000.0);
    let source = pure(2.5, 0.0, 400.0, 300_000.0);
    let law = interface(gas.model(), support(), 1.0);
    let gases = BTreeMap::from([(0, gas)]);
    let vents = [outlet(&law, AREA)];
    let ports = [port(&source, &law)];
    let result = open::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &vents,
        &ports,
        1.0,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let delta = result
        .step()
        .exchange
        .boundaries
        .iter()
        .map(|r| r.transfer.net_mass_into_gas_kg)
        .sum::<f64>();
    assert!(delta > 0.0);
    let expected = 400.0 * (1.0 + delta / (3.0 * (gases[&0].mass_kg() + delta)));
    near(
        result.parcels()[&0].state().state().temperature,
        expected,
        2e-8,
    );
    let mut restrictive = controls();
    restrictive.exchange.maximum_species_fraction = 0.02;
    assert!(matches!(
        open::step_candidate(
            &gases,
            GasTransportPaths::default(),
            &vents,
            &ports,
            1.0,
            &endpoints(&gases),
            restrictive,
            || false
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        ))
    ));
    let legacy = step_candidate(
        &gases,
        GasTransportPaths::default(),
        &vents,
        0.1,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let combined = open::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &vents,
        &[],
        0.1,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(legacy.parcels(), combined.parcels());
    assert_eq!(legacy.step().outlets, combined.step().vacuum_outlets);
    assert_eq!(legacy.step().exchange, combined.step().exchange);
}

#[test]
fn source_can_add_an_original_species_absent_from_the_finite_vessel() {
    let gas = pure(2.5, 0.0, 400.0, 100_000.0);
    let supply_model = model(&[("A", M, 2.5), ("B", 0.004, 2.5)], vec![0.5, 0.5], 0.0);
    let source = ClosedGasParcel::new(supply_model, point(500.0, 200_000.0), VOLUME).unwrap();
    let saved = source.clone();
    let law = interface(source.model(), support(), 1.0);
    let gases = BTreeMap::from([(0, gas)]);
    let result = open::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[],
        &[port(&source, &law)],
        0.1,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let b = result
        .step()
        .exchange
        .boundaries
        .iter()
        .find(|r| r.transfer.species.as_str() == "B")
        .unwrap();
    assert!(b.transfer.incoming_mass_kg > 0.0);
    assert_eq!(b.transfer.outgoing_mass_kg, 0.0);
    let mixture = result.parcels()[&0].model().thermodynamics();
    assert_eq!(mixture.components().len(), 2);
    assert!(mixture.mole_composition().fractions()[1] > 0.0);
    let before = gases[&0].mass_kg() * gases[&0].state().specific_internal_energy_j_kg();
    let after = result.parcels()[&0].mass_kg()
        * result.parcels()[&0].state().specific_internal_energy_j_kg();
    let supplied = result
        .step()
        .exchange
        .boundaries
        .iter()
        .map(|r| r.energy_into_gas_j)
        .sum::<f64>();
    near(after - before, supplied, 2e-8);
    assert_eq!(source, saved);
}

#[test]
fn source_endpoint_identity_limits_and_cancellation_refuse_atomically() {
    let gases = BTreeMap::from([(0, pure(2.5, 0.0, 400.0, 100_000.0))]);
    let saved = gases.clone();
    let source = pure(2.5, 0.0, 200.0, 100_000.0);
    let narrow = interface(
        gases[&0].model(),
        EffusiveSupport {
            pressure_pa: (99_999.0, 100_001.0),
            ..support()
        },
        1.0,
    );
    assert!(
        open::step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[],
            &[port(&source, &narrow)],
            0.1,
            &endpoints(&gases),
            controls(),
            || false
        )
        .is_err()
    );
    let law = interface(gases[&0].model(), support(), 1.0);
    let foreign = pure(3.5, 0.0, 400.0, 100_000.0);
    assert!(
        open::step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[],
            &[port(&foreign, &law)],
            0.1,
            &endpoints(&gases),
            controls(),
            || false
        )
        .is_err()
    );
    let closed = open::EffusionReservoirLink {
        area_m2: 0.0,
        ..port(&foreign, &law)
    };
    let idle = open::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[],
        &[closed],
        0.1,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(idle.parcels(), &gases);
    let duplicate = open::EffusionReservoirLink {
        id: 7,
        ..port(&source, &law)
    };
    assert!(
        open::step_candidate(
            &gases,
            GasTransportPaths::default(),
            &[outlet(&law, AREA)],
            &[duplicate],
            0.1,
            &endpoints(&gases),
            controls(),
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    open::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[],
        &[port(&source, &law)],
        0.1,
        &endpoints(&gases),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls - 1, polls] {
        let mut count = 0;
        assert!(
            open::step_candidate(
                &gases,
                GasTransportPaths::default(),
                &[],
                &[port(&source, &law)],
                0.1,
                &endpoints(&gases),
                controls(),
                || {
                    count += 1;
                    count == stop
                }
            )
            .is_err()
        );
        assert_eq!(gases, saved);
    }
}
