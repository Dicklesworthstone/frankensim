//! End-to-end source-backed rates, multiport reconstruction and time refinement.
use super::*;
use crate::thermochemical_gas::tests::{model, point};
use crate::thermochemical_gas::energy::parcel::transfer::{GasMassTransferSolve,
    nozzle::GasNozzleStepSolve};

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn controls() -> GasSnapshotSolve {
    GasSnapshotSolve { maximum_withdrawn_fraction: 0.2, mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }
}
fn solves(p: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    p.keys().map(|&id| (id, energy())).collect()
}
fn gas(x: f64, t: f64, p: f64, slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], slope), point(t, p), 0.001).unwrap()
}
fn link(id: u64, a: u64, b: u64, area_m2: f64) -> SnapshotNozzleLink {
    SnapshotNozzleLink { id, parcels: [a, b], nozzle: GasNozzle { area_m2, discharge_coefficient: 0.8 },
        nozzle_solve: energy() }
}
fn initial() -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([(0, gas(0.75, 500.0, 300000.0, 0.001)),
        (1, gas(0.5, 450.0, 200000.0, 0.001)), (2, gas(0.25, 400.0, 100000.0, 0.001))])
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn actual_rates_are_from_the_original_snapshot_and_link_ids_do_not_change_physics() {
    let parcels = initial();
    let links = [link(10, 0, 1, 1e-6), link(20, 1, 2, 1e-6), link(30, 0, 2, 0.5e-6)];
    let trial = step_candidate(&parcels, &links, 0.0002, &solves(&parcels), controls(), || false).unwrap();
    for report in &trial.step().links {
        let source_link = links.iter().find(|l| l.id == report.id).unwrap();
        let original = parcels[&report.donor].nozzle_flow_to(parcels[&report.receiver].state().state().pressure,
            source_link.nozzle, energy(), || false).unwrap();
        assert_eq!(report.flow, original);
        assert_eq!(report.transferred_mass_kg, original.mass_flow_kg_s*0.0002);
    }
    let renamed = [link(2, 2, 0, 0.5e-6), link(1, 2, 1, 1e-6), link(3, 1, 0, 1e-6)];
    let other = step_candidate(&parcels, &renamed, 0.0002, &solves(&parcels), controls(), || false).unwrap();
    assert_eq!(trial.parcels(), other.parcels());
    assert_eq!(trial.step().exchange, other.step().exchange);
}

#[test]
fn equal_fanout_receivers_stay_symmetric_without_sequential_pressure_bias() {
    let receiver = gas(0.25, 400.0, 100000.0, 0.001);
    let parcels = BTreeMap::from([(0, gas(0.75, 500.0, 300000.0, 0.001)),
        (1, receiver.clone()), (2, receiver)]);
    let links = [link(1, 0, 1, 1e-6), link(2, 0, 2, 1e-6)];
    let trial = step_candidate(&parcels, &links, 0.001, &solves(&parcels), controls(), || false).unwrap();
    assert_eq!(trial.step().links[0].transferred_mass_kg, trial.step().links[1].transferred_mass_kg);
    assert_eq!(trial.parcels()[&1], trial.parcels()[&2]);
    let fraction = trial.step().exchange.parcels[0].outgoing_mass_kg / parcels[&0].mass_kg();
    let strict = GasSnapshotSolve { maximum_withdrawn_fraction: 0.75*fraction, ..controls() };
    assert!(matches!(step_candidate(&parcels, &links, 0.001, &solves(&parcels), strict, || false),
        Err(Error::Transfer(GasSnapshotError::Withdrawal { id: 0, .. }))));
}

#[test]
fn joint_pressure_crossing_refuses_even_when_each_isolated_pair_would_pass() {
    let donor = gas(0.75, 400.0, 110000.0, 0.0);
    let parcels = BTreeMap::from([(0, donor.clone()), (1, donor), (2, gas(0.75, 400.0, 100000.0, 0.0))]);
    let links = [link(1, 0, 2, 1e-6), link(2, 1, 2, 1e-6)];
    let g = parcels[&0].state().state();
    // Each donor drops 4 kPa: an isolated pair ends at 106/104 kPa,
    // but the combined receiver rises to 108 kPa, reversing both paths.
    let dm = 4000.0*0.001/(g.gamma*g.specific_gas_constant*g.temperature);
    let flow = parcels[&0].nozzle_flow_to(100000.0, links[0].nozzle, energy(), || false).unwrap();
    let dt = dm/flow.mass_flow_kg_s;
    let pair = GasNozzleStepSolve { nozzle: energy(), maximum_donor_fraction: 0.2,
        transfer: GasMassTransferSolve { donor: energy(), receiver: energy(),
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } };
    assert!(parcels[&0].through_nozzle_to(&parcels[&2], links[0].nozzle, dt, pair, || false).is_ok());
    let saved = parcels.clone();
    assert!(matches!(step_candidate(&parcels, &links, dt, &solves(&parcels), controls(), || false),
        Err(Error::EquilibriumCrossed { .. })));
    assert_eq!(parcels, saved);
}

#[test]
fn source_data_changes_flow_and_updated_acoustics_without_changing_graph_code() {
    let mut outcomes = Vec::new();
    for x in [0.75, 0.25] {
        let parcels = BTreeMap::from([(0, gas(x, 500.0, 300000.0, 0.001)),
            (1, gas(0.5, 400.0, 100000.0, 0.001))]);
        let trial = step_candidate(&parcels, &[link(7, 0, 1, 1e-6)], 0.001,
            &solves(&parcels), controls(), || false).unwrap();
        assert!(trial.step().links[0].flow.choked);
        let state = trial.parcels()[&1].state().state();
        near(state.sound_speed.powi(2), state.gamma*state.pressure/state.density, 1e-12);
        outcomes.push((trial.step().links[0].flow.mass_flow_kg_s, state.sound_speed));
    }
    assert_ne!(outcomes[0].0, outcomes[1].0);
    assert_ne!(outcomes[0].1, outcomes[1].1);
}

#[test]
fn graph_trajectory_refines_with_global_conservation() {
    let original = initial();
    let links = [link(10, 0, 1, 1e-6), link(20, 1, 2, 1e-6)];
    let mut pressures = Vec::new();
    for count in [32, 64, 128] {
        let mut parcels = original.clone();
        for _ in 0..count {
            let trial = step_candidate(&parcels, &links, 0.02/f64::from(count),
                &solves(&parcels), controls(), || false).unwrap();
            assert!(trial.step().exchange.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
            assert!(trial.step().exchange.maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
            parcels = trial.into_parts().0;
        }
        pressures.push(parcels[&2].state().state().pressure);
        let drift: f64 = parcels.iter().map(|(id,p)| super::super::super::super::stored_change(&original[id],p)).sum();
        assert!(drift.abs() < 1e-6);
    }
    assert!((pressures[0]-pressures[1]).abs() > 1.8*(pressures[1]-pressures[2]).abs());
}

#[test]
fn closed_paths_preflight_failures_and_late_cancellation_return_no_partial_graph() {
    let parcels = initial();
    let saved = parcels.clone();
    let closed = [link(1, 0, 1, 0.0), link(2, 1, 2, 0.0)];
    let trial = step_candidate(&parcels, &closed, 0.001, &solves(&parcels), controls(), || false).unwrap();
    assert_eq!(trial.parcels(), &parcels);
    assert!(trial.step().links.iter().all(|s| s.transferred_mass_kg == 0.0));
    for links in [[link(1, 0, 1, 1e-6), link(1, 1, 2, 1e-6)],
        [link(1, 0, 9, 1e-6), link(2, 1, 2, 1e-6)],
        [link(1, 0, 0, 1e-6), link(2, 1, 2, 1e-6)]] {
        assert!(step_candidate(&parcels, &links, 0.001, &solves(&parcels), controls(), || false).is_err());
    }
    let links = [link(1, 0, 1, 1e-6), link(2, 1, 2, 1e-6)];
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, 1000.0] {
        assert!(step_candidate(&parcels, &links, dt, &solves(&parcels), controls(), || false).is_err());
    }
    let mut total = 0;
    step_candidate(&parcels, &links, 0.001, &solves(&parcels), controls(), || { total+=1; false }).unwrap();
    for stop in [1, total/2, total] {
        let mut calls = 0;
        assert_eq!(step_candidate(&parcels, &links, 0.001, &solves(&parcels), controls(),
            || { calls+=1; calls==stop }), Err(Error::Cancelled));
        assert_eq!(parcels, saved);
    }
}
