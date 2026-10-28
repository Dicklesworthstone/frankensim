//! Real source-resolved endpoints; fixtures are synthetic, not measured gases.
use super::*;
use super::super::{SpeciesId, stored_change};
use crate::thermochemical_gas::tests::{model, point, species, transport};
use crate::thermochemical_gas::ThermochemicalGasModel;
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn controls() -> GasSnapshotSolve {
    GasSnapshotSolve { maximum_withdrawn_fraction: 0.25,
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }
}
fn solves(parcels: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    parcels.keys().map(|&id| (id, energy())).collect()
}
fn gas(x: f64, t: f64, p: f64, slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], slope),
        point(t, p).with("local-source-axis", p).unwrap(), 0.001).unwrap()
}
fn pure(name: &str, molar_mass: f64, cp_r: f64, t: f64, p: f64) -> ClosedGasParcel {
    let source = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species(name, molar_mass, cp_r, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport(name, 1.7e-5, 0.025)], None).unwrap();
    ClosedGasParcel::new(source, point(t, p), 0.001).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1.0), "{a} != {b}");
}
fn row(donor: u64, receiver: u64, mass_kg: f64) -> GasSnapshotTransfer {
    GasSnapshotTransfer { donor, receiver, mass_kg }
}

#[test]
fn g3_bulk_network_preserves_nonlinear_heating_under_independent_species_datums() {
    use super::super::tests::reference_parcel;
    let graph = |reference| BTreeMap::from([
        (0, reference_parcel(0.75, 0.001, 500.0, 250000.0, reference)),
        (1, reference_parcel(0.25, 0.001, 350.0, 100000.0, reference)),
        (2, reference_parcel(0.5, 0.001, 400.0, 150000.0, reference)),
    ]);
    let ordinary = graph([0.0; 2]);
    let shifted = graph([1e24, -2e24]);
    let dm = 0.02 * ordinary[&1].mass_kg();
    let rows = [row(0, 1, dm), row(1, 2, dm), row(2, 0, dm)];
    let expected = transfer_candidate(&ordinary, &rows, &solves(&ordinary), controls(), || false).unwrap();
    let actual = transfer_candidate(&shifted, &rows, &solves(&shifted), controls(), || false).unwrap();
    for id in [0, 1, 2] {
        assert_eq!(actual.parcels()[&id].state().state(), expected.parcels()[&id].state().state());
        assert_eq!(actual.parcels()[&id].mass_kg(), expected.parcels()[&id].mass_kg());
        assert_ne!(actual.parcels()[&id].state().state().temperature,
            shifted[&id].state().state().temperature);
        assert!(actual.step().parcels[id as usize].energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
    assert_eq!(actual.step().total_energy_residual_j, expected.step().total_energy_residual_j);
    assert!(actual.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    let mut reversed = rows; reversed.reverse();
    assert_eq!(actual, transfer_candidate(&shifted, &reversed, &solves(&shifted), controls(), || false).unwrap());
}

#[test]
fn bulk_edge_capacity_is_not_reduced_by_expanding_multiple_species() {
    let source = gas(0.75, 400.0, 100000.0, 0.001);
    let parcels: BTreeMap<_, _> = (0..65).map(|id| (id, source.clone())).collect();
    let dm = source.mass_kg() * 1e-4;
    let mut rows = Vec::new();
    for donor in 0..65 {
        for offset in 1..=33 { rows.push(row(donor, (donor + offset) % 65, dm)); }
    }
    assert!(rows.len() < MAX_SNAPSHOT_TRANSFERS);
    assert!(2 * rows.len() > MAX_SNAPSHOT_TRANSFERS);
    let candidate = transfer_candidate(&parcels, &rows, &solves(&parcels), controls(), || false).unwrap();
    assert_eq!(candidate.parcels(), &parcels);
    assert_eq!(candidate.step().fluxes.len(), rows.len());
    assert!(candidate.step().parcels.iter().all(|p| p.iterations == 0));
}

#[test]
fn one_link_reuses_the_original_pair_material_and_energy_owners() {
    for slope in [0.0, 0.001] {
        let parcels = BTreeMap::from([(3, gas(0.75, 500.0, 200000.0, slope)),
            (9, gas(0.25, 400.0, 100000.0, slope))]);
        let dm = 0.05*parcels[&3].mass_kg();
        let original = parcels[&3].transfer_mass_to(&parcels[&9], dm, super::super::GasMassTransferSolve {
            donor: energy(), receiver: energy(), mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8,
        }, || false).unwrap();
        let trial = transfer_candidate(&parcels, &[row(3, 9, dm)], &solves(&parcels), controls(), || false).unwrap();
        for (id, reference) in [(3, original.donor()), (9, original.receiver())] {
            let actual = &trial.parcels()[&id];
            assert_eq!(actual.model(), reference.model());
            near(actual.mass_kg(), reference.mass_kg(), 1e-14);
            near(actual.state().state().temperature, reference.state().state().temperature, 1e-10);
            near(actual.state().state().pressure, reference.state().state().pressure, 1e-10);
            assert_eq!(actual.point().axes().get("local-source-axis"), parcels[&id].point().axes().get("local-source-axis"));
        }
        near(trial.step().fluxes[0].enthalpy_j, original.step().transferred_enthalpy_j, 1e-12);
        assert!(trial.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
}

#[test]
fn incoming_species_do_not_jump_through_an_outlet_in_the_same_interval() {
    let parcels = BTreeMap::from([(0, pure("A", 0.028, 3.5, 500.0, 250000.0)),
        (1, pure("B", 0.004, 2.5, 450.0, 180000.0)),
        (2, pure("B", 0.004, 2.5, 400.0, 100000.0))]);
    let dm = 0.02*parcels[&1].mass_kg();
    let transfers = [row(0, 1, dm), row(1, 2, dm)];
    let trial = transfer_candidate(&parcels, &transfers, &solves(&parcels), controls(), || false).unwrap();
    assert_eq!(trial.parcels()[&1].model().thermodynamics().components().len(), 2);
    assert_eq!(trial.parcels()[&2].model().thermodynamics().components().len(), 1);
    assert_eq!(trial.step().fluxes[1].species[0].species, SpeciesId::new("B").unwrap());
    let inv = inventory(trial.parcels()[&1].model(), trial.parcels()[&1].mass_kg()).unwrap();
    near(inv[&SpeciesId::new("A").unwrap()], dm, 1e-14);
    near(inv[&SpeciesId::new("B").unwrap()], parcels[&1].mass_kg()-dm, 1e-14);
    // The second physical interval CAN carry A, using its then-accepted source.
    let next = transfer_candidate(trial.parcels(), &transfers, &solves(trial.parcels()), controls(), || false).unwrap();
    assert_eq!(next.parcels()[&2].model().thermodynamics().components().len(), 2);
    // Falsifier: immediate sequential pair updates do carry A across both links.
    let first = parcels[&0].transfer_mass_to(&parcels[&1], dm, super::super::GasMassTransferSolve {
        donor: energy(), receiver: energy(), mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8,
    }, || false).unwrap();
    let split = first.receiver().transfer_mass_to(&parcels[&2], dm, super::super::GasMassTransferSolve {
        donor: energy(), receiver: energy(), mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8,
    }, || false).unwrap();
    assert_eq!(split.receiver().model().thermodynamics().components().len(), 2);
}

#[test]
fn multiple_inlets_outlets_and_parallel_rows_are_order_independent() {
    let parcels = BTreeMap::from([(0, gas(0.75, 500.0, 250000.0, 0.001)),
        (1, gas(0.25, 450.0, 180000.0, 0.001)), (2, gas(0.5, 400.0, 120000.0, 0.001)),
        (3, gas(0.5, 350.0, 90000.0, 0.001))]);
    let dm = 0.01*parcels[&3].mass_kg();
    let rows = vec![row(0, 2, dm), row(1, 2, 2.0*dm), row(2, 3, dm),
        row(0, 2, 0.5*dm), row(3, 1, dm), row(0, 3, 0.0)];
    let expected = transfer_candidate(&parcels, &rows, &solves(&parcels), controls(), || false).unwrap();
    for rotation in 0..rows.len() {
        let mut reordered = rows.clone(); reordered.rotate_left(rotation); reordered.reverse();
        assert_eq!(transfer_candidate(&parcels, &reordered, &solves(&parcels), controls(), || false).unwrap(), expected);
    }
    assert_eq!(expected.step().fluxes.len(), 4);
    for account in &expected.step().parcels {
        assert!(account.energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(account.maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    }
    // Continued thermal/mechanical consumers use the updated, source-bearing medium.
    let before = &expected.parcels()[&2];
    let moved = before.midpoint_heat_volume(0.00101, 2.0, energy(), || false).unwrap();
    let load = moved.discrete_mechanical_load(0.00101, &[0.01], 100000.0).unwrap();
    near(load.generalized_forces().unwrap()[0], 0.01*(moved.midpoint().state().pressure-100000.0), 1e-12);
    near(stored_change(before, moved.endpoint())+moved.step().work_by_gas_j, 2.0, 1e-8);
}

#[test]
fn counterflow_can_change_energy_and_composition_without_changing_net_mass() {
    let parcels = BTreeMap::from([(0, gas(0.75, 500.0, 150000.0, 0.001)),
        (1, gas(0.25, 350.0, 100000.0, 0.001))]);
    let dm = 0.01*parcels[&1].mass_kg();
    let trial = transfer_candidate(&parcels, &[row(0, 1, dm), row(1, 0, dm)],
        &solves(&parcels), controls(), || false).unwrap();
    let exchanged = dm*(parcels[&1].state().specific_enthalpy_j_kg()-parcels[&0].state().specific_enthalpy_j_kg());
    near(trial.parcels()[&0].mass_kg(), parcels[&0].mass_kg(), 1e-14);
    near(stored_change(&parcels[&0], &trial.parcels()[&0]), exchanged, 1e-8);
    assert_ne!(trial.parcels()[&0].model().thermodynamics().mole_composition(), parcels[&0].model().thermodynamics().mole_composition());
    assert_ne!(trial.parcels()[&0].state().state().sound_speed, parcels[&0].state().state().sound_speed);
}

#[test]
fn aggregate_withdrawal_is_bounded_before_incoming_mass_can_hide_it() {
    let parcels = BTreeMap::from([(0, gas(0.75, 400.0, 100000.0, 0.001)),
        (1, gas(0.75, 400.0, 100000.0, 0.001)), (2, gas(0.75, 400.0, 100000.0, 0.001))]);
    let dm = 0.2*parcels[&0].mass_kg();
    let rows = [row(0, 1, dm), row(0, 2, dm), row(1, 0, dm)];
    assert!(matches!(transfer_candidate(&parcels, &rows, &solves(&parcels), controls(), || false),
        Err(Error::Withdrawal { id: 0, .. })));
}

#[test]
fn idle_paths_remain_bitwise_and_malformed_or_conflicting_sources_refuse() {
    let parcels = BTreeMap::from([(0, gas(0.75, 500.0, 200000.0, 0.001)),
        (1, gas(0.25, 400.0, 100000.0, 0.002))]);
    let noop = transfer_candidate(&parcels, &[row(0, 1, 0.0)], &solves(&parcels), controls(), || false).unwrap();
    assert_eq!(noop.parcels(), &parcels);
    assert!(noop.step().parcels.iter().all(|r| r.iterations == 0));
    let dm = 0.01*parcels[&0].mass_kg();
    for rows in [vec![row(0, 1, dm)], vec![row(0, 9, dm)], vec![row(0, 0, dm)],
        vec![row(0, 1, -1.0)], vec![row(0, 1, f64::NAN)], vec![row(0, 1, f64::INFINITY)],
        vec![row(0, 1, 0.0); MAX_SNAPSHOT_TRANSFERS+1]] {
        assert!(transfer_candidate(&parcels, &rows, &solves(&parcels), controls(), || false).is_err());
    }
    let mut missing = solves(&parcels); missing.remove(&1);
    assert!(transfer_candidate(&parcels, &[], &missing, controls(), || false).is_err());
}

#[test]
fn endpoint_refusal_and_cancellation_leave_the_complete_snapshot_unchanged() {
    let parcels = BTreeMap::from([(0, gas(0.75, 500.0, 200000.0, 0.001)),
        (1, gas(0.25, 400.0, 100000.0, 0.001))]);
    let original = parcels.clone();
    let rows = [row(0, 1, 0.01*parcels[&0].mass_kg())];
    let mut narrow = solves(&parcels); narrow.get_mut(&1).unwrap().maximum_temperature_k = 250.0;
    assert!(matches!(transfer_candidate(&parcels, &rows, &narrow, controls(), || false),
        Err(Error::Parcel { id: 1, .. })));
    let mut total = 0;
    transfer_candidate(&parcels, &rows, &solves(&parcels), controls(), || { total += 1; false }).unwrap();
    for stop in [1, total/2, total] {
        let mut calls = 0;
        assert_eq!(transfer_candidate(&parcels, &rows, &solves(&parcels), controls(),
            || { calls += 1; calls == stop }), Err(Error::Cancelled));
        assert_eq!(parcels, original);
    }
}
