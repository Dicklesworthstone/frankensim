//! Independent open storage and bath entropy balances on original public states.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::DiffusionNetworkError;
use transport::{GasReservoirExchange, ReservoirSpeciesTransfer, audit_reservoir_window};

fn species_changes(old: &ClosedGasParcel, new: &ClosedGasParcel) -> Vec<ReservoirSpeciesTransfer> {
    ["A", "B"]
        .into_iter()
        .map(|name| ReservoirSpeciesTransfer {
            species: id(name),
            mass_into_gas_kg: mass(new, name) - mass(old, name),
        })
        .collect()
}

#[test]
fn reservoir_window_matches_independent_mixing_entropy_with_an_initially_absent_species() {
    let run = |formation| {
        let (bank, _) = fixture(MASS, formation, 0.0);
        let old = BTreeMap::from([(1, parcel(&bank, 1.0, 400.0, 1e5))]);
        let new = BTreeMap::from([(1, parcel(&bank, 0.5, 400.0, 1e5))]);
        let bath = bank.evaluate(&point(400.0, 1e5)).unwrap();
        let rows = species_changes(&old[&1], &new[&1]);
        let exchanges = [GasReservoirExchange {
            parcel: 1,
            catalog: &bank,
            state: &bath,
            species: &rows,
            reduced_heat_into_gas_j: 0.0,
        }];
        let report =
            audit_reservoir_window(&old, &new, &exchanges, controls().physical, || false).unwrap();
        // Both synthetic species have the same Cv and entropy datum. The bath
        // exchanges equal amounts at equal partial entropies, so its entropy
        // flux cancels and the gas gains N R ln(2), independently of formation h.
        let expected = old[&1].mass_kg() / MASS * R * 2.0_f64.ln();
        close(report.entropy_production_j_k, expected, 2e-13);
        close(report.energy_residual_j, 0.0, 2e-10);
        close(report.total_mass_residual_kg, 0.0, 1e-18);
        let reversed: Vec<_> = rows
            .iter()
            .map(|r| ReservoirSpeciesTransfer {
                species: r.species.clone(),
                mass_into_gas_kg: -r.mass_into_gas_kg,
            })
            .collect();
        let reverse = [GasReservoirExchange {
            species: &reversed,
            ..exchanges[0]
        }];
        assert!(matches!(
            audit_reservoir_window(&new, &old, &reverse, controls().physical, || false),
            Err(GasTransportError::Diffusion(
                DiffusionNetworkError::EntropyDecrease { .. }
            ))
        ));
        report
    };
    let original = run(0.0);
    let shifted = run(1e24);
    close(shifted.energy_residual_j, original.energy_residual_j, 1e-12);
    close(
        shifted.entropy_production_j_k,
        original.entropy_production_j_k,
        1e-14,
    );
}

#[test]
fn reservoir_window_allows_gas_entropy_loss_and_detects_hidden_energy_defects() {
    let (bank, _) = fixture(MASS, 1e24, 0.0);
    let old = BTreeMap::from([(1, parcel(&bank, 0.5, 420.0, 1e5))]);
    let new = BTreeMap::from([(1, parcel(&bank, 0.5, 400.0, 1e5 * 400.0 / 420.0))]);
    let bath = bank.evaluate(&point(300.0, 1e5)).unwrap();
    let capacity = old[&1].mass_kg() / MASS * 3.0 * R;
    let heat = capacity * (400.0 - 420.0);
    let exchange = GasReservoirExchange {
        parcel: 1,
        catalog: &bank,
        state: &bath,
        species: &[],
        reduced_heat_into_gas_j: heat,
    };
    let report =
        audit_reservoir_window(&old, &new, &[exchange], controls().physical, || false).unwrap();
    let gas_entropy = capacity * (400.0_f64 / 420.0).ln();
    assert!(gas_entropy < 0.0);
    close(
        report.entropy_production_j_k,
        gas_entropy - heat / 300.0,
        1e-12,
    );
    assert!(report.entropy_production_j_k > 0.0);
    close(report.energy_residual_j, 0.0, 1e-10);
    let wrong_heat = GasReservoirExchange {
        reduced_heat_into_gas_j: heat + 1e-4,
        ..exchange
    };
    assert!(matches!(
        audit_reservoir_window(&old, &new, &[wrong_heat], controls().physical, || false),
        Err(GasTransportError::Invalid(_))
    ));
}

#[test]
fn reservoir_window_refuses_changed_sources_even_when_species_move_between_parcels() {
    let (bank, _) = fixture(MASS, 0.0, 0.0);
    let old = BTreeMap::from([
        (1, parcel(&bank, 1.0, 400.0, 1e5)),
        (2, parcel(&bank, 0.0, 400.0, 1e5)),
    ]);
    let valid = BTreeMap::from([(1, old[&2].clone()), (2, old[&1].clone())]);
    audit_reservoir_window(&old, &valid, &[], controls().physical, || false).unwrap();
    let mut data = bank.transport_materials().to_vec();
    data[0].material = card(
        "substituted transport source with identical values",
        [
            scalar("dynamic_viscosity", DynViscosity::DIMS, 1e-5),
            scalar(
                "thermal_conductivity",
                ThermalConductivityQuantityV1::DIMS,
                CONDUCTIVITY,
            ),
        ],
    );
    let substitute =
        ThermochemicalGasModel::new(bank.thermodynamics().clone(), data, None).unwrap();
    let forged = BTreeMap::from([
        (1, old[&2].clone()),
        (2, parcel(&substitute, 1.0, 400.0, 1e5)),
    ]);
    assert!(matches!(
        audit_reservoir_window(&old, &forged, &[], controls().physical, || false),
        Err(GasTransportError::Invalid(_))
    ));
    let mut loose_mass_budget = controls().physical;
    loose_mass_budget.exchange.mass_tolerance_kg =
        2.0 * old.values().map(|p| p.mass_kg()).sum::<f64>();
    assert!(matches!(
        audit_reservoir_window(&old, &forged, &[], loose_mass_budget, || false),
        Err(GasTransportError::Invalid(_))
    ));
    let actual_bath = bank.evaluate(&point(400.0, 1e5)).unwrap();
    let mismatched = GasReservoirExchange {
        parcel: 1,
        catalog: &substitute,
        state: &actual_bath,
        species: &[],
        reduced_heat_into_gas_j: 0.0,
    };
    assert!(matches!(
        audit_reservoir_window(&old, &old, &[mismatched], controls().physical, || false),
        Err(GasTransportError::Invalid(_))
    ));
}

#[test]
fn reservoir_window_keeps_fixed_species_budgets_zero_activity_and_late_cancellation() {
    let (bank, _) = fixture(MASS, 0.0, 0.0);
    let old = BTreeMap::from([(1, parcel(&bank, 1.0, 400.0, 1e5))]);
    let new = BTreeMap::from([(1, parcel(&bank, 0.5, 400.0, 1e5))]);
    let bath = bank.evaluate(&point(400.0, 1e5)).unwrap();
    let rows = species_changes(&old[&1], &new[&1]);
    let exchange = GasReservoirExchange {
        parcel: 1,
        catalog: &bank,
        state: &bath,
        species: &rows,
        reduced_heat_into_gas_j: 0.0,
    };
    let mut wrong = rows.clone();
    wrong[0].mass_into_gas_kg += 1e-7;
    wrong[1].mass_into_gas_kg -= 1e-7;
    assert!(matches!(
        audit_reservoir_window(
            &old,
            &new,
            &[GasReservoirExchange {
                species: &wrong,
                ..exchange
            }],
            controls().physical,
            || false
        ),
        Err(GasTransportError::Invalid(_))
    ));
    let pure = old[&1].model();
    let pure_state = pure.evaluate(&point(400.0, 1e5)).unwrap();
    assert!(matches!(
        audit_reservoir_window(
            &old,
            &new,
            &[GasReservoirExchange {
                catalog: pure,
                state: &pure_state,
                ..exchange
            }],
            controls().physical,
            || false
        ),
        Err(GasTransportError::Invalid(_))
    ));
    let calls = std::cell::Cell::new(0_usize);
    audit_reservoir_window(&old, &new, &[exchange], controls().physical, || {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let cutoff = calls.get() - 1;
    calls.set(0);
    assert!(matches!(
        audit_reservoir_window(&old, &new, &[exchange], controls().physical, || {
            calls.set(calls.get() + 1);
            calls.get() >= cutoff
        }),
        Err(GasTransportError::Cancelled)
    ));
    assert_eq!(old[&1].state().state().temperature, 400.0);
    assert_eq!(new[&1].state().state().temperature, 400.0);
}
