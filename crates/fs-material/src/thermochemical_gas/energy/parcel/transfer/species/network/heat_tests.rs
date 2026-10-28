//! Original source/EOS/caloric owners with synthetic ideal-gas fixtures.
//! These tests verify integrated network energy accounting, not an effusion law.

use super::*;
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermochem::{
    Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const MOLAR_MASS: f64 = 0.028;
const CV: f64 = 2.5 * R / MOLAR_MASS;

fn parcel(id: u64, temperature: f64, reference: f64) -> ClosedGasParcel {
    let species = species_with_energy_reference("A", MOLAR_MASS, 3.5, 0.0, reference);
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 2e-5, 0.03)],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(
        model,
        point(temperature, 100_000.0)
            .with("local-lot", id as f64)
            .unwrap(),
        0.001,
    )
    .unwrap()
}
fn gases(reference: f64) -> BTreeMap<u64, ClosedGasParcel> {
    [
        (1, parcel(1, 320.0, reference)),
        (2, parcel(2, 400.0, reference)),
        (3, parcel(3, 500.0, reference)),
    ]
    .into()
}
fn endpoints(gas: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gas.keys()
        .map(|&id| {
            (
                id,
                GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
            )
        })
        .collect()
}
fn controls() -> SpeciesNetworkSolve {
    SpeciesNetworkSolve {
        maximum_species_fraction: 0.2,
        mass_tolerance_kg: 1e-13,
        energy_tolerance_j: 1e-7,
    }
}
fn mass(donor: u64, receiver: u64, mass_kg: f64) -> SpeciesNetworkTransfer {
    SpeciesNetworkTransfer {
        donor,
        receiver,
        species: SpeciesId::new("A").unwrap(),
        mass_kg,
    }
}
fn heat(left: u64, right: u64, joules: f64) -> SpeciesNetworkHeatTransfer {
    SpeciesNetworkHeatTransfer {
        parcels: [left, right],
        heat_left_to_right_j: joules,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (allowance {tolerance:.3e})"
    );
}
fn account(candidate: &SpeciesNetworkCandidate, id: u64) -> &SpeciesNetworkParcelStep {
    candidate
        .step()
        .parcels
        .iter()
        .find(|row| row.id == id)
        .unwrap()
}

#[test]
fn heat_only_nodes_use_the_original_constant_cv_owner_and_idle_nodes_are_unchanged() {
    let gas = gases(0.0);
    let saved = gas.clone();
    let result = exchange_candidate_with_heat(
        &gas,
        &[],
        &[heat(1, 2, 0.25)],
        &endpoints(&gas),
        controls(),
        || false,
    )
    .unwrap();
    assert!(result.step().fluxes.is_empty());
    assert_eq!(result.step().heat, vec![heat(1, 2, 0.25)]);
    for (id, q) in [(1, -0.25), (2, 0.25)] {
        let old = &gas[&id];
        let next = &result.parcels()[&id];
        near(
            next.state().state().temperature,
            old.state().state().temperature + q / (old.mass_kg() * CV),
            2e-8,
        );
        assert_eq!(next.mass_kg(), old.mass_kg());
        assert_eq!(next.model(), old.model());
        assert_eq!(next.volume_m3(), old.volume_m3());
        assert_eq!(next.point().axes().get("local-lot"), Some(&(id as f64)));
        let row = account(&result, id);
        assert_eq!(row.net_enthalpy_j, 0.0);
        assert_eq!(row.net_additional_heat_j, q);
        assert_eq!(row.incoming_mass_kg, 0.0);
        assert_eq!(row.outgoing_mass_kg, 0.0);
        near(row.energy_change_j, q, controls().energy_tolerance_j);
    }
    assert_eq!(&result.parcels()[&3], &gas[&3]);
    assert_eq!(account(&result, 3).iterations, 0);
    assert_eq!(gas, saved);
    assert!(result.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert_eq!(result.step().total_mass_residual_kg, 0.0);
}

#[test]
fn species_enthalpy_and_additional_heat_share_one_complete_endpoint_energy_target() {
    let gas = gases(0.0);
    let dm = 2e-6;
    let q = -0.25; // A physical stream correction may oppose the mass direction.
    let transfer = mass(1, 2, dm);
    let result = exchange_candidate_with_heat(
        &gas,
        &[transfer.clone()],
        &[heat(1, 2, q)],
        &endpoints(&gas),
        controls(),
        || false,
    )
    .unwrap();
    let carried = dm * 3.5 * R * gas[&1].state().state().temperature / MOLAR_MASS;
    for (id, sign) in [(1, -1.0), (2, 1.0)] {
        let old = &gas[&id];
        let next = &result.parcels()[&id];
        let target_mass = old.mass_kg() + sign * dm;
        let target_energy =
            old.mass_kg() * CV * old.state().state().temperature + sign * (carried + q);
        near(next.mass_kg(), target_mass, controls().mass_tolerance_kg);
        near(
            next.state().state().temperature,
            target_energy / (target_mass * CV),
            2e-8,
        );
        let row = account(&result, id);
        near(row.net_enthalpy_j, sign * carried, 1e-12);
        assert_eq!(row.net_additional_heat_j, sign * q);
        near(
            row.energy_change_j,
            sign * (carried + q),
            controls().energy_tolerance_j,
        );
        near(
            row.energy_residual_j,
            row.energy_change_j - row.net_enthalpy_j - row.net_additional_heat_j,
            1e-12,
        );
    }
    assert_eq!(result.step().fluxes.len(), 1);
    assert_eq!(result.step().fluxes[0].transfer, transfer);
    assert!(result.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    assert!(result.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
}

#[test]
fn signed_heat_edges_are_canonical_and_the_correction_is_independent_of_energy_reference() {
    let transfers = [mass(1, 2, 1e-6), mass(3, 2, 2e-6)];
    let heat_rows = [
        heat(2, 1, -0.3),
        heat(1, 2, 0.1),
        heat(3, 2, 0.2),
        heat(2, 3, 0.05),
        heat(1, 3, -0.0),
    ];
    let gas = gases(0.0);
    let result = exchange_candidate_with_heat(
        &gas,
        &transfers,
        &heat_rows,
        &endpoints(&gas),
        controls(),
        || false,
    )
    .unwrap();
    let reversed = heat_rows
        .iter()
        .rev()
        .map(|row| heat(row.parcels[1], row.parcels[0], -row.heat_left_to_right_j))
        .collect::<Vec<_>>();
    let reordered = exchange_candidate_with_heat(
        &gas,
        &[transfers[1].clone(), transfers[0].clone()],
        &reversed,
        &endpoints(&gas),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result, reordered);
    assert_eq!(result.step().heat.len(), 2);
    assert_eq!(result.step().heat[0].parcels, [1, 2]);
    assert_eq!(result.step().heat[1].parcels, [2, 3]);
    near(result.step().heat[0].heat_left_to_right_j, 0.4, 1e-15);
    near(result.step().heat[1].heat_left_to_right_j, -0.15, 1e-15);

    let reference = 1e8;
    let shifted = gases(reference);
    let shifted_result = exchange_candidate_with_heat(
        &shifted,
        &transfers,
        &heat_rows,
        &endpoints(&shifted),
        controls(),
        || false,
    )
    .unwrap();
    for id in [1, 2, 3] {
        let base = account(&result, id);
        let actual = account(&shifted_result, id);
        near(
            shifted_result.parcels()[&id].state().state().temperature,
            result.parcels()[&id].state().state().temperature,
            2e-7,
        );
        near(
            actual.net_enthalpy_j - base.net_enthalpy_j,
            reference * (base.incoming_mass_kg - base.outgoing_mass_kg),
            1e-7,
        );
        assert_eq!(actual.net_additional_heat_j, base.net_additional_heat_j);
        assert!(actual.energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
    assert_eq!(
        exchange_candidate(&gas, &transfers, &endpoints(&gas), controls(), || false).unwrap(),
        exchange_candidate_with_heat(&gas, &transfers, &[], &endpoints(&gas), controls(), || {
            false
        })
        .unwrap()
    );
    let cancelled_heat = exchange_candidate_with_heat(
        &gas,
        &[],
        &[heat(1, 2, 0.5), heat(2, 1, 0.5)],
        &endpoints(&gas),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(cancelled_heat.parcels(), &gas);
    assert!(cancelled_heat.step().heat.is_empty());
}

#[test]
fn malformed_overlarge_withdrawal_and_cancelled_work_preserve_the_original_snapshot() {
    let gas = gases(0.0);
    let saved = gas.clone();
    let solve = endpoints(&gas);
    for bad in [
        heat(1, 1, 1.0),
        heat(1, 99, 1.0),
        heat(1, 2, f64::NAN),
        heat(1, 2, f64::INFINITY),
        heat(1, 2, 1e9),
    ] {
        assert!(
            exchange_candidate_with_heat(&gas, &[], &[bad], &solve, controls(), || false).is_err()
        );
    }
    assert!(
        exchange_candidate_with_heat(
            &gas,
            &[],
            &[heat(1, 2, f64::MAX), heat(1, 2, f64::MAX)],
            &solve,
            controls(),
            || false
        )
        .is_err()
    );
    let too_many = vec![heat(1, 2, 0.0); MAX_SNAPSHOT_TRANSFERS];
    assert!(matches!(
        exchange_candidate_with_heat(
            &gas,
            &[mass(1, 2, 1e-6)],
            &too_many,
            &solve,
            controls(),
            || false,
        ),
        Err(SpeciesNetworkError::Invalid(_))
    ));
    let excessive = mass(1, 2, 0.3 * gas[&1].mass_kg());
    assert!(matches!(
        exchange_candidate_with_heat(
            &gas,
            &[excessive],
            &[heat(2, 1, 10.0)],
            &solve,
            controls(),
            || false,
        ),
        Err(SpeciesNetworkError::Withdrawal { id: 1, .. })
    ));

    let transfer = [mass(1, 2, 1e-6)];
    let heat_rows = [heat(1, 2, 0.1)];
    let mut calls = 0;
    exchange_candidate_with_heat(&gas, &transfer, &heat_rows, &solve, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut seen = 0;
        assert!(matches!(
            exchange_candidate_with_heat(&gas, &transfer, &heat_rows, &solve, controls(), || {
                seen += 1;
                seen == stop
            },),
            Err(SpeciesNetworkError::Cancelled)
        ));
    }
    assert_eq!(gas, saved);
}

#[test]
fn large_common_energy_datum_does_not_change_heat_only_temperatures_or_joule_closure() {
    let mut solve = controls();
    solve.energy_tolerance_j = 1e-10;
    let make = |reference| {
        let gas = gases(reference);
        let mut endpoint_controls = endpoints(&gas);
        for control in endpoint_controls.values_mut() {
            control.energy_tolerance_j_kg = 1e-9;
        }
        let result = exchange_candidate_with_heat(
            &gas,
            &[],
            &[heat(1, 2, 0.25)],
            &endpoint_controls,
            solve,
            || false,
        )
        .unwrap();
        (gas, result)
    };
    let (base, expected) = make(0.0);
    let (shifted, actual) = make(1e15);
    for id in [1, 2, 3] {
        near(
            actual.parcels()[&id].state().state().temperature,
            expected.parcels()[&id].state().state().temperature,
            1e-10,
        );
        near(
            account(&actual, id).energy_change_j,
            account(&expected, id).energy_change_j,
            1e-12,
        );
        assert_eq!(actual.parcels()[&id].model(), shifted[&id].model());
        assert_eq!(actual.parcels()[&id].mass_kg(), base[&id].mass_kg());
    }
    assert!(actual.step().total_energy_residual_j.abs() <= solve.energy_tolerance_j);
}

#[test]
fn equal_gross_counterstreams_retain_small_real_heat_at_a_huge_formation_datum() {
    let gas: BTreeMap<u64, ClosedGasParcel> =
        [(1, parcel(1, 300.00001, 1e15)), (2, parcel(2, 300.0, 1e15))].into();
    let dm = 1e-6;
    let delta_t = gas[&1].state().state().temperature - gas[&2].state().state().temperature;
    // Effusive beam energy h-RT/2: equal opposing masses have zero net mass,
    // but carry dm*(Cp-R/2)*deltaT, independently of any common energy datum.
    let original_enthalpy = dm * 3.5 * R / MOLAR_MASS * delta_t;
    let correction = -dm * 0.5 * R / MOLAR_MASS * delta_t;
    let complete_energy = original_enthalpy + correction;
    let mut limits = controls();
    limits.energy_tolerance_j = 1e-11;
    let mut endpoint_controls = endpoints(&gas);
    for control in endpoint_controls.values_mut() {
        control.energy_tolerance_j_kg = 1e-10;
    }
    let result = exchange_candidate_with_heat(
        &gas,
        &[mass(1, 2, dm), mass(2, 1, dm)],
        &[heat(1, 2, correction)],
        &endpoint_controls,
        limits,
        || false,
    )
    .unwrap();
    assert_eq!(result.step().fluxes.len(), 2); // Preserve both gross withdrawals.
    for (id, sign) in [(1, -1.0), (2, 1.0)] {
        let old = &gas[&id];
        let next = &result.parcels()[&id];
        assert_eq!(next.mass_kg(), old.mass_kg());
        near(
            next.state().state().temperature,
            old.state().state().temperature + sign * complete_energy / (old.mass_kg() * CV),
            1e-10,
        );
        near(
            account(&result, id).net_enthalpy_j,
            sign * original_enthalpy,
            1e-17,
        );
        near(
            account(&result, id).energy_change_j,
            sign * complete_energy,
            limits.energy_tolerance_j,
        );
        near(account(&result, id).incoming_mass_kg, dm, 0.0);
        near(account(&result, id).outgoing_mass_kg, dm, 0.0);
    }
    assert!(result.parcels()[&1].state().state().temperature < gas[&1].state().state().temperature);
    assert!(result.parcels()[&2].state().state().temperature > gas[&2].state().state().temperature);
    assert!(result.step().total_energy_residual_j.abs() <= limits.energy_tolerance_j);
}

#[test]
fn changing_mixture_cancels_each_distinct_species_datum_before_the_single_caloric_solve() {
    let make = |references: [f64; 2]| -> BTreeMap<u64, ClosedGasParcel> {
        let model = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![
                    species_with_energy_reference("A", MOLAR_MASS, 3.5, 0.0, references[0]),
                    species_with_energy_reference("B", 0.004, 2.5, 0.0, references[1]),
                ],
                Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
            )
            .unwrap(),
            vec![transport("A", 2e-5, 0.03), transport("B", 2e-5, 0.03)],
            None,
        )
        .unwrap();
        [
            (
                1,
                ClosedGasParcel::new(model.clone(), point(500.0, 100_000.0), 0.001).unwrap(),
            ),
            (
                2,
                ClosedGasParcel::new(model, point(300.0, 100_000.0), 0.001).unwrap(),
            ),
        ]
        .into()
    };
    let moles = 1e-5;
    let streams = [
        mass(1, 2, moles * MOLAR_MASS),
        SpeciesNetworkTransfer {
            donor: 2,
            receiver: 1,
            species: SpeciesId::new("B").unwrap(),
            mass_kg: moles * 0.004,
        },
    ];
    let q = 0.005;
    let base = make([0.0, 0.0]);
    let shifted = make([1e15, -2e15]);
    let mut limits = controls();
    limits.energy_tolerance_j = 1e-9;
    let solve = |gas: &BTreeMap<u64, ClosedGasParcel>| {
        let mut endpoint_controls = endpoints(gas);
        for control in endpoint_controls.values_mut() {
            control.energy_tolerance_j_kg = 1e-9;
        }
        exchange_candidate_with_heat(
            gas,
            &streams,
            &[heat(1, 2, q)],
            &endpoint_controls,
            limits,
            || false,
        )
        .unwrap()
    };
    let expected = solve(&base);
    let result = solve(&shifted);
    let id_a = SpeciesId::new("A").unwrap();
    let id_b = SpeciesId::new("B").unwrap();
    let partial_stream_energy = moles * R * (3.5 * 500.0 - 2.5 * 300.0) + q;
    for (id, sign) in [(1, -1.0), (2, 1.0)] {
        let old = &base[&id];
        let a = old.species_mass_kg(&id_a).unwrap() / MOLAR_MASS;
        let b = old.species_mass_kg(&id_b).unwrap() / 0.004;
        let next_a = a + sign * moles;
        let next_b = b - sign * moles;
        let caloric_cv = R * (2.5 * next_a + 1.5 * next_b);
        let initial_sensible = R * old.state().state().temperature * (2.5 * a + 1.5 * b);
        let target_t = (initial_sensible + sign * partial_stream_energy) / caloric_cv;
        near(
            result.parcels()[&id].state().state().temperature,
            target_t,
            1e-8,
        );
        near(
            result.parcels()[&id].state().state().temperature,
            expected.parcels()[&id].state().state().temperature,
            1e-9,
        );
        near(
            result.parcels()[&id].species_mass_kg(&id_a).unwrap(),
            next_a * MOLAR_MASS,
            limits.mass_tolerance_kg,
        );
        near(
            result.parcels()[&id].species_mass_kg(&id_b).unwrap(),
            next_b * 0.004,
            limits.mass_tolerance_kg,
        );
        assert!(account(&result, id).energy_residual_j.abs() <= limits.energy_tolerance_j);
        assert!(account(&result, id).maximum_species_mass_residual_kg <= limits.mass_tolerance_kg);
    }
    assert!(result.step().total_energy_residual_j.abs() <= limits.energy_tolerance_j);
}
