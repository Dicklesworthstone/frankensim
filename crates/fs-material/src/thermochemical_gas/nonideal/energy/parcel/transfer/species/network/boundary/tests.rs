//! Open PR mass/energy balances with explicit external phase source charts.
use super::super::super::tests::{
    at, close, energy, pair, partial_enthalpy, real_model, source, species_mass,
};
use super::super::{
    GasEnergySolve, SpeciesNetworkExternalHeat, SpeciesNetworkHeatTransfer, SpeciesNetworkTransfer,
    exchange_candidate, exchange_candidate_with_boundaries,
    exchange_candidate_with_boundaries_and_external_heat,
};
use super::*;
use fs_qty::Amount;

fn controls() -> SpeciesNetworkSolve {
    SpeciesNetworkSolve {
        maximum_species_fraction: 0.25,
        mass_tolerance_kg: 1e-14,
        energy_tolerance_j: 1e-7,
    }
}
fn endpoints(nodes: &BTreeMap<u64, IsochoricRealGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    nodes.keys().map(|&id| (id, energy())).collect()
}
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn internal(donor: u64, receiver: u64, name: &str, mass_kg: f64) -> SpeciesNetworkTransfer {
    SpeciesNetworkTransfer {
        donor,
        receiver,
        species: id(name),
        mass_kg,
    }
}
fn stored(gas: &IsochoricRealGasParcel) -> f64 {
    gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
}
fn ideal_h(gas: &IsochoricRealGasParcel, species: &SpeciesId) -> f64 {
    let source = gas
        .model()
        .source()
        .thermodynamics()
        .components()
        .iter()
        .find(|s| s.species() == species)
        .unwrap();
    source
        .evaluate_properties(Temperature::new(gas.state().temperature_k()))
        .unwrap()
        .h()
        .value()
        / source.molar_mass().value()
}

#[test]
fn actual_gas_boundary_streams_recover_internal_partial_enthalpy_transport() {
    let [left, right] = pair([0.0; 2]);
    let parcels = BTreeMap::from([(10, left), (20, right)]);
    let saved = parcels.clone();
    let mass = 2e-4;
    let incoming =
        RealGasBoundaryStream::inflow(&parcels[&20], &parcels[&10], &id("A"), mass).unwrap();
    let outgoing = RealGasBoundaryStream::outflow(&parcels[&20], &id("A"), mass).unwrap();
    assert!(std::ptr::eq(incoming.catalog(), parcels[&20].model()));
    let boundaries = [incoming.boundary(7, 10), outgoing.boundary(3, 20)];
    let open = exchange_candidate_with_boundaries(
        &parcels,
        &[],
        &[],
        &boundaries,
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    let closed = exchange_candidate(
        &parcels,
        &[internal(20, 10, "A", mass)],
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    let carried = mass * partial_enthalpy(&parcels[&20], &id("A"));
    assert_eq!(open.step().boundaries[0].id, 3);
    close(open.step().boundaries[0].energy_into_gas_j, -carried, 1e-10);
    close(open.step().boundaries[1].energy_into_gas_j, carried, 1e-10);
    for (&node, original) in &parcels {
        let actual = &open.parcels()[&node];
        let expected = &closed.parcels()[&node];
        close(
            actual.state().temperature_k(),
            expected.state().temperature_k(),
            1e-10,
        );
        close(
            actual.state().pressure_pa(),
            expected.state().pressure_pa(),
            1e-5,
        );
        close(
            stored(actual) - stored(original),
            if node == 10 { carried } else { -carried },
            1e-7,
        );
        assert_eq!(actual.mass_kg(), expected.mass_kg());
        assert_eq!(actual.volume_m3(), original.volume_m3());
    }
    assert!(open.step().total_mass_residual_kg.abs() < 1e-14);
    assert!(open.step().maximum_species_mass_residual_kg < 1e-14);
    assert!(open.step().total_energy_residual_j.abs() < 1e-7);
    assert_eq!(parcels, saved);
}

#[test]
fn authored_phase_energy_and_joint_heat_enter_the_same_open_endpoint_once() {
    let [left, right] = pair([0.0; 2]);
    let parcels = BTreeMap::from([
        (10, left),
        (20, right),
        (
            99,
            at(real_model(source(0.6, [1e9, -2e9]), 0.04), 320.0, 2e6, 99.0),
        ),
    ]);
    let rows = [
        BoundarySpeciesTransfer {
            species: id("A"),
            incoming_mass_kg: 2e-4,
            outgoing_mass_kg: 0.0,
            net_mass_into_gas_kg: 2e-4,
            relative_energy_into_gas_j: 1.75,
        },
        // A deliberately authored uptake law. A PR gas departure must NOT be
        // added to this explicit external-phase energy a second time.
        BoundarySpeciesTransfer {
            species: id("B"),
            incoming_mass_kg: 0.0,
            outgoing_mass_kg: 1e-5,
            net_mass_into_gas_kg: -1e-5,
            relative_energy_into_gas_j: 0.5,
        },
    ];
    let boundaries = [RealGasSpeciesNetworkBoundary {
        id: 1,
        parcel: 10,
        catalog: parcels[&10].model(),
        species: &rows,
    }];
    let streams = [internal(10, 20, "B", 3e-5)];
    let heat = [SpeciesNetworkHeatTransfer {
        parcels: [10, 20],
        heat_left_to_right_j: 0.75,
    }];
    let external = [SpeciesNetworkExternalHeat {
        id: 1,
        parcel: 20,
        heat_into_gas_j: 2.0,
    }];
    let candidate = exchange_candidate_with_boundaries_and_external_heat(
        &parcels,
        &streams,
        &heat,
        &boundaries,
        &external,
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    let boundary_energy: f64 = rows
        .iter()
        .map(|r| {
            r.net_mass_into_gas_kg * ideal_h(&parcels[&10], &r.species)
                + r.relative_energy_into_gas_j
        })
        .sum();
    let carried = 3e-5 * partial_enthalpy(&parcels[&10], &id("B"));
    close(
        stored(&candidate.parcels()[&10]) - stored(&parcels[&10]),
        boundary_energy - carried - 0.75,
        1e-7,
    );
    close(
        stored(&candidate.parcels()[&20]) - stored(&parcels[&20]),
        carried + 2.75,
        1e-7,
    );
    close(
        candidate
            .parcels()
            .values()
            .map(|p| p.mass_kg())
            .sum::<f64>()
            - parcels.values().map(|p| p.mass_kg()).sum::<f64>(),
        1.9e-4,
        1e-14,
    );
    close(
        candidate
            .step()
            .boundaries
            .iter()
            .map(|row| row.energy_into_gas_j)
            .sum(),
        boundary_energy,
        1e-10,
    );
    assert_eq!(candidate.parcels()[&99], parcels[&99]);
    assert!(candidate.step().total_mass_residual_kg.abs() < 1e-14);
    assert!(candidate.step().maximum_species_mass_residual_kg < 1e-14);
    assert!(candidate.step().total_energy_residual_j.abs() < 1e-7);
}

#[test]
fn independent_large_species_datums_preserve_counterflow_heat_and_density() {
    let run = |references| {
        let [receiver, donor] = pair(references);
        let parcels = BTreeMap::from([(10, receiver)]);
        let inlet = RealGasBoundaryStream::inflow(&donor, &parcels[&10], &id("A"), 2e-4).unwrap();
        let outlet = RealGasBoundaryStream::outflow(&parcels[&10], &id("A"), 2e-4).unwrap();
        let rows = [inlet.boundary(3, 10), outlet.boundary(5, 10)];
        let answer = exchange_candidate_with_boundaries(
            &parcels,
            &[],
            &[],
            &rows,
            &endpoints(&parcels),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(answer.parcels()[&10].mass_kg(), parcels[&10].mass_kg());
        assert_eq!(
            answer.parcels()[&10].density_kg_m3().to_bits(),
            parcels[&10].density_kg_m3().to_bits()
        );
        assert!(
            answer.parcels()[&10].state().temperature_k() > parcels[&10].state().temperature_k()
        );
        (
            answer,
            inlet.transfer().relative_energy_into_gas_j,
            outlet.transfer().relative_energy_into_gas_j,
        )
    };
    let ordinary = run([0.0; 2]);
    let shifted = run([1e16, -2e16]);
    assert_eq!(ordinary.1, shifted.1);
    assert_eq!(ordinary.2, shifted.2);
    close(
        ordinary.0.parcels()[&10].state().temperature_k(),
        shifted.0.parcels()[&10].state().temperature_k(),
        1e-10,
    );
    close(
        ordinary.0.parcels()[&10].state().pressure_pa(),
        shifted.0.parcels()[&10].state().pressure_pa(),
        1e-5,
    );
    close(
        ordinary.0.step().parcels[0].net_enthalpy_j,
        shifted.0.step().parcels[0].net_enthalpy_j,
        1e-10,
    );
}

#[test]
fn incoming_species_require_original_complete_compatible_pr_catalogs() {
    let base = source(0.75, [0.0; 2]);
    let pure_a = base
        .with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
        .unwrap();
    let pure_b = base
        .with_mole_amounts(&[Amount::new(0.0), Amount::new(1.0)])
        .unwrap();
    let receiver = at(real_model(pure_a, 0.03), 300.0, 5e6, 10.0);
    let parcels = BTreeMap::from([(10, receiver)]);
    let donor = at(real_model(base, 0.03), 350.0, 3.8e6, 20.0);
    let incoming = RealGasBoundaryStream::inflow(&donor, &parcels[&10], &id("B"), 2e-5).unwrap();
    let candidate = exchange_candidate_with_boundaries(
        &parcels,
        &[],
        &[],
        &[incoming.boundary(1, 10)],
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    close(
        species_mass(&candidate.parcels()[&10], &id("B")),
        2e-5,
        1e-14,
    );
    assert_eq!(
        candidate.parcels()[&10]
            .model()
            .mechanics()
            .binary_interactions(),
        &[0.0, 0.03, 0.03, 0.0]
    );
    let pure_donor = at(real_model(pure_b, 0.03), 350.0, 3.8e6, 20.0);
    let missing =
        RealGasBoundaryStream::inflow(&pure_donor, &parcels[&10], &id("B"), 2e-5).unwrap();
    assert!(
        exchange_candidate_with_boundaries(
            &parcels,
            &[],
            &[],
            &[missing.boundary(1, 10)],
            &endpoints(&parcels),
            controls(),
            || false
        )
        .is_err()
    );
    let [mixed, _] = pair([0.0; 2]);
    let mixed = BTreeMap::from([(10, mixed)]);
    for incompatible in [
        real_model(source(0.35, [10.0, 0.0]), 0.03),
        real_model(source(0.35, [0.0; 2]), 0.04),
    ] {
        let donor = at(incompatible, 350.0, 3.8e6, 20.0);
        let incoming = RealGasBoundaryStream::inflow(&donor, &mixed[&10], &id("A"), 2e-5).unwrap();
        assert!(
            exchange_candidate_with_boundaries(
                &mixed,
                &[],
                &[],
                &[incoming.boundary(1, 10)],
                &endpoints(&mixed),
                controls(),
                || false
            )
            .is_err()
        );
    }
}

#[test]
fn boundary_and_internal_gross_withdrawal_share_original_inventory() {
    let [left, right] = pair([0.0; 2]);
    let parcels = BTreeMap::from([(10, left), (20, right)]);
    let mass = 0.15 * species_mass(&parcels[&10], &id("A"));
    let outgoing = RealGasBoundaryStream::outflow(&parcels[&10], &id("A"), mass).unwrap();
    let incoming =
        RealGasBoundaryStream::inflow(&parcels[&20], &parcels[&10], &id("A"), mass).unwrap();
    assert!(matches!(
        exchange_candidate_with_boundaries(
            &parcels,
            &[internal(10, 20, "A", mass)],
            &[],
            &[incoming.boundary(1, 10), outgoing.boundary(2, 10)],
            &endpoints(&parcels),
            controls(),
            || false
        ),
        Err(Error::Withdrawal { id: 10, .. })
    ));
    let pure = source(0.75, [0.0; 2])
        .with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
        .unwrap();
    let pure = BTreeMap::from([(10, at(real_model(pure, 0.03), 300.0, 5e6, 10.0))]);
    let absent = [BoundarySpeciesTransfer {
        species: id("B"),
        incoming_mass_kg: 2e-5,
        outgoing_mass_kg: 1e-5,
        net_mass_into_gas_kg: 1e-5,
        relative_energy_into_gas_j: 0.0,
    }];
    assert!(matches!(
        exchange_candidate_with_boundaries(
            &pure,
            &[],
            &[],
            &[RealGasSpeciesNetworkBoundary {
                id: 1,
                parcel: 10,
                catalog: parcels[&20].model(),
                species: &absent
            }],
            &endpoints(&pure),
            controls(),
            || false
        ),
        Err(Error::Withdrawal { id: 10, .. })
    ));
}

#[test]
fn boundary_order_zero_rows_and_late_refusal_preserve_transactionality() {
    let [left, right] = pair([0.0; 2]);
    let parcels = BTreeMap::from([(10, left), (20, right)]);
    let saved = parcels.clone();
    let inlet =
        RealGasBoundaryStream::inflow(&parcels[&20], &parcels[&10], &id("A"), 2e-5).unwrap();
    let outlet = RealGasBoundaryStream::outflow(&parcels[&10], &id("B"), 1e-5).unwrap();
    let rows = [inlet.boundary(4, 10), outlet.boundary(2, 10)];
    let mut calls = 0;
    let expected = exchange_candidate_with_boundaries(
        &parcels,
        &[],
        &[],
        &rows,
        &endpoints(&parcels),
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let reversed = [rows[1], rows[0]];
    assert_eq!(
        expected,
        exchange_candidate_with_boundaries(
            &parcels,
            &[],
            &[],
            &reversed,
            &endpoints(&parcels),
            controls(),
            || false
        )
        .unwrap()
    );
    let mut count = 0;
    assert!(matches!(
        exchange_candidate_with_boundaries(
            &parcels,
            &[],
            &[],
            &rows,
            &endpoints(&parcels),
            controls(),
            || {
                count += 1;
                count == calls
            }
        ),
        Err(Error::Cancelled)
    ));
    let huge = [SpeciesNetworkExternalHeat {
        id: 9,
        parcel: 10,
        heat_into_gas_j: 1e8,
    }];
    assert!(
        exchange_candidate_with_boundaries_and_external_heat(
            &parcels,
            &[],
            &[],
            &rows,
            &huge,
            &endpoints(&parcels),
            controls(),
            || false
        )
        .is_err()
    );
    let zero = [BoundarySpeciesTransfer {
        species: id("A"),
        incoming_mass_kg: 0.0,
        outgoing_mass_kg: 0.0,
        net_mass_into_gas_kg: 0.0,
        relative_energy_into_gas_j: 0.0,
    }];
    let incompatible = real_model(source(0.35, [1e9, 2e9]), 0.04);
    let idle = exchange_candidate_with_boundaries(
        &parcels,
        &[],
        &[],
        &[RealGasSpeciesNetworkBoundary {
            id: 0,
            parcel: 10,
            catalog: &incompatible,
            species: &zero,
        }],
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(idle.parcels(), &parcels);
    assert!(idle.step().boundaries.is_empty());
    assert_eq!(parcels, saved);
}

#[test]
fn boundary_schema_and_joint_submitted_work_bound_are_enforced() {
    let [left, right] = pair([0.0; 2]);
    let parcels = BTreeMap::from([(10, left), (20, right)]);
    let outgoing = RealGasBoundaryStream::outflow(&parcels[&10], &id("A"), 1e-7).unwrap();
    let owner = outgoing.boundary(1, 10);
    assert!(
        exchange_candidate_with_boundaries(
            &parcels,
            &[],
            &[],
            &[owner, owner],
            &endpoints(&parcels),
            controls(),
            || false
        )
        .is_err()
    );
    let row = outgoing.transfer().clone();
    for rows in [
        vec![row.clone(), row.clone()],
        vec![BoundarySpeciesTransfer {
            incoming_mass_kg: 0.0,
            outgoing_mass_kg: 0.0,
            ..row.clone()
        }],
        vec![BoundarySpeciesTransfer {
            relative_energy_into_gas_j: f64::NAN,
            ..row.clone()
        }],
        vec![BoundarySpeciesTransfer {
            net_mass_into_gas_kg: -1e-30,
            outgoing_mass_kg: 1e-30,
            ..row.clone()
        }],
    ] {
        assert!(
            exchange_candidate_with_boundaries(
                &parcels,
                &[],
                &[],
                &[RealGasSpeciesNetworkBoundary {
                    species: &rows,
                    ..owner
                }],
                &endpoints(&parcels),
                controls(),
                || false
            )
            .is_err()
        );
    }
    let streams = vec![internal(10, 20, "A", 1e-9); 4094];
    exchange_candidate_with_boundaries(
        &parcels,
        &streams,
        &[],
        &[owner],
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    let heat = [SpeciesNetworkExternalHeat {
        id: 1,
        parcel: 10,
        heat_into_gas_j: 0.0,
    }];
    assert!(
        exchange_candidate_with_boundaries_and_external_heat(
            &parcels,
            &streams,
            &[],
            &[owner],
            &heat,
            &endpoints(&parcels),
            controls(),
            || false
        )
        .is_err()
    );
}
