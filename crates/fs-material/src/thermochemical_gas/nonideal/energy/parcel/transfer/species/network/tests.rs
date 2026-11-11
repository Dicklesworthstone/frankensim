//! Original snapshot balances on branches/cycles and source-preserving refusal.
use super::super::tests::{
    at, close, energy, pair, partial_enthalpy, real_model, source, species_mass,
};
use super::*;
use crate::fluid_transport::FluidTransport;
use fs_qty::Amount;
use fs_thermochem::SpeciesId;

fn controls() -> SpeciesNetworkSolve {
    SpeciesNetworkSolve {
        maximum_species_fraction: 0.25,
        mass_tolerance_kg: 1e-14,
        energy_tolerance_j: 1e-7,
    }
}
fn nodes(references: [f64; 2]) -> BTreeMap<u64, IsochoricRealGasParcel> {
    let [left, right] = pair(references);
    BTreeMap::from([
        (10, left),
        (20, right),
        (
            30,
            at(real_model(source(0.55, references), 0.03), 325.0, 4e6, 33.0),
        ),
        // Same IDs with incompatible reference data are harmless in an idle,
        // disconnected owner. A global fictitious catalog would reject it.
        (
            99,
            at(
                real_model(source(0.65, [1e6, -2e6]), 0.04),
                310.0,
                2e6,
                99.0,
            ),
        ),
    ])
}
fn endpoints(parcels: &BTreeMap<u64, IsochoricRealGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    parcels.keys().map(|&id| (id, energy())).collect()
}
fn stream(donor: u64, receiver: u64, name: &str, mass: f64) -> SpeciesNetworkTransfer {
    SpeciesNetworkTransfer {
        donor,
        receiver,
        species: SpeciesId::new(name).unwrap(),
        mass_kg: mass,
    }
}
fn streams() -> Vec<SpeciesNetworkTransfer> {
    vec![
        stream(10, 20, "A", 2e-4),
        stream(10, 30, "B", 5e-5),
        stream(20, 30, "A", 1e-4),
        stream(30, 10, "B", 5e-5),
    ]
}

#[test]
fn branches_and_cycles_close_each_original_species_and_absolute_energy_once() {
    let parcels = nodes([0.0; 2]);
    let saved = parcels.clone();
    let rows = streams();
    let heat = [SpeciesNetworkHeatTransfer {
        parcels: [10, 30],
        heat_left_to_right_j: 1.25,
    }];
    let external = [
        SpeciesNetworkExternalHeat {
            id: 7,
            parcel: 20,
            heat_into_gas_j: 3.0,
        },
        SpeciesNetworkExternalHeat {
            id: 8,
            parcel: 30,
            heat_into_gas_j: -0.75,
        },
    ];
    let candidate = exchange_candidate_with_external_heat(
        &parcels,
        &rows,
        &heat,
        &external,
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    let mut total_energy = 0.0;
    for (&id, before) in &parcels {
        let after = &candidate.parcels()[&id];
        let mut carried = 0.0;
        for row in &rows {
            let h = partial_enthalpy(&parcels[&row.donor], &row.species);
            if id == row.receiver {
                carried += row.mass_kg * h;
            }
            if id == row.donor {
                carried -= row.mass_kg * h;
            }
        }
        let extra = match id {
            10 => -1.25,
            20 => 3.0,
            30 => 0.5,
            _ => 0.0,
        };
        let delta_u = after.mass_kg() * after.state().specific_internal_energy_j_kg()
            - before.mass_kg() * before.state().specific_internal_energy_j_kg();
        close(delta_u, carried + extra, controls().energy_tolerance_j);
        total_energy += delta_u;
        for name in ["A", "B"] {
            let species = SpeciesId::new(name).unwrap();
            let intended: f64 = rows
                .iter()
                .filter(|row| row.species == species)
                .map(|row| {
                    if row.receiver == id {
                        row.mass_kg
                    } else if row.donor == id {
                        -row.mass_kg
                    } else {
                        0.0
                    }
                })
                .sum();
            close(
                species_mass(after, &species) - species_mass(before, &species),
                intended,
                controls().mass_tolerance_kg,
            );
        }
        assert_eq!(after.volume_m3(), before.volume_m3());
        assert_eq!(after.state().density_kg_m3(), after.density_kg_m3());
        assert_eq!(
            after.state().transport_state().unwrap().density_kg_m3(),
            after.density_kg_m3()
        );
        if id != 99 {
            let report = candidate
                .step()
                .parcels
                .iter()
                .find(|report| report.id == id)
                .unwrap();
            close(report.net_enthalpy_j, carried, 1e-8);
            close(report.net_additional_heat_j, extra, 1e-12);
            close(
                report.energy_change_j,
                delta_u,
                controls().energy_tolerance_j,
            );
            assert!(report.energy_residual_j.abs() <= controls().energy_tolerance_j);
        }
    }
    close(total_energy, 2.25, controls().energy_tolerance_j);
    assert_eq!(candidate.parcels()[&99], parcels[&99]);
    assert_eq!(parcels, saved);
    assert!(candidate.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(candidate.step().total_mass_residual_kg.abs() <= controls().mass_tolerance_kg);
    assert!(candidate.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
}

#[test]
fn canonical_parallel_rows_and_heat_orientation_preserve_the_entire_candidate() {
    let parcels = nodes([0.0; 2]);
    let original = streams();
    let a = exchange_candidate_with_heat(
        &parcels,
        &original,
        &[SpeciesNetworkHeatTransfer {
            parcels: [10, 30],
            heat_left_to_right_j: 2.0,
        }],
        &BTreeMap::new(),
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    let mut permuted = original.clone();
    permuted[0].mass_kg *= 0.5;
    permuted.push(permuted[0].clone());
    permuted.reverse();
    let b = exchange_candidate_with_heat(
        &parcels,
        &permuted,
        &[
            SpeciesNetworkHeatTransfer {
                parcels: [30, 10],
                heat_left_to_right_j: -0.75,
            },
            SpeciesNetworkHeatTransfer {
                parcels: [10, 30],
                heat_left_to_right_j: 1.25,
            },
        ],
        &BTreeMap::new(),
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(a, b);
    let zero =
        exchange_candidate(&parcels, &[], &endpoints(&parcels), controls(), || false).unwrap();
    assert_eq!(zero.parcels(), &parcels);
    assert!(zero.step().fluxes.is_empty());
    assert_eq!(zero.step().total_energy_residual_j, 0.0);
    assert_eq!(zero.into_parts().0, parcels);
}

#[test]
fn reference_shifts_and_balanced_throughflow_preserve_real_gas_physics() {
    let base = nodes([0.0; 2]);
    let expected =
        exchange_candidate(&base, &streams(), &endpoints(&base), controls(), || false).unwrap();
    for references in [[1e24, -2e24], [-3e24, 4e24]] {
        let parcels = nodes(references);
        let actual = exchange_candidate(
            &parcels,
            &streams(),
            &endpoints(&parcels),
            controls(),
            || false,
        )
        .unwrap();
        for (id, before) in expected.parcels() {
            let after = &actual.parcels()[id];
            assert_eq!(after.mass_kg(), before.mass_kg());
            assert_eq!(after.density_kg_m3(), before.density_kg_m3());
            close(
                after.state().temperature_k(),
                before.state().temperature_k(),
                2e-8,
            );
            close(
                after.state().pressure_pa(),
                before.state().pressure_pa(),
                5e-4,
            );
        }
        assert!(actual.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
    let through = [
        stream(10, 20, "A", 1e-4),
        stream(20, 30, "A", 1e-4),
        stream(30, 10, "A", 1e-4),
    ];
    let candidate =
        exchange_candidate(&base, &through, &endpoints(&base), controls(), || false).unwrap();
    for (&id, before) in &base {
        let after = &candidate.parcels()[&id];
        assert_eq!(after.mass_kg().to_bits(), before.mass_kg().to_bits());
        assert_eq!(
            after.density_kg_m3().to_bits(),
            before.density_kg_m3().to_bits()
        );
        if id != 99 {
            assert!((after.state().temperature_k() - before.state().temperature_k()).abs() > 1e-3);
        }
    }
    // Heat-only exchange uses the same endpoint owner; it retains the original
    // material and its density rather than joining unrelated catalogs.
    let heat = exchange_candidate_with_external_heat(
        &base,
        &[],
        &[],
        &[SpeciesNetworkExternalHeat {
            id: 1,
            parcel: 99,
            heat_into_gas_j: 2.0,
        }],
        &endpoints(&base),
        controls(),
        || false,
    )
    .unwrap();
    let mut expected = base[&99].clone();
    expected.add_heat(2.0, energy(), || false).unwrap();
    assert_eq!(&heat.parcels()[&99], &expected);
}

#[test]
fn gross_original_withdrawal_cannot_be_hidden_by_parallel_incoming_streams() {
    let mut parcels = nodes([0.0; 2]);
    let available = species_mass(&parcels[&10], &SpeciesId::new("A").unwrap());
    let rows = [
        stream(10, 20, "A", 0.15 * available),
        stream(10, 30, "A", 0.15 * available),
        stream(20, 10, "A", 0.05 * available),
    ];
    assert!(matches!(
        exchange_candidate(&parcels, &rows, &endpoints(&parcels), controls(), || false),
        Err(SpeciesNetworkError::Withdrawal { id: 10, .. })
    ));
    let pure_a = source(0.75, [0.0; 2])
        .with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
        .unwrap();
    parcels.insert(20, at(real_model(pure_a, 0.03), 350.0, 3.8e6, 22.0));
    let saved = parcels.clone();
    // Even a larger arrival cannot be forwarded through an initially absent B
    // inventory during this original-state interval.
    let forward = [stream(10, 20, "B", 2e-5), stream(20, 30, "B", 1e-5)];
    assert!(
        exchange_candidate(&parcels, &forward, &endpoints(&parcels), controls(), || {
            false
        })
        .is_err()
    );
    assert_eq!(parcels, saved);
}

#[test]
fn source_and_late_endpoint_refusals_or_cancellation_publish_no_network() {
    let mut parcels = nodes([0.0; 2]);
    let saved = parcels.clone();
    let mut polls = 0;
    exchange_candidate(
        &parcels,
        &streams(),
        &endpoints(&parcels),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        exchange_candidate(
            &parcels,
            &streams(),
            &endpoints(&parcels),
            controls(),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(SpeciesNetworkError::Cancelled)
    ));
    let too_hot = [SpeciesNetworkExternalHeat {
        id: 7,
        parcel: 30,
        heat_into_gas_j: 1e6,
    }];
    assert!(matches!(
        exchange_candidate_with_external_heat(
            &parcels,
            &streams(),
            &[],
            &too_hot,
            &endpoints(&parcels),
            controls(),
            || false
        ),
        Err(SpeciesNetworkError::Parcel { id: 30, .. })
    ));
    assert_eq!(parcels, saved);
    parcels.insert(
        30,
        at(real_model(source(0.55, [0.0; 2]), 0.04), 325.0, 4e6, 33.0),
    );
    assert!(
        exchange_candidate(
            &parcels,
            &streams(),
            &endpoints(&parcels),
            controls(),
            || false
        )
        .is_err()
    );
    let duplicated = [
        SpeciesNetworkExternalHeat {
            id: 7,
            parcel: 20,
            heat_into_gas_j: 1.0,
        },
        SpeciesNetworkExternalHeat {
            id: 7,
            parcel: 30,
            heat_into_gas_j: 2.0,
        },
    ];
    assert!(
        exchange_candidate_with_external_heat(
            &saved,
            &[],
            &[],
            &duplicated,
            &endpoints(&saved),
            controls(),
            || false
        )
        .is_err()
    );
    for rows in [
        vec![stream(10, 10, "A", 1e-5)],
        vec![stream(10, 88, "A", 1e-5)],
        vec![stream(10, 20, "absent", 1e-5)],
        vec![stream(10, 20, "A", f64::NAN)],
        vec![stream(10, 20, "A", 1e-30)],
    ] {
        assert!(
            exchange_candidate(&saved, &rows, &endpoints(&saved), controls(), || false).is_err()
        );
    }
}

#[test]
fn bounded_parallel_expansion_is_combined_before_the_original_energy_solve() {
    let parcels = nodes([0.0; 2]);
    let rows = vec![stream(10, 20, "A", 1e-9); 4096];
    let candidate =
        exchange_candidate(&parcels, &rows, &endpoints(&parcels), controls(), || false).unwrap();
    assert_eq!(candidate.step().fluxes.len(), 1);
    close(candidate.step().fluxes[0].transfer.mass_kg, 4096e-9, 1e-19);
    assert!(
        exchange_candidate(
            &parcels,
            &vec![rows[0].clone(); 4097],
            &endpoints(&parcels),
            controls(),
            || false
        )
        .is_err()
    );
    let mut missing = endpoints(&parcels);
    missing.remove(&99);
    assert!(exchange_candidate(&parcels, &[], &missing, controls(), || false).is_err());
}

#[test]
fn a_counterstream_within_its_mass_allowance_keeps_both_parcels_exactly() {
    // Planted case (b2). Two species cross one face in opposite directions and
    // net less than half an ulp of either parcel's mass. Within an allowance
    // equal to that net both parcels keep their mass and density bit for bit
    // while their compositions change. With half of it, or none, the net is an
    // increment the parcel cannot represent and is refused as it always was.
    let parcels = nodes([0.0; 2]);
    let (a, b) = (1e-6, 1e-6 - 1e-21);
    let net = a - b;
    assert!(net > 0.0);
    assert_eq!(parcels[&10].mass_kg() + net, parcels[&10].mass_kg());
    assert_eq!(parcels[&20].mass_kg() - net, parcels[&20].mass_kg());
    let rows = [stream(10, 20, "A", a), stream(20, 10, "B", b)];
    let run = |allowance: &BTreeMap<u64, f64>| {
        exchange_candidate_with_heat(
            &parcels,
            &rows,
            &[],
            allowance,
            &endpoints(&parcels),
            controls(),
            || false,
        )
    };
    let candidate = run(&BTreeMap::from([(10, net), (20, net)])).unwrap();
    for id in [10, 20] {
        let (old, new) = (&parcels[&id], &candidate.parcels()[&id]);
        assert_eq!(new.mass_kg().to_bits(), old.mass_kg().to_bits());
        assert_eq!(new.density_kg_m3().to_bits(), old.density_kg_m3().to_bits());
    }
    let species_a = SpeciesId::new("A").unwrap();
    close(
        species_mass(&parcels[&10], &species_a)
            - species_mass(&candidate.parcels()[&10], &species_a),
        a,
        controls().mass_tolerance_kg,
    );
    let refused =
        SpeciesNetworkError::Invalid("real-gas network mass/density increment is unrepresentable");
    let half = BTreeMap::from([(10, 0.5 * net), (20, 0.5 * net)]);
    assert_eq!(run(&half).unwrap_err(), refused);
    assert_eq!(run(&BTreeMap::new()).unwrap_err(), refused);
    // An allowance must name an original parcel and be finite and nonnegative.
    for invalid in [(42, net), (10, -net), (10, f64::NAN)] {
        assert!(run(&BTreeMap::from([invalid])).is_err());
    }
}
