//! Pressure-driven graph composition against original single-path physics.
use super::super::tests::{bound, close, controls as flow_controls, energy, nozzle, parcel};
use super::*;
use crate::fluid_transport::FluidTransport;
use crate::thermochemical_gas::nonideal::energy::tests::bound_with_reference;
use fs_thermochem::SpeciesId;

fn controls() -> SpeciesNetworkSolve {
    SpeciesNetworkSolve {
        maximum_species_fraction: 0.25,
        mass_tolerance_kg: 1e-14,
        energy_tolerance_j: 1e-7,
    }
}
fn nodes(reference: f64) -> BTreeMap<u64, IsochoricRealGasParcel> {
    let model = bound_with_reference(reference);
    BTreeMap::from([
        (10, parcel(model.clone(), 5e6, 0.01)),
        (20, parcel(model.clone(), 4.8e6, 0.01)),
        (30, parcel(model, 3.8e6, 0.01)),
    ])
}
fn endpoints(parcels: &BTreeMap<u64, IsochoricRealGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    parcels.keys().map(|&id| (id, energy())).collect()
}
fn link(id: u64, a: u64, b: u64) -> RealGasNozzleNetworkLink {
    RealGasNozzleNetworkLink {
        id,
        parcels: [a, b],
        nozzle: nozzle(),
        solve: flow_controls(),
    }
}
fn links() -> Vec<RealGasNozzleNetworkLink> {
    vec![link(8, 30, 10), link(4, 20, 30), link(1, 20, 10)]
}

#[test]
fn one_pressure_path_recovers_original_real_nozzle_and_bulk_enthalpy_endpoints() {
    let mut parcels = nodes(0.0);
    parcels.remove(&30);
    let duration = 2e-5;
    let candidate = advance(
        &parcels,
        &[link(7, 20, 10)],
        duration,
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    let (expected, flow) = parcels[&10]
        .through_nozzle_to(&parcels[&20], nozzle(), duration, flow_controls(), || false)
        .unwrap();
    let report = &candidate.step().links[0];
    assert_eq!((report.id, report.donor, report.receiver), (7, 10, 20));
    assert_eq!(report.flow, flow);
    close(
        report.flow.mass_flow_kg_s,
        0.003_429_609_230_876_759_5,
        2e-11,
    );
    close(
        report.flow.throat().temperature_k(),
        296.327_413_254_839_35,
        2e-8,
    );
    assert_eq!(report.transferred_mass_kg, flow.mass_flow_kg_s * duration);
    for (actual, expected) in [&candidate.parcels()[&10], &candidate.parcels()[&20]]
        .into_iter()
        .zip([expected.donor(), expected.receiver()])
    {
        close(actual.mass_kg(), expected.mass_kg(), 1e-14);
        close(
            actual.state().temperature_k(),
            expected.state().temperature_k(),
            2e-8,
        );
        close(
            actual.state().pressure_pa(),
            expected.state().pressure_pa(),
            5e-4,
        );
    }
    close(
        candidate
            .step()
            .exchange
            .fluxes
            .iter()
            .map(|row| row.carried_enthalpy_j)
            .sum(),
        expected.step().transferred_enthalpy_j,
        1e-9,
    );
    assert!(
        candidate.step().exchange.total_energy_residual_j.abs() <= controls().energy_tolerance_j
    );
}

#[test]
fn pressure_branches_selective_counterflow_and_heat_share_one_original_snapshot() {
    let parcels = nodes(0.0);
    let saved = parcels.clone();
    let paths = links();
    let extra = [SpeciesNetworkTransfer {
        donor: 30,
        receiver: 10,
        species: SpeciesId::new("B").unwrap(),
        mass_kg: 2e-8,
    }];
    let heat = [SpeciesNetworkHeatTransfer {
        parcels: [10, 30],
        heat_left_to_right_j: 0.05,
    }];
    let external = [SpeciesNetworkExternalHeat {
        id: 77,
        parcel: 20,
        heat_into_gas_j: 0.01,
    }];
    let candidate = advance_with_exchange(
        &parcels,
        &paths,
        2e-5,
        &extra,
        &heat,
        &external,
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(
        candidate
            .step()
            .links
            .iter()
            .map(|row| row.id)
            .collect::<Vec<_>>(),
        vec![1, 4, 8]
    );
    let mut expected_mass = BTreeMap::from([(10, 0.0), (20, 0.0), (30, 0.0)]);
    let mut expected_energy = BTreeMap::from([(10, -0.05), (20, 0.01), (30, 0.05)]);
    for row in &candidate.step().links {
        let direct = parcels[&row.donor]
            .nozzle_flow_to(
                parcels[&row.receiver].state().pressure_pa(),
                nozzle(),
                flow_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(direct, row.flow);
        let h = row.transferred_mass_kg * parcels[&row.donor].state().specific_enthalpy_j_kg();
        *expected_mass.get_mut(&row.donor).unwrap() -= row.transferred_mass_kg;
        *expected_mass.get_mut(&row.receiver).unwrap() += row.transferred_mass_kg;
        *expected_energy.get_mut(&row.donor).unwrap() -= h;
        *expected_energy.get_mut(&row.receiver).unwrap() += h;
    }
    let partial = parcels[&30]
        .model()
        .partial_calorics_at_density(
            parcels[&30].density_kg_m3(),
            parcels[&30].state().temperature_k(),
        )
        .unwrap();
    let h_b = partial.partial_specific_enthalpies_j_kg()[1];
    *expected_mass.get_mut(&30).unwrap() -= extra[0].mass_kg;
    *expected_mass.get_mut(&10).unwrap() += extra[0].mass_kg;
    *expected_energy.get_mut(&30).unwrap() -= extra[0].mass_kg * h_b;
    *expected_energy.get_mut(&10).unwrap() += extra[0].mass_kg * h_b;
    let mut total_u = 0.0;
    for (&id, before) in &parcels {
        let after = &candidate.parcels()[&id];
        close(
            after.mass_kg() - before.mass_kg(),
            expected_mass[&id],
            1e-14,
        );
        let du = after.mass_kg() * after.state().specific_internal_energy_j_kg()
            - before.mass_kg() * before.state().specific_internal_energy_j_kg();
        close(du, expected_energy[&id], 1e-7);
        total_u += du;
        assert_eq!(after.volume_m3(), before.volume_m3());
        assert_eq!(
            after.state().transport_state().unwrap().density_kg_m3(),
            after.density_kg_m3()
        );
    }
    close(total_u, 0.01, 1e-7);
    let mut reversed = paths;
    reversed.reverse();
    for path in &mut reversed {
        path.parcels.swap(0, 1);
    }
    let repeated = advance_with_exchange(
        &parcels,
        &reversed,
        2e-5,
        &extra,
        &heat,
        &external,
        &endpoints(&parcels),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(candidate, repeated);
    assert_eq!(parcels, saved);
}

#[test]
fn all_outgoing_paths_share_the_strictest_original_donor_inventory_allowance() {
    let parcels = BTreeMap::from([
        (10, parcel(bound(), 5e6, 0.01)),
        (20, parcel(bound(), 1e6, 1.0)),
        (30, parcel(bound(), 1e6, 1.0)),
    ]);
    let mass = parcels[&10].mass_kg();
    let rate = parcels[&10]
        .nozzle_flow_to(1e6, nozzle(), flow_controls(), || false)
        .unwrap()
        .mass_flow_kg_s;
    let duration = 0.06 * mass / rate;
    // Each path alone would withdraw 6%, but the shared 10% original-donor limit
    // must reject their combined 12% before mixing or endpoint publication.
    assert!(
        matches!(advance(&parcels, &[link(1, 10, 20), link(2, 10, 30)], duration,
        &endpoints(&parcels), controls(), || false),
        Err(Error::Withdrawal { donor: 10, fraction, maximum }) if fraction > 0.11 && maximum == 0.1)
    );
    let selective = [SpeciesNetworkTransfer {
        donor: 10,
        receiver: 30,
        species: SpeciesId::new("A").unwrap(),
        mass_kg: 0.05 * mass,
    }];
    assert!(
        matches!(advance_with_exchange(&parcels, &[link(1, 10, 20)], duration,
        &selective, &[], &[], &endpoints(&parcels), controls(), || false),
        Err(Error::Withdrawal { donor: 10, fraction, .. }) if fraction > 0.1)
    );
}

#[test]
fn pressure_crossing_and_late_cancellation_return_no_physical_endpoint() {
    let parcels = BTreeMap::from([
        (10, parcel(bound(), 5e6, 0.001)),
        (20, parcel(bound(), 4.8e6, 0.001)),
    ]);
    let saved = parcels.clone();
    let paths = [link(7, 10, 20)];
    assert!(matches!(
        advance(
            &parcels,
            &paths,
            0.4,
            &endpoints(&parcels),
            controls(),
            || false
        ),
        Err(Error::EquilibriumCrossing {
            link: 7,
            donor: 10,
            receiver: 20
        })
    ));
    let mut polls = 0;
    let good = advance(
        &parcels,
        &paths,
        2e-5,
        &endpoints(&parcels),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        advance(
            &parcels,
            &paths,
            2e-5,
            &endpoints(&parcels),
            controls(),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(Error::Cancelled)
    ));
    let mut unsupported = paths;
    unsupported[0].solve.entropy.minimum_temperature_k = 299.0;
    assert!(matches!(
        advance(
            &parcels,
            &unsupported,
            2e-5,
            &endpoints(&parcels),
            controls(),
            || false
        ),
        Err(Error::Link { id: 7, .. })
    ));
    assert_eq!(parcels, saved);
    assert_eq!(
        advance(
            &parcels,
            &paths,
            2e-5,
            &endpoints(&parcels),
            controls(),
            || false
        )
        .unwrap(),
        good
    );
}

#[test]
fn huge_enthalpy_reference_retains_rates_endpoints_and_closed_path_identity() {
    let base = nodes(0.0);
    let paths = links();
    let expected = advance(&base, &paths, 2e-5, &endpoints(&base), controls(), || false).unwrap();
    let shifted = nodes(1e24);
    let actual = advance(
        &shifted,
        &paths,
        2e-5,
        &endpoints(&shifted),
        controls(),
        || false,
    )
    .unwrap();
    for (a, b) in actual.step().links.iter().zip(&expected.step().links) {
        close(a.flow.mass_flow_kg_s, b.flow.mass_flow_kg_s, 1e-13);
        close(a.flow.velocity_m_s, b.flow.velocity_m_s, 1e-9);
        close(a.transferred_mass_kg, b.transferred_mass_kg, 1e-18);
    }
    for (id, before) in expected.parcels() {
        let after = &actual.parcels()[id];
        close(after.mass_kg(), before.mass_kg(), 1e-14);
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
    let mut closed = paths;
    for path in &mut closed {
        path.nozzle.area_m2 = 0.0;
    }
    let zero = advance(&base, &closed, 1.0, &endpoints(&base), controls(), || false).unwrap();
    assert_eq!(zero.parcels(), &base);
    assert!(
        zero.step()
            .links
            .iter()
            .all(|row| row.transferred_mass_kg == 0.0)
    );
    assert!(zero.step().exchange.fluxes.is_empty());
    assert_eq!(zero.into_parts().0, base);
}

#[test]
fn malformed_paths_and_lowered_species_overflow_refuse_before_publication() {
    let parcels = nodes(0.0);
    for paths in [
        vec![link(1, 10, 10)],
        vec![link(1, 10, 88)],
        vec![link(1, 10, 20), link(1, 20, 30)],
    ] {
        assert!(
            advance(
                &parcels,
                &paths,
                2e-5,
                &endpoints(&parcels),
                controls(),
                || false
            )
            .is_err()
        );
    }
    for duration in [0.0, -1.0, f64::NAN] {
        assert!(
            advance(
                &parcels,
                &links(),
                duration,
                &endpoints(&parcels),
                controls(),
                || false
            )
            .is_err()
        );
    }
    let zeros = vec![
        SpeciesNetworkTransfer {
            donor: 10,
            receiver: 20,
            species: SpeciesId::new("A").unwrap(),
            mass_kg: 0.0
        };
        4095
    ];
    assert!(matches!(
        advance_with_exchange(
            &parcels,
            &[link(1, 10, 20)],
            2e-5,
            &zeros,
            &[],
            &[],
            &endpoints(&parcels),
            controls(),
            || false
        ),
        Err(Error::Invalid(
            "real-gas nozzle species expansion exceeds the shared snapshot bound"
        ))
    ));
}
