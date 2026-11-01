//! Simultaneous finite continuum diffusion and explicitly held bath exchange.

use super::*;
use fs_gas::material::transport::parcel::reservoir::{
    MolecularReservoir, ReservoirMolecularCandidate, ReservoirMolecularLink,
};

fn bath(bank: &ThermochemicalGasModel, x: f64, t: f64) -> MolecularReservoir {
    MolecularReservoir::new(
        bank.with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        point(t, 1e5),
        || false,
    )
    .unwrap()
}
fn boundary(bath: &MolecularReservoir) -> ReservoirMolecularLink<'_> {
    ReservoirMolecularLink {
        id: 8,
        parcel: 1,
        reservoir: bath,
        area_m2: 0.01,
        distance_m: 0.1,
    }
}
fn network(
    op: &mut Operator,
    gas: &BTreeMap<u64, ClosedGasParcel>,
    internal: &[ParcelMolecularLink],
    baths: &[ReservoirMolecularLink<'_>],
    dt: f64,
) -> Result<ReservoirMolecularCandidate, ParcelMolecularError> {
    let endpoints = gas.keys().map(|&id| (id, endpoint())).collect();
    op.network_candidate(gas, internal, baths, dt, &endpoints, controls(), || false)
}

#[test]
fn one_original_heat_snapshot_matches_two_finite_cells_and_held_bath() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 1.0, 390.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 1.0, 400.0, 1e5)),
        (2, parcel(&bank, 1.0, 420.0, 1.05e5)),
    ]);
    let next = network(&mut op, &original, &[link()], &[boundary(&bath)], 0.1).unwrap();
    let conductance = CONDUCTIVITY * 0.01 / 0.1;
    let internal_q = conductance * (400.0 - 420.0) * 0.1;
    let bath_q = conductance * (390.0 - 400.0) * 0.1;
    let c1 = original[&1].mass_kg() * 3.0 * R / MASS;
    let c2 = original[&2].mass_kg() * 3.0 * R / MASS;
    let t1 = 400.0 + (bath_q - internal_q) / c1;
    let t2 = 420.0 + internal_q / c2;
    close(next.parcels()[&1].state().state().temperature, t1, 2e-10);
    close(next.parcels()[&2].state().state().temperature, t2, 2e-10);
    close(
        next.step().internal_links[0].heat_correction_left_to_right_j,
        internal_q,
        1e-14,
    );
    close(next.step().links[0].reduced_heat_into_gas_j, bath_q, 1e-14);
    close(
        next.step().internal_links[0].face.enthalpy_temperature_k,
        410.0,
        1e-12,
    );
    close(
        next.step().links[0].face.enthalpy_temperature_k,
        390.0,
        1e-12,
    );
    assert!(next.step().exchange.fluxes.is_empty());
    assert!(next.step().exchange.boundaries.is_empty());
    assert_eq!(next.step().exchange.heat.len(), 1);
    assert_eq!(next.step().exchange.external_heat.len(), 1);
    let entropy = c1 * (t1 / 400.0).ln() + c2 * (t2 / 420.0).ln() - bath_q / 390.0;
    close(next.step().audit.entropy_production_j_k, entropy, 2e-14);
    assert!(entropy > 0.0);
    assert!(next.step().audit.energy_residual_j.abs() < 1e-9);
}

#[test]
fn internal_and_bath_species_and_dufour_heat_share_original_compositions() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 0.9, 400.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.7, 400.0, 1e5)),
        (2, parcel(&bank, 0.3, 400.0, 1e5)),
    ]);
    let next = network(&mut op, &original, &[link()], &[boundary(&bath)], 0.1).unwrap();
    let (ji, qi) = analytic(400.0, 1e5, 0.5, -4.0, 0.0);
    let (jb, qb) = analytic(400.0, 1e5, 0.8, -2.0, 0.0);
    let (mi, mb, hi, hb) = (ji * 0.001, jb * 0.001, qi * 0.001, qb * 0.001);
    close(
        next.step().internal_links[0].mass_left_to_right_kg[0],
        mi,
        1e-18,
    );
    close(next.step().links[0].mass_into_gas_kg[0], mb, 1e-18);
    close(
        mass(&next.parcels()[&1], "A") - mass(&original[&1], "A"),
        mb - mi,
        1e-18,
    );
    close(
        mass(&next.parcels()[&2], "A") - mass(&original[&2], "A"),
        mi,
        1e-18,
    );
    for (id, heat) in [(1, hb - hi), (2, hi)] {
        let capacity = original[&id].mass_kg() * 3.0 * R / MASS;
        close(
            next.parcels()[&id].state().state().temperature,
            400.0 + heat / capacity,
            2e-10,
        );
        close(
            next.parcels()[&id].mass_kg(),
            original[&id].mass_kg(),
            1e-18,
        );
    }
    assert_eq!(next.step().exchange.fluxes.len(), 2);
    assert_eq!(next.step().exchange.boundaries.len(), 2);
    assert!(next.step().audit.maximum_species_mass_residual_kg < 1e-13);
    assert!(next.step().audit.energy_residual_j.abs() < 1e-9);
    assert!(next.step().audit.entropy_production_j_k > 0.0);
}

#[test]
fn empty_face_families_preserve_original_closed_and_reservoir_candidates() {
    let (bank, mut op) = fixture(2.0 * MASS, 1e14, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.7, 390.0, 1e5)),
        (2, parcel(&bank, 0.3, 410.0, 1e5)),
    ]);
    let closed = solve(&mut op, &original, &[link()], 0.1).unwrap();
    let combined = network(&mut op, &original, &[link()], &[], 0.1).unwrap();
    assert_eq!(closed.parcels(), combined.parcels());
    assert_eq!(closed.step().links, combined.step().internal_links);
    assert_eq!(closed.step().transport.exchange, combined.step().exchange);
    assert!(combined.step().links.is_empty());
    close(
        closed.step().transport.entropy_change_j_k,
        combined.step().audit.entropy_production_j_k,
        1e-12,
    );
    let bath = bath(&bank, 0.4, 400.0);
    let endpoints = original.keys().map(|&id| (id, endpoint())).collect();
    let only_bath = op
        .reservoir_candidate(
            &original,
            &[boundary(&bath)],
            0.1,
            &endpoints,
            controls(),
            || false,
        )
        .unwrap();
    let combined = network(&mut op, &original, &[], &[boundary(&bath)], 0.1).unwrap();
    assert_eq!(only_bath, combined);
    assert!(combined.step().internal_links.is_empty());
}

#[test]
fn internal_and_external_withdrawals_compete_for_one_original_species_allowance() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 0.2, 400.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.8, 400.0, 1e5)),
        (2, parcel(&bank, 0.2, 400.0, 1e5)),
    ]);
    let retained = original.clone();
    // Each exposure withdraws 22.5% of original A and is admissible alone.
    solve(&mut op, &original, &[link()], 30.0).unwrap();
    network(&mut op, &original, &[], &[boundary(&bath)], 30.0).unwrap();
    // Together the gross 45% exceeds the shared 40% allowance.
    assert!(matches!(
        network(&mut op, &original, &[link()], &[boundary(&bath)], 30.0),
        Err(ParcelMolecularError::Transport(
            GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { id: 1, .. })
        ))
    ));
    assert_eq!(original, retained);
}

#[test]
fn internal_only_endpoints_are_source_qualified_and_late_cancellation_is_atomic() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 1.0, 649.9);
    let original = BTreeMap::from([
        (1, parcel(&bank, 1.0, 649.9, 1e5 * 649.9 / 649.0)),
        (2, parcel(&bank, 1.0, 649.0, 1e5)),
    ]);
    let retained = original.clone();
    network(&mut op, &original, &[link()], &[boundary(&bath)], 0.1).unwrap();
    // Gas 2 has no bath face. Its overlarge internal Fourier increment reaches
    // about 650.17 K, outside the actual heat source but inside calorics.
    assert!(matches!(
        network(&mut op, &original, &[link()], &[boundary(&bath)], 120.0),
        Err(ParcelMolecularError::Face(DissipationError::Material(_)))
    ));
    let endpoints = original.keys().map(|&id| (id, endpoint())).collect();
    let mut polls = 0;
    op.network_candidate(
        &original,
        &[link()],
        &[boundary(&bath)],
        0.1,
        &endpoints,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    let result = op.network_candidate(
        &original,
        &[link()],
        &[boundary(&bath)],
        0.1,
        &endpoints,
        controls(),
        || {
            count += 1;
            count == polls
        },
    );
    assert_eq!(result, Err(ParcelMolecularError::Cancelled));
    assert_eq!(original, retained);
}

#[test]
fn both_face_families_share_physical_ids_and_one_bounded_work_allowance() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 0.5, 400.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.5, 400.0, 1e5)),
        (2, parcel(&bank, 0.5, 400.0, 1e5)),
    ]);
    let duplicate = ReservoirMolecularLink {
        id: link().id,
        ..boundary(&bath)
    };
    assert!(matches!(
        network(&mut op, &original, &[link()], &[duplicate], 0.1),
        Err(ParcelMolecularError::Invalid(
            "continuum and reservoir faces need globally unique IDs and positive geometry"
        ))
    ));
    // Each family fits by itself (2100 and 2400 entries); together 4500 exceeds
    // the existing 4096 bound. Refuse before any rate allocation/reconstruction.
    let internal: Vec<_> = (0..700)
        .map(|id| ParcelMolecularLink { id, ..link() })
        .collect();
    let boundaries: Vec<_> = (10000..10600)
        .map(|id| ReservoirMolecularLink {
            id,
            ..boundary(&bath)
        })
        .collect();
    assert!(matches!(
        network(&mut op, &original, &internal, &boundaries, 0.1),
        Err(ParcelMolecularError::Invalid(
            "bounded reservoir graph, exact endpoint map and positive duration required"
        ))
    ));
}
