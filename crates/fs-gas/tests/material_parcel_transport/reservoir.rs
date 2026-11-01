//! Original public finite-gas/bath faces, including open energy and entropy.

use super::*;
use fs_gas::material::transport::parcel::reservoir::{
    MolecularReservoir, ReservoirMolecularCandidate, ReservoirMolecularLink,
};

fn bath(bank: &ThermochemicalGasModel, x: f64, t: f64, p: f64) -> MolecularReservoir {
    MolecularReservoir::new(
        bank.with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        point(t, p),
        || false,
    )
    .unwrap()
}
fn face(bath: &MolecularReservoir) -> ReservoirMolecularLink<'_> {
    ReservoirMolecularLink {
        id: 7,
        parcel: 1,
        reservoir: bath,
        area_m2: 0.01,
        distance_m: 0.1,
    }
}
fn candidate(
    op: &mut Operator,
    gas: &BTreeMap<u64, ClosedGasParcel>,
    faces: &[ReservoirMolecularLink<'_>],
    dt: f64,
) -> Result<ReservoirMolecularCandidate, ParcelMolecularError> {
    let endpoints = gas.keys().map(|&id| (id, endpoint())).collect();
    op.reservoir_candidate(gas, faces, dt, &endpoints, controls(), || false)
}

#[test]
fn source_soret_dufour_bath_cools_finite_gas_with_positive_total_entropy() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 0.5, 390.0, 1e5);
    let retained = bath.clone();
    let original = BTreeMap::from([(1, parcel(&bank, 0.5, 410.0, 1e5))]);
    let next = candidate(&mut op, &original, &[face(&bath)], 0.1).unwrap();
    let (j, q) = analytic(400.0, 1e5, 0.5, 0.0, 200.0);
    let link = &next.step().links[0];
    assert!(j < 0.0 && q < 0.0);
    close(link.mass_into_gas_kg[0], j * 0.001, 1e-18);
    close(link.mass_into_gas_kg[1], -j * 0.001, 1e-18);
    close(link.reduced_heat_into_gas_j, q * 0.001, 1e-13);
    close(link.face.enthalpy_temperature_k, 390.0, 1e-12);
    assert_eq!(next.step().species, vec![id("A"), id("B")]);
    assert_eq!(&link.reservoir, bath.state());
    assert_eq!(bath, retained);
    let old = &original[&1];
    let new = &next.parcels()[&1];
    let amount = old.mass_kg() / MASS;
    let expected_t = 410.0 + q * 0.001 / (amount * 3.0 * R);
    close(new.state().state().temperature, expected_t, 2e-10);
    close(new.mass_kg(), old.mass_kg(), 1e-18);
    close(mass(new, "A") - mass(old, "A"), j * 0.001, 1e-18);
    // Equal masses and equal Cp make the bath species entropy stream vanish.
    // Independently evaluate finite gas entropy and the actual Q/T_bath term.
    let x = mass(new, "A") / new.mass_kg();
    let mixing = x * x.ln() + (1.0 - x) * (1.0 - x).ln() + 2.0_f64.ln();
    let gas_entropy = amount * R * (3.0 * (expected_t / 410.0).ln() - mixing);
    assert!(gas_entropy < 0.0);
    let expected_entropy = gas_entropy - q * 0.001 / 390.0;
    assert!(expected_entropy > 0.0);
    close(
        next.step().audit.entropy_production_j_k,
        expected_entropy,
        2e-14,
    );
    assert!(next.step().audit.energy_residual_j.abs() < 1e-9);
    assert_eq!(next.step().exchange.boundaries.len(), 2);
    assert_eq!(next.step().exchange.external_heat.len(), 1);
}

#[test]
fn bath_enthalpy_not_midpoint_or_donor_enthalpy_sets_unequal_mass_temperature() {
    let run = |formation, heat_shift| {
        let (bank, mut op) = fixture(2.0 * MASS, formation, heat_shift);
        let bath = bath(&bank, 0.7, 420.0, 1e5);
        let original = BTreeMap::from([(1, parcel(&bank, 0.3, 380.0, 1e5))]);
        let next = candidate(&mut op, &original, &[face(&bath)], 0.1).unwrap();
        let mean_mass = 1.5 * MASS;
        let centered_a_heat = (2.0 / 3.0) * 100000.0;
        let driving = -4.0 + 0.5 * MASS * centered_a_heat / (R * 400.0 * 400.0) * -400.0;
        let j = -1e5 / (R * 400.0) * DIFFUSIVITY * (MASS * 2.0 * MASS / mean_mass) * driving;
        // The original material law uses WMS mixture conductivity. Equal pure
        // conductivities/viscosities do not cancel the unequal-mass Wilke
        // weights. Independently evaluate them at the 50/50 midpoint.
        let phi_ab = (1.0 + 2.0_f64.sqrt().sqrt()).powi(2) / 12.0_f64.sqrt();
        let phi_ba = (1.0 + 0.5_f64.sqrt().sqrt()).powi(2) / 24.0_f64.sqrt();
        let conductivity = CONDUCTIVITY / (1.0 + phi_ab) + CONDUCTIVITY / (1.0 + phi_ba);
        let q = conductivity * 400.0 + 100000.0 * j;
        let dm = j * 0.001;
        let dn = dm / MASS - dm / (2.0 * MASS);
        let n0 = 1e5 * 1e-3 / (R * 380.0);
        // Both species have Cv=3R and Cp=4R per mole. Bath h is evaluated at
        // 420 K for BOTH signs; using donor or midpoint h gives a different T.
        let expected_t =
            (3.0 * R * n0 * 380.0 + 4.0 * R * 420.0 * dn + q * 0.001) / (3.0 * R * (n0 + dn));
        let report = &next.step().links[0];
        close(report.mass_into_gas_kg[0], dm, 2e-18);
        close(report.reduced_heat_into_gas_j, q * 0.001, 1e-13);
        close(report.face.enthalpy_temperature_k, 420.0, 1e-12);
        close(
            next.parcels()[&1].state().state().temperature,
            expected_t,
            2e-10,
        );
        for boundary in &next.step().exchange.boundaries {
            let mm = if boundary.transfer.species == id("A") {
                MASS
            } else {
                2.0 * MASS
            };
            close(
                boundary.transfer.relative_energy_into_gas_j,
                boundary.transfer.net_mass_into_gas_kg * 4.0 * R * 40.0 / mm,
                1e-14,
            );
        }
        assert!(next.step().audit.entropy_production_j_k > 0.0);
        assert!(next.step().audit.energy_residual_j.abs() < 1e-9);
        (
            next.parcels()[&1].state().state().temperature,
            report.mass_into_gas_kg.clone(),
            next.step().audit,
        )
    };
    let plain = run(0.0, 0.0);
    let shifted = run(1e14, 1e12);
    close(plain.0, shifted.0, 1e-10);
    for (a, b) in plain.1.into_iter().zip(shifted.1) {
        close(a, b, 1e-18);
    }
    close(
        plain.2.energy_residual_j,
        shifted.2.energy_residual_j,
        1e-10,
    );
    close(
        plain.2.entropy_production_j_k,
        shifted.2.entropy_production_j_k,
        1e-12,
    );
}

#[test]
fn pure_species_bath_delivers_heat_without_an_artificial_species_inventory() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 1.0, 440.0, 1e5);
    let original = BTreeMap::from([(1, parcel(&bank, 1.0, 400.0, 1e5))]);
    let next = candidate(&mut op, &original, &[face(&bath)], 0.1).unwrap();
    let q = CONDUCTIVITY * 400.0 * 0.001;
    let capacity = original[&1].mass_kg() * 3.0 * R / MASS;
    assert!(next.step().exchange.boundaries.is_empty());
    assert_eq!(next.step().links[0].mass_into_gas_kg, vec![0.0, 0.0]);
    close(next.step().links[0].reduced_heat_into_gas_j, q, 1e-14);
    close(
        next.parcels()[&1].state().state().temperature,
        400.0 + q / capacity,
        2e-10,
    );
    close(
        mass(&next.parcels()[&1], "A"),
        original[&1].mass_kg(),
        1e-18,
    );
    assert_eq!(mass(&next.parcels()[&1], "B"), 0.0);
    assert!(next.step().audit.entropy_production_j_k > 0.0);
}

#[test]
fn opposing_baths_share_original_gross_withdrawal_and_are_order_independent() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let rich = bath(&bank, 0.9, 400.0, 1e5);
    let lean = bath(&bank, 0.1, 400.0, 1e5);
    let links = [
        face(&rich),
        ReservoirMolecularLink {
            id: 8,
            ..face(&lean)
        },
    ];
    let original = BTreeMap::from([(1, parcel(&bank, 0.5, 400.0, 1e5))]);
    let retained = original.clone();
    let first = candidate(&mut op, &original, &links, 0.1).unwrap();
    let second = candidate(&mut op, &original, &[links[1], links[0]], 0.1).unwrap();
    assert_eq!(first, second);
    close(
        mass(&first.parcels()[&1], "A"),
        mass(&original[&1], "A"),
        1e-18,
    );
    close(
        first.parcels()[&1].state().state().temperature,
        400.0,
        1e-10,
    );
    assert_eq!(first.step().exchange.boundaries.len(), 4);
    assert!(first.step().audit.entropy_production_j_k > 0.0);
    // Net mass and heat almost cancel, but each gross outflow would withdraw
    // 64% of one ORIGINAL species. The opposite bath cannot replenish it.
    assert!(matches!(
        candidate(&mut op, &original, &links, 80.0),
        Err(ParcelMolecularError::Transport(
            GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { .. })
        ))
    ));
    assert_eq!(original, retained);
}

#[test]
fn absent_finite_species_can_enter_but_nonzero_stream_at_zero_bath_activity_refuses() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let mixed = bath(&bank, 0.5, 400.0, 1e5);
    let pure = BTreeMap::from([(1, parcel(&bank, 1.0, 400.0, 1e5))]);
    let next = candidate(&mut op, &pure, &[face(&mixed)], 0.01).unwrap();
    assert!(mass(&next.parcels()[&1], "B") > 0.0);
    assert!(next.step().audit.entropy_production_j_k > 0.0);
    let pure_bath = bath(&bank, 1.0, 400.0, 1e5);
    let gas = BTreeMap::from([(1, parcel(&bank, 0.5, 400.0, 1e5))]);
    assert!(matches!(
        candidate(&mut op, &gas, &[face(&pure_bath)], 0.01),
        Err(ParcelMolecularError::Invalid(
            "nonzero reservoir stream requires positive bath species activity"
        ))
    ));
    assert_eq!(pure[&1].model().thermodynamics().components().len(), 1);
}

#[test]
fn reservoir_initial_and_final_sources_and_late_cancellation_are_atomic() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([(1, parcel(&bank, 0.5, 400.0, 1e5))]);
    let retained = original.clone();
    // The midpoint 550 K is supported; the held 700 K bath is outside the
    // authored heat-of-transport domain and must be refused independently.
    let unsupported = bath(&bank, 0.5, 700.0, 1e5);
    assert!(matches!(
        candidate(&mut op, &original, &[face(&unsupported)], 0.1),
        Err(ParcelMolecularError::Face(DissipationError::Material(_)))
    ));
    // Both original traces and midpoint are supported. The explicit endpoint
    // overshoots the 650 K source boundary while remaining inside calorics.
    let edge = bath(&bank, 0.5, 649.9, 1e5);
    let near_edge = BTreeMap::from([(1, parcel(&bank, 0.5, 649.0, 1e5))]);
    candidate(&mut op, &near_edge, &[face(&edge)], 0.1).unwrap();
    assert!(matches!(
        candidate(&mut op, &near_edge, &[face(&edge)], 120.0),
        Err(ParcelMolecularError::Face(DissipationError::Material(_)))
    ));
    let supported = bath(&bank, 0.5, 420.0, 1e5);
    let endpoints = BTreeMap::from([(1, endpoint())]);
    let mut polls = 0;
    op.reservoir_candidate(
        &original,
        &[face(&supported)],
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
    let result = op.reservoir_candidate(
        &original,
        &[face(&supported)],
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
fn supported_heat_overshoot_refuses_actual_negative_gas_plus_bath_entropy() {
    use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::DiffusionNetworkError;

    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let bath = bath(&bank, 1.0, 420.0, 1e5);
    let original = BTreeMap::from([(1, parcel(&bank, 1.0, 400.0, 1e5))]);
    let retained = original.clone();
    candidate(&mut op, &original, &[face(&bath)], 0.1).unwrap();
    // Constant Fourier conductance G=lambda*A/d=0.005 W/K delivers 100 J.
    // The resulting ~533 K is within every source/caloric domain, but exceeds
    // the supported physical entropy budget for this explicit time interval.
    let q = CONDUCTIVITY * 0.01 / 0.1 * 20.0 * 1000.0;
    let capacity = original[&1].mass_kg() * 3.0 * R / MASS;
    let final_t = 400.0 + q / capacity;
    assert!((500.0..650.0).contains(&final_t));
    let expected = capacity * (final_t / 400.0).ln() - q / 420.0;
    assert!(expected < -controls().physical.entropy_tolerance_j_k);
    match candidate(&mut op, &original, &[face(&bath)], 1000.0) {
        Err(ParcelMolecularError::Transport(GasTransportError::Diffusion(
            DiffusionNetworkError::EntropyDecrease {
                change_j_k,
                tolerance_j_k,
            },
        ))) => {
            close(change_j_k, expected, 1e-12);
            assert_eq!(tolerance_j_k, controls().physical.entropy_tolerance_j_k);
        }
        other => panic!("expected actual open-system entropy refusal, got {other:?}"),
    }
    assert_eq!(original, retained);
}
