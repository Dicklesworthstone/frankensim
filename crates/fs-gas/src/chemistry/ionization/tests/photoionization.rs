//! G0/G1/G3: analytic photon inventories, source energy and composition semantics.
//! Cross sections here are manufactured constants, not measured gas data.
use super::*;
use crate::chemistry::photoionization::{
    PhotoionizationControls, PhotoionizationError, PhotoionizationExposure,
    advance_monochromatic_photoionization,
};
use fs_qty::{Area, Energy};
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1, PhotonNumberFluxV1,
};

const SIGMA: f64 = 1e-22;

fn photon_source(channel: AtomicIonizationV1) -> AtomicPhotoionizationV1 {
    AtomicPhotoionizationV1::new(channel, Provenance {
        source: "constant cross-section analytic photon-budget fixture".into(),
        license: "analytic test fixture".into(), artifact: None,
    }, vec![
        PhotoionizationCrossSectionPointV1 {
            photon_energy: Energy::new(14.0 * EV), cross_section: Area::new(SIGMA),
        },
        PhotoionizationCrossSectionPointV1 {
            photon_energy: Energy::new(24.0 * EV), cross_section: Area::new(SIGMA),
        },
    ]).unwrap()
}
fn photo_controls() -> PhotoionizationControls {
    PhotoionizationControls { relative_charge_tolerance: 1e-12, energy_tolerance_j: 1e-7 }
}
fn exposure(neutral: f64, tau: f64, dose: f64) -> PhotoionizationExposure {
    PhotoionizationExposure {
        beam: PhotonNumberFluxV1 {
            photon_energy: Energy::new(16.0 * EV), number_flux_per_m2_s: dose / SIGMA,
        },
        illuminated_area: Area::new(SIGMA * NA * neutral / tau),
        duration: Time::new(1.0),
    }
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e}");
}

#[test]
fn photons_create_first_electrons_and_close_original_source_energy() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic);
    let result = advance_monochromatic_photoionization(
        &initial, &source, exposure(1.0, 1.0, 1.0), photo_controls(), || false,
    ).unwrap();
    // Direct closed form at tau0=dose=1, independent of the stable production branches.
    let remaining = (2.0 - (-1.0_f64).exp()).ln();
    let extent = 1.0 - remaining;
    let end = result.state();
    relative(result.reaction_extent().value(), extent, 2e-14);
    relative(end.amounts()[0].value(), remaining, 2e-14);
    relative(end.amounts()[1].value(), extent, 2e-14);
    assert_eq!(end.amounts()[1], end.amounts()[2]);
    relative(end.mass().value(), initial.mass().value(), 2e-14);
    assert_eq!(end.heavy_temperature(), initial.heavy_temperature());
    let born = (16.0 * EV * NA - CHI) / (1.5 * R);
    relative(end.electron_temperature().value(), born, 2e-14);
    let b = result.balance();
    relative(b.absorbed_photons, NA * extent, 2e-14);
    relative(b.absorbed_photons + b.transmitted_photons, b.incident_photons, 2e-14);
    relative(b.chemical_storage_change_j, CHI * extent, 2e-14);
    relative(b.electron_sensible_energy_change_j, 1.5 * R * extent * born, 2e-14);
    assert!((b.absorbed_energy_j - b.chemical_storage_change_j
        - b.electron_sensible_energy_change_j).abs() < 1e-7);
    assert!((end.internal_energy().value() - initial.internal_energy().value()
        - b.absorbed_energy_j).abs() < 1e-7);
    assert_eq!(end.thermodynamics(), initial.thermodynamics());
}

#[test]
fn thin_column_depletes_exponentially_and_thick_column_is_photon_limited() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic);
    let thin = advance_monochromatic_photoionization(
        &initial, &source, exposure(1.0, 1e-12, 0.25), photo_controls(), || false,
    ).unwrap();
    relative(thin.state().amounts()[0].value(), (-0.25_f64).exp(), 1e-12);
    relative(thin.reaction_extent().value(), -(-0.25_f64).exp_m1(), 1e-12);
    let thick = advance_monochromatic_photoionization(
        &initial, &source, exposure(1.0, 1e6, 100.0), photo_controls(), || false,
    ).unwrap();
    relative(thick.reaction_extent().value(), 1e-4, 2e-14);
    relative(thick.balance().absorbed_photons, thick.balance().incident_photons, 2e-14);
    assert!(thick.balance().transmitted_photons >= 0.0);
}

#[test]
fn exact_depletion_preserves_small_events_and_small_surviving_populations() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic);
    let small = advance_monochromatic_photoionization(
        &initial, &source, exposure(1.0, 1.0, 1e-12), photo_controls(), || false,
    ).unwrap();
    relative(small.reaction_extent().value(), -(-1.0_f64).exp_m1() * 1e-12, 1e-12);
    let depleted = advance_monochromatic_photoionization(
        &initial, &source, exposure(1.0, 1.0, 100.0), photo_controls(), || false,
    ).unwrap();
    let surviving = 1.0_f64.exp_m1() * (-100.0_f64).exp();
    relative(depleted.state().amounts()[0].value(), surviving, 5e-14);
    assert!(depleted.state().amounts()[0].value() > 0.0);
    assert!(depleted.balance().absorbed_photons <= depleted.balance().incident_photons);
}

#[test]
fn repeated_exposure_has_the_exact_semigroup_and_retains_temperature_dilution() {
    let (initial, atomic) = initial(0.1, 30000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic);
    let full_exposure = exposure(0.9, 3.0, 2.0);
    let full = advance_monochromatic_photoionization(
        &initial, &source, full_exposure, photo_controls(), || false,
    ).unwrap();
    let mut half_exposure = full_exposure;
    half_exposure.duration = Time::new(0.5);
    let first = advance_monochromatic_photoionization(
        &initial, &source, half_exposure, photo_controls(), || false,
    ).unwrap();
    let second = advance_monochromatic_photoionization(
        first.state(), &source, half_exposure, photo_controls(), || false,
    ).unwrap();
    for (a, b) in full.state().amounts().iter().zip(second.state().amounts()) {
        relative(a.value(), b.value(), 5e-14);
    }
    relative(full.state().electron_temperature().value(),
        second.state().electron_temperature().value(), 5e-14);
    relative(full.balance().absorbed_photons,
        first.balance().absorbed_photons + second.balance().absorbed_photons, 5e-14);
    // A positive absorbed energy can lower Te: the new electrons are cooler.
    assert!(full.balance().absorbed_energy_j > 0.0);
    assert!(full.state().electron_temperature().value() < initial.electron_temperature().value());
}

#[test]
fn photon_energy_account_is_invariant_under_conserved_source_reference_shifts() {
    let mut outputs = Vec::new();
    for (elemental, charge) in [(0.0, 0.0), (1e12, 1e8)] {
        let (initial, atomic) = initial(0.05, 20000.0, 4000.0, elemental, charge);
        let source = photon_source(atomic);
        let result = advance_monochromatic_photoionization(
            &initial, &source, exposure(0.95, 2.0, 0.5), photo_controls(), || false,
        ).unwrap();
        assert!(result.balance().energy_residual_j.abs() <= 1e-7);
        outputs.push(result);
    }
    assert_eq!(outputs[0].state().amounts(), outputs[1].state().amounts());
    relative(outputs[0].state().electron_temperature().value(),
        outputs[1].state().electron_temperature().value(), 1e-9);
}

#[test]
fn below_threshold_zero_flux_and_zero_time_preserve_state_and_photon_accounts() {
    let (initial, atomic) = initial(0.1, 20000.0, 4000.0, 0.0, 0.0);
    let source = photon_source(atomic);
    let original = exposure(0.9, 1.0, 1.0);
    let mut below = original;
    below.beam.photon_energy = Energy::new(13.0 * EV);
    let mut dark = original;
    dark.beam.number_flux_per_m2_s = 0.0;
    let mut zero_time = original;
    zero_time.duration = Time::new(0.0);
    for illumination in [below, dark, zero_time] {
        let result = advance_monochromatic_photoionization(
            &initial, &source, illumination, photo_controls(), || false,
        ).unwrap();
        assert_eq!(result.state(), &initial);
        assert_eq!(result.reaction_extent().value(), 0.0);
        assert_eq!(result.balance().absorbed_energy_j, 0.0);
        assert_eq!(result.balance().absorbed_photons, 0.0);
        assert_eq!(result.balance().incident_photons, result.balance().transmitted_photons);
    }
}

#[test]
fn source_binding_temperature_domains_and_cancellation_refuse_atomically() {
    let (initial, atomic) = initial(0.0, 12000.0, 4000.0, 0.0, 0.0);
    let before = initial.clone();
    let source = photon_source(atomic);
    let beam = exposure(1.0, 1.0, 1.0);
    assert_eq!(advance_monochromatic_photoionization(
        &initial, &source, beam, photo_controls(), || true,
    ), Err(PhotoionizationError::Cancelled));
    let mut calls = 0;
    assert_eq!(advance_monochromatic_photoionization(
        &initial, &source, beam, photo_controls(), || { calls += 1; calls >= 3 },
    ), Err(PhotoionizationError::Cancelled));
    let (_, other_atomic) = initial_with_shift();
    assert!(matches!(advance_monochromatic_photoionization(
        &initial, &photon_source(other_atomic), beam, photo_controls(), || false,
    ), Err(PhotoionizationError::Invalid(_))));
    let mut unsupported = beam;
    unsupported.beam.photon_energy = Energy::new(25.0 * EV);
    assert!(matches!(advance_monochromatic_photoionization(
        &initial, &source, unsupported, photo_controls(), || false,
    ), Err(PhotoionizationError::Source(_))));
    let mut hot = beam;
    hot.beam.photon_energy = Energy::new(22.0 * EV);
    // 8.4 eV of excess per newborn electron exceeds this fixture's 50,000 K source.
    assert!(matches!(advance_monochromatic_photoionization(
        &initial, &source, hot, photo_controls(), || false,
    ), Err(PhotoionizationError::State(_))));
    assert_eq!(initial, before);
}
fn initial_with_shift() -> (FrozenIonizedGasStateV1, AtomicIonizationV1) {
    initial(0.0, 12000.0, 4000.0, 1e5, 0.0)
}
