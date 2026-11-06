use super::{channel, close};
use crate::electron_collision::{AVOGADRO_PER_MOL, ELECTRONVOLT_J};
use crate::kinetics::atomic_ionization::photoionization::*;
use fs_matdb::Provenance;
use fs_qty::{Area, Energy};

fn source() -> AtomicPhotoionizationV1 {
    let atomic = channel();
    let gap = atomic.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    AtomicPhotoionizationV1::new(
        atomic,
        Provenance {
            source: "synthetic linear event table for independent algebra tests".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap),
                cross_section: Area::new(2e-22),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(2.0 * gap),
                cross_section: Area::new(6e-22),
            },
        ],
    ).unwrap()
}

#[test]
fn photoedge_is_nonzero_below_edge_is_zero_and_upper_support_refuses() {
    let source = source();
    let gap = source.threshold().value();
    assert_eq!(source.channel(), &channel());
    assert_eq!(source.cross_section(Energy::new(0.5 * gap)).unwrap().value(), 0.0);
    assert_eq!(source.cross_section(source.threshold()).unwrap().value(), 2e-22);
    close(source.cross_section(Energy::new(1.5 * gap)).unwrap().value(),
        4e-22, 2e-15, 0.0);
    assert_eq!(source.cross_section(source.maximum_energy()).unwrap().value(), 6e-22);
    assert_eq!(source.cross_section(Energy::new(2.01 * gap)),
        Err(PhotoionizationErrorV1::OutsideEnergyDomain));
}

#[test]
fn discrete_number_flux_moments_close_original_chemical_and_electron_energy() {
    let source = source();
    let gap = source.threshold().value();
    let result = source.evaluate(&[
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(gap), number_flux_per_m2_s: 3e20,
        },
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(1.5 * gap), number_flux_per_m2_s: 2e20,
        },
    ], || false).unwrap();
    // Independent event algebra: 0.06 + 0.08 events/(atom*s).
    close(result.rate_per_neutral_per_s, 0.14, 3e-15, 0.0);
    close(result.absorbed_power_per_neutral_w, 0.18 * gap, 3e-15, 0.0);
    close(result.chemical_power_per_neutral_w, 0.14 * gap, 3e-15, 0.0);
    close(result.electron_heating_per_neutral_w, 0.04 * gap, 3e-15, 0.0);
    close(result.absorbed_power_per_neutral_w,
        result.chemical_power_per_neutral_w + result.electron_heating_per_neutral_w,
        3e-15, 0.0);
    assert_eq!(source.evaluate(&[], || false).unwrap().rate_per_neutral_per_s, 0.0);
}

#[test]
fn source_gap_order_and_local_field_admission_refuse() {
    let source = source();
    let mut points = source.points().unwrap().to_vec();
    points[0].photon_energy = Energy::new(0.9 * source.threshold().value());
    assert!(AtomicPhotoionizationV1::new(
        channel(), source.provenance().clone(), points,
    ).is_err());
    let mut points = source.points().unwrap().to_vec();
    points.swap(0, 1);
    assert!(AtomicPhotoionizationV1::new(
        channel(), source.provenance().clone(), points,
    ).is_err());
    for flux in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(source.evaluate(&[PhotonNumberFluxV1 {
            photon_energy: source.threshold(), number_flux_per_m2_s: flux,
        }], || false).is_err());
    }
    assert_eq!(source.evaluate(&[PhotonNumberFluxV1 {
        photon_energy: Energy::new(source.maximum_energy().value() * 2.0),
        number_flux_per_m2_s: 0.0,
    }], || false), Err(PhotoionizationErrorV1::OutsideEnergyDomain));
}

#[test]
fn cancellation_is_observed_inside_group_work_without_a_partial_result() {
    let source = source();
    let groups = [PhotonNumberFluxV1 {
        photon_energy: source.threshold(), number_flux_per_m2_s: 1e20,
    }; 65];
    let mut polls = 0;
    let result = source.evaluate(&groups, || { polls += 1; polls == 2 });
    assert_eq!(result, Err(PhotoionizationErrorV1::Cancelled));
    assert_eq!(polls, 2);
}

#[test]
fn verner_hydrogen_matches_independent_published_parameter_anchors() {
    let source = AtomicPhotoionizationV1::verner_hydrogen(channel()).unwrap();
    // Direct-power evaluation of the literal Table 1 H I row, independent of
    // the production logarithmic operation tree; these are numerical anchors,
    // not measured-data uncertainties or a bound on the fit's physical error.
    for (energy_ev, expected_m2) in [
        (13.6, 6.346_296_358_990_503e-22),
        (20.0, 2.211_102_984_039_216_3e-22),
        (100.0, 1.939_792_233_567_902_5e-24),
        (1000.0, 1.140_867_331_110_972e-27),
        (50_000.0, 1.715_120_863_209_725_6e-33),
    ] {
        close(
            source.cross_section(Energy::new(energy_ev * ELECTRONVOLT_J)).unwrap().value(),
            expected_m2, 1e-13, 0.0,
        );
    }
}

#[test]
fn verner_hydrogen_keeps_published_edge_and_upper_energy_domain() {
    let source = AtomicPhotoionizationV1::verner_hydrogen(channel()).unwrap();
    assert!(source.points().is_none());
    let edge = source.threshold().value();
    let maximum = source.maximum_energy().value();
    assert_eq!(edge, 13.60 * ELECTRONVOLT_J);
    assert_eq!(maximum, 50_000.0 * ELECTRONVOLT_J);
    assert_eq!(
        source.cross_section(Energy::new(f64::from_bits(edge.to_bits() - 1))).unwrap().value(),
        0.0,
    );
    assert!(source.cross_section(Energy::new(edge)).unwrap().value() > 0.0);
    assert!(source.cross_section(Energy::new(maximum)).unwrap().value() > 0.0);
    assert_eq!(
        source.cross_section(Energy::new(f64::from_bits(maximum.to_bits() + 1))),
        Err(PhotoionizationErrorV1::OutsideEnergyDomain),
    );
}

#[test]
fn verner_hydrogen_refuses_wrong_element_mass_gap_and_statistics() {
    use super::{
        AtomicIonizationV1, AtomicStatisticalWeightsV1, ChargeVector, ElementId,
        ElementalMatrix, IdealGasEquilibriumV1, ReactionId,
    };
    let models = super::models();
    let ids = models.iter().map(|m| m.species().clone()).collect::<Vec<_>>();
    // Species still named H/H+/e, but their declared nuclear row is helium.
    let wrong_element = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("He").unwrap()], ids.clone(), vec![vec![1, 1, 0]],
        ).unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        models,
    ).unwrap();
    let atomic = AtomicIonizationV1::new(
        wrong_element, ReactionId::new("wrong-element").unwrap(), super::id("H+"),
        super::impact(), super::WEIGHTS, 1e-12,
    ).unwrap();
    assert!(AtomicPhotoionizationV1::verner_hydrogen(atomic).is_err());

    let mut models = super::models();
    models[0] = super::source(
        "H", 2.0 * super::MP + super::ELECTRON_MASS_KG, 2.0, 0.0, 2.5, 200_000.0,
    );
    models[1] = super::source("H+", 2.0 * super::MP, 1.0, super::CHI, 2.5, 200_000.0);
    assert!(AtomicPhotoionizationV1::verner_hydrogen(super::bind(models).unwrap()).is_err());

    let mut models = super::models();
    models[1] = super::source("H+", super::MP, 1.0, 0.95 * super::CHI, 2.5, 200_000.0);
    assert!(AtomicPhotoionizationV1::verner_hydrogen(super::bind(models).unwrap()).is_err());

    let mut models = super::models();
    models[0] = super::source("H", super::MA, 4.0, 0.0, 2.5, 200_000.0);
    let atomic = AtomicIonizationV1::new(
        super::bank(models, 1), ReactionId::new("wrong-statistics").unwrap(), super::id("H+"),
        super::impact(), AtomicStatisticalWeightsV1 { neutral: 4, ion: 1 }, 1e-12,
    ).unwrap();
    assert!(AtomicPhotoionizationV1::verner_hydrogen(atomic).is_err());
}

#[test]
fn verner_heating_uses_original_gap_instead_of_rounded_fit_threshold() {
    let atomic = channel();
    let gap = atomic.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    let source = AtomicPhotoionizationV1::verner_hydrogen(atomic.clone()).unwrap();
    assert_eq!(source.channel(), &atomic);
    let energy = source.threshold().value();
    let result = source.evaluate(&[PhotonNumberFluxV1 {
        photon_energy: Energy::new(energy), number_flux_per_m2_s: 1e20,
    }], || false).unwrap();
    assert!(energy > gap);
    assert!(result.electron_heating_per_neutral_w > 0.0);
    close(result.chemical_power_per_neutral_w,
        result.rate_per_neutral_per_s * gap, 2e-15, 0.0);
    close(result.electron_heating_per_neutral_w,
        result.rate_per_neutral_per_s * (energy - gap), 2e-15, 0.0);
    close(result.absorbed_power_per_neutral_w,
        result.chemical_power_per_neutral_w + result.electron_heating_per_neutral_w,
        2e-15, 0.0);
}
