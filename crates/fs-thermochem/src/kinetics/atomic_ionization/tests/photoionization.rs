use super::{channel, close};
use crate::electron_collision::AVOGADRO_PER_MOL;
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
    let mut points = source.points().to_vec();
    points[0].photon_energy = Energy::new(0.9 * source.threshold().value());
    assert!(AtomicPhotoionizationV1::new(
        channel(), source.provenance().clone(), points,
    ).is_err());
    let mut points = source.points().to_vec();
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
