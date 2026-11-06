//! Analytic Maxwellian and independent published-source checks for Milne capture.

use super::{channel, close};
use crate::electron_collision::{
    AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K as KB, MaxwellianRateControlsV1,
};
use crate::kinetics::atomic_ionization::{
    photoionization::{AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1},
    radiative_recombination::*,
};
use fs_matdb::Provenance;
use fs_qty::{Area, Energy, Temperature};

fn provenance() -> Provenance {
    Provenance {
        source: "explicit uniform all-energy envelope for analytic/source tests".into(),
        license: "analytic test fixture".into(),
        artifact: None,
    }
}
fn controls() -> MaxwellianRateControlsV1 {
    MaxwellianRateControlsV1 {
        relative_tolerance: 1e-12,
        maximum_subdivisions: 4096,
    }
}
fn table(points: &[(f64, f64)], bound: f64) -> AtomicRadiativeRecombinationV1 {
    let photo = AtomicPhotoionizationV1::new(
        channel(),
        provenance(),
        points
            .iter()
            .map(|&(energy, sigma)| PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(energy),
                cross_section: Area::new(sigma),
            })
            .collect(),
    )
    .unwrap();
    AtomicRadiativeRecombinationV1::new(photo, Area::new(bound), provenance()).unwrap()
}
fn constant(temperature: f64, shift: f64, coverage: f64) -> AtomicRadiativeRecombinationV1 {
    let gap = channel().ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    let thermal = KB * temperature;
    let onset = gap + shift * thermal;
    table(
        &[(onset, 1e-22), (onset + coverage * thermal, 1e-22)],
        1e-22,
    )
}

#[test]
fn constant_photo_cross_section_matches_independent_complete_maxwellian_moments() {
    // 75-digit decimal evaluation of the closed integrals:
    // alpha=C*sigma/sqrt(theta)*(chi²+2chi*theta+2theta²)
    // L=C*sigma*sqrt(theta)*(chi²+4chi*theta+6theta²).
    let source = constant(10_000.0, 0.0, 100.0);
    let r = source
        .evaluate(Temperature::new(10_000.0), controls(), || false)
        .unwrap();
    close(r.coefficient_m3_s, 2.960_469_529_028_938e-20, 2e-12, 0.0);
    close(
        r.captured_kinetic_energy_j_m3_s,
        4.601_736_332_018_349e-39,
        2e-12,
        0.0,
    );
    close(
        r.mean_captured_energy_j().unwrap(),
        1.554_394_087_456_715_6e-19,
        2e-12,
        0.0,
    );
    assert!((r.mean_captured_energy_j().unwrap() / (1.5 * KB * 10_000.0) - 1.0).abs() > 0.2);
    let gap = source.source().channel().ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    close(
        r.emitted_photon_energy_j_m3_s,
        gap * r.coefficient_m3_s + r.captured_kinetic_energy_j_m3_s,
        2e-15,
        0.0,
    );
    assert!(r.relative_quadrature_error_estimate <= controls().relative_tolerance * 0.25);
    assert!(r.relative_tail_bound <= controls().relative_tolerance * 0.25);
    assert!(r.relative_tail_bound > 0.0);
}

#[test]
fn explicit_onset_shift_retains_boltzmann_suppression_and_capture_energy() {
    // The complete constant-sigma integral starts at epsilon=5*kT.
    let source = constant(10_000.0, 5.0, 100.0);
    let r = source
        .evaluate(Temperature::new(10_000.0), controls(), || false)
        .unwrap();
    close(r.coefficient_m3_s, 3.355_761_106_194_091e-22, 2e-12, 0.0);
    close(
        r.captured_kinetic_energy_j_m3_s,
        2.824_280_809_695_143e-40,
        2e-12,
        0.0,
    );
    close(
        r.mean_captured_energy_j().unwrap(),
        8.416_215_339_292_366e-19,
        2e-12,
        0.0,
    );
    assert!(r.mean_captured_energy_j().unwrap() > 5.0 * KB * 10_000.0);
}

#[test]
fn narrow_source_knots_and_million_thermal_energy_coverage_do_not_hide_the_peak() {
    let gap = channel().ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    let theta = KB * 10_000.0;
    let source = table(
        &[
            (gap, 0.0),
            (gap + 0.99 * theta, 0.0),
            (gap + theta, 1e-22),
            (gap + 1.01 * theta, 0.0),
            (gap + 1e6 * theta, 0.0),
        ],
        1e-22,
    );
    let r = source
        .evaluate(Temperature::new(10_000.0), controls(), || false)
        .unwrap();
    // Independent 75-digit integration of each polynomial using integer
    // incomplete gamma functions, not sampling a second numerical quadrature.
    close(r.coefficient_m3_s, 1.085_248_704_843_075_6e-22, 2e-11, 0.0);
    close(
        r.captured_kinetic_energy_j_m3_s,
        1.498_325_543_061_298_2e-41,
        2e-11,
        0.0,
    );
}

#[test]
fn original_verner_hydrogen_milne_moments_match_independent_simpson_oracles() {
    let photo = AtomicPhotoionizationV1::verner_hydrogen(channel()).unwrap();
    let source =
        AtomicRadiativeRecombinationV1::new(photo.clone(), Area::new(6.4e-22), provenance())
            .unwrap();
    assert_eq!(source.source(), &photo);
    // Composite Simpson with Richardson extrapolation, 2048..16384 intervals
    // over 80 thermal energies; final changes <1.2e-13 relative for both moments.
    // Verner's original 13.60 eV onset and original NIST thermochemical gap differ.
    for (t, alpha, kinetic) in [
        (1000.0, 5.120_408_283_482_090e-19, 7.167_514_554_236_459e-39),
        (
            10_000.0,
            1.586_467_063_069_108_3e-19,
            2.109_062_355_317_087_3e-38,
        ),
        (
            100_000.0,
            3.885_042_045_874_216e-20,
            4.176_374_693_530_577e-38,
        ),
    ] {
        let r = source
            .evaluate(Temperature::new(t), controls(), || false)
            .unwrap();
        close(r.coefficient_m3_s, alpha, 2e-11, 0.0);
        close(r.captured_kinetic_energy_j_m3_s, kinetic, 2e-11, 0.0);
    }
}

#[test]
fn finite_source_coverage_must_resolve_energy_as_well_as_event_tail() {
    let source = constant(10_000.0, 0.0, 0.2);
    assert!(matches!(
        source.evaluate(Temperature::new(10_000.0), controls(), || false),
        Err(RadiativeRecombinationErrorV1::UnresolvedTail { .. })
    ));
    // At 22*kT the event tail passes tol/4 but the kinetic-energy tail does not.
    let source = constant(10_000.0, 0.0, 22.0);
    let control = MaxwellianRateControlsV1 {
        relative_tolerance: 1e-7,
        maximum_subdivisions: 4096,
    };
    let a = source.source().threshold().value() / (KB * 10_000.0);
    let x = a + 22.0;
    let rate_tail = (-22.0_f64).exp() * (x * x + 2.0 * x + 2.0);
    let energy_tail = (-22.0_f64).exp() * (22.0 * (x * x + 2.0 * x + 2.0) + x * x + 4.0 * x + 6.0);
    assert!(rate_tail / (a * a + 2.0 * a + 2.0 - rate_tail) < control.relative_tolerance * 0.25);
    assert!(
        energy_tail / (a * a + 4.0 * a + 6.0 - energy_tail) > control.relative_tolerance * 0.25
    );
    assert!(matches!(
        source.evaluate(Temperature::new(10_000.0), control, || false),
        Err(RadiativeRecombinationErrorV1::UnresolvedTail { .. })
    ));
}

#[test]
fn global_envelope_and_zero_channel_are_explicit() {
    let source = constant(10_000.0, 0.0, 100.0);
    for bound in [-1.0, f64::NAN, f64::INFINITY, 0.9e-22] {
        assert!(
            AtomicRadiativeRecombinationV1::new(
                source.source().clone(),
                Area::new(bound),
                provenance(),
            )
            .is_err()
        );
    }
    let mut missing = provenance();
    missing.source.clear();
    assert!(
        AtomicRadiativeRecombinationV1::new(source.source().clone(), Area::new(1e-22), missing,)
            .is_err()
    );
    let gap = source.source().channel().ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    let dark = table(&[(gap, 0.0), (2.0 * gap, 0.0)], 0.0);
    let r = dark
        .evaluate(Temperature::new(10_000.0), controls(), || false)
        .unwrap();
    assert_eq!(r.coefficient_m3_s, 0.0);
    assert_eq!(r.captured_kinetic_energy_j_m3_s, 0.0);
    assert_eq!(r.emitted_photon_energy_j_m3_s, 0.0);
    assert_eq!(r.mean_captured_energy_j(), None);
    let unknown =
        AtomicRadiativeRecombinationV1::new(dark.source().clone(), Area::new(1e-22), provenance())
            .unwrap();
    assert!(matches!(
        unknown.evaluate(Temperature::new(10_000.0), controls(), || false),
        Err(RadiativeRecombinationErrorV1::UnresolvedTail { .. })
    ));
}

#[test]
fn work_limits_cancellation_and_positive_underflow_refuse_without_result() {
    let source = constant(10_000.0, 0.0, 100.0);
    let mut control = controls();
    control.maximum_subdivisions = 1;
    assert_eq!(
        source.evaluate(Temperature::new(10_000.0), control, || false),
        Err(RadiativeRecombinationErrorV1::QuadratureBudget),
    );
    let mut polls = 0;
    assert_eq!(
        source.evaluate(Temperature::new(10_000.0), controls(), || {
            polls += 1;
            polls == 3
        }),
        Err(RadiativeRecombinationErrorV1::Cancelled),
    );
    assert_eq!(polls, 3);
    let inaccessible = constant(10_000.0, 1000.0, 100.0);
    assert!(matches!(
        inaccessible.evaluate(Temperature::new(10_000.0), controls(), || false),
        Err(RadiativeRecombinationErrorV1::Arithmetic(_)),
    ));
    assert!(matches!(
        source.evaluate(Temperature::new(999.0), controls(), || false),
        Err(RadiativeRecombinationErrorV1::Atomic(_)),
    ));
}

#[test]
fn arbitrarily_hot_authored_calorics_cannot_bypass_nonrelativistic_admission() {
    let models = vec![
        super::source("H", super::MA, 2.0, 0.0, 2.5, 1e9),
        super::source("H+", super::MP, 1.0, super::CHI, 2.5, 1e9),
        super::source("e", super::ELECTRON_MASS_KG, 2.0, 0.0, 2.5, 1e9),
    ];
    let atomic = super::bind(models).unwrap();
    let gap = atomic.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    let photo = AtomicPhotoionizationV1::new(
        atomic,
        provenance(),
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap),
                cross_section: Area::new(1e-22),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap + 100.0 * KB * 1e7),
                cross_section: Area::new(1e-22),
            },
        ],
    )
    .unwrap();
    let source =
        AtomicRadiativeRecombinationV1::new(photo, Area::new(1e-22), provenance()).unwrap();
    assert!(matches!(
        source.evaluate(Temperature::new(1e7), controls(), || false),
        Err(RadiativeRecombinationErrorV1::Invalid(_)),
    ));
}
