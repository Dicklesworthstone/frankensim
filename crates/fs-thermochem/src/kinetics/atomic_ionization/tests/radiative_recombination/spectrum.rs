use super::{constant, controls, table};
use crate::electron_collision::{
    AVOGADRO_PER_MOL as NA, BOLTZMANN_J_PER_K as KB, ELECTRON_MASS_KG as ME,
};
use crate::kinetics::atomic_ionization::radiative_recombination::RadiativeRecombinationErrorV1;
use fs_qty::{Energy, Temperature};

fn relative(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual / expected - 1.0).abs() < allowance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_spectral_density_and_envelope_match_independent_constant_source_integrals() {
    let source = constant(10_000.0, 0.0, 100.0);
    let spectrum = source.spectrum(Temperature::new(10_000.0)).unwrap();
    let theta = KB * 10_000.0;
    let chi = source.source().channel().ionization_energy_j_mol() / NA;
    let factor =
        2.0 * (2.0 / std::f64::consts::PI).sqrt() / (ME.powf(1.5) * 299_792_458.0_f64.powi(2));
    for x in [0.0, 0.25, 1.0, 4.0, 20.0, 80.0] {
        let energy = chi + x * theta;
        let kinetic = energy - chi;
        let actual = spectrum.evaluate(Energy::new(energy)).unwrap();
        let density = factor * 1e-22 * energy.powi(2) / theta.powf(1.5) * (-kinetic / theta).exp();
        relative(actual.coefficient_m3_s_per_j, density, 5e-14);
        if kinetic == 0.0 {
            assert_eq!(actual.captured_kinetic_energy_m3_s, 0.0);
        } else {
            relative(
                actual.captured_kinetic_energy_m3_s,
                kinetic * density,
                5e-14,
            );
        }
        let bound = spectrum.tail_bound(Energy::new(energy)).unwrap();
        let polynomial = energy * energy + 2.0 * energy * theta + 2.0 * theta * theta;
        let event = factor * 1e-22 / theta.sqrt() * (-kinetic / theta).exp() * polynomial;
        let captured = factor * 1e-22 / theta.sqrt()
            * (-kinetic / theta).exp()
            * (kinetic * polynomial
                + theta * (energy * energy + 4.0 * energy * theta + 6.0 * theta * theta));
        relative(bound.coefficient_m3_s, event, 5e-12);
        relative(bound.captured_kinetic_energy_j_m3_s, captured, 5e-12);
        assert!(bound.coefficient_m3_s >= event * (1.0 - 8e-15));
        assert!(bound.captured_kinetic_energy_j_m3_s >= captured * (1.0 - 8e-15));
    }
    let complete = source
        .evaluate(Temperature::new(10_000.0), controls(), || false)
        .unwrap();
    let bound = spectrum.tail_bound(Energy::new(0.0)).unwrap();
    relative(bound.coefficient_m3_s, complete.coefficient_m3_s, 3e-12);
    relative(
        bound.captured_kinetic_energy_j_m3_s,
        complete.captured_kinetic_energy_j_m3_s,
        3e-12,
    );
}

#[test]
fn g1_original_onset_and_differential_quadrature_recover_shifted_milne_moments() {
    let source = constant(10_000.0, 5.0, 100.0);
    let spectrum = source.spectrum(Temperature::new(10_000.0)).unwrap();
    let onset = source.source().threshold().value();
    let theta = KB * 10_000.0;
    assert_eq!(
        spectrum
            .evaluate(Energy::new(onset - theta))
            .unwrap()
            .coefficient_m3_s_per_j,
        0.0
    );
    // Independent composite Simpson integration, with a second grid to check
    // its own convergence; no production Gauss or tail routine is reused.
    let mut previous = [0.0; 2];
    for n in [4096, 8192] {
        let h = 60.0 * theta / n as f64;
        let mut sum = [0.0; 2];
        for i in 0..=n {
            let value = spectrum
                .evaluate(Energy::new(onset + i as f64 * h))
                .unwrap();
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum[0] += weight * value.coefficient_m3_s_per_j;
            sum[1] += weight * value.captured_kinetic_energy_m3_s;
        }
        let current = sum.map(|x| x * h / 3.0);
        if n == 8192 {
            for j in 0..2 {
                relative(current[j], previous[j], 2e-9);
            }
            let total = source
                .evaluate(Temperature::new(10_000.0), controls(), || false)
                .unwrap();
            relative(current[0], total.coefficient_m3_s, 2e-10);
            relative(current[1], total.captured_kinetic_energy_j_m3_s, 2e-10);
        }
        previous = current;
    }
}

#[test]
fn g0_spectral_source_domains_zero_laws_and_tail_underflow_remain_explicit() {
    let source = constant(10_000.0, 0.0, 100.0);
    let spectrum = source.spectrum(Temperature::new(10_000.0)).unwrap();
    assert!(matches!(
        spectrum.evaluate(Energy::new(source.source().maximum_energy().value() * 2.0)),
        Err(RadiativeRecombinationErrorV1::Photoionization(_))
    ));
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(spectrum.tail_bound(Energy::new(value)).is_err());
        assert!(source.spectrum(Temperature::new(value)).is_err());
    }
    let tiny = spectrum
        .tail_bound(Energy::new(1e6 * spectrum.thermal_energy().value()))
        .unwrap();
    assert!(tiny.coefficient_m3_s > 0.0);
    assert!(tiny.captured_kinetic_energy_j_m3_s > 0.0);
    let onset = source.source().threshold().value();
    let zero = table(&[(onset, 0.0), (onset * 2.0, 0.0)], 0.0);
    let dark = zero.spectrum(Temperature::new(10_000.0)).unwrap();
    assert_eq!(
        dark.evaluate(Energy::new(onset))
            .unwrap()
            .coefficient_m3_s_per_j,
        0.0
    );
    assert_eq!(
        dark.tail_bound(Energy::new(onset))
            .unwrap()
            .emitted_photon_energy_j_m3_s,
        0.0
    );
}
