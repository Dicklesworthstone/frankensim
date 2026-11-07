use super::*;
use crate::gas::{GasSpec, GasState};
use crate::relaxation::EvaluatedProcess;
use fs_thermochem::SpeciesId;

fn medium(terms: &[(f64, f64)]) -> RelaxationState {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    RelaxationState {
        gas,
        cv_frozen: gas.specific_heat_cp
            - gas.specific_gas_constant
            - terms.iter().map(|(capacity, _)| capacity).sum::<f64>(),
        processes: terms
            .iter()
            .enumerate()
            .map(|(index, &(capacity, tau))| EvaluatedProcess {
                species: SpeciesId::new("laplace-oracle").unwrap(),
                index,
                heat_capacity_j_kg_k: capacity,
                relaxation_time_s: tau,
            })
            .collect(),
    }
}

fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-30),
        "{a:?} != {b:?}"
    );
}

#[test]
fn single_original_reservoir_matches_independent_elimination_and_propagation() {
    let material = medium(&[(100.0, 0.002)]);
    let gas = material.gas();
    let frozen = material.frozen_cv_j_kg_k();
    let equilibrium = frozen + 100.0;
    for s in [
        C64::ZERO,
        C64::new(3.0, 40.0),
        C64::new(500.0, -800.0),
        C64::new(1e5, 2e4),
    ] {
        let expected = (C64::from_re(equilibrium + gas.specific_gas_constant)
            + s.scale(0.002 * (frozen + gas.specific_gas_constant)))
            / (C64::from_re(equilibrium) + s.scale(0.002 * frozen));
        let actual = material.laplace_response(s).unwrap();
        close(actual.bulk_modulus_pa, expected.scale(gas.pressure), 4e-15);
        close(
            actual.slowness_s_m * actual.slowness_s_m,
            actual.bulk_modulus_pa.recip().scale(gas.density),
            6e-15,
        );
        if s != C64::ZERO {
            close(
                actual.excess_propagation_rate_per_m,
                s * (actual.slowness_s_m - C64::from_re(material.frozen_sound_speed_m_s().recip())),
                2e-12,
            );
        }
        assert!(actual.excess_propagation_rate_per_m.re >= 0.0);
    }
}

#[test]
fn right_half_plane_agrees_with_original_harmonic_and_positive_bulk_realization() {
    let material = medium(&[(40.0, 0.0002), (60.0, 0.003), (20.0, 0.03)]);
    let bulk = material.bulk_relaxation().unwrap();
    for omega in [0.0, 1e-6, 1.0, 200.0, 2e4, 1e10] {
        let boundary = material.laplace_response(C64::new(0.0, -omega)).unwrap();
        let harmonic = material.caloric_response(omega).unwrap();
        close(boundary.cv_j_kg_k, harmonic.cv_j_kg_k, 4e-15);
        close(boundary.bulk_modulus_pa, harmonic.bulk_modulus_pa, 4e-15);
        if omega != 0.0 {
            close(
                boundary.slowness_s_m.scale(omega),
                material.relaxing_wavenumber(omega).unwrap(),
                8e-15,
            );
        }
        for sigma in [1.0, 300.0, 3e4] {
            let s = C64::new(sigma, omega);
            let actual = material.laplace_response(s).unwrap();
            let mut expected = C64::from_re(bulk.equilibrium_bulk_modulus_pa());
            for arm in bulk.arms() {
                let z = s.scale(arm.relaxation_time_s);
                expected = expected + (z / (C64::ONE + z)).scale(arm.delta_bulk_modulus_pa);
            }
            close(actual.bulk_modulus_pa, expected, 3e-13);
            let reflected = material.laplace_response(s.conj()).unwrap();
            assert_eq!(reflected.bulk_modulus_pa, actual.bulk_modulus_pa.conj());
            assert_eq!(
                reflected.excess_propagation_rate_per_m,
                actual.excess_propagation_rate_per_m.conj()
            );
            assert!(actual.excess_propagation_rate_per_m.re >= 0.0);
        }
    }
}

#[test]
fn high_frequency_excess_survives_loss_of_the_slowness_difference() {
    let material = medium(&[(40.0, 0.002)]);
    let s = C64::from_re(1e24);
    let response = material.laplace_response(s).unwrap();
    assert_eq!(
        response.slowness_s_m.re,
        material.frozen_sound_speed_m_s().recip()
    );
    let cv = material.frozen_cv_j_kg_k();
    let r = material.gas().specific_gas_constant;
    let expected = r * 40.0 / (2.0 * material.frozen_sound_speed_m_s() * cv * (cv + r) * 0.002);
    assert!(expected > 0.0);
    close(
        response.excess_propagation_rate_per_m,
        C64::from_re(expected),
        4e-15,
    );
}

#[test]
fn empty_reservoir_limit_extreme_debye_scaling_and_domain_refusals() {
    for terms in [vec![], vec![(0.0, 0.002)], vec![(40.0, 0.002)]] {
        let material = medium(&terms);
        for s in [
            C64::ZERO,
            C64::from_re(f64::MAX),
            C64::new(0.0, f64::MAX),
            C64::new(f64::MAX, f64::MAX),
        ] {
            let value = material.laplace_response(s).unwrap();
            assert!(value.bulk_modulus_pa.re > 0.0);
            assert!(value.excess_propagation_rate_per_m.re >= 0.0);
            if terms.iter().all(|(capacity, _)| *capacity == 0.0) {
                assert_eq!(value.excess_propagation_rate_per_m, C64::ZERO);
                assert_eq!(
                    value.slowness_s_m,
                    C64::from_re(material.frozen_sound_speed_m_s().recip())
                );
            }
        }
        for s in [
            C64::new(-1e-30, 0.0),
            C64::new(f64::NAN, 1.0),
            C64::new(0.0, f64::INFINITY),
        ] {
            assert!(material.laplace_response(s).is_err());
        }
    }
}
