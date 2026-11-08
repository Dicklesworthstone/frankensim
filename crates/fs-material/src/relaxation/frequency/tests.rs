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
                species: SpeciesId::new("frequency-oracle").unwrap(),
                index,
                heat_capacity_j_kg_k: capacity,
                relaxation_time_s: tau,
            })
            .collect(),
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn frequency_values_retain_original_bits_and_empty_reservoirs_recover_constant_slowness() {
    for terms in [vec![], vec![(0.0, 0.002)], vec![(80.0, 0.001), (30.0, 0.02)]] {
        let material = medium(&terms);
        for omega in [1e-8, 1.0, 300.0, 3000.0, 1e8] {
            let response = material.relaxing_frequency_response(omega).unwrap();
            assert_eq!(response.wavenumber_per_m, material.relaxing_wavenumber(omega).unwrap());
            if terms.iter().all(|(capacity, _)| *capacity == 0.0) {
                let slope = material.relaxing_wavenumber(1.0).unwrap();
                assert_eq!(response.wavenumber_derivative_s_m, slope);
                assert_eq!(response.wavenumber_derivative_s_m.im, 0.0);
            }
        }
    }
}

#[test]
fn original_reservoir_derivative_matches_independent_rational_elimination() {
    let (capacity, tau) = (100.0, 0.002);
    let material = medium(&[(capacity, tau)]);
    let gas = material.gas();
    let frozen = material.frozen_cv_j_kg_k();
    let r = gas.specific_gas_constant;
    for omega in [0.01, 10.0, 300.0, 1000.0, 3e4] {
        // Eliminate the one reservoir first: slowness^2=(rho/p)*N/D.
        let n = C64::new(frozen + capacity, -omega * tau * frozen);
        let d = C64::new(frozen + capacity + r, -omega * tau * (frozen + r));
        let slowness = (n / d).scale(gas.density / gas.pressure).sqrt();
        let expected = slowness
            * (C64::ONE + C64::new(0.0, 0.5 * omega * tau * r * capacity) / (n * d));
        let actual = material.relaxing_frequency_response(omega).unwrap();
        close(actual.wavenumber_derivative_s_m, expected, 5e-15);
    }
}

#[test]
fn multiple_reservoir_slopes_match_five_point_original_values() {
    let material = medium(&[(80.0, 0.001), (30.0, 0.02)]);
    for omega in [20.0, 300.0, 2000.0, 2e4] {
        let h = omega * 0.001;
        let value = |offset: f64| material.relaxing_wavenumber(omega + offset * h).unwrap();
        let expected = (value(-2.0) - value(-1.0).scale(8.0)
            + value(1.0).scale(8.0) - value(2.0)).scale(1.0 / (12.0 * h));
        let actual = material.relaxing_frequency_response(omega).unwrap().wavenumber_derivative_s_m;
        close(actual, expected, 2e-10);
        assert!((actual.im - expected.im).abs() < 2e-8 * expected.im.abs());
    }
    let response = material.relaxing_frequency_response(600.0).unwrap();
    assert!((response.wavenumber_derivative_s_m.re - 1.0 / material.gas().sound_speed).abs() > 1e-5);
    assert!((response.wavenumber_derivative_s_m.re - response.wavenumber_per_m.re / 600.0).abs() > 1e-5);
}

#[test]
fn frozen_tail_retains_the_small_attenuation_slope_without_cancellation() {
    let terms = [(80.0, 0.001), (30.0, 0.02)];
    let material = medium(&terms);
    let frozen = material.frozen_cv_j_kg_k();
    let r = material.gas().specific_gas_constant;
    let cp_inf = frozen + r;
    let a1: f64 = terms.iter().map(|(c, t)| c / t).sum();
    let a2: f64 = terms.iter().map(|(c, t)| c / (t * t)).sum();
    let a3: f64 = terms.iter().map(|(c, t)| c / (t * t * t)).sum();
    let scale = r / (frozen * cp_inf);
    let b1 = scale * a1;
    let b2 = scale * (a2 + a1 * a1 / cp_inf);
    let b3 = scale * (a3 + 2.0 * a1 * a2 / cp_inf + a1 * a1 * a1 / (cp_inf * cp_inf));
    // Independent large-frequency expansion of sqrt(rho*cv/(p*cp)).
    let coefficient = (b3 + 0.5 * b1 * b2 + b1 * b1 * b1 / 8.0)
        / material.frozen_sound_speed_m_s();
    for omega in [1e12, 1e16, 1e50] {
        let actual = material.relaxing_frequency_response(omega).unwrap();
        let expected = coefficient / (omega * omega * omega);
        assert!(actual.wavenumber_derivative_s_m.im > 0.0);
        assert!((actual.wavenumber_derivative_s_m.im / expected - 1.0).abs() < 2e-13);
    }
}

#[test]
fn extreme_debye_scales_remain_finite_and_invalid_frequencies_refuse() {
    let material = medium(&[(80.0, 0.001)]);
    for omega in [1e-250, 1e250] {
        let actual = material.relaxing_frequency_response(omega).unwrap();
        assert!(actual.wavenumber_per_m.re > 0.0);
        assert!(actual.wavenumber_derivative_s_m.re > 0.0);
    }
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(material.relaxing_frequency_response(omega).is_err());
    }
}
