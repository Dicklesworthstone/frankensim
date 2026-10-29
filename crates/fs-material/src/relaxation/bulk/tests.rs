use super::*;
use crate::gas::{GasSpec, GasState};
use crate::relaxation::EvaluatedProcess;
use fs_thermochem::SpeciesId;

fn source(pressure: f64, terms: &[(f64, f64)]) -> RelaxationState {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    RelaxationState {
        gas,
        cv_frozen: gas.specific_heat_cp
            - gas.specific_gas_constant
            - terms.iter().map(|(c, _)| c).sum::<f64>(),
        processes: terms
            .iter()
            .enumerate()
            .map(|(index, &(capacity, tau))| EvaluatedProcess {
                species: SpeciesId::new("fixture").unwrap(),
                index,
                heat_capacity_j_kg_k: capacity,
                relaxation_time_s: tau,
            })
            .collect(),
    }
}

#[test]
fn single_reservoir_bulk_memory_matches_analytic_elimination() {
    // G1: independently eliminate one internal temperature, with no eigensolve.
    let source = source(101325.0, &[(100.0, 0.002)]);
    let law = source.bulk_relaxation().unwrap();
    let cv = source.frozen_cv_j_kg_k();
    let gas = source.gas();
    let expected_tau = 0.002 * cv / (cv + 100.0);
    let expected_delta = gas.pressure * gas.specific_gas_constant * 100.0 / (cv * (cv + 100.0));
    assert_eq!(law.arms().len(), 1);
    assert!((law.arms()[0].relaxation_time_s / expected_tau - 1.0).abs() < 2e-15);
    assert!((law.arms()[0].delta_bulk_modulus_pa / expected_delta - 1.0).abs() < 2e-15);
    assert_eq!(law.source(), &source);
    for omega in [0.0, 1e-8, 1.0, 500.0, 1e5, 1e20, f64::MAX] {
        let actual = law.response(omega).unwrap();
        let exact = source.caloric_response(omega).unwrap().bulk_modulus_pa;
        assert!((actual - exact).abs() < 1e-14 * exact.abs());
        assert!(actual.im <= 0.0);
    }
}

#[test]
fn coupled_reservoir_spectrum_preserves_the_original_caloric_law() {
    // G1/G3: disparate/equal molecular times, zero strengths and no reservoirs.
    for terms in [
        vec![],
        vec![(0.0, 0.001)],
        vec![(40.0, 0.002), (60.0, 0.002)],
        vec![(10.0, 1e-6), (40.0, 2e-4), (80.0, 3e-2)],
    ] {
        let source = source(85000.0, &terms);
        let law = source.bulk_relaxation().unwrap();
        if terms.len() == 2 {
            assert_eq!(law.arms().len(), 1);
        }
        assert!(
            law.arms()
                .iter()
                .all(|a| a.delta_bulk_modulus_pa > 0.0 && a.relaxation_time_s > 0.0)
        );
        let frozen = law.equilibrium_bulk_modulus_pa()
            + law
                .arms()
                .iter()
                .map(|a| a.delta_bulk_modulus_pa)
                .sum::<f64>();
        assert!((frozen / law.frozen_bulk_modulus_pa() - 1.0).abs() < 1e-12);
        for i in 0..=80 {
            let omega = det::pow(10.0, -3.0 + 0.15 * i as f64);
            let exact = source.caloric_response(omega).unwrap().bulk_modulus_pa;
            let actual = law.response(omega).unwrap();
            assert!(
                (actual - exact).abs() < 2e-12 * exact.abs(),
                "{omega}: {actual:?} != {exact:?}"
            );
        }
        for omega in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(law.response(omega).is_err());
        }
    }
}

#[test]
fn source_pressure_and_collision_time_change_the_same_bulk_memory() {
    // G3: a pressure-independent p*tau gives proportional modulus and rates.
    let first = source(50000.0, &[(80.0, 0.004)]).bulk_relaxation().unwrap();
    let second = source(100000.0, &[(80.0, 0.002)])
        .bulk_relaxation()
        .unwrap();
    assert_eq!(first.arms().len(), second.arms().len());
    for (a, b) in first.arms().iter().zip(second.arms()) {
        assert_eq!(2.0 * a.delta_bulk_modulus_pa, b.delta_bulk_modulus_pa);
        assert_eq!(a.relaxation_time_s, 2.0 * b.relaxation_time_s);
    }
    let a = first.response(500.0).unwrap();
    let b = second.response(1000.0).unwrap();
    assert!((b - a.scale(2.0)).abs() < 1e-10);
}
