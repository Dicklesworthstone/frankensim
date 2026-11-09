//! G1/G3 exact original caloric response versus independent pressure memories.
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
fn one_reservoir_has_the_independent_pressure_driven_time_and_strength() {
    let source = source(101325.0, &[(100.0, 0.002)]);
    let law = source.compressibility_relaxation().unwrap();
    let cp = source.frozen_cv_j_kg_k() + source.gas().specific_gas_constant;
    let expected_tau = 0.002 * cp / (cp + 100.0);
    let expected_strength =
        source.gas().specific_gas_constant * 100.0 / (source.gas().pressure * cp * (cp + 100.0));
    assert_eq!(law.arms().len(), 1);
    assert!((law.arms()[0].relaxation_time_s / expected_tau - 1.0).abs() < 3e-15);
    assert!((law.arms()[0].delta_compressibility_pa_inv / expected_strength - 1.0).abs() < 3e-15);
    assert_eq!(law.source(), &source);
    // Cp relaxation at fixed pressure has a different time than Cv-driven bulk memory.
    let bulk = source.bulk_relaxation().unwrap();
    assert!(
        (law.arms()[0].relaxation_time_s / bulk.arms()[0].relaxation_time_s - 1.0).abs() > 0.01
    );
}
#[test]
fn multiple_original_reservoirs_match_reciprocal_caloric_modulus_over_the_band() {
    let mut worst = 0.0_f64;
    for terms in [
        vec![],
        vec![(0.0, 0.001)],
        vec![(40.0, 0.002), (60.0, 0.002)],
        vec![(10.0, 1e-6), (40.0, 2e-4), (80.0, 3e-2)],
    ] {
        let source = source(85000.0, &terms);
        let law = source.compressibility_relaxation().unwrap();
        if terms.len() == 2 {
            assert_eq!(law.arms().len(), 1);
        }
        let sum: f64 = law
            .arms()
            .iter()
            .map(|a| a.delta_compressibility_pa_inv)
            .sum();
        assert!(
            (law.frozen_compressibility_pa_inv() + sum - law.equilibrium_compressibility_pa_inv())
                .abs()
                < 2e-14 * law.equilibrium_compressibility_pa_inv()
        );
        for omega in [0.0, 1e-8, 0.1, 10.0, 500.0, 1e4, 1e6, 1e12, 1e100, f64::MAX] {
            let actual = law.response(omega).unwrap();
            let expected = source
                .caloric_response(omega)
                .unwrap()
                .bulk_modulus_pa
                .recip();
            let error = (actual - expected).abs() / expected.abs();
            worst = worst.max(error);
            assert!(error < 5e-12, "omega={omega:e}, error={error:e}");
            assert!(actual.im >= 0.0);
        }
    }
    eprintln!("original molecular compressibility maximum relative error {worst:e}");
}
#[test]
fn pressure_scale_preserves_original_reservoir_shape_and_storage_units() {
    let a = source(50000.0, &[(20.0, 0.002), (90.0, 0.03)])
        .compressibility_relaxation()
        .unwrap();
    let b = source(200000.0, &[(20.0, 0.002), (90.0, 0.03)])
        .compressibility_relaxation()
        .unwrap();
    assert_eq!(a.arms().len(), b.arms().len());
    for (x, y) in a.arms().iter().zip(b.arms()) {
        assert_eq!(x.relaxation_time_s, y.relaxation_time_s);
        assert_eq!(
            x.delta_compressibility_pa_inv,
            4.0 * y.delta_compressibility_pa_inv
        );
    }
    for omega in [0.0, 7.0, 1000.0] {
        let x = a.response(omega).unwrap();
        let y = b.response(omega).unwrap();
        assert!((x - y.scale(4.0)).abs() < 1e-18);
    }
}
#[test]
fn original_material_coherence_and_public_arm_scope_are_enforced() {
    let mut stale = source(101325.0, &[(100.0, 0.002)]);
    stale.gas.density *= 1.1;
    assert!(stale.compressibility_relaxation().is_err());
    for (strength, tau) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, 0.0),
        (1.0, f64::INFINITY),
        (f64::NAN, 1.0),
    ] {
        assert!(
            CompressibilityArm {
                delta_compressibility_pa_inv: strength,
                relaxation_time_s: tau
            }
            .response(1.0)
            .is_err()
        );
    }
    let law = source(101325.0, &[]).compressibility_relaxation().unwrap();
    assert!(law.arms().is_empty());
    assert!(law.response(-1.0).is_err());
    assert!(law.response(f64::NAN).is_err());
}
