//! Synthetic NASA-9 coefficients provide an independent caloric oracle.
use super::*;
use crate::thermochemical_gas::tests::{model, point};

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 280.0, maximum_temperature_k: 360.0,
        energy_tolerance_j_kg: 1e-8, maximum_iterations: 64 }
}
fn mole(fractions: Vec<f64>) -> Composition {
    Composition::new(CompositionBasis::MoleFraction, fractions).unwrap()
}

#[test]
fn pressure_enthalpy_inverse_uses_nasa_curvature_and_returns_one_coherent_medium() {
    let source = model(vec![0.75, 0.25], 0.001);
    let law = GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 140000.0), controls()).unwrap();
    let t: f64 = 337.0;
    let molar_mass = 0.75 * 0.028 + 0.25 * 0.004;
    // NASA-9 fixture h/R = x_A*(3.5*T+0.001*T^2/2)+x_B*2.5*T.
    let expected_h = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (3.25 * t + 0.75 * 0.001 * t.powi(2) / 2.0) / molar_mass;
    let solved = law.resolve_enthalpy(expected_h, || false).unwrap();
    assert!((solved.state().state().temperature - t).abs() < 1e-9);
    assert!(solved.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
    assert_eq!(solved.state().state().pressure, 140000.0);
    solved.state().state().validate_coherence().unwrap();
    let expected_density = 140000.0 * molar_mass / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t);
    assert!((solved.state().state().density / expected_density - 1.0).abs() < 1e-11);
    for t in [280.0, 300.0, 337.0, 360.0] {
        let h = law.enthalpy_at_temperature(t).unwrap();
        assert!((law.temperature_k(h).unwrap() - t).abs() < 1e-9);
        let exact_cp = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            * (3.25 + 0.75 * 0.001 * t) / molar_mass;
        assert!(law.maximum_temperature_slope() >= 1.0 / exact_cp);
    }
    // A frozen initial cp would miss the quadratic enthalpy by a finite amount.
    let h300 = law.enthalpy_at_temperature(300.0).unwrap();
    let cp300 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (3.25 + 0.00075 * 300.0) / molar_mass;
    assert!((h300 + cp300 * (t - 300.0) - expected_h).abs() > 100.0);
}

#[test]
fn zero_species_and_shared_reference_survive_caloric_basis_conversion() {
    let source = model(vec![0.75, 0.25], 0.0);
    let a = GasEnthalpyLaw::from_mole(&source, mole(vec![1.0, 0.0]),
        point(300.0, 101325.0), controls()).unwrap();
    let b = GasEnthalpyLaw::from_mole(&source, mole(vec![0.0, 1.0]),
        point(300.0, 101325.0), controls()).unwrap();
    assert_eq!(a.mole_composition().fractions(), &[1.0, 0.0]);
    assert_eq!(a.advected_composition().unwrap().mass_fractions, &[1.0, 0.0]);
    assert_eq!(a.advected_composition().unwrap().species, b.advected_composition().unwrap().species);
    assert_eq!(a.advected_composition().unwrap().enthalpy_reference,
        b.advected_composition().unwrap().enthalpy_reference);
    let h = a.enthalpy_at_temperature(320.0).unwrap();
    assert!(b.temperature_k(h).is_err());
    let mixed = GasEnthalpyLaw::new(&source,
        Composition::new(CompositionBasis::MassFraction, vec![0.25, 0.75]).unwrap(),
        point(300.0, 101325.0), controls()).unwrap();
    let actual = mixed.model().evaluate(&point(320.0, 101325.0)).unwrap();
    let expected_h = 320.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (0.25 * 3.5 / 0.028 + 0.75 * 2.5 / 0.004);
    assert!((actual.specific_enthalpy_j_kg() - expected_h).abs() < 1e-8);
}

#[test]
fn caloric_domain_work_and_cancellation_refuse_without_publishing_a_state() {
    let source = model(vec![0.75, 0.25], 0.001);
    let law = GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 101325.0), controls()).unwrap();
    let h = law.enthalpy_at_temperature(337.0).unwrap();
    assert!(matches!(law.resolve_enthalpy(h, || true), Err(GasEnergyError::Cancelled)));
    assert!(law.resolve_enthalpy(law.enthalpy_at_temperature(360.0).unwrap() + 1.0,
        || false).is_err());
    assert!(GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 101325.0),
        GasEnergySolve { minimum_temperature_k: 199.0, ..controls() }).is_err());
    assert!(GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 101325.0),
        GasEnergySolve { maximum_iterations: 0, ..controls() }).is_err());
    assert!(law.enthalpy_at_temperature(279.0).is_err());
}
