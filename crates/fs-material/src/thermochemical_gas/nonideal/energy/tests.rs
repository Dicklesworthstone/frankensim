//! G0/G1 conservative-state inversion; continuing-vessel tests live in parcel.
use super::*;
use crate::thermochemical_gas::nonideal::PengRobinsonSpeciesBinding;
use crate::thermochemical_gas::tests::{model, point};
use fs_qty::Pressure;
use fs_thermochem::SpeciesId;

fn bound() -> ThermochemicalPengRobinsonModel {
    ThermochemicalPengRobinsonModel::new(model(vec![0.75, 0.25], 0.001), vec![
        PengRobinsonSpeciesBinding { species: SpeciesId::new("A").unwrap(),
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6), acentric_factor: 0.01142 },
        PengRobinsonSpeciesBinding { species: SpeciesId::new("B").unwrap(),
            critical_temperature: Temperature::new(150.687),
            critical_pressure: Pressure::new(4.863e6), acentric_factor: -0.00219 },
    ], vec![0.0, 0.03, 0.03, 0.0]).unwrap()
}
fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 280.0, maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1.0e-7, maximum_iterations: 128 }
}

#[test]
fn g1_absolute_and_reference_energy_recover_one_conserved_density_state() {
    let model = bound();
    let density = 47.889_003_167_633_395;
    let context = point(300.0, 5.0e6).with("lot", 7.0).unwrap();
    let final_point = point(360.0, 5.0e6).with("lot", 7.0).unwrap();
    let expected = model.evaluate_density_temperature(&final_point, density).unwrap();
    let initial = model.evaluate_density_temperature(&context, density).unwrap();
    let absolute = model.reconstruct_density_energy(&context, density,
        expected.specific_internal_energy_j_kg(), controls(), || false).unwrap();
    let relative = model.reconstruct_density_energy_change(&context, density,
        initial.specific_internal_energy_j_kg(),
        expected.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg(),
        controls(), || false).unwrap();
    for solved in [absolute, relative] {
        assert_eq!(solved.state().density_kg_m3().to_bits(), density.to_bits());
        assert!((solved.state().temperature_k() - 360.0).abs() < 2.0e-9);
        assert!((solved.state().pressure_pa() - 6_226_134.941_447_334).abs() < 1.0e-4);
        assert!(solved.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
        assert_eq!(solved.state().source().components()[0].query().axes().get("lot"), Some(&7.0));
    }
}

#[test]
fn g0_invalid_density_unbracketed_energy_and_cancellation_do_not_return_states() {
    let model = bound();
    let context = point(300.0, 5.0e6);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.reconstruct_density_energy(&context, bad, 100_000.0, controls(), || false).is_err());
    }
    assert!(matches!(model.reconstruct_density_energy(&context, 50.0, 1.0e12,
        controls(), || false), Err(GasEnergyError::OutsideBracket { .. })));
    assert!(matches!(model.reconstruct_density_energy(&context, 50.0, 100_000.0,
        controls(), || true), Err(GasEnergyError::Cancelled)));
    assert!(model.reconstruct_density_internal_energy_increment(&point(200.0, 5.0e6),
        50.0, 1.0, controls(), || false).is_err());
    let zero = model.reconstruct_density_internal_energy_increment(&context, 50.0,
        0.0, controls(), || false).unwrap();
    assert_eq!(zero.state().temperature_k(), 300.0);
    assert_eq!(zero.state().density_kg_m3(), 50.0);
    assert_eq!(zero.iterations(), 0);
}

#[test]
fn g1_joule_expansion_retains_original_density_departure_and_matches_oracle() {
    let model = bound();
    let rho0 = 47.889_003_167_633_395;
    let rho1 = rho0 / 1.25;
    let context = point(300.0, 5.0e6);
    // Independent 75-digit Decimal PR energy balance, retaining the original
    // residual energy at rho0: Delta u_ideal + uR(rho1,T1)-uR(rho0,T0)=0.
    // This is a constitutive endpoint; the parcel owns free-expansion admission.
    let solved = model.reconstruct_density_internal_energy_increment_between_densities(
        &context, rho0, rho1, 0.0, controls(), || false).unwrap();
    let state = solved.state();
    assert_eq!(state.density_kg_m3().to_bits(), rho1.to_bits());
    assert!((state.temperature_k() - 294.872_058_857_231_76).abs() < 2.0e-9);
    assert!((state.pressure_pa() - 3_975_006.042_948_880_2).abs() < 1.0e-4);
    assert!(solved.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
    let entropy = model.specific_entropy_change_between_densities(rho0, 300.0,
        rho1, state.temperature_k()).unwrap();
    assert!((entropy - 78.087_722_432_685_2).abs() < 1.0e-8);
    // Fixed-u expansion is not fixed-T for this nonideal material.
    assert!((state.temperature_k() - 300.0).abs() > 5.0);
    let returned = model.reconstruct_density_internal_energy_increment_between_densities(
        &point(state.temperature_k(), state.pressure_pa()), rho1, rho0,
        0.0, controls(), || false).unwrap();
    assert!((returned.state().temperature_k() - 300.0).abs() < 2.0e-9);
}

#[test]
fn g1_general_energy_increment_keeps_both_density_endpoints_and_heat_direction() {
    let model = bound();
    let (rho0, rho1) = (40.0, 60.0);
    let context = point(300.0, 5.0e6);
    let isothermal_change = model.specific_internal_energy_change_between_densities(
        rho0, 300.0, rho1, 300.0).unwrap();
    assert!(isothermal_change.abs() > 1000.0);
    let isothermal = model.reconstruct_density_internal_energy_increment_between_densities(
        &context, rho0, rho1, isothermal_change, controls(), || false).unwrap();
    assert!((isothermal.state().temperature_k() - 300.0).abs() < 2.0e-9);
    let target = model.specific_internal_energy_change_between_densities(
        rho0, 300.0, rho1, 350.0).unwrap();
    let returned = model.specific_internal_energy_change_between_densities(
        rho1, 350.0, rho0, 300.0).unwrap();
    assert_eq!(returned, -target);
    let solved = model.reconstruct_density_internal_energy_increment_between_densities(
        &context, rho0, rho1, target, controls(), || false).unwrap();
    assert!((solved.state().temperature_k() - 350.0).abs() < 2.0e-9);
    // Existing fixed-density interface delegates to the same conserved law.
    assert_eq!(model.specific_internal_energy_change_at_density(rho0, 300.0, 350.0).unwrap(),
        model.specific_internal_energy_change_between_densities(rho0, 300.0, rho0, 350.0).unwrap());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.reconstruct_density_internal_energy_increment_between_densities(
            &context, bad, rho1, target, controls(), || false).is_err());
        assert!(model.reconstruct_density_internal_energy_increment_between_densities(
            &context, rho0, bad, target, controls(), || false).is_err());
    }
}
