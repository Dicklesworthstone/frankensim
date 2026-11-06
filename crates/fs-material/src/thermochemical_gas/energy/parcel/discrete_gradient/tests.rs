//! G1/G3/G4: source calorics, conjugate work, entropy and atomic proposals.
use super::*;
use crate::thermochemical_gas::energy::entropy::GasEntropySolve;
use crate::thermochemical_gas::tests::{
    model, point, species_with_energy_reference, transport,
};
use fs_thermochem::{Composition, CompositionBasis};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn parcel(slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![0.75, 0.25], slope),
        point(400.0, 100000.0).with("material-lot", 42.0).unwrap(), 0.001).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance * b.abs().max(1.0), "{a} != {b}");
}
fn entropy_change(old: &ClosedGasParcel, end: &ClosedGasParcel) -> f64 {
    old.mass_kg() * old.model().specific_entropy_change_j_kg_k(
        old.point(), end.point()).unwrap()
}

#[test]
fn g1_constant_capacity_adiabatic_motion_matches_the_exact_isentrope() {
    let old = parcel(0.0);
    let gas = old.state().state();
    let cv = old.state().thermodynamics().properties().mass_specific_cv().value();
    for volume in [0.00055, 0.0007, 0.001, 0.0014, 0.002] {
        let result = old.discrete_gradient_heat_volume(volume, 0.0, controls(), || false).unwrap();
        let expected = gas.temperature * (old.volume_m3()/volume).powf(gas.specific_gas_constant/cv);
        near(result.endpoint().state().state().temperature, expected, 2e-11);
        near(result.step().stored_energy_change_j, old.mass_kg()*cv*(expected-gas.temperature), 2e-9);
        near(result.step().work_by_gas_j, -result.step().stored_energy_change_j, 2e-9);
        assert!(result.step().entropy_residual_j_k.abs()
            <= old.mass_kg()*controls().energy_tolerance_j_kg/result.thermodynamic_mean().state().temperature);
        assert!(entropy_change(&old, result.endpoint()).abs() < 1e-10);
        assert_eq!(result.endpoint().mass_kg(), old.mass_kg());
    }
}

#[test]
fn g1_nonlinear_capacity_adiabatic_endpoint_matches_original_entropy_inversion() {
    let old = parcel(0.002);
    for volume in [0.0006, 0.0009, 0.0015] {
        let result = old.discrete_gradient_heat_volume(volume, 0.0, controls(), || false).unwrap();
        let mut reference = old.clone();
        reference.set_volume_isentropic(volume, GasEntropySolve {
            minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
            entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128,
        }, || false).unwrap();
        near(result.endpoint().state().state().temperature,
            reference.state().state().temperature, 2e-10);
        near(result.endpoint().state().state().pressure,
            reference.state().state().pressure, 2e-10);
        near(result.step().work_by_gas_j, result.thermodynamic_mean().state().pressure
            *(volume-old.volume_m3()), 1e-14);
        assert_eq!(result.endpoint().point().axes().get("material-lot"), Some(&42.0));
    }
}

#[test]
fn g1_combined_heat_and_work_follow_the_constant_capacity_discrete_gibbs_equation() {
    let old = parcel(0.0);
    let r = old.state().state().specific_gas_constant;
    let cv = old.state().thermodynamics().properties().mass_specific_cv().value();
    let final_temperature = 470.0;
    let theta = (final_temperature-400.0)/(final_temperature/400.0_f64).ln();
    for volume in [0.0007, 0.001, 0.0015] {
        let heat = old.mass_kg()*(cv*(final_temperature-400.0)
            +r*theta*(volume/old.volume_m3()).ln());
        let result = old.discrete_gradient_heat_volume(volume, heat, controls(), || false).unwrap();
        near(result.endpoint().state().state().temperature, final_temperature, 2e-11);
        near(result.thermodynamic_mean().state().temperature, theta, 2e-11);
        near(result.step().entropy_change_j_k, heat/theta, 2e-10);
        assert!(result.step().energy_residual_j.abs() <= old.mass_kg()*controls().energy_tolerance_j_kg);
        assert!(result.step().gibbs_residual_j.abs() <= old.mass_kg()*controls().energy_tolerance_j_kg);
    }
}

#[test]
fn g1_fixed_volume_heat_uses_original_storage_and_conjugate_temperature() {
    let old = parcel(0.002);
    let result = old.discrete_gradient_heat_volume(old.volume_m3(), 50.0, controls(), || false).unwrap();
    let mut reference = old.clone();
    reference.add_heat(50.0, controls(), || false).unwrap();
    near(result.endpoint().state().state().temperature, reference.state().state().temperature, 2e-11);
    assert_eq!(result.step().work_by_gas_j, 0.0);
    let theta = result.thermodynamic_mean().state().temperature;
    let arithmetic = 0.5*(old.state().state().temperature+result.endpoint().state().state().temperature);
    assert!((theta-arithmetic).abs() > 1e-4);
    near(entropy_change(&old, result.endpoint()), 50.0/theta, 1e-11);
    assert!(result.endpoint().state().state().sound_speed > old.state().state().sound_speed);
    assert!(result.endpoint().state().state().thermal_conductivity > old.state().state().thermal_conductivity);
}

#[test]
fn g1_implicit_pressure_sensitivities_match_independent_endpoint_perturbations() {
    for slope in [0.0, 0.002] {
        let old = parcel(slope);
        for (volume, heat) in [(0.0007, 20.0), (0.001, 0.0), (0.001, 30.0), (0.0015, -10.0)] {
            let result = old.discrete_gradient_heat_volume(volume, heat, controls(), || false).unwrap();
            let pressure = |v, q| old.discrete_gradient_heat_volume(v, q, controls(), || false)
                .unwrap().thermodynamic_mean().state().pressure;
            let dv = 1e-8;
            let dq = 0.001;
            near(result.step().pressure_final_volume_derivative_pa_m3,
                (pressure(volume+dv, heat)-pressure(volume-dv, heat))/(2.0*dv), 2e-7);
            near(result.step().pressure_heat_derivative_pa_j,
                (pressure(volume, heat+dq)-pressure(volume, heat-dq))/(2.0*dq), 2e-7);
        }
    }
}

#[test]
fn g1_close_temperature_mean_derivative_retains_the_logarithmic_limit() {
    let old = parcel(0.0);
    for delta in [-1e-3, -1e-6, -1e-9, 1e-9, 1e-6, 1e-3] {
        let t = 400.0+delta;
        let c = caloric_mean(old.model(), 400.0, t, &mut || false).unwrap();
        let expected_derivative = 0.5-(t-400.0)/(6.0*400.0);
        near(c.derivative, expected_derivative, 2e-12);
        near(c.temperature, 400.0+0.5*(t-400.0), 1e-11);
    }
}

#[test]
fn g3_reverse_heat_volume_steps_restore_the_initial_thermodynamic_state() {
    let old = parcel(0.002);
    for (volume, heat) in [(0.0007, 20.0), (0.0015, -10.0)] {
        let forward = old.discrete_gradient_heat_volume(volume, heat, controls(), || false).unwrap();
        let reverse = forward.endpoint().discrete_gradient_heat_volume(
            old.volume_m3(), -heat, controls(), || false).unwrap();
        near(reverse.endpoint().state().state().temperature, old.state().state().temperature, 3e-10);
        near(reverse.endpoint().state().state().pressure, old.state().state().pressure, 3e-10);
        near(reverse.thermodynamic_mean().state().pressure, forward.thermodynamic_mean().state().pressure, 3e-10);
        near(reverse.step().work_by_gas_j, -forward.step().work_by_gas_j, 3e-10);
        near(reverse.step().entropy_change_j_k, -forward.step().entropy_change_j_k, 2e-10);
    }
}

fn offset_parcel(reference: f64) -> ClosedGasParcel {
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A", 0.028, 3.5, 0.002, reference),
            species_with_energy_reference("B", 0.004, 2.5, 0.0, -reference)],
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    ).unwrap();
    let model = ThermochemicalGasModel::new(mixture,
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap();
    ClosedGasParcel::new(model, point(400.0, 100000.0), 0.001).unwrap()
}

#[test]
fn g3_formation_references_do_not_change_work_entropy_or_discrete_tangents() {
    let base = offset_parcel(0.0);
    for (volume, heat) in [(0.0007, 20.0), (0.001, 1e-8), (0.0015, 0.0)] {
        let expected = base.discrete_gradient_heat_volume(volume, heat, controls(), || false).unwrap();
        for reference in [-1e24, 1e24] {
            let old = offset_parcel(reference);
            let result = old.discrete_gradient_heat_volume(volume, heat, controls(), || false).unwrap();
            assert_eq!(result.step(), expected.step());
            assert_eq!(result.thermodynamic_mean().state(), expected.thermodynamic_mean().state());
            assert_eq!(result.endpoint().state().state(), expected.endpoint().state().state());
            assert!(result.endpoint().state().specific_internal_energy_j_kg().abs() > 1e20);
        }
    }
}

#[test]
fn g1_tiny_heating_and_isothermal_work_do_not_disappear() {
    let old = parcel(0.002);
    let heat = 1e-8;
    let heated = old.discrete_gradient_heat_volume(old.volume_m3(), heat, controls(), || false).unwrap();
    assert!(heated.endpoint().state().state().temperature > old.state().state().temperature);
    assert!(heated.step().stored_energy_change_j > 0.0);
    near(heated.step().stored_energy_change_j, heat, 1e-10);
    for volume in [0.0008, 0.0014, old.volume_m3()*(1.0+1e-10)] {
        let heat = old.mass_kg()*old.state().state().specific_gas_constant*400.0
            * log_ratio(volume, old.volume_m3());
        let result = old.discrete_gradient_heat_volume(volume, heat, controls(), || false).unwrap();
        near(result.endpoint().state().state().temperature, 400.0, 1e-11);
        near(result.step().work_by_gas_j, heat, 1e-10);
        assert_ne!(result.step().work_by_gas_j, 0.0);
    }
}

#[test]
fn g4_noop_and_refusals_preserve_the_original_parcel() {
    let old = parcel(0.002);
    let saved = old.clone();
    let result = old.discrete_gradient_heat_volume(old.volume_m3(), 0.0, controls(), || false).unwrap();
    assert_eq!(result.endpoint(), &old);
    assert_eq!(result.thermodynamic_mean(), old.state());
    assert_eq!(result.step().entropy_change_j_k, 0.0);
    assert_eq!(result.step().energy_residual_j, 0.0);
    for (volume, heat) in [(0.0, 0.0), (0.00001, 0.0), (0.001, 1e9), (0.001, f64::NAN)] {
        assert!(old.discrete_gradient_heat_volume(volume, heat, controls(), || false).is_err());
        assert_eq!(old, saved);
    }
    let mut calls = 0;
    old.discrete_gradient_heat_volume(0.001, 1e-5, controls(), || {
        calls += 1; false
    }).unwrap();
    for cancelled_at in [1, calls/2, calls] {
        let mut current = 0;
        assert!(matches!(old.discrete_gradient_heat_volume(0.001, 1e-5, controls(), || {
            current += 1;
            current == cancelled_at
        }), Err(GasEnergyError::Cancelled)));
        assert_eq!(old, saved);
    }
}
