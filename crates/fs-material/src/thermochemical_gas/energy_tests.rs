//! Numerical closure tests against independent synthetic caloric formulas.
use super::*;
use super::super::tests::{close, model, point};
use crate::fluid_transport::FluidTransport;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn specific_u(t: f64, a_fraction: f64, slope: f64) -> f64 {
    let mass = a_fraction * 0.028 + (1.0 - a_fraction) * 0.004;
    let a = a_fraction * 3.5 + (1.0 - a_fraction) * 2.5;
    let b = a_fraction * slope;
    (R / mass) * ((a - 1.0) * t + 0.5 * b * t * t)
}

#[test]
fn fixed_volume_heat_changes_pressure_and_the_complete_shared_medium() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.0), point(300.0, 100000.0), 0.01).unwrap();
    let before = parcel.state().clone();
    let mass = parcel.mass_kg();
    let q = mass * (specific_u(600.0, 0.75, 0.0) - specific_u(300.0, 0.75, 0.0));
    let report = parcel.add_heat(q, controls(), || false).unwrap();
    close(parcel.state().state().temperature, 600.0);
    close(parcel.state().state().pressure, 200000.0);
    close(parcel.state().state().density, before.state().density);
    close(parcel.state().state().sound_speed / before.state().sound_speed, 2.0f64.sqrt());
    assert!(parcel.state().state().dynamic_viscosity > before.state().dynamic_viscosity);
    assert!(parcel.state().transport_state().unwrap().thermal_diffusivity_m2_s() >
        before.transport_state().unwrap().thermal_diffusivity_m2_s());
    assert_eq!(parcel.mass_kg(), mass);
    assert_eq!(parcel.volume_m3(), 0.01);
    assert!(report.residual_j.abs() <= mass * controls().energy_tolerance_j_kg);
    close(report.stored_energy_change_j, q);
    close(*parcel.point().axes().get("temperature").unwrap(), 600.0);
}

#[test]
fn nonlinear_caloric_source_changes_heat_response_without_a_new_solver() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 100000.0), 0.01).unwrap();
    let initial = parcel.state().clone();
    let du = specific_u(700.0, 0.75, 0.001) - specific_u(300.0, 0.75, 0.001);
    let q = parcel.mass_kg() * du;
    let cv0 = initial.thermodynamics().properties().mass_specific_cv().value();
    assert!((300.0 + du / cv0 - 700.0).abs() > 10.0);
    let result = parcel.add_heat(q, controls(), || false).unwrap();
    close(parcel.state().state().temperature, 700.0);
    close(parcel.state().state().pressure, 100000.0 * 700.0 / 300.0);
    assert!(parcel.state().state().gamma < initial.state().gamma);
    assert!(result.iterations > 1);
    let inverse = parcel.add_heat(-q, controls(), || false).unwrap();
    close(parcel.state().state().temperature, 300.0);
    assert!(inverse.residual_j.abs() <= parcel.mass_kg() * controls().energy_tolerance_j_kg);
}

#[test]
fn conserved_density_energy_reconstruction_roundtrips_and_preserves_other_conditions() {
    let model = model(vec![0.25, 0.75], 0.001);
    let context = point(300.0, 100000.0).with("material-lot", 12.0).unwrap();
    let density = model.evaluate(&context).unwrap().state().density;
    let result = model.reconstruct_density_energy(&context, density,
        specific_u(725.0, 0.25, 0.001), controls(), || false).unwrap();
    close(result.state().state().temperature, 725.0);
    close(result.state().state().density, density);
    let mut parcel = IsochoricGasParcel::new(model, context, 0.01).unwrap();
    parcel.add_heat(100.0, controls(), || false).unwrap();
    assert_eq!(parcel.point().axes().get("material-lot"), Some(&12.0));
}

#[test]
fn failure_cancellation_and_work_exhaustion_are_atomic() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 100000.0), 0.01).unwrap();
    let saved = parcel.clone();
    assert!(matches!(parcel.add_heat(10.0, controls(), || true), Err(GasEnergyError::Cancelled)));
    assert_eq!(parcel, saved);
    assert!(parcel.add_heat(f64::NAN, controls(), || false).is_err());
    assert_eq!(parcel, saved);
    assert!(matches!(parcel.add_heat(1e9, controls(), || false), Err(GasEnergyError::OutsideBracket { .. })));
    assert_eq!(parcel, saved);
    let limited = GasEnergySolve { maximum_iterations: 1, ..controls() };
    assert!(matches!(parcel.add_heat(100.0, limited, || false), Err(GasEnergyError::NotConverged { .. })));
    assert_eq!(parcel, saved);
    let mut polls = 0;
    assert!(matches!(parcel.add_heat(100.0, controls(), || { polls += 1; polls == 4 }), Err(GasEnergyError::Cancelled)));
    assert_eq!(parcel, saved);
    // Target is inside the thermodynamic bracket but pressure crosses the
    // shared gas admission limit. No tentative state may become observable.
    let mut dense = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 8e6), 0.01).unwrap();
    let saved_dense = dense.clone();
    let q = dense.mass_kg() * (specific_u(600.0, 0.75, 0.001) - specific_u(300.0, 0.75, 0.001));
    assert!(dense.add_heat(q, controls(), || false).is_err());
    assert_eq!(dense, saved_dense);
}

#[test]
fn zero_heat_is_bitwise_unchanged_and_bad_controls_are_not_ignored() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 100000.0), 0.01).unwrap();
    let saved = parcel.clone();
    let step = parcel.add_heat(0.0, controls(), || false).unwrap();
    assert_eq!(parcel, saved);
    assert_eq!(step.stored_energy_change_j, 0.0);
    for invalid in [GasEnergySolve { maximum_iterations: 0, ..controls() },
        GasEnergySolve { energy_tolerance_j_kg: 0.0, ..controls() },
        GasEnergySolve { minimum_temperature_k: 1600.0, ..controls() }] {
        assert!(parcel.add_heat(0.0, invalid, || false).is_err());
        assert_eq!(parcel, saved);
    }
}
