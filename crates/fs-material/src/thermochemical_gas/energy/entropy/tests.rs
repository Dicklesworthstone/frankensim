//! Independent synthetic caloric identities, not measured gas qualification.
use super::*;
use crate::thermochemical_gas::tests::{close, model, point};

fn controls() -> GasEntropySolve {
    GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128 }
}

#[test]
fn density_entropy_roundtrip_preserves_complete_material_evaluation() {
    let model = model(vec![0.75, 0.25], 0.001);
    let context = point(600.0, 200000.0).with("material-lot", 12.0).unwrap();
    let original = model.evaluate(&context).unwrap();
    let result = model.reconstruct_density_entropy(&context, original.state().density,
        original.specific_entropy_j_kg_k().unwrap(), controls(), || false).unwrap();
    close(result.state().state().temperature, 600.0);
    close(result.state().state().pressure, 200000.0);
    close(result.state().state().dynamic_viscosity, original.state().dynamic_viscosity);
    close(result.state().state().thermal_conductivity, original.state().thermal_conductivity);
    assert!(result.residual_j_kg_k().abs() <= controls().entropy_tolerance_j_kg_k);
    assert_eq!(result.state().components().len(), 2);
    let energy = model.reconstruct_density_energy(&context, result.state().state().density,
        result.state().specific_internal_energy_j_kg(), super::super::GasEnergySolve {
            minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
            energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
        }, || false).unwrap();
    close(energy.state().specific_entropy_j_kg_k().unwrap(), original.specific_entropy_j_kg_k().unwrap());
}

#[test]
fn variable_heat_capacity_uses_integrated_entropy_not_fixed_gamma() {
    let model = model(vec![0.75, 0.25], 0.001);
    let context = point(300.0, 100000.0);
    let initial = model.evaluate(&context).unwrap();
    // cp/R = a + b*T. Independently integrating ds at fixed composition gives
    // ln(rho1/rho0) = (a-1)*ln(T1/T0) + b*(T1-T0).
    let (a, b) = (3.25, 0.00075);
    let density_ratio = ((a - 1.0) * (600.0_f64/300.0).ln() + b * 300.0).exp();
    let result = model.reconstruct_density_entropy(&context, initial.state().density*density_ratio,
        initial.specific_entropy_j_kg_k().unwrap(), controls(), || false).unwrap();
    close(result.state().state().temperature, 600.0);
    close(result.state().state().pressure, 100000.0*density_ratio*2.0);
    assert!((300.0*density_ratio.powf(initial.state().gamma - 1.0) - 600.0).abs() > 10.0);
    assert!(result.state().state().gamma < initial.state().gamma);
}

#[test]
fn isentropic_pressure_derivative_matches_shared_acoustic_sound_speed() {
    let model = model(vec![0.25, 0.75], 0.001);
    let context = point(500.0, 100000.0);
    let initial = model.evaluate(&context).unwrap();
    let rho = initial.state().density;
    let epsilon = 1e-4;
    let solve = |density| model.reconstruct_density_entropy(&context, density,
        initial.specific_entropy_j_kg_k().unwrap(), controls(), || false).unwrap();
    let lo = solve(rho*(1.0-epsilon));
    let hi = solve(rho*(1.0+epsilon));
    let dp_drho = (hi.state().state().pressure-lo.state().state().pressure)/(2.0*rho*epsilon);
    let c2 = initial.state().sound_speed.powi(2);
    assert!((dp_drho/c2 - 1.0).abs() < 1e-7);
}

#[test]
fn bad_states_brackets_work_limits_and_cancellation_refuse() {
    let model = model(vec![0.75, 0.25], 0.001);
    let context = point(300.0, 100000.0);
    let initial = model.evaluate(&context).unwrap();
    let rho = initial.state().density;
    let s = initial.specific_entropy_j_kg_k().unwrap();
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s, controls(), || true),
        Err(GasEntropyError::Cancelled)));
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s+1e9, controls(), || false),
        Err(GasEntropyError::OutsideBracket { .. })));
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s,
        GasEntropySolve { maximum_iterations: 1, ..controls() }, || false),
        Err(GasEntropyError::NotConverged { .. })));
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.reconstruct_density_entropy(&context, density, s, controls(), || false).is_err());
    }
    for bad in [GasEntropySolve { maximum_iterations: 0, ..controls() },
        GasEntropySolve { entropy_tolerance_j_kg_k: 0.0, ..controls() },
        GasEntropySolve { minimum_temperature_k: 1600.0, ..controls() }] {
        assert!(model.reconstruct_density_entropy(&context, rho, s, bad, || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(model.reconstruct_density_entropy(&context, rho, s, controls(), || {
        polls += 1; polls == 4
    }), Err(GasEntropyError::Cancelled)));
}
