//! Shared numerical budgets for native and finite-heat electrode circuit tests.
use fs_em::spatial::vector_potential::induction::terminals::circuit::{
    SpatialTerminalCircuitConfig, SpatialTerminalCircuitInitial,
    SpatialTerminalCircuitNonlinearControls,
};
use fs_solver::NewtonKrylovConfig;

pub fn config(reference: f64) -> SpatialTerminalCircuitConfig {
    SpatialTerminalCircuitConfig {
        reference_resistance_ohm: reference,
        maximum_native_states: 128,
        maximum_interface_scalars: 100_000,
    }
}
pub fn controls<F>(field: F) -> SpatialTerminalCircuitNonlinearControls<F> {
    SpatialTerminalCircuitNonlinearControls {
        field,
        difference_step_v: 1e-4,
        maximum_field_trials: 300,
        maximum_newton_iterations: 40,
        voltage_tolerance_v: 2e-8,
        current_tolerance_a: 2e-8,
        energy_tolerance_j: 2e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
pub fn initial(
    fluxes: Vec<f64>,
    charges: Vec<f64>,
    inputs: Vec<f64>,
) -> SpatialTerminalCircuitInitial {
    SpatialTerminalCircuitInitial {
        fluxes_wb: fluxes,
        charges_c: charges,
        external_inputs: inputs,
    }
}

/// Scalar oracle only: no production matrix, constitutive solver or trajectory.
pub fn bracket(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..100 {
        let mid = f64::midpoint(lo, hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    f64::midpoint(lo, hi)
}
