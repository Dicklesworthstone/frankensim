//! The coupled stiffness screen refuses a network above the Nyquist guard with both figures (G0).
use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use fs_math::c64::C64;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 4, max_connections: 4, max_setup_terms: 100_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9 }
}

// One 20 rad/s mode at a 1 kHz clock, joined to the world by a unit-participation spring of the given stiffness.
fn network(stiffness_n_m: f64) -> Result<CoupledModalSystem, ModalCouplingError> {
    let gate = CancelGate::new_clock_free();
    let model = ModalAcousticTimeModel::try_new(1000,
        vec![ModalAcousticMode { angular_frequency_rad_s: 20.0, damping_ratio: 0.01,
            pressure_per_modal_velocity: C64::ZERO }],
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    CoupledModalSystem::new_with_columns(vec![model], vec![ModalColumnConnection {
        column: vec![1.0], stiffness_n_m, damping_n_s_m: 0.0, rest_extension_m: 0.0 }], config(), &gate)
}

// The screen is lambda_max <= max(omega^2) + sum k ||b||^2, held to pi/dt x nyquist_guard_fraction (dt = 1 ms here).
#[test]
fn a_stiffness_bound_above_the_nyquist_guard_is_refused_with_its_value_and_the_guard() {
    let guard = core::f64::consts::PI * 1000.0 * 0.9;
    // k ||b||^2 = 1e8: the bound is sqrt(20^2 + 1e8), about 1.0e4 rad/s, above the 2.83e3 rad/s guard.
    let error = network(1.0e8).err().expect("a stiffness bound above the guard was admitted");
    let ModalCouplingError::Budget { what, value, limit } = error else {
        panic!("a bound above the guard was not refused with both figures: {error:?}");
    };
    assert!(what.contains("Nyquist guard"), "{what}");
    assert!((value - (400.0_f64 + 1.0e8).sqrt()).abs() <= 1e-9 * value, "bound {value:e}");
    assert!((limit - guard).abs() <= 1e-12 * guard, "guard {limit:e}");
    let text = ModalCouplingError::Budget { what, value, limit }.to_string();
    assert!(text.contains("exceeds") && text.contains(&format!("{value:e}")) && text.contains(&format!("{limit:e}")),
        "{text}");
    // Under the guard the same network is admitted: sqrt(20^2 + 2e6) is about 1.41e3 rad/s.
    assert!(network(2.0e6).is_ok());
}
