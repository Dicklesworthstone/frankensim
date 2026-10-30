#![allow(dead_code)]
//! Reuse manufactured original material claims; no planar geometry enters the
//! circular bridge or its numerical verification.
#[path = "../contact_line_bridge/common.rs"]
mod material_fixture;
use fs_contact::capillary_bridge::axisymmetric::contact_line::*;
use fs_contact::capillary_bridge::axisymmetric::*;

pub const VOLUME: f64 = material_fixture::VOLUME;
pub const SIGMA: f64 = material_fixture::SIGMA;
pub const INITIAL_Q: f64 = 2.0;

pub fn numerical() -> AxisymmetricNumerics {
    AxisymmetricNumerics {
        quadrature_relative_tolerance: 2e-13,
        gap_relative_tolerance: 2e-13,
        maximum_quadrature_panels: 4096,
        maximum_root_iterations: 48,
    }
}
pub fn config() -> AxisymmetricContactLineConfig {
    AxisymmetricContactLineConfig {
        shape: numerical(),
        minimum_angle_radians: 0.7,
        maximum_angle_radians: 1.4,
        minimum_neck_radius_m: 1e-6,
        maximum_gap_increment_m: 2e-5,
        maximum_radius_increment_m: 2e-5,
        maximum_step_s: 0.05,
        maximum_contact_line_speed_m_s: 0.02,
        maximum_angle_departure_radians: 0.35,
        maximum_capillary_number: 0.005,
        maximum_weber_number: 0.005,
        maximum_iterations: 24,
        maximum_backtracks: 20,
        maximum_shape_evaluations: 8192,
        geometry_absolute_tolerance_m: 1e-18,
        force_absolute_tolerance_n: 2e-12,
        force_relative_tolerance: 2e-9,
        energy_absolute_tolerance_j: 2e-18,
        energy_relative_tolerance: 2e-10,
    }
}
pub fn original(
    theta: f64,
    friction: [f64; 2],
    receding: f64,
    advancing: f64,
) -> SourceAxisymmetricCapillaryBridge {
    let material = material_fixture::original(theta, friction, receding, advancing);
    SourceAxisymmetricCapillaryBridge::new(
        material.liquid().clone(),
        material.wetting().clone(),
        VOLUME,
        AxisymmetricBranch {
            minimum_suction_parameter: 0.0,
            maximum_suction_parameter: 32.0,
            minimum_stability_margin: 0.01,
        },
        AxisymmetricLimits {
            minimum_gap_m: 1e-6,
            maximum_gap_m: 0.01,
            maximum_wetted_radius_m: 0.01,
            maximum_suction_pa: 1e7,
            maximum_abs_force_n: 1.0,
            maximum_abs_surface_energy_j: 0.01,
            maximum_abs_step_work_j: 0.001,
            gravity_m_s2: 0.0,
            maximum_bond_number: 0.01,
        },
    )
    .unwrap()
}
pub fn source(
    theta: f64,
    friction: f64,
    receding: f64,
    advancing: f64,
) -> SourceAxisymmetricContactLineBridge {
    SourceAxisymmetricContactLineBridge::new(
        original(theta, [friction; 2], receding, advancing),
        config(),
    )
    .unwrap()
}
pub fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}
