//! Browser boundary for the source-bounded US 194,047 Otto topology.
//!
//! Generic revolute/prismatic joints and slider-crank closure remain in
//! `fs-mbd`; this crate performs browser admission and stable envelope
//! serialization only.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_mbd::otto::{OttoTopologyError, OttoTopologyParams, step_otto_topology};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Admit and close one source-order US 194,047 mechanism pose.
///
/// Dimensions are finite, positive caller-declared display geometry in a
/// consistent unit. `engine_rpm` is admitted in the museum range `[0, 600]`.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
pub fn otto_topology_step(
    crank_angle_rad: f64,
    crank_radius: f64,
    connecting_rod_length: f64,
    engine_rpm: f64,
) -> String {
    let result = match step_otto_topology(OttoTopologyParams {
        crank_angle_rad,
        crank_radius,
        connecting_rod_length,
        engine_rpm,
    }) {
        Ok(result) => result,
        Err(OttoTopologyError::InvalidInput) => {
            return refusal_json(
                "input-outside-domain",
                "angle and geometry must be finite, geometry positive, and engine speed in [0, 600] rpm",
                "Use finite display geometry and a speed inside the registered museum control domain",
            );
        }
        Err(OttoTopologyError::ImpossibleLinkage) => {
            return refusal_json(
                "impossible-linkage",
                "the connecting rod must be longer than the crank radius",
                "Increase rod length or reduce crank radius before stepping",
            );
        }
        Err(OttoTopologyError::Multibody(_)) => {
            return refusal_json(
                "multibody-refusal",
                "the generic revolute/prismatic joint owner refused the composition",
                "Inspect the fs-mbd Otto topology contract",
            );
        }
    };

    let mut output = String::with_capacity(1_500);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"scalar_joint_coordinates\":{},\
         \"independent_drive_dofs\":{},\
         \"crank_axis\":[{},{},{}],\
         \"piston_axis\":[{},{},{}],\
         \"side_shaft_axis\":[{},{},{}],\
         \"slide_valve_axis\":[{},{},{}],\
         \"exhaust_valve_axis\":[{},{},{}],\
         \"governor_axis\":[{},{},{}],\
         \"cycle_angle_rad\":{},\
         \"crank_pin_x\":{},\
         \"crank_pin_y\":{},\
         \"piston_pin_x\":{},\
         \"piston_pin_y\":{},\
         \"connecting_rod_angle_rad\":{},\
         \"connecting_rod_span\":{},\
         \"side_shaft_angle_rad\":{},\
         \"slide_valve_normalized\":{},\
         \"exhaust_lift_normalized\":{},\
         \"governor_spread_normalized\":{},\
         \"cycle_phase\":\"{}\"\
         }}}}",
        result.scalar_joint_coordinates,
        result.independent_drive_dofs,
        result.crank_axis[0],
        result.crank_axis[1],
        result.crank_axis[2],
        result.piston_axis[0],
        result.piston_axis[1],
        result.piston_axis[2],
        result.side_shaft_axis[0],
        result.side_shaft_axis[1],
        result.side_shaft_axis[2],
        result.slide_valve_axis[0],
        result.slide_valve_axis[1],
        result.slide_valve_axis[2],
        result.exhaust_valve_axis[0],
        result.exhaust_valve_axis[1],
        result.exhaust_valve_axis[2],
        result.governor_axis[0],
        result.governor_axis[1],
        result.governor_axis[2],
        result.cycle_angle_rad,
        result.crank_pin_x,
        result.crank_pin_y,
        result.piston_pin_x,
        result.piston_pin_y,
        result.connecting_rod_angle_rad,
        result.connecting_rod_span,
        result.side_shaft_angle_rad,
        result.slide_valve_normalized,
        result.exhaust_lift_normalized,
        result.governor_spread_normalized,
        result.cycle_phase,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_closed_generic_joint_pose() {
        let output = otto_topology_step(core::f64::consts::FRAC_PI_2, 0.65, 2.4, 180.0);
        assert!(output.contains("\"scalar_joint_coordinates\":8"));
        assert!(output.contains("\"independent_drive_dofs\":1"));
        assert!(output.contains("\"piston_axis\":[1,0,0]"));
        assert!(output.contains("\"connecting_rod_span\":2.4"));
        assert!(output.contains("\"cycle_phase\":\"intake\""));
    }

    #[test]
    fn emits_exhaust_lift_from_half_speed_source_order() {
        let output = otto_topology_step(3.5 * core::f64::consts::PI, 0.65, 2.4, 180.0);
        assert!(output.contains("\"cycle_phase\":\"exhaust\""));
        assert!(output.contains("\"exhaust_lift_normalized\":1"));
    }

    #[test]
    fn invalid_geometry_refuses_without_partial_success() {
        let refusal = otto_topology_step(0.0, 0.65, 0.65, 180.0);
        assert!(refusal.contains("\"code\":\"impossible-linkage\""));
        assert!(!refusal.contains("\"ok\""));
    }
}
