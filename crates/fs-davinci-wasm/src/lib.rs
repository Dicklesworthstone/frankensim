//! Browser boundary for the source-bounded US 6,331,181 manipulator topology.
//!
//! The generic joint law remains in `fs-mbd`; this crate performs browser
//! input admission and stable envelope serialization only.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_mbd::davinci::{DaVinciTopologyError, DaVinciTopologyParams, step_davinci_topology};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Admit one normalized tool-manipulator topology state.
///
/// Angles are radians. `insertion_normalized` is dimensionless because the
/// grant supplies no dimensioned manipulator travel.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn davinci_topology_step(
    base_yaw_rad: f64,
    carriage_pitch_rad: f64,
    distal_pitch_rad: f64,
    distal_yaw_rad: f64,
    tool_roll_rad: f64,
    insertion_normalized: f64,
    compatibility_identifier_present: bool,
) -> String {
    let result = match step_davinci_topology(DaVinciTopologyParams {
        base_yaw_rad,
        carriage_pitch_rad,
        distal_pitch_rad,
        distal_yaw_rad,
        tool_roll_rad,
        insertion_normalized,
        compatibility_identifier_present,
    }) {
        Ok(result) => result,
        Err(DaVinciTopologyError::InvalidInput) => {
            return refusal_json(
                "input-outside-domain",
                "angles must be finite and normalized insertion must be in [-1, 1]",
                "Use finite radians and clamp the dimensionless insertion coordinate",
            );
        }
        Err(DaVinciTopologyError::Multibody(_)) => {
            return refusal_json(
                "multibody-refusal",
                "the generic revolute/prismatic joint owner refused the composition",
                "Inspect the fs-mbd da Vinci topology contract",
            );
        }
    };

    let mut output = String::with_capacity(760);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"joint_dofs\":{},\
         \"base_yaw_axis\":[{},{},{}],\
         \"carriage_pitch_axis\":[{},{},{}],\
         \"insertion_axis\":[{},{},{}],\
         \"distal_pitch_axis\":[{},{},{}],\
         \"distal_yaw_axis\":[{},{},{}],\
         \"tool_roll_axis\":[{},{},{}],\
         \"base_yaw_rad\":{},\
         \"carriage_pitch_rad\":{},\
         \"distal_pitch_rad\":{},\
         \"distal_yaw_rad\":{},\
         \"tool_roll_rad\":{},\
         \"insertion_normalized\":{},\
         \"compatibility_identifier_present\":{}\
         }}}}",
        result.joint_dofs,
        result.base_yaw_axis[0],
        result.base_yaw_axis[1],
        result.base_yaw_axis[2],
        result.carriage_pitch_axis[0],
        result.carriage_pitch_axis[1],
        result.carriage_pitch_axis[2],
        result.insertion_axis[0],
        result.insertion_axis[1],
        result.insertion_axis[2],
        result.distal_pitch_axis[0],
        result.distal_pitch_axis[1],
        result.distal_pitch_axis[2],
        result.distal_yaw_axis[0],
        result.distal_yaw_axis[1],
        result.distal_yaw_axis[2],
        result.tool_roll_axis[0],
        result.tool_roll_axis[1],
        result.tool_roll_axis[2],
        result.base_yaw_rad,
        result.carriage_pitch_rad,
        result.distal_pitch_rad,
        result.distal_yaw_rad,
        result.tool_roll_rad,
        result.insertion_normalized,
        result.compatibility_identifier_present,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_generic_joint_axes_and_coordinates() {
        let output = davinci_topology_step(0.2, -0.3, 0.4, -0.1, 1.2, 0.25, true);
        assert!(output.contains("\"joint_dofs\":6"));
        assert!(output.contains("\"base_yaw_axis\":[0,1,0]"));
        assert!(output.contains("\"insertion_axis\":[0,-1,0]"));
        assert!(output.contains("\"compatibility_identifier_present\":true"));
    }

    #[test]
    fn absence_of_identifier_is_reported_without_rewriting_joint_state() {
        let output = davinci_topology_step(0.2, -0.3, 0.4, -0.1, 1.2, -0.4, false);
        assert!(output.contains("\"base_yaw_rad\":0.2"));
        assert!(output.contains("\"insertion_normalized\":-0.4"));
        assert!(output.contains("\"compatibility_identifier_present\":false"));
    }

    #[test]
    fn invalid_input_refuses() {
        let refusal = davinci_topology_step(0.0, 0.0, 0.0, 0.0, 0.0, 1.2, true);
        assert!(refusal.contains("\"code\":\"input-outside-domain\""));
        assert!(!refusal.contains("\"ok\""));
    }
}
