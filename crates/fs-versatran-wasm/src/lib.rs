//! Browser boundary for the source-bounded US 3,212,649 Versatran topology.
//!
//! `fs-mbd` owns the generic revolute/prismatic composition. This crate admits
//! browser inputs and serializes a stable receipt without adding geometry,
//! inverse kinematics, dynamics, hydraulics, or program execution.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_mbd::versatran::{VersatranTopologyError, VersatranTopologyParams, step_versatran_topology};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Admit one source-bounded AMF Versatran topology state.
///
/// Rotations are radians. The vertical-arm, horizontal-arm, and internal
/// work-member-rack coordinates are dimensionless normalized presentation
/// coordinates because US 3,212,649 does not print their travel dimensions.
/// Selecting `automatic_program_mode_selected` reports the patent's mode;
/// this function does not execute or synthesize a program.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn versatran_topology_step(
    column_rotation_rad: f64,
    arm_vertical_normalized: f64,
    arm_horizontal_normalized: f64,
    wrist_rotation_rad: f64,
    wrist_swing_rad: f64,
    work_member_rack_normalized: f64,
    automatic_program_mode_selected: bool,
) -> String {
    let result = match step_versatran_topology(VersatranTopologyParams {
        column_rotation_rad,
        arm_vertical_normalized,
        arm_horizontal_normalized,
        wrist_rotation_rad,
        wrist_swing_rad,
        work_member_rack_normalized,
        automatic_program_mode_selected,
    }) {
        Ok(result) => result,
        Err(VersatranTopologyError::InvalidInput) => {
            return refusal_json(
                "input-outside-domain",
                "rotations must be finite and normalized arm/rack coordinates must be in [0, 1]",
                "Use finite radians and clamp each dimensionless coordinate to its declared interval",
            );
        }
        Err(VersatranTopologyError::Multibody(_)) => {
            return refusal_json(
                "multibody-refusal",
                "the generic revolute/prismatic joint owner refused the composition",
                "Inspect the fs-mbd Versatran topology contract",
            );
        }
    };

    let mut output = String::with_capacity(1_200);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"scalar_channels\":{},\
         \"geometric_motion_joints\":{},\
         \"work_member_operation_channels\":{},\
         \"revolute_joint_count\":{},\
         \"prismatic_joint_count\":{},\
         \"column_rotation_axis\":[{},{},{}],\
         \"arm_vertical_axis\":[{},{},{}],\
         \"arm_horizontal_axis\":[{},{},{}],\
         \"wrist_rotation_axis\":[{},{},{}],\
         \"wrist_swing_axis\":[{},{},{}],\
         \"work_member_rack_axis\":[{},{},{}],\
         \"column_rotation_rad\":{},\
         \"arm_vertical_normalized\":{},\
         \"arm_horizontal_normalized\":{},\
         \"wrist_rotation_rad\":{},\
         \"wrist_swing_rad\":{},\
         \"work_member_rack_normalized\":{},\
         \"automatic_program_mode_selected\":{},\
         \"historical_geometry_available\":{},\
         \"historical_dynamics_available\":{}\
         }}}}",
        result.scalar_channels,
        result.geometric_motion_joints,
        result.work_member_operation_channels,
        result.revolute_joint_count,
        result.prismatic_joint_count,
        result.column_rotation_axis[0],
        result.column_rotation_axis[1],
        result.column_rotation_axis[2],
        result.arm_vertical_axis[0],
        result.arm_vertical_axis[1],
        result.arm_vertical_axis[2],
        result.arm_horizontal_axis[0],
        result.arm_horizontal_axis[1],
        result.arm_horizontal_axis[2],
        result.wrist_rotation_axis[0],
        result.wrist_rotation_axis[1],
        result.wrist_rotation_axis[2],
        result.wrist_swing_axis[0],
        result.wrist_swing_axis[1],
        result.wrist_swing_axis[2],
        result.work_member_rack_axis[0],
        result.work_member_rack_axis[1],
        result.work_member_rack_axis[2],
        result.column_rotation_rad,
        result.arm_vertical_normalized,
        result.arm_horizontal_normalized,
        result.wrist_rotation_rad,
        result.wrist_swing_rad,
        result.work_member_rack_normalized,
        result.automatic_program_mode_selected,
        result.historical_geometry_available,
        result.historical_dynamics_available,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_the_bounded_generic_joint_receipt() {
        let output = versatran_topology_step(0.2, 0.35, 0.65, -0.4, 0.3, 0.8, true);
        assert!(output.contains("\"scalar_channels\":6"));
        assert!(output.contains("\"geometric_motion_joints\":5"));
        assert!(output.contains("\"work_member_operation_channels\":1"));
        assert!(output.contains("\"revolute_joint_count\":3"));
        assert!(output.contains("\"prismatic_joint_count\":3"));
        assert!(output.contains("\"column_rotation_axis\":[0,1,0]"));
        assert!(output.contains("\"work_member_rack_axis\":[1,0,0]"));
        assert!(output.contains("\"historical_geometry_available\":false"));
        assert!(output.contains("\"historical_dynamics_available\":false"));
    }

    #[test]
    fn program_mode_is_reported_without_rewriting_joint_state() {
        let automatic = versatran_topology_step(0.2, 0.35, 0.65, -0.4, 0.3, 0.8, true);
        let manual = versatran_topology_step(0.2, 0.35, 0.65, -0.4, 0.3, 0.8, false);
        assert!(automatic.contains("\"automatic_program_mode_selected\":true"));
        assert!(manual.contains("\"automatic_program_mode_selected\":false"));
        assert!(manual.contains("\"work_member_rack_normalized\":0.8"));
    }

    #[test]
    fn invalid_input_refuses() {
        let out_of_domain = versatran_topology_step(0.0, 1.2, 0.0, 0.0, 0.0, 0.0, false);
        assert!(out_of_domain.contains("\"code\":\"input-outside-domain\""));
        assert!(!out_of_domain.contains("\"ok\""));

        let non_finite = versatran_topology_step(f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, false);
        assert!(non_finite.contains("\"code\":\"input-outside-domain\""));
    }
}
