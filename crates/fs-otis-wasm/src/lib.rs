//! Browser boundary for the source-bounded US 31,128 hoisting topology.
//!
//! Generic revolute and prismatic joint ownership remains in `fs-mbd`; this
//! crate performs browser admission and stable envelope serialization only.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_mbd::otis::{OtisTopologyError, OtisTopologyParams, step_otis_topology};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Admit one source-order US 31,128 topology state.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn otis_topology_step(
    platform_position_normalized: f64,
    drive_phase_rad: f64,
    drive_command: i8,
    rope_g_intact: bool,
    stop_rope_pulled: bool,
    claim_1_hook_lock_enabled: bool,
    claim_3_brake_interlock_enabled: bool,
    claim_4_counterpoise_enabled: bool,
) -> String {
    let result = match step_otis_topology(OtisTopologyParams {
        platform_position_normalized,
        drive_phase_rad,
        drive_command,
        rope_g_intact,
        stop_rope_pulled,
        claim_1_hook_lock_enabled,
        claim_3_brake_interlock_enabled,
        claim_4_counterpoise_enabled,
    }) {
        Ok(result) => result,
        Err(OtisTopologyError::InvalidInput) => {
            return refusal_json(
                "input-outside-domain",
                "platform position must be in [0,1], phase finite, and drive command -1, 0, or 1",
                "Clamp platform travel, use a finite phase, and choose a declared drive command",
            );
        }
        Err(OtisTopologyError::Multibody(_)) => {
            return refusal_json(
                "multibody-refusal",
                "the generic revolute/prismatic joint owner refused the composition",
                "Inspect the fs-mbd Otis topology contract",
            );
        }
    };

    let mut output = String::with_capacity(2_000);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"scalar_joint_coordinates\":{},\
         \"independent_drive_dofs\":{},\
         \"platform_axis\":[{},{},{}],\
         \"safety_bar_axis\":[{},{},{}],\
         \"safety_lever_axis\":[{},{},{}],\
         \"winding_drum_axis\":[{},{},{}],\
         \"shipper_axis\":[{},{},{}],\
         \"brake_axis\":[{},{},{}],\
         \"counterpoise_axis\":[{},{},{}],\
         \"platform_position_normalized\":{},\
         \"counterpoise_position_normalized\":{},\
         \"drive_phase_rad\":{},\
         \"requested_drive_direction\":{},\
         \"platform_motion_direction\":{},\
         \"shipper_position_normalized\":{},\
         \"straight_belt_o_working\":{},\
         \"cross_belt_p_working\":{},\
         \"both_belts_idle\":{},\
         \"brake_z_engaged\":{},\
         \"stop_rope_geometry_active\":{},\
         \"lower_limit_stop_active\":{},\
         \"rope_g_taut\":{},\
         \"safety_bar_release_normalized\":{},\
         \"safety_lever_rotation_normalized\":{},\
         \"pawls_f_engaged\":{},\
         \"claim_1_hook_lock_satisfied\":{},\
         \"free_fall_counterfactual\":{},\
         \"claim_3_stop_interlock_satisfied\":{},\
         \"claim_4_counterpoise_topology_satisfied\":{},\
         \"mechanism_mode\":\"{}\"\
         }}}}",
        result.scalar_joint_coordinates,
        result.independent_drive_dofs,
        result.platform_axis[0],
        result.platform_axis[1],
        result.platform_axis[2],
        result.safety_bar_axis[0],
        result.safety_bar_axis[1],
        result.safety_bar_axis[2],
        result.safety_lever_axis[0],
        result.safety_lever_axis[1],
        result.safety_lever_axis[2],
        result.winding_drum_axis[0],
        result.winding_drum_axis[1],
        result.winding_drum_axis[2],
        result.shipper_axis[0],
        result.shipper_axis[1],
        result.shipper_axis[2],
        result.brake_axis[0],
        result.brake_axis[1],
        result.brake_axis[2],
        result.counterpoise_axis[0],
        result.counterpoise_axis[1],
        result.counterpoise_axis[2],
        result.platform_position_normalized,
        result.counterpoise_position_normalized,
        result.drive_phase_rad,
        result.requested_drive_direction,
        result.platform_motion_direction,
        result.shipper_position_normalized,
        result.straight_belt_o_working,
        result.cross_belt_p_working,
        result.both_belts_idle,
        result.brake_z_engaged,
        result.stop_rope_geometry_active,
        result.lower_limit_stop_active,
        result.rope_g_taut,
        result.safety_bar_release_normalized,
        result.safety_lever_rotation_normalized,
        result.pawls_f_engaged,
        result.claim_1_hook_lock_satisfied,
        result.free_fall_counterfactual,
        result.claim_3_stop_interlock_satisfied,
        result.claim_4_counterpoise_topology_satisfied,
        result.mechanism_mode,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_the_complete_generic_joint_topology() {
        let output = otis_topology_step(0.55, 0.0, 1, true, false, true, true, true);
        assert!(output.contains("\"scalar_joint_coordinates\":12"));
        assert!(output.contains("\"platform_axis\":[0,1,0]"));
        assert!(output.contains("\"straight_belt_o_working\":true"));
        assert!(output.contains("\"mechanism_mode\":\"raise\""));
    }

    #[test]
    fn rope_failure_proves_the_hook_lock_and_counterfactual() {
        let caught = otis_topology_step(0.55, 0.0, 1, false, false, true, true, true);
        assert!(caught.contains("\"pawls_f_engaged\":true"));
        assert!(caught.contains("\"platform_motion_direction\":0"));

        let removed = otis_topology_step(0.55, 0.0, 1, false, false, false, true, true);
        assert!(removed.contains("\"free_fall_counterfactual\":true"));
        assert!(removed.contains("\"platform_motion_direction\":-1"));
    }

    #[test]
    fn invalid_travel_and_command_refuse() {
        let bad_travel = otis_topology_step(1.2, 0.0, 1, true, false, true, true, true);
        assert!(bad_travel.contains("\"code\":\"input-outside-domain\""));
        let bad_command = otis_topology_step(0.5, 0.0, 2, true, false, true, true, true);
        assert!(bad_command.contains("\"code\":\"input-outside-domain\""));
    }
}
