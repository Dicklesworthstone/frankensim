//! Browser boundary for the source-bounded US 4,750 sewing-machine topology.
//!
//! Generic revolute and prismatic joint ownership remains in `fs-mbd`; this
//! crate performs browser admission and stable envelope serialization only.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_mbd::howe::{HoweTopologyError, HoweTopologyParams, step_howe_topology};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Admit one source-order US 4,750 topology state.
///
/// `crank_angle_rad` is a finite prescribed coordinate and
/// `loop_slack_normalized` is a dimensionless display input in `[0, 1]`.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
pub fn howe_topology_step(
    crank_angle_rad: f64,
    loop_slack_normalized: f64,
    claim_1_interlock_enabled: bool,
) -> String {
    let result = match step_howe_topology(HoweTopologyParams {
        crank_angle_rad,
        loop_slack_normalized,
        claim_1_interlock_enabled,
    }) {
        Ok(result) => result,
        Err(HoweTopologyError::InvalidInput) => {
            return refusal_json(
                "input-outside-domain",
                "main-shaft angle must be finite and normalized loop slack must be in [0, 1]",
                "Use a finite radian coordinate and clamp displayed loop slack to [0, 1]",
            );
        }
        Err(HoweTopologyError::Multibody(_)) => {
            return refusal_json(
                "multibody-refusal",
                "the generic revolute/prismatic joint owner refused the composition",
                "Inspect the fs-mbd Howe topology contract",
            );
        }
    };

    let mut output = String::with_capacity(1_500);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"scalar_joint_coordinates\":{},\
         \"independent_drive_dofs\":{},\
         \"main_shaft_axis\":[{},{},{}],\
         \"needle_arm_axis\":[{},{},{}],\
         \"shuttle_axis\":[{},{},{}],\
         \"lifting_rod_axis\":[{},{},{}],\
         \"baster_feed_axis\":[{},{},{}],\
         \"crank_angle_rad\":{},\
         \"needle_penetration_normalized\":{},\
         \"needle_arm_angle_rad\":{},\
         \"needle_retracting\":{},\
         \"shuttle_travel_normalized\":{},\
         \"loop_open_fraction\":{},\
         \"loop_open\":{},\
         \"shuttle_passes_loop\":{},\
         \"shuttle_track_offset_normalized\":{},\
         \"picker_left_normalized\":{},\
         \"picker_right_normalized\":{},\
         \"lifting_rod_normalized\":{},\
         \"feed_advance_fraction\":{},\
         \"thread_clamp_engaged\":{},\
         \"claim_1_interlock_satisfied\":{},\
         \"cycle_phase\":\"{}\",\
         \"needle_eye_offset_in\":{},\
         \"baster_point_pitch_in\":{}\
         }}}}",
        result.scalar_joint_coordinates,
        result.independent_drive_dofs,
        result.main_shaft_axis[0],
        result.main_shaft_axis[1],
        result.main_shaft_axis[2],
        result.needle_arm_axis[0],
        result.needle_arm_axis[1],
        result.needle_arm_axis[2],
        result.shuttle_axis[0],
        result.shuttle_axis[1],
        result.shuttle_axis[2],
        result.lifting_rod_axis[0],
        result.lifting_rod_axis[1],
        result.lifting_rod_axis[2],
        result.baster_feed_axis[0],
        result.baster_feed_axis[1],
        result.baster_feed_axis[2],
        result.crank_angle_rad,
        result.needle_penetration_normalized,
        result.needle_arm_angle_rad,
        result.needle_retracting,
        result.shuttle_travel_normalized,
        result.loop_open_fraction,
        result.loop_open,
        result.shuttle_passes_loop,
        result.shuttle_track_offset_normalized,
        result.picker_left_normalized,
        result.picker_right_normalized,
        result.lifting_rod_normalized,
        result.feed_advance_fraction,
        result.thread_clamp_engaged,
        result.claim_1_interlock_satisfied,
        result.cycle_phase,
        result.needle_eye_offset_in,
        result.baster_point_pitch_in,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_generic_axes_and_source_order_pass_state() {
        let output = howe_topology_step(1.5 * core::f64::consts::PI, 0.65, true);
        assert!(output.contains("\"scalar_joint_coordinates\":7"));
        assert!(output.contains("\"independent_drive_dofs\":1"));
        assert!(output.contains("\"main_shaft_axis\":[0,0,1]"));
        assert!(output.contains("\"shuttle_axis\":[1,0,0]"));
        assert!(output.contains("\"shuttle_passes_loop\":true"));
        assert!(output.contains("\"cycle_phase\":\"shuttle-pass\""));
    }

    #[test]
    fn removed_claim_reports_guided_counterfactual_without_an_interlock() {
        let output = howe_topology_step(1.5 * core::f64::consts::PI, 0.65, false);
        assert!(output.contains("\"shuttle_passes_loop\":false"));
        assert!(output.contains("\"claim_1_interlock_satisfied\":false"));
        assert!(output.contains("\"shuttle_track_offset_normalized\":0.55"));
    }

    #[test]
    fn invalid_slack_refuses() {
        let refusal = howe_topology_step(0.0, 1.2, true);
        assert!(refusal.contains("\"code\":\"input-outside-domain\""));
        assert!(!refusal.contains("\"ok\""));
    }
}
