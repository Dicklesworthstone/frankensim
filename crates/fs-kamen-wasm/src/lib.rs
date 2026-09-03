//! Browser boundary for the source-dimensioned US 5,701,965 tri-wheel states.
//!
//! Table 1 owns the nominal dimensions and Figures 39--42 own the discrete
//! teaching poses. Generic rigid transforms and horizontal-support gap checks
//! remain owned by `fs-mbd`; this crate performs source mapping and envelope
//! serialization only.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::f64::consts::{PI, TAU};
use core::fmt::Write as _;
use fs_mbd::tri_wheel_cluster::{
    TriWheelStairError, TriWheelStairInput, step_tri_wheel_stair_contact,
};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

const INCH_M: f64 = 0.0254;
const SYSTEM_CENTRE_OFFSET_M: f64 = 21.0 * INCH_M;
const CLUSTER_RADIUS_M: f64 = 5.581 * INCH_M;
const ADJACENT_WHEEL_CENTRE_DISTANCE_M: f64 = 9.667 * INCH_M;
const STAIR_TREAD_M: f64 = 10.9 * INCH_M;
const STAIR_RISE_M: f64 = 6.85 * INCH_M;
const RISER_TO_LOWER_CONTACT_M: f64 = 3.011 * INCH_M;
const WHEEL_RADIUS_M: f64 = 3.81 * INCH_M;
const CONTACT_TOLERANCE_M: f64 = 1.0e-8;
const GENERIC_OWNER: &str = "fs-mbd::tri_wheel_cluster::step_tri_wheel_stair_contact";
const MODEL_BOUNDARY: &str = "rigid-planar-three-equal-wheels-horizontal-tread-contact-no-force-friction-compliance-impact-or-riser-side-contact";
const SOURCE_RECEIPT: &str = "us-5701965-table-1-figures-39-through-42";

#[derive(Clone, Copy)]
struct SourcePose {
    name: &'static str,
    figure: &'static str,
    axle_x_m: f64,
    axle_y_m: f64,
    carrier_rotation_rad: f64,
    chassis_pitch_rad: f64,
    stair_active: bool,
}

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

fn source_start_rotation_rad() -> f64 {
    let adjacent_distance_from_radius_m = 3.0_f64.sqrt() * CLUSTER_RADIUS_M;
    -(PI / 3.0 - (STAIR_RISE_M / adjacent_distance_from_radius_m).asin())
}

fn source_pose(state_index: u8) -> Option<SourcePose> {
    let start_rotation = source_start_rotation_rad();
    let wheel_a_angle = -PI / 2.0 + start_rotation;
    let start_axle_x = -RISER_TO_LOWER_CONTACT_M - CLUSTER_RADIUS_M * wheel_a_angle.cos();
    let start_axle_y = WHEEL_RADIUS_M - CLUSTER_RADIUS_M * wheel_a_angle.sin();
    let start_pitch = wheel_a_angle + (TAU - 2.814) - PI / 2.0;

    let wheel_b_angle = PI / 6.0 + start_rotation;
    let transfer_pitch = wheel_b_angle + (TAU - 5.236) - PI / 2.0;

    let climb_rotation = start_rotation - 2.0 * PI / 3.0;
    let climb_wheel_b_angle = PI / 6.0 + climb_rotation;
    let climb_axle_x =
        STAIR_TREAD_M - RISER_TO_LOWER_CONTACT_M - CLUSTER_RADIUS_M * climb_wheel_b_angle.cos();
    let climb_axle_y = STAIR_RISE_M + WHEEL_RADIUS_M - CLUSTER_RADIUS_M * climb_wheel_b_angle.sin();
    let climb_pitch = climb_wheel_b_angle + (TAU - 2.814) - PI / 2.0;

    match state_index {
        0 => Some(SourcePose {
            name: "ground_support",
            figure: "claim-16-four-ground-wheels-comparison",
            axle_x_m: 0.0,
            axle_y_m: WHEEL_RADIUS_M + CLUSTER_RADIUS_M / 2.0,
            carrier_rotation_rad: -PI / 3.0,
            chassis_pitch_rad: 0.0,
            stair_active: false,
        }),
        1 => Some(SourcePose {
            name: "balance",
            figure: "figure-39a",
            axle_x_m: 0.0,
            axle_y_m: WHEEL_RADIUS_M + CLUSTER_RADIUS_M,
            carrier_rotation_rad: 0.0,
            chassis_pitch_rad: 0.0,
            stair_active: false,
        }),
        2 => Some(SourcePose {
            name: "stair_start",
            figure: "figure-39b",
            axle_x_m: start_axle_x,
            axle_y_m: start_axle_y,
            carrier_rotation_rad: start_rotation,
            chassis_pitch_rad: start_pitch,
            stair_active: true,
        }),
        3 => Some(SourcePose {
            name: "weight_transfer",
            figure: "figure-41b",
            axle_x_m: start_axle_x,
            axle_y_m: start_axle_y,
            carrier_rotation_rad: start_rotation,
            chassis_pitch_rad: transfer_pitch,
            stair_active: true,
        }),
        4 => Some(SourcePose {
            name: "climb",
            figure: "figure-42c",
            axle_x_m: climb_axle_x,
            axle_y_m: climb_axle_y,
            carrier_rotation_rad: climb_rotation,
            chassis_pitch_rad: climb_pitch,
            stair_active: true,
        }),
        5 => Some(SourcePose {
            name: "transition",
            figure: "figure-38-zero-crossing-on-upper-tread",
            axle_x_m: 1.5 * STAIR_TREAD_M,
            axle_y_m: 2.0 * STAIR_RISE_M + WHEEL_RADIUS_M + CLUSTER_RADIUS_M,
            carrier_rotation_rad: -4.0 * PI / 3.0,
            chassis_pitch_rad: 0.0,
            stair_active: true,
        }),
        _ => None,
    }
}

fn contact_error_json(error: TriWheelStairError) -> String {
    match error {
        TriWheelStairError::InvalidInput(_) => refusal_json(
            "invalid-source-pose",
            "the source-dimensioned tri-wheel pose contains an invalid scalar",
            "Inspect the Table 1 dimension mapping",
        ),
        TriWheelStairError::PenetratingSupport { .. } => refusal_json(
            "support-penetration",
            "the rigid wheel geometry penetrates a horizontal support",
            "Correct the source-pose transform before rendering",
        ),
        TriWheelStairError::Unsupported { .. } => refusal_json(
            "unsupported-pose",
            "no wheel touches a horizontal support in the selected state",
            "Correct the source-pose transform before rendering",
        ),
    }
}

/// Resolve one source-dimensioned tri-wheel support state.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
pub fn kamen_cluster_step(state_index: u8) -> String {
    let Some(pose) = source_pose(state_index) else {
        return refusal_json(
            "state-outside-domain",
            "state index must be in the inclusive range 0 through 5",
            "Choose one of the six declared source-reading states",
        );
    };

    let result = match step_tri_wheel_stair_contact(TriWheelStairInput {
        cluster_radius_m: CLUSTER_RADIUS_M,
        wheel_radius_m: WHEEL_RADIUS_M,
        axle_x_m: pose.axle_x_m,
        axle_y_m: pose.axle_y_m,
        carrier_rotation_rad: pose.carrier_rotation_rad,
        stair_rise_m: STAIR_RISE_M,
        stair_tread_m: STAIR_TREAD_M,
        stair_active: pose.stair_active,
        contact_tolerance_m: CONTACT_TOLERANCE_M,
    }) {
        Ok(result) => result,
        Err(error) => return contact_error_json(error),
    };

    let mut output = String::with_capacity(1_800);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"owner\":\"{}\",\
         \"boundary\":\"{}\",\
         \"source_receipt\":\"{}\",\
         \"state\":\"{}\",\
         \"source_figure\":\"{}\",\
         \"system_centre_offset_m\":{},\
         \"cluster_radius_m\":{},\
         \"adjacent_wheel_centre_distance_m\":{},\
         \"wheel_radius_m\":{},\
         \"stair_rise_m\":{},\
         \"stair_tread_m\":{},\
         \"riser_to_lower_contact_m\":{},\
         \"axle_x_m\":{},\
         \"axle_y_m\":{},\
         \"carrier_rotation_rad\":{},\
         \"chassis_pitch_rad\":{},\
         \"stair_active\":{},\
         \"wheel_centres_m\":[[{},{}],[{},{}],[{},{}]],\
         \"signed_vertical_gaps_m\":[{},{},{}],\
         \"contact_mask\":[{},{},{}],\
         \"contact_count\":{},\
         \"minimum_gap_m\":{}\
         }}}}",
        GENERIC_OWNER,
        MODEL_BOUNDARY,
        SOURCE_RECEIPT,
        pose.name,
        pose.figure,
        SYSTEM_CENTRE_OFFSET_M,
        CLUSTER_RADIUS_M,
        ADJACENT_WHEEL_CENTRE_DISTANCE_M,
        WHEEL_RADIUS_M,
        STAIR_RISE_M,
        STAIR_TREAD_M,
        RISER_TO_LOWER_CONTACT_M,
        pose.axle_x_m,
        pose.axle_y_m,
        pose.carrier_rotation_rad,
        pose.chassis_pitch_rad,
        pose.stair_active,
        result.wheel_centres_m[0][0],
        result.wheel_centres_m[0][1],
        result.wheel_centres_m[1][0],
        result.wheel_centres_m[1][1],
        result.wheel_centres_m[2][0],
        result.wheel_centres_m[2][1],
        result.signed_vertical_gaps_m[0],
        result.signed_vertical_gaps_m[1],
        result.signed_vertical_gaps_m[2],
        result.contact_mask[0],
        result.contact_mask[1],
        result.contact_mask[2],
        result.contact_count,
        result.minimum_gap_m,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn balance_has_one_ground_contact() {
        let output = kamen_cluster_step(1);
        assert!(output.contains("\"state\":\"balance\""));
        assert!(output.contains("\"contact_mask\":[true,false,false]"));
        assert!(output.contains("\"contact_count\":1"));
    }

    #[test]
    fn start_and_climb_touch_successive_horizontal_levels() {
        let start = kamen_cluster_step(2);
        assert!(start.contains("\"source_figure\":\"figure-39b\""));
        assert!(start.contains("\"contact_mask\":[true,true,false]"));

        let climb = kamen_cluster_step(4);
        assert!(climb.contains("\"source_figure\":\"figure-42c\""));
        assert!(climb.contains("\"contact_mask\":[false,true,true]"));
    }

    #[test]
    fn rejects_unknown_state_indices() {
        let output = kamen_cluster_step(6);
        assert!(output.contains("\"code\":\"state-outside-domain\""));
    }
}
