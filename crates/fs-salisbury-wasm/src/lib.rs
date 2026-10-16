//! Browser boundary for the source-bounded US 4,921,293 tendon law.
//!
//! `fs-mbd` owns the generic revolute joints and static torque equations. This
//! crate admits browser inputs and serializes the resulting receipt without
//! adding a dynamic or contact model.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_mbd::salisbury::{SalisburyHandError, SalisburyHandParams, step_salisbury_hand};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Evaluate the source-printed four-tension/three-torque relation.
///
/// Tensions are the four cable ends of one representative digit, in newtons.
/// The returned topology separately reports all twelve cable ends. The
/// `radius_scale_m` argument is the visitor-declared R2 radius in metres and
/// does not claim to reproduce an unpublished historic dimension.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn salisbury_hand_step(
    tension_t1_n: f64,
    tension_t2_n: f64,
    tension_t3_n: f64,
    tension_t4_n: f64,
    radius_scale_m: f64,
    first_idler_fixed: bool,
) -> String {
    let result = match step_salisbury_hand(SalisburyHandParams {
        tension_t1_n,
        tension_t2_n,
        tension_t3_n,
        tension_t4_n,
        radius_scale_m,
        first_idler_fixed,
    }) {
        Ok(result) => result,
        Err(SalisburyHandError::InvalidInput) => {
            return refusal_json(
                "input-outside-domain",
                "cable tensions must be finite and non-negative and the declared radius scale must be finite and positive",
                "Use non-negative newtons and a positive radius in metres",
            );
        }
        Err(SalisburyHandError::Multibody(_)) => {
            return refusal_json(
                "multibody-refusal",
                "the generic revolute-joint owner refused the hand topology",
                "Inspect the fs-mbd Salisbury hand contract",
            );
        }
    };

    let mut output = String::with_capacity(720);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"scalar_joint_coordinates\":{},\
         \"digit_count\":{},\
         \"palm_root_present\":{},\
         \"joint_parent_coordinates\":[{},{},{},{},{},{},{},{},{}],\
         \"cable_end_count\":{},\
         \"axis_1\":[{},{},{}],\
         \"axis_2\":[{},{},{}],\
         \"axis_3\":[{},{},{}],\
         \"tendon_tensions_n\":[{},{},{},{}],\
         \"pulley_radii_m\":[{},{},{}],\
         \"joint_torques_nm\":[{},{},{}],\
         \"claim_1_routing_present\":{},\
         \"claim_2_first_idler_fixed\":{},\
         \"historical_dynamics_available\":{}\
        }}}}",
        result.scalar_joint_coordinates,
        result.digit_count,
        result.palm_root_present,
        result.joint_parent_coordinates[0],
        result.joint_parent_coordinates[1],
        result.joint_parent_coordinates[2],
        result.joint_parent_coordinates[3],
        result.joint_parent_coordinates[4],
        result.joint_parent_coordinates[5],
        result.joint_parent_coordinates[6],
        result.joint_parent_coordinates[7],
        result.joint_parent_coordinates[8],
        result.cable_end_count,
        result.axis_1[0],
        result.axis_1[1],
        result.axis_1[2],
        result.axis_2[0],
        result.axis_2[1],
        result.axis_2[2],
        result.axis_3[0],
        result.axis_3[1],
        result.axis_3[2],
        result.tendon_tensions_n[0],
        result.tendon_tensions_n[1],
        result.tendon_tensions_n[2],
        result.tendon_tensions_n[3],
        result.pulley_radii_m[0],
        result.pulley_radii_m[1],
        result.pulley_radii_m[2],
        result.joint_torques_nm[0],
        result.joint_torques_nm[1],
        result.joint_torques_nm[2],
        result.claim_1_routing_present,
        result.claim_2_first_idler_fixed,
        result.historical_dynamics_available,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_the_generic_topology_and_source_torques() {
        let output = salisbury_hand_step(20.0, 15.0, 5.0, 10.0, 0.01, true);
        assert!(output.contains("\"scalar_joint_coordinates\":9"));
        assert!(output.contains("\"digit_count\":3"));
        assert!(output.contains("\"palm_root_present\":true"));
        assert!(output.contains("\"joint_parent_coordinates\":[-1,0,1,-1,3,4,-1,6,7]"));
        assert!(output.contains("\"cable_end_count\":12"));
        assert!(
            output.contains("\"joint_torques_nm\":[-0.15999999999999998,0.24,0.09999999999999999]")
        );
        assert!(output.contains("\"historical_dynamics_available\":false"));
    }

    #[test]
    fn idler_probe_is_reported_without_rewriting_torque() {
        let fixed = salisbury_hand_step(20.0, 15.0, 5.0, 10.0, 0.01, true);
        let free = salisbury_hand_step(20.0, 15.0, 5.0, 10.0, 0.01, false);
        assert!(fixed.contains("\"claim_2_first_idler_fixed\":true"));
        assert!(free.contains("\"claim_2_first_idler_fixed\":false"));
        assert!(
            free.contains("\"joint_torques_nm\":[-0.15999999999999998,0.24,0.09999999999999999]")
        );
    }

    #[test]
    fn invalid_physical_input_refuses() {
        let negative = salisbury_hand_step(20.0, 15.0, -1.0, 10.0, 0.01, true);
        assert!(negative.contains("\"code\":\"input-outside-domain\""));
        assert!(!negative.contains("\"ok\""));

        let zero_radius = salisbury_hand_step(20.0, 15.0, 5.0, 10.0, 0.0, true);
        assert!(zero_radius.contains("\"code\":\"input-outside-domain\""));
    }
}
