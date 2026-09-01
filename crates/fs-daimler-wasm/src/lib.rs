//! Browser boundary for the source-bounded US 361,931 marine apparatus.
//!
//! The generic multibody law remains in `fs-mbd`; this crate performs only
//! browser input admission and stable envelope serialization.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_mbd::daimler::{DaimlerMarineError, DaimlerMarineParams, step_daimler_marine};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Step the normalized source topology for Daimler's US 361,931 installation.
///
/// `shaft_selection` is exactly `-1` astern, `0` neutral, or `1` ahead.
/// Distances are normalized display coordinates because the grant supplies no
/// physical shaft travel.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
pub fn daimler_marine_step(shaft_selection: i32, cooling_pump_enabled: bool) -> String {
    if !matches!(shaft_selection, -1..=1) {
        return refusal_json(
            "input-outside-domain",
            "shaft_selection must be -1 (astern), 0 (neutral), or 1 (ahead)",
            "Choose one of the three source-facing drive states",
        );
    }

    let result = match step_daimler_marine(DaimlerMarineParams {
        shaft_selection: shaft_selection as i8,
        cooling_pump_enabled,
    }) {
        Ok(result) => result,
        Err(DaimlerMarineError::InvalidShaftSelection) => {
            return refusal_json(
                "input-outside-domain",
                "the multibody owner refused the shaft selector",
                "Choose -1, 0, or 1",
            );
        }
        Err(DaimlerMarineError::Multibody(_)) => {
            return refusal_json(
                "multibody-refusal",
                "the generic prismatic-joint owner refused the composition",
                "Inspect the fs-mbd Daimler composition contract",
            );
        }
    };

    let mut output = String::with_capacity(640);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"shaft_translation_along_axis_normalized\":{},\
         \"shaft_axis\":[{},{},{}],\
         \"shaft_joint_dofs\":{},\
         \"motor_rotation_sign\":{},\
         \"propeller_rotation_sign\":{},\
         \"ahead_coupling_engaged\":{},\
         \"astern_gearing_engaged\":{},\
         \"neutral\":{},\
         \"thrust_can_maintain_ahead_contact\":{},\
         \"passive_fore_aft_cooling_path_present\":{},\
         \"cooling_pump_active\":{}\
         }}}}",
        result.shaft_translation_along_axis_normalized,
        result.shaft_axis[0],
        result.shaft_axis[1],
        result.shaft_axis[2],
        result.shaft_joint_dofs,
        result.motor_rotation_sign,
        result.propeller_rotation_sign,
        result.ahead_coupling_engaged,
        result.astern_gearing_engaged,
        result.neutral,
        result.thrust_can_maintain_ahead_contact,
        result.passive_fore_aft_cooling_path_present,
        result.cooling_pump_active,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ahead_and_astern_emit_opposite_exclusive_topologies() {
        let ahead = daimler_marine_step(1, false);
        assert!(ahead.contains("\"shaft_translation_along_axis_normalized\":-1"));
        assert!(ahead.contains("\"ahead_coupling_engaged\":true"));
        assert!(ahead.contains("\"astern_gearing_engaged\":false"));
        assert!(ahead.contains("\"propeller_rotation_sign\":1"));

        let astern = daimler_marine_step(-1, false);
        assert!(astern.contains("\"shaft_translation_along_axis_normalized\":1"));
        assert!(astern.contains("\"ahead_coupling_engaged\":false"));
        assert!(astern.contains("\"astern_gearing_engaged\":true"));
        assert!(astern.contains("\"propeller_rotation_sign\":-1"));
    }

    #[test]
    fn neutral_opens_both_drive_paths() {
        let neutral = daimler_marine_step(0, false);
        assert!(neutral.contains("\"neutral\":true"));
        assert!(neutral.contains("\"ahead_coupling_engaged\":false"));
        assert!(neutral.contains("\"astern_gearing_engaged\":false"));
        assert!(neutral.contains("\"motor_rotation_sign\":1"));
        assert!(neutral.contains("\"propeller_rotation_sign\":0"));
    }

    #[test]
    fn pump_is_additive_to_the_passive_path() {
        let pumped = daimler_marine_step(0, true);
        assert!(pumped.contains("\"passive_fore_aft_cooling_path_present\":true"));
        assert!(pumped.contains("\"cooling_pump_active\":true"));
    }

    #[test]
    fn invalid_selector_refuses() {
        let refusal = daimler_marine_step(2, false);
        assert!(refusal.contains("\"code\":\"input-outside-domain\""));
        assert!(!refusal.contains("\"ok\""));
    }
}
