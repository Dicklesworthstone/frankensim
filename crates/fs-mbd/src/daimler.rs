//! Source-bounded multibody composition for Daimler's US 361,931 marine installation.
//!
//! The grant supplies a one-axis longitudinal propeller-shaft motion and the
//! contact topology selected by that motion. It supplies no travel distance,
//! friction coefficient, normal load, shaft speed, cooling flow, or power.
//! This module therefore composes the generic [`crate::articulated::JointModel`]
//! prismatic owner, publishes normalized state and exact topology predicates,
//! and deliberately refuses to manufacture quantitative performance.

use crate::articulated::{ArticulatedError, JointKind, JointModel};
use fs_ga::Vec3;

/// Source-state input for the US 361,931 marine installation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DaimlerMarineParams {
    /// Reader-selected drive state: `-1` astern, `0` neutral, `1` ahead.
    pub shaft_selection: i8,
    /// Whether the optional centrifugal pump `u` is active alongside the
    /// always-present fore/aft pipe path.
    pub cooling_pump_enabled: bool,
}

/// Typed refusal from the source-bounded Daimler composition.
#[derive(Debug, Clone, PartialEq)]
pub enum DaimlerMarineError {
    /// The selector was not one of the three source-facing discrete states.
    InvalidShaftSelection,
    /// The generic multibody owner refused the prismatic joint composition.
    Multibody(ArticulatedError),
}

impl From<ArticulatedError> for DaimlerMarineError {
    fn from(value: ArticulatedError) -> Self {
        Self::Multibody(value)
    }
}

/// Source-bounded contact and cooling state for one reader selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DaimlerMarineResult {
    /// Normalized translation along the vessel's positive sternward axis.
    /// Ahead is negative because the shaft moves toward the motor; the source
    /// supplies no dimensional travel.
    pub shaft_translation_along_axis_normalized: f64,
    /// Normalized axis owned by the generic prismatic joint.
    pub shaft_axis: [f64; 3],
    /// Scalar degrees of freedom in the composed shaft joint.
    pub shaft_joint_dofs: usize,
    /// The continuously one-direction motor-shaft sign.
    pub motor_rotation_sign: i8,
    /// Propeller sign selected by the contact topology (`-1`, `0`, or `1`).
    pub propeller_rotation_sign: i8,
    /// Whether half-couplings `a` and `a²` are in ahead contact.
    pub ahead_coupling_engaged: bool,
    /// Whether reverse disks `e¹/e²` contact `a²/c` for astern drive.
    pub astern_gearing_engaged: bool,
    /// Whether the selected state is neutral with both drive paths open.
    pub neutral: bool,
    /// Whether the source-stated propeller thrust can maintain ahead contact.
    pub thrust_can_maintain_ahead_contact: bool,
    /// Claims 7–9 retain the fore/aft outside-water path in every pump state.
    pub passive_fore_aft_cooling_path_present: bool,
    /// Whether optional centrifugal pump `u` is active.
    pub cooling_pump_active: bool,
}

/// Compose the printed prismatic shaft and mutually exclusive contact states.
///
/// # Errors
/// Refuses a selector outside `-1..=1` or a refusal from the generic
/// multibody joint owner.
pub fn step_daimler_marine(
    params: DaimlerMarineParams,
) -> Result<DaimlerMarineResult, DaimlerMarineError> {
    if !matches!(params.shaft_selection, -1..=1) {
        return Err(DaimlerMarineError::InvalidShaftSelection);
    }

    let shaft_joint = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    debug_assert_eq!(shaft_joint.kind(), JointKind::Prismatic);
    let axis = shaft_joint.motion_subspace().linear;
    let ahead = params.shaft_selection == 1;
    let astern = params.shaft_selection == -1;

    Ok(DaimlerMarineResult {
        shaft_translation_along_axis_normalized: -f64::from(params.shaft_selection),
        shaft_axis: [axis.x, axis.y, axis.z],
        shaft_joint_dofs: shaft_joint.dof_count(),
        motor_rotation_sign: 1,
        propeller_rotation_sign: params.shaft_selection,
        ahead_coupling_engaged: ahead,
        astern_gearing_engaged: astern,
        neutral: params.shaft_selection == 0,
        thrust_can_maintain_ahead_contact: ahead,
        passive_fore_aft_cooling_path_present: true,
        cooling_pump_active: params.cooling_pump_enabled,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ahead_moves_toward_motor_and_closes_only_ahead_contact() {
        let result = step_daimler_marine(DaimlerMarineParams {
            shaft_selection: 1,
            cooling_pump_enabled: false,
        })
        .expect("valid source state");

        assert_eq!(result.shaft_translation_along_axis_normalized, -1.0);
        assert_eq!(result.shaft_axis, [1.0, 0.0, 0.0]);
        assert_eq!(result.shaft_joint_dofs, 1);
        assert_eq!(result.motor_rotation_sign, 1);
        assert_eq!(result.propeller_rotation_sign, 1);
        assert!(result.ahead_coupling_engaged);
        assert!(!result.astern_gearing_engaged);
        assert!(!result.neutral);
        assert!(result.thrust_can_maintain_ahead_contact);
    }

    #[test]
    fn neutral_and_astern_are_mutually_exclusive_source_states() {
        let neutral = step_daimler_marine(DaimlerMarineParams {
            shaft_selection: 0,
            cooling_pump_enabled: false,
        })
        .expect("neutral source state");
        assert!(neutral.neutral);
        assert!(!neutral.ahead_coupling_engaged);
        assert!(!neutral.astern_gearing_engaged);
        assert_eq!(neutral.propeller_rotation_sign, 0);

        let astern = step_daimler_marine(DaimlerMarineParams {
            shaft_selection: -1,
            cooling_pump_enabled: true,
        })
        .expect("astern source state");
        assert_eq!(astern.shaft_translation_along_axis_normalized, 1.0);
        assert!(!astern.ahead_coupling_engaged);
        assert!(astern.astern_gearing_engaged);
        assert!(!astern.neutral);
        assert_eq!(astern.propeller_rotation_sign, -1);
    }

    #[test]
    fn optional_pump_never_erases_the_printed_passive_pipe_path() {
        for cooling_pump_enabled in [false, true] {
            let result = step_daimler_marine(DaimlerMarineParams {
                shaft_selection: 0,
                cooling_pump_enabled,
            })
            .expect("valid cooling state");
            assert!(result.passive_fore_aft_cooling_path_present);
            assert_eq!(result.cooling_pump_active, cooling_pump_enabled);
        }
    }

    #[test]
    fn selector_outside_the_three_source_states_refuses() {
        assert_eq!(
            step_daimler_marine(DaimlerMarineParams {
                shaft_selection: 2,
                cooling_pump_enabled: false,
            }),
            Err(DaimlerMarineError::InvalidShaftSelection)
        );
    }
}
