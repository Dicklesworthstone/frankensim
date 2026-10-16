//! Source-bounded generic-joint topology for AMF's US 3,212,649 Versatran.
//!
//! Claim 1 identifies a column that rotates about a vertical axis, a horizontal
//! arm with vertical and horizontal movement, a wrist with rotation about a
//! horizontal axis and swing about a central vertical axis, and operation of a
//! work-manipulating member. The specification also identifies a reciprocating
//! rack/sleeve path for the gripper mechanism. It does not disclose link
//! dimensions, zero transforms, masses, inertias, actuator stroke, pressure,
//! flow, payload, timing, or a program trajectory. This module therefore owns
//! only the six generic scalar channels: five member-motion joints plus one
//! internal work-member operating coordinate.

use crate::articulated::{ArticulatedError, JointKind, JointModel};
use fs_ga::Vec3;

/// Number of scalar channels disclosed by the bounded composition.
pub const VERSATRAN_SCALAR_CHANNELS: usize = 6;
/// Number of channels that move the named column, arm, or wrist members.
pub const VERSATRAN_GEOMETRIC_MOTION_JOINTS: usize = 5;
/// Number of internal work-member operating channels.
pub const VERSATRAN_WORK_MEMBER_OPERATION_CHANNELS: usize = 1;

/// Typed refusal from the source-bounded Versatran composition.
#[derive(Debug, Clone, PartialEq)]
pub enum VersatranTopologyError {
    /// A coordinate was non-finite or a normalized coordinate was outside its
    /// declared presentation interval.
    InvalidInput,
    /// The generic multibody joint owner refused the source composition.
    Multibody(ArticulatedError),
}

impl From<ArticulatedError> for VersatranTopologyError {
    fn from(value: ArticulatedError) -> Self {
        Self::Multibody(value)
    }
}

/// Browser-facing coordinates for one source-bounded Versatran topology query.
///
/// Angles are radians. The three normalized coordinates have no length unit:
/// US 3,212,649 identifies their motions but prints no travel dimensions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VersatranTopologyParams {
    /// Column rotation about the declared central vertical axis, in radians.
    pub column_rotation_rad: f64,
    /// Horizontal arm vertical movement, normalized in `[0, 1]`.
    pub arm_vertical_normalized: f64,
    /// Horizontal arm horizontal movement, normalized in `[0, 1]`.
    pub arm_horizontal_normalized: f64,
    /// Wrist rotation about the declared horizontal axis, in radians.
    pub wrist_rotation_rad: f64,
    /// Wrist swing about the declared central vertical axis, in radians.
    pub wrist_swing_rad: f64,
    /// Reciprocating work-member rack/sleeve coordinate, normalized in `[0, 1]`.
    ///
    /// This is an internal operating channel for the work-manipulating member,
    /// not a sixth rigid-pose degree of freedom.
    pub work_member_rack_normalized: f64,
    /// Whether the patent's manually programmed automatic-operation mode is
    /// selected. This boolean does not execute or synthesize a program.
    pub automatic_program_mode_selected: bool,
}

/// Deterministic generic-joint receipt consumed by browser renderers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VersatranTopologyStep {
    /// Total scalar coordinates in the composed generic-joint topology.
    pub scalar_channels: usize,
    /// Coordinates that move the named column, arm, or wrist members.
    pub geometric_motion_joints: usize,
    /// Internal coordinate that operates the work-manipulating member.
    pub work_member_operation_channels: usize,
    /// Revolute generic joints in the six-channel composition.
    pub revolute_joint_count: usize,
    /// Prismatic generic joints in the six-channel composition.
    pub prismatic_joint_count: usize,
    /// Column rotation axis in the bounded museum convention.
    pub column_rotation_axis: [f64; 3],
    /// Horizontal-arm vertical-motion axis in the bounded museum convention.
    pub arm_vertical_axis: [f64; 3],
    /// Horizontal-arm horizontal-motion axis in the bounded museum convention.
    pub arm_horizontal_axis: [f64; 3],
    /// Wrist horizontal rotation axis in the bounded museum convention.
    pub wrist_rotation_axis: [f64; 3],
    /// Wrist central vertical swing axis in the bounded museum convention.
    pub wrist_swing_axis: [f64; 3],
    /// Internal reciprocal rack/sleeve axis in the bounded museum convention.
    pub work_member_rack_axis: [f64; 3],
    /// Admitted column rotation, in radians.
    pub column_rotation_rad: f64,
    /// Admitted normalized vertical arm movement.
    pub arm_vertical_normalized: f64,
    /// Admitted normalized horizontal arm movement.
    pub arm_horizontal_normalized: f64,
    /// Admitted wrist horizontal-axis rotation, in radians.
    pub wrist_rotation_rad: f64,
    /// Admitted wrist central-vertical-axis swing, in radians.
    pub wrist_swing_rad: f64,
    /// Admitted normalized internal work-member rack coordinate.
    pub work_member_rack_normalized: f64,
    /// Reported automatic-operation mode selection; no program executor is
    /// included in this topology query.
    pub automatic_program_mode_selected: bool,
    /// The historic grant does not provide enough data for a calibrated rigid
    /// geometry or forward-kinematics model.
    pub historical_geometry_available: bool,
    /// The historic grant does not provide enough data for a dynamics or
    /// actuator-control model.
    pub historical_dynamics_available: bool,
}

fn revolute_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().angular;
    [axis.x, axis.y, axis.z]
}

fn prismatic_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().linear;
    [axis.x, axis.y, axis.z]
}

/// Compose the six source-bounded scalar channels from generic joints.
///
/// The returned axes use a named museum normalization convention only:
/// `+Y` is the claim's vertical column direction and `+X` is the horizontal
/// arm/rack direction. They are not dimensions, zero transforms, or an
/// assertion of an original AMF global coordinate frame.
///
/// # Errors
/// Refuses non-finite inputs; normalized arm/rack coordinates outside `[0, 1]`;
/// or an unexpected refusal from the generic multibody joint owner.
pub fn step_versatran_topology(
    params: VersatranTopologyParams,
) -> Result<VersatranTopologyStep, VersatranTopologyError> {
    let values = [
        params.column_rotation_rad,
        params.arm_vertical_normalized,
        params.arm_horizontal_normalized,
        params.wrist_rotation_rad,
        params.wrist_swing_rad,
        params.work_member_rack_normalized,
    ];
    if values.iter().any(|value| !value.is_finite())
        || !(0.0..=1.0).contains(&params.arm_vertical_normalized)
        || !(0.0..=1.0).contains(&params.arm_horizontal_normalized)
        || !(0.0..=1.0).contains(&params.work_member_rack_normalized)
    {
        return Err(VersatranTopologyError::InvalidInput);
    }

    // This is a source-bounded axis convention, not a fabricated link frame:
    // Y is the claim's vertical direction and X is the horizontal arm/rack
    // direction. The wrist's two stated axes are correspondingly X and Y.
    let column_rotation = JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None)?;
    let arm_vertical = JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None)?;
    let arm_horizontal = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    let wrist_rotation = JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None)?;
    let wrist_swing = JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None)?;
    let work_member_rack = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    let geometric_joints = [
        column_rotation,
        arm_vertical,
        arm_horizontal,
        wrist_rotation,
        wrist_swing,
    ];
    let joints = [
        column_rotation,
        arm_vertical,
        arm_horizontal,
        wrist_rotation,
        wrist_swing,
        work_member_rack,
    ];
    debug_assert!(joints.iter().all(|joint| joint.dof_count() == 1));
    debug_assert_eq!(column_rotation.kind(), JointKind::Revolute);
    debug_assert_eq!(arm_vertical.kind(), JointKind::Prismatic);
    debug_assert_eq!(arm_horizontal.kind(), JointKind::Prismatic);
    debug_assert_eq!(wrist_rotation.kind(), JointKind::Revolute);
    debug_assert_eq!(wrist_swing.kind(), JointKind::Revolute);
    debug_assert_eq!(work_member_rack.kind(), JointKind::Prismatic);

    Ok(VersatranTopologyStep {
        scalar_channels: joints.iter().map(|joint| joint.dof_count()).sum(),
        geometric_motion_joints: geometric_joints.iter().map(|joint| joint.dof_count()).sum(),
        work_member_operation_channels: work_member_rack.dof_count(),
        revolute_joint_count: [column_rotation, wrist_rotation, wrist_swing]
            .iter()
            .map(|joint| joint.dof_count())
            .sum(),
        prismatic_joint_count: [arm_vertical, arm_horizontal, work_member_rack]
            .iter()
            .map(|joint| joint.dof_count())
            .sum(),
        column_rotation_axis: revolute_axis(&column_rotation),
        arm_vertical_axis: prismatic_axis(&arm_vertical),
        arm_horizontal_axis: prismatic_axis(&arm_horizontal),
        wrist_rotation_axis: revolute_axis(&wrist_rotation),
        wrist_swing_axis: revolute_axis(&wrist_swing),
        work_member_rack_axis: prismatic_axis(&work_member_rack),
        column_rotation_rad: params.column_rotation_rad,
        arm_vertical_normalized: params.arm_vertical_normalized,
        arm_horizontal_normalized: params.arm_horizontal_normalized,
        wrist_rotation_rad: params.wrist_rotation_rad,
        wrist_swing_rad: params.wrist_swing_rad,
        work_member_rack_normalized: params.work_member_rack_normalized,
        automatic_program_mode_selected: params.automatic_program_mode_selected,
        historical_geometry_available: false,
        historical_dynamics_available: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn baseline() -> VersatranTopologyParams {
        VersatranTopologyParams {
            column_rotation_rad: 0.2,
            arm_vertical_normalized: 0.35,
            arm_horizontal_normalized: 0.65,
            wrist_rotation_rad: -0.4,
            wrist_swing_rad: 0.3,
            work_member_rack_normalized: 0.8,
            automatic_program_mode_selected: true,
        }
    }

    #[test]
    fn composes_three_revolute_and_three_prismatic_source_channels() {
        let step = step_versatran_topology(baseline()).expect("valid source topology");
        assert_eq!(step.scalar_channels, VERSATRAN_SCALAR_CHANNELS);
        assert_eq!(
            step.geometric_motion_joints,
            VERSATRAN_GEOMETRIC_MOTION_JOINTS
        );
        assert_eq!(
            step.work_member_operation_channels,
            VERSATRAN_WORK_MEMBER_OPERATION_CHANNELS
        );
        assert_eq!(step.revolute_joint_count, 3);
        assert_eq!(step.prismatic_joint_count, 3);
        assert_eq!(step.column_rotation_axis, [0.0, 1.0, 0.0]);
        assert_eq!(step.arm_vertical_axis, [0.0, 1.0, 0.0]);
        assert_eq!(step.arm_horizontal_axis, [1.0, 0.0, 0.0]);
        assert_eq!(step.wrist_rotation_axis, [1.0, 0.0, 0.0]);
        assert_eq!(step.wrist_swing_axis, [0.0, 1.0, 0.0]);
        assert_eq!(step.work_member_rack_axis, [1.0, 0.0, 0.0]);
        assert!(!step.historical_geometry_available);
        assert!(!step.historical_dynamics_available);
    }

    #[test]
    fn program_selection_is_reported_without_executing_or_rewriting_coordinates() {
        let selected = step_versatran_topology(baseline()).expect("selected mode");
        let not_selected = step_versatran_topology(VersatranTopologyParams {
            automatic_program_mode_selected: false,
            ..baseline()
        })
        .expect("manual mode");
        assert!(selected.automatic_program_mode_selected);
        assert!(!not_selected.automatic_program_mode_selected);
        assert_eq!(
            selected.work_member_rack_normalized,
            not_selected.work_member_rack_normalized
        );
        assert_eq!(selected.wrist_swing_rad, not_selected.wrist_swing_rad);
    }

    #[test]
    fn rejects_non_finite_or_out_of_domain_presentation_coordinates() {
        assert_eq!(
            step_versatran_topology(VersatranTopologyParams {
                arm_vertical_normalized: 1.01,
                ..baseline()
            }),
            Err(VersatranTopologyError::InvalidInput)
        );
        assert_eq!(
            step_versatran_topology(VersatranTopologyParams {
                work_member_rack_normalized: -0.01,
                ..baseline()
            }),
            Err(VersatranTopologyError::InvalidInput)
        );
        assert_eq!(
            step_versatran_topology(VersatranTopologyParams {
                wrist_rotation_rad: f64::NAN,
                ..baseline()
            }),
            Err(VersatranTopologyError::InvalidInput)
        );
    }
}
