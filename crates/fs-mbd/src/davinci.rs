//! Source-bounded joint topology for the US 6,331,181 tool manipulator.
//!
//! Figures 2 and 2A print a supported manipulator, insertion along the tool
//! axis, rotations about intersecting axes, a releasable tool holder, and a
//! distal tool. They do not print link dimensions, motor data, or a universal
//! commercial degree-of-freedom count. This module therefore composes only the
//! generic joint kinds and normalized presentation coordinate used at the
//! museum boundary.

use crate::articulated::{ArticulatedError, JointKind, JointModel};
use fs_ga::Vec3;

/// Refusal channel for the source-bounded manipulator topology.
#[derive(Debug, Clone, PartialEq)]
pub enum DaVinciTopologyError {
    /// A presentation coordinate was non-finite or outside its declared
    /// normalized interval.
    InvalidInput,
    /// The generic multibody joint owner refused the composition.
    Multibody(ArticulatedError),
}

impl From<ArticulatedError> for DaVinciTopologyError {
    fn from(value: ArticulatedError) -> Self {
        Self::Multibody(value)
    }
}

/// Normalized browser inputs. Angles are radians; insertion has no length unit
/// because the grant supplies no dimensioned travel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DaVinciTopologyParams {
    /// Rotation about the upright support/tool plane.
    pub base_yaw_rad: f64,
    /// Rotation of the supported tool carriage.
    pub carriage_pitch_rad: f64,
    /// Rotation of the distal pitch joint.
    pub distal_pitch_rad: f64,
    /// Rotation of the distal yaw joint.
    pub distal_yaw_rad: f64,
    /// Rotation about the tool shaft.
    pub tool_roll_rad: f64,
    /// Normalized insertion coordinate in `[-1, 1]`.
    pub insertion_normalized: f64,
    /// Whether the claimed compatibility identifier is present.
    pub compatibility_identifier_present: bool,
}

/// Generic-joint receipt consumed by the browser model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DaVinciTopologyStep {
    /// Number of composed scalar joints.
    pub joint_dofs: usize,
    /// Unit axis for support yaw.
    pub base_yaw_axis: [f64; 3],
    /// Unit axis for carriage pitch.
    pub carriage_pitch_axis: [f64; 3],
    /// Unit axis for insertion.
    pub insertion_axis: [f64; 3],
    /// Unit axis for distal pitch.
    pub distal_pitch_axis: [f64; 3],
    /// Unit axis for distal yaw.
    pub distal_yaw_axis: [f64; 3],
    /// Unit axis for tool roll.
    pub tool_roll_axis: [f64; 3],
    /// Admitted support yaw coordinate.
    pub base_yaw_rad: f64,
    /// Admitted carriage pitch coordinate.
    pub carriage_pitch_rad: f64,
    /// Admitted distal pitch coordinate.
    pub distal_pitch_rad: f64,
    /// Admitted distal yaw coordinate.
    pub distal_yaw_rad: f64,
    /// Admitted roll coordinate.
    pub tool_roll_rad: f64,
    /// Admitted normalized insertion coordinate.
    pub insertion_normalized: f64,
    /// Source-facing compatibility predicate; this is not a fabricated motor
    /// interlock or safety certification.
    pub compatibility_identifier_present: bool,
}

/// Compose the generic revolute/prismatic joints printed schematically in
/// Figs. 2 and 2A and return their stable axes and admitted coordinates.
///
/// # Errors
/// Refuses non-finite inputs, an insertion coordinate outside `[-1, 1]`, or
/// any refusal from the generic joint owner.
pub fn step_davinci_topology(
    params: DaVinciTopologyParams,
) -> Result<DaVinciTopologyStep, DaVinciTopologyError> {
    let values = [
        params.base_yaw_rad,
        params.carriage_pitch_rad,
        params.distal_pitch_rad,
        params.distal_yaw_rad,
        params.tool_roll_rad,
        params.insertion_normalized,
    ];
    if values.iter().any(|value| !value.is_finite())
        || !(-1.0..=1.0).contains(&params.insertion_normalized)
    {
        return Err(DaVinciTopologyError::InvalidInput);
    }

    let base_yaw = JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None)?;
    let carriage_pitch = JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None)?;
    let insertion = JointModel::prismatic(Vec3::new(0.0, -1.0, 0.0), None)?;
    let distal_pitch = JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None)?;
    let distal_yaw = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let tool_roll = JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None)?;
    let joints = [
        base_yaw,
        carriage_pitch,
        insertion,
        distal_pitch,
        distal_yaw,
        tool_roll,
    ];
    debug_assert!(joints.iter().all(|joint| joint.dof_count() == 1));
    debug_assert_eq!(insertion.kind(), JointKind::Prismatic);
    debug_assert!(
        joints
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != 2)
            .all(|(_, joint)| joint.kind() == JointKind::Revolute)
    );

    Ok(DaVinciTopologyStep {
        joint_dofs: joints.iter().map(|joint| joint.dof_count()).sum(),
        base_yaw_axis: [0.0, 1.0, 0.0],
        carriage_pitch_axis: [1.0, 0.0, 0.0],
        insertion_axis: [0.0, -1.0, 0.0],
        distal_pitch_axis: [1.0, 0.0, 0.0],
        distal_yaw_axis: [0.0, 0.0, 1.0],
        tool_roll_axis: [0.0, 1.0, 0.0],
        base_yaw_rad: params.base_yaw_rad,
        carriage_pitch_rad: params.carriage_pitch_rad,
        distal_pitch_rad: params.distal_pitch_rad,
        distal_yaw_rad: params.distal_yaw_rad,
        tool_roll_rad: params.tool_roll_rad,
        insertion_normalized: params.insertion_normalized,
        compatibility_identifier_present: params.compatibility_identifier_present,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composes_five_revolute_axes_and_one_normalized_insertion_axis() {
        let step = step_davinci_topology(DaVinciTopologyParams {
            base_yaw_rad: 0.2,
            carriage_pitch_rad: -0.3,
            distal_pitch_rad: 0.4,
            distal_yaw_rad: -0.1,
            tool_roll_rad: 1.2,
            insertion_normalized: 0.25,
            compatibility_identifier_present: true,
        })
        .unwrap();
        assert_eq!(step.joint_dofs, 6);
        assert_eq!(step.base_yaw_axis, [0.0, 1.0, 0.0]);
        assert_eq!(step.carriage_pitch_axis, [1.0, 0.0, 0.0]);
        assert_eq!(step.insertion_axis, [0.0, -1.0, 0.0]);
        assert_eq!(step.distal_pitch_axis, [1.0, 0.0, 0.0]);
        assert_eq!(step.distal_yaw_axis, [0.0, 0.0, 1.0]);
        assert_eq!(step.tool_roll_axis, [0.0, 1.0, 0.0]);
        assert!(step.compatibility_identifier_present);
    }

    #[test]
    fn compatibility_is_a_reported_identifier_predicate_not_a_motion_rewrite() {
        let present = step_davinci_topology(DaVinciTopologyParams {
            base_yaw_rad: 0.2,
            carriage_pitch_rad: -0.3,
            distal_pitch_rad: 0.4,
            distal_yaw_rad: -0.1,
            tool_roll_rad: 1.2,
            insertion_normalized: -0.4,
            compatibility_identifier_present: true,
        })
        .unwrap();
        let absent = step_davinci_topology(DaVinciTopologyParams {
            compatibility_identifier_present: false,
            ..DaVinciTopologyParams {
                base_yaw_rad: 0.2,
                carriage_pitch_rad: -0.3,
                distal_pitch_rad: 0.4,
                distal_yaw_rad: -0.1,
                tool_roll_rad: 1.2,
                insertion_normalized: -0.4,
                compatibility_identifier_present: true,
            }
        })
        .unwrap();
        assert_eq!(present.base_yaw_rad, absent.base_yaw_rad);
        assert_eq!(present.insertion_normalized, absent.insertion_normalized);
        assert!(!absent.compatibility_identifier_present);
    }

    #[test]
    fn rejects_non_finite_or_out_of_domain_presentation_coordinates() {
        let baseline = DaVinciTopologyParams {
            base_yaw_rad: 0.0,
            carriage_pitch_rad: 0.0,
            distal_pitch_rad: 0.0,
            distal_yaw_rad: 0.0,
            tool_roll_rad: 0.0,
            insertion_normalized: 0.0,
            compatibility_identifier_present: true,
        };
        assert_eq!(
            step_davinci_topology(DaVinciTopologyParams {
                insertion_normalized: 1.01,
                ..baseline
            }),
            Err(DaVinciTopologyError::InvalidInput)
        );
        assert_eq!(
            step_davinci_topology(DaVinciTopologyParams {
                distal_yaw_rad: f64::NAN,
                ..baseline
            }),
            Err(DaVinciTopologyError::InvalidInput)
        );
    }
}
