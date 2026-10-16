//! Source-bounded multibody composition for Elias Howe's US 4,750 sewing machine.
//!
//! The grant prints one main shaft carrying the needle, shuttle, and feed cams;
//! a curved eye-pointed needle fixed to a vibrating arm; shuttle K constrained
//! to trough I and driven by picker-staves J; lifting rod W; and the pinned,
//! rack-holed baster plate H. It supplies no mass, inertia, torque, force,
//! friction, speed, or overall machine dimensions, so this module owns only
//! joint topology, normalized coordinates, sequence predicates, and the two
//! local dimensions actually printed by the specification.

use crate::articulated::{ArticulatedError, JointKind, JointModel};
use fs_ga::Vec3;

/// Approximate distance from the curved needle point to its eye, in inches.
pub const NEEDLE_EYE_OFFSET_IN: f64 = 1.0 / 8.0;
/// Approximate pitch between baster-plate points, in inches.
pub const BASTER_POINT_PITCH_IN: f64 = 3.0 / 4.0;
/// Declared presentation-domain loop-clearance boundary, not a grant dimension.
pub const MINIMUM_LOOP_SLACK_NORMALIZED: f64 = 0.38;

/// Inputs to one deterministic source-order mechanism query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HoweTopologyParams {
    /// Main-shaft coordinate in radians; finite values are wrapped to one turn.
    pub crank_angle_rad: f64,
    /// Normalized displayed loop slack in `[0, 1]`.
    pub loop_slack_normalized: f64,
    /// Whether the Claim 1 needle-and-shuttle combination remains present.
    pub claim_1_interlock_enabled: bool,
}

/// Typed refusal from the source-bounded Howe composition.
#[derive(Debug, Clone, PartialEq)]
pub enum HoweTopologyError {
    /// A coordinate was non-finite or outside its declared interval.
    InvalidInput,
    /// The generic multibody joint owner refused the composition.
    Multibody(ArticulatedError),
}

impl From<ArticulatedError> for HoweTopologyError {
    fn from(value: ArticulatedError) -> Self {
        Self::Multibody(value)
    }
}

/// Generic-joint and source-sequence state consumed by browser renderers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HoweTopologyStep {
    /// Scalar coordinates in the composed joint topology.
    pub scalar_joint_coordinates: usize,
    /// Independent prescribed drives after applying the source linkages.
    pub independent_drive_dofs: usize,
    /// Main-shaft revolute axis.
    pub main_shaft_axis: [f64; 3],
    /// Needle-arm revolute axis.
    pub needle_arm_axis: [f64; 3],
    /// Shuttle prismatic axis in trough I.
    pub shuttle_axis: [f64; 3],
    /// Lifting-rod W prismatic axis.
    pub lifting_rod_axis: [f64; 3],
    /// Baster-plate H feed axis.
    pub baster_feed_axis: [f64; 3],
    /// Wrapped main-shaft angle in radians.
    pub crank_angle_rad: f64,
    /// Normalized depth of the curved needle below its top position.
    pub needle_penetration_normalized: f64,
    /// Needle-arm coordinate in radians.
    pub needle_arm_angle_rad: f64,
    /// Whether the needle is on its retracting half-cycle.
    pub needle_retracting: bool,
    /// Normalized shuttle coordinate along trough I.
    pub shuttle_travel_normalized: f64,
    /// Instantaneous normalized upper-thread loop opening.
    pub loop_open_fraction: f64,
    /// Whether the displayed loop clears the declared shuttle section.
    pub loop_open: bool,
    /// Whether shuttle K currently passes through that loop.
    pub shuttle_passes_loop: bool,
    /// Counterfactual normalized track offset when Claim 1 is disabled.
    pub shuttle_track_offset_normalized: f64,
    /// Normalized left picker-stave engagement.
    pub picker_left_normalized: f64,
    /// Normalized right picker-stave engagement.
    pub picker_right_normalized: f64,
    /// Normalized lifting-rod W coordinate.
    pub lifting_rod_normalized: f64,
    /// Normalized intermittent baster-plate advance.
    pub feed_advance_fraction: f64,
    /// Whether the thread clamp is closed at the cycle boundary.
    pub thread_clamp_engaged: bool,
    /// Whether the selected topology and slack admit Claim 1's combination.
    pub claim_1_interlock_satisfied: bool,
    /// Stable source-order phase label.
    pub cycle_phase: &'static str,
    /// Printed approximate needle-eye offset in inches.
    pub needle_eye_offset_in: f64,
    /// Printed approximate baster-point pitch in inches.
    pub baster_point_pitch_in: f64,
}

fn smoothstep(value: f64) -> f64 {
    let t = value.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

fn revolute_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().angular;
    [axis.x, axis.y, axis.z]
}

fn prismatic_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().linear;
    [axis.x, axis.y, axis.z]
}

/// Compose the source-printed joints and evaluate one shaft-driven cycle.
///
/// # Errors
/// Refuses non-finite angle/slack, slack outside `[0, 1]`, or a generic
/// multibody-joint refusal.
pub fn step_howe_topology(
    params: HoweTopologyParams,
) -> Result<HoweTopologyStep, HoweTopologyError> {
    if !params.crank_angle_rad.is_finite()
        || !params.loop_slack_normalized.is_finite()
        || !(0.0..=1.0).contains(&params.loop_slack_normalized)
    {
        return Err(HoweTopologyError::InvalidInput);
    }

    let main_shaft = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let needle_arm = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let shuttle = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    let left_picker = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let right_picker = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None)?;
    let lifting_rod = JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None)?;
    let baster_feed = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None)?;
    let joints = [
        main_shaft,
        needle_arm,
        shuttle,
        left_picker,
        right_picker,
        lifting_rod,
        baster_feed,
    ];
    debug_assert_eq!(main_shaft.kind(), JointKind::Revolute);
    debug_assert_eq!(needle_arm.kind(), JointKind::Revolute);
    debug_assert_eq!(shuttle.kind(), JointKind::Prismatic);
    debug_assert_eq!(lifting_rod.kind(), JointKind::Prismatic);
    debug_assert_eq!(baster_feed.kind(), JointKind::Prismatic);

    let angle = params.crank_angle_rad.rem_euclid(core::f64::consts::TAU);
    let angle_deg = angle.to_degrees();
    let needle_penetration = 0.5 * (1.0 - angle.cos());
    let needle_retracting = angle > core::f64::consts::PI;
    let shuttle_travel = -angle.cos();
    let loop_envelope = if (180.0..330.0).contains(&angle_deg) {
        (((angle_deg - 180.0) / 150.0) * core::f64::consts::PI)
            .sin()
            .max(0.0)
    } else {
        0.0
    };
    let loop_open_fraction = loop_envelope * params.loop_slack_normalized;
    let loop_open =
        params.claim_1_interlock_enabled && loop_open_fraction >= MINIMUM_LOOP_SLACK_NORMALIZED;
    let shuttle_at_needle_plane =
        (210.0..320.0).contains(&angle_deg) && shuttle_travel.abs() < 0.22;
    let feed_advance_fraction = smoothstep((angle_deg - 315.0) / 45.0);
    let cycle_phase = if angle_deg < 180.0 {
        "penetrate"
    } else if angle_deg < 235.0 {
        "retract-and-open-loop"
    } else if angle_deg < 315.0 {
        "shuttle-pass"
    } else {
        "feed"
    };

    Ok(HoweTopologyStep {
        scalar_joint_coordinates: joints.iter().map(|joint| joint.dof_count()).sum(),
        independent_drive_dofs: 1,
        main_shaft_axis: revolute_axis(&main_shaft),
        needle_arm_axis: revolute_axis(&needle_arm),
        shuttle_axis: prismatic_axis(&shuttle),
        lifting_rod_axis: prismatic_axis(&lifting_rod),
        baster_feed_axis: prismatic_axis(&baster_feed),
        crank_angle_rad: angle,
        needle_penetration_normalized: needle_penetration,
        needle_arm_angle_rad: 0.12 - 0.24 * needle_penetration,
        needle_retracting,
        shuttle_travel_normalized: shuttle_travel,
        loop_open_fraction,
        loop_open,
        shuttle_passes_loop: loop_open && shuttle_at_needle_plane,
        shuttle_track_offset_normalized: if params.claim_1_interlock_enabled {
            0.0
        } else {
            0.55
        },
        picker_left_normalized: (-shuttle_travel).max(0.0),
        picker_right_normalized: shuttle_travel.max(0.0),
        lifting_rod_normalized: loop_open_fraction,
        feed_advance_fraction,
        thread_clamp_engaged: angle_deg >= 320.0 || angle_deg <= 35.0,
        claim_1_interlock_satisfied: params.claim_1_interlock_enabled
            && params.loop_slack_normalized >= MINIMUM_LOOP_SLACK_NORMALIZED,
        cycle_phase,
        needle_eye_offset_in: NEEDLE_EYE_OFFSET_IN,
        baster_point_pitch_in: BASTER_POINT_PITCH_IN,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composes_one_driven_source_topology_from_generic_joints() {
        let step = step_howe_topology(HoweTopologyParams {
            crank_angle_rad: 1.5 * core::f64::consts::PI,
            loop_slack_normalized: 0.65,
            claim_1_interlock_enabled: true,
        })
        .expect("valid source topology");
        assert_eq!(step.scalar_joint_coordinates, 7);
        assert_eq!(step.independent_drive_dofs, 1);
        assert_eq!(step.main_shaft_axis, [0.0, 0.0, 1.0]);
        assert_eq!(step.shuttle_axis, [1.0, 0.0, 0.0]);
        assert_eq!(step.lifting_rod_axis, [0.0, 1.0, 0.0]);
        assert!(step.needle_retracting);
        assert!(step.loop_open);
        assert!(step.shuttle_passes_loop);
        assert_eq!(step.cycle_phase, "shuttle-pass");
    }

    #[test]
    fn insufficient_slack_and_removed_claim_never_report_an_interlock() {
        let low_slack = step_howe_topology(HoweTopologyParams {
            crank_angle_rad: 1.5 * core::f64::consts::PI,
            loop_slack_normalized: 0.2,
            claim_1_interlock_enabled: true,
        })
        .unwrap();
        assert!(!low_slack.loop_open);
        assert!(!low_slack.shuttle_passes_loop);

        let removed = step_howe_topology(HoweTopologyParams {
            crank_angle_rad: 1.5 * core::f64::consts::PI,
            loop_slack_normalized: 0.65,
            claim_1_interlock_enabled: false,
        })
        .unwrap();
        assert!(!removed.loop_open);
        assert!(!removed.shuttle_passes_loop);
        assert_eq!(removed.shuttle_track_offset_normalized, 0.55);
    }

    #[test]
    fn source_dimensions_are_reported_without_inventing_machine_dimensions() {
        let step = step_howe_topology(HoweTopologyParams {
            crank_angle_rad: 0.0,
            loop_slack_normalized: 0.65,
            claim_1_interlock_enabled: true,
        })
        .unwrap();
        assert_eq!(step.needle_eye_offset_in, 0.125);
        assert_eq!(step.baster_point_pitch_in, 0.75);
        assert_eq!(step.needle_penetration_normalized, 0.0);
        assert_eq!(step.shuttle_travel_normalized, -1.0);
    }

    #[test]
    fn invalid_slack_and_non_finite_angle_refuse() {
        let baseline = HoweTopologyParams {
            crank_angle_rad: 0.0,
            loop_slack_normalized: 0.65,
            claim_1_interlock_enabled: true,
        };
        assert_eq!(
            step_howe_topology(HoweTopologyParams {
                loop_slack_normalized: 1.01,
                ..baseline
            }),
            Err(HoweTopologyError::InvalidInput)
        );
        assert_eq!(
            step_howe_topology(HoweTopologyParams {
                crank_angle_rad: f64::NAN,
                ..baseline
            }),
            Err(HoweTopologyError::InvalidInput)
        );
    }
}
