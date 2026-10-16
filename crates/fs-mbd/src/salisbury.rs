//! Source-bounded tendon topology and torque map for US 4,921,293.
//!
//! The grant prints three joint-torque equations for the Figure 3 routing and
//! identifies three revolute joints per digit, four cable ends per digit, and
//! three digits. It does not print pulley dimensions, link geometry, inertia,
//! damping, motor limits, friction, or contact properties. This composition
//! therefore owns the connected generic joint topology and exact static torque
//! law only; all radii and tensions are caller-declared SI study inputs.

use crate::articulated::{ArticulatedError, JointModel};
use fs_ga::Vec3;

/// Failure to admit the source-bounded hand composition.
#[derive(Clone, Debug, PartialEq)]
pub enum SalisburyHandError {
    /// A cable tension was negative/non-finite or a radius scale was not positive and finite.
    InvalidInput,
    /// The generic revolute-joint owner refused a source axis.
    Multibody(ArticulatedError),
}

impl From<ArticulatedError> for SalisburyHandError {
    fn from(value: ArticulatedError) -> Self {
        Self::Multibody(value)
    }
}

/// Caller-declared SI inputs to the Figure 3 tendon law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SalisburyHandParams {
    /// Cable-end tension T1 in newtons.
    pub tension_t1_n: f64,
    /// Cable-end tension T2 in newtons.
    pub tension_t2_n: f64,
    /// Cable-end tension T3 in newtons.
    pub tension_t3_n: f64,
    /// Cable-end tension T4 in newtons.
    pub tension_t4_n: f64,
    /// Visitor-declared R2 study radius in metres.
    pub radius_scale_m: f64,
    /// Whether the first idler is held for the Claim 2 teaching probe.
    pub first_idler_fixed: bool,
}

/// One deterministic source-topology and torque-law receipt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SalisburyHandStep {
    /// Three articulated digits are attached to the palm.
    pub digit_count: usize,
    /// The palm is the common fixed root for all three serial digit chains.
    pub palm_root_present: bool,
    /// Three revolute coordinates per digit across three digits.
    pub scalar_joint_coordinates: usize,
    /// Parent coordinate for each joint; `-1` denotes the common palm root.
    pub joint_parent_coordinates: [i8; 9],
    /// Four cable ends per digit across three digits.
    pub cable_end_count: usize,
    /// Axis 1 of each digit in the museum model frame.
    pub axis_1: [f64; 3],
    /// Axis 2 of each digit, perpendicular to Axis 1.
    pub axis_2: [f64; 3],
    /// Axis 3 of each digit, parallel to Axis 2 in the preferred embodiment.
    pub axis_3: [f64; 3],
    /// One representative digit's admitted tensions `[T1, T2, T3, T4]` in newtons.
    pub tendon_tensions_n: [f64; 4],
    /// Declared illustrative radii `[R1, R2, R3]` in metres.
    pub pulley_radii_m: [f64; 3],
    /// Printed Figure 3 outputs `[Torque1, Torque2, Torque3]` in newton metres.
    pub joint_torques_nm: [f64; 3],
    /// The complete four-cable/three-joint source topology is present.
    pub claim_1_routing_present: bool,
    /// The first idler is fixed for the Claim 2 teaching predicate.
    pub claim_2_first_idler_fixed: bool,
    /// Historic dynamic performance remains unavailable from this grant.
    pub historical_dynamics_available: bool,
}

fn revolute_axis(joint: &JointModel) -> [f64; 3] {
    let axis = joint.motion_subspace().angular;
    [axis.x, axis.y, axis.z]
}

fn valid_tension(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

/// Compose the nine generic revolute joints and evaluate the printed SI law.
///
/// The caller declares `radius_scale_m = R2`. The museum study uses
/// `R1 = 1.2 R2` and `R3 = 1.4 R2` solely to preserve Figure 3's depicted
/// ordering `R3 > R1 > R2`; those ratios are not claimed historical dimensions.
pub fn step_salisbury_hand(
    params: SalisburyHandParams,
) -> Result<SalisburyHandStep, SalisburyHandError> {
    // Figure 3 prints this law for one digit. The topology receipt separately
    // proves three four-cable digit chains; these four values are not the full
    // hand's twelve independent cable tensions.
    let tensions = [
        params.tension_t1_n,
        params.tension_t2_n,
        params.tension_t3_n,
        params.tension_t4_n,
    ];
    if tensions.into_iter().any(|value| !valid_tension(value))
        || !params.radius_scale_m.is_finite()
        || params.radius_scale_m <= 0.0
    {
        return Err(SalisburyHandError::InvalidInput);
    }

    let axis_1_joint = JointModel::revolute(Vec3::new(0.0, 1.0, 0.0), None)?;
    let axis_2_joint = JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None)?;
    let axis_3_joint = JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None)?;
    let digit_joints = [axis_1_joint, axis_2_joint, axis_3_joint];
    let hand_joints = [digit_joints, digit_joints, digit_joints];
    let scalar_joint_coordinates = hand_joints
        .iter()
        .flatten()
        .map(|joint| joint.dof_count())
        .sum::<usize>();
    // Each first axis is palm-anchored; the other two are serial children.
    let joint_parent_coordinates = [-1, 0, 1, -1, 3, 4, -1, 6, 7];

    let r2 = params.radius_scale_m;
    let r1 = 1.2 * r2;
    let r3 = 1.4 * r2;
    let [t1, t2, t3, t4] = tensions;

    // Exact signs and products printed in the Figure 3 description.
    let torque_1 = -t1 * r1 + t2 * r2 + t3 * r2 - t4 * r1;
    let torque_2 = t1 * r3 + t2 * r2 - t3 * r2 - t4 * r3;
    let torque_3 = t2 * r2 - t3 * r2;
    let joint_torques_nm = [torque_1, torque_2, torque_3];
    if joint_torques_nm.into_iter().any(|value| !value.is_finite()) {
        return Err(SalisburyHandError::InvalidInput);
    }

    Ok(SalisburyHandStep {
        digit_count: hand_joints.len(),
        palm_root_present: true,
        scalar_joint_coordinates,
        joint_parent_coordinates,
        cable_end_count: 12,
        axis_1: revolute_axis(&axis_1_joint),
        axis_2: revolute_axis(&axis_2_joint),
        axis_3: revolute_axis(&axis_3_joint),
        tendon_tensions_n: tensions,
        pulley_radii_m: [r1, r2, r3],
        joint_torques_nm,
        claim_1_routing_present: true,
        claim_2_first_idler_fixed: params.first_idler_fixed,
        historical_dynamics_available: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn example() -> SalisburyHandParams {
        SalisburyHandParams {
            tension_t1_n: 20.0,
            tension_t2_n: 15.0,
            tension_t3_n: 5.0,
            tension_t4_n: 10.0,
            radius_scale_m: 0.01,
            first_idler_fixed: true,
        }
    }

    #[test]
    fn composes_three_connected_revolute_axes_for_each_of_three_digits() {
        let step = step_salisbury_hand(example()).expect("valid source topology");
        assert_eq!(step.scalar_joint_coordinates, 9);
        assert_eq!(step.cable_end_count, 12);
        assert_eq!(step.digit_count, 3);
        assert!(step.palm_root_present);
        assert_eq!(
            step.joint_parent_coordinates,
            [-1, 0, 1, -1, 3, 4, -1, 6, 7]
        );
        assert_eq!(step.axis_1, [0.0, 1.0, 0.0]);
        assert_eq!(step.axis_2, [1.0, 0.0, 0.0]);
        assert_eq!(step.axis_3, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn evaluates_the_three_source_printed_torque_equations() {
        let step = step_salisbury_hand(example()).expect("valid torque study");
        assert!((step.pulley_radii_m[0] - 0.012).abs() < 1.0e-15);
        assert!((step.pulley_radii_m[1] - 0.010).abs() < 1.0e-15);
        assert!((step.pulley_radii_m[2] - 0.014).abs() < 1.0e-15);
        assert!((step.joint_torques_nm[0] - -0.16).abs() < 1.0e-14);
        assert!((step.joint_torques_nm[1] - 0.24).abs() < 1.0e-14);
        assert!((step.joint_torques_nm[2] - 0.10).abs() < 1.0e-14);
    }

    #[test]
    fn refuses_nonphysical_tension_and_radius_inputs() {
        let mut negative = example();
        negative.tension_t3_n = -1.0;
        assert_eq!(
            step_salisbury_hand(negative),
            Err(SalisburyHandError::InvalidInput)
        );

        let mut zero_radius = example();
        zero_radius.radius_scale_m = 0.0;
        assert_eq!(
            step_salisbury_hand(zero_radius),
            Err(SalisburyHandError::InvalidInput)
        );
    }

    #[test]
    fn idler_probe_does_not_rewrite_the_printed_torque_law() {
        let fixed = step_salisbury_hand(example()).expect("fixed-idler study");
        let mut free_params = example();
        free_params.first_idler_fixed = false;
        let free = step_salisbury_hand(free_params).expect("free-idler study");
        assert!(fixed.claim_2_first_idler_fixed);
        assert!(!free.claim_2_first_idler_fixed);
        assert_eq!(fixed.joint_torques_nm, free.joint_torques_nm);
        assert!(!free.historical_dynamics_available);
    }
}
