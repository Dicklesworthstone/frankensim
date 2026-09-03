//! Rigid planar kinematics for an equal-radius three-wheel cluster on level
//! ground or an ideal two-riser stair profile.
//!
//! This module owns geometry and vertical contact-gap evaluation only. It does
//! not infer normal force, friction, tire compliance, motor torque, controller
//! response, impact, or side contact with a stair riser.

use core::f64::consts::{PI, TAU};
use core::fmt;

/// Number of equally spaced wheels in the supported cluster.
pub const TRI_WHEEL_COUNT: usize = 3;

/// Complete input for one rigid tri-wheel support/contact evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriWheelStairInput {
    /// Distance from the carrier axis to every wheel axis [m].
    pub cluster_radius_m: f64,
    /// Outside rolling radius of every wheel [m].
    pub wheel_radius_m: f64,
    /// Carrier-axis horizontal coordinate [m].
    pub axle_x_m: f64,
    /// Carrier-axis vertical coordinate above the ground datum [m].
    pub axle_y_m: f64,
    /// Counter-clockwise carrier rotation from wheel A vertically downward [rad].
    pub carrier_rotation_rad: f64,
    /// Height of each ideal stair riser [m].
    pub stair_rise_m: f64,
    /// Depth of each ideal stair tread [m].
    pub stair_tread_m: f64,
    /// Whether to evaluate the two-riser stair rather than level ground.
    pub stair_active: bool,
    /// Absolute vertical gap admitted as touching [m].
    pub contact_tolerance_m: f64,
}

/// One admitted rigid support configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriWheelStairResult {
    /// Wheel-axis coordinates, ordered A/B/C [m].
    pub wheel_centres_m: [[f64; 2]; TRI_WHEEL_COUNT],
    /// Wheel-bottom minus horizontal support height, ordered A/B/C [m].
    pub signed_vertical_gaps_m: [f64; TRI_WHEEL_COUNT],
    /// Whether each vertical gap lies within the declared tolerance.
    pub contact_mask: [bool; TRI_WHEEL_COUNT],
    /// Number of admitted horizontal support contacts.
    pub contact_count: u8,
    /// Smallest signed vertical gap [m].
    pub minimum_gap_m: f64,
}

/// Typed refusal from the rigid tri-wheel support evaluator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TriWheelStairError {
    /// A named input is non-finite or outside its physical domain.
    InvalidInput(&'static str),
    /// A wheel penetrates a horizontal support beyond the declared tolerance.
    PenetratingSupport {
        /// Zero-based A/B/C wheel index.
        wheel_index: usize,
        /// Signed vertical gap [m].
        gap_m: f64,
    },
    /// No wheel touches a horizontal support within the declared tolerance.
    Unsupported {
        /// Smallest absolute vertical gap [m].
        nearest_gap_m: f64,
    },
}

impl fmt::Display for TriWheelStairError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(field) => write!(formatter, "invalid tri-wheel input: {field}"),
            Self::PenetratingSupport { wheel_index, gap_m } => write!(
                formatter,
                "tri-wheel {wheel_index} penetrates its horizontal support by signed gap {gap_m:.9e} m"
            ),
            Self::Unsupported { nearest_gap_m } => write!(
                formatter,
                "tri-wheel cluster has no horizontal support contact; nearest gap is {nearest_gap_m:.9e} m"
            ),
        }
    }
}

impl std::error::Error for TriWheelStairError {}

/// Height of the ideal horizontal support beneath one wheel axis [m].
///
/// With stairs active, the profile has level ground for `x < 0`, the first
/// tread over `[0, tread)`, and a second-tread plateau for `x >= tread`.
#[must_use]
pub fn horizontal_support_height_m(
    x_m: f64,
    stair_rise_m: f64,
    stair_tread_m: f64,
    stair_active: bool,
) -> f64 {
    if !stair_active || x_m < 0.0 {
        0.0
    } else if x_m < stair_tread_m {
        stair_rise_m
    } else {
        2.0 * stair_rise_m
    }
}

/// Evaluate wheel-axis coordinates and horizontal-support gaps for one rigid
/// tri-wheel cluster.
///
/// Wheel A begins vertically below the carrier axis. Wheels B and C are spaced
/// by `2π/3`; the input carrier rotation is positive counter-clockwise.
pub fn step_tri_wheel_stair_contact(
    input: TriWheelStairInput,
) -> Result<TriWheelStairResult, TriWheelStairError> {
    validate_input(input)?;

    let mut wheel_centres_m = [[0.0; 2]; TRI_WHEEL_COUNT];
    let mut signed_vertical_gaps_m = [0.0; TRI_WHEEL_COUNT];
    let mut contact_mask = [false; TRI_WHEEL_COUNT];
    let mut contact_count = 0_u8;
    let mut minimum_gap_m = f64::INFINITY;
    let mut nearest_gap_m = f64::INFINITY;

    for wheel_index in 0..TRI_WHEEL_COUNT {
        let phase_rad = -PI / 2.0 + wheel_index as f64 * TAU / TRI_WHEEL_COUNT as f64;
        let angle_rad = phase_rad + input.carrier_rotation_rad;
        let x_m = input
            .cluster_radius_m
            .mul_add(angle_rad.cos(), input.axle_x_m);
        let y_m = input
            .cluster_radius_m
            .mul_add(angle_rad.sin(), input.axle_y_m);
        let support_y_m = horizontal_support_height_m(
            x_m,
            input.stair_rise_m,
            input.stair_tread_m,
            input.stair_active,
        );
        let gap_m = y_m - input.wheel_radius_m - support_y_m;
        if gap_m < -input.contact_tolerance_m {
            return Err(TriWheelStairError::PenetratingSupport { wheel_index, gap_m });
        }

        let touching = gap_m.abs() <= input.contact_tolerance_m;
        wheel_centres_m[wheel_index] = [x_m, y_m];
        signed_vertical_gaps_m[wheel_index] = gap_m;
        contact_mask[wheel_index] = touching;
        contact_count += u8::from(touching);
        minimum_gap_m = minimum_gap_m.min(gap_m);
        nearest_gap_m = nearest_gap_m.min(gap_m.abs());
    }

    if contact_count == 0 {
        return Err(TriWheelStairError::Unsupported { nearest_gap_m });
    }

    Ok(TriWheelStairResult {
        wheel_centres_m,
        signed_vertical_gaps_m,
        contact_mask,
        contact_count,
        minimum_gap_m,
    })
}

fn validate_input(input: TriWheelStairInput) -> Result<(), TriWheelStairError> {
    for (field, value) in [
        ("cluster_radius_m", input.cluster_radius_m),
        ("wheel_radius_m", input.wheel_radius_m),
        ("axle_x_m", input.axle_x_m),
        ("axle_y_m", input.axle_y_m),
        ("carrier_rotation_rad", input.carrier_rotation_rad),
        ("stair_rise_m", input.stair_rise_m),
        ("stair_tread_m", input.stair_tread_m),
        ("contact_tolerance_m", input.contact_tolerance_m),
    ] {
        if !value.is_finite() {
            return Err(TriWheelStairError::InvalidInput(field));
        }
    }
    if input.cluster_radius_m <= 0.0 {
        return Err(TriWheelStairError::InvalidInput("cluster_radius_m"));
    }
    if input.wheel_radius_m <= 0.0 {
        return Err(TriWheelStairError::InvalidInput("wheel_radius_m"));
    }
    if input.stair_rise_m <= 0.0 {
        return Err(TriWheelStairError::InvalidInput("stair_rise_m"));
    }
    if input.stair_tread_m <= 0.0 {
        return Err(TriWheelStairError::InvalidInput("stair_tread_m"));
    }
    if input.contact_tolerance_m <= 0.0 || input.contact_tolerance_m > input.wheel_radius_m {
        return Err(TriWheelStairError::InvalidInput("contact_tolerance_m"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const INCH_M: f64 = 0.0254;
    const CLUSTER_RADIUS_M: f64 = 5.581 * INCH_M;
    const WHEEL_RADIUS_M: f64 = 3.81 * INCH_M;
    const STAIR_RISE_M: f64 = 6.85 * INCH_M;
    const STAIR_TREAD_M: f64 = 10.9 * INCH_M;
    const RISER_OFFSET_M: f64 = 3.011 * INCH_M;
    const TOLERANCE_M: f64 = 1.0e-9;

    fn start_rotation_rad() -> f64 {
        PI / 3.0 - (STAIR_RISE_M / (3.0_f64.sqrt() * CLUSTER_RADIUS_M)).asin()
    }

    #[test]
    fn one_wheel_supports_the_vertical_balance_pose() {
        let result = step_tri_wheel_stair_contact(TriWheelStairInput {
            cluster_radius_m: CLUSTER_RADIUS_M,
            wheel_radius_m: WHEEL_RADIUS_M,
            axle_x_m: 0.0,
            axle_y_m: WHEEL_RADIUS_M + CLUSTER_RADIUS_M,
            carrier_rotation_rad: 0.0,
            stair_rise_m: STAIR_RISE_M,
            stair_tread_m: STAIR_TREAD_M,
            stair_active: false,
            contact_tolerance_m: TOLERANCE_M,
        })
        .unwrap();

        assert_eq!(result.contact_mask, [true, false, false]);
        assert_eq!(result.contact_count, 1);
        assert!(result.minimum_gap_m.abs() <= TOLERANCE_M);
    }

    #[test]
    fn source_dimensioned_start_pose_touches_ground_and_first_tread() {
        let rotation = -start_rotation_rad();
        let wheel_a_relative_x = CLUSTER_RADIUS_M * (-PI / 2.0 + rotation).cos();
        let wheel_a_relative_y = CLUSTER_RADIUS_M * (-PI / 2.0 + rotation).sin();
        let result = step_tri_wheel_stair_contact(TriWheelStairInput {
            cluster_radius_m: CLUSTER_RADIUS_M,
            wheel_radius_m: WHEEL_RADIUS_M,
            axle_x_m: -RISER_OFFSET_M - wheel_a_relative_x,
            axle_y_m: WHEEL_RADIUS_M - wheel_a_relative_y,
            carrier_rotation_rad: rotation,
            stair_rise_m: STAIR_RISE_M,
            stair_tread_m: STAIR_TREAD_M,
            stair_active: true,
            contact_tolerance_m: TOLERANCE_M,
        })
        .unwrap();

        assert_eq!(result.contact_mask, [true, true, false]);
        assert_eq!(result.contact_count, 2);
    }

    #[test]
    fn source_dimensioned_climb_pose_touches_successive_treads() {
        let rotation = -start_rotation_rad() - 2.0 * PI / 3.0;
        let wheel_b_angle = PI / 6.0 + rotation;
        let wheel_b_relative_x = CLUSTER_RADIUS_M * wheel_b_angle.cos();
        let wheel_b_relative_y = CLUSTER_RADIUS_M * wheel_b_angle.sin();
        let result = step_tri_wheel_stair_contact(TriWheelStairInput {
            cluster_radius_m: CLUSTER_RADIUS_M,
            wheel_radius_m: WHEEL_RADIUS_M,
            axle_x_m: STAIR_TREAD_M - RISER_OFFSET_M - wheel_b_relative_x,
            axle_y_m: STAIR_RISE_M + WHEEL_RADIUS_M - wheel_b_relative_y,
            carrier_rotation_rad: rotation,
            stair_rise_m: STAIR_RISE_M,
            stair_tread_m: STAIR_TREAD_M,
            stair_active: true,
            contact_tolerance_m: TOLERANCE_M,
        })
        .unwrap();

        assert_eq!(result.contact_mask, [false, true, true]);
        assert_eq!(result.contact_count, 2);
    }

    #[test]
    fn refuses_penetrating_and_floating_clusters() {
        let base = TriWheelStairInput {
            cluster_radius_m: CLUSTER_RADIUS_M,
            wheel_radius_m: WHEEL_RADIUS_M,
            axle_x_m: 0.0,
            axle_y_m: WHEEL_RADIUS_M + CLUSTER_RADIUS_M,
            carrier_rotation_rad: 0.0,
            stair_rise_m: STAIR_RISE_M,
            stair_tread_m: STAIR_TREAD_M,
            stair_active: false,
            contact_tolerance_m: TOLERANCE_M,
        };

        assert!(matches!(
            step_tri_wheel_stair_contact(TriWheelStairInput {
                axle_y_m: base.axle_y_m - 0.01,
                ..base
            }),
            Err(TriWheelStairError::PenetratingSupport { .. })
        ));
        assert!(matches!(
            step_tri_wheel_stair_contact(TriWheelStairInput {
                axle_y_m: base.axle_y_m + 0.01,
                ..base
            }),
            Err(TriWheelStairError::Unsupported { .. })
        ));
    }
}
