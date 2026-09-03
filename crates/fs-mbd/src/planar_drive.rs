//! Deterministic planar differential-drive kinematics.
//!
//! This module owns the constant-twist SE(2) update for a two-wheel drive.
//! It is deliberately narrower than a tire/contact model: wheel speeds are
//! prescribed, the ground is an ideal kinematic reference, and no force,
//! traction, slip, friction, or motor performance is inferred.

use core::fmt;

/// Admitted pose and wheel coordinates for a planar differential drive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlanarDriveState {
    /// World x coordinate of the axle midpoint, in metres.
    pub x_m: f64,
    /// World y coordinate of the axle midpoint, in metres.
    pub y_m: f64,
    /// Chassis heading from world +x toward world +y, in radians.
    pub heading_rad: f64,
    /// Left wheel rotation coordinate, in radians.
    pub left_wheel_angle_rad: f64,
    /// Right wheel rotation coordinate, in radians.
    pub right_wheel_angle_rad: f64,
}

/// Prescribed wheel speeds and geometry for one fixed step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DifferentialDriveStep {
    /// Left wheel tangential speed at the ideal ground reference, in m/s.
    pub left_speed_mps: f64,
    /// Right wheel tangential speed at the ideal ground reference, in m/s.
    pub right_speed_mps: f64,
    /// Distance between the wheel center planes, in metres.
    pub track_width_m: f64,
    /// Effective rolling radius used only for wheel-angle display, in metres.
    pub wheel_radius_m: f64,
    /// Fixed step duration, in seconds.
    pub dt_s: f64,
}

/// Typed refusal from the planar differential-drive boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanarDriveError {
    /// A named input was NaN or infinite.
    NonFinite(&'static str),
    /// Track width must be finite and strictly positive.
    InvalidTrackWidth,
    /// Wheel radius must be finite and strictly positive.
    InvalidWheelRadius,
    /// The fixed step must be finite, positive, and no greater than 0.25 s.
    InvalidStepDuration,
    /// Finite inputs produced an unrepresentable output.
    UnrepresentableOutput,
}

impl fmt::Display for PlanarDriveError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(field) => write!(formatter, "{field} must be finite"),
            Self::InvalidTrackWidth => {
                formatter.write_str("track width must be finite and positive")
            }
            Self::InvalidWheelRadius => {
                formatter.write_str("wheel radius must be finite and positive")
            }
            Self::InvalidStepDuration => formatter
                .write_str("step duration must be finite, positive, and no greater than 0.25 s"),
            Self::UnrepresentableOutput => {
                formatter.write_str("planar-drive output is not representable as finite f64")
            }
        }
    }
}

impl std::error::Error for PlanarDriveError {}

fn require_finite(value: f64, field: &'static str) -> Result<(), PlanarDriveError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(PlanarDriveError::NonFinite(field))
    }
}

/// Advance a prescribed differential drive through one constant-twist step.
///
/// The chassis update is the exact SE(2) exponential for constant left and
/// right wheel speeds over `dt_s`; it does not use heading-first Euler drift.
///
/// # Errors
/// Refuses non-finite state/speed input, non-positive geometry, a step outside
/// `(0, 0.25]` seconds, or a non-finite derived coordinate.
pub fn step_differential_drive(
    state: PlanarDriveState,
    step: DifferentialDriveStep,
) -> Result<PlanarDriveState, PlanarDriveError> {
    for (value, field) in [
        (state.x_m, "state.x_m"),
        (state.y_m, "state.y_m"),
        (state.heading_rad, "state.heading_rad"),
        (state.left_wheel_angle_rad, "state.left_wheel_angle_rad"),
        (state.right_wheel_angle_rad, "state.right_wheel_angle_rad"),
        (step.left_speed_mps, "step.left_speed_mps"),
        (step.right_speed_mps, "step.right_speed_mps"),
    ] {
        require_finite(value, field)?;
    }
    if !step.track_width_m.is_finite() || step.track_width_m <= 0.0 {
        return Err(PlanarDriveError::InvalidTrackWidth);
    }
    if !step.wheel_radius_m.is_finite() || step.wheel_radius_m <= 0.0 {
        return Err(PlanarDriveError::InvalidWheelRadius);
    }
    if !step.dt_s.is_finite() || step.dt_s <= 0.0 || step.dt_s > 0.25 {
        return Err(PlanarDriveError::InvalidStepDuration);
    }

    let linear_speed_mps = 0.5 * (step.left_speed_mps + step.right_speed_mps);
    let angular_speed_rad_s = (step.right_speed_mps - step.left_speed_mps) / step.track_width_m;
    let next_heading = angular_speed_rad_s.mul_add(step.dt_s, state.heading_rad);

    let (next_x, next_y) = if angular_speed_rad_s.abs() <= 1.0e-12 {
        let distance = linear_speed_mps * step.dt_s;
        (
            distance.mul_add(state.heading_rad.cos(), state.x_m),
            distance.mul_add(state.heading_rad.sin(), state.y_m),
        )
    } else {
        let radius = linear_speed_mps / angular_speed_rad_s;
        (
            radius.mul_add(next_heading.sin() - state.heading_rad.sin(), state.x_m),
            (-radius).mul_add(next_heading.cos() - state.heading_rad.cos(), state.y_m),
        )
    };

    let next = PlanarDriveState {
        x_m: next_x,
        y_m: next_y,
        heading_rad: next_heading,
        left_wheel_angle_rad: (step.left_speed_mps / step.wheel_radius_m)
            .mul_add(step.dt_s, state.left_wheel_angle_rad),
        right_wheel_angle_rad: (step.right_speed_mps / step.wheel_radius_m)
            .mul_add(step.dt_s, state.right_wheel_angle_rad),
    };
    if [
        next.x_m,
        next.y_m,
        next.heading_rad,
        next.left_wheel_angle_rad,
        next.right_wheel_angle_rad,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        Ok(next)
    } else {
        Err(PlanarDriveError::UnrepresentableOutput)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state() -> PlanarDriveState {
        PlanarDriveState {
            x_m: 0.0,
            y_m: 0.0,
            heading_rad: 0.0,
            left_wheel_angle_rad: 0.0,
            right_wheel_angle_rad: 0.0,
        }
    }

    #[test]
    fn equal_wheel_speeds_translate_without_yaw() {
        let next = step_differential_drive(
            state(),
            DifferentialDriveStep {
                left_speed_mps: 0.4,
                right_speed_mps: 0.4,
                track_width_m: 0.24,
                wheel_radius_m: 0.04,
                dt_s: 0.1,
            },
        )
        .expect("valid straight step");
        assert!((next.x_m - 0.04).abs() < 1.0e-14);
        assert_eq!(next.y_m, 0.0);
        assert_eq!(next.heading_rad, 0.0);
        assert!((next.left_wheel_angle_rad - 1.0).abs() < 1.0e-14);
        assert_eq!(next.left_wheel_angle_rad, next.right_wheel_angle_rad);
    }

    #[test]
    fn opposite_wheel_speeds_spin_about_the_axle_midpoint() {
        let next = step_differential_drive(
            state(),
            DifferentialDriveStep {
                left_speed_mps: -0.12,
                right_speed_mps: 0.12,
                track_width_m: 0.24,
                wheel_radius_m: 0.04,
                dt_s: 0.1,
            },
        )
        .expect("valid spin step");
        assert_eq!(next.x_m, 0.0);
        assert_eq!(next.y_m, 0.0);
        assert!((next.heading_rad - 0.1).abs() < 1.0e-14);
        assert_eq!(next.left_wheel_angle_rad, -next.right_wheel_angle_rad);
    }

    #[test]
    fn constant_twist_arc_matches_the_closed_form_oracle() {
        let next = step_differential_drive(
            state(),
            DifferentialDriveStep {
                left_speed_mps: 0.1,
                right_speed_mps: 0.3,
                track_width_m: 0.2,
                wheel_radius_m: 0.05,
                dt_s: 0.2,
            },
        )
        .expect("valid arc step");
        assert!((next.heading_rad - 0.2).abs() < 1.0e-14);
        assert!((next.x_m - 0.2 * 0.2_f64.sin()).abs() < 1.0e-14);
        assert!((next.y_m - 0.2 * (1.0 - 0.2_f64.cos())).abs() < 1.0e-14);
    }

    #[test]
    fn invalid_geometry_time_and_non_finite_state_refuse() {
        let base = DifferentialDriveStep {
            left_speed_mps: 0.1,
            right_speed_mps: 0.1,
            track_width_m: 0.24,
            wheel_radius_m: 0.04,
            dt_s: 0.1,
        };
        assert_eq!(
            step_differential_drive(
                state(),
                DifferentialDriveStep {
                    track_width_m: 0.0,
                    ..base
                }
            ),
            Err(PlanarDriveError::InvalidTrackWidth)
        );
        assert_eq!(
            step_differential_drive(state(), DifferentialDriveStep { dt_s: 0.5, ..base }),
            Err(PlanarDriveError::InvalidStepDuration)
        );
        assert!(matches!(
            step_differential_drive(
                PlanarDriveState {
                    x_m: f64::NAN,
                    ..state()
                },
                base,
            ),
            Err(PlanarDriveError::NonFinite("state.x_m"))
        ));
    }
}
