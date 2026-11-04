//! Prescribed relative rotation through the original finite-body graph owners.

use super::{
    Dd, DistanceGraphControls, DistanceGraphError, DistanceStepError, MassProperties,
    RigidBodyState, Row, RowKind, Vec3, driven, positive, prismatic,
};
use crate::UnitQuaternion;

const AXIS: Vec3 = Vec3::new(0.0, 0.0, 1.0);

/// Held relative angular acceleration of a finite revolute actuator.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AngularDrive {
    /// Initial spin of body two relative to body one about joint-local Z [rad/s].
    pub angular_rate_rad_s: f64,
    /// Constant signed relative spin acceleration for this interval [rad/s²].
    pub angular_acceleration_rad_s2: f64,
}

/// A finite revolute joint whose remaining axial angle is prescribed.
///
/// Both material points coincide. The local joint frames initially coincide
/// in the world; the commanded spin subsequently turns the first joint frame
/// about its own Z axis. Body two follows that moving frame, while both bodies
/// retain finite inertia and recoil. Frames map joint coordinates into native
/// COM body coordinates and need not follow either inertia's principal axes.
///
/// Six independent rows (three point closures and three orientation equations)
/// enter the existing simultaneous graph. This is a prescribed ideal actuator,
/// with signed source work, not a finite-power motor or a fixed support.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RevoluteMotorJoint {
    bodies: [usize; 2],
    points_body_m: [Vec3; 2],
    frames_body: [UnitQuaternion; 2],
    angular_length_m: f64,
}
impl RevoluteMotorJoint {
    /// Prepare the current commanded joint frames and their coincident points.
    ///
    /// `angular_length_m = L` gives the same full-frame sine residual, local
    /// chart and complete-couple bound as [`super::PrismaticJoint`]. Thus
    /// angular gap/rate tolerances are scaled by L; L must exceed the graph's
    /// position tolerance. The three point impulses have their separate
    /// complete-vector bound. The graph admits actual body indices and initial
    /// constraints, including the commanded relative angular rate.
    pub fn new(
        bodies: [usize; 2],
        points_body_m: [Vec3; 2],
        frames_body: [UnitQuaternion; 2],
        angular_length_m: f64,
    ) -> Result<Self, DistanceStepError> {
        if bodies[0] == bodies[1]
            || points_body_m.iter().any(|point| !point.is_finite())
            || !positive(angular_length_m)
        {
            return Err(DistanceStepError::Invalid("finite revolute motor geometry"));
        }
        Ok(Self {
            bodies,
            points_body_m,
            frames_body,
            angular_length_m,
        })
    }

    /// First reacting body, then the body driven in the positive joint-Z sense.
    #[must_use]
    pub const fn bodies(self) -> [usize; 2] {
        self.bodies
    }
    /// Coincident material points measured from each native COM [m].
    #[must_use]
    pub const fn points_body_m(self) -> [Vec3; 2] {
        self.points_body_m
    }
    /// Current commanded joint-to-body frames, retained in the continuation.
    #[must_use]
    pub const fn frames_body(self) -> [UnitQuaternion; 2] {
        self.frames_body
    }
    /// Dimensional scale of the angular constraints and couple bound [m].
    #[must_use]
    pub const fn angular_length_m(self) -> f64 {
        self.angular_length_m
    }

    pub(super) fn write_rows(self, rows: &mut [Row]) {
        for component in 0..3 {
            rows[component] = Row {
                bodies: self.bodies,
                points_body_m: self.points_body_m,
                kind: RowKind::Cartesian(component),
            };
            rows[3 + component] = Row {
                bodies: self.bodies,
                points_body_m: [Vec3::ZERO; 2],
                kind: RowKind::FrameLock(prismatic::AngularConstraint {
                    frames_body: self.frames_body,
                    length_m: self.angular_length_m,
                    component,
                }),
            };
        }
    }
}

/// Measured representation of one compensated quadratic angular command.
///
/// The original polynomial is evaluated before binary64 rounding. Quaternion
/// frame preparation uses the native exponential map. Axial phase and unwanted
/// transverse tilt are reconstructed from the actual prepared frame pair with
/// compensated products; no unwrapped absolute angle is inferred from a pose.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AngularDriveStep {
    /// Represented interval displacement `rate*h + acceleration*h²/2` [rad].
    pub displacement_rad: f64,
    /// Prepared frame's axial angle minus the compensated displacement [rad].
    pub displacement_representation_residual_rad: f64,
    /// Transverse tilt introduced while preparing the advanced frame [rad].
    pub transverse_frame_residual_rad: f64,
    /// Represented continuation rate minus its compensated command [rad/s].
    pub rate_representation_residual_rad_s: f64,
    /// Largest absolute command excursion from its initial angle, including
    /// an interior quadratic turning point [rad].
    pub maximum_excursion_rad: f64,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct AngularTarget {
    pub frame_after: UnitQuaternion,
    pub axis_first_body: Vec3,
    pub report: AngularDriveStep,
}

pub(super) fn prepare(
    index: usize,
    joint: RevoluteMotorJoint,
    drive: AngularDrive,
    duration: f64,
    controls: DistanceGraphControls,
) -> Result<driven::Target, DistanceGraphError> {
    let error = |reason| DistanceGraphError::PrescribedMotion {
        constraint: index,
        reason,
    };
    if !drive.angular_rate_rad_s.is_finite()
        || !drive.angular_acceleration_rad_s2.is_finite()
        || joint.angular_length_m <= controls.position_tolerance_m
    {
        return Err(error(
            "finite angular command and resolved angular scale required",
        ));
    }
    let (displacement, final_rate) = driven::motion_increment(
        drive.angular_rate_rad_s,
        drive.angular_acceleration_rad_s2,
        duration,
    );
    let represented = displacement.to_f64();
    let rate_residual = Dd::from_f64(final_rate.to_f64()) - final_rate;
    if !driven::finite(displacement)
        || !driven::finite(final_rate)
        || !represented.is_finite()
        || !final_rate.to_f64().is_finite()
        || !driven::finite(rate_residual)
    {
        return Err(error("angular endpoint prescription is unrepresentable"));
    }
    let mut excursion = displacement.abs();
    let rate = Dd::from_f64(drive.angular_rate_rad_s);
    let acceleration = Dd::from_f64(drive.angular_acceleration_rad_s2);
    if (rate.lt(Dd::ZERO) && Dd::ZERO.lt(final_rate))
        || (Dd::ZERO.lt(rate) && final_rate.lt(Dd::ZERO))
    {
        let turn = -rate / acceleration;
        let peak = (turn * (rate + Dd::from_f64(0.5) * acceleration * turn)).abs();
        if !driven::finite(peak) {
            return Err(error("angular interior command is unrepresentable"));
        }
        if excursion.lt(peak) {
            excursion = peak;
        }
    }
    // Two finite bodies can each travel at most the native angular cap. This
    // independent command check also rejects hidden winding/out-and-back paths
    // whose endpoint quaternion alone would look harmless.
    let maximum = Dd::from_f64(2.0 * controls.step.maximum_rotation_increment_rad);
    if !driven::finite(excursion) || maximum.lt(excursion) {
        return Err(error(
            "angular command exceeds the local interval excursion",
        ));
    }
    let frame_before = joint.frames_body[0];
    let frame_after = if displacement == Dd::ZERO {
        frame_before
    } else {
        frame_before.right_exp(AXIS.scale(represented))?
    };
    if displacement != Dd::ZERO && frame_after == frame_before {
        return Err(error("nonzero angular displacement vanished in its frame"));
    }
    let mut relative = relative_frame(frame_before, frame_after);
    // The stored native frames are canonical individually. Their relative
    // quaternion can change sign while crossing pi in the accumulated frame;
    // the admitted small command has the positive scalar representative.
    if relative[0].lt(Dd::ZERO) {
        relative = relative.map(|part| -part);
    }
    let phase = 2.0 * relative[3].to_f64().atan2(relative[0].to_f64());
    let tilt = 2.0
        * relative[1]
            .to_f64()
            .hypot(relative[2].to_f64())
            .atan2(relative[0].to_f64().hypot(relative[3].to_f64()));
    let phase_residual = (Dd::from_f64(phase) - displacement).to_f64();
    let report = AngularDriveStep {
        displacement_rad: represented,
        displacement_representation_residual_rad: phase_residual,
        transverse_frame_residual_rad: tilt,
        rate_representation_residual_rad_s: rate_residual.to_f64(),
        maximum_excursion_rad: excursion.to_f64(),
    };
    let preparation = representation_error_m(report, joint.angular_length_m);
    let rate_error = (rate_residual.abs() * Dd::from_f64(joint.angular_length_m)).to_f64();
    if !preparation.is_finite()
        || !rate_error.is_finite()
        || preparation > controls.position_tolerance_m
        || rate_error > controls.velocity_tolerance_mps
    {
        return Err(error(
            "angular continuation exceeds its representation tolerances",
        ));
    }
    Ok(driven::Target {
        initial_rate: drive.angular_rate_rad_s,
        final_rate,
        angular: Some(AngularTarget {
            frame_after,
            axis_first_body: frame_before.rotate_body_to_world(AXIS),
            report,
        }),
        ..driven::Target::default()
    })
}

pub(super) fn continuation(
    joint: RevoluteMotorJoint,
    drive: AngularDrive,
    target: driven::Target,
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    controls: DistanceGraphControls,
) -> Result<(RevoluteMotorJoint, AngularDrive, AngularDriveStep), DistanceGraphError> {
    let angular = target
        .angular
        .ok_or(DistanceStepError::Invalid("motor target missing"))?;
    let joint_after = RevoluteMotorJoint {
        frames_body: [angular.frame_after, joint.frames_body[1]],
        ..joint
    };
    let drive_after = AngularDrive {
        angular_rate_rad_s: target.final_rate.to_f64(),
        ..drive
    };
    let spin = Some((
        joint_after.frames_body[0].rotate_body_to_world(AXIS),
        Dd::from_f64(drive_after.angular_rate_rad_s),
    ));
    let (mut gap, mut rate) = (0.0_f64, 0.0_f64);
    for component in 0..3 {
        let g = prismatic::geometry_with_spin(
            states,
            properties,
            prismatic::AngularConstraint {
                frames_body: joint_after.frames_body,
                length_m: joint.angular_length_m,
                component,
            },
            spin,
        )?;
        gap = gap.hypot(g.gap);
        rate = rate.hypot(g.rate);
    }
    // Admit the actual continuation frame/rate as well as the original command.
    // The separate measured preparation error is never projected out of a body.
    if gap + representation_error_m(angular.report, joint.angular_length_m)
        > controls.position_tolerance_m
        || rate > controls.velocity_tolerance_mps
    {
        return Err(DistanceStepError::NotConverged.into());
    }
    Ok((joint_after, drive_after, angular.report))
}

pub(super) fn representation_error_m(report: AngularDriveStep, length: f64) -> f64 {
    ((Dd::from_f64(report.displacement_representation_residual_rad).abs()
        + Dd::from_f64(report.transverse_frame_residual_rad))
        * Dd::from_f64(length))
    .to_f64()
}

// Exact represented-input products retain the tiny preparation rotation beside
// the unit diagonal; this does not reuse the rounded production frame product.
fn relative_frame(a: UnitQuaternion, b: UnitQuaternion) -> [Dd; 4] {
    let [w, x, y, z] = a.components().map(Dd::from_f64);
    let [s, u, v, t] = b.components().map(Dd::from_f64);
    [
        w * s + x * u + y * v + z * t,
        w * u - x * s - y * t + z * v,
        w * v + x * t - y * s - z * u,
        w * t - x * v + y * u - z * s,
    ]
}
