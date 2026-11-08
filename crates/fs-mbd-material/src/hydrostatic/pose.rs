//! Shared gravity-aligned coordinates for floating equilibria and stiffness.

use fs_mbd::{Gravity, Pose, RigidBodyState, Vec3};

use super::{HydrostaticError, gravity_direction};

/// World axes remain independent of arbitrary principal-axis ordering.
pub(super) struct HydrostaticPoseAxes {
    pub up: Vec3,
    pub roll: Vec3,
    pub pitch: Vec3,
}

impl HydrostaticPoseAxes {
    pub fn new(gravity: Gravity, roll_direction_world: Vec3) -> Result<Self, HydrostaticError> {
        let (down, _) = gravity_direction(gravity)?;
        let up = Vec3::new(-down.x, -down.y, -down.z);
        let scale = roll_direction_world
            .x
            .abs()
            .max(roll_direction_world.y.abs())
            .max(roll_direction_world.z.abs());
        if !roll_direction_world.is_finite() || scale == 0.0 {
            return Err(HydrostaticError::Invalid(
                "finite nonzero roll direction required",
            ));
        }
        // Component division also admits very small finite headings without
        // overflowing a reciprocal before the normalization.
        let heading = Vec3::new(
            roll_direction_world.x / scale,
            roll_direction_world.y / scale,
            roll_direction_world.z / scale,
        );
        // Remove the parallel component through cross products, then build the
        // roll axis back from a unit transverse vector. A single subtraction
        // loses orthogonality for an almost gravity-parallel input heading.
        let transverse = up.cross(heading);
        let horizontal = transverse.cross(up);
        let length = horizontal.norm_squared().sqrt();
        if !length.is_finite() || length <= 64.0 * f64::EPSILON {
            return Err(HydrostaticError::Invalid(
                "roll direction is parallel to gravity",
            ));
        }
        let roll = horizontal.scale(length.recip());
        let pitch = up.cross(roll);
        let pitch = pitch.scale(pitch.norm_squared().sqrt().recip());
        Ok(Self {
            up,
            roll: pitch.cross(up),
            pitch,
        })
    }

    /// Rotate about the current COM with a world exponential, then heave.
    /// No momentum is inferred from a static pose search or perturbation.
    pub fn displaced_state(
        &self,
        reference: RigidBodyState,
        offset_m: f64,
        tilt_radians: [f64; 2],
    ) -> Result<RigidBodyState, HydrostaticError> {
        if !offset_m.is_finite() || tilt_radians.iter().any(|x| !x.is_finite()) {
            return Err(HydrostaticError::Invalid(
                "finite hydrostatic pose offsets required",
            ));
        }
        let original = reference.pose();
        let rotation_world = self
            .roll
            .scale(tilt_radians[0])
            .add(self.pitch.scale(tilt_radians[1]));
        let orientation = if tilt_radians == [0.0, 0.0] {
            original.orientation()
        } else {
            // R Exp(R^-1 delta_world) = Exp(delta_world) R. Keep rotation
            // arithmetic with the native quaternion owner.
            original
                .orientation()
                .right_exp(original.orientation().rotate_world_to_body(rotation_world))?
        };
        Ok(RigidBodyState::new(
            Pose::new(
                original.position_world().add(self.up.scale(offset_m)),
                orientation,
            )?,
            Vec3::ZERO,
            Vec3::ZERO,
        )?)
    }
}
