//! Held world angular velocity and a translating pivot for the existing wake.

use fs_ga::{So3, So3Tangent, Vec3};

use super::{ImagePlane, Refusal, refusal};

/// Prescribed rigid lifting-line motion over one shedding interval.
///
/// All vectors use fixed world coordinates. The pivot starts at P, translates
/// at the held velocity V, and the line rotates about it at the held angular
/// velocity omega: x(t) = P + V*t + Exp(omega*t) * (x(0) - P).
/// Positive rotation follows the right-hand rule; elapsed t is in seconds.
/// The initial line is the owner's currently accepted geometry, not a second
/// pose state. This supports finite angles and multiple turns without chord
/// interpolation, shortening the span, or choosing a shortest rotation arc.
///
/// This is a prescribed kinematic law, not a force/torque-driven body solver.
/// The material point velocity is V + omega cross (x(t) - P - V*t).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RigidLineMotion {
    /// Pivot position at the start of this interval [m].
    pub pivot_start_m: [f64; 3],
    /// Constant world velocity of the pivot [m/s].
    pub pivot_velocity_mps: [f64; 3],
    /// Constant world angular velocity [rad/s], not a finite endpoint angle.
    pub angular_velocity_world_rad_s: [f64; 3],
}

impl RigidLineMotion {
    /// Sample a material point's world position [m] and velocity [m/s].
    ///
    /// The input point is its accepted position at this interval's start.
    /// Use the same elapsed seconds as the wake's incident-field callback.
    /// A later shedding interval starts from its newly accepted line and a
    /// correspondingly advanced pivot; no cumulative time is stored here.
    ///
    /// # Errors
    /// Non-finite input, negative/non-finite elapsed time, or an unrepresentable
    /// rotation, position or material velocity.
    pub fn point_and_velocity(
        self,
        start_m: [f64; 3],
        elapsed_s: f64,
    ) -> Result<([f64; 3], [f64; 3]), Refusal> {
        let stage = self.stage(elapsed_s)?;
        let (point, arm) = stage.position_and_arm(start_m)?;
        let velocity = if let Some(arm) = arm {
            let spin = vector(self.angular_velocity_world_rad_s).cross(arm);
            [
                self.pivot_velocity_mps[0] + spin.x,
                self.pivot_velocity_mps[1] + spin.y,
                self.pivot_velocity_mps[2] + spin.z,
            ]
        } else {
            self.pivot_velocity_mps
        };
        if velocity.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable rigid-line material velocity"));
        }
        Ok((point, velocity))
    }

    fn stage(self, elapsed_s: f64) -> Result<RigidStage, Refusal> {
        if !self
            .pivot_start_m
            .into_iter()
            .chain(self.pivot_velocity_mps)
            .chain(self.angular_velocity_world_rad_s)
            .all(f64::is_finite)
            || !elapsed_s.is_finite()
            || elapsed_s < 0.0
        {
            return Err(invalid("finite rigid motion and nonnegative elapsed seconds required"));
        }
        let translation = self.pivot_velocity_mps.map(|v| v * elapsed_s);
        if translation.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable rigid-line translation"));
        }
        let rotation = if self.angular_velocity_world_rad_s == [0.0; 3] {
            None
        } else {
            let angle = self.angular_velocity_world_rad_s.map(|v| v * elapsed_s);
            Some(
                So3::exp(So3Tangent::new(vector(angle)))
                    .map_err(|_| invalid("unrepresentable rigid-line rotation"))?,
            )
        };
        Ok(RigidStage {
            pivot: self.pivot_start_m,
            elapsed_s,
            translation,
            rotation,
        })
    }
}

/// Rigid lifting-line prescription with the existing stationary slip plane.
/// Free wake vertices remain material fluid points in world coordinates.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RigidWakeBoundary {
    /// Finite rotation and pivot translation during this interval.
    pub line_motion: RigidLineMotion,
    /// Optional stationary horizontal slip plane. As in WakeBoundary,
    /// clearance is checked at stages and the endpoint, not continuously.
    pub image_plane: Option<ImagePlane>,
}

pub(super) struct PreparedRigidMotion {
    motion: RigidLineMotion,
    dt_s: f64,
    pub(super) travel_m: f64,
}

impl PreparedRigidMotion {
    pub(super) fn new(
        motion: RigidLineMotion,
        dt_s: f64,
        start: &[[f64; 3]],
        plane: Option<ImagePlane>,
    ) -> Result<Self, Refusal> {
        let endpoint = motion.stage(dt_s)?;
        let angle = vector(motion.angular_velocity_world_rad_s.map(|v| v * dt_s));
        let mut rotation_travel = 0.0_f64;
        for &point in start {
            let (end, _) = endpoint.position_and_arm(point)?;
            if plane.is_some_and(|p| end[2] <= p.z0_m) {
                return Err(invalid("rigid-line endpoint crosses the image plane"));
            }
            if endpoint.rotation.is_some() {
                let arm = vector(std::array::from_fn(|i| point[i] - motion.pivot_start_m[i]));
                // |omega cross arm| is constant under this fixed-axis rotation.
                // Scale by dt before the cross product so the guard describes
                // travel, not a possibly overflowing instantaneous speed.
                let arc = angle.cross(arm);
                let distance = arc.x.hypot(arc.y).hypot(arc.z);
                if !distance.is_finite() {
                    return Err(invalid("unrepresentable rigid-line arc travel"));
                }
                rotation_travel = rotation_travel.max(distance);
            }
        }
        let shift = endpoint.translation;
        let travel_m = shift[0].hypot(shift[1]).hypot(shift[2]) + rotation_travel;
        if !travel_m.is_finite() {
            return Err(invalid("unrepresentable rigid-line path travel"));
        }
        Ok(Self { motion, dt_s, travel_m })
    }

    pub(super) fn stage(&self, fraction: f64) -> Result<RigidStage, Refusal> {
        self.motion.stage(fraction * self.dt_s)
    }
}

pub(super) struct RigidStage {
    pivot: [f64; 3],
    elapsed_s: f64,
    translation: [f64; 3],
    rotation: Option<So3>,
}

impl RigidStage {
    pub(super) fn position(&self, start: [f64; 3]) -> Result<[f64; 3], Refusal> {
        self.position_and_arm(start).map(|value| value.0)
    }

    fn position_and_arm(
        &self,
        start: [f64; 3],
    ) -> Result<([f64; 3], Option<Vec3>), Refusal> {
        if start.iter().any(|v| !v.is_finite()) {
            return Err(invalid("non-finite rigid-line initial point"));
        }
        let (displacement, arm) = if let Some(rotation) = self.rotation {
            let original = vector(std::array::from_fn(|i| start[i] - self.pivot[i]));
            let rotated = rotation
                .rotate(original)
                .map_err(|_| invalid("unrepresentable rigid-line pivot arm"))?;
            let delta = [rotated.x - original.x, rotated.y - original.y, rotated.z - original.z];
            (delta, Some(rotated))
        } else {
            ([0.0; 3], None)
        };
        // Apply the rotational increment to the accepted world point, avoiding
        // a subtract/add of a distant pivot in the zero-motion limit.
        let point = if self.elapsed_s == 0.0 {
            start
        } else {
            std::array::from_fn(|i| {
                let moved = if displacement[i] == 0.0 {
                    start[i]
                } else {
                    start[i] + displacement[i]
                };
                if self.translation[i] == 0.0 {
                    moved
                } else {
                    moved + self.translation[i]
                }
            })
        };
        if point.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable rigid-line world position"));
        }
        Ok((point, arm))
    }
}

fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn invalid(message: &str) -> Refusal {
    refusal("wake-advection-line-invalid", message)
}
