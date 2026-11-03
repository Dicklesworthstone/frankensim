//! Full-inertia RATTLE projection for one rigid material point on a fixed plane.
//!
//! The projection uses the body's actual kinetic metric rather than translating
//! the reference point by hand. For the scalar constraint
//! `g = n · (x + R r - p) = 0`, the body-frame Jacobian is
//! `J = [r × n_b, n_b]`. Corrections are therefore
//! `delta = M^-1 J^T lambda`, retaining products of inertia and the exact
//! translational/rotational coupling of an offset point.
//!
//! This module implements the `fs_time::se3::RattleProjection` seam. It is a
//! projection owner, not collision detection, contact response, or a penalty
//! force. Position projection uses bounded Newton corrections of COM position
//! and body attitude; the COM correction follows the fixed plane normal.
//! The composed candidate projection carries world linear momentum and central
//! spin across the configuration repair before applying the normal mass-metric
//! impulse. This preserves the plane's translation/rotation symmetries without
//! treating a changed body basis as a physical impulse. It is not a complete
//! variational RATTLE time integrator.

/// Constrained point/plane evolution using native physical impulse responses.
#[path = "rigid_rattle/step.rs"]
pub mod step;

/// Rigid distance between material points on two finite native rigid bodies.
#[path = "rigid_rattle/distance.rs"]
pub mod distance;

#[path = "rigid_rattle/pair_step.rs"]
mod pair_step;

/// A material point sliding on a plane carried by another finite rigid body.
#[path = "rigid_rattle/finite_plane.rs"]
pub mod finite_plane;

use crate::{MassProperties, Vec3 as MbdVec3};
use fs_ga::{Se3, So3Tangent, Twist, Vec3};
use fs_time::se3::{RattleProjection, Se3Error};

/// Bounded controls for one point/plane projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointPlaneRattleControls {
    /// Maximum absolute signed gap after configuration projection [m].
    pub position_tolerance_m: f64,
    /// Maximum absolute normal material-point speed after tangent projection [m/s].
    pub velocity_tolerance_mps: f64,
    /// Maximum nonlinear configuration corrections.
    pub maximum_position_iterations: u32,
}

impl Default for PointPlaneRattleControls {
    fn default() -> Self {
        Self {
            position_tolerance_m: 1e-12,
            velocity_tolerance_mps: 1e-12,
            maximum_position_iterations: 12,
        }
    }
}

/// A rigid material point constrained to one fixed world plane.
///
/// The body reference point must be its center of mass, as required by
/// [`MassProperties`]. The plane normal is normalized once at construction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointPlaneRattle {
    properties: MassProperties,
    point_body_m: Vec3,
    plane_point_world_m: Vec3,
    normal_world: Vec3,
    controls: PointPlaneRattleControls,
}

/// Actual projection diagnostics for a candidate pose/twist pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointPlaneProjectionReceipt {
    /// Magnitude removed by the position projection [m].
    pub position_violation_removed_m: f64,
    /// Magnitude removed by the tangent projection [m/s].
    pub velocity_violation_removed_mps: f64,
    /// Final signed plane gap [m].
    pub final_gap_m: f64,
    /// Final signed normal material-point speed [m/s].
    pub final_normal_speed_mps: f64,
}

impl PointPlaneRattle {
    /// Prepare a full-inertia point/plane projector.
    ///
    /// # Errors
    /// Refuses non-finite geometry, a zero normal, invalid tolerances, or a
    /// zero iteration budget.
    pub fn new(
        properties: MassProperties,
        point_body_m: Vec3,
        plane_point_world_m: Vec3,
        normal_world: Vec3,
        controls: PointPlaneRattleControls,
    ) -> Result<Self, Se3Error> {
        if !finite_vec(point_body_m)
            || !finite_vec(plane_point_world_m)
            || !finite_vec(normal_world)
            || !controls.position_tolerance_m.is_finite()
            || controls.position_tolerance_m < 0.0
            || !controls.velocity_tolerance_mps.is_finite()
            || controls.velocity_tolerance_mps < 0.0
            || controls.maximum_position_iterations == 0
        {
            return Err(invalid("rigid point/plane RATTLE controls"));
        }
        let normal_world = normalize(normal_world)?;
        Ok(Self {
            properties,
            point_body_m,
            plane_point_world_m,
            normal_world,
            controls,
        })
    }

    /// Signed world-plane gap of the constrained material point.
    pub fn gap(&self, pose: Se3) -> Result<f64, Se3Error> {
        let point = pose.transform_point(self.point_body_m)?;
        let gap = self.normal_world.dot(point - self.plane_point_world_m);
        if gap.is_finite() {
            Ok(gap)
        } else {
            Err(invalid("rigid point/plane gap"))
        }
    }

    /// Signed normal speed of the constrained material point.
    pub fn normal_speed(&self, pose: Se3, twist_body: Twist) -> Result<f64, Se3Error> {
        let (angular_row, linear_row) = self.body_row(pose)?;
        let speed = angular_row.dot(twist_body.angular) + linear_row.dot(twist_body.linear);
        if speed.is_finite() {
            Ok(speed)
        } else {
            Err(invalid("rigid point/plane normal speed"))
        }
    }

    /// Scalar free mobility `J M^-1 J^T` at the supplied pose.
    pub fn normal_mobility(&self, pose: Se3) -> Result<f64, Se3Error> {
        let (angular_row, linear_row) = self.body_row(pose)?;
        let inverse_angular = self.inverse_inertia(angular_row)?;
        let value =
            angular_row.dot(inverse_angular) + linear_row.dot(linear_row) / self.properties.mass();
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(invalid("rigid point/plane mobility"))
        }
    }

    /// Project a complete candidate and return measured residuals.
    ///
    /// The configuration repair keeps world linear momentum and world central
    /// spin fixed, then the tangent projection applies one normal point impulse.
    /// COM translation during repair is normal to the plane. Consequently the
    /// plane-tangent linear and plane-normal total angular momenta are retained
    /// to roundoff. Other orbital components change with the explicit COM
    /// relocation; that change must not be converted into an artificial spin.
    /// The configuration repair can change kinetic energy for anisotropic
    /// inertia and does not claim a symplectic or energy-preserving time step.
    /// Refusal leaves both caller-owned inputs unchanged.
    pub fn project_candidate(
        &self,
        pose: &mut Se3,
        twist_body: &mut Twist,
    ) -> Result<PointPlaneProjectionReceipt, Se3Error> {
        let mut candidate_pose = *pose;
        let mut candidate_twist = *twist_body;
        let position_violation_removed_m = self.project_position(&mut candidate_pose)?;
        if candidate_pose.rotation() != pose.rotation() {
            let velocity_world = pose.rotation().rotate(twist_body.linear)?;
            let spin_body = self
                .properties
                .angular_momentum_body_checked(MbdVec3::new(
                    twist_body.angular.x,
                    twist_body.angular.y,
                    twist_body.angular.z,
                ))
                .map_err(|_| invalid("rigid point/plane angular momentum"))?;
            let spin_world =
                pose.rotation()
                    .rotate(Vec3::new(spin_body.x, spin_body.y, spin_body.z))?;
            let world_to_body = candidate_pose.rotation().inverse();
            candidate_twist = Twist::new(
                self.inverse_inertia(world_to_body.rotate(spin_world)?)?,
                world_to_body.rotate(velocity_world)?,
            );
        }
        let velocity_violation_removed_mps =
            self.project_velocity(&candidate_pose, &mut candidate_twist)?;
        let receipt = PointPlaneProjectionReceipt {
            position_violation_removed_m,
            velocity_violation_removed_mps,
            final_gap_m: self.gap(candidate_pose)?,
            final_normal_speed_mps: self.normal_speed(candidate_pose, candidate_twist)?,
        };
        *pose = candidate_pose;
        *twist_body = candidate_twist;
        Ok(receipt)
    }

    fn body_row(&self, pose: Se3) -> Result<(Vec3, Vec3), Se3Error> {
        let normal_body = pose.rotation().inverse().rotate(self.normal_world)?;
        Ok((self.point_body_m.cross(normal_body), normal_body))
    }

    fn inverse_inertia(&self, body: Vec3) -> Result<Vec3, Se3Error> {
        let result = self
            .properties
            .angular_velocity_body_checked(MbdVec3::new(body.x, body.y, body.z))
            .map_err(|_| invalid("rigid point/plane inertia response"))?;
        let result = Vec3::new(result.x, result.y, result.z);
        if finite_vec(result) {
            Ok(result)
        } else {
            Err(invalid("rigid point/plane inertia response"))
        }
    }

    // The scalar accompanies the body-coordinate correction so a position
    // repair can move the COM strictly along the world normal. A coupled SE(3)
    // screw exponential would add a tangential translation at finite rotation.
    fn correction(&self, pose: Se3, residual: f64) -> Result<(Twist, f64), Se3Error> {
        let (angular_row, linear_row) = self.body_row(pose)?;
        let inverse_angular = self.inverse_inertia(angular_row)?;
        let mobility =
            angular_row.dot(inverse_angular) + linear_row.dot(linear_row) / self.properties.mass();
        if !mobility.is_finite() || mobility <= 0.0 {
            return Err(invalid("rigid point/plane mobility"));
        }
        let multiplier = -residual / mobility;
        if !multiplier.is_finite() {
            return Err(invalid("rigid point/plane multiplier"));
        }
        let linear_amount = multiplier / self.properties.mass();
        Ok((
            Twist::new(
                inverse_angular.scale(multiplier),
                linear_row.scale(linear_amount),
            ),
            linear_amount,
        ))
    }
}

impl RattleProjection for PointPlaneRattle {
    fn project_position(&self, pose: &mut Se3) -> Result<f64, Se3Error> {
        let initial = self.gap(*pose)?;
        if initial.abs() <= self.controls.position_tolerance_m {
            return Ok(0.0);
        }
        let mut candidate = *pose;
        let mut residual = initial;
        for iteration in 1..=self.controls.maximum_position_iterations {
            let (correction, linear_amount) = self.correction(candidate, residual)?;
            candidate = Se3::from_parts(
                candidate
                    .rotation()
                    .body_plus(So3Tangent::new(correction.angular))?,
                candidate.translation() + self.normal_world.scale(linear_amount),
            )?;
            residual = self.gap(candidate)?;
            if residual.abs() <= self.controls.position_tolerance_m {
                *pose = candidate;
                return Ok(initial.abs());
            }
            if iteration == self.controls.maximum_position_iterations {
                return Err(Se3Error::SolverDiverged {
                    iters: iteration,
                    residual: residual.abs(),
                });
            }
        }
        unreachable!("positive iteration cap returns inside the loop")
    }

    fn project_velocity(&self, pose: &Se3, twist: &mut Twist) -> Result<f64, Se3Error> {
        let initial = self.normal_speed(*pose, *twist)?;
        if initial.abs() <= self.controls.velocity_tolerance_mps {
            return Ok(0.0);
        }
        let candidate = twist.plus(self.correction(*pose, initial)?.0);
        let residual = self.normal_speed(*pose, candidate)?;
        if residual.abs() > self.controls.velocity_tolerance_mps {
            return Err(Se3Error::SolverDiverged {
                iters: 1,
                residual: residual.abs(),
            });
        }
        *twist = candidate;
        Ok(initial.abs())
    }
}

fn finite_vec(value: Vec3) -> bool {
    value.x.is_finite() && value.y.is_finite() && value.z.is_finite()
}

fn normalize(value: Vec3) -> Result<Vec3, Se3Error> {
    let scale = value.x.abs().max(value.y.abs()).max(value.z.abs());
    if !scale.is_finite() || scale == 0.0 {
        return Err(invalid("rigid point/plane normal"));
    }
    let scaled = value.scale(1.0 / scale);
    let norm = scaled.dot(scaled).sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return Err(invalid("rigid point/plane normal"));
    }
    Ok(scaled.scale(1.0 / norm))
}

fn invalid(context: &'static str) -> Se3Error {
    Se3Error::InvalidParameter { context }
}
