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
//! force. Position projection uses bounded Newton corrections on SE(3);
//! velocity projection is the exact mass-metric tangent projection at the
//! supplied pose.

use crate::{MassProperties, Vec3 as MbdVec3};
use fs_ga::{Se3, Twist, Vec3};
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
        let value = angular_row.dot(inverse_angular)
            + linear_row.dot(linear_row) / self.properties.mass();
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(invalid("rigid point/plane mobility"))
        }
    }

    /// Project a complete candidate and return measured residuals.
    pub fn project_candidate(
        &self,
        pose: &mut Se3,
        twist_body: &mut Twist,
    ) -> Result<PointPlaneProjectionReceipt, Se3Error> {
        let position_violation_removed_m = self.project_position(pose)?;
        let velocity_violation_removed_mps = self.project_velocity(pose, twist_body)?;
        Ok(PointPlaneProjectionReceipt {
            position_violation_removed_m,
            velocity_violation_removed_mps,
            final_gap_m: self.gap(*pose)?,
            final_normal_speed_mps: self.normal_speed(*pose, *twist_body)?,
        })
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

    fn correction(&self, pose: Se3, residual: f64) -> Result<Twist, Se3Error> {
        let (angular_row, linear_row) = self.body_row(pose)?;
        let inverse_angular = self.inverse_inertia(angular_row)?;
        let mobility = angular_row.dot(inverse_angular)
            + linear_row.dot(linear_row) / self.properties.mass();
        if !mobility.is_finite() || mobility <= 0.0 {
            return Err(invalid("rigid point/plane mobility"));
        }
        let multiplier = -residual / mobility;
        if !multiplier.is_finite() {
            return Err(invalid("rigid point/plane multiplier"));
        }
        Ok(Twist::new(
            inverse_angular.scale(multiplier),
            linear_row.scale(multiplier / self.properties.mass()),
        ))
    }
}

impl RattleProjection for PointPlaneRattle {
    fn project_position(&self, pose: &mut Se3) -> Result<f64, Se3Error> {
        let initial = self.gap(*pose)?;
        if initial.abs() <= self.controls.position_tolerance_m {
            return Ok(0.0);
        }
        let mut residual = initial;
        for iteration in 1..=self.controls.maximum_position_iterations {
            *pose = pose.body_plus(self.correction(*pose, residual)?)?;
            residual = self.gap(*pose)?;
            if residual.abs() <= self.controls.position_tolerance_m {
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
        *twist = twist.plus(self.correction(*pose, initial)?);
        let residual = self.normal_speed(*pose, *twist)?;
        if residual.abs() > self.controls.velocity_tolerance_mps {
            return Err(Se3Error::SolverDiverged {
                iters: 1,
                residual: residual.abs(),
            });
        }
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
