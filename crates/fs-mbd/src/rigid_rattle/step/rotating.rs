//! Prescribed rigid support motion for the native bilateral point/plane solve.

use super::*;
use crate::{PointKinematics, UnitQuaternion};

/// One held world-frame rigid velocity field for an infinite plane.
/// The pivot can lie anywhere, including off the initial plane. It translates
/// at the held velocity while the plane rotates about it at the held angular
/// velocity. Tangential surface speed does no work on a frictionless reaction.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PlaneRigidMotion {
    /// Initial world position of the translating rotation pivot [m].
    pub pivot_world_m: Vec3,
    /// Constant world velocity of that pivot [m/s].
    pub pivot_velocity_world_mps: Vec3,
    /// Constant world angular velocity [rad/s].
    pub angular_velocity_world_rad_s: Vec3,
}

/// Native constrained motion plus the prescribed support's impulse/work account.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotatingPlaneStepReceipt {
    /// Original impulse/free-step/impulse map. The two reactions use the
    /// initial and final plane normals, respectively, at actual body points.
    pub step: PointPlaneStepReceipt,
    /// Represented endpoint plane, ready for continuation.
    pub plane_after: PointPlaneRattle,
    /// Same held velocity field, with the pivot advanced to the endpoint.
    pub motion_after: PlaneRigidMotion,
    /// Equal-and-opposite sum of both body reaction impulses [N s].
    pub support_impulse_world_ns: Vec3,
    /// Reciprocal angular impulse about the FIXED initial pivot [N m s].
    /// Both endpoint lever arms use that same origin, including pivot travel.
    pub support_angular_impulse_about_initial_pivot_world_nms: Vec3,
    /// Sum of each body impulse dotted with its physical support velocity [J].
    /// This is signed external work, not dissipated heat.
    pub boundary_work_j: f64,
    /// Impulse work minus support work, paired directly in the relative frame.
    /// This measured time-discretization defect is not projected away.
    pub constraint_work_defect_j: f64,
    /// Endpoint gap evaluated before rounding the prescribed plane geometry [m].
    /// Both this and the represented-plane residual must pass the original cap.
    pub prescribed_final_gap_m: f64,
    /// Endpoint relative speed using the unrounded prescribed normal/pivot [m/s].
    pub prescribed_final_normal_speed_mps: f64,
}

impl PointPlaneRattle {
    /// Advance a material point on a prescribed translating and rotating plane.
    ///
    /// With held `U`, `omega` and initial pivot `c`, the prescribed geometry is
    /// `c(t)=c+U*t`, `n(t)=R(omega*t)*n`, and
    /// `p(t)=c(t)+R(omega*t)*(p-c)`. The velocity constraint is
    /// `n(t)·(v_point-U-omega×(x_point-c(t)))=0`. Merely replacing the plane's
    /// endpoint normal would omit the rotating support's normal velocity.
    ///
    /// Each physical impulse uses its own endpoint normal; the bounded position
    /// solve still calls the original free rigid-body integrator. Both body and
    /// support angular travel must fit the existing local rotation cap. The
    /// source has prescribed motion, with no inferred finite inertia. The
    /// smooth bilateral model permits tensile reactions; no collision, friction
    /// or contact-mode selection is added. No allocation or pose repair occurs.
    ///
    /// The original scalar preconditioner remains local: failure of decreasing
    /// residuals requires a smaller step rather than a search for a remote root.
    /// Zero support velocity recovers the original stationary native step.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_rotating_plane(
        &self,
        state: RigidBodyState,
        integrator: RigidBodyIntegrator,
        wrench: Wrench,
        duration_s: f64,
        motion: PlaneRigidMotion,
        controls: PointPlaneStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotatingPlaneStepReceipt, PointPlaneStepError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !motion.pivot_world_m.is_finite()
            || !motion.pivot_velocity_world_mps.is_finite()
            || !motion.angular_velocity_world_rad_s.is_finite()
        {
            return Err(PointPlaneStepError::Invalid(
                "finite rigid plane motion required",
            ));
        }
        let omega = motion.angular_velocity_world_rad_s;
        let rate = omega.x.hypot(omega.y).hypot(omega.z);
        let angle = rate * duration_s;
        if !angle.is_finite()
            || (rate != 0.0 && angle == 0.0)
            || angle > controls.maximum_rotation_increment_rad
        {
            return Err(PointPlaneStepError::LocalBranch);
        }
        let rotation = if rate == 0.0 {
            UnitQuaternion::IDENTITY
        } else {
            UnitQuaternion::from_axis_angle(omega, angle)?
        };
        let native = |v: fs_ga::Vec3| Vec3::new(v.x, v.y, v.z);
        let pivot = dd_vec(motion.pivot_world_m);
        let velocity = dd_vec(motion.pivot_velocity_world_mps);
        let pivot_after = add(pivot, scale(velocity, Dd::from_f64(duration_s)));
        // Retain small rotations and translations until physical gap/speed
        // admission. Computing p+delta in f64 first could silently swallow a
        // commanded movement at a distant world origin.
        let prescribed_normal = rotate(rotation, dd_vec(native(self.normal_world)));
        let prescribed_point = add(
            add(
                dd_vec(native(self.plane_point_world_m)),
                scale(velocity, Dd::from_f64(duration_s)),
            ),
            rotation_increment(
                rotation,
                sub(dd_vec(native(self.plane_point_world_m)), pivot),
            ),
        );
        let end_pivot = value(pivot_after);
        let end_point = value(prescribed_point);
        let end_normal = value(prescribed_normal);
        if !end_pivot.is_finite() || !end_point.is_finite() || !end_normal.is_finite() {
            return Err(PointPlaneStepError::Invalid(
                "rigid plane endpoint is unrepresentable",
            ));
        }
        // A returned pivot must describe the same source on the next interval.
        // This bound is geometric, independent of an arbitrary energy datum.
        let pivot_error = value(sub(dd_vec(end_pivot), pivot_after));
        if pivot_error.x.hypot(pivot_error.y).hypot(pivot_error.z)
            > self.controls.position_tolerance_m
        {
            return Err(PointPlaneStepError::Invalid(
                "rigid plane pivot travel is unrepresentable",
            ));
        }
        let mut plane_after = *self;
        plane_after.plane_point_world_m = fs_ga::Vec3::new(end_point.x, end_point.y, end_point.z);
        plane_after.normal_world = fs_ga::Vec3::new(end_normal.x, end_normal.y, end_normal.z);
        let speed = |pivot| PlaneSpeed::Rigid {
            pivot,
            velocity: motion.pivot_velocity_world_mps,
            angular_velocity: omega,
        };
        let initial_speed = speed(pivot);
        let final_speed = speed(pivot_after);
        let stationary = motion.pivot_velocity_world_mps == Vec3::ZERO && omega == Vec3::ZERO;
        let step = if stationary {
            // In particular, an arbitrary stationary pivot cannot change a
            // pre-existing plane through an unnecessary subtraction roundtrip.
            plane_after = *self;
            self.advance(
                state,
                integrator,
                wrench,
                duration_s,
                controls,
                &mut cancelled,
            )?
        } else {
            self.advance_between(
                state,
                integrator,
                wrench,
                duration_s,
                controls,
                plane_after.plane_point_world_m,
                plane_after.normal_world,
                [initial_speed, final_speed],
                &mut cancelled,
            )?
        };
        let endpoint = step.velocity_impulse.point_kinematics_after;
        let prescribed_final_gap_m = gap(
            prescribed_normal,
            step.state_after,
            endpoint.arm_world,
            prescribed_point,
        );
        let prescribed_final_normal_speed_mps = relative_speed(
            prescribed_normal,
            point_velocity(endpoint),
            final_speed.velocity_at(step.state_after, endpoint),
        );
        if !prescribed_final_gap_m.is_finite()
            || !prescribed_final_normal_speed_mps.is_finite()
            || prescribed_final_gap_m.abs() > self.controls.position_tolerance_m
            || prescribed_final_normal_speed_mps.abs() > self.controls.velocity_tolerance_mps
        {
            return Err(PointPlaneStepError::NotConverged);
        }
        let mut impulse = [Dd::ZERO; 3];
        let mut angular_impulse = [Dd::ZERO; 3];
        let mut work = Dd::ZERO;
        let mut defect = Dd::ZERO;
        for (event, support) in [
            (step.position_impulse, initial_speed),
            (step.velocity_impulse, final_speed),
        ] {
            let j = dd_vec(event.impulse_world);
            let support_velocity =
                support.velocity_at(event.state_before, event.point_kinematics_before);
            let midpoint = scale(
                add(
                    point_velocity(event.point_kinematics_before),
                    point_velocity(event.point_kinematics_after),
                ),
                Dd::from_f64(0.5),
            );
            impulse = sub(impulse, j);
            angular_impulse = sub(
                angular_impulse,
                cross(
                    point_offset(
                        event.state_before,
                        event.point_kinematics_before.arm_world,
                        pivot,
                    ),
                    j,
                ),
            );
            work = work + dot(j, support_velocity);
            defect = defect + dot(j, sub(midpoint, support_velocity));
        }
        let support_impulse_world_ns = value(impulse);
        let support_angular_impulse_about_initial_pivot_world_nms = value(angular_impulse);
        let boundary_work_j = work.to_f64();
        let constraint_work_defect_j = defect.to_f64();
        if !support_impulse_world_ns.is_finite()
            || !support_angular_impulse_about_initial_pivot_world_nms.is_finite()
            || !boundary_work_j.is_finite()
            || !constraint_work_defect_j.is_finite()
        {
            return Err(PointPlaneStepError::Invalid(
                "rotating support account is unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RotatingPlaneStepReceipt {
            step,
            plane_after,
            motion_after: PlaneRigidMotion {
                pivot_world_m: end_pivot,
                ..motion
            },
            support_impulse_world_ns,
            support_angular_impulse_about_initial_pivot_world_nms,
            boundary_work_j,
            constraint_work_defect_j,
            prescribed_final_gap_m,
            prescribed_final_normal_speed_mps,
        })
    }
}

impl PlaneSpeed {
    pub(super) fn velocity_at(self, state: RigidBodyState, point: PointKinematics) -> [Dd; 3] {
        let Self::Rigid {
            pivot,
            velocity,
            angular_velocity,
        } = self
        else {
            unreachable!("a physical velocity field is requested only for rigid plane motion")
        };
        add(
            dd_vec(velocity),
            cross(
                dd_vec(angular_velocity),
                point_offset(state, point.arm_world, pivot),
            ),
        )
    }
}

pub(super) fn dd_vec(v: Vec3) -> [Dd; 3] {
    [v.x, v.y, v.z].map(Dd::from_f64)
}
fn value(v: [Dd; 3]) -> Vec3 {
    Vec3::new(v[0].to_f64(), v[1].to_f64(), v[2].to_f64())
}
fn add(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|k| a[k] + b[k])
}
fn sub(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|k| a[k] - b[k])
}
fn scale(a: [Dd; 3], s: Dd) -> [Dd; 3] {
    a.map(|v| v * s)
}
fn dot(a: [Dd; 3], b: [Dd; 3]) -> Dd {
    (0..3).fold(Dd::ZERO, |s, k| s + a[k] * b[k])
}
fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(rotation: UnitQuaternion, v: [Dd; 3]) -> [Dd; 3] {
    add(v, rotation_increment(rotation, v))
}
fn rotation_increment(rotation: UnitQuaternion, v: [Dd; 3]) -> [Dd; 3] {
    let [w, x, y, z] = rotation.components();
    let q = [x, y, z].map(Dd::from_f64);
    let t = scale(cross(q, v), Dd::from_f64(2.0));
    add(scale(t, Dd::from_f64(w)), cross(q, t))
}
pub(super) fn relative_speed(normal: [Dd; 3], velocity: [Dd; 3], support: [Dd; 3]) -> f64 {
    dot(normal, sub(velocity, support)).to_f64()
}

// Form relative geometry before adding a small body arm to a remote world
// origin. The original PointKinematics supplies these native components even
// when its convenient absolute point/velocity values have rounded them away.
fn point_offset(state: RigidBodyState, arm: Vec3, origin: [Dd; 3]) -> [Dd; 3] {
    add(
        sub(dd_vec(state.pose().position_world()), origin),
        dd_vec(arm),
    )
}
pub(super) fn point_velocity(point: PointKinematics) -> [Dd; 3] {
    add(
        dd_vec(point.center_of_mass_velocity_world),
        cross(
            dd_vec(point.angular_velocity_world),
            dd_vec(point.arm_world),
        ),
    )
}
pub(super) fn gap(normal: [Dd; 3], state: RigidBodyState, arm: Vec3, origin: [Dd; 3]) -> f64 {
    dot(normal, point_offset(state, arm, origin)).to_f64()
}
