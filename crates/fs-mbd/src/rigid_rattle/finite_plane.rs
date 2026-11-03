//! Bilateral point/plane motion with a finite, freely recoiling support body.
//!
//! The material point belongs to body 0. Body 1 carries an infinite plane;
//! its reaction acts at the current projected point on that plane, rather
//! than at a fixed support anchor. This retains the moving normal's angular
//! Jacobian and reciprocal torque while permitting frictionless sliding.

use super::distance::{DistanceStepControls, DistanceStepError};
use super::pair_step::{self, Constraint, Geometry};
use crate::{
    ActionReactionImpulseReceipt, MassProperties, PointKinematics, RigidBodyIntegrator,
    RigidBodyState, StepReceipt, Vec3, Wrench,
};
use fs_math::dd::Dd;

/// A material point of body 0 constrained to a plane carried by finite body 1.
/// Both body frames are centered at their actual centers of mass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinitePointPlaneConstraint {
    properties: [MassProperties; 2],
    point_body_m: Vec3,
    plane_offset_body_m: Dd,
    plane_normal_body: Vec3,
    position_tolerance_m: f64,
    velocity_tolerance_mps: f64,
}

/// Both native bodies and their internal constraint impulse/work account.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinitePlaneStepReceipt<FreeStep = StepReceipt> {
    /// Point body followed by the finite plane carrier.
    pub states_after: [RigidBodyState; 2],
    pub position_impulse: ActionReactionImpulseReceipt,
    pub free_steps: [FreeStep; 2],
    pub velocity_impulse: ActionReactionImpulseReceipt,
    /// Signed constraint impulse work, paired from relative native COM/spin
    /// velocities. This numerical integration defect is neither heat nor an
    /// external support-work supply. No body energy is changed to erase it.
    pub constraint_work_j: f64,
    /// Signed endpoint gap, positive on the plane normal's side [m].
    pub final_gap_m: f64,
    /// Endpoint material-point speed relative to the rotating carrier [m/s].
    pub final_normal_speed_mps: f64,
    pub position_iterations: u32,
}

impl FinitePointPlaneConstraint {
    /// Declare both finite inertias and body-local geometry. A nonzero normal
    /// is normalized once; positive absolute gap/speed tolerances are explicit.
    pub fn new(
        properties: [MassProperties; 2],
        point_body_m: Vec3,
        plane_point_body_m: Vec3,
        plane_normal_body: Vec3,
        position_tolerance_m: f64,
        velocity_tolerance_mps: f64,
    ) -> Result<Self, DistanceStepError> {
        if !point_body_m.is_finite()
            || !plane_point_body_m.is_finite()
            || !position_tolerance_m.is_finite()
            || position_tolerance_m <= 0.0
            || !velocity_tolerance_mps.is_finite()
            || velocity_tolerance_mps <= 0.0
        {
            return Err(DistanceStepError::Invalid(
                "finite point/plane geometry or tolerances",
            ));
        }
        let plane_normal_body = plane_normal_body.normalized("plane normal")?;
        let plane_offset_body_m = dot(vector(plane_normal_body), vector(plane_point_body_m));
        if !plane_offset_body_m.to_f64().is_finite() {
            return Err(DistanceStepError::Invalid("finite plane normal offset"));
        }
        Ok(Self {
            properties,
            point_body_m,
            plane_offset_body_m,
            plane_normal_body,
            position_tolerance_m,
            velocity_tolerance_mps,
        })
    }

    /// Advance both bodies using the original scalar pair constraint engine.
    ///
    /// The support's translation and rotation follow its own full finite
    /// inertia, external wrench and native free step. Normal reactions at the
    /// current projected point transfer reciprocal impulse and torque; no
    /// prescribed pivot, infinite-mass substitution or pose repair is used.
    /// The returned states can be passed directly to the next interval.
    ///
    /// Both initial constraints must already hold. Body rotation, plane-normal
    /// travel, impulses and local backtracking share `DistanceStepControls`.
    /// This is a smooth bilateral joint, permitting tensile reactions and
    /// arbitrary in-plane sliding. It adds no unilateral collision, friction,
    /// finite plane extent, or exact-energy/symplectic claim. No heap storage is
    /// allocated; cancellation or any late refusal publishes neither body.
    pub fn advance(
        &self,
        states: [RigidBodyState; 2],
        integrator: RigidBodyIntegrator,
        wrenches: [Wrench; 2],
        duration_s: f64,
        controls: DistanceStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FinitePlaneStepReceipt, DistanceStepError> {
        let step = pair_step::advance(
            self,
            states,
            integrator,
            wrenches,
            duration_s,
            controls,
            &mut cancelled,
        )?;
        finish(step, &mut cancelled)
    }

    /// Evolve a loaded finite carrier and sliding body with spatial moments.
    ///
    /// World-held point forces follow their actual rotating material arms in
    /// every trial; body and world free couples retain their declared frames.
    /// Both original spatial-load receipts expose external impulse and work
    /// separately from the reciprocal sliding reaction. The finite carrier,
    /// local branch, tolerances and cancellation contract are those of
    /// [`Self::advance`]; loads are held over this interval, with no extra heap
    /// buffers or alternate free-motion solver.
    pub fn advance_spatial(
        &self,
        states: [RigidBodyState; 2],
        integrator: RigidBodyIntegrator,
        loads: [crate::spatial_loads::PreparedSpatialLoads; 2],
        duration_s: f64,
        controls: DistanceStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FinitePlaneStepReceipt<crate::spatial_loads::SpatialStepReceipt>, DistanceStepError>
    {
        let step = pair_step::advance_spatial(
            self,
            states,
            integrator,
            loads,
            duration_s,
            controls,
            &mut cancelled,
        )?;
        finish(step, &mut cancelled)
    }
}

fn finish<FreeStep>(
    step: super::distance::DistanceStepReceipt<FreeStep>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FinitePlaneStepReceipt<FreeStep>, DistanceStepError> {
    // Internal work must cancel common translation before local spin is
    // added: a rounded laboratory point velocity can lose that spin.
    let constraint_work_j =
        (relative_work(step.position_impulse) + relative_work(step.velocity_impulse)).to_f64();
    if !constraint_work_j.is_finite() {
        return Err(DistanceStepError::Invalid("finite plane constraint work"));
    }
    if cancelled() {
        return Err(DistanceStepError::Cancelled);
    }
    Ok(FinitePlaneStepReceipt {
        states_after: step.states_after,
        position_impulse: step.position_impulse,
        free_steps: step.free_steps,
        velocity_impulse: step.velocity_impulse,
        constraint_work_j,
        final_gap_m: step.final_gap_m,
        final_normal_speed_mps: step.final_range_rate_mps,
        position_iterations: step.position_iterations,
    })
}

impl Constraint for FinitePointPlaneConstraint {
    fn properties(&self) -> [MassProperties; 2] {
        self.properties
    }
    fn tolerances(&self) -> [f64; 2] {
        [self.position_tolerance_m, self.velocity_tolerance_mps]
    }

    fn geometry(&self, states: [RigidBodyState; 2]) -> Result<Geometry, DistanceStepError> {
        let point = states[0].point_kinematics(self.properties[0], self.point_body_m)?;
        let carrier = states[1].point_kinematics(self.properties[1], Vec3::ZERO)?;
        let normal = states[1]
            .pose()
            .orientation()
            .rotate_body_to_world(self.plane_normal_body);
        let n = vector(normal);
        // Keep the physical material arms at distant world origins.
        let from_carrier = add(
            sub(
                vector(states[0].pose().position_world()),
                vector(states[1].pose().position_world()),
            ),
            vector(point.arm_world),
        );
        // An infinite plane depends only on its body-normal offset. A huge
        // tangential anchor must not create a false offset through rounding
        // when it is rotated into world coordinates, or overflow a query for
        // a physically irrelevant material point far along the same plane.
        let gap = dot(n, from_carrier) - self.plane_offset_body_m;
        let projected = sub(from_carrier, scale(n, gap));
        let projected_value = value(projected);
        if !projected_value.is_finite() || !normal.is_finite() || !gap.to_f64().is_finite() {
            return Err(DistanceStepError::Invalid("finite plane relative geometry"));
        }
        let carrier_arm = states[1]
            .pose()
            .orientation()
            .rotate_world_to_body(projected_value);
        let represented_arm = states[1]
            .pose()
            .orientation()
            .rotate_body_to_world(carrier_arm);
        let arm_error = value(sub(vector(represented_arm), projected));
        if !carrier_arm.is_finite()
            || !arm_error.is_finite()
            || arm_error.x.hypot(arm_error.y).hypot(arm_error.z) > self.position_tolerance_m
        {
            return Err(DistanceStepError::Invalid(
                "finite plane reaction arm is unrepresentable",
            ));
        }
        // The plane's angular velocity acts on the CURRENT sliding point.
        // Using the fixed plane anchor here would omit its rotational Jacobian.
        let relative = add(
            sub(
                vector(point.center_of_mass_velocity_world),
                vector(carrier.center_of_mass_velocity_world),
            ),
            sub(
                cross(
                    vector(point.angular_velocity_world),
                    vector(point.arm_world),
                ),
                cross(vector(carrier.angular_velocity_world), from_carrier),
            ),
        );
        let rate = dot(n, relative).to_f64();
        if !rate.is_finite() {
            return Err(DistanceStepError::Invalid("finite plane relative speed"));
        }
        Ok(Geometry {
            // Shared convention: normal impulse on body 1 increases the gap.
            normal: normal.scale(-1.0),
            arms_body_m: [self.point_body_m, carrier_arm],
            gap: gap.to_f64(),
            rate,
        })
    }
}

fn relative_work(event: ActionReactionImpulseReceipt) -> Dd {
    let before = relative_velocity(
        event.first.point_kinematics_before,
        event.second.point_kinematics_before,
    );
    let after = relative_velocity(
        event.first.point_kinematics_after,
        event.second.point_kinematics_after,
    );
    dot(
        vector(event.second.impulse_world),
        scale(add(before, after), Dd::from_f64(0.5)),
    )
}
fn relative_velocity(a: PointKinematics, b: PointKinematics) -> [Dd; 3] {
    add(
        sub(
            vector(b.center_of_mass_velocity_world),
            vector(a.center_of_mass_velocity_world),
        ),
        sub(
            cross(vector(b.angular_velocity_world), vector(b.arm_world)),
            cross(vector(a.angular_velocity_world), vector(a.arm_world)),
        ),
    )
}
fn vector(v: Vec3) -> [Dd; 3] {
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
