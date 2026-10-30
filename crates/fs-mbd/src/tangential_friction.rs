//! Finite-patch tangential impulses on the original rigid-body state.
//!
//! Masses, full body inertias, poses and actual contact arms produce the full
//! shared-body mobility. No surrogate effective mass is accepted. This is a
//! frozen-geometry impulse phase: poses stay fixed, and the surrounding motion
//! integrator owns configuration advance and normal contact/impact admission.

/// Material-derived rolling loss with finite-body recoil and feasible no slip.
pub mod rolling;

/// Joint Hertz sliding/spin friction with full two-body recoil.
pub mod hertz;

/// Simultaneous Hertz shear/spin and material-derived rolling on one native pair.
pub mod hertz_rolling;

use crate::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_contact::tangential::nonsmooth::{
    NonsmoothAudit, NonsmoothContact, NonsmoothContactResult, NonsmoothControls,
    NonsmoothTangentialError, NonsmoothTangentialState, TangentialMobility,
};
use fs_exec::CancelGate;

/// Physical body pair corresponding to one supplied normal/tangential patch.
#[derive(Debug, Clone, Copy)]
pub struct RigidContactBinding {
    pub body_a: usize,
    pub body_b: usize,
    /// Both arms are relative to their body's COM in its declared body frame.
    pub arm_a_body: Vec3,
    pub arm_b_body: Vec3,
}

/// Independent dimensional admission and whole-body balance gates.
#[derive(Debug, Clone, Copy)]
pub struct RigidFrictionControls {
    pub maximum_bodies: usize,
    pub maximum_contacts: usize,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_mps: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
    pub contact: NonsmoothControls,
}

/// Actual body observables, in addition to the constitutive/solver receipt.
#[derive(Debug, Clone)]
pub struct RigidFrictionReport {
    pub contacts: Vec<NonsmoothContactResult>,
    pub friction: NonsmoothAudit,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub endpoint_velocity_residual_mps: f64,
}

/// Every refusal preserves both the original body slice and friction history.
#[derive(Debug)]
pub enum RigidFrictionError {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Contact(NonsmoothTangentialError),
    Cancelled,
    Audit(RigidFrictionReport),
}
impl From<DynamicsError> for RigidFrictionError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<NonsmoothTangentialError> for RigidFrictionError {
    fn from(error: NonsmoothTangentialError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for RigidFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid-body tangential friction: {self:?}")
    }
}
impl std::error::Error for RigidFrictionError {}

/// Solve coupled tangential motion and accept actual body momenta/history together.
///
/// Contacts/bindings correspond in input order; shared bodies are fully coupled.
/// The supplied points must coincide within the explicit geometry tolerance;
/// they are never moved to a fabricated common point. Surface normal patches
/// remain the normal owner's responsibility. Common material spin is derived
/// from the accepted body's actual angular velocities, replacing the supplied
/// contact's generic common-spin input for this frozen-geometry phase.
///
/// No normal impulses, position repair, restitution, or pose integration occur.
/// A changing-coefficient driver must separately balance the reported signed
/// coefficient-work port; that port and numerical losses are not friction heat.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_friction(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut NonsmoothTangentialState,
    contacts: &[NonsmoothContact<'_>],
    bindings: &[RigidContactBinding],
    dt_s: f64,
    controls: RigidFrictionControls,
    gate: &CancelGate,
) -> Result<RigidFrictionReport, RigidFrictionError> {
    poll(gate)?;
    let dofs = bodies
        .len()
        .checked_mul(6)
        .ok_or(RigidFrictionError::Invalid("body dimension overflow"))?;
    let rows = contacts
        .len()
        .checked_mul(3)
        .ok_or(RigidFrictionError::Invalid("contact dimension overflow"))?;
    let entries = rows
        .checked_mul(dofs)
        .ok_or(RigidFrictionError::Invalid("mobility dimension overflow"))?;
    if bodies.len() != properties.len()
        || contacts.len() != bindings.len()
        || contacts.len() > controls.maximum_contacts
        || bodies.len() > controls.maximum_bodies
        || dofs > controls.contact.maximum_mechanical_dofs
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_mps,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(RigidFrictionError::Invalid(
            "matching bounded bodies/contacts and positive finite gates required",
        ));
    }
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_by(|&a, &b| {
        contacts[a]
            .patch
            .patch_id()
            .cmp(contacts[b].patch.patch_id())
    });
    let mut factor = vec![0.0; entries];
    let mut free_velocity = Vec::with_capacity(contacts.len());
    let mut resolved = Vec::with_capacity(contacts.len());
    for (slot, &index) in order.iter().enumerate() {
        poll(gate)?;
        let binding = bindings[index];
        let contact = &contacts[index];
        if binding.body_a >= bodies.len()
            || binding.body_b >= bodies.len()
            || binding.body_a == binding.body_b
        {
            return Err(RigidFrictionError::Invalid(
                "two distinct valid body indices required",
            ));
        }
        let a = bodies[binding.body_a]
            .point_kinematics(properties[binding.body_a], binding.arm_a_body)?;
        let b = bodies[binding.body_b]
            .point_kinematics(properties[binding.body_b], binding.arm_b_body)?;
        if magnitude(a.point_world.sub(b.point_world)) > controls.geometry_tolerance_m {
            return Err(RigidFrictionError::Invalid(
                "contact points do not coincide",
            ));
        }
        let axes = [
            vector(contact.frame.longitudinal()),
            vector(contact.frame.lateral()),
            vector(contact.frame.normal()),
        ];
        let relative = a.point_velocity_world.sub(b.point_velocity_world);
        let relative_spin = a.angular_velocity_world.sub(b.angular_velocity_world);
        free_velocity.push([
            relative.dot(axes[0]),
            relative.dot(axes[1]),
            contact.torsional_length_m * relative_spin.dot(axes[2]),
        ]);
        resolved.push(NonsmoothContact {
            law: contact.law,
            patch: contact.patch,
            interface: contact.interface,
            frame: contact.frame,
            common_spin_rad_per_s: a
                .angular_velocity_world
                .scale(0.5)
                .add(b.angular_velocity_world.scale(0.5))
                .dot(axes[2]),
            torsional_length_m: contact.torsional_length_m,
            creepage_reference_speed_mps: contact.creepage_reference_speed_mps,
            work: contact.work.clone(),
        });
        for (body, arm, sign) in [
            (binding.body_a, binding.arm_a_body, 1.0),
            (binding.body_b, binding.arm_b_body, -1.0),
        ] {
            let inverse_mass_root = properties[body].mass().sqrt().recip();
            let orientation = bodies[body].pose().orientation();
            for axis in 0..3 {
                let direction = orientation.rotate_world_to_body(axes[axis]);
                let translation = if axis < 2 {
                    components(axes[axis].scale(inverse_mass_root * sign))
                } else {
                    [0.0; 3]
                };
                let rotation = if axis < 2 {
                    arm.cross(direction).scale(sign)
                } else {
                    direction.scale(sign * contact.torsional_length_m)
                };
                let rotation = components(properties[body].angular_mobility_factor(rotation));
                let offset = (3 * slot + axis) * dofs + 6 * body;
                for j in 0..3 {
                    factor[offset + j] = translation[j];
                    factor[offset + 3 + j] = rotation[j];
                }
            }
        }
    }
    let trial = history.prepare(
        &resolved,
        &free_velocity,
        TangentialMobility {
            mechanical_dofs: dofs,
            factor: &factor,
        },
        dt_s,
        controls.contact,
        gate,
    )?;
    let mut next = bodies.to_vec();
    for (result, &index) in trial.contacts().iter().zip(&order) {
        poll(gate)?;
        let binding = bindings[index];
        let impulse = vector(result.action.force_n).scale(dt_s);
        let couple = vector(result.action.torque_nm).scale(dt_s);
        for (body, arm, sign) in [
            (binding.body_a, binding.arm_a_body, 1.0),
            (binding.body_b, binding.arm_b_body, -1.0),
        ] {
            let after_force = next[body]
                .apply_impulse_at_body_point(properties[body], arm, impulse.scale(sign))?
                .state_after;
            let angular = after_force
                .pose()
                .orientation()
                .rotate_world_to_body(couple.scale(sign));
            next[body] = RigidBodyState::new(
                after_force.pose(),
                after_force.linear_momentum_world(),
                after_force.angular_momentum_body().add(angular),
            )?;
        }
    }
    let (energy_before, linear_before, angular_before) = totals(bodies, properties, gate)?;
    let (energy_after, linear_after, angular_after) = totals(&next, properties, gate)?;
    let mut midpoint_work = 0.0;
    let mut endpoint_residual = 0.0_f64;
    for ((result, contact), &index) in trial.contacts().iter().zip(&resolved).zip(&order) {
        poll(gate)?;
        let binding = bindings[index];
        let before_a = bodies[binding.body_a]
            .point_kinematics(properties[binding.body_a], binding.arm_a_body)?;
        let before_b = bodies[binding.body_b]
            .point_kinematics(properties[binding.body_b], binding.arm_b_body)?;
        let after_a = next[binding.body_a]
            .point_kinematics(properties[binding.body_a], binding.arm_a_body)?;
        let after_b = next[binding.body_b]
            .point_kinematics(properties[binding.body_b], binding.arm_b_body)?;
        let before_relative = before_a
            .point_velocity_world
            .sub(before_b.point_velocity_world);
        let after_relative = after_a
            .point_velocity_world
            .sub(after_b.point_velocity_world);
        let before_spin = before_a
            .angular_velocity_world
            .sub(before_b.angular_velocity_world);
        let after_spin = after_a
            .angular_velocity_world
            .sub(after_b.angular_velocity_world);
        midpoint_work += vector(result.action.force_n)
            .scale(dt_s)
            .dot(before_relative.add(after_relative).scale(0.5))
            + vector(result.action.torque_nm)
                .scale(dt_s)
                .dot(before_spin.add(after_spin).scale(0.5));
        let actual = [
            after_relative.dot(vector(contact.frame.longitudinal())),
            after_relative.dot(vector(contact.frame.lateral())),
            contact.torsional_length_m * after_spin.dot(vector(contact.frame.normal())),
        ];
        for j in 0..3 {
            endpoint_residual =
                endpoint_residual.max((actual[j] - result.endpoint_velocity_mps[j]).abs());
        }
    }
    let friction = trial.audit();
    let kinetic_change = energy_after - energy_before;
    let report = RigidFrictionReport {
        contacts: trial.contacts().to_vec(),
        friction,
        kinetic_before_j: energy_before,
        kinetic_after_j: energy_after,
        midpoint_impulse_work_j: midpoint_work,
        work_energy_residual_j: kinetic_change - midpoint_work,
        total_energy_residual_j: kinetic_change
            + friction.storage_change_j
            + friction.physical_heat_j
            + friction.elastic_numerical_loss_j
            + friction.mechanical_numerical_loss_j
            - friction.coefficient_work_j,
        linear_momentum_change_ns: linear_after.sub(linear_before),
        angular_momentum_change_nms: angular_after.sub(angular_before),
        endpoint_velocity_residual_mps: endpoint_residual,
    };
    if [
        energy_before,
        energy_after,
        midpoint_work,
        report.work_energy_residual_j,
        report.total_energy_residual_j,
        endpoint_residual,
        magnitude(report.linear_momentum_change_ns),
        magnitude(report.angular_momentum_change_nms),
    ]
    .iter()
    .any(|x| !x.is_finite())
        || endpoint_residual > controls.velocity_tolerance_mps
        || magnitude(report.linear_momentum_change_ns) > controls.linear_momentum_tolerance_ns
        || magnitude(report.angular_momentum_change_nms) > controls.angular_momentum_tolerance_nms
        || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
        || (kinetic_change - friction.mechanical_work_j).abs() > controls.energy_tolerance_j
        || report.total_energy_residual_j.abs() > controls.energy_tolerance_j
    {
        return Err(RigidFrictionError::Audit(report));
    }
    let mut next_history = history.clone();
    next_history.commit(trial)?;
    poll(gate)?;
    // No fallible operation follows this boundary.
    bodies.copy_from_slice(&next);
    *history = next_history;
    Ok(report)
}

fn totals(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    gate: &CancelGate,
) -> Result<(f64, Vec3, Vec3), RigidFrictionError> {
    let mut energy = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    for (&body, &property) in bodies.iter().zip(properties) {
        poll(gate)?;
        energy += body.kinetic_energy(property)?;
        linear = linear.add(body.linear_momentum_world());
        angular = angular
            .add(
                body.pose()
                    .position_world()
                    .cross(body.linear_momentum_world()),
            )
            .add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(body.angular_momentum_body()),
            );
    }
    Ok((energy, linear, angular))
}
fn vector(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn poll(gate: &CancelGate) -> Result<(), RigidFrictionError> {
    if gate.is_requested() {
        Err(RigidFrictionError::Cancelled)
    } else {
        Ok(())
    }
}
