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

/// Sourced real-area adhesive sliding and twisting with finite-body recoil.
pub mod adhesive;

// Law-neutral native mechanics shared by simultaneous circular traction owners.
mod island;

/// Simultaneous Hertz shear/spin and material-derived rolling on one native pair.
pub mod hertz_rolling;

/// Source-backed changing-load Hertz--Mindlin contact with native normal/shear recoil.
pub mod mindlin;

/// Source-capable directional Coulomb friction with full two-body mobility.
pub mod directional;

use crate::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_contact::tangential::nonsmooth::{
    NonsmoothAudit, NonsmoothContact, NonsmoothContactResult, NonsmoothControls,
    NonsmoothTangentialError, NonsmoothTangentialState, TangentialMobility,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

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
    /// Actual native momentum increments contracted with midpoint velocities.
    /// Retains small losses beside unchanged common translational motion.
    pub kinetic_change_j: f64,
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
    let mut resultants = vec![ResultantImpulse::default(); bodies.len()];
    for (result, &index) in trial.contacts().iter().zip(&order) {
        poll(gate)?;
        let binding = bindings[index];
        let impulse = vector(result.action.force_n).scale(dt_s);
        let couple = vector(result.action.torque_nm).scale(dt_s);
        for (body, arm, sign) in [
            (binding.body_a, binding.arm_a_body, 1.0),
            (binding.body_b, binding.arm_b_body, -1.0),
        ] {
            resultants[body].add(bodies[body], arm, impulse.scale(sign), couple.scale(sign));
        }
    }
    let native = apply_resultants(bodies, properties, &resultants, gate)?;
    let next = &native.bodies_after;
    let energy_before = native.kinetic_before_j;
    let energy_after = native.kinetic_after_j;
    let mut midpoint_work = Dd::from_f64(0.0);
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
        midpoint_work = midpoint_work
            + Dd::from_f64(
                vector(result.action.force_n)
                    .scale(dt_s)
                    .dot(before_relative.add(after_relative).scale(0.5)),
            )
            + Dd::from_f64(
                vector(result.action.torque_nm)
                    .scale(dt_s)
                    .dot(before_spin.add(after_spin).scale(0.5)),
            );
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
    let midpoint_work = midpoint_work.to_f64();
    let friction = trial.audit();
    let kinetic_change = native.kinetic_change_j;
    let report = RigidFrictionReport {
        contacts: trial.contacts().to_vec(),
        friction,
        kinetic_before_j: energy_before,
        kinetic_after_j: energy_after,
        kinetic_change_j: kinetic_change,
        midpoint_impulse_work_j: midpoint_work,
        work_energy_residual_j: kinetic_change - midpoint_work,
        total_energy_residual_j: kinetic_change
            + friction.storage_change_j
            + friction.physical_heat_j
            + friction.elastic_numerical_loss_j
            + friction.mechanical_numerical_loss_j
            - friction.coefficient_work_j,
        linear_momentum_change_ns: native.linear_momentum_change_ns,
        angular_momentum_change_nms: native.angular_momentum_change_nms,
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
    bodies.copy_from_slice(next);
    *history = next_history;
    Ok(report)
}

/// All simultaneous contact covectors are summed before rounding native state.
/// Couples retain their original body-frame central moment independently of
/// the point-force arm; no contact point is relocated.
#[derive(Clone, Copy, Default)]
struct ResultantImpulse {
    linear: [Dd; 3],
    angular: [Dd; 3],
}

impl ResultantImpulse {
    fn add(&mut self, body: RigidBodyState, arm: Vec3, impulse: Vec3, couple: Vec3) {
        let orientation = body.pose().orientation();
        add_components(&mut self.linear, impulse);
        add_components(
            &mut self.angular,
            arm.cross(orientation.rotate_world_to_body(impulse)),
        );
        add_components(&mut self.angular, orientation.rotate_world_to_body(couple));
    }
}

struct NativeImpulsePhase {
    bodies_after: Vec<RigidBodyState>,
    kinetic_before_j: f64,
    kinetic_after_j: f64,
    kinetic_change_j: f64,
    linear_momentum_change_ns: Vec3,
    angular_momentum_change_nms: Vec3,
}

fn apply_resultants(
    before: &[RigidBodyState],
    properties: &[MassProperties],
    resultants: &[ResultantImpulse],
    gate: &CancelGate,
) -> Result<NativeImpulsePhase, RigidFrictionError> {
    let mut after = Vec::with_capacity(before.len());
    let mut energy_before = Dd::from_f64(0.0);
    let mut energy_after = Dd::from_f64(0.0);
    let mut kinetic_change = Dd::from_f64(0.0);
    let mut linear = [Dd::from_f64(0.0); 3];
    let mut angular = [Dd::from_f64(0.0); 3];
    for ((&body, &property), impulse) in before.iter().zip(properties).zip(resultants) {
        poll(gate)?;
        let receipt = body.apply_momentum_impulse(
            property,
            rounded_components(impulse.linear),
            rounded_components(impulse.angular),
        )?;
        energy_before = energy_before + Dd::from_f64(receipt.work.kinetic_energy_before_j);
        energy_after = energy_after + Dd::from_f64(receipt.work.kinetic_energy_after_j);
        kinetic_change = kinetic_change + Dd::from_f64(receipt.work.kinetic_energy_change_j);
        add_components(&mut linear, receipt.linear_increment_world);
        add_components(
            &mut angular,
            body.pose()
                .position_world()
                .cross(receipt.linear_increment_world),
        );
        add_components(
            &mut angular,
            body.pose()
                .orientation()
                .rotate_body_to_world(receipt.angular_increment_body),
        );
        after.push(receipt.state_after);
    }
    Ok(NativeImpulsePhase {
        bodies_after: after,
        kinetic_before_j: energy_before.to_f64(),
        kinetic_after_j: energy_after.to_f64(),
        kinetic_change_j: kinetic_change.to_f64(),
        linear_momentum_change_ns: rounded_components(linear),
        angular_momentum_change_nms: rounded_components(angular),
    })
}

fn add_components(sum: &mut [Dd; 3], value: Vec3) {
    for (entry, component) in sum.iter_mut().zip(components(value)) {
        *entry = *entry + Dd::from_f64(component);
    }
}

fn rounded_components(sum: [Dd; 3]) -> Vec3 {
    Vec3::new(sum[0].to_f64(), sum[1].to_f64(), sum[2].to_f64())
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

#[cfg(test)]
mod impulse_tests {
    use super::*;
    use crate::Pose;

    #[test]
    fn simultaneous_opposing_wrenches_preserve_small_momentum_and_net_load() {
        let property = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 2.0)).unwrap();
        let before = RigidBodyState::new(
            Pose::identity(),
            Vec3::new(0.125, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.125),
        )
        .unwrap();
        let large = 2.0_f64.powi(54);
        // Three simultaneous contacts at the same physical arm: the two large
        // ones cancel exactly, leaving the 1/4 Ns load. Sequentially updating
        // native momenta loses the original 1/8 Ns; summing in f64 loses 1/4 Ns.
        for loads in [[large, 0.25, -large], [-large, 0.25, large]] {
            let mut impulses = [ResultantImpulse::default()];
            for load in loads {
                impulses[0].add(
                    before,
                    Vec3::new(0.0, -1.0, 0.0),
                    Vec3::new(load, 0.0, 0.0),
                    Vec3::ZERO,
                );
            }
            let result =
                apply_resultants(&[before], &[property], &impulses, &CancelGate::new()).unwrap();
            let after = result.bodies_after[0];
            assert_eq!(after.linear_momentum_world(), Vec3::new(0.375, 0.0, 0.0));
            assert_eq!(after.angular_momentum_body(), Vec3::new(0.0, 0.0, 0.375));
            assert_eq!(result.linear_momentum_change_ns, Vec3::new(0.25, 0.0, 0.0));
            assert_eq!(
                result.angular_momentum_change_nms,
                Vec3::new(0.0, 0.0, 0.25)
            );
            // Independent quadratic kinetic-energy difference, for m=Izz=2.
            assert_eq!(result.kinetic_change_j, 0.0625);
        }
    }
}
