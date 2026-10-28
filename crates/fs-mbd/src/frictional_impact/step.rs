//! Configuration-admitted first-order force/contact/drift on native rigid bodies.
//!
//! The collision owner must admit the entire proposed geometric path and feature
//! coverage. Independent initial/endpoint gap checks do not constitute CCD.
//! Forces kick momenta, simultaneous normal/Coulomb contact resolves velocities,
//! and the original torque-free midpoint integrator advances the real poses.
//! No gap projection, fabricated contact, or impact-to-heat allocation occurs.

use super::{
    RigidFrictionalContact, RigidFrictionalImpactError, RigidFrictionalImpactReport,
    RigidImpactControls, RigidImpactCounterpart, apply_rigid_frictional_impacts,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};
use fs_exec::CancelGate;

/// Independent phase and complete-step acceptance gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidContactStepControls {
    pub impact: RigidImpactControls,
    pub drift_energy_tolerance_j: f64,
    pub drift_angular_momentum_tolerance_nms: f64,
    pub work_energy_tolerance_j: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
}

/// A signed configuration gap from the actual material points, ordered by ID.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidContactGap {
    pub id: u64,
    pub signed_gap_m: f64,
}

/// The ORIGINAL proposed drift path, supplied to the required geometry owner.
///
/// The callback must certify that no omitted feature collides and that supplied
/// witnesses/normals remain applicable over the entire interval. Checking only a
/// finite set of samples cannot in general establish this. For planar translation
/// the caller can prove a linear gap bound; curved/rotating features need their
/// original collision/continuous-path authority.
pub struct RigidContactPath<'a> {
    pub initial: &'a [RigidBodyState],
    pub post_contact: &'a [RigidBodyState],
    pub endpoint: &'a [RigidBodyState],
    pub properties: &'a [MassProperties],
    pub contacts: &'a [RigidFrictionalContact],
    pub duration_s: f64,
    gate: &'a CancelGate,
}
impl RigidContactPath<'_> {
    /// Evaluate the original zero-force midpoint drift at a fraction in [0,1].
    /// This defines the candidate numerical path, not exact Euler-top dynamics.
    /// Sampling is a tool for the geometry owner, not a coverage certificate.
    pub fn sample(&self, fraction: f64) -> Result<Vec<RigidBodyState>, RigidFrictionalImpactError> {
        super::poll(self.gate)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(RigidFrictionalImpactError::Invalid(
                "path fraction outside [0,1]",
            ));
        }
        if fraction == 0.0 {
            return Ok(self.post_contact.to_vec());
        }
        if fraction == 1.0 {
            return Ok(self.endpoint.to_vec());
        }
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        self.post_contact
            .iter()
            .zip(self.properties)
            .map(|(&body, &property)| {
                super::poll(self.gate)?;
                Ok(integrator
                    .step(body, property, Wrench::ZERO, fraction * self.duration_s)?
                    .state_after)
            })
            .collect()
    }
}

/// No refusal publishes momenta or poses, including a late geometry refusal.
#[derive(Debug)]
pub enum RigidContactStepError<E> {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Impact(RigidFrictionalImpactError),
    /// Original geometry owner rejected the whole path or its feature coverage.
    Path(E),
    Penetration {
        contact_id: u64,
        endpoint: bool,
        signed_gap_m: f64,
    },
    Cancelled,
    Audit(&'static str),
}
impl<E> From<DynamicsError> for RigidContactStepError<E> {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl<E> From<RigidFrictionalImpactError> for RigidContactStepError<E> {
    fn from(error: RigidFrictionalImpactError) -> Self {
        Self::Impact(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidContactStepError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid contact time step: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidContactStepError<E> {}

/// Work and configuration from the original accepted native states.
pub struct RigidContactStepReport<P> {
    pub duration_s: f64,
    pub impact: RigidFrictionalImpactReport,
    pub path_admission: P,
    pub initial_gaps: Vec<RigidContactGap>,
    pub endpoint_gaps: Vec<RigidContactGap>,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Force impulse dotted with actual COM midpoint velocity, including gravity.
    /// This is the declared kick's work, not force times the later displacement.
    pub force_impulse_work_j: f64,
    /// Exact work of the constant total world force along accepted COM displacement.
    pub force_displacement_work_j: f64,
    /// Signed kick work minus force-displacement work. This first-order splitting
    /// defect is reported for refinement; it is not physical friction heat.
    pub force_splitting_work_j: f64,
    /// Applied body torque impulse dotted with actual body midpoint angular speed.
    pub torque_impulse_work_j: f64,
    pub kick_work_energy_residual_j: f64,
    /// Measured kinetic change during the original zero-force midpoint drift.
    pub drift_energy_defect_j: f64,
    pub drift_angular_momentum_defect_nms: Vec3,
    /// Final kinetic change minus kick/support work plus total impact loss.
    /// Integration and normal-impact losses are never called friction heat.
    pub total_energy_residual_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// Advance actual native poses/momenta under a caller-admitted contact path.
///
/// Each wrench is held constant in its declared world-force/body-torque frame
/// for the initial kick; uniform gravity contributes its original mass times g.
/// The supplied normal law remains explicit. A sustained resting contact
/// ordinarily uses zero restitution/target, rather than creating repeated bounce.
///
/// Plane witnesses translate at their supplied constant velocity. Dynamic
/// witnesses and normals are frozen feature descriptions for this interval;
/// the geometry callback must prove their continued applicability. The step
/// itself rejects every negative supplied initial/endpoint gap, with no repair.
///
/// This is a first-order kick/contact/drift composition. Its impact solve keeps
/// its explicit impulse/velocity determinacy scope and total unallocated loss. It is
/// not general collision discovery, a full Moreau-Jean method, or a guarantee
/// of rigid-contact existence. Timestep refinement remains necessary.
///
/// The mandatory admission callback runs on staged states and must not publish
/// external history. Its receipt publishes only with the accepted native bodies.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn step_rigid_frictional_contacts<P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    wrenches: &[Wrench],
    gravity: Gravity,
    contacts: &[RigidFrictionalContact],
    dt_s: f64,
    controls: RigidContactStepControls,
    admit_path: impl FnOnce(&RigidContactPath<'_>) -> Result<P, E>,
    gate: &CancelGate,
) -> Result<RigidContactStepReport<P>, RigidContactStepError<E>> {
    poll(gate)?;
    if bodies.len() != properties.len()
        || bodies.len() != wrenches.len()
        || bodies.len() > controls.impact.maximum_bodies
        || contacts.len() > controls.impact.maximum_contacts
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || [
            controls.drift_energy_tolerance_j,
            controls.drift_angular_momentum_tolerance_nms,
            controls.work_energy_tolerance_j,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(RigidContactStepError::Invalid(
            "matching bounded native bodies and positive finite step gates required",
        ));
    }
    let initial_gaps = gaps(bodies, properties, contacts, 0.0, false, gate)?;
    let (before_energy, before_linear, before_angular) = super::totals(bodies, properties, gate)?;
    let mut kicked = Vec::with_capacity(bodies.len());
    let mut force_work = 0.0;
    let mut torque_work = 0.0;
    let mut applied_linear = Vec3::ZERO;
    let mut applied_angular = Vec3::ZERO;
    for ((&body, &property), &wrench) in bodies.iter().zip(properties).zip(wrenches) {
        poll(gate)?;
        wrench.validate()?;
        let force = wrench
            .force_world
            .add(gravity.acceleration_world().scale(property.mass()));
        // Original native force-impulse owner, at the real COM.
        let receipt = body.apply_force_at_body_point(property, Vec3::ZERO, force, dt_s)?;
        let after_force = receipt.impulse.state_after;
        let torque_impulse = wrench.torque_body.scale(dt_s);
        let next = RigidBodyState::new(
            after_force.pose(),
            after_force.linear_momentum_world(),
            after_force.angular_momentum_body().add(torque_impulse),
        )?;
        let omega_before =
            property.angular_velocity_body_checked(after_force.angular_momentum_body())?;
        let omega_after = property.angular_velocity_body_checked(next.angular_momentum_body())?;
        force_work += receipt.impulse.work.impulse_work_j;
        torque_work += torque_impulse.dot(omega_before.add(omega_after).scale(0.5));
        applied_linear = applied_linear.add(receipt.impulse.impulse_world);
        applied_angular = applied_angular
            .add(
                body.pose()
                    .position_world()
                    .cross(receipt.impulse.impulse_world),
            )
            .add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(torque_impulse),
            );
        kicked.push(next);
    }
    let (kick_energy, _, _) = super::totals(&kicked, properties, gate)?;
    let kick_residual = kick_energy - before_energy - force_work - torque_work;
    if !kick_residual.is_finite() || kick_residual.abs() > controls.work_energy_tolerance_j {
        return Err(RigidContactStepError::Audit(
            "native force/torque impulse work failed",
        ));
    }
    let impact =
        apply_rigid_frictional_impacts(&mut kicked, properties, contacts, controls.impact, gate)?;
    let (_, drift_linear_before, drift_angular_before) = super::totals(&kicked, properties, gate)?;
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut endpoint = Vec::with_capacity(kicked.len());
    for (&body, &property) in kicked.iter().zip(properties) {
        poll(gate)?;
        endpoint.push(
            integrator
                .step(body, property, Wrench::ZERO, dt_s)?
                .state_after,
        );
    }
    let endpoint_gaps = gaps(&endpoint, properties, contacts, dt_s, true, gate)?;
    let (after_energy, after_linear, after_angular) = super::totals(&endpoint, properties, gate)?;
    let force_displacement_work: f64 = bodies
        .iter()
        .zip(&endpoint)
        .zip(properties)
        .zip(wrenches)
        .map(|(((&before, &after), &property), &wrench)| {
            wrench
                .force_world
                .add(gravity.acceleration_world().scale(property.mass()))
                .dot(
                    after
                        .pose()
                        .position_world()
                        .sub(before.pose().position_world()),
                )
        })
        .sum();
    let force_splitting_work = force_work - force_displacement_work;
    let drift_energy = after_energy - impact.kinetic_after_j;
    let drift_angular = after_angular.sub(drift_angular_before);
    let total_energy =
        after_energy - before_energy - force_work - torque_work - impact.support_work_j
            + impact.contact.dissipated_energy;
    let linear_residual = after_linear
        .sub(before_linear)
        .sub(applied_linear)
        .sub(impact.support_impulse_ns);
    let angular_residual = after_angular
        .sub(before_angular)
        .sub(applied_angular)
        .sub(impact.support_angular_impulse_nms);
    if [
        before_energy,
        after_energy,
        force_work,
        force_displacement_work,
        force_splitting_work,
        torque_work,
        drift_energy,
        total_energy,
        super::magnitude(drift_angular),
        super::magnitude(after_linear.sub(drift_linear_before)),
        super::magnitude(linear_residual),
        super::magnitude(angular_residual),
    ]
    .iter()
    .any(|x| !x.is_finite())
        || drift_energy.abs() > controls.drift_energy_tolerance_j
        || super::magnitude(drift_angular) > controls.drift_angular_momentum_tolerance_nms
        || super::magnitude(after_linear.sub(drift_linear_before))
            > controls.linear_momentum_tolerance_ns
        || total_energy.abs() > controls.work_energy_tolerance_j
        || super::magnitude(linear_residual) > controls.linear_momentum_tolerance_ns
        || super::magnitude(angular_residual) > controls.angular_momentum_tolerance_nms
    {
        return Err(RigidContactStepError::Audit(
            "native drift or complete force/contact/drift work failed",
        ));
    }
    let path_admission = admit_path(&RigidContactPath {
        initial: bodies,
        post_contact: &kicked,
        endpoint: &endpoint,
        properties,
        contacts,
        duration_s: dt_s,
        gate,
    })
    .map_err(RigidContactStepError::Path)?;
    poll(gate)?;
    // No fallible work or callback follows this publication boundary.
    bodies.copy_from_slice(&endpoint);
    Ok(RigidContactStepReport {
        duration_s: dt_s,
        impact,
        path_admission,
        initial_gaps,
        endpoint_gaps,
        kinetic_before_j: before_energy,
        kinetic_after_j: after_energy,
        force_impulse_work_j: force_work,
        force_displacement_work_j: force_displacement_work,
        force_splitting_work_j: force_splitting_work,
        torque_impulse_work_j: torque_work,
        kick_work_energy_residual_j: kick_residual,
        drift_energy_defect_j: drift_energy,
        drift_angular_momentum_defect_nms: drift_angular,
        total_energy_residual_j: total_energy,
        linear_momentum_residual_ns: linear_residual,
        angular_momentum_residual_nms: angular_residual,
    })
}

fn gaps<E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact],
    elapsed_s: f64,
    endpoint: bool,
    gate: &CancelGate,
) -> Result<Vec<RigidContactGap>, RigidContactStepError<E>> {
    let mut result = Vec::with_capacity(contacts.len());
    for contact in contacts {
        poll(gate)?;
        if contact.point.body >= bodies.len() {
            return Err(RigidContactStepError::Invalid(
                "invalid positive witness body",
            ));
        }
        let normal = contact.normal_world.normalized("step.contact.normal")?;
        let a = bodies[contact.point.body]
            .point_kinematics(properties[contact.point.body], contact.point.arm_body)?
            .point_world;
        let b = match contact.other {
            RigidImpactCounterpart::Body(point) => {
                if point.body >= bodies.len() || point.body == contact.point.body {
                    return Err(RigidContactStepError::Invalid(
                        "invalid counterpart witness body",
                    ));
                }
                bodies[point.body]
                    .point_kinematics(properties[point.body], point.arm_body)?
                    .point_world
            }
            RigidImpactCounterpart::PrescribedPlane {
                point_world,
                velocity_world,
            } => {
                if !point_world.is_finite() || !velocity_world.is_finite() {
                    return Err(RigidContactStepError::Invalid(
                        "nonfinite prescribed plane path",
                    ));
                }
                point_world.add(velocity_world.scale(elapsed_s))
            }
        };
        let gap = normal.dot(a.sub(b));
        if !gap.is_finite() {
            return Err(RigidContactStepError::Invalid(
                "nonfinite native contact gap",
            ));
        }
        if gap < 0.0 {
            return Err(RigidContactStepError::Penetration {
                contact_id: contact.id,
                endpoint,
                signed_gap_m: gap,
            });
        }
        result.push(RigidContactGap {
            id: contact.id,
            signed_gap_m: gap,
        });
    }
    result.sort_by_key(|gap| gap.id);
    Ok(result)
}
fn poll<E>(gate: &CancelGate) -> Result<(), RigidContactStepError<E>> {
    if gate.is_requested() {
        Err(RigidContactStepError::Cancelled)
    } else {
        Ok(())
    }
}
