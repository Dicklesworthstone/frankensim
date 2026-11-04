//! Configuration-admitted first-order force/contact/drift on native rigid bodies.
//!
//! The collision owner must admit the entire proposed geometric path and feature
//! coverage. Independent initial/endpoint gap checks do not constitute CCD.
//! Forces kick momenta, simultaneous normal/Coulomb contact resolves velocities,
//! and the original torque-free midpoint integrator advances the real poses.
//! No gap projection, fabricated contact, or impact-to-heat allocation occurs.

/// Directional coefficients and source histories through the same native step.
pub mod elliptic;

use super::{
    ImpactResponse, RigidFrictionalContact, RigidFrictionalImpactError,
    RigidFrictionalImpactReport, RigidImpactControls, RigidImpactCounterpart,
    apply_rigid_frictional_impacts,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};
use fs_contact::frictional_multicontact::{FrictionalContactError, FrictionalContactResult};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

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
pub struct RigidContactPath<'a, F = f64> {
    pub initial: &'a [RigidBodyState],
    pub post_contact: &'a [RigidBodyState],
    pub endpoint: &'a [RigidBodyState],
    pub properties: &'a [MassProperties],
    pub contacts: &'a [RigidFrictionalContact<F>],
    pub duration_s: f64,
    gate: &'a CancelGate,
}
impl<F> RigidContactPath<'_, F> {
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
pub enum RigidContactStepError<E, R = FrictionalContactResult, I = FrictionalContactError> {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Impact(RigidFrictionalImpactError<R, I>),
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
impl<E, R, I> From<DynamicsError> for RigidContactStepError<E, R, I> {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl<E, R, I> From<RigidFrictionalImpactError<R, I>> for RigidContactStepError<E, R, I> {
    fn from(error: RigidFrictionalImpactError<R, I>) -> Self {
        Self::Impact(error)
    }
}
impl<E: core::fmt::Debug, R: core::fmt::Debug, I: core::fmt::Debug> core::fmt::Display
    for RigidContactStepError<E, R, I>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid contact time step: {self:?}")
    }
}
impl<E: std::error::Error + 'static, R: core::fmt::Debug, I: core::fmt::Debug> std::error::Error
    for RigidContactStepError<E, R, I>
{
}

/// Work and configuration from the original accepted native states.
pub struct RigidContactStepReport<P, R = FrictionalContactResult> {
    pub duration_s: f64,
    pub impact: RigidFrictionalImpactReport<R>,
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
/// Momentum and energy gates compose actual native increments from each phase;
/// large common translation or orbital momentum cannot swallow a small contact
/// response by subtraction of scene totals. Unrepresentable applied impulses
/// and drift errors remain subject to the original absolute tolerances.
#[allow(clippy::too_many_arguments)]
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
    step_with(
        bodies,
        properties,
        wrenches,
        gravity,
        contacts,
        dt_s,
        controls,
        apply_rigid_frictional_impacts,
        admit_path,
        gate,
    )
}

// Constitutive dispatch changes only the simultaneous impulse owner. All
// force work, native drift, geometric coverage and final publication stay here.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn step_with<F: Copy, R: ImpactResponse, I, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    wrenches: &[Wrench],
    gravity: Gravity,
    contacts: &[RigidFrictionalContact<F>],
    dt_s: f64,
    controls: RigidContactStepControls,
    apply_impact: impl FnOnce(
        &mut [RigidBodyState],
        &[MassProperties],
        &[RigidFrictionalContact<F>],
        RigidImpactControls,
        &CancelGate,
    ) -> Result<
        RigidFrictionalImpactReport<R>,
        RigidFrictionalImpactError<R, I>,
    >,
    admit_path: impl FnOnce(&RigidContactPath<'_, F>) -> Result<P, E>,
    gate: &CancelGate,
) -> Result<RigidContactStepReport<P, R>, RigidContactStepError<E, R, I>> {
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
    let (before_energy, _, _) = super::totals_with::<R, I>(bodies, properties, gate)?;
    let mut kicked = Vec::with_capacity(bodies.len());
    let mut force_work = Dd::ZERO;
    let mut torque_work = Dd::ZERO;
    let mut kick_change = Dd::ZERO;
    let mut linear_residual = [Dd::ZERO; 3];
    let mut angular_residual = [Dd::ZERO; 3];
    for ((&body, &property), &wrench) in bodies.iter().zip(properties).zip(wrenches) {
        poll(gate)?;
        wrench.validate()?;
        let force = wrench
            .force_world
            .add(gravity.acceleration_world().scale(property.mass()));
        // Apply both resultants through the original atomic native owner. Its
        // increment/impulse residual exposes kicks lost by rounded state addition
        // so the original acceptance gates can refuse an unresolved update.
        let force_impulse = force.scale(dt_s);
        let torque_impulse = wrench.torque_body.scale(dt_s);
        let receipt = body.apply_momentum_impulse(property, force_impulse, torque_impulse)?;
        let next = receipt.state_after;
        let velocity_mid = body
            .center_of_mass_velocity_world(property)?
            .scale(0.5)
            .add(next.center_of_mass_velocity_world(property)?.scale(0.5));
        let omega_before = property.angular_velocity_body_checked(body.angular_momentum_body())?;
        let omega_after = property.angular_velocity_body_checked(next.angular_momentum_body())?;
        force_work = force_work + Dd::from_f64(force_impulse.dot(velocity_mid));
        torque_work = torque_work
            + Dd::from_f64(torque_impulse.dot(omega_before.scale(0.5).add(omega_after.scale(0.5))));
        kick_change = kick_change + Dd::from_f64(receipt.work.kinetic_energy_change_j);
        super::add_vector(&mut linear_residual, receipt.linear_residual_world);
        super::add_vector(
            &mut angular_residual,
            body.pose()
                .position_world()
                .cross(receipt.linear_residual_world),
        );
        super::add_vector(
            &mut angular_residual,
            body.pose()
                .orientation()
                .rotate_body_to_world(receipt.angular_residual_body),
        );
        kicked.push(next);
    }
    let kick_residual = (kick_change - force_work - torque_work).to_f64();
    if !kick_residual.is_finite() || kick_residual.abs() > controls.work_energy_tolerance_j {
        return Err(RigidContactStepError::Audit(
            "native force/torque impulse work failed",
        ));
    }
    let impact = apply_impact(&mut kicked, properties, contacts, controls.impact, gate)?;
    super::add_vector(&mut linear_residual, impact.linear_momentum_residual_ns);
    super::add_vector(&mut angular_residual, impact.angular_momentum_residual_nms);
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut endpoint = Vec::with_capacity(kicked.len());
    let mut drift_energy = Dd::ZERO;
    let mut drift_linear = [Dd::ZERO; 3];
    let mut drift_angular = [Dd::ZERO; 3];
    for (&body, &property) in kicked.iter().zip(properties) {
        poll(gate)?;
        let after = integrator
            .step(body, property, Wrench::ZERO, dt_s)?
            .state_after;
        let dp = after
            .linear_momentum_world()
            .sub(body.linear_momentum_world());
        let dl = after
            .angular_momentum_body()
            .sub(body.angular_momentum_body());
        let velocity_mid = body
            .center_of_mass_velocity_world(property)?
            .scale(0.5)
            .add(after.center_of_mass_velocity_world(property)?.scale(0.5));
        let omega_mid = property
            .angular_velocity_body_checked(body.angular_momentum_body())?
            .scale(0.5)
            .add(
                property
                    .angular_velocity_body_checked(after.angular_momentum_body())?
                    .scale(0.5),
            );
        drift_energy =
            drift_energy + Dd::from_f64(dp.dot(velocity_mid)) + Dd::from_f64(dl.dot(omega_mid));
        super::add_vector(&mut drift_linear, dp);
        // Delta(x cross p) = x_before cross delta_p + delta_x cross p_after.
        // Unlike subtracting two orbital totals, this retains small exchange
        // alongside a large unchanged origin/translation contribution.
        let dx = after
            .pose()
            .position_world()
            .sub(body.pose().position_world());
        super::add_vector(&mut drift_angular, body.pose().position_world().cross(dp));
        super::add_vector(&mut drift_angular, dx.cross(after.linear_momentum_world()));
        super::add_vector(
            &mut drift_angular,
            after
                .pose()
                .orientation()
                .rotate_body_to_world(after.angular_momentum_body())
                .sub(
                    body.pose()
                        .orientation()
                        .rotate_body_to_world(body.angular_momentum_body()),
                ),
        );
        endpoint.push(after);
    }
    let endpoint_gaps = gaps(&endpoint, properties, contacts, dt_s, true, gate)?;
    let (after_energy, _, _) = super::totals_with::<R, I>(&endpoint, properties, gate)?;
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
    let total_energy = (kick_change + Dd::from_f64(impact.kinetic_change_j) + drift_energy
        - force_work
        - torque_work
        - Dd::from_f64(impact.support_work_j)
        + Dd::from_f64(impact.contact.physical().dissipated_energy))
    .to_f64();
    for axis in 0..3 {
        linear_residual[axis] = linear_residual[axis] + drift_linear[axis];
        angular_residual[axis] = angular_residual[axis] + drift_angular[axis];
    }
    let linear_residual = super::rounded_vector(linear_residual);
    let angular_residual = super::rounded_vector(angular_residual);
    let drift_linear = super::rounded_vector(drift_linear);
    let drift_angular = super::rounded_vector(drift_angular);
    let drift_energy = drift_energy.to_f64();
    let force_work = force_work.to_f64();
    let torque_work = torque_work.to_f64();
    let force_splitting_work = force_work - force_displacement_work;
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
        super::magnitude(drift_linear),
        super::magnitude(linear_residual),
        super::magnitude(angular_residual),
    ]
    .iter()
    .any(|x| !x.is_finite())
        || drift_energy.abs() > controls.drift_energy_tolerance_j
        || super::magnitude(drift_angular) > controls.drift_angular_momentum_tolerance_nms
        || super::magnitude(drift_linear) > controls.linear_momentum_tolerance_ns
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

fn gaps<F, E, R, I>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact<F>],
    elapsed_s: f64,
    endpoint: bool,
    gate: &CancelGate,
) -> Result<Vec<RigidContactGap>, RigidContactStepError<E, R, I>> {
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
            .point_kinematics(properties[contact.point.body], contact.point.arm_body)?;
        let (centre_b, arm_b) = match contact.other {
            RigidImpactCounterpart::Body(point) => {
                if point.body >= bodies.len() || point.body == contact.point.body {
                    return Err(RigidContactStepError::Invalid(
                        "invalid counterpart witness body",
                    ));
                }
                let b =
                    bodies[point.body].point_kinematics(properties[point.body], point.arm_body)?;
                (bodies[point.body].pose().position_world(), b.arm_world)
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
                (point_world, velocity_world.scale(elapsed_s))
            }
        };
        let gap = normal.dot(super::material_separation(
            bodies[contact.point.body].pose().position_world(),
            a.arm_world,
            centre_b,
            arm_b,
        ));
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
fn poll<E, R, I>(gate: &CancelGate) -> Result<(), RigidContactStepError<E, R, I>> {
    if gate.is_requested() {
        Err(RigidContactStepError::Cancelled)
    } else {
        Ok(())
    }
}
