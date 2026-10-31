//! Geometry-admitted first-order force/contact/drift on the complete sphere scene.
//!
//! Original solid forces and native rigid force/torque impulses act first. ONE
//! original contact solve then couples every active rigid/rigid and rigid/solid
//! row. Original kinetic subflows advance the actual native configurations.
//! Restitution targets use the incoming PRE-FORCE velocities, so a held load
//! does not masquerade as a new impact at every resting-contact step.
//!
//! Every retry restarts the entire kick/contact/drift map from the unchanged
//! accepted state. The shared geometry owner admits all affine sphere/triangle,
//! sphere/sphere and tetrahedron-orientation paths; endpoints are never repaired.
//! COM-centred spheres make their geometric centre path affine despite body
//! rotation. This is the declared numerical path, not an exact-flow certificate.
//!
//! This is a bounded first-order split, not a general Moreau-Jean method or an
//! existence/uniqueness proof for arbitrary rigid contact. Native determinacy
//! refusals survive. Normal/impact loss and splitting defects remain unallocated;
//! neither is silently called heat. Explicit elastic and force splitting must
//! pass noncancelling refinement gates, and finite budgets can refuse progress.

use super::free::{rotation_difference, rotational_change};
use super::{
    FiniteSolidSpherePathRefusal, FiniteSolidSpherePathReport, FiniteSolidSphereScene,
    FiniteSolidSphereSceneKey, FiniteSolidSphereTrajectoryControls,
    FiniteSolidSphereTrajectoryError, admit_path, current_contacts, norm, validate, vector,
};
use crate::finite_solid::contact::mixed::{FiniteSolidContactCounterpart, FiniteSolidMixedContact};
use crate::finite_solid::contact::surface::FiniteSolidContactSurface;
use crate::finite_solid::contact::{
    FiniteSolidRigidContactError, FiniteSolidRigidContactReport, contact_frame,
};
use crate::finite_solid::impulse::{FiniteSolidImpulse, FiniteSolidMomentum};
use crate::finite_solid::split::{FiniteSolidContactDrift, FiniteSolidForceKick};
use crate::finite_solid::{FiniteSolidBody, FiniteSolidError};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::Cx;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, MomentumImpulseReceipt, RigidBodyIntegrator,
    RigidBodyState, StepReceipt, Vec3, Wrench,
};
use fs_solid::hyper3::HyperTetError;

/// Compact kick receipt; private continuation vectors are not retained per step.
#[derive(Clone, Debug)]
pub struct FiniteSolidSustainedSolidKick {
    pub impulse: FiniteSolidImpulse,
    pub source_impulse: FiniteSolidMomentum,
    pub internal_impulse: FiniteSolidMomentum,
    pub source_impulse_work_j: f64,
    pub internal_impulse_work_j: f64,
    pub work_decomposition_residual_j: f64,
    pub linear_decomposition_residual_kg_m_s: [f64; 3],
    pub angular_decomposition_residual_kg_m2_s: [f64; 3],
}
impl From<FiniteSolidForceKick<'_>> for FiniteSolidSustainedSolidKick {
    fn from(kick: FiniteSolidForceKick<'_>) -> Self {
        Self {
            impulse: kick.impulse,
            source_impulse: kick.source_impulse,
            internal_impulse: kick.internal_impulse,
            source_impulse_work_j: kick.source_impulse_work_j,
            internal_impulse_work_j: kick.internal_impulse_work_j,
            work_decomposition_residual_j: kick.work_decomposition_residual_j,
            linear_decomposition_residual_kg_m_s: kick.linear_decomposition_residual_kg_m_s,
            angular_decomposition_residual_kg_m2_s: kick.angular_decomposition_residual_kg_m2_s,
        }
    }
}

/// Actual original full-inertia kick and zero-force rotational drift.
#[derive(Clone, Debug)]
pub struct FiniteSolidSustainedRigidStep {
    pub kick: MomentumImpulseReceipt,
    pub drift: StepReceipt,
    pub force_impulse_work_j: f64,
    /// Body-frame torque impulse dotted with its actual midpoint angular speed.
    /// This is the declared discrete torque work, not an exact forced-top flow.
    pub torque_impulse_work_j: f64,
    pub force_displacement_work_j: f64,
    pub force_splitting_work_j: f64,
    pub kinetic_drift_defect_j: f64,
    pub drift_linear_increment_kg_m_s: Vec3,
    pub drift_angular_increment_kg_m2_s: Vec3,
    pub kick_angular_residual_kg_m2_s: Vec3,
    /// Sum of separate kick orbital/spin residual norms and free-drift
    /// orbital/spin increment norms. No exchange exists between those channels
    /// in these phases, so opposing errors cannot cancel the admission gate.
    pub absolute_angular_residual_kg_m2_s: f64,
}

#[derive(Clone, Debug)]
pub struct FiniteSolidSustainedTrajectoryStep {
    pub elapsed_s: f64,
    pub duration_s: f64,
    pub maximum_contact_gap_m: f64,
    pub contact_keys: Vec<FiniteSolidSphereSceneKey>,
    /// Same scene-row order as keys, before this interval's force kick.
    pub incoming_normal_speeds_m_s: Vec<f64>,
    pub normal_targets_m_s: Vec<f64>,
    pub solid_kick: FiniteSolidSustainedSolidKick,
    pub solid_drift: FiniteSolidContactDrift,
    pub rigid: Vec<FiniteSolidSustainedRigidStep>,
    pub impact: Option<FiniteSolidRigidContactReport>,
    pub path: FiniteSolidSpherePathReport,
}

/// Whole-call storage, source/support transfer and separately visible split errors.
#[derive(Clone, Debug)]
pub struct FiniteSolidSustainedTrajectoryReport {
    pub duration_s: f64,
    pub rejected_steps: usize,
    pub steps: Vec<FiniteSolidSustainedTrajectoryStep>,
    pub kinetic_change_j: f64,
    pub strain_change_j: f64,
    /// Solid/rigid force displacement work and the declared body-torque kick work.
    pub source_work_j: f64,
    /// Native impact-map loss, including its normal and first-order splitting
    /// contribution. This is not an authoritative friction-heat partition.
    pub unallocated_contact_loss_j: f64,
    pub solid_source_splitting_work_j: f64,
    pub solid_elastic_splitting_work_j: f64,
    pub rigid_force_splitting_work_j: f64,
    pub rigid_drift_energy_defect_j: f64,
    pub integration_defect_j: f64,
    pub absolute_integration_defect_j: f64,
    pub work_energy_residual_j: f64,
    pub linear_momentum_change_kg_m_s: Vec3,
    pub angular_momentum_change_kg_m2_s: Vec3,
    pub external_impulse_kg_m_s: Vec3,
    pub external_angular_impulse_kg_m2_s: Vec3,
    pub support_impulse_kg_m_s: Vec3,
    pub support_angular_impulse_kg_m2_s: Vec3,
    pub linear_momentum_residual_kg_m_s: Vec3,
    pub angular_momentum_residual_kg_m2_s: Vec3,
    pub accumulated_linear_residual_kg_m_s: f64,
    pub accumulated_angular_residual_kg_m2_s: f64,
    pub accumulated_energy_residual_j: f64,
}

#[derive(Clone, Debug)]
pub enum FiniteSolidSustainedRefusal {
    Path(FiniteSolidSpherePathRefusal),
    Integration {
        defect_j: f64,
        allowance_j: f64,
    },
    /// Original elastic owner refused a trial endpoint; no accepted state changed.
    SolidTrial(HyperTetError),
    /// Original bounded torque-free solve requests a shorter interval.
    RotationalSolve {
        rigid_body: usize,
    },
}

/// Every refusal leaves accepted q/v/a, native clock/history and rigid states intact.
#[derive(Debug)]
pub enum FiniteSolidSustainedTrajectoryError {
    Invalid(&'static str),
    Trajectory(FiniteSolidSphereTrajectoryError),
    Solid(FiniteSolidError),
    Rigid(DynamicsError),
    Impact(FiniteSolidRigidContactError),
    RejectionBudget {
        elapsed_s: f64,
        last: FiniteSolidSustainedRefusal,
    },
    PhaseAudit {
        rigid_body: usize,
        linear_residual_kg_m_s: f64,
        angular_residual_kg_m2_s: f64,
        energy_residual_j: f64,
    },
    Audit(Box<FiniteSolidSustainedTrajectoryReport>),
}
impl From<FiniteSolidSphereTrajectoryError> for FiniteSolidSustainedTrajectoryError {
    fn from(e: FiniteSolidSphereTrajectoryError) -> Self {
        Self::Trajectory(e)
    }
}
impl From<FiniteSolidError> for FiniteSolidSustainedTrajectoryError {
    fn from(e: FiniteSolidError) -> Self {
        Self::Solid(e)
    }
}
impl From<DynamicsError> for FiniteSolidSustainedTrajectoryError {
    fn from(e: DynamicsError) -> Self {
        Self::Rigid(e)
    }
}
impl From<FiniteSolidRigidContactError> for FiniteSolidSustainedTrajectoryError {
    fn from(e: FiniteSolidRigidContactError) -> Self {
        Self::Impact(e)
    }
}
impl core::fmt::Display for FiniteSolidSustainedTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native sustained rigid/finite-solid contact: {self:?}")
    }
}
impl std::error::Error for FiniteSolidSustainedTrajectoryError {}

impl<'a> FiniteSolidBody<'a> {
    /// Advance the complete declared sphere scene through loaded sustained contact.
    ///
    /// Forces use the same conventions as `rigid_free_flight_candidate`: rigid
    /// gravity is separate; any solid gravity belongs in explicit solid loads.
    /// Solid history records actual kinetic-subflow solves, with original smooth
    /// endpoint acceleration. Newton targets use pre-force contact velocity;
    /// explicit target-speed laws remain explicit. Geometry, source coefficients
    /// and complete-set contact determinacy retain their original refusals.
    ///
    /// The integration tolerance bounds a sum of ABSOLUTE phase defects. Each
    /// interval receives its duration fraction of that budget; failing force or
    /// elastic splitting halves the entire interval and recomputes every phase.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn sustained_sphere_trajectory_candidate(
        &self,
        surface: &FiniteSolidContactSurface,
        scene: &FiniteSolidSphereScene,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        wrenches: &[Wrench],
        gravity: Gravity,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FiniteSolidSphereTrajectoryControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            Vec<RigidBodyState>,
            FiniteSolidSustainedTrajectoryReport,
        ),
        FiniteSolidSustainedTrajectoryError,
    > {
        super::poll(cx)?;
        validate(
            bodies,
            properties,
            scene.spheres(),
            wrenches,
            duration_s,
            controls,
        )?;
        if !vector(origin_m).is_finite()
            || [
                controls.free.linear_tolerance_kg_m_s,
                controls.free.angular_tolerance_kg_m2_s,
                controls.free.energy_tolerance_j,
                controls.free.solid.energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(FiniteSolidSustainedTrajectoryError::Invalid(
                "finite origin and positive phase gates required",
            ));
        }
        let start_time = self.native_state().t;
        let end_time = start_time + duration_s;
        if !end_time.is_finite() || end_time <= start_time {
            return Err(FiniteSolidSphereTrajectoryError::NoProgress { elapsed_s: 0.0 }.into());
        }
        current_contacts(
            surface,
            Some(scene),
            self,
            bodies,
            scene.spheres(),
            controls.surface,
            cx,
        )?;
        let mut solid = self.clone();
        let mut rigid = bodies.to_vec();
        let mut report = empty_report(duration_s);
        let mut elapsed = 0.0;
        let mut proposed_step_s = controls.maximum_step_s;
        let clock_tolerance = 8.0 * f64::EPSILON * duration_s;
        while elapsed < duration_s {
            super::poll(cx)?;
            if report.steps.len() == controls.maximum_steps {
                return Err(
                    FiniteSolidSphereTrajectoryError::StepBudget { elapsed_s: elapsed }.into(),
                );
            }
            let (contacts, keys, maximum_gap) = current_contacts(
                surface,
                Some(scene),
                &solid,
                &rigid,
                scene.spheres(),
                controls.surface,
                cx,
            )?;
            let (contacts, incoming, targets) =
                incoming_laws(&solid, &rigid, properties, &contacts, cx)?;
            let mut dt = proposed_step_s.min(duration_s - elapsed);
            let (next_solid, next_rigid, mut step) = loop {
                super::poll(cx)?;
                if !dt.is_finite()
                    || dt <= 0.0
                    || elapsed + dt <= elapsed
                    || solid.native_state().t + dt <= solid.native_state().t
                {
                    return Err(FiniteSolidSphereTrajectoryError::NoProgress {
                        elapsed_s: elapsed,
                    }
                    .into());
                }
                let allowance = controls.integration_tolerance_j * (dt / duration_s);
                if !allowance.is_finite() || allowance <= 0.0 {
                    return Err(FiniteSolidSphereTrajectoryError::NoProgress {
                        elapsed_s: elapsed,
                    }
                    .into());
                }
                let attempt = interval(
                    &solid,
                    surface,
                    scene,
                    &rigid,
                    properties,
                    &contacts,
                    nodal_forces_n,
                    body_force_n_m3,
                    wrenches,
                    gravity,
                    dt,
                    origin_m,
                    controls,
                    allowance,
                    clock_tolerance,
                    cx,
                )?;
                match attempt {
                    Ok(result) => break result,
                    Err(last) => {
                        if report.rejected_steps == controls.maximum_rejections {
                            return Err(FiniteSolidSustainedTrajectoryError::RejectionBudget {
                                elapsed_s: elapsed,
                                last,
                            });
                        }
                        report.rejected_steps += 1;
                        let fraction = match last {
                            FiniteSolidSustainedRefusal::Path(p) => p.start(),
                            FiniteSolidSustainedRefusal::Integration { .. }
                            | FiniteSolidSustainedRefusal::SolidTrial(_)
                            | FiniteSolidSustainedRefusal::RotationalSolve { .. } => 0.0,
                        };
                        dt *= if fraction > 0.0 { 0.9 * fraction } else { 0.5 };
                    }
                }
            };
            let next_elapsed = elapsed + dt;
            if next_elapsed <= elapsed || next_elapsed > duration_s {
                return Err(
                    FiniteSolidSphereTrajectoryError::NoProgress { elapsed_s: elapsed }.into(),
                );
            }
            elapsed = next_elapsed;
            if (next_solid.native_state().t - (start_time + elapsed)).abs() > clock_tolerance {
                return Err(FiniteSolidSphereTrajectoryError::TimeResolution {
                    requested_step_s: elapsed,
                    represented_step_s: next_solid.native_state().t - start_time,
                }
                .into());
            }
            step.elapsed_s = elapsed;
            step.maximum_contact_gap_m = maximum_gap;
            step.contact_keys = keys;
            step.incoming_normal_speeds_m_s = incoming;
            step.normal_targets_m_s = targets;
            // Retain an admitted scale rather than repeatedly paying the same
            // force/contact solves to rediscover it. Grow only with a factor-four
            // work margin; geometry admission still checks every new proposal.
            let split_defect = step.solid_drift.absolute_integration_defect_j
                + step
                    .rigid
                    .iter()
                    .map(|body| {
                        body.force_splitting_work_j.abs() + body.kinetic_drift_defect_j.abs()
                    })
                    .sum::<f64>();
            proposed_step_s = if split_defect
                <= 0.25 * controls.integration_tolerance_j * (dt / duration_s)
                && step.solid_drift.absolute_integration_defect_j
                    <= 0.25 * controls.free.solid.energy_tolerance_j
            {
                (2.0 * dt).min(controls.maximum_step_s)
            } else {
                dt
            };
            accumulate(&mut report, &step, vector(origin_m));
            report.steps.push(step);
            solid = next_solid;
            rigid = next_rigid;
        }
        if (solid.native_state().t - end_time).abs() > clock_tolerance {
            return Err(FiniteSolidSphereTrajectoryError::TimeResolution {
                requested_step_s: duration_s,
                represented_step_s: solid.native_state().t - start_time,
            }
            .into());
        }
        finish(&mut report);
        if !admitted(&report, controls) {
            return Err(FiniteSolidSustainedTrajectoryError::Audit(Box::new(report)));
        }
        super::poll(cx)?;
        Ok((solid, rigid, report))
    }
}

fn incoming_laws(
    solid: &FiniteSolidBody<'_>,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[FiniteSolidMixedContact],
    cx: &Cx<'_>,
) -> Result<(Vec<FiniteSolidMixedContact>, Vec<f64>, Vec<f64>), FiniteSolidSustainedTrajectoryError>
{
    let mut result = contacts.to_vec();
    let mut incoming = Vec::with_capacity(contacts.len());
    let mut targets = Vec::with_capacity(contacts.len());
    for contact in &mut result {
        super::poll(cx)?;
        let axes = contact_frame(contact.normal_world, contact.tangent_world)
            .map_err(FiniteSolidSustainedTrajectoryError::Invalid)?;
        let velocity = bodies[contact.point.body]
            .point_kinematics(properties[contact.point.body], contact.point.arm_body)?
            .point_velocity_world;
        let other = match contact.counterpart {
            FiniteSolidContactCounterpart::Solid(p) => vector(solid.point_velocity_m_s(p, cx)?),
            FiniteSolidContactCounterpart::Rigid(p) => {
                bodies[p.body]
                    .point_kinematics(properties[p.body], p.arm_body)?
                    .point_velocity_world
            }
        };
        let speed = axes[0].dot(velocity.sub(other));
        let target = match contact.normal_law {
            NormalVelocityLaw::NewtonRestitution { coefficient }
                if coefficient.is_finite() && (0.0..=1.0).contains(&coefficient) =>
            {
                -coefficient * speed.min(0.0)
            }
            NormalVelocityLaw::Target { speed } if speed.is_finite() && speed >= 0.0 => speed,
            _ => {
                return Err(FiniteSolidSustainedTrajectoryError::Invalid(
                    "valid incoming normal law required",
                ));
            }
        };
        if !speed.is_finite() || !target.is_finite() {
            return Err(FiniteSolidSustainedTrajectoryError::Invalid(
                "unrepresentable incoming contact speed",
            ));
        }
        incoming.push(speed);
        targets.push(target);
        contact.normal_law = NormalVelocityLaw::Target { speed: target };
    }
    Ok((result, incoming, targets))
}

type Attempt<'a> = Result<
    (
        FiniteSolidBody<'a>,
        Vec<RigidBodyState>,
        FiniteSolidSustainedTrajectoryStep,
    ),
    FiniteSolidSustainedRefusal,
>;

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn interval<'a>(
    solid: &FiniteSolidBody<'a>,
    surface: &FiniteSolidContactSurface,
    scene: &FiniteSolidSphereScene,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[FiniteSolidMixedContact],
    nodal_forces_n: &[(usize, f64)],
    body_force_n_m3: [f64; 3],
    wrenches: &[Wrench],
    gravity: Gravity,
    dt: f64,
    origin_m: [f64; 3],
    controls: FiniteSolidSphereTrajectoryControls,
    allowance: f64,
    clock_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<Attempt<'a>, FiniteSolidSustainedTrajectoryError> {
    let (kicked_solid, solid_kick) = solid.force_kick_candidate(
        dt,
        nodal_forces_n,
        body_force_n_m3,
        origin_m,
        controls.impact.solid_impulse,
        cx,
    )?;
    let mut kicks = Vec::with_capacity(bodies.len());
    let mut kicked_rigid = Vec::with_capacity(bodies.len());
    for ((&body, &property), &wrench) in bodies.iter().zip(properties).zip(wrenches) {
        super::poll(cx)?;
        let force = wrench
            .force_world
            .add(gravity.acceleration_world().scale(property.mass()));
        let kick =
            body.apply_momentum_impulse(property, force.scale(dt), wrench.torque_body.scale(dt))?;
        kicked_rigid.push(kick.state_after);
        kicks.push(kick);
    }
    let (contacted_solid, contacted_rigid, impact) = if contacts.is_empty() {
        (kicked_solid, kicked_rigid, None)
    } else {
        let (s, r, i) = kicked_solid.mixed_contact_candidate(
            &kicked_rigid,
            properties,
            contacts,
            origin_m,
            controls.impact,
            cx,
        )?;
        (s, r, Some(i))
    };
    let mut solid_controls = controls.free.solid;
    solid_controls.energy_tolerance_j = solid_controls.energy_tolerance_j.min(allowance);
    let (next_solid, solid_drift) =
        match contacted_solid.contact_drift_candidate(&solid_kick, solid_controls, cx) {
            Ok(value) => value,
            Err(FiniteSolidError::Energy {
                defect_j,
                allowance_j,
            }) => {
                return Ok(Err(FiniteSolidSustainedRefusal::Integration {
                    defect_j,
                    allowance_j,
                }));
            }
            Err(FiniteSolidError::Solid(error @ HyperTetError::MaterialRefused { .. })) => {
                return Ok(Err(FiniteSolidSustainedRefusal::SolidTrial(error)));
            }
            Err(error) => return Err(error.into()),
        };
    let represented = next_solid.native_state().t - solid.native_state().t;
    if (represented - dt).abs() > clock_tolerance {
        return Err(FiniteSolidSphereTrajectoryError::TimeResolution {
            requested_step_s: dt,
            represented_step_s: represented,
        }
        .into());
    }
    let origin = vector(origin_m);
    let mut next_rigid = Vec::with_capacity(bodies.len());
    let mut rigid = Vec::with_capacity(bodies.len());
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut absolute_defect = solid_drift.absolute_integration_defect_j;
    for (index, (((&contacted, &property), &wrench), kick)) in contacted_rigid
        .iter()
        .zip(properties)
        .zip(wrenches)
        .zip(kicks)
        .enumerate()
    {
        super::poll(cx)?;
        let drift = match integrator.step(contacted, property, Wrench::ZERO, dt) {
            Ok(drift) => drift,
            Err(DynamicsError::RotationalSolveFailed) => {
                return Ok(Err(FiniteSolidSustainedRefusal::RotationalSolve {
                    rigid_body: index,
                }));
            }
            Err(error) => return Err(error.into()),
        };
        let next = drift.state_after;
        let dp = next
            .linear_momentum_world()
            .sub(contacted.linear_momentum_world());
        let p_average = contacted
            .linear_momentum_world()
            .scale(0.5)
            .add(next.linear_momentum_world().scale(0.5));
        let dx = next
            .pose()
            .position_world()
            .sub(contacted.pose().position_world());
        let r_average = contacted
            .pose()
            .position_world()
            .sub(origin)
            .scale(0.5)
            .add(next.pose().position_world().sub(origin).scale(0.5));
        let dl = next
            .angular_momentum_body()
            .sub(contacted.angular_momentum_body());
        let drift_orbital = r_average.cross(dp).add(dx.cross(p_average));
        let drift_spin =
            contacted
                .pose()
                .orientation()
                .rotate_body_to_world(dl)
                .add(rotation_difference(
                    contacted.pose().orientation(),
                    next.pose().orientation(),
                    next.angular_momentum_body(),
                ));
        let drift_angular = drift_orbital.add(drift_spin);
        let kinetic_drift_defect = p_average.scale(property.mass().recip()).dot(dp)
            + rotational_change(
                property,
                contacted.angular_momentum_body(),
                next.angular_momentum_body(),
            )?;
        let force = wrench
            .force_world
            .add(gravity.acceleration_world().scale(property.mass()));
        let force_displacement_work = force.dot(
            next.pose()
                .position_world()
                .sub(kick.state_before.pose().position_world()),
        );
        let velocity_average = kick
            .state_before
            .center_of_mass_velocity_world(property)?
            .scale(0.5)
            .add(
                kick.state_after
                    .center_of_mass_velocity_world(property)?
                    .scale(0.5),
            );
        let omega_average = property
            .angular_velocity_body_checked(kick.state_before.angular_momentum_body())?
            .scale(0.5)
            .add(
                property
                    .angular_velocity_body_checked(kick.state_after.angular_momentum_body())?
                    .scale(0.5),
            );
        let force_kick_work = kick.linear_impulse_world.dot(velocity_average);
        let torque_kick_work = kick.angular_impulse_body.dot(omega_average);
        let force_splitting = force_kick_work - force_displacement_work;
        let kick_orbital_residual = kick
            .state_before
            .pose()
            .position_world()
            .sub(origin)
            .cross(kick.linear_residual_world);
        let kick_spin_residual = kick
            .state_before
            .pose()
            .orientation()
            .rotate_body_to_world(kick.angular_residual_body);
        let kick_angular_residual = kick_orbital_residual.add(kick_spin_residual);
        let linear_residual = norm(kick.linear_residual_world) + norm(dp);
        let angular_residual = norm(kick_orbital_residual)
            + norm(kick_spin_residual)
            + norm(drift_orbital)
            + norm(drift_spin);
        let energy_residual = kick.work.work_energy_residual_j.abs()
            + (kick.work.impulse_work_j - (force_kick_work + torque_kick_work)).abs();
        if [
            linear_residual,
            angular_residual,
            energy_residual,
            kinetic_drift_defect,
            force_displacement_work,
            force_splitting,
            torque_kick_work,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || linear_residual > controls.free.linear_tolerance_kg_m_s
            || angular_residual > controls.free.angular_tolerance_kg_m2_s
            || energy_residual > controls.free.energy_tolerance_j
        {
            return Err(FiniteSolidSustainedTrajectoryError::PhaseAudit {
                rigid_body: index,
                linear_residual_kg_m_s: linear_residual,
                angular_residual_kg_m2_s: angular_residual,
                energy_residual_j: energy_residual,
            });
        }
        absolute_defect += force_splitting.abs() + kinetic_drift_defect.abs();
        next_rigid.push(next);
        rigid.push(FiniteSolidSustainedRigidStep {
            kick,
            drift,
            force_impulse_work_j: force_kick_work,
            torque_impulse_work_j: torque_kick_work,
            force_displacement_work_j: force_displacement_work,
            force_splitting_work_j: force_splitting,
            kinetic_drift_defect_j: kinetic_drift_defect,
            drift_linear_increment_kg_m_s: dp,
            drift_angular_increment_kg_m2_s: drift_angular,
            kick_angular_residual_kg_m2_s: kick_angular_residual,
            absolute_angular_residual_kg_m2_s: angular_residual,
        });
    }
    if !absolute_defect.is_finite() {
        return Err(FiniteSolidSustainedTrajectoryError::Invalid(
            "unrepresentable absolute split defect",
        ));
    }
    if absolute_defect > allowance {
        return Ok(Err(FiniteSolidSustainedRefusal::Integration {
            defect_j: absolute_defect,
            allowance_j: allowance,
        }));
    }
    let path = match admit_path(
        surface,
        &contacted_solid,
        &next_solid,
        &contacted_rigid,
        &next_rigid,
        scene.spheres(),
        true,
        controls,
        cx,
    )? {
        Ok(path) => path,
        Err(refusal) => return Ok(Err(FiniteSolidSustainedRefusal::Path(refusal))),
    };
    current_contacts(
        surface,
        Some(scene),
        &next_solid,
        &next_rigid,
        scene.spheres(),
        controls.surface,
        cx,
    )?;
    Ok(Ok((
        next_solid,
        next_rigid,
        FiniteSolidSustainedTrajectoryStep {
            elapsed_s: 0.0,
            duration_s: dt,
            maximum_contact_gap_m: 0.0,
            contact_keys: Vec::new(),
            incoming_normal_speeds_m_s: Vec::new(),
            normal_targets_m_s: Vec::new(),
            solid_kick: solid_kick.into(),
            solid_drift,
            rigid,
            impact,
            path,
        },
    )))
}

fn empty_report(duration_s: f64) -> FiniteSolidSustainedTrajectoryReport {
    FiniteSolidSustainedTrajectoryReport {
        duration_s,
        rejected_steps: 0,
        steps: Vec::new(),
        kinetic_change_j: 0.0,
        strain_change_j: 0.0,
        source_work_j: 0.0,
        unallocated_contact_loss_j: 0.0,
        solid_source_splitting_work_j: 0.0,
        solid_elastic_splitting_work_j: 0.0,
        rigid_force_splitting_work_j: 0.0,
        rigid_drift_energy_defect_j: 0.0,
        integration_defect_j: 0.0,
        absolute_integration_defect_j: 0.0,
        work_energy_residual_j: 0.0,
        linear_momentum_change_kg_m_s: Vec3::ZERO,
        angular_momentum_change_kg_m2_s: Vec3::ZERO,
        external_impulse_kg_m_s: Vec3::ZERO,
        external_angular_impulse_kg_m2_s: Vec3::ZERO,
        support_impulse_kg_m_s: Vec3::ZERO,
        support_angular_impulse_kg_m2_s: Vec3::ZERO,
        linear_momentum_residual_kg_m_s: Vec3::ZERO,
        angular_momentum_residual_kg_m2_s: Vec3::ZERO,
        accumulated_linear_residual_kg_m_s: 0.0,
        accumulated_angular_residual_kg_m2_s: 0.0,
        accumulated_energy_residual_j: 0.0,
    }
}

fn accumulate(
    r: &mut FiniteSolidSustainedTrajectoryReport,
    s: &FiniteSolidSustainedTrajectoryStep,
    origin: Vec3,
) {
    let k = &s.solid_kick;
    let d = &s.solid_drift;
    let native = &k.impulse;
    r.kinetic_change_j += native.kinetic_change_j + d.kinetic_change_j;
    r.strain_change_j += d.strain_after_j - d.strain_before_j;
    r.source_work_j += d.source_displacement_work_j;
    r.solid_source_splitting_work_j += d.source_splitting_work_j;
    r.solid_elastic_splitting_work_j += d.elastic_splitting_work_j;
    r.integration_defect_j += d.integration_defect_j;
    r.absolute_integration_defect_j += d.absolute_integration_defect_j;
    // Native fixed-configuration residual is actual change minus applied minus
    // support. Retain measured internal-force moment, never count it as source.
    let delta_p = vector(native.applied.linear_kg_m_s)
        .add(vector(native.support.linear_kg_m_s))
        .add(vector(native.linear_residual_kg_m_s));
    let delta_l = vector(native.applied.angular_kg_m2_s)
        .add(vector(native.support.angular_kg_m2_s))
        .add(vector(native.angular_residual_kg_m2_s));
    r.linear_momentum_change_kg_m_s = r
        .linear_momentum_change_kg_m_s
        .add(delta_p)
        .add(vector(d.linear_residual_kg_m_s));
    r.angular_momentum_change_kg_m2_s = r
        .angular_momentum_change_kg_m2_s
        .add(delta_l)
        .add(vector(d.angular_residual_kg_m2_s));
    r.external_impulse_kg_m_s = r
        .external_impulse_kg_m_s
        .add(vector(k.source_impulse.linear_kg_m_s));
    r.external_angular_impulse_kg_m2_s = r
        .external_angular_impulse_kg_m2_s
        .add(vector(k.source_impulse.angular_kg_m2_s));
    r.support_impulse_kg_m_s = r
        .support_impulse_kg_m_s
        .add(vector(native.support.linear_kg_m_s));
    r.support_angular_impulse_kg_m2_s = r
        .support_angular_impulse_kg_m2_s
        .add(vector(native.support.angular_kg_m2_s));
    r.accumulated_linear_residual_kg_m_s += norm(vector(native.linear_residual_kg_m_s))
        + norm(vector(k.linear_decomposition_residual_kg_m_s))
        + norm(vector(k.internal_impulse.linear_kg_m_s))
        + norm(vector(d.linear_residual_kg_m_s));
    r.accumulated_angular_residual_kg_m2_s += norm(vector(native.angular_residual_kg_m2_s))
        + norm(vector(k.angular_decomposition_residual_kg_m2_s))
        + norm(vector(k.internal_impulse.angular_kg_m2_s))
        + norm(vector(d.angular_residual_kg_m2_s));
    r.accumulated_energy_residual_j +=
        native.work_residual_j.abs() + k.work_decomposition_residual_j.abs();
    if let Some(i) = &s.impact {
        r.kinetic_change_j += i.kinetic_change_j;
        r.unallocated_contact_loss_j += i.contact.dissipated_energy;
        let support = vector(i.solid_impulse.support.linear_kg_m_s);
        let moment = vector(i.solid_impulse.support.angular_kg_m2_s);
        r.linear_momentum_change_kg_m_s = r
            .linear_momentum_change_kg_m_s
            .add(support)
            .add(i.linear_momentum_residual_kg_m_s);
        r.angular_momentum_change_kg_m2_s = r
            .angular_momentum_change_kg_m2_s
            .add(moment)
            .add(i.angular_momentum_residual_kg_m2_s);
        r.support_impulse_kg_m_s = r.support_impulse_kg_m_s.add(support);
        r.support_angular_impulse_kg_m2_s = r.support_angular_impulse_kg_m2_s.add(moment);
        let (linear, angular, energy) = super::impact_residuals(i);
        r.accumulated_linear_residual_kg_m_s += linear;
        r.accumulated_angular_residual_kg_m2_s += angular;
        r.accumulated_energy_residual_j += energy;
    }
    for b in &s.rigid {
        let old = b.kick.state_before;
        let impulse_moment = old
            .pose()
            .position_world()
            .sub(origin)
            .cross(b.kick.linear_impulse_world)
            .add(
                old.pose()
                    .orientation()
                    .rotate_body_to_world(b.kick.angular_impulse_body),
            );
        let delta_moment = old
            .pose()
            .position_world()
            .sub(origin)
            .cross(b.kick.linear_increment_world)
            .add(
                old.pose()
                    .orientation()
                    .rotate_body_to_world(b.kick.angular_increment_body),
            );
        r.kinetic_change_j += b.kick.work.kinetic_energy_change_j + b.kinetic_drift_defect_j;
        r.source_work_j += b.force_displacement_work_j + b.torque_impulse_work_j;
        r.rigid_force_splitting_work_j += b.force_splitting_work_j;
        r.rigid_drift_energy_defect_j += b.kinetic_drift_defect_j;
        r.integration_defect_j += b.force_splitting_work_j + b.kinetic_drift_defect_j;
        r.absolute_integration_defect_j +=
            b.force_splitting_work_j.abs() + b.kinetic_drift_defect_j.abs();
        r.linear_momentum_change_kg_m_s = r
            .linear_momentum_change_kg_m_s
            .add(b.kick.linear_increment_world)
            .add(b.drift_linear_increment_kg_m_s);
        r.angular_momentum_change_kg_m2_s = r
            .angular_momentum_change_kg_m2_s
            .add(delta_moment)
            .add(b.drift_angular_increment_kg_m2_s);
        r.external_impulse_kg_m_s = r.external_impulse_kg_m_s.add(b.kick.linear_impulse_world);
        r.external_angular_impulse_kg_m2_s = r.external_angular_impulse_kg_m2_s.add(impulse_moment);
        r.accumulated_linear_residual_kg_m_s +=
            norm(b.kick.linear_residual_world) + norm(b.drift_linear_increment_kg_m_s);
        r.accumulated_angular_residual_kg_m2_s += b.absolute_angular_residual_kg_m2_s;
        r.accumulated_energy_residual_j += b.kick.work.work_energy_residual_j.abs()
            + (b.kick.work.impulse_work_j - (b.force_impulse_work_j + b.torque_impulse_work_j))
                .abs();
    }
}

fn finish(r: &mut FiniteSolidSustainedTrajectoryReport) {
    r.work_energy_residual_j = r.kinetic_change_j + r.strain_change_j - r.source_work_j
        + r.unallocated_contact_loss_j
        - r.integration_defect_j;
    r.linear_momentum_residual_kg_m_s = r
        .linear_momentum_change_kg_m_s
        .sub(r.external_impulse_kg_m_s)
        .sub(r.support_impulse_kg_m_s);
    r.angular_momentum_residual_kg_m2_s = r
        .angular_momentum_change_kg_m2_s
        .sub(r.external_angular_impulse_kg_m2_s)
        .sub(r.support_angular_impulse_kg_m2_s);
}

fn admitted(
    r: &FiniteSolidSustainedTrajectoryReport,
    c: FiniteSolidSphereTrajectoryControls,
) -> bool {
    [
        r.kinetic_change_j,
        r.strain_change_j,
        r.source_work_j,
        r.unallocated_contact_loss_j,
        r.solid_source_splitting_work_j,
        r.solid_elastic_splitting_work_j,
        r.rigid_force_splitting_work_j,
        r.rigid_drift_energy_defect_j,
        r.integration_defect_j,
        r.absolute_integration_defect_j,
        r.work_energy_residual_j,
        r.accumulated_linear_residual_kg_m_s,
        r.accumulated_angular_residual_kg_m2_s,
        r.accumulated_energy_residual_j,
    ]
    .iter()
    .all(|x| x.is_finite())
        && [
            r.linear_momentum_change_kg_m_s,
            r.angular_momentum_change_kg_m2_s,
            r.external_impulse_kg_m_s,
            r.external_angular_impulse_kg_m2_s,
            r.support_impulse_kg_m_s,
            r.support_angular_impulse_kg_m2_s,
            r.linear_momentum_residual_kg_m_s,
            r.angular_momentum_residual_kg_m2_s,
        ]
        .iter()
        .all(|v| v.is_finite())
        && r.accumulated_linear_residual_kg_m_s <= c.linear_tolerance_kg_m_s
        && r.accumulated_angular_residual_kg_m2_s <= c.angular_tolerance_kg_m2_s
        && r.accumulated_energy_residual_j <= c.energy_tolerance_j
        && norm(r.linear_momentum_residual_kg_m_s) <= c.linear_tolerance_kg_m_s
        && norm(r.angular_momentum_residual_kg_m2_s) <= c.angular_tolerance_kg_m2_s
        && r.work_energy_residual_j.abs() <= c.energy_tolerance_j
        && r.absolute_integration_defect_j <= c.integration_tolerance_j
}
