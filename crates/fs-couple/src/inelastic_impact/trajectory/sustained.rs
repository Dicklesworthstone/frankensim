//! Native force/contact/drift with retained Maxwell memory and a complete sphere scene.
//!
//! Restitution uses PRE-FORCE velocities. Every retry repeats the complete native
//! split from the same accepted state; rejected attempts advance no memory/time.
//! Original interval sweeps admit the declared affine coordinate path, including
//! rigid pairs, not the unknown continuum flow. Numerical splitting losses are
//! explicit and require refinement; none becomes an inferred heat deposition.

use super::{norm, vector};
use crate::finite_solid::contact::{
    contact_frame,
    step::{
        FiniteSolidSpherePathControls, FiniteSolidSpherePathRefusal, FiniteSolidSpherePathReport,
        FiniteSolidSphereTrajectoryError, admit_path_at_positions,
        sustained::{
            FiniteSolidSustainedRigidStep,
            rigid::{rigid_contact_drift, rigid_force_kick},
        },
    },
};
use crate::inelastic_impact::{
    multiple::{
        FiniteSolidContactCounterpart, InelasticMixedContact, InelasticMixedContactControls,
        InelasticMixedContactError, InelasticMixedContactReport, visco_mixed_contact_candidate,
    },
    surface::{
        FiniteSolidContactSurface, FiniteSolidSurfaceControls, FiniteSolidSurfaceError,
        scene::{
            FiniteSolidSphereScene, FiniteSolidSphereSceneError, FiniteSolidSphereSceneKey,
            InelasticSphereSceneQuery,
        },
    },
};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_mbd::{DynamicsError, Gravity, MassProperties, RigidBodyState, Vec3, Wrench};
use fs_solid::hyper3::{
    HyperTetError,
    visco::{
        FiniteViscoTetError,
        dynamics::{
            FiniteViscoContactDrift, FiniteViscoContactDriftControls, FiniteViscoForceKick,
            FiniteViscoTetMotion, SolidImpulseReport, SolidMomentum,
        },
    },
};

/// Original physical phase gates and bounded whole-call refinement allowances.
#[derive(Clone, Copy, Debug)]
pub struct ViscoSustainedSphereTrajectoryControls {
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    pub path: FiniteSolidSpherePathControls,
    pub surface: FiniteSolidSurfaceControls,
    /// Also supplies independent per-phase rigid and native force-kick gates.
    pub impact: InelasticMixedContactControls,
    pub solid: FiniteViscoContactDriftControls,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
    /// Noncancelling whole-call budget; each trial receives dt/duration of it.
    pub integration_tolerance_j: f64,
}

/// Compact native force account; private continuation bindings are not retained per step.
#[derive(Clone, Debug)]
pub struct ViscoSustainedSolidKick {
    pub impulse: SolidImpulseReport,
    pub source_impulse: SolidMomentum,
    pub internal_impulse: SolidMomentum,
    pub source_impulse_work_j: f64,
    pub internal_impulse_work_j: f64,
    pub work_decomposition_residual_j: f64,
    pub linear_decomposition_residual_kg_m_s: [f64; 3],
    pub angular_decomposition_residual_kg_m2_s: [f64; 3],
}
impl From<FiniteViscoForceKick<'_>> for ViscoSustainedSolidKick {
    fn from(k: FiniteViscoForceKick<'_>) -> Self {
        Self {
            impulse: k.impulse,
            source_impulse: k.source_impulse,
            internal_impulse: k.internal_impulse,
            source_impulse_work_j: k.source_impulse_work_j,
            internal_impulse_work_j: k.internal_impulse_work_j,
            work_decomposition_residual_j: k.work_decomposition_residual_j,
            linear_decomposition_residual_kg_m_s: k.linear_decomposition_residual_kg_m_s,
            angular_decomposition_residual_kg_m2_s: k.angular_decomposition_residual_kg_m2_s,
        }
    }
}

/// One accepted native split; durations and original scene keys permit replay.
#[derive(Clone, Debug)]
pub struct SustainedSphereTrajectoryStep<P> {
    pub elapsed_s: f64,
    pub duration_s: f64,
    pub maximum_contact_gap_m: f64,
    pub contact_keys: Vec<FiniteSolidSphereSceneKey>,
    pub incoming_normal_speeds_m_s: Vec<f64>,
    pub normal_targets_m_s: Vec<f64>,
    pub solid_kick: ViscoSustainedSolidKick,
    pub solid_drift: FiniteViscoContactDrift,
    pub rigid: Vec<FiniteSolidSustainedRigidStep>,
    pub impact: Option<InelasticMixedContactReport>,
    pub path: P,
}

/// Original physical storage/loss and separate signed/noncancelling numerical errors.
#[derive(Clone, Debug)]
pub struct SustainedSphereTrajectoryReport<P> {
    pub duration_s: f64,
    pub rejected_steps: usize,
    pub steps: Vec<SustainedSphereTrajectoryStep<P>>,
    pub kinetic_change_j: f64,
    pub stored_change_j: f64,
    pub source_work_j: f64,
    pub support_work_j: f64,
    pub relaxation_loss_j: f64,
    /// Impulsive-map loss includes refinable support splitting, not proven friction heat.
    pub unallocated_contact_loss_j: f64,
    pub solid_source_splitting_work_j: f64,
    pub solid_force_splitting_work_j: f64,
    pub material_endpoint_defect_j: f64,
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

/// Original Tet4 receipt names retain their concrete geometric path type.
pub type ViscoSustainedSphereTrajectoryStep =
    SustainedSphereTrajectoryStep<FiniteSolidSpherePathReport>;
pub type ViscoSustainedSphereTrajectoryReport =
    SustainedSphereTrajectoryReport<FiniteSolidSpherePathReport>;

impl<P> Default for SustainedSphereTrajectoryReport<P> {
    fn default() -> Self {
        Self {
            duration_s: 0.0,
            rejected_steps: 0,
            steps: Vec::new(),
            kinetic_change_j: 0.0,
            stored_change_j: 0.0,
            source_work_j: 0.0,
            support_work_j: 0.0,
            relaxation_loss_j: 0.0,
            unallocated_contact_loss_j: 0.0,
            solid_source_splitting_work_j: 0.0,
            solid_force_splitting_work_j: 0.0,
            material_endpoint_defect_j: 0.0,
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
}

#[derive(Clone, Copy)]
pub(super) struct SustainedTrajectoryAuditControls {
    pub(super) linear_tolerance_kg_m_s: f64,
    pub(super) angular_tolerance_kg_m2_s: f64,
    pub(super) energy_tolerance_j: f64,
    pub(super) integration_tolerance_j: f64,
}

impl From<ViscoSustainedSphereTrajectoryControls> for SustainedTrajectoryAuditControls {
    fn from(c: ViscoSustainedSphereTrajectoryControls) -> Self {
        Self {
            linear_tolerance_kg_m_s: c.linear_tolerance_kg_m_s,
            angular_tolerance_kg_m2_s: c.angular_tolerance_kg_m2_s,
            energy_tolerance_j: c.energy_tolerance_j,
            integration_tolerance_j: c.integration_tolerance_j,
        }
    }
}

#[derive(Debug)]
pub struct ViscoSustainedSphereTrajectoryCandidate<'a> {
    pub solid: FiniteViscoTetMotion<'a>,
    pub rigid_bodies: Vec<RigidBodyState>,
    pub report: ViscoSustainedSphereTrajectoryReport,
}

#[derive(Debug)]
pub enum ViscoSustainedSphereRefusal {
    Path(FiniteSolidSpherePathRefusal),
    Integration { defect_j: f64, allowance_j: f64 },
    SolidTrial(FiniteViscoTetError),
    RotationalSolve { rigid_body: usize },
}

/// Failures retain original source/solver details and publish no partial motion.
#[derive(Debug)]
pub enum ViscoSustainedSphereTrajectoryError {
    Invalid(&'static str),
    Solid(FiniteViscoTetError),
    Surface(FiniteSolidSurfaceError),
    Scene(FiniteSolidSphereSceneError),
    Rigid(DynamicsError),
    Impact(InelasticMixedContactError<FiniteViscoTetError>),
    Path(FiniteSolidSphereTrajectoryError),
    RejectionBudget {
        elapsed_s: f64,
        last: ViscoSustainedSphereRefusal,
    },
    PhaseAudit {
        rigid_body: usize,
        linear_residual_kg_m_s: f64,
        angular_residual_kg_m2_s: f64,
        energy_residual_j: f64,
    },
    Cancelled,
    Audit(Box<ViscoSustainedSphereTrajectoryReport>),
}
impl core::fmt::Display for ViscoSustainedSphereTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native sustained viscoelastic sphere contact: {self:?}")
    }
}
impl std::error::Error for ViscoSustainedSphereTrajectoryError {}
type Error = ViscoSustainedSphereTrajectoryError;
type Refusal = ViscoSustainedSphereRefusal;

/// Advance loaded contact with actual Piola forces, native memory and full rigid inertia.
///
/// Supports hold their current displacement and must initially have zero velocity.
/// Dead nodal/reference-volume loads stay fixed; rigid gravity does not also load
/// the solid. Each interval applies original forces, one simultaneous mixed contact
/// solve, then original kinetic drift and Maxwell evolution. Restitution targets
/// are fixed from the pre-force accepted velocities throughout all retries.
///
/// Complete affine-path admission covers both kinds of contact and tetrahedron
/// orientation. Positive activation gaps stay visible. The final endpoint is
/// nonpenetrating; a newly active endpoint row is resolved on continuation.
/// This is a first-order numerical split, not true-flow CCD, general Moreau--Jean
/// dynamics, a thermal model or a guarantee of progress at grazing/Zeno contact.
/// Unresolved geometry, native determinacy, phase errors and budgets never repair
/// a state or publish a partial trajectory. Physical relaxation, impact-map loss
/// and all split defects remain separate, without inferred heat allocation.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn visco_sustained_sphere_trajectory_candidate<'a>(
    solid: &FiniteViscoTetMotion<'a>,
    surface: &FiniteSolidContactSurface,
    scene: &FiniteSolidSphereScene,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    held_velocity_dofs: &[usize],
    nodal_forces_n: &[(usize, f64)],
    body_force_n_m3: [f64; 3],
    wrenches: &[Wrench],
    gravity: Gravity,
    duration_s: f64,
    origin_m: [f64; 3],
    controls: ViscoSustainedSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<ViscoSustainedSphereTrajectoryCandidate<'a>, Error> {
    poll(cx)?;
    validate(
        scene, bodies, properties, wrenches, duration_s, origin_m, controls,
    )?;
    let ndofs = solid.velocity_mps().len() * 3;
    if held_velocity_dofs.len() > ndofs {
        return Err(Error::Invalid("too many held solid DOFs"));
    }
    let mut held = held_velocity_dofs.to_vec();
    held.sort_unstable();
    if held.windows(2).any(|w| w[0] == w[1])
        || held
            .iter()
            .any(|&d| d >= ndofs || solid.velocity_mps()[d / 3][d % 3] != 0.0)
    {
        return Err(Error::Invalid(
            "distinct valid stationary held velocities required",
        ));
    }
    let start_time = solid.solid().time_s();
    let end_time = start_time + duration_s;
    if !end_time.is_finite() || end_time <= start_time {
        return Err(no_progress(0.0));
    }
    let clock_tolerance = 8.0 * f64::EPSILON * start_time.abs().max(end_time.abs()).max(duration_s);
    let mut current: Option<FiniteViscoTetMotion<'a>> = None;
    let mut rigid = bodies.to_vec();
    let mut geometry = scene
        .visco_contacts(surface, solid, &rigid, controls.surface, cx)
        .map_err(Error::Scene)?;
    let mut report = ViscoSustainedSphereTrajectoryReport {
        duration_s,
        ..Default::default()
    };
    let mut elapsed = 0.0;
    let mut proposed = controls.maximum_step_s;
    while elapsed < duration_s {
        poll(cx)?;
        if report.steps.len() == controls.maximum_steps {
            return Err(Error::Path(FiniteSolidSphereTrajectoryError::StepBudget {
                elapsed_s: elapsed,
            }));
        }
        let accepted = current.as_ref().unwrap_or(solid);
        let (contacts, incoming, targets) =
            incoming_laws(accepted, &rigid, properties, &geometry.contacts, cx)?;
        let start = positions(surface, accepted, cx)?;
        let mut dt = proposed.min(duration_s - elapsed);
        let (next, next_rigid, next_geometry, mut step) = loop {
            poll(cx)?;
            if !dt.is_finite()
                || dt <= 0.0
                || elapsed + dt <= elapsed
                || accepted.solid().time_s() + dt <= accepted.solid().time_s()
            {
                return Err(no_progress(elapsed));
            }
            let allowance = controls.integration_tolerance_j * (dt / duration_s);
            if !allowance.is_finite() || allowance <= 0.0 {
                return Err(no_progress(elapsed));
            }
            match interval(
                accepted,
                surface,
                scene,
                &start,
                &rigid,
                properties,
                &held,
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
            )? {
                Ok(value) => break value,
                Err(last) => {
                    if report.rejected_steps == controls.maximum_rejections {
                        return Err(Error::RejectionBudget {
                            elapsed_s: elapsed,
                            last,
                        });
                    }
                    report.rejected_steps += 1;
                    let fraction = match &last {
                        Refusal::Path(p) => p.start(),
                        _ => 0.0,
                    };
                    dt *= if fraction > 0.0 { 0.9 * fraction } else { 0.5 };
                }
            }
        };
        let next_elapsed = elapsed + dt;
        if next_elapsed <= elapsed || next_elapsed > duration_s {
            return Err(no_progress(elapsed));
        }
        elapsed = next_elapsed;
        check_clock(elapsed, next.solid().time_s() - start_time, clock_tolerance)?;
        step.elapsed_s = elapsed;
        step.maximum_contact_gap_m = geometry.maximum_contact_gap_m;
        step.contact_keys = geometry.keys;
        step.incoming_normal_speeds_m_s = incoming;
        step.normal_targets_m_s = targets;
        let defect = step.solid_drift.absolute_integration_defect_j
            + step
                .rigid
                .iter()
                .map(|r| r.force_splitting_work_j.abs() + r.kinetic_drift_defect_j.abs())
                .sum::<f64>();
        proposed = if defect <= 0.25 * controls.integration_tolerance_j * (dt / duration_s)
            && step.solid_drift.absolute_integration_defect_j
                <= 0.25 * controls.solid.maximum_integration_defect_j
        {
            (2.0 * dt).min(controls.maximum_step_s)
        } else {
            dt
        };
        accumulate(&mut report, &step, vector(origin_m));
        report.steps.push(step);
        finish(&mut report);
        if !admitted(&report, controls.into()) {
            return Err(Error::Audit(Box::new(report)));
        }
        current = Some(next);
        rigid = next_rigid;
        geometry = next_geometry;
    }
    let solid = current.expect("positive interval accepted a complete native split");
    check_clock(
        duration_s,
        solid.solid().time_s() - start_time,
        clock_tolerance,
    )?;
    poll(cx)?;
    Ok(ViscoSustainedSphereTrajectoryCandidate {
        solid,
        rigid_bodies: rigid,
        report,
    })
}

type Attempt<'a> = Result<
    (
        FiniteViscoTetMotion<'a>,
        Vec<RigidBodyState>,
        InelasticSphereSceneQuery,
        ViscoSustainedSphereTrajectoryStep,
    ),
    Refusal,
>;

#[allow(clippy::too_many_arguments)]
fn interval<'a>(
    solid: &FiniteViscoTetMotion<'a>,
    surface: &FiniteSolidContactSurface,
    scene: &FiniteSolidSphereScene,
    start: &[[f64; 3]],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    held: &[usize],
    contacts: &[InelasticMixedContact],
    nodal_forces: &[(usize, f64)],
    body_force: [f64; 3],
    wrenches: &[Wrench],
    gravity: Gravity,
    dt: f64,
    origin: [f64; 3],
    c: ViscoSustainedSphereTrajectoryControls,
    allowance: f64,
    clock_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<Attempt<'a>, Error> {
    let (kicked, kick) = solid
        .force_kick_candidate(
            dt,
            held,
            nodal_forces,
            body_force,
            origin,
            c.impact.solid_impulse,
            cx,
        )
        .map_err(Error::Solid)?;
    let mut kicks = Vec::with_capacity(bodies.len());
    let mut kicked_bodies = Vec::with_capacity(bodies.len());
    for ((&body, &property), &wrench) in bodies.iter().zip(properties).zip(wrenches) {
        poll(cx)?;
        let k = rigid_force_kick(body, property, wrench, gravity, dt).map_err(Error::Rigid)?;
        kicked_bodies.push(k.state_after);
        kicks.push(k);
    }
    let (contacted, contacted_bodies, impact) = if contacts.is_empty() {
        (kicked, kicked_bodies, None)
    } else {
        let event = visco_mixed_contact_candidate(
            &kicked,
            held,
            &kicked_bodies,
            properties,
            contacts,
            origin,
            c.impact,
            cx,
        )
        .map_err(Error::Impact)?;
        (event.solid, event.rigid_bodies, Some(event.report))
    };
    let mut drift_controls = c.solid;
    drift_controls.maximum_integration_defect_j =
        drift_controls.maximum_integration_defect_j.min(allowance);
    let (next, solid_drift) = match contacted.contact_drift_candidate(&kick, drift_controls, cx) {
        Ok(value) => value,
        Err(error) if refinable(&error) => return Ok(Err(Refusal::SolidTrial(error))),
        Err(error) => return Err(Error::Solid(error)),
    };
    check_clock(
        dt,
        next.solid().time_s() - solid.solid().time_s(),
        clock_tolerance,
    )?;
    let mut next_bodies = Vec::with_capacity(bodies.len());
    let mut rigid = Vec::with_capacity(bodies.len());
    let mut absolute_defect = solid_drift.absolute_integration_defect_j;
    for (index, (((&body, &property), &wrench), kick)) in contacted_bodies
        .iter()
        .zip(properties)
        .zip(wrenches)
        .zip(kicks)
        .enumerate()
    {
        poll(cx)?;
        let step =
            match rigid_contact_drift(kick, body, property, wrench, gravity, dt, vector(origin)) {
                Ok(value) => value,
                Err(DynamicsError::RotationalSolveFailed) => {
                    return Ok(Err(Refusal::RotationalSolve { rigid_body: index }));
                }
                Err(error) => return Err(Error::Rigid(error)),
            };
        let (linear, angular, energy) = step.phase_residuals();
        if [
            linear,
            angular,
            energy,
            step.kinetic_drift_defect_j,
            step.force_displacement_work_j,
            step.force_splitting_work_j,
            step.torque_impulse_work_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || linear > c.impact.linear_tolerance_kg_m_s
            || angular > c.impact.angular_tolerance_kg_m2_s
            || energy > c.impact.energy_tolerance_j
        {
            return Err(Error::PhaseAudit {
                rigid_body: index,
                linear_residual_kg_m_s: linear,
                angular_residual_kg_m2_s: angular,
                energy_residual_j: energy,
            });
        }
        absolute_defect += step.force_splitting_work_j.abs() + step.kinetic_drift_defect_j.abs();
        next_bodies.push(step.drift.state_after);
        rigid.push(step);
    }
    if !absolute_defect.is_finite() {
        return Err(Error::Invalid("unrepresentable absolute split defect"));
    }
    if absolute_defect > allowance {
        return Ok(Err(Refusal::Integration {
            defect_j: absolute_defect,
            allowance_j: allowance,
        }));
    }
    let end = positions(surface, &next, cx)?;
    let path = match admit_path_at_positions(
        surface,
        solid.solid().reference_tetrahedra(),
        start,
        &end,
        &contacted_bodies,
        &next_bodies,
        scene.spheres(),
        true,
        c.path,
        cx,
    )
    .map_err(Error::Path)?
    {
        Ok(path) => path,
        Err(reason) => return Ok(Err(Refusal::Path(reason))),
    };
    let geometry = scene
        .visco_contacts(surface, &next, &next_bodies, c.surface, cx)
        .map_err(Error::Scene)?;
    poll(cx)?;
    Ok(Ok((
        next,
        next_bodies,
        geometry,
        ViscoSustainedSphereTrajectoryStep {
            elapsed_s: 0.0,
            duration_s: dt,
            maximum_contact_gap_m: 0.0,
            contact_keys: Vec::new(),
            incoming_normal_speeds_m_s: Vec::new(),
            normal_targets_m_s: Vec::new(),
            solid_kick: kick.into(),
            solid_drift,
            rigid,
            impact,
            path,
        },
    )))
}

fn incoming_laws(
    solid: &FiniteViscoTetMotion<'_>,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact],
    cx: &Cx<'_>,
) -> Result<(Vec<InelasticMixedContact>, Vec<f64>, Vec<f64>), Error> {
    let mut mapped = contacts.to_vec();
    let mut incoming = Vec::with_capacity(contacts.len());
    let mut targets = Vec::with_capacity(contacts.len());
    for c in &mut mapped {
        poll(cx)?;
        let axes = contact_frame(c.normal_world, c.tangent_world).map_err(Error::Invalid)?;
        let first = bodies[c.point.body]
            .point_kinematics(properties[c.point.body], c.point.arm_body)
            .map_err(Error::Rigid)?
            .point_velocity_world;
        let second = match c.counterpart {
            FiniteSolidContactCounterpart::Solid(p) => vector(
                solid
                    .point_kinematics(p, cx)
                    .map_err(Error::Solid)?
                    .velocity_mps,
            ),
            FiniteSolidContactCounterpart::Rigid(p) => {
                bodies[p.body]
                    .point_kinematics(properties[p.body], p.arm_body)
                    .map_err(Error::Rigid)?
                    .point_velocity_world
            }
        };
        let speed = axes[0].dot(first.sub(second));
        let target = match c.normal_law {
            NormalVelocityLaw::NewtonRestitution { coefficient }
                if coefficient.is_finite() && (0.0..=1.0).contains(&coefficient) =>
            {
                -coefficient * speed.min(0.0)
            }
            NormalVelocityLaw::Target { speed } if speed.is_finite() && speed >= 0.0 => speed,
            _ => return Err(Error::Invalid("valid incoming normal law required")),
        };
        if !speed.is_finite() || !target.is_finite() {
            return Err(Error::Invalid("unrepresentable pre-force contact speed"));
        }
        incoming.push(speed);
        targets.push(target);
        c.normal_law = NormalVelocityLaw::Target { speed: target };
    }
    Ok((mapped, incoming, targets))
}

pub(super) fn accumulate<P>(
    r: &mut SustainedSphereTrajectoryReport<P>,
    s: &SustainedSphereTrajectoryStep<P>,
    origin: Vec3,
) {
    let k = &s.solid_kick;
    let d = &s.solid_drift;
    let native = &k.impulse;
    r.kinetic_change_j += native.kinetic_change_j + d.kinetic_change_j;
    r.stored_change_j += d.stored_after_j - d.stored_before_j;
    r.relaxation_loss_j += d.dissipated_j;
    r.source_work_j += d.source_displacement_work_j;
    r.support_work_j += native.support_work_j;
    r.solid_source_splitting_work_j += d.source_splitting_work_j;
    r.solid_force_splitting_work_j += d.force_splitting_work_j;
    r.material_endpoint_defect_j += d.material_endpoint_defect_j;
    r.integration_defect_j += d.integration_defect_j;
    r.absolute_integration_defect_j += d.absolute_integration_defect_j;
    add(
        &mut r.linear_momentum_change_kg_m_s,
        vector(native.applied.linear_kg_m_s)
            .add(vector(native.support.linear_kg_m_s))
            .add(vector(native.linear_residual_kg_m_s))
            .add(vector(d.linear_residual_kg_m_s)),
    );
    add(
        &mut r.angular_momentum_change_kg_m2_s,
        vector(native.applied.angular_kg_m2_s)
            .add(vector(native.support.angular_kg_m2_s))
            .add(vector(native.angular_residual_kg_m2_s))
            .add(vector(d.angular_residual_kg_m2_s)),
    );
    add(
        &mut r.external_impulse_kg_m_s,
        vector(k.source_impulse.linear_kg_m_s),
    );
    add(
        &mut r.external_angular_impulse_kg_m2_s,
        vector(k.source_impulse.angular_kg_m2_s),
    );
    add(
        &mut r.support_impulse_kg_m_s,
        vector(native.support.linear_kg_m_s),
    );
    add(
        &mut r.support_angular_impulse_kg_m2_s,
        vector(native.support.angular_kg_m2_s),
    );
    r.accumulated_linear_residual_kg_m_s += norm(vector(native.linear_residual_kg_m_s))
        + norm(vector(k.linear_decomposition_residual_kg_m_s))
        + norm(vector(k.internal_impulse.linear_kg_m_s))
        + norm(vector(d.linear_residual_kg_m_s));
    r.accumulated_angular_residual_kg_m2_s += norm(vector(native.angular_residual_kg_m2_s))
        + norm(vector(k.angular_decomposition_residual_kg_m2_s))
        + norm(vector(k.internal_impulse.angular_kg_m2_s))
        + norm(vector(d.angular_residual_kg_m2_s));
    r.accumulated_energy_residual_j += native.work_residual_j.abs()
        + k.work_decomposition_residual_j.abs()
        + d.material_balance_error_j.abs();
    if let Some(i) = &s.impact {
        r.kinetic_change_j += i.kinetic_change_j;
        r.support_work_j += i.support_work_j;
        r.unallocated_contact_loss_j += i.contact.dissipated_energy;
        let support = vector(i.solid_impulse.support.linear_kg_m_s);
        let moment = vector(i.solid_impulse.support.angular_kg_m2_s);
        add(
            &mut r.linear_momentum_change_kg_m_s,
            support.add(i.linear_momentum_residual_kg_m_s),
        );
        add(
            &mut r.angular_momentum_change_kg_m2_s,
            moment.add(i.angular_momentum_residual_kg_m2_s),
        );
        add(&mut r.support_impulse_kg_m_s, support);
        add(&mut r.support_angular_impulse_kg_m2_s, moment);
        r.accumulated_linear_residual_kg_m_s +=
            i.absolute_linear_residual_kg_m_s + norm(i.linear_momentum_residual_kg_m_s);
        r.accumulated_angular_residual_kg_m2_s +=
            i.absolute_angular_residual_kg_m2_s + norm(i.angular_momentum_residual_kg_m2_s);
        r.accumulated_energy_residual_j += i.absolute_work_residual_j
            + i.work_energy_residual_j.abs()
            + i.dissipation_residual_j.abs()
            + (i.relative_impulse_work_j - i.contact.relative_work).abs();
    }
    for b in &s.rigid {
        let old = b.kick.state_before;
        let arm = old.pose().position_world().sub(origin);
        let rotation = old.pose().orientation();
        let applied_moment = arm
            .cross(b.kick.linear_impulse_world)
            .add(rotation.rotate_body_to_world(b.kick.angular_impulse_body));
        let delta_moment = arm
            .cross(b.kick.linear_increment_world)
            .add(rotation.rotate_body_to_world(b.kick.angular_increment_body));
        r.kinetic_change_j += b.kick.work.kinetic_energy_change_j + b.kinetic_drift_defect_j;
        r.source_work_j += b.force_displacement_work_j + b.torque_impulse_work_j;
        r.rigid_force_splitting_work_j += b.force_splitting_work_j;
        r.rigid_drift_energy_defect_j += b.kinetic_drift_defect_j;
        r.integration_defect_j += b.force_splitting_work_j + b.kinetic_drift_defect_j;
        r.absolute_integration_defect_j +=
            b.force_splitting_work_j.abs() + b.kinetic_drift_defect_j.abs();
        add(
            &mut r.linear_momentum_change_kg_m_s,
            b.kick
                .linear_increment_world
                .add(b.drift_linear_increment_kg_m_s),
        );
        add(
            &mut r.angular_momentum_change_kg_m2_s,
            delta_moment.add(b.drift_angular_increment_kg_m2_s),
        );
        add(&mut r.external_impulse_kg_m_s, b.kick.linear_impulse_world);
        add(&mut r.external_angular_impulse_kg_m2_s, applied_moment);
        let (linear, angular, energy) = b.phase_residuals();
        r.accumulated_linear_residual_kg_m_s += linear;
        r.accumulated_angular_residual_kg_m2_s += angular;
        r.accumulated_energy_residual_j += energy;
    }
}

pub(super) fn finish<P>(r: &mut SustainedSphereTrajectoryReport<P>) {
    r.work_energy_residual_j =
        r.kinetic_change_j + r.stored_change_j + r.relaxation_loss_j + r.unallocated_contact_loss_j
            - r.source_work_j
            - r.support_work_j
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
pub(super) fn admitted<P>(
    r: &SustainedSphereTrajectoryReport<P>,
    c: SustainedTrajectoryAuditControls,
) -> bool {
    [
        r.kinetic_change_j,
        r.stored_change_j,
        r.source_work_j,
        r.support_work_j,
        r.relaxation_loss_j,
        r.unallocated_contact_loss_j,
        r.solid_source_splitting_work_j,
        r.solid_force_splitting_work_j,
        r.material_endpoint_defect_j,
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

#[allow(clippy::too_many_arguments)]
fn validate(
    scene: &FiniteSolidSphereScene,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    wrenches: &[Wrench],
    duration: f64,
    origin: [f64; 3],
    c: ViscoSustainedSphereTrajectoryControls,
) -> Result<(), Error> {
    if bodies.len() != properties.len()
        || bodies.len() != wrenches.len()
        || bodies.len() != scene.spheres().len()
        || bodies.len() > c.surface.maximum_spheres
        || c.impact.maximum_rigid_bodies == 0
        || bodies.len() > c.impact.maximum_rigid_bodies
        || !(1..=21).contains(&c.impact.contact.maximum_contacts)
        || !(1..=65_536).contains(&c.maximum_steps)
        || c.maximum_rejections > 65_536
        || !c.path.valid()
        || origin.iter().any(|v| !v.is_finite())
        || c.surface.activation_tolerance_m > c.impact.geometry_tolerance_m
        || [
            duration,
            c.maximum_step_s,
            c.surface.activation_tolerance_m,
            c.impact.geometry_tolerance_m,
            c.impact.linear_tolerance_kg_m_s,
            c.impact.angular_tolerance_kg_m2_s,
            c.impact.energy_tolerance_j,
            c.linear_tolerance_kg_m_s,
            c.angular_tolerance_kg_m2_s,
            c.energy_tolerance_j,
            c.integration_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || [
            c.solid.material_balance_tolerance_j,
            c.solid.maximum_integration_defect_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(Error::Invalid(
            "bounded matching scene and finite positive physical/time gates required",
        ));
    }
    Ok(())
}
fn refinable(error: &FiniteViscoTetError) -> bool {
    matches!(
        error,
        FiniteViscoTetError::EnergyGate { .. }
            | FiniteViscoTetError::Material {
                source: MaterialError::State { .. },
                ..
            }
            | FiniteViscoTetError::Native(HyperTetError::MaterialRefused { .. })
    )
}
fn positions(
    surface: &FiniteSolidContactSurface,
    solid: &FiniteViscoTetMotion<'_>,
    cx: &Cx<'_>,
) -> Result<Vec<[f64; 3]>, Error> {
    let state = solid.solid();
    surface
        .current_positions_from_reference(
            state.reference_nodes_m(),
            state.reference_tetrahedra(),
            state.displacement_m(),
            cx,
        )
        .map_err(Error::Surface)
}
fn check_clock(requested: f64, represented: f64, tolerance: f64) -> Result<(), Error> {
    if (represented - requested).abs() > tolerance {
        Err(Error::Path(
            FiniteSolidSphereTrajectoryError::TimeResolution {
                requested_step_s: requested,
                represented_step_s: represented,
            },
        ))
    } else {
        Ok(())
    }
}
fn no_progress(elapsed_s: f64) -> Error {
    Error::Path(FiniteSolidSphereTrajectoryError::NoProgress { elapsed_s })
}
fn add(target: &mut Vec3, value: Vec3) {
    *target = target.add(value);
}
fn poll(cx: &Cx<'_>) -> Result<(), Error> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
