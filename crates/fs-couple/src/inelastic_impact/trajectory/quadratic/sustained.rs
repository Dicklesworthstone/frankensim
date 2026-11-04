//! Loaded curved contact through original quadratic Maxwell force and drift owners.
//!
//! Incoming restitution is fixed before external forces act. Each retry repeats
//! the entire force/contact/drift map from the accepted state. Pressure solves,
//! quadrature histories and physical storage remain with the native mixed owner;
//! geometric and splitting refusals publish neither partial motion nor memory.

use crate::finite_solid::contact::{contact_frame, step::FiniteSolidSphereTrajectoryError};
use crate::inelastic_impact::{
    multiple::{
        FiniteSolidContactCounterpart, InelasticMixedContact, InelasticMixedContactControls,
        InelasticMixedContactError,
    },
    surface::{
        quadratic::{
            QuadraticSolidContactSurface, QuadraticSolidSurfaceError,
            path::QuadraticSolidPathControls,
            scene::{
                QuadraticSphereSceneError, QuadraticSphereScenePathReport,
                QuadraticSphereScenePathVerdict, QuadraticSphereSceneQuery,
            },
            volume::{
                PreparedQuadraticSolidVolume, QuadraticSolidVolumeError,
                QuadraticSolidVolumePathReport, QuadraticSolidVolumePathVerdict,
            },
        },
        scene::FiniteSolidSphereScene,
    },
    trajectory::{
        check_clock,
        sustained::{
            SustainedSphereTrajectoryReport, SustainedSphereTrajectoryStep,
            SustainedTrajectoryAuditControls, ViscoSustainedSolidKick, accumulate, admitted,
            finish,
        },
        vector,
    },
};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::Cx;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_mbd::{DynamicsError, Gravity, MassProperties, RigidBodyState, Wrench};
use fs_query::tetrahedron::quadratic::QuadraticTetPathControls;
use fs_solid::hyper3::mixed::{
    MixedUpSettings,
    visco::{
        MixedViscoTetError,
        dynamics::{
            FiniteViscoContactDriftControls, MixedViscoForceKick, MixedViscoTetMotion,
            Tet10ImpulsePoint,
        },
    },
};

mod interval;

/// Native phase gates, complete geometric paths and whole-call refinement limits.
#[derive(Clone, Copy, Debug)]
pub struct MixedViscoSustainedSphereTrajectoryControls {
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    /// Curved boundary and rigid-pair sweeps share their authored work limits.
    pub path: QuadraticSolidPathControls,
    /// Closed material-domain admission across all original TET10 elements.
    pub volume: QuadraticTetPathControls,
    /// Also supplies independent native force-kick and rigid phase gates.
    pub impact: InelasticMixedContactControls,
    /// The original finite-compliance pressure solve at fixed drift coordinates.
    pub settings: MixedUpSettings,
    pub solid: FiniteViscoContactDriftControls,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
    /// Noncancelling whole-call allowance; each trial receives dt/duration of it.
    pub integration_tolerance_j: f64,
}

/// Complete material and contact path evidence for one accepted native interval.
#[derive(Clone, Debug)]
pub struct MixedViscoSustainedSpherePathReport {
    pub volume: QuadraticSolidVolumePathReport,
    pub contact: QuadraticSphereScenePathReport,
}

pub type MixedViscoSustainedSphereTrajectoryStep =
    SustainedSphereTrajectoryStep<MixedViscoSustainedSpherePathReport>;
pub type MixedViscoSustainedSphereTrajectoryReport =
    SustainedSphereTrajectoryReport<MixedViscoSustainedSpherePathReport>;

/// Adopt every owner together only after the complete requested interval passes.
#[derive(Debug)]
pub struct MixedViscoSustainedSphereTrajectoryCandidate<'a> {
    pub solid: MixedViscoTetMotion<'a>,
    pub rigid_bodies: Vec<RigidBodyState>,
    pub reference_volume: QuadraticSolidVolumePathReport,
    pub initial_volume: QuadraticSolidVolumePathReport,
    pub report: MixedViscoSustainedSphereTrajectoryReport,
}

/// Failed private trial, retained when the authored refinement limit is reached.
#[derive(Debug)]
pub enum MixedViscoSustainedSphereRefusal {
    Volume(QuadraticSolidVolumePathReport),
    Path(QuadraticSphereScenePathReport),
    Endpoint(QuadraticSphereSceneError),
    Integration { defect_j: f64, allowance_j: f64 },
    SolidTrial(MixedViscoTetError),
    RotationalSolve { rigid_body: usize },
}

/// Every refusal leaves caller-owned mechanics, pressure and histories unchanged.
#[derive(Debug)]
pub enum MixedViscoSustainedSphereTrajectoryError {
    Invalid(&'static str),
    Solid(MixedViscoTetError),
    Surface(QuadraticSolidSurfaceError),
    Scene(QuadraticSphereSceneError),
    Volume(QuadraticSolidVolumeError),
    InitialVolumeUnresolved {
        report: QuadraticSolidVolumePathReport,
    },
    Rigid(DynamicsError),
    Impact(InelasticMixedContactError<MixedViscoTetError>),
    Path(FiniteSolidSphereTrajectoryError),
    RejectionBudget {
        elapsed_s: f64,
        last: Box<MixedViscoSustainedSphereRefusal>,
    },
    PhaseAudit {
        rigid_body: usize,
        linear_residual_kg_m_s: f64,
        angular_residual_kg_m2_s: f64,
        energy_residual_j: f64,
    },
    Cancelled,
    Audit(Box<MixedViscoSustainedSphereTrajectoryReport>),
}
impl core::fmt::Display for MixedViscoSustainedSphereTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native sustained quadratic Maxwell contact: {self:?}")
    }
}
impl std::error::Error for MixedViscoSustainedSphereTrajectoryError {}
type Error = MixedViscoSustainedSphereTrajectoryError;
type Refusal = MixedViscoSustainedSphereRefusal;
type Contact = InelasticMixedContact<Tet10ImpulsePoint>;

struct AcceptedInterval<'a> {
    solid: MixedViscoTetMotion<'a>,
    rigid: Vec<RigidBodyState>,
    positions: Vec<[f64; 3]>,
    geometry: QuadraticSphereSceneQuery,
    step: MixedViscoSustainedSphereTrajectoryStep,
}

/// Advance a complete loaded sphere scene against the original curved P2/P1 solid.
///
/// Held displacement DOFs must initially be stationary. Authored nodal and
/// reference-volume loads, rigid wrenches and rigid gravity remain fixed for the
/// call. Solid gravity is an explicit native load. All active curved-solid and
/// declared rigid-pair rows enter one simultaneous native contact solve, with
/// restitution targets from the accepted PRE-FORCE velocities.
///
/// The original reference and current volume maps are admitted before forces or
/// impulses. Every accepted interval then admits the complete declared affine
/// nodal/centre paths, including volume, curved boundary and sphere pairs. The
/// original pressure-only endpoint solve updates finite-compliance pressure at
/// the drift coordinates; its coordinate-fixing reactions are not physical
/// supports. Only successful intervals advance native material history or time.
///
/// Numerical source/force/material/rigid split defects have independent signed
/// and noncancelling accounts. Impact-map loss may contain refinable support
/// splitting; it is not inferred friction heat. This is a first-order split
/// requiring refinement, with local regularity of its declared numerical path,
/// not continuum-flow CCD or global self-intersection exclusion. Exact
/// incompressibility requires a different compatible kinetic subflow and refuses.
/// Grazing/Zeno motion, unresolved geometry and exhausted allowances may refuse.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn mixed_visco_sustained_sphere_trajectory_candidate<'a>(
    solid: &MixedViscoTetMotion<'a>,
    surface: &QuadraticSolidContactSurface,
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
    controls: MixedViscoSustainedSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<MixedViscoSustainedSphereTrajectoryCandidate<'a>, Error> {
    poll(cx)?;
    validate(
        scene, bodies, properties, wrenches, duration_s, origin_m, controls,
    )?;
    if solid.volumetric_treatment() == VolumetricTreatment::ExactIncompressible {
        return Err(Error::Solid(MixedViscoTetError::ExactIncompressibleSplit));
    }
    let ndofs = solid.velocity_mps().len() * 3;
    if held_velocity_dofs.len() > ndofs {
        return Err(Error::Invalid("too many held quadratic solid DOFs"));
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
    let prepared_volume = surface
        .prepare_volume_paths(controls.volume, cx)
        .map_err(Error::Volume)?;
    let reference_volume = *prepared_volume.reference_admission();
    let mut coordinates = positions(surface, solid, cx)?;
    let initial_volume = prepared_volume
        .affine_path_at_positions(&coordinates, &coordinates, controls.volume, cx)
        .map_err(Error::Volume)?;
    if matches!(
        initial_volume.verdict,
        QuadraticSolidVolumePathVerdict::Unresolved { .. }
    ) {
        return Err(Error::InitialVolumeUnresolved {
            report: initial_volume,
        });
    }
    let mut current: Option<MixedViscoTetMotion<'a>> = None;
    let mut rigid = bodies.to_vec();
    let mut geometry = scene
        .quadratic_contacts_at_positions(surface, &coordinates, &rigid, controls.path.initial, cx)
        .map_err(Error::Scene)?;
    let mut report = MixedViscoSustainedSphereTrajectoryReport {
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
        let mut dt = proposed.min(duration_s - elapsed);
        let mut next = loop {
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
            match interval::advance(
                accepted,
                surface,
                &prepared_volume,
                scene,
                &coordinates,
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
                Ok(next) => break next,
                Err(last) => {
                    if report.rejected_steps == controls.maximum_rejections {
                        return Err(Error::RejectionBudget {
                            elapsed_s: elapsed,
                            last: Box::new(last),
                        });
                    }
                    report.rejected_steps += 1;
                    let fraction = match &last {
                        Refusal::Volume(QuadraticSolidVolumePathReport {
                            verdict:
                                QuadraticSolidVolumePathVerdict::Unresolved { fraction_start, .. },
                            ..
                        }) => *fraction_start,
                        Refusal::Path(QuadraticSphereScenePathReport {
                            verdict:
                                QuadraticSphereScenePathVerdict::PossibleContact {
                                    fraction_start, ..
                                },
                            ..
                        }) => *fraction_start,
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
        check_clock(
            elapsed,
            next.solid.solid().time_s() - start_time,
            clock_tolerance,
        )
        .map_err(Error::Path)?;
        next.step.elapsed_s = elapsed;
        next.step.maximum_contact_gap_m = geometry.maximum_contact_gap_m;
        next.step.contact_keys = geometry.keys;
        next.step.incoming_normal_speeds_m_s = incoming;
        next.step.normal_targets_m_s = targets;
        let defect = next.step.solid_drift.absolute_integration_defect_j
            + next
                .step
                .rigid
                .iter()
                .map(|r| r.force_splitting_work_j.abs() + r.kinetic_drift_defect_j.abs())
                .sum::<f64>();
        proposed = if defect <= 0.25 * controls.integration_tolerance_j * (dt / duration_s)
            && next.step.solid_drift.absolute_integration_defect_j
                <= 0.25 * controls.solid.maximum_integration_defect_j
        {
            (2.0 * dt).min(controls.maximum_step_s)
        } else {
            dt
        };
        accumulate(&mut report, &next.step, vector(origin_m));
        report.steps.push(next.step);
        finish(&mut report);
        if !admitted(&report, controls.into()) {
            return Err(Error::Audit(Box::new(report)));
        }
        current = Some(next.solid);
        rigid = next.rigid;
        coordinates = next.positions;
        geometry = next.geometry;
    }
    let solid = current.expect("positive interval accepted a complete native quadratic split");
    check_clock(
        duration_s,
        solid.solid().time_s() - start_time,
        clock_tolerance,
    )
    .map_err(Error::Path)?;
    poll(cx)?;
    Ok(MixedViscoSustainedSphereTrajectoryCandidate {
        solid,
        rigid_bodies: rigid,
        reference_volume,
        initial_volume,
        report,
    })
}

impl From<MixedViscoForceKick<'_>> for ViscoSustainedSolidKick {
    fn from(k: MixedViscoForceKick<'_>) -> Self {
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

impl From<MixedViscoSustainedSphereTrajectoryControls> for SustainedTrajectoryAuditControls {
    fn from(c: MixedViscoSustainedSphereTrajectoryControls) -> Self {
        Self {
            linear_tolerance_kg_m_s: c.linear_tolerance_kg_m_s,
            angular_tolerance_kg_m2_s: c.angular_tolerance_kg_m2_s,
            energy_tolerance_j: c.energy_tolerance_j,
            integration_tolerance_j: c.integration_tolerance_j,
        }
    }
}

fn incoming_laws(
    solid: &MixedViscoTetMotion<'_>,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[Contact],
    cx: &Cx<'_>,
) -> Result<(Vec<Contact>, Vec<f64>, Vec<f64>), Error> {
    let mut mapped = contacts.to_vec();
    let mut incoming = Vec::with_capacity(contacts.len());
    let mut targets = Vec::with_capacity(contacts.len());
    for contact in &mut mapped {
        poll(cx)?;
        let axes =
            contact_frame(contact.normal_world, contact.tangent_world).map_err(Error::Invalid)?;
        let first = bodies[contact.point.body]
            .point_kinematics(properties[contact.point.body], contact.point.arm_body)
            .map_err(Error::Rigid)?
            .point_velocity_world;
        let second = match contact.counterpart {
            FiniteSolidContactCounterpart::Solid(point) => vector(
                solid
                    .point_kinematics(point, cx)
                    .map_err(Error::Solid)?
                    .velocity_mps,
            ),
            FiniteSolidContactCounterpart::Rigid(point) => {
                bodies[point.body]
                    .point_kinematics(properties[point.body], point.arm_body)
                    .map_err(Error::Rigid)?
                    .point_velocity_world
            }
        };
        let speed = axes[0].dot(first.sub(second));
        let target = match contact.normal_law {
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
        contact.normal_law = NormalVelocityLaw::Target { speed: target };
    }
    Ok((mapped, incoming, targets))
}

#[allow(clippy::too_many_arguments)]
fn validate(
    scene: &FiniteSolidSphereScene,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    wrenches: &[Wrench],
    duration: f64,
    origin: [f64; 3],
    c: MixedViscoSustainedSphereTrajectoryControls,
) -> Result<(), Error> {
    let surface = c.path.initial.surface;
    if bodies.len() != properties.len()
        || bodies.len() != wrenches.len()
        || bodies.len() != scene.spheres().len()
        || bodies.len() > surface.maximum_spheres
        || c.impact.maximum_rigid_bodies == 0
        || bodies.len() > c.impact.maximum_rigid_bodies
        || !(1..=21).contains(&c.impact.contact.maximum_contacts)
        || !(1..=65_536).contains(&c.maximum_steps)
        || c.maximum_rejections > 65_536
        || origin.iter().any(|v| !v.is_finite())
        || surface.activation_tolerance_m > c.impact.geometry_tolerance_m
        || [
            duration,
            c.maximum_step_s,
            surface.activation_tolerance_m,
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
    c.path
        .sweep
        .validate()
        .map_err(|e| Error::Surface(QuadraticSolidSurfaceError::Query(e)))?;
    Ok(())
}

fn positions(
    surface: &QuadraticSolidContactSurface,
    solid: &MixedViscoTetMotion<'_>,
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
fn no_progress(elapsed_s: f64) -> Error {
    Error::Path(FiniteSolidSphereTrajectoryError::NoProgress { elapsed_s })
}
fn poll(cx: &Cx<'_>) -> Result<(), Error> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
