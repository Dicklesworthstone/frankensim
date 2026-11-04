//! Native P2/P1 Maxwell motion with curved impacts and whole-element path admission.
//!
//! Each accepted interval retains the original pressure, quadrature histories,
//! consistent mass, rigid inertia and simultaneous contact solve. The declared
//! collision path interpolates all ten nodal coordinates and sphere centres;
//! it is not an enclosure of the unknown continuous material trajectory.

/// Loaded contact through original quadratic force, contact and material substeps.
pub mod sustained;

use super::{
    SphereTrajectoryReport, SphereTrajectoryStep, TrajectoryAuditControls, accumulate, admitted,
    check_clock, finish, vector,
};
use crate::finite_solid::contact::step::{
    FiniteSolidSphereTrajectoryError,
    free::{FreeFlightError, rigid_free_flight},
};
use crate::inelastic_impact::{
    multiple::{
        InelasticMixedContactControls, InelasticMixedContactError,
        mixed_visco_mixed_contact_candidate,
    },
    surface::{
        FiniteSolidSphere,
        quadratic::{
            QuadraticSolidContactSurface, QuadraticSolidSurfaceError,
            path::{
                QuadraticSolidPathControls, QuadraticSolidPathReport, QuadraticSolidPathVerdict,
            },
            scene::{
                QuadraticSphereSceneError, QuadraticSphereScenePathReport,
                QuadraticSphereScenePathVerdict,
            },
            volume::{
                QuadraticSolidVolumeError, QuadraticSolidVolumePathReport,
                QuadraticSolidVolumePathVerdict,
            },
        },
        scene::{FiniteSolidSphereScene, FiniteSolidSphereSceneKey},
    },
};
use fs_exec::Cx;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_mbd::{Gravity, MassProperties, RigidBodyState, Vec3, Wrench};
use fs_query::tetrahedron::quadratic::QuadraticTetPathControls;
use fs_solid::hyper3::{
    mixed::{
        MixedUpSettings,
        visco::{MixedViscoTetEnergy, MixedViscoTetError, dynamics::MixedViscoTetMotion},
    },
    visco::dynamics::FiniteViscoMotionControls,
};

mod geometry;
mod paths;
use geometry::TrajectoryScene;

/// Explicit native solver, geometric work and complete-call physical allowances.
#[derive(Clone, Copy, Debug)]
pub struct MixedViscoSphereTrajectoryControls {
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    /// Initial/endpoint geometry and complete curved boundary paths. Pair tests
    /// share the query-visit and time-window caps with the boundary path.
    pub path: QuadraticSolidPathControls,
    /// Aggregate closed material-domain Jacobian work over all original elements.
    pub volume: QuadraticTetPathControls,
    pub impact: InelasticMixedContactControls,
    pub settings: MixedUpSettings,
    pub solid: FiniteViscoMotionControls,
    /// Gates apply to summed residual magnitudes as well as signed final balances.
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
    /// Sum of numerical kinetic loss and absolute material/rigid work defects.
    pub integration_tolerance_j: f64,
    /// Sum of native BE angular-defect magnitudes, kept separate from torque.
    pub angular_integration_tolerance_kg_m2_s: f64,
}

/// Complete-volume, curved-boundary and rigid-pair path checks.
#[derive(Clone, Debug, PartialEq)]
pub struct MixedViscoSpherePathReport {
    /// Original active features in the same order as this interval's impulses.
    pub contact_keys: Vec<FiniteSolidSphereSceneKey>,
    pub volume: QuadraticSolidVolumePathReport,
    pub boundary: QuadraticSolidPathReport,
    pub sphere_pair_comparisons: usize,
    pub sphere_pair_candidates: usize,
    pub sphere_pairs_checked: usize,
    pub sphere_pair_windows: usize,
    /// No pair lower bound exists when there are fewer than two spheres.
    pub minimum_pair_clearance_m: Option<f64>,
}

pub type MixedViscoSphereTrajectoryStep =
    SphereTrajectoryStep<MixedViscoTetEnergy, MixedViscoSpherePathReport>;
pub type MixedViscoSphereTrajectoryReport =
    SphereTrajectoryReport<MixedViscoTetEnergy, MixedViscoSpherePathReport>;

/// All owners can be adopted together only after the complete interval passes.
#[derive(Debug)]
pub struct MixedViscoSphereTrajectoryCandidate<'a> {
    pub solid: MixedViscoTetMotion<'a>,
    pub rigid_bodies: Vec<RigidBodyState>,
    /// The complete original reference map was proved once before advancing.
    pub reference_volume: QuadraticSolidVolumePathReport,
    /// The complete accepted current map was admitted before its first impulse.
    pub initial_volume: QuadraticSolidVolumePathReport,
    pub report: MixedViscoSphereTrajectoryReport,
}

/// Conservative retry evidence. Unresolved geometry is never repaired.
#[derive(Debug)]
pub enum MixedViscoSpherePathRefusal {
    Volume(QuadraticSolidVolumePathReport),
    Boundary(QuadraticSolidPathReport),
    ScenePath(QuadraticSphereScenePathReport),
    Endpoint(QuadraticSolidSurfaceError),
}

/// Every refusal leaves caller-owned mechanics, pressure and histories unchanged.
#[derive(Debug)]
pub enum MixedViscoSphereTrajectoryError {
    Invalid(&'static str),
    /// The native exact-volume BE chord is not a guaranteed tangent velocity;
    /// a pressure-constrained impulse requires a different compatible integrator.
    ExactIncompressibleChordVelocity,
    Solid(MixedViscoTetError),
    Surface(QuadraticSolidSurfaceError),
    Scene(QuadraticSphereSceneError),
    Volume(QuadraticSolidVolumeError),
    InitialVolumeUnresolved {
        report: QuadraticSolidVolumePathReport,
    },
    Impact(InelasticMixedContactError<MixedViscoTetError>),
    Rigid(FreeFlightError),
    Path(FiniteSolidSphereTrajectoryError),
    RejectionBudget {
        elapsed_s: f64,
        last: Box<MixedViscoSpherePathRefusal>,
    },
    Cancelled,
    Audit(Box<MixedViscoSphereTrajectoryReport>),
}
impl core::fmt::Display for MixedViscoSphereTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "native quadratic viscoelastic sphere trajectory: {self:?}"
        )
    }
}
impl std::error::Error for MixedViscoSphereTrajectoryError {}

/// Advance finite-compliance P2/P1 Maxwell material and COM-centred rigid spheres.
///
/// Supports hold their initially stationary displacement DOFs. Dead nodal and
/// reference-volume forces, rigid wrenches and rigid gravity stay fixed during
/// this call. Solid gravity must enter its original load projection separately.
/// The original pressure gauge, Newton controls, mass, material and rigid tensor
/// are passed unchanged to their native owners.
///
/// Positive activation localizes transient impacts. All resolved active rows
/// enter one simultaneous native impulse. Rejected paths shorten the interval
/// and recompute from the same accepted post-impact Maxwell state; no pressure,
/// position, velocity or history is interpolated or projected. A contact first
/// activating at the last endpoint resolves on the next continuation.
///
/// Each trial must prove positive Jacobians throughout every closed TET10
/// material domain and its complete affine time path before boundary admission.
/// This local regularity does not certify global absence of solid self-contact;
/// the original mesh must represent a nonoverlapping, consistently oriented body.
/// No sphere-pair law or rotating offset-centre path is inferred. Persistent
/// contact, unresolved grazing/Zeno motion and exhausted work allowances refuse.
/// Physical bulk/remainder storage, relaxation, support/source work and contact
/// loss use the same complete-call accounting as the original Tet4 driver;
/// numerical kinetic, angular and endpoint-work defects retain separate gates.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn mixed_visco_sphere_trajectory_candidate<'a>(
    solid: &MixedViscoTetMotion<'a>,
    surface: &QuadraticSolidContactSurface,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[FiniteSolidSphere],
    held_velocity_dofs: &[usize],
    nodal_forces_n: &[(usize, f64)],
    body_force_n_m3: [f64; 3],
    pressure_gauge: Option<(usize, f64)>,
    wrenches: &[Wrench],
    gravity: Gravity,
    duration_s: f64,
    origin_m: [f64; 3],
    controls: MixedViscoSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<MixedViscoSphereTrajectoryCandidate<'a>, MixedViscoSphereTrajectoryError> {
    advance_scene(
        solid,
        surface,
        bodies,
        properties,
        TrajectoryScene::Surface(spheres),
        held_velocity_dofs,
        nodal_forces_n,
        body_force_n_m3,
        pressure_gauge,
        wrenches,
        gravity,
        duration_s,
        origin_m,
        controls,
        cx,
    )
}

/// Advance a complete prepared sphere scene through one simultaneous native
/// impact system, including rigid pairs and curved quadratic solid contacts.
///
/// Every rigid pair keeps its explicitly declared law and canonical scene key.
/// Original native mass, pressure, histories, rigid inertia and all physical
/// and numerical gates are identical to the surface-only entry. Possible pair
/// encounters shorten the native interval before an active row can be resolved;
/// no pair is ignored or assigned a sphere's different solid-interface law.
/// This is a transient-impact driver; loaded dwell uses [`sustained`].
#[allow(clippy::too_many_arguments)]
pub fn mixed_visco_sphere_scene_trajectory_candidate<'a>(
    solid: &MixedViscoTetMotion<'a>,
    surface: &QuadraticSolidContactSurface,
    scene: &FiniteSolidSphereScene,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    held_velocity_dofs: &[usize],
    nodal_forces_n: &[(usize, f64)],
    body_force_n_m3: [f64; 3],
    pressure_gauge: Option<(usize, f64)>,
    wrenches: &[Wrench],
    gravity: Gravity,
    duration_s: f64,
    origin_m: [f64; 3],
    controls: MixedViscoSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<MixedViscoSphereTrajectoryCandidate<'a>, MixedViscoSphereTrajectoryError> {
    advance_scene(
        solid,
        surface,
        bodies,
        properties,
        TrajectoryScene::Complete(scene),
        held_velocity_dofs,
        nodal_forces_n,
        body_force_n_m3,
        pressure_gauge,
        wrenches,
        gravity,
        duration_s,
        origin_m,
        controls,
        cx,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn advance_scene<'a>(
    solid: &MixedViscoTetMotion<'a>,
    surface: &QuadraticSolidContactSurface,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    scene: TrajectoryScene<'_>,
    held_velocity_dofs: &[usize],
    nodal_forces_n: &[(usize, f64)],
    body_force_n_m3: [f64; 3],
    pressure_gauge: Option<(usize, f64)>,
    wrenches: &[Wrench],
    gravity: Gravity,
    duration_s: f64,
    origin_m: [f64; 3],
    controls: MixedViscoSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<MixedViscoSphereTrajectoryCandidate<'a>, MixedViscoSphereTrajectoryError> {
    use MixedViscoSphereTrajectoryError as E;
    poll(cx)?;
    let spheres = scene.spheres();
    validate(
        bodies, properties, spheres, wrenches, duration_s, origin_m, controls,
    )?;
    if solid.volumetric_treatment() == VolumetricTreatment::ExactIncompressible {
        return Err(E::ExactIncompressibleChordVelocity);
    }
    let ndofs = solid.velocity_mps().len() * 3;
    if held_velocity_dofs.len() > ndofs {
        return Err(E::Invalid("held velocities exceed the original solid DOFs"));
    }
    let mut held = held_velocity_dofs.to_vec();
    held.sort_unstable();
    if held.windows(2).any(|w| w[0] == w[1])
        || held
            .iter()
            .any(|&d| d >= ndofs || solid.velocity_mps()[d / 3][d % 3] != 0.0)
    {
        return Err(E::Invalid(
            "distinct valid stationary support DOFs required",
        ));
    }
    let prescribed: Vec<_> = held
        .iter()
        .map(|&d| (d, solid.solid().displacement_m()[d / 3][d % 3]))
        .collect();
    let mut spheres = spheres.to_vec();
    spheres.sort_unstable_by_key(|s| s.id);
    let start_time = solid.solid().time_s();
    let end_time = start_time + duration_s;
    if !end_time.is_finite() || end_time <= start_time {
        return Err(E::Path(FiniteSolidSphereTrajectoryError::NoProgress {
            elapsed_s: 0.0,
        }));
    }
    let prepared_volume = surface
        .prepare_volume_paths(controls.volume, cx)
        .map_err(E::Volume)?;
    let reference_volume = *prepared_volume.reference_admission();
    let mut coordinates = positions(surface, solid, cx)?;
    let initial_volume = prepared_volume
        .affine_path_at_positions(&coordinates, &coordinates, controls.volume, cx)
        .map_err(E::Volume)?;
    if matches!(
        initial_volume.verdict,
        QuadraticSolidVolumePathVerdict::Unresolved { .. }
    ) {
        return Err(E::InitialVolumeUnresolved {
            report: initial_volume,
        });
    }
    let mut current: Option<MixedViscoTetMotion<'a>> = None;
    let mut rigid = bodies.to_vec();
    let mut report = MixedViscoSphereTrajectoryReport {
        duration_s,
        ..Default::default()
    };
    let mut elapsed = 0.0;
    let clock_tolerance = 8.0 * f64::EPSILON * start_time.abs().max(end_time.abs()).max(duration_s);
    let mut geometry = scene.contacts(surface, &coordinates, &rigid, controls.path.initial, cx)?;
    while elapsed < duration_s {
        poll(cx)?;
        if report.steps.len() == controls.maximum_steps {
            return Err(E::Path(FiniteSolidSphereTrajectoryError::StepBudget {
                elapsed_s: elapsed,
            }));
        }
        let accepted = current.as_ref().unwrap_or(solid);
        let maximum_gap = geometry.maximum_contact_gap_m;
        let event = if geometry.contacts.is_empty() {
            None
        } else {
            Some(
                mixed_visco_mixed_contact_candidate(
                    accepted,
                    &held,
                    &rigid,
                    properties,
                    &geometry.contacts,
                    origin_m,
                    controls.impact,
                    cx,
                )
                .map_err(E::Impact)?,
            )
        };
        let post = event.as_ref().map_or(accepted, |e| &e.solid);
        let post_rigid = event
            .as_ref()
            .map_or(rigid.as_slice(), |e| e.rigid_bodies.as_slice());
        // Original impulses change velocity only. Reuse the accepted coordinates
        // and endpoint geometry instead of rebuilding their material traces.
        let start = &coordinates;
        let mut dt = controls.maximum_step_s.min(duration_s - elapsed);
        let (step, momentum, flights, rigid_next, mut path, end, geometry_next) = loop {
            poll(cx)?;
            if !dt.is_finite()
                || dt <= 0.0
                || elapsed + dt <= elapsed
                || post.solid().time_s() + dt <= post.solid().time_s()
            {
                return Err(E::Path(FiniteSolidSphereTrajectoryError::NoProgress {
                    elapsed_s: elapsed,
                }));
            }
            let (step, momentum) = post
                .advance_with_momentum(
                    dt,
                    &prescribed,
                    nodal_forces_n,
                    body_force_n_m3,
                    pressure_gauge,
                    controls.settings,
                    controls.solid,
                    origin_m,
                    cx,
                )
                .map_err(E::Solid)?;
            check_clock(
                dt,
                step.motion.solid().time_s() - post.solid().time_s(),
                clock_tolerance,
            )
            .map_err(E::Path)?;
            let end = positions(surface, &step.motion, cx)?;
            let volume = prepared_volume
                .affine_path_at_positions(start, &end, controls.volume, cx)
                .map_err(E::Volume)?;
            if matches!(
                volume.verdict,
                QuadraticSolidVolumePathVerdict::Unresolved { .. }
            ) {
                dt = shortened(
                    dt,
                    elapsed,
                    &mut report,
                    controls,
                    MixedViscoSpherePathRefusal::Volume(volume),
                )?;
                continue;
            }
            let flights: Vec<_> = post_rigid
                .iter()
                .zip(properties)
                .zip(wrenches)
                .map(|((&old, &property), &wrench)| {
                    rigid_free_flight(old, property, wrench, gravity, dt, vector(origin_m), cx)
                })
                .collect::<Result<_, _>>()
                .map_err(E::Rigid)?;
            let next: Vec<_> = flights.iter().map(|f| f.receipt.state_after).collect();
            let path = match paths::admit_path(
                scene,
                surface,
                volume,
                start,
                &end,
                post_rigid,
                &next,
                &spheres,
                controls.path,
                cx,
            )? {
                Ok(path) => path,
                Err(refusal) => {
                    dt = shortened(dt, elapsed, &mut report, controls, refusal)?;
                    continue;
                }
            };
            let geometry = match scene.contacts(surface, &end, &next, controls.path.initial, cx) {
                Ok(geometry) => geometry,
                Err(error) => match geometry::refinable_endpoint(error) {
                    Ok(refusal) => {
                        dt = shortened(dt, elapsed, &mut report, controls, refusal)?;
                        continue;
                    }
                    Err(error) => return Err(error),
                },
            };
            break (step, momentum, flights, next, path, end, geometry);
        };
        path.contact_keys = core::mem::take(&mut geometry.keys);
        let impact = event.map(|e| e.report);
        accumulate(
            &mut report,
            (&step.material_energy).into(),
            step.inertia,
            momentum,
            &flights,
            impact.as_ref(),
        );
        elapsed += dt;
        if elapsed > duration_s {
            return Err(E::Path(FiniteSolidSphereTrajectoryError::NoProgress {
                elapsed_s: elapsed,
            }));
        }
        check_clock(
            elapsed,
            step.motion.solid().time_s() - start_time,
            clock_tolerance,
        )
        .map_err(E::Path)?;
        report.steps.push(MixedViscoSphereTrajectoryStep {
            duration_s: dt,
            elapsed_s: elapsed,
            maximum_contact_gap_m: maximum_gap,
            impact,
            material: step.material_energy,
            inertia: step.inertia,
            momentum,
            rigid: flights.into_iter().map(|f| f.receipt).collect(),
            path,
        });
        current = Some(step.motion);
        coordinates = end;
        rigid = rigid_next;
        geometry = geometry_next;
        finish(&mut report);
        if !admitted(&report, controls.into()) {
            return Err(E::Audit(Box::new(report)));
        }
    }
    let solid = current.expect("positive interval accepted a native material step");
    check_clock(
        duration_s,
        solid.solid().time_s() - start_time,
        clock_tolerance,
    )
    .map_err(E::Path)?;
    poll(cx)?;
    Ok(MixedViscoSphereTrajectoryCandidate {
        solid,
        rigid_bodies: rigid,
        reference_volume,
        initial_volume,
        report,
    })
}

fn shortened(
    dt: f64,
    elapsed_s: f64,
    report: &mut MixedViscoSphereTrajectoryReport,
    controls: MixedViscoSphereTrajectoryControls,
    refusal: MixedViscoSpherePathRefusal,
) -> Result<f64, MixedViscoSphereTrajectoryError> {
    if report.rejected_steps == controls.maximum_rejections {
        return Err(MixedViscoSphereTrajectoryError::RejectionBudget {
            elapsed_s,
            last: Box::new(refusal),
        });
    }
    report.rejected_steps += 1;
    let fraction = match refusal {
        MixedViscoSpherePathRefusal::Volume(QuadraticSolidVolumePathReport {
            verdict: QuadraticSolidVolumePathVerdict::Unresolved { fraction_start, .. },
            ..
        })
        | MixedViscoSpherePathRefusal::Boundary(QuadraticSolidPathReport {
            verdict: QuadraticSolidPathVerdict::PossibleContact { fraction_start, .. },
            ..
        })
        | MixedViscoSpherePathRefusal::ScenePath(QuadraticSphereScenePathReport {
            verdict: QuadraticSphereScenePathVerdict::PossibleContact { fraction_start, .. },
            ..
        }) => fraction_start,
        _ => 0.0,
    };
    Ok(dt * if fraction > 0.0 { 0.9 * fraction } else { 0.5 })
}

impl From<MixedViscoSphereTrajectoryControls> for TrajectoryAuditControls {
    fn from(c: MixedViscoSphereTrajectoryControls) -> Self {
        Self {
            linear_tolerance_kg_m_s: c.linear_tolerance_kg_m_s,
            angular_tolerance_kg_m2_s: c.angular_tolerance_kg_m2_s,
            energy_tolerance_j: c.energy_tolerance_j,
            integration_tolerance_j: c.integration_tolerance_j,
            angular_integration_tolerance_kg_m2_s: c.angular_integration_tolerance_kg_m2_s,
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn validate(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[FiniteSolidSphere],
    wrenches: &[Wrench],
    duration: f64,
    origin: [f64; 3],
    c: MixedViscoSphereTrajectoryControls,
) -> Result<(), MixedViscoSphereTrajectoryError> {
    use MixedViscoSphereTrajectoryError as E;
    let surface = c.path.initial.surface;
    if bodies.len() != properties.len()
        || bodies.len() != wrenches.len()
        || bodies.len() != spheres.len()
        || bodies.len() > surface.maximum_spheres
        || c.impact.maximum_rigid_bodies == 0
        || bodies.len() > c.impact.maximum_rigid_bodies
        || !(1..=21).contains(&c.impact.contact.maximum_contacts)
        || !(1..=65_536).contains(&c.maximum_steps)
        || c.maximum_rejections > 65_536
        || origin.iter().any(|v| !v.is_finite())
        || !surface.activation_tolerance_m.is_finite()
        || surface.activation_tolerance_m <= 0.0
        || !c.impact.geometry_tolerance_m.is_finite()
        || surface.activation_tolerance_m > c.impact.geometry_tolerance_m
        || [
            duration,
            c.maximum_step_s,
            c.linear_tolerance_kg_m_s,
            c.angular_tolerance_kg_m2_s,
            c.energy_tolerance_j,
            c.integration_tolerance_j,
            c.angular_integration_tolerance_kg_m2_s,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(E::Invalid(
            "bounded matching sphere scene and positive physical/time gates required",
        ));
    }
    c.path
        .sweep
        .validate()
        .map_err(|e| E::Surface(QuadraticSolidSurfaceError::Query(e)))?;
    let mut occupied = vec![false; bodies.len()];
    for sphere in spheres {
        if sphere.rigid_body >= bodies.len()
            || occupied[sphere.rigid_body]
            || sphere.center_body_m != Vec3::ZERO
        {
            return Err(E::Invalid(
                "exactly one COM-centred sphere per rigid body required",
            ));
        }
        occupied[sphere.rigid_body] = true;
    }
    Ok(())
}

fn positions(
    surface: &QuadraticSolidContactSurface,
    solid: &MixedViscoTetMotion<'_>,
    cx: &Cx<'_>,
) -> Result<Vec<[f64; 3]>, MixedViscoSphereTrajectoryError> {
    let state = solid.solid();
    surface
        .current_positions_from_reference(
            state.reference_nodes_m(),
            state.reference_tetrahedra(),
            state.displacement_m(),
            cx,
        )
        .map_err(MixedViscoSphereTrajectoryError::Surface)
}

fn poll(cx: &Cx<'_>) -> Result<(), MixedViscoSphereTrajectoryError> {
    cx.checkpoint()
        .map_err(|_| MixedViscoSphereTrajectoryError::Cancelled)
}
