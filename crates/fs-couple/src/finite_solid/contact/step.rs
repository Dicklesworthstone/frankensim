//! Native rigid/deformable free motion and geometry-admitted impact trajectories.
//!
//! Automatic current-surface witnesses, the original simultaneous impulse solve,
//! and the original smooth mechanical integrators advance one private candidate.
//! A cached BVH and interval separating-plane bounds admit the ENTIRE declared
//! straight coordinate interpolation between native endpoints. An unresolved
//! window shortens the native step; every retry recomputes the native trajectory.
//! No endpoint is projected and no sampled gap is called a path certificate.
//!
//! The smooth-flight entry points refine transient impacts; the `sustained`
//! child advances loaded contact with an explicitly first-order force/contact/
//! kinetic-drift split. Both retain positive contact activation distance and
//! share the complete affine-path authority. Neither is a general Moreau--Jean
//! method, exact-flow CCD, or a guarantee against Zeno behavior. Unresolved
//! contact, geometry or time progress returns no candidate.

/// Original smooth mechanical owners and their complete momentum/work account.
pub mod free;
/// Original force/contact/kinetic-drift composition for sustained loaded contact.
pub mod sustained;

use super::mixed::FiniteSolidMixedContact;
use super::scene::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneError, FiniteSolidSphereSceneKey,
};
use super::surface::{
    FiniteSolidContactSurface, FiniteSolidSphere, FiniteSolidSurfaceControls,
    FiniteSolidSurfaceError,
};
use super::{
    FiniteSolidRigidContactControls, FiniteSolidRigidContactError, FiniteSolidRigidContactReport,
};
use crate::finite_solid::FiniteSolidBody;
use free::{FreeFlightControls, FreeFlightError, FreeFlightReport};
use fs_exec::Cx;
use fs_mbd::{Gravity, MassProperties, RigidBodyState, Vec3, Wrench};
use fs_query::TriangleQueryError;
use fs_query::triangle::sweep::{
    LinearSweepControls, LinearSweepVerdict, linear_sphere_triangle, linear_spheres,
    linear_tetrahedron_orientation,
};

/// Explicit time, geometric work, contact and complete-trajectory balance gates.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidSphereTrajectoryControls {
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    /// Aggregate interval windows per proposed path, including element quality.
    pub maximum_path_windows: usize,
    /// BVH nodes/triangle tests plus rigid-sphere pair tests per proposed path.
    pub maximum_path_visits: usize,
    pub maximum_path_candidates: usize,
    pub sweep: LinearSweepControls,
    pub surface: FiniteSolidSurfaceControls,
    pub impact: FiniteSolidRigidContactControls,
    pub free: FreeFlightControls,
    /// Gates on sums of phase residual magnitudes, not cancelling signed errors.
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
    /// Sum of absolute native smooth integration defects across the trajectory.
    pub integration_tolerance_j: f64,
}
impl Default for FiniteSolidSphereTrajectoryControls {
    fn default() -> Self {
        Self {
            maximum_step_s: 1e-3,
            maximum_steps: 4096,
            maximum_rejections: 4096,
            maximum_path_windows: 1 << 20,
            maximum_path_visits: 1 << 20,
            maximum_path_candidates: 1 << 16,
            sweep: LinearSweepControls::default(),
            surface: FiniteSolidSurfaceControls::default(),
            impact: FiniteSolidRigidContactControls::default(),
            free: FreeFlightControls::default(),
            linear_tolerance_kg_m_s: 1e-7,
            angular_tolerance_kg_m2_s: 1e-7,
            energy_tolerance_j: 1e-7,
            integration_tolerance_j: 1e-6,
        }
    }
}

/// Complete path admission for accepted native endpoints.
#[derive(Clone, Debug)]
pub struct FiniteSolidSpherePathReport {
    pub triangle_candidates: usize,
    pub query_visits: usize,
    pub interval_windows: usize,
    pub tetrahedra_checked: usize,
    pub sphere_pairs_checked: usize,
    /// Lower bound for interval-checked sphere/triangle and sphere-pair paths; broad-phase
    /// exclusions are clear. Infinity means no exact pair needed checking.
    pub minimum_clearance_m: f64,
}

/// One private accepted interval and its actual physical event, if any.
#[derive(Clone, Debug)]
pub struct FiniteSolidSphereTrajectoryStep {
    pub elapsed_s: f64,
    /// Largest actual activation gap; witnesses have never been moved together.
    pub maximum_contact_gap_m: f64,
    /// Original scene features in the same order as the submitted contact rows.
    pub contact_keys: Vec<FiniteSolidSphereSceneKey>,
    pub impact: Option<FiniteSolidRigidContactReport>,
    pub free: FreeFlightReport,
    pub path: FiniteSolidSpherePathReport,
}

/// One whole-call candidate, including all unallocated impact loss.
#[derive(Clone, Debug)]
pub struct FiniteSolidSphereTrajectoryReport {
    pub duration_s: f64,
    pub rejected_steps: usize,
    pub steps: Vec<FiniteSolidSphereTrajectoryStep>,
    pub kinetic_change_j: f64,
    pub strain_change_j: f64,
    pub source_work_j: f64,
    pub unallocated_contact_loss_j: f64,
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

#[derive(Clone, Copy, Debug)]
pub enum FiniteSolidSpherePathRefusal {
    /// Retry hint for a possible encounter; not a claimed first impact time.
    PossibleContact {
        sphere_id: u64,
        face: usize,
        fraction_start: f64,
        fraction_end: f64,
    },
    PossibleElementCollapse {
        tetrahedron: usize,
        fraction_start: f64,
        fraction_end: f64,
    },
    PossibleSpherePair {
        first: u64,
        second: u64,
        fraction_start: f64,
        fraction_end: f64,
    },
}
impl FiniteSolidSpherePathRefusal {
    fn start(self) -> f64 {
        match self {
            Self::PossibleContact { fraction_start, .. }
            | Self::PossibleElementCollapse { fraction_start, .. }
            | Self::PossibleSpherePair { fraction_start, .. } => fraction_start,
        }
    }
}

/// Every error leaves all caller-owned states and native histories unchanged.
#[derive(Debug)]
pub enum FiniteSolidSphereTrajectoryError {
    Invalid(&'static str),
    Surface(FiniteSolidSurfaceError),
    Scene(FiniteSolidSphereSceneError),
    Query(TriangleQueryError),
    /// The legacy surface-only entry point has no explicit sphere-pair laws.
    /// Use `mixed_sphere_trajectory_candidate` for a complete declared scene.
    UnsupportedSpherePair {
        first: u64,
        second: u64,
        fraction: f64,
    },
    Impact(FiniteSolidRigidContactError),
    Free(FreeFlightError),
    PathBudget,
    StepBudget {
        elapsed_s: f64,
    },
    RejectionBudget {
        elapsed_s: f64,
        last: FiniteSolidSpherePathRefusal,
    },
    NoProgress {
        elapsed_s: f64,
    },
    TimeResolution {
        requested_step_s: f64,
        represented_step_s: f64,
    },
    Cancelled,
    Audit(Box<FiniteSolidSphereTrajectoryReport>),
}
impl From<FiniteSolidSurfaceError> for FiniteSolidSphereTrajectoryError {
    fn from(e: FiniteSolidSurfaceError) -> Self {
        Self::Surface(e)
    }
}
impl From<FiniteSolidSphereSceneError> for FiniteSolidSphereTrajectoryError {
    fn from(e: FiniteSolidSphereSceneError) -> Self {
        Self::Scene(e)
    }
}
impl From<TriangleQueryError> for FiniteSolidSphereTrajectoryError {
    fn from(e: TriangleQueryError) -> Self {
        Self::Query(e)
    }
}
impl From<FiniteSolidRigidContactError> for FiniteSolidSphereTrajectoryError {
    fn from(e: FiniteSolidRigidContactError) -> Self {
        Self::Impact(e)
    }
}
impl From<FreeFlightError> for FiniteSolidSphereTrajectoryError {
    fn from(e: FreeFlightError) -> Self {
        Self::Free(e)
    }
}
impl core::fmt::Display for FiniteSolidSphereTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid/finite-solid sphere trajectory: {self:?}")
    }
}
impl std::error::Error for FiniteSolidSphereTrajectoryError {}

impl<'a> FiniteSolidBody<'a> {
    /// Advance automatically detected sphere/finite-solid transient impacts.
    ///
    /// Every rigid body has exactly one COM-centred sphere, with held declared
    /// friction and normal law. The rigid inertia remains its full native tensor.
    /// Offset sphere geometry is supported by discrete `sphere_contacts`, but
    /// rotation makes its centre path curved and is refused by this straight
    /// interpolation driver. This surface-only compatibility entry point has no
    /// sphere-pair laws and refuses rigid/rigid encounters. The mixed entry point
    /// below resolves both kinds together through the same trajectory owner.
    ///
    /// An impact acts at the actual current configuration within the positive
    /// activation distance. Native smooth motion then advances both owners.
    /// Original surface geometry and interval bounds admit every interpolated
    /// path, including triangles omitted by discrete proximity. Any possible
    /// crossing retries a shorter ORIGINAL mechanical step. Endpoint positions
    /// are never repaired or replaced with interpolated event samples.
    ///
    /// Straight interpolation defines the collision path of this numerical
    /// composition; it is not the true nonlinear flow between endpoints. Smooth
    /// native order does not establish a global impact-trajectory order. Spatial
    /// activation and time limits require refinement for the intended observable.
    /// Sustained loaded contact, unresolved grazing or accumulating impacts may
    /// exhaust progress/work bounds; no artificial resting/runout state is made.
    /// Pressure-bound solids retain the explicit-pressure native refusal.
    ///
    /// Gravity applies to rigid bodies. Supply the solid's resolved gravity in
    /// its held nodal/body loads. Total impact loss remains unallocated; neither
    /// it nor the reported integration defect is automatically converted to heat.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn sphere_trajectory_candidate(
        &self,
        surface: &FiniteSolidContactSurface,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        spheres: &[FiniteSolidSphere],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        wrenches: &[Wrench],
        gravity: Gravity,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FiniteSolidSphereTrajectoryControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, Vec<RigidBodyState>, FiniteSolidSphereTrajectoryReport),
        FiniteSolidSphereTrajectoryError,
    > {
        self.sphere_trajectory_with(
            surface,
            None,
            bodies,
            properties,
            spheres,
            nodal_forces_n,
            body_force_n_m3,
            wrenches,
            gravity,
            duration_s,
            origin_m,
            controls,
            cx,
        )
    }

    /// Advance automatically detected rigid/rigid AND rigid/solid impacts.
    ///
    /// Every active row enters one simultaneous solve, including cross-coupling
    /// through a body shared by both interfaces. The prepared scene owns explicit
    /// pair laws and stable feature identities. Original full rigid inertia and
    /// consistent solid mass remain in force. All admission, refinement, native
    /// free flight and whole-trajectory audits are shared with the surface lane.
    /// This transient lane retains that lane's declared straight-path/spatial
    /// activation scope; sustained loaded contact requires the force/contact
    /// split rather than repeatedly applying gravity after an impact.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_sphere_trajectory_candidate(
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
        (Self, Vec<RigidBodyState>, FiniteSolidSphereTrajectoryReport),
        FiniteSolidSphereTrajectoryError,
    > {
        self.sphere_trajectory_with(
            surface,
            Some(scene),
            bodies,
            properties,
            scene.spheres(),
            nodal_forces_n,
            body_force_n_m3,
            wrenches,
            gravity,
            duration_s,
            origin_m,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn sphere_trajectory_with(
        &self,
        surface: &FiniteSolidContactSurface,
        scene: Option<&FiniteSolidSphereScene>,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        spheres: &[FiniteSolidSphere],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        wrenches: &[Wrench],
        gravity: Gravity,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FiniteSolidSphereTrajectoryControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, Vec<RigidBodyState>, FiniteSolidSphereTrajectoryReport),
        FiniteSolidSphereTrajectoryError,
    > {
        poll(cx)?;
        validate(bodies, properties, spheres, wrenches, duration_s, controls)?;
        let start_time = self.native_state().t;
        let end_time = start_time + duration_s;
        if !end_time.is_finite() || end_time <= self.native_state().t {
            return Err(FiniteSolidSphereTrajectoryError::NoProgress { elapsed_s: 0.0 });
        }
        let mut sphere_order = spheres.to_vec();
        sphere_order.sort_unstable_by_key(|s| s.id);
        // Initial admission precedes any candidate mechanical work.
        current_contacts(
            surface,
            scene,
            self,
            bodies,
            &sphere_order,
            controls.surface,
            cx,
        )?;
        let mut solid = self.clone();
        let mut rigid = bodies.to_vec();
        let mut report = empty_report(duration_s);
        let mut elapsed = 0.0;
        let clock_tolerance = 8.0 * f64::EPSILON * duration_s;
        while elapsed < duration_s {
            poll(cx)?;
            if report.steps.len() == controls.maximum_steps {
                return Err(FiniteSolidSphereTrajectoryError::StepBudget { elapsed_s: elapsed });
            }
            let (contacts, contact_keys, maximum_gap) = current_contacts(
                surface,
                scene,
                &solid,
                &rigid,
                &sphere_order,
                controls.surface,
                cx,
            )?;
            let (post_solid, post_rigid, impact) = if contacts.is_empty() {
                (solid.clone(), rigid.clone(), None)
            } else {
                let (s, r, impact) = solid.mixed_contact_candidate(
                    &rigid,
                    properties,
                    &contacts,
                    origin_m,
                    controls.impact,
                    cx,
                )?;
                (s, r, Some(impact))
            };
            let mut dt = controls.maximum_step_s.min(duration_s - elapsed);
            let (next_solid, next_rigid, free, path) = loop {
                poll(cx)?;
                if !dt.is_finite()
                    || dt <= 0.0
                    || elapsed + dt <= elapsed
                    || solid.native_state().t + dt <= solid.native_state().t
                {
                    return Err(FiniteSolidSphereTrajectoryError::NoProgress {
                        elapsed_s: elapsed,
                    });
                }
                let (candidate, rigid_candidate, free) = post_solid.rigid_free_flight_candidate(
                    &post_rigid,
                    properties,
                    nodal_forces_n,
                    body_force_n_m3,
                    wrenches,
                    gravity,
                    dt,
                    origin_m,
                    controls.free,
                    cx,
                )?;
                let represented = candidate.native_state().t - post_solid.native_state().t;
                if (represented - dt).abs() > clock_tolerance {
                    return Err(FiniteSolidSphereTrajectoryError::TimeResolution {
                        requested_step_s: dt,
                        represented_step_s: represented,
                    });
                }
                match admit_path(
                    surface,
                    &post_solid,
                    &candidate,
                    &post_rigid,
                    &rigid_candidate,
                    &sphere_order,
                    scene.is_some(),
                    controls,
                    cx,
                )? {
                    Ok(path) => {
                        // Recheck actual reconstructed material witnesses, not
                        // just the coordinate interpolation's analytic support.
                        current_contacts(
                            surface,
                            scene,
                            &candidate,
                            &rigid_candidate,
                            &sphere_order,
                            controls.surface,
                            cx,
                        )?;
                        break (candidate, rigid_candidate, free, path);
                    }
                    Err(refusal) => {
                        if report.rejected_steps == controls.maximum_rejections {
                            return Err(FiniteSolidSphereTrajectoryError::RejectionBudget {
                                elapsed_s: elapsed,
                                last: refusal,
                            });
                        }
                        report.rejected_steps += 1;
                        let fraction = refusal.start();
                        dt *= if fraction > 0.0 { 0.9 * fraction } else { 0.5 };
                    }
                }
            };
            accumulate(&mut report, &free, impact.as_ref());
            let next_elapsed = elapsed + dt;
            if next_elapsed > duration_s || next_elapsed <= elapsed {
                return Err(FiniteSolidSphereTrajectoryError::NoProgress { elapsed_s: elapsed });
            }
            elapsed = next_elapsed;
            solid = next_solid;
            rigid = next_rigid;
            if (solid.native_state().t - (start_time + elapsed)).abs() > clock_tolerance {
                return Err(FiniteSolidSphereTrajectoryError::TimeResolution {
                    requested_step_s: elapsed,
                    represented_step_s: solid.native_state().t - start_time,
                });
            }
            report.steps.push(FiniteSolidSphereTrajectoryStep {
                elapsed_s: elapsed,
                maximum_contact_gap_m: maximum_gap,
                contact_keys,
                impact,
                free,
                path,
            });
        }
        if (solid.native_state().t - end_time).abs() > clock_tolerance {
            return Err(FiniteSolidSphereTrajectoryError::TimeResolution {
                requested_step_s: duration_s,
                represented_step_s: solid.native_state().t - start_time,
            });
        }
        finish(&mut report);
        let scalars = [
            report.kinetic_change_j,
            report.strain_change_j,
            report.source_work_j,
            report.unallocated_contact_loss_j,
            report.integration_defect_j,
            report.absolute_integration_defect_j,
            report.work_energy_residual_j,
            report.accumulated_linear_residual_kg_m_s,
            report.accumulated_angular_residual_kg_m2_s,
            report.accumulated_energy_residual_j,
        ];
        let vectors = [
            report.linear_momentum_change_kg_m_s,
            report.angular_momentum_change_kg_m2_s,
            report.external_impulse_kg_m_s,
            report.external_angular_impulse_kg_m2_s,
            report.support_impulse_kg_m_s,
            report.support_angular_impulse_kg_m2_s,
            report.linear_momentum_residual_kg_m_s,
            report.angular_momentum_residual_kg_m2_s,
        ];
        if scalars.into_iter().any(|x| !x.is_finite())
            || vectors.into_iter().any(|v| !v.is_finite())
            || report.accumulated_linear_residual_kg_m_s > controls.linear_tolerance_kg_m_s
            || report.accumulated_angular_residual_kg_m2_s > controls.angular_tolerance_kg_m2_s
            || report.accumulated_energy_residual_j > controls.energy_tolerance_j
            || norm(report.linear_momentum_residual_kg_m_s) > controls.linear_tolerance_kg_m_s
            || norm(report.angular_momentum_residual_kg_m2_s) > controls.angular_tolerance_kg_m2_s
            || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
            || report.absolute_integration_defect_j > controls.integration_tolerance_j
        {
            return Err(FiniteSolidSphereTrajectoryError::Audit(Box::new(report)));
        }
        poll(cx)?;
        Ok((solid, rigid, report))
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn current_contacts(
    surface: &FiniteSolidContactSurface,
    scene: Option<&FiniteSolidSphereScene>,
    solid: &FiniteSolidBody<'_>,
    bodies: &[RigidBodyState],
    spheres: &[FiniteSolidSphere],
    controls: FiniteSolidSurfaceControls,
    cx: &Cx<'_>,
) -> Result<
    (
        Vec<FiniteSolidMixedContact>,
        Vec<FiniteSolidSphereSceneKey>,
        f64,
    ),
    FiniteSolidSphereTrajectoryError,
> {
    if let Some(scene) = scene {
        let query = scene.contacts(surface, solid, bodies, controls, cx)?;
        return Ok((query.contacts, query.keys, query.maximum_contact_gap_m));
    }
    let query = surface.sphere_separations(solid, bodies, spheres, controls, cx)?;
    let contacts = query
        .admitted_contacts()?
        .into_iter()
        .map(Into::into)
        .collect();
    let keys = query
        .contacts
        .iter()
        .map(|w| FiniteSolidSphereSceneKey::Solid {
            sphere: w.sphere_id,
            contact_id: w.contact.id,
            face: w.face,
        })
        .collect();
    let gap = query.contacts.iter().map(|w| w.gap_m).fold(0.0, f64::max);
    Ok((contacts, keys, gap))
}

#[allow(clippy::too_many_arguments)]
fn admit_path(
    surface: &FiniteSolidContactSurface,
    before: &FiniteSolidBody<'_>,
    after: &FiniteSolidBody<'_>,
    rigid_before: &[RigidBodyState],
    rigid_after: &[RigidBodyState],
    spheres: &[FiniteSolidSphere],
    allow_sphere_pairs: bool,
    controls: FiniteSolidSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<
    Result<FiniteSolidSpherePathReport, FiniteSolidSpherePathRefusal>,
    FiniteSolidSphereTrajectoryError,
> {
    let start = surface.current_positions(before, cx)?;
    let end = surface.current_positions(after, cx)?;
    let mut report = FiniteSolidSpherePathReport {
        triangle_candidates: 0,
        query_visits: 0,
        interval_windows: 0,
        tetrahedra_checked: 0,
        sphere_pairs_checked: 0,
        minimum_clearance_m: f64::INFINITY,
    };
    for (i, tet) in before.problem().tetrahedra.iter().enumerate() {
        poll(cx)?;
        let sweep = linear_tetrahedron_orientation(
            tet.map(|n| start[n]),
            tet.map(|n| end[n]),
            surface.reference_orientations()[i],
            window_controls(controls, &report)?,
            cx,
        )?;
        report.interval_windows += sweep.examined_windows;
        report.tetrahedra_checked += 1;
        if let LinearSweepVerdict::PossibleContact {
            fraction_start,
            fraction_end,
        } = sweep.verdict
        {
            return Ok(Err(FiniteSolidSpherePathRefusal::PossibleElementCollapse {
                tetrahedron: i,
                fraction_start,
                fraction_end,
            }));
        }
    }
    // Do not silently omit dynamic-body pairs from a purported complete scene.
    for (i, a) in spheres.iter().enumerate() {
        for b in &spheres[i + 1..] {
            poll(cx)?;
            if report.query_visits == controls.maximum_path_visits {
                return Err(FiniteSolidSphereTrajectoryError::PathBudget);
            }
            report.query_visits += 1;
            report.sphere_pairs_checked += 1;
            let pair = linear_spheres(
                array(rigid_before[a.rigid_body].pose().position_world()),
                array(rigid_after[a.rigid_body].pose().position_world()),
                a.radius_m,
                array(rigid_before[b.rigid_body].pose().position_world()),
                array(rigid_after[b.rigid_body].pose().position_world()),
                b.radius_m,
                window_controls(controls, &report)?,
                cx,
            )?;
            report.interval_windows += pair.examined_windows;
            match pair.verdict {
                LinearSweepVerdict::PossibleContact {
                    fraction_start,
                    fraction_end,
                } => {
                    if allow_sphere_pairs {
                        return Ok(Err(FiniteSolidSpherePathRefusal::PossibleSpherePair {
                            first: a.id,
                            second: b.id,
                            fraction_start,
                            fraction_end,
                        }));
                    }
                    return Err(FiniteSolidSphereTrajectoryError::UnsupportedSpherePair {
                        first: a.id,
                        second: b.id,
                        fraction: fraction_start,
                    });
                }
                LinearSweepVerdict::Clear { minimum_clearance } => {
                    report.minimum_clearance_m = report.minimum_clearance_m.min(minimum_clearance);
                }
            }
        }
    }
    let swept = surface.triangle_surface().refit_swept(&start, &end, cx)?;
    for sphere in spheres {
        poll(cx)?;
        let c0 = array(rigid_before[sphere.rigid_body].pose().position_world());
        let c1 = array(rigid_after[sphere.rigid_body].pose().position_world());
        let candidates = swept.sphere_candidates(
            c0,
            c1,
            sphere.radius_m,
            controls.maximum_path_visits - report.query_visits,
            controls.maximum_path_candidates - report.triangle_candidates,
            cx,
        )?;
        report.query_visits += candidates.visited_nodes + candidates.tested_triangles;
        report.triangle_candidates += candidates.triangles.len();
        for face in candidates.triangles {
            let nodes = surface.faces()[face];
            let sweep = linear_sphere_triangle(
                c0,
                c1,
                sphere.radius_m,
                nodes.map(|n| start[n]),
                nodes.map(|n| end[n]),
                window_controls(controls, &report)?,
                cx,
            )?;
            report.interval_windows += sweep.examined_windows;
            match sweep.verdict {
                LinearSweepVerdict::Clear { minimum_clearance } => {
                    report.minimum_clearance_m = report.minimum_clearance_m.min(minimum_clearance);
                }
                LinearSweepVerdict::PossibleContact {
                    fraction_start,
                    fraction_end,
                } => {
                    return Ok(Err(FiniteSolidSpherePathRefusal::PossibleContact {
                        sphere_id: sphere.id,
                        face,
                        fraction_start,
                        fraction_end,
                    }));
                }
            }
        }
    }
    Ok(Ok(report))
}

fn window_controls(
    c: FiniteSolidSphereTrajectoryControls,
    r: &FiniteSolidSpherePathReport,
) -> Result<LinearSweepControls, FiniteSolidSphereTrajectoryError> {
    if r.interval_windows == c.maximum_path_windows {
        return Err(FiniteSolidSphereTrajectoryError::PathBudget);
    }
    Ok(LinearSweepControls {
        maximum_windows: c
            .sweep
            .maximum_windows
            .min(c.maximum_path_windows - r.interval_windows),
        fraction_tolerance: c.sweep.fraction_tolerance,
    })
}

fn validate(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[FiniteSolidSphere],
    wrenches: &[Wrench],
    duration: f64,
    c: FiniteSolidSphereTrajectoryControls,
) -> Result<(), FiniteSolidSphereTrajectoryError> {
    if bodies.len() != properties.len()
        || bodies.len() != spheres.len()
        || bodies.len() != wrenches.len()
        || bodies.len() > c.surface.maximum_spheres
        || bodies.len() > c.impact.maximum_rigid_bodies
        || bodies.len() > c.free.maximum_rigid_bodies
        || c.maximum_steps == 0
        || c.maximum_steps > 65_536
        || c.maximum_rejections > 65_536
        || c.maximum_path_windows == 0
        || c.maximum_path_windows > 1 << 24
        || c.maximum_path_visits == 0
        || c.maximum_path_visits > 1 << 26
        || c.maximum_path_candidates == 0
        || c.maximum_path_candidates > 1 << 24
        || c.sweep.maximum_windows == 0
        || c.sweep.maximum_windows > 1 << 20
        || !c.sweep.fraction_tolerance.is_finite()
        || c.sweep.fraction_tolerance <= 0.0
        || c.sweep.fraction_tolerance >= 1.0
        || c.surface.activation_tolerance_m <= 0.0
        || c.surface.activation_tolerance_m > c.impact.geometry_tolerance_m
        || [
            duration,
            c.maximum_step_s,
            c.linear_tolerance_kg_m_s,
            c.angular_tolerance_kg_m2_s,
            c.energy_tolerance_j,
            c.integration_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(FiniteSolidSphereTrajectoryError::Invalid(
            "matching bounded sphere scene and positive finite time/physical gates required",
        ));
    }
    let mut occupied = vec![false; bodies.len()];
    for sphere in spheres {
        if sphere.rigid_body >= bodies.len()
            || occupied[sphere.rigid_body]
            || sphere.center_body_m != Vec3::ZERO
        {
            return Err(FiniteSolidSphereTrajectoryError::Invalid(
                "exactly one COM-centred sphere per rigid body required",
            ));
        }
        occupied[sphere.rigid_body] = true;
    }
    Ok(())
}

fn empty_report(duration_s: f64) -> FiniteSolidSphereTrajectoryReport {
    FiniteSolidSphereTrajectoryReport {
        duration_s,
        rejected_steps: 0,
        steps: Vec::new(),
        kinetic_change_j: 0.0,
        strain_change_j: 0.0,
        source_work_j: 0.0,
        unallocated_contact_loss_j: 0.0,
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
    r: &mut FiniteSolidSphereTrajectoryReport,
    f: &FreeFlightReport,
    impact: Option<&FiniteSolidRigidContactReport>,
) {
    r.kinetic_change_j += f.kinetic_change_j;
    r.strain_change_j += f.strain_change_j;
    r.source_work_j += f.source_work_j;
    r.integration_defect_j += f.integration_defect_j;
    r.absolute_integration_defect_j += f.absolute_integration_defect_j;
    r.linear_momentum_change_kg_m_s = r
        .linear_momentum_change_kg_m_s
        .add(f.linear_momentum_change_kg_m_s);
    r.angular_momentum_change_kg_m2_s = r
        .angular_momentum_change_kg_m2_s
        .add(f.angular_momentum_change_kg_m2_s);
    r.external_impulse_kg_m_s = r.external_impulse_kg_m_s.add(f.external_impulse_kg_m_s);
    r.external_angular_impulse_kg_m2_s = r
        .external_angular_impulse_kg_m2_s
        .add(f.external_angular_impulse_kg_m2_s);
    r.support_impulse_kg_m_s = r.support_impulse_kg_m_s.add(f.support_impulse_kg_m_s);
    r.support_angular_impulse_kg_m2_s = r
        .support_angular_impulse_kg_m2_s
        .add(f.support_angular_impulse_kg_m2_s);
    r.accumulated_linear_residual_kg_m_s += f.absolute_linear_residual_kg_m_s;
    r.accumulated_angular_residual_kg_m2_s += f.absolute_angular_residual_kg_m2_s;
    r.accumulated_energy_residual_j += f.absolute_work_residual_j;
    if let Some(i) = impact {
        let support = vector(i.solid_impulse.support.linear_kg_m_s);
        let moment = vector(i.solid_impulse.support.angular_kg_m2_s);
        r.kinetic_change_j += i.kinetic_change_j;
        r.unallocated_contact_loss_j += i.contact.dissipated_energy;
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
        let (linear, angular, energy) = impact_residuals(i);
        r.accumulated_linear_residual_kg_m_s += linear;
        r.accumulated_angular_residual_kg_m2_s += angular;
        r.accumulated_energy_residual_j += energy;
    }
}

/// Include both native representability errors and the complete event's
/// geometry/law residuals. In particular, real finite-gap tangential couples
/// must not cancel between accepted events or hide behind exact native updates.
fn impact_residuals(i: &FiniteSolidRigidContactReport) -> (f64, f64, f64) {
    (
        i.absolute_linear_residual_kg_m_s + norm(i.linear_momentum_residual_kg_m_s),
        i.absolute_angular_residual_kg_m2_s + norm(i.angular_momentum_residual_kg_m2_s),
        i.absolute_work_residual_j
            + i.work_energy_residual_j.abs()
            + i.dissipation_residual_j.abs()
            + (i.relative_impulse_work_j - i.contact.relative_work).abs(),
    )
}

fn finish(r: &mut FiniteSolidSphereTrajectoryReport) {
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
fn vector(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn poll(cx: &Cx<'_>) -> Result<(), FiniteSolidSphereTrajectoryError> {
    cx.checkpoint()
        .map_err(|_| FiniteSolidSphereTrajectoryError::Cancelled)
}
