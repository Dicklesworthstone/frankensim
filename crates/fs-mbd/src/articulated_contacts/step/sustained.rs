//! Sustained plane forces coupled to original free/prescribed articulated motion.
//!
//! The contact owner solves force/displacement complementarity. Every response
//! restarts from the same borrowed native state, with gravity, actuator effort,
//! body loads and held contact forces present together at both original stages.
//! There is no force kick, restitution solve, proxy inertia or contact spring.

use super::mixed::{native_input, totals};
use super::{
    ArticulatedContactStepError, com_displacement, finite_vector, impulse_pairing,
    relative_kinematics,
};
use crate::articulated::contact_geometry::{ContactGeometryError, prepare_contact_frame};
use crate::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBaseMotionStep,
    PrescribedBasePath,
};
use crate::articulated::motion::spatial_loads::prescribed::step_prescribed_base_spatial_with_effort;
use crate::articulated::motion::spatial_loads::{
    MaterialPoint, SpatialMotionLoads, WorldPointLoad, step_free_floating_spatial_with_effort,
};
use crate::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionStep, MAX_FREE_FLOATING_STEP_LINKS,
};
use crate::articulated::{ArticulatedError, BaseState, Kinematics, forward_kinematics};
use crate::articulated_contacts::mixed::ContactSystem;
use crate::articulated_contacts::trajectory::mixed::{
    MixedArticulatedTrajectoryState as State, MixedArticulatedTrajectorySystem as System,
};
use crate::articulated_contacts::trajectory::{ArticulatedTrajectoryMomentum, Impulses};
use crate::articulated_contacts::{ContactCounterpart, LinkContactPoint};
use fs_contact::sustained::{
    NaturalMapScales,
    solve::{self, SolveControls, SolveError, SolveReport},
};
use fs_ga::Vec3;
use fs_math::dd::Dd;

/// Original material point against a uniformly translating, fixed-normal plane.
#[derive(Clone, Copy, Debug)]
pub struct ArticulatedSustainedPlaneContact {
    /// Unique stable identity; force and path reports use ascending ID order.
    pub id: u64,
    /// Original material point receiving the force.
    pub point: LinkContactPoint,
    /// Any point on the plane at the common initial clock, in metres.
    pub plane_point_world: Vec3,
    /// Constant plane MATERIAL velocity in world axes, m/s, including belt slip.
    pub plane_velocity_world: Vec3,
    /// Separating normal and tangent hint, normalized by the native geometry owner.
    pub normal_world: Vec3,
    pub tangent_world: Vec3,
    /// Finite nonnegative isotropic Coulomb coefficient.
    pub friction_coefficient: f64,
    /// Numerical conditioning in N/m; these scales add no physical compliance.
    pub scales: NaturalMapScales,
}

/// Separate numerical, geometric and physical allowances in SI units.
#[derive(Clone, Copy, Debug)]
pub struct ArticulatedSustainedControls {
    pub duration: f64,
    pub angular_origin_world: Vec3,
    /// Original bounded Newton--Krylov controls, with a fresh native replay.
    pub solve: SolveControls,
    /// Maximum separation of a loaded witness; negative gaps always refuse.
    pub loaded_gap_tolerance_m: f64,
    /// Cone and unloaded-force resolution, N. Normal attraction always refuses.
    pub force_tolerance_n: f64,
    /// Normal dwell-speed resolution, m/s. Closing impacts need an impact owner.
    pub normal_velocity_tolerance_m_s: f64,
    /// Maximum normal complementarity and tangent duality-gap residual, J.
    pub complementarity_tolerance_j: f64,
    /// Actual kinetic change minus independently evaluated work, J.
    pub energy_tolerance_j: f64,
    /// Actual native momentum changes minus independent native source impulses.
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
}

/// Full native motion evidence for one retained tree.
#[derive(Clone, Debug)]
pub enum ArticulatedSustainedMotion {
    Free(FreeFloatingMotionStep),
    Prescribed(PrescribedBaseMotionStep),
}

/// Original candidate path and held law supplied to mandatory physical admission.
///
/// The caller must establish complete-path separation, omitted-feature coverage,
/// applicability of these fixed material points/normals, and continuous contact
/// mode: loaded points stay in contact; sticking has zero slip; sliding retains
/// the opposing Coulomb direction without an interior reversal. Stage samples
/// alone do not establish these claims. The caller also owns any additional
/// branch/determinacy evidence needed to adopt this local numerical candidate.
/// Callback-owned history must remain staged until the complete call succeeds.
pub struct ArticulatedSustainedPath<'a> {
    pub systems: &'a [System<'a>],
    /// Native midpoint PREDICTORS used for forces/work, not inferred interpolation.
    pub midpoint: &'a [State],
    pub endpoint: &'a [State],
    pub support_paths: &'a [Option<PrescribedBasePath>],
    /// Contacts and forces in the same canonical ID order.
    pub contacts: &'a [ArticulatedSustainedPlaneContact],
    pub axes_world: &'a [[Vec3; 3]],
    pub forces_n: &'a [[f64; 3]],
    pub time_before_s: f64,
    pub duration: f64,
}

/// Physically audited LOCAL candidate, without an existence/uniqueness claim.
///
/// Neither Newton convergence nor the physical residuals prove a unique force,
/// motion or contact mode. The supplied path admission remains part of this
/// conditional result. No Painlevé selection, collision discovery or generic
/// Moreau--Jean method is inferred from this sustained-force entry point.
pub struct ArticulatedSustainedStepCandidate<P> {
    pub time_before_s: f64,
    pub time_after_s: f64,
    pub states_after: Vec<State>,
    pub motion: Vec<ArticulatedSustainedMotion>,
    pub support_paths: Vec<Option<PrescribedBasePath>>,
    pub contact_ids: Vec<u64>,
    /// Unmodified final force solve and freshly recomputed endpoint gap/slip, SI units.
    pub solve: SolveReport,
    /// Total original response calls, including an optional inactive-mode restart.
    pub total_response_evaluations: usize,
    /// Total original outer attempts across both solves, within the caller cap.
    pub total_nonlinear_iterations: usize,
    pub path_admission: P,
    pub momentum_before: ArticulatedTrajectoryMomentum,
    pub momentum_after: ArticulatedTrajectoryMomentum,
    /// Actual per-native-stage increments, accumulated before rounding.
    pub kinetic_energy_change_j: f64,
    pub linear_momentum_change_world: Vec3,
    pub angular_momentum_change_about_origin_world: Vec3,
    /// Includes gravity, body loads, plane forces and prescribed support reaction.
    pub external_linear_impulse_world: Vec3,
    pub external_angular_impulse_about_origin_world: Vec3,
    /// Held effort and body-wrench work at the actual native midpoint stage, J.
    pub actuator_work_j: f64,
    pub body_wrench_work_j: f64,
    /// Uniform-gravity work from actual original COM displacements, J.
    pub gravity_work_j: f64,
    /// Original independently recovered prescribed support work, J.
    pub support_work_j: f64,
    /// Signed work supplied by translating plane material, J.
    pub plane_work_j: f64,
    /// Signed normal force times actual relative normal displacement, J.
    /// Finite configuration error is visible here; it is never frictional heat.
    pub normal_contact_work_j: f64,
    /// Nonnegative -force_t dot actual relative slip, J; no thermal assignment.
    pub frictional_dissipation_j: f64,
    /// Native midpoint contact work minus force times actual point displacement.
    /// This is a work-quadrature defect, not physical dissipation, J.
    pub contact_work_quadrature_defect_j: f64,
    /// Original native integration work defect; no contact loss is added to it.
    pub native_kinetic_work_residual_j: f64,
    /// Actual kinetic increment minus held/support/plane/normal work plus friction loss.
    pub kinetic_work_residual_j: f64,
    pub linear_momentum_residual_world: Vec3,
    pub angular_momentum_residual_about_origin_world: Vec3,
}

/// No native state publishes on any source, law, geometry, budget or balance refusal.
#[derive(Debug)]
pub enum ArticulatedSustainedError<E> {
    Invalid(&'static str),
    Geometry(ContactGeometryError),
    Dynamics {
        system_index: usize,
        error: ArticulatedError,
    },
    FreeMotion {
        system_index: usize,
        error: FreeFloatingMotionError,
    },
    PrescribedMotion {
        system_index: usize,
        error: PrescribedBaseMotionError,
    },
    /// Original local solver refusal, preserving its first mechanical source error.
    Solve(Box<SolveError<ArticulatedSustainedError<core::convert::Infallible>>>),
    /// Initial touching/near-contact motion requires a separate impact event.
    ImpactRequired {
        contact_id: u64,
        relative_normal_speed_m_s: f64,
    },
    Penetration {
        contact_id: u64,
        gap_m: f64,
    },
    ContactLaw {
        contact_id: u64,
        reason: &'static str,
    },
    Balance(&'static str),
    Path(E),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for ArticulatedSustainedError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated sustained contact: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ArticulatedSustainedError<E> {}

/// Solve held sustained-plane forces against original articulated trajectories.
///
/// All forces act simultaneously at both native stages, preserving full ABA
/// coupling, finite inertia, base recoil and driven support reactions. A resting
/// gravity equilibrium therefore needs no per-step velocity reset or artificial
/// kick loss. Tangent slip is actual material displacement relative to plane
/// material, including a translating belt's signed work.
///
/// Supply one initial force triple per contact in input order. Support samples
/// are requested once per prescribed tree and reused unchanged in every native
/// response. The original contact owner bounds all trial/replay evaluations;
/// total tree size is at most 4096 links in 128 systems and contact count at most
/// 21. Native mechanical passes remain O(links); no dense tree mass matrix is
/// introduced. Original stage/joint limits and source refusals are retained.
///
/// Trials may have negative gaps/normal forces because the nonlinear law must
/// see those residuals. Only the final replay undergoes physical adoption gates:
/// no penetration or attractive normal force, native force/cone/work/momentum
/// admission, and the required complete-path/mode callback. No state, geometry
/// or force is projected to make a failed candidate pass.
///
/// A clearly separated contact with a tiny attractive numerical force may seed
/// one fresh solve with that contact's exact inactive force triple. The original
/// native response and contact law must independently accept the new candidate;
/// both solves share the caller's response and outer-iteration budgets.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn step_sustained_articulated_planes<P, E>(
    systems: &[System<'_>],
    contacts: &[ArticulatedSustainedPlaneContact],
    initial_forces_n: &[[f64; 3]],
    start_time_s: f64,
    controls: ArticulatedSustainedControls,
    mut support_path: impl FnMut(
        usize,
        &PrescribedBaseMotionState,
        f64,
    ) -> Result<PrescribedBasePath, PrescribedBaseMotionError>,
    admit_path: impl FnOnce(&ArticulatedSustainedPath<'_>) -> Result<P, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedSustainedStepCandidate<P>, ArticulatedSustainedError<E>> {
    poll(&mut cancelled)?;
    let end_time = admit(systems, contacts, initial_forces_n, start_time_s, controls)?;
    let initial: Vec<_> = systems.iter().map(|s| s.state.clone()).collect();
    let before = frames(systems, &initial, &mut cancelled).map_err(lift)?;
    let borrowed: Vec<_> = before.iter().map(|f| &f.world).collect();
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    if order
        .windows(2)
        .any(|w| contacts[w[0]].id == contacts[w[1]].id)
    {
        return Err(ArticulatedSustainedError::Invalid(
            "duplicate contact identity",
        ));
    }
    let sorted: Vec<_> = order.iter().map(|&i| contacts[i]).collect();
    let mut prepared = Vec::with_capacity(contacts.len());
    for (index, contact) in sorted.iter().enumerate() {
        poll(&mut cancelled)?;
        let geometry = prepare_contact_frame(
            contact.point,
            ContactCounterpart::PrescribedPlane {
                point_world: contact.plane_point_world,
                velocity_world: contact.plane_velocity_world,
            },
            contact.normal_world,
            contact.tangent_world,
            index,
            &borrowed,
        )
        .map_err(ArticulatedSustainedError::Geometry)?;
        let axes = geometry.rows.map(|row| row.axis_world);
        let parts = point_parts(&before[contact.point.system_index], contact.point)
            .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
        let gap = point_gap(parts, contact.plane_point_world, axes[0]);
        let speed = relative_speed(&before[contact.point.system_index], *contact, axes)
            .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
        check_gap(contact.id, gap.to_f64())?;
        if gap.to_f64() <= controls.loaded_gap_tolerance_m
            && speed[0] < -controls.normal_velocity_tolerance_m_s
        {
            return Err(ArticulatedSustainedError::ImpactRequired {
                contact_id: contact.id,
                relative_normal_speed_m_s: speed[0],
            });
        }
        prepared.push(Prepared {
            contact: *contact,
            axes,
            initial_gap: gap,
            initial_parts: parts,
            initial_speed: speed,
        });
    }
    let mut support_paths = Vec::with_capacity(systems.len());
    for (index, system) in systems.iter().enumerate() {
        poll(&mut cancelled)?;
        support_paths.push(match system.state {
            State::Free(_) => None,
            State::Prescribed(state) => Some(
                support_path(index, state, controls.duration)
                    .map_err(|error| prescribed_error(index, error))?,
            ),
        });
    }
    let friction: Vec<_> = sorted.iter().map(|c| c.friction_coefficient).collect();
    let scales: Vec<_> = sorted.iter().map(|c| c.scales).collect();
    let mut warm: Vec<_> = order.iter().map(|&i| initial_forces_n[i]).collect();
    let mut solve_controls = controls.solve;
    let mut total_response_evaluations = 0;
    let mut total_nonlinear_iterations = 0;
    let mut restarted = false;
    let (report, forces, native) = loop {
        let mut replay = None;
        let report = solve::solve(
            &friction,
            &scales,
            &warm,
            solve_controls,
            &mut |forces, stop| {
                let native = response(
                    systems,
                    &before,
                    &prepared,
                    &support_paths,
                    forces,
                    controls,
                    stop,
                )?;
                let paths = native.paths.clone();
                replay = Some((forces.to_vec(), native));
                Ok(paths)
            },
            &mut cancelled,
        )
        .map_err(|error| {
            let error = match error {
                SolveError::ResponseBudget { evaluations, .. } => SolveError::ResponseBudget {
                    evaluations: total_response_evaluations + evaluations,
                    maximum: controls.solve.maximum_response_evaluations,
                },
                error => error,
            };
            ArticulatedSustainedError::Solve(Box::new(error))
        })?;
        total_response_evaluations += report.response_evaluations;
        total_nonlinear_iterations += report.nonlinear.iterations;
        let Some((forces, native)) = replay else {
            return Err(ArticulatedSustainedError::Invalid(
                "missing native final replay",
            ));
        };
        if forces != report.forces_n || native.paths != report.paths_m {
            return Err(ArticulatedSustainedError::Invalid(
                "native replay differs from the original solver candidate",
            ));
        }

        // This changes only a warm start. Nothing from this discarded response
        // is adopted, and every force/path/state is recomputed by the same owners.
        let mut inactive = Vec::new();
        if !restarted {
            for (index, (force, path)) in forces.iter().zip(&native.paths).enumerate() {
                if force[0] < 0.0
                    && force.iter().all(|f| f.abs() <= controls.force_tolerance_n)
                    && path[0] > controls.loaded_gap_tolerance_m
                {
                    inactive.push(index);
                }
            }
        }
        if inactive.is_empty() {
            break (report, forces, native);
        }
        if total_nonlinear_iterations >= controls.solve.maximum_iterations {
            return Err(ArticulatedSustainedError::ContactLaw {
                contact_id: sorted[inactive[0]].id,
                reason: "attractive normal force at exhausted outer iteration budget",
            });
        }
        if total_response_evaluations >= controls.solve.maximum_response_evaluations {
            return Err(ArticulatedSustainedError::Solve(Box::new(
                SolveError::ResponseBudget {
                    evaluations: total_response_evaluations,
                    maximum: controls.solve.maximum_response_evaluations,
                },
            )));
        }
        warm = forces;
        for index in inactive {
            warm[index] = [0.0; 3];
        }
        solve_controls.maximum_response_evaluations =
            controls.solve.maximum_response_evaluations - total_response_evaluations;
        solve_controls.maximum_iterations =
            controls.solve.maximum_iterations - total_nonlinear_iterations;
        restarted = true;
    };
    let mut plane_work = Dd::ZERO;
    let mut normal_work = Dd::ZERO;
    let mut tangential_work = Dd::ZERO;
    for (index, prepared) in prepared.iter().enumerate() {
        poll(&mut cancelled)?;
        let force = forces[index];
        let path = native.paths[index];
        check_contact(prepared, force, path, native.midpoint_gaps[index], controls)?;
        let relative = native.relative_displacement[index];
        normal_work = normal_work + Dd::from_f64(force[0]) * relative[0];
        let tangent = Dd::from_f64(force[1]) * relative[1] + Dd::from_f64(force[2]) * relative[2];
        if tangent.to_f64() > 0.0 {
            return Err(ArticulatedSustainedError::ContactLaw {
                contact_id: prepared.contact.id,
                reason: "positive relative tangential work",
            });
        }
        tangential_work = tangential_work + tangent;
        for (&force, axis) in force.iter().zip(prepared.axes) {
            plane_work = plane_work
                + Dd::from_f64(force)
                    * dot(axis, prepared.contact.plane_velocity_world)
                    * Dd::from_f64(controls.duration);
        }
    }
    let a = &native.accounts;
    let exact_contact_work = plane_work + normal_work + tangential_work;
    let work_residual =
        (a.kinetic - a.actuator - a.body - a.gravity - a.support - exact_contact_work).to_f64();
    let quadrature_defect = (a.contact_quadrature - exact_contact_work).to_f64();
    let (linear_residual, angular_residual) = a.actual.difference(&a.sources);
    let values = [
        a.kinetic,
        a.actuator,
        a.body,
        a.gravity,
        a.support,
        plane_work,
        normal_work,
        -tangential_work,
        a.native_defect,
    ]
    .map(Dd::to_f64);
    if values.iter().any(|v| !v.is_finite())
        || !work_residual.is_finite()
        || !quadrature_defect.is_finite()
        || ![
            a.actual.linear,
            a.actual.angular,
            a.sources.linear,
            a.sources.angular,
            linear_residual,
            angular_residual,
        ]
        .into_iter()
        .all(finite_vector)
    {
        return Err(ArticulatedSustainedError::Invalid(
            "nonfinite native sustained balance",
        ));
    }
    if work_residual.abs() > controls.energy_tolerance_j
        || linear_residual.norm() > controls.linear_momentum_tolerance_ns
        || angular_residual.norm() > controls.angular_momentum_tolerance_nms
    {
        return Err(ArticulatedSustainedError::Balance(
            "actual native work or momentum exceeds caller gate",
        ));
    }
    let momentum_before = totals(
        systems,
        &initial,
        controls.angular_origin_world,
        &mut cancelled,
    )
    .map_err(step_error)?;
    let momentum_after = totals(
        systems,
        &native.states,
        controls.angular_origin_world,
        &mut cancelled,
    )
    .map_err(step_error)?;
    let axes: Vec<_> = prepared.iter().map(|p| p.axes).collect();
    poll(&mut cancelled)?;
    let path_admission = admit_path(&ArticulatedSustainedPath {
        systems,
        midpoint: &native.midpoints,
        endpoint: &native.states,
        support_paths: &support_paths,
        contacts: &sorted,
        axes_world: &axes,
        forces_n: &forces,
        time_before_s: start_time_s,
        duration: controls.duration,
    })
    .map_err(ArticulatedSustainedError::Path)?;
    poll(&mut cancelled)?;
    Ok(ArticulatedSustainedStepCandidate {
        time_before_s: start_time_s,
        time_after_s: end_time,
        states_after: native.states,
        motion: native.motion,
        support_paths,
        contact_ids: sorted.iter().map(|c| c.id).collect(),
        solve: report,
        total_response_evaluations,
        total_nonlinear_iterations,
        path_admission,
        momentum_before,
        momentum_after,
        kinetic_energy_change_j: values[0],
        linear_momentum_change_world: a.actual.linear,
        angular_momentum_change_about_origin_world: a.actual.angular,
        external_linear_impulse_world: a.sources.linear,
        external_angular_impulse_about_origin_world: a.sources.angular,
        actuator_work_j: values[1],
        body_wrench_work_j: values[2],
        gravity_work_j: values[3],
        support_work_j: values[4],
        plane_work_j: values[5],
        normal_contact_work_j: values[6],
        frictional_dissipation_j: values[7],
        contact_work_quadrature_defect_j: quadrature_defect,
        native_kinetic_work_residual_j: values[8],
        kinetic_work_residual_j: work_residual,
        linear_momentum_residual_world: linear_residual,
        angular_momentum_residual_about_origin_world: angular_residual,
    })
}

struct Frame {
    state: State,
    base: BaseState,
    world: Kinematics,
    shape: Kinematics,
}
struct Prepared {
    contact: ArticulatedSustainedPlaneContact,
    axes: [Vec3; 3],
    initial_gap: Dd,
    initial_parts: [Vec3; 3],
    initial_speed: [f64; 3],
}
#[derive(Default)]
struct Accounts {
    actual: Impulses,
    sources: Impulses,
    kinetic: Dd,
    actuator: Dd,
    body: Dd,
    gravity: Dd,
    support: Dd,
    contact_quadrature: Dd,
    native_defect: Dd,
}
struct Response {
    states: Vec<State>,
    midpoints: Vec<State>,
    motion: Vec<ArticulatedSustainedMotion>,
    paths: Vec<[f64; 3]>,
    relative_displacement: Vec<[Dd; 3]>,
    midpoint_gaps: Vec<f64>,
    accounts: Accounts,
}
type NativeError = ArticulatedSustainedError<core::convert::Infallible>;

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn response(
    systems: &[System<'_>],
    before: &[Frame],
    prepared: &[Prepared],
    supports: &[Option<PrescribedBasePath>],
    forces: &[[f64; 3]],
    controls: ArticulatedSustainedControls,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Response, NativeError> {
    let mut points = vec![Vec::new(); systems.len()];
    for (p, force) in prepared.iter().zip(forces) {
        poll(cancelled)?;
        let coordinates: [Dd; 3] = core::array::from_fn(|axis| {
            p.axes.iter().zip(force).fold(Dd::ZERO, |s, (a, f)| {
                s + Dd::from_f64(components(*a)[axis]) * Dd::from_f64(*f)
            })
        });
        let force_world = vector(coordinates);
        if !finite_vector(force_world)
            || (force.iter().any(|&f| f != 0.0) && force_world == Vec3::default())
        {
            return Err(ArticulatedSustainedError::Invalid(
                "unrepresentable world contact force",
            ));
        }
        points[p.contact.point.system_index].push(WorldPointLoad {
            point: MaterialPoint {
                link: p.contact.point.link_index,
                point_body: p.contact.point.point_body,
            },
            force_world,
            couple_world: Vec3::default(),
        });
    }
    let mut states = Vec::with_capacity(systems.len());
    let mut midpoints = Vec::with_capacity(systems.len());
    let mut motion = Vec::with_capacity(systems.len());
    let mut accounts = Accounts::default();
    for (index, system) in systems.iter().enumerate() {
        poll(cancelled)?;
        let loads = SpatialMotionLoads {
            held: system.loads,
            world_points: &points[index],
            connections: &[],
        };
        let mut midpoint = None;
        let record = match system.state {
            State::Free(state) => {
                let step = step_free_floating_spatial_with_effort(
                    system.model,
                    state,
                    loads,
                    controls.duration,
                    |time, state, _efforts| {
                        if time > 0.0 {
                            midpoint = Some(State::Free(state.clone()));
                        }
                        Ok(())
                    },
                    &mut *cancelled,
                )
                .map_err(|error| free_error(index, error))?
                .motion;
                accounts.actual.add(
                    step.linear_momentum_change_world,
                    step.angular_momentum_change_world,
                    controls.angular_origin_world,
                );
                accounts.sources.add(
                    step.external_linear_impulse_world,
                    step.external_angular_impulse_world,
                    controls.angular_origin_world,
                );
                accounts.kinetic = accounts.kinetic + Dd::from_f64(step.kinetic_energy_change);
                accounts.native_defect =
                    accounts.native_defect + Dd::from_f64(step.kinetic_work_residual);
                states.push(State::Free(step.state_after.clone()));
                ArticulatedSustainedMotion::Free(step)
            }
            State::Prescribed(state) => {
                let step = step_prescribed_base_spatial_with_effort(
                    system.model,
                    state,
                    supports[index].ok_or(ArticulatedSustainedError::Invalid(
                        "missing prescribed path",
                    ))?,
                    loads,
                    controls.duration,
                    |time, state, _efforts| {
                        if time > 0.0 {
                            midpoint = Some(State::Prescribed(state.clone()));
                        }
                        Ok(())
                    },
                    &mut *cancelled,
                )
                .map_err(|error| prescribed_error(index, error))?
                .motion;
                if step.state_after.time_s != state.time_s + controls.duration {
                    return Err(ArticulatedSustainedError::Invalid(
                        "prescribed response changed clock",
                    ));
                }
                accounts.actual.add(
                    step.linear_momentum_change_world,
                    step.angular_momentum_change_world,
                    controls.angular_origin_world,
                );
                accounts.sources.add(
                    step.external_linear_impulse_world,
                    step.external_angular_impulse_world,
                    controls.angular_origin_world,
                );
                accounts.sources.add(
                    step.support_linear_impulse_world,
                    step.support_angular_impulse_world,
                    controls.angular_origin_world,
                );
                accounts.kinetic = accounts.kinetic + Dd::from_f64(step.kinetic_energy_change);
                accounts.support = accounts.support + Dd::from_f64(step.support_work);
                accounts.native_defect =
                    accounts.native_defect + Dd::from_f64(step.kinetic_work_residual);
                states.push(State::Prescribed(step.state_after.clone()));
                ArticulatedSustainedMotion::Prescribed(step)
            }
        };
        midpoints.push(midpoint.ok_or(ArticulatedSustainedError::Invalid(
            "native owner omitted midpoint stage",
        ))?);
        motion.push(record);
    }
    let middle = frames(systems, &midpoints, cancelled)?;
    let after = frames(systems, &states, cancelled)?;
    for (index, system) in systems.iter().enumerate() {
        let input = native_input(system.model, &middle[index].state);
        for (&effort, &velocity) in system
            .loads
            .generalized_force
            .iter()
            .zip(input.joint_velocity)
        {
            accounts.actuator = accounts.actuator
                + Dd::from_f64(effort) * Dd::from_f64(velocity) * Dd::from_f64(controls.duration);
        }
        for (link, (description, body)) in system
            .model
            .links()
            .iter()
            .zip(system.loads.external_body_wrench)
            .enumerate()
        {
            if link % 64 == 0 {
                poll(cancelled)?;
            }
            let velocity = middle[index].world.body_twist[link];
            accounts.body = accounts.body
                + impulse_pairing(*body, velocity, velocity) * Dd::from_f64(controls.duration);
            if system.loads.gravity_world != Vec3::default() {
                let displacement = com_displacement(
                    before[index].base.world_from_base,
                    after[index].base.world_from_base,
                    before[index].shape.world_from_link[link],
                    after[index].shape.world_from_link[link],
                    description.inertia().center_of_mass(),
                )
                .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
                for (delta, gravity) in displacement
                    .into_iter()
                    .zip(components(system.loads.gravity_world))
                {
                    accounts.gravity = accounts.gravity
                        + Dd::from_f64(description.inertia().mass())
                            * Dd::from_f64(gravity)
                            * delta;
                }
            }
        }
    }
    let mut paths = Vec::with_capacity(prepared.len());
    let mut slips = Vec::with_capacity(prepared.len());
    let mut midpoint_gaps = Vec::with_capacity(prepared.len());
    for (p, force) in prepared.iter().zip(forces) {
        poll(cancelled)?;
        let index = p.contact.point.system_index;
        let end_parts = point_parts(&after[index], p.contact.point)
            .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
        let mid_parts = point_parts(&middle[index], p.contact.point)
            .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
        let displacement = difference(end_parts, p.initial_parts);
        let midpoint_displacement = difference(mid_parts, p.initial_parts);
        let relative = p.axes.map(|axis| {
            dot_dd(axis, displacement)
                - dot(axis, p.contact.plane_velocity_world) * Dd::from_f64(controls.duration)
        });
        paths.push([
            (p.initial_gap + relative[0]).to_f64(),
            relative[1].to_f64(),
            relative[2].to_f64(),
        ]);
        slips.push(relative);
        midpoint_gaps.push(
            (p.initial_gap + dot_dd(p.axes[0], midpoint_displacement)
                - dot(p.axes[0], p.contact.plane_velocity_world)
                    * Dd::from_f64(0.5 * controls.duration))
            .to_f64(),
        );
        let velocity = velocity_parts(&middle[index], p.contact.point)
            .map_err(|error| ArticulatedSustainedError::Geometry(error.into()))?;
        for (&force, axis) in force.iter().zip(p.axes) {
            accounts.contact_quadrature = accounts.contact_quadrature
                + Dd::from_f64(force)
                    * (dot(axis, velocity[0]) + dot(axis, velocity[1]))
                    * Dd::from_f64(controls.duration);
        }
    }
    Ok(Response {
        states,
        midpoints,
        motion,
        paths,
        relative_displacement: slips,
        midpoint_gaps,
        accounts,
    })
}

fn frames(
    systems: &[System<'_>],
    states: &[State],
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Vec<Frame>, NativeError> {
    systems
        .iter()
        .zip(states)
        .enumerate()
        .map(|(index, (system, state))| {
            poll(cancelled)?;
            let input = native_input(system.model, state);
            let view = input.system_view();
            let world = forward_kinematics(
                system.model,
                view.base,
                input.joint_position,
                input.joint_velocity,
            )
            .map_err(|error| ArticulatedSustainedError::Dynamics {
                system_index: index,
                error,
            })?;
            poll(cancelled)?;
            let shape =
                relative_kinematics(system.model, input.joint_position).map_err(|error| {
                    ArticulatedSustainedError::Dynamics {
                        system_index: index,
                        error,
                    }
                })?;
            Ok(Frame {
                state: state.clone(),
                base: view.base,
                world,
                shape,
            })
        })
        .collect()
}
fn point_parts(frame: &Frame, point: LinkContactPoint) -> Result<[Vec3; 3], fs_ga::GaError> {
    let relative = frame.shape.world_from_link[point.link_index];
    let rotation = frame.base.world_from_base.rotation();
    Ok([
        frame.base.world_from_base.translation(),
        rotation.rotate(relative.translation())?,
        rotation.rotate(relative.rotation().rotate(point.point_body)?)?,
    ])
}
fn velocity_parts(frame: &Frame, point: LinkContactPoint) -> Result<[Vec3; 2], fs_ga::GaError> {
    let rotation = frame.world.world_from_link[point.link_index].rotation();
    let velocity = frame.world.body_twist[point.link_index];
    Ok([
        rotation.rotate(velocity.linear)?,
        rotation.rotate(velocity.angular.cross(point.point_body))?,
    ])
}
fn relative_speed(
    frame: &Frame,
    contact: ArticulatedSustainedPlaneContact,
    axes: [Vec3; 3],
) -> Result<[f64; 3], fs_ga::GaError> {
    let parts = velocity_parts(frame, contact.point)?;
    Ok(axes.map(|axis| {
        (dot(axis, parts[0]) + dot(axis, parts[1]) - dot(axis, contact.plane_velocity_world))
            .to_f64()
    }))
}
fn point_gap(parts: [Vec3; 3], plane: Vec3, normal: Vec3) -> Dd {
    dot(normal, parts[0]) - dot(normal, plane) + dot(normal, parts[1]) + dot(normal, parts[2])
}
fn difference(last: [Vec3; 3], first: [Vec3; 3]) -> [Dd; 3] {
    core::array::from_fn(|axis| {
        last.iter().zip(first).fold(Dd::ZERO, |sum, (b, a)| {
            sum + Dd::from_f64(components(*b)[axis]) - Dd::from_f64(components(a)[axis])
        })
    })
}
fn dot(axis: Vec3, value: Vec3) -> Dd {
    dot_dd(axis, components(value).map(Dd::from_f64))
}
fn dot_dd(axis: Vec3, value: [Dd; 3]) -> Dd {
    components(axis)
        .into_iter()
        .zip(value)
        .fold(Dd::ZERO, |s, (a, b)| s + Dd::from_f64(a) * b)
}
fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn vector(value: [Dd; 3]) -> Vec3 {
    Vec3::new(value[0].to_f64(), value[1].to_f64(), value[2].to_f64())
}

fn check_gap<E>(contact_id: u64, gap: f64) -> Result<(), ArticulatedSustainedError<E>> {
    if !gap.is_finite() {
        return Err(ArticulatedSustainedError::Invalid("nonfinite contact gap"));
    }
    if gap < 0.0 {
        return Err(ArticulatedSustainedError::Penetration {
            contact_id,
            gap_m: gap,
        });
    }
    Ok(())
}
fn check_contact<E>(
    p: &Prepared,
    force: [f64; 3],
    path: [f64; 3],
    midpoint_gap: f64,
    c: ArticulatedSustainedControls,
) -> Result<(), ArticulatedSustainedError<E>> {
    check_gap(p.contact.id, path[0])?;
    check_gap(p.contact.id, midpoint_gap)?;
    let fail = |reason| ArticulatedSustainedError::ContactLaw {
        contact_id: p.contact.id,
        reason,
    };
    if force[0] < 0.0 {
        return Err(fail("attractive normal force"));
    }
    let capacity = p.contact.friction_coefficient * force[0];
    let tangent = force[1].hypot(force[2]);
    let work = Dd::from_f64(force[1]) * Dd::from_f64(path[1])
        + Dd::from_f64(force[2]) * Dd::from_f64(path[2]);
    let duality = Dd::from_f64(capacity) * Dd::from_f64(path[1].hypot(path[2])) + work;
    let normal_complementarity = (Dd::from_f64(force[0]) * Dd::from_f64(path[0])).to_f64();
    if !capacity.is_finite()
        || !tangent.is_finite()
        || !duality.to_f64().is_finite()
        || !normal_complementarity.is_finite()
        || tangent - capacity > c.force_tolerance_n
        || duality.to_f64().abs() > c.complementarity_tolerance_j
        || normal_complementarity.abs() > c.complementarity_tolerance_j
    {
        return Err(fail("force cone or configuration/slip complementarity"));
    }
    if force[0] > c.force_tolerance_n
        && ([p.initial_gap.to_f64(), midpoint_gap, path[0]]
            .iter()
            .any(|&gap| gap > c.loaded_gap_tolerance_m)
            || p.initial_speed[0].abs() > c.normal_velocity_tolerance_m_s)
    {
        return Err(fail("loaded point is outside the sustained contact regime"));
    }
    Ok(())
}

fn admit<E>(
    systems: &[System<'_>],
    contacts: &[ArticulatedSustainedPlaneContact],
    forces: &[[f64; 3]],
    time: f64,
    c: ArticulatedSustainedControls,
) -> Result<f64, ArticulatedSustainedError<E>> {
    let end = time + c.duration;
    if systems.is_empty()
        || systems.len() > 128
        || contacts.is_empty()
        || contacts.len() > solve::MAX_CONTACTS
        || contacts.len() > c.solve.maximum_contacts
        || forces.len() != contacts.len()
        || !time.is_finite()
        || !end.is_finite()
        || end <= time
        || !finite_vector(c.angular_origin_world)
        || [
            c.duration,
            c.loaded_gap_tolerance_m,
            c.force_tolerance_n,
            c.normal_velocity_tolerance_m_s,
            c.complementarity_tolerance_j,
            c.energy_tolerance_j,
            c.linear_momentum_tolerance_ns,
            c.angular_momentum_tolerance_nms,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(ArticulatedSustainedError::Invalid(
            "finite advancing clock, bounded complete inputs and positive independent gates required",
        ));
    }
    let mut links = 0usize;
    for system in systems {
        links = links
            .checked_add(system.model.link_count())
            .ok_or(ArticulatedSustainedError::Invalid("link count overflow"))?;
        if links > MAX_FREE_FLOATING_STEP_LINKS {
            return Err(ArticulatedSustainedError::Invalid(
                "native total link cap exceeded",
            ));
        }
        if let State::Prescribed(state) = system.state {
            if state.time_s != time {
                return Err(ArticulatedSustainedError::Invalid(
                    "prescribed clock differs from common clock",
                ));
            }
        }
    }
    Ok(end)
}
fn poll<E>(
    cancelled: &mut (impl FnMut() -> bool + ?Sized),
) -> Result<(), ArticulatedSustainedError<E>> {
    if cancelled() {
        Err(ArticulatedSustainedError::Cancelled)
    } else {
        Ok(())
    }
}
fn free_error<E>(
    system_index: usize,
    error: FreeFloatingMotionError,
) -> ArticulatedSustainedError<E> {
    if matches!(error, FreeFloatingMotionError::Cancelled) {
        ArticulatedSustainedError::Cancelled
    } else {
        ArticulatedSustainedError::FreeMotion {
            system_index,
            error,
        }
    }
}
fn prescribed_error<E>(
    system_index: usize,
    error: PrescribedBaseMotionError,
) -> ArticulatedSustainedError<E> {
    if matches!(error, PrescribedBaseMotionError::Cancelled) {
        ArticulatedSustainedError::Cancelled
    } else {
        ArticulatedSustainedError::PrescribedMotion {
            system_index,
            error,
        }
    }
}
fn step_error<E>(
    error: ArticulatedContactStepError<core::convert::Infallible>,
) -> ArticulatedSustainedError<E> {
    match error {
        ArticulatedContactStepError::Cancelled => ArticulatedSustainedError::Cancelled,
        ArticulatedContactStepError::Dynamics {
            system_index,
            error,
        } => ArticulatedSustainedError::Dynamics {
            system_index,
            error,
        },
        _ => ArticulatedSustainedError::Invalid("native momentum total refused"),
    }
}
fn lift<E>(error: NativeError) -> ArticulatedSustainedError<E> {
    match error {
        ArticulatedSustainedError::Invalid(s) => ArticulatedSustainedError::Invalid(s),
        ArticulatedSustainedError::Geometry(e) => ArticulatedSustainedError::Geometry(e),
        ArticulatedSustainedError::Dynamics {
            system_index,
            error,
        } => ArticulatedSustainedError::Dynamics {
            system_index,
            error,
        },
        ArticulatedSustainedError::FreeMotion {
            system_index,
            error,
        } => free_error(system_index, error),
        ArticulatedSustainedError::PrescribedMotion {
            system_index,
            error,
        } => prescribed_error(system_index, error),
        ArticulatedSustainedError::Cancelled => ArticulatedSustainedError::Cancelled,
        _ => ArticulatedSustainedError::Invalid("unexpected native frame refusal"),
    }
}
