//! Automatic encounters of free and supported articulated spherical carriers.
//!
//! Coincident link origins/COMs and fixed/revolute joints make each sphere
//! centre independent of internal rotation. Free centres follow uniform gravity;
//! prescribed centres have constant world acceleration and their support frames
//! have constant world angular velocity. Every relative centre path is quadratic.
//! The original contact owner checks its complete path, including all interior
//! distance extrema. This is floating analytic geometry, not an interval tube
//! or an arbitrary articulated collision certificate.

pub mod planes;
pub mod translating;

use super::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectoryError, ArticulatedTrajectoryMomentum,
    Impulses, MixedArticulatedTrajectoryEventResult, MixedArticulatedTrajectoryResult,
    MixedArticulatedTrajectoryState as State, MixedArticulatedTrajectorySystem as System,
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
    advance_mixed_articulated_events, finite_vector, poll,
};
use crate::articulated::contact_geometry::{ContactGeometryError, prepare_contact_frame};
use crate::articulated::{BaseState, JointKind};
use crate::articulated_contacts::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
};
use crate::articulated_contacts::mixed::{
    MixedArticulatedImpulseResponse, resolve_mixed_articulated_frictional_contacts,
};
use crate::articulated_contacts::{
    ArticulatedContactError, ContactCounterpart, LinkContactPoint, unit_normal,
};
use crate::frictional_impact::spheres::{ContactSphere, SpherePairLaw};
use fs_contact::frictional_multicontact::{FrictionalContactResult, NormalVelocityLaw};
use fs_contact::swept_sphere::accelerated::{AcceleratedSphereInput, sphere_path};
use fs_contact::swept_sphere::{SweptSphereError, SweptSphereInput};
use fs_ga::{GaError, Se3, So3, So3Tangent, Twist, Vec3};
use fs_math::dd::Dd;
use planes::{
    AcceleratingPlane, MixedSphereSceneContactKey, MixedSphereSceneEventGeometry,
    MixedSphereSceneResult, PreparedPlanes, SpherePlaneLaw,
};
use translating::{CentrePath, PreparedGeometry};

type Impact =
    ArticulatedFrictionalContactResult<FrictionalContactResult, MixedArticulatedImpulseResponse>;

/// Whole-window numerical and geometric bounds.
#[derive(Debug, Clone, Copy)]
pub struct MixedSphereControls {
    /// Native smooth accuracy, instantaneous law controls and momentum origin.
    /// The interval budget counts every attempted common smooth step, including
    /// rejected floating-time neighbours.
    pub trajectory: ArticulatedTrajectoryConfig,
    /// Positive nonpenetrating activation resolution, in metres. All pairs
    /// inside this band enter the same event, including separating neighbours.
    pub event_gap_tolerance_m: f64,
    /// Complete unordered-pair count, bounded by 21 (at most seven spheres).
    pub maximum_pairs: usize,
    /// At most 64 simultaneous events. Zero permits only event-free windows.
    pub maximum_events: usize,
    /// Complete-path and material-witness observations, in 1..=1_000_000.
    pub maximum_geometry_queries: usize,
}

/// Geometry accompanying one original simultaneous event.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedSphereEventGeometry {
    /// Stable sphere IDs in canonical pair order, matching the native contact IDs.
    pub pairs: Vec<[u64; 2]>,
    /// Largest actual centre/material activation gap at this event, in metres.
    pub maximum_gap_m: f64,
}

/// Complete native trajectory and independently checked sphere geometry.
#[derive(Debug, Clone)]
pub struct MixedSphereResult {
    /// Original native endpoints, reactions, actual increments and signed work.
    pub trajectory: MixedArticulatedTrajectoryResult,
    /// One geometry entry per trajectory event, in the same chronological order.
    pub event_geometry: Vec<MixedSphereEventGeometry>,
    pub smooth_steps_attempted: usize,
    pub geometry_queries: usize,
    /// Minimum complete accepted-path and actual surface-witness separation.
    /// This floating observation is not an interval enclosure.
    pub minimum_gap_m: f64,
}

/// Atomic refusal of the complete window.
#[derive(Debug)]
pub enum MixedSphereError {
    Invalid(&'static str),
    Trajectory(ArticulatedTrajectoryError),
    Geometry(GaError),
    Path {
        pair: [u64; 2],
        time: f64,
        error: SweptSphereError,
    },
    PlanePath {
        sphere: u64,
        plane: u64,
        time: f64,
        error: fs_contact::ballistic::BallisticContactError,
    },
    /// Near-rest inward contact needs a sustained-contact owner.
    UnresolvedContact {
        pair: [u64; 2],
        time: f64,
    },
    /// Near-rest inward plane contact needs a sustained-contact owner.
    UnresolvedPlaneContact {
        sphere: u64,
        plane: u64,
        time: f64,
    },
    /// Original Cartesian acceleration cannot resolve the admitted affine
    /// native joint path. Residual and allowance have units of m/s².
    UnresolvedLinkPath {
        system_index: usize,
        link_index: usize,
        time: f64,
        acceleration_residual: f64,
        allowed_acceleration_residual: f64,
    },
    /// No nonpenetrating native endpoint/witness was found in 64 time neighbours.
    UnresolvedGeometry {
        time: f64,
    },
    NoEventProgress,
    SmoothBudget,
    EventBudget,
    GeometryBudget,
}
impl core::fmt::Display for MixedSphereError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "mixed articulated spheres: {self:?}")
    }
}
impl std::error::Error for MixedSphereError {}
impl From<ArticulatedTrajectoryError> for MixedSphereError {
    fn from(error: ArticulatedTrajectoryError) -> Self {
        Self::Trajectory(error)
    }
}
impl From<GaError> for MixedSphereError {
    fn from(error: GaError) -> Self {
        Self::Geometry(error)
    }
}

#[derive(Clone, Copy)]
struct CentreCommand {
    acceleration_world: Vec3,
    angular_velocity_world: Vec3,
}
#[derive(Clone, Copy)]
struct Pair {
    first: usize,
    second: usize,
    id: u64,
    law: SpherePairLaw,
}
impl Pair {
    fn ids(self) -> [u64; 2] {
        [self.law.first, self.law.second]
    }
}
#[derive(Default)]
struct Work {
    smooth: usize,
    accepted: usize,
    geometry: usize,
}
impl Work {
    fn query(&mut self, controls: MixedSphereControls) -> Result<(), MixedSphereError> {
        if self.geometry == controls.maximum_geometry_queries {
            return Err(MixedSphereError::GeometryBudget);
        }
        self.geometry += 1;
        Ok(())
    }
}

/// Advance an automatically located sequence of mixed-boundary sphere impacts.
///
/// Spheres correspond to systems in input order. Every unordered stable-ID pair
/// requires exactly one SpherePairLaw; pair IDs need no caller enumeration.
/// Each declared root-attached sphere encloses its complete original tree.
/// Link origins and COMs must coincide, with fixed/revolute joints only. Mass
/// moments are checked against the radius; they do not prove unknown shapes fit.
///
/// Free systems retain native base recoil, torque loads and actuator dynamics;
/// only gravity may supply their translation. A prescribed system may have a
/// movable revolute root. Its initial base data define constant WORLD angular
/// velocity and constant WORLD centre acceleration. Initial angular acceleration
/// must be zero. Native support samples use a_body = R^T*a_world - omega x v_body.
/// This command is evaluated directly, not inferred by interpolating samples.
///
/// Every complete quadratic pair path uses the existing accelerated-sphere
/// owner. Actual native endpoints and regenerated root-link material points
/// must also remain nonpenetrating. At most 64 earlier floating-time neighbours
/// are tried; no pose, arm, speed or inertia is projected or repaired. One native
/// simultaneous Coulomb event updates every participating tree before discovery
/// restarts from the accepted state.
///
/// All input states remain immutable on failure or cancellation. Numerical
/// smooth balance defects remain explicit, separate from impact loss. This lane
/// does not resolve sustained contact, general offset articulated geometry,
/// interval-certified CCD, arbitrary support commands or heat partition.
pub fn advance_mixed_sphere_contacts(
    systems: &[System<'_>],
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    start_time_s: f64,
    controls: MixedSphereControls,
    cancelled: impl FnMut() -> bool,
) -> Result<MixedSphereResult, MixedSphereError> {
    let scene = advance_scene(
        systems,
        spheres,
        laws,
        &[],
        &[],
        start_time_s,
        controls,
        cancelled,
    )?;
    let event_geometry = scene
        .event_geometry
        .into_iter()
        .map(|event| {
            let pairs = event
                .contacts
                .into_iter()
                .map(|key| match key {
                    MixedSphereSceneContactKey::SpherePair(pair) => Ok(pair),
                    MixedSphereSceneContactKey::SpherePlane { .. } => Err(
                        MixedSphereError::Invalid("plane appeared in a pair-only scene"),
                    ),
                })
                .collect::<Result<_, _>>()?;
            Ok(MixedSphereEventGeometry {
                pairs,
                maximum_gap_m: event.maximum_gap_m,
            })
        })
        .collect::<Result<_, MixedSphereError>>()?;
    Ok(MixedSphereResult {
        trajectory: scene.trajectory,
        event_geometry,
        smooth_steps_attempted: scene.smooth_steps_attempted,
        geometry_queries: scene.geometry_queries,
        minimum_gap_m: scene.minimum_gap_m,
    })
}

#[allow(clippy::too_many_arguments)]
fn advance_scene(
    systems: &[System<'_>],
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    planes: &[AcceleratingPlane],
    plane_laws: &[SpherePlaneLaw],
    start_time_s: f64,
    controls: MixedSphereControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixedSphereSceneResult, MixedSphereError> {
    poll(&mut cancelled)?;
    let geometry = PreparedGeometry::carriers(
        systems,
        spheres,
        laws,
        start_time_s,
        controls,
        !planes.is_empty(),
    )?;
    advance_prepared(
        systems,
        geometry,
        planes,
        plane_laws,
        start_time_s,
        controls,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_prepared(
    systems: &[System<'_>],
    geometry: PreparedGeometry,
    planes: &[AcceleratingPlane],
    plane_laws: &[SpherePlaneLaw],
    start_time_s: f64,
    controls: MixedSphereControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixedSphereSceneResult, MixedSphereError> {
    let spheres = &geometry.spheres;
    let pairs = &geometry.pairs;
    let commands = &geometry.commands;
    let plane_geometry = PreparedPlanes::admit(spheres, planes, plane_laws, pairs.len())?;
    let mut states: Vec<_> = systems.iter().map(|system| system.state.clone()).collect();
    let origin = controls.trajectory.angular_origin_world;
    let before = super::momentum(systems, &states, origin, 0.0, &mut cancelled)?;
    let mut totals = Totals::default();
    let mut work = Work::default();
    let mut events = Vec::<MixedArticulatedTrajectoryEventResult>::new();
    let mut event_geometry = Vec::new();
    let mut minimum_gap = f64::INFINITY;
    let mut elapsed = 0.0;
    let mut clock = start_time_s;
    while elapsed < controls.trajectory.duration {
        poll(&mut cancelled)?;
        let paths = geometry.paths(systems, &states, elapsed, &mut cancelled)?;
        let horizon =
            (controls.trajectory.duration - elapsed).min(controls.trajectory.maximum_smooth_step);
        let mut encounter: Option<(f64, u64)> = None;
        for &pair in pairs {
            poll(&mut cancelled)?;
            work.query(controls)?;
            let path =
                sphere_path(path_input(pair, spheres, &paths, horizon)?).map_err(|error| {
                    MixedSphereError::Path {
                        pair: pair.ids(),
                        time: elapsed,
                        error,
                    }
                })?;
            if let Some(hit) = path.first_entering_contact {
                if hit.relative_normal_speed.abs() <= controls.trajectory.contact.velocity_tolerance
                {
                    return Err(MixedSphereError::UnresolvedContact {
                        pair: pair.ids(),
                        time: elapsed + hit.time_s,
                    });
                }
                if encounter.is_none_or(|(earliest, _)| hit.time_s < earliest) {
                    encounter = Some((hit.time_s, pair.id));
                }
            }
        }
        if let Some(hit) = plane_geometry.encounter(
            spheres,
            &paths,
            elapsed,
            horizon,
            controls,
            &mut work,
            &mut cancelled,
        )? {
            if encounter.is_none_or(|(earliest, _)| hit.0 < earliest) {
                encounter = Some(hit);
            }
        }
        let mut duration = encounter.map_or(horizon, |(time, _)| time);
        let mut accepted = None;
        for _ in 0..=64 {
            poll(&mut cancelled)?;
            let end = elapsed + duration;
            if !end.is_finite() || (duration > 0.0 && end <= elapsed) {
                return Err(MixedSphereError::UnresolvedGeometry { time: end });
            }
            let step = if duration > 0.0 {
                if work.smooth == controls.trajectory.maximum_smooth_intervals {
                    return Err(MixedSphereError::SmoothBudget);
                }
                work.smooth += 1;
                let inputs: Vec<_> = systems
                    .iter()
                    .zip(&states)
                    .map(|(system, state)| System {
                        model: system.model,
                        state,
                        loads: system.loads,
                    })
                    .collect();
                Some(advance_mixed_articulated_events(
                    &inputs,
                    &[],
                    clock,
                    ArticulatedTrajectoryConfig {
                        duration,
                        maximum_smooth_intervals: 1,
                        ..controls.trajectory
                    },
                    |index, current, h| support_path(current, commands[index], h),
                    &mut cancelled,
                )?)
            } else {
                None
            };
            let candidate = step.as_ref().map_or(states.as_slice(), |s| &s.states_after);
            let observation = observe(
                systems,
                &geometry,
                &plane_geometry,
                &paths,
                candidate,
                duration,
                elapsed,
                encounter.map(|(_, id)| id),
                controls,
                &mut work,
                &mut cancelled,
            )?;
            if let Some(observation) = observation {
                accepted = Some((step, end, observation));
                break;
            }
            if duration == 0.0 {
                break;
            }
            duration = duration.next_down();
        }
        let (step, end, observation) = accepted.ok_or(MixedSphereError::UnresolvedGeometry {
            time: elapsed + duration,
        })?;
        if let Some(step) = step {
            totals.smooth(&step);
            clock = step.time_after_s;
            states = step.states_after;
            work.accepted += step.smooth_intervals;
        }
        elapsed = end;
        minimum_gap = minimum_gap.min(observation.minimum_gap);
        if !observation.contacts.is_empty() {
            if events.len() == controls.maximum_events {
                return Err(MixedSphereError::EventBudget);
            }
            if events.last().is_some_and(|event| event.time == elapsed) {
                return Err(MixedSphereError::NoEventProgress);
            }
            let inputs: Vec<_> = systems
                .iter()
                .zip(&states)
                .map(|(system, state)| state.contact_system(system.model))
                .collect();
            let impact = resolve_mixed_articulated_frictional_contacts(
                &inputs,
                &observation.contacts,
                controls.trajectory.contact,
                &mut cancelled,
            )
            .map_err(|error| {
                if matches!(error, ArticulatedContactError::Cancelled) {
                    ArticulatedTrajectoryError::Cancelled
                } else {
                    ArticulatedTrajectoryError::Contact {
                        event_index: events.len(),
                        time: elapsed,
                        error,
                    }
                }
            })?;
            totals.impact(&impact, origin);
            for (state, response) in states.iter_mut().zip(&impact.responses) {
                match (state, response) {
                    (State::Free(state), MixedArticulatedImpulseResponse::Free(response)) => {
                        state.base = response.base_after;
                        state
                            .joint_velocity
                            .clone_from(&response.joint_velocity_after);
                    }
                    (
                        State::Prescribed(state),
                        MixedArticulatedImpulseResponse::Prescribed(response),
                    ) => {
                        state.base = response.base_after;
                        state
                            .joint_velocity
                            .clone_from(&response.joint_velocity_after);
                    }
                    _ => return Err(MixedSphereError::Invalid("native boundary kind changed")),
                }
            }
            events.push(MixedArticulatedTrajectoryEventResult {
                time: elapsed,
                time_s: clock,
                impact,
            });
            event_geometry.push(MixedSphereSceneEventGeometry {
                contacts: observation.keys,
                maximum_gap_m: observation.maximum_gap,
            });
        }
    }
    let after = super::momentum(systems, &states, origin, elapsed, &mut cancelled)?;
    let trajectory = totals.finish(
        before,
        after,
        states,
        events,
        start_time_s,
        clock,
        work.accepted,
        origin,
    )?;
    let planes_after = plane_geometry.endpoints(elapsed)?;
    poll(&mut cancelled)?;
    Ok(MixedSphereSceneResult {
        trajectory,
        event_geometry,
        planes_after,
        smooth_steps_attempted: work.smooth,
        geometry_queries: work.geometry,
        minimum_gap_m: minimum_gap,
    })
}

fn admit(
    systems: &[System<'_>],
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    clock: f64,
    controls: MixedSphereControls,
    allow_single: bool,
) -> Result<(Vec<Pair>, Vec<CentreCommand>), MixedSphereError> {
    admit_window(systems, clock, controls)?;
    if (systems.len() < 2 && !allow_single) || systems.len() != spheres.len() {
        return Err(MixedSphereError::Invalid(
            "one enclosing carrier sphere per system required",
        ));
    }
    let pairs = admit_pairs(spheres, laws, controls)?;
    let commands = systems
        .iter()
        .zip(spheres)
        .map(|(system, sphere)| admit_carrier(system, *sphere))
        .collect::<Result<_, _>>()?;
    Ok((pairs, commands))
}

fn admit_window(
    systems: &[System<'_>],
    clock: f64,
    controls: MixedSphereControls,
) -> Result<(), MixedSphereError> {
    super::event_schedule(systems.len(), &[], controls.trajectory)?;
    super::admit_states(systems, clock)?;
    if systems.is_empty()
        || !(1..=21).contains(&controls.maximum_pairs)
        || controls.maximum_events > 64
        || !(1..=1_000_000).contains(&controls.maximum_geometry_queries)
        || !controls.event_gap_tolerance_m.is_finite()
        || controls.event_gap_tolerance_m <= 0.0
        || !controls.trajectory.contact.velocity_tolerance.is_finite()
        || controls.trajectory.contact.velocity_tolerance <= 0.0
    {
        return Err(MixedSphereError::Invalid(
            "invalid sphere scene or work bounds",
        ));
    }
    Ok(())
}

fn admit_carrier(
    system: &System<'_>,
    sphere: ContactSphere,
) -> Result<CentreCommand, MixedSphereError> {
    if !sphere.radius_m.is_finite() || sphere.radius_m <= 0.0 {
        return Err(MixedSphereError::Invalid(
            "positive finite sphere radius required",
        ));
    }
    for link in system.model.links() {
        if link.inertia().center_of_mass() != Vec3::default()
            || link.parent_from_child_zero().translation() != Vec3::default()
            || !matches!(link.joint().kind(), JointKind::Fixed | JointKind::Revolute)
        {
            return Err(MixedSphereError::Invalid(
                "coincident carrier origins/COMs and fixed/revolute joints required",
            ));
        }
        let inertia = link.inertia();
        let matrix = inertia.inertia_com().m;
        let scale = matrix[0].max(matrix[4]).max(matrix[8]);
        let rms = (scale.sqrt() / inertia.mass().sqrt())
            * (0.5 * (matrix[0] / scale + matrix[4] / scale + matrix[8] / scale)).sqrt();
        if !rms.is_finite() || rms > sphere.radius_m {
            return Err(MixedSphereError::Invalid(
                "inertia exceeds spherical enclosure",
            ));
        }
    }
    if system.loads.external_body_wrench.len() != system.model.link_count()
        || system.loads.generalized_force.len() != system.model.dof_count()
    {
        return Err(MixedSphereError::Invalid(
            "native load dimensions must match model",
        ));
    }
    let command = match system.state {
        State::Free(_) => {
            if system.model.links()[0].joint().kind() != JointKind::Fixed
                || system
                    .loads
                    .external_body_wrench
                    .iter()
                    .any(|w| w.force != Vec3::default())
            {
                return Err(MixedSphereError::Invalid(
                    "free carrier requires a fixed root and gravity-only translation",
                ));
            }
            CentreCommand {
                acceleration_world: system.loads.gravity_world,
                angular_velocity_world: Vec3::default(),
            }
        }
        State::Prescribed(state) => {
            if state.base.acceleration_body.angular != Vec3::default() {
                return Err(MixedSphereError::Invalid(
                    "prescribed carrier requires constant world angular velocity",
                ));
            }
            let base = state.base;
            CentreCommand {
                acceleration_world: base.world_from_base.rotation().rotate(
                    base.acceleration_body.linear
                        + base.twist_body.angular.cross(base.twist_body.linear),
                )?,
                angular_velocity_world: base
                    .world_from_base
                    .rotation()
                    .rotate(base.twist_body.angular)?,
            }
        }
    };
    if !finite_vector(command.acceleration_world) || !finite_vector(command.angular_velocity_world)
    {
        return Err(MixedSphereError::Invalid(
            "nonfinite centre/support command",
        ));
    }
    Ok(command)
}

fn admit_pairs(
    spheres: &[ContactSphere],
    laws: &[SpherePairLaw],
    controls: MixedSphereControls,
) -> Result<Vec<Pair>, MixedSphereError> {
    if spheres.is_empty() || spheres.len() > 7 {
        return Err(MixedSphereError::Invalid(
            "one to seven declared spheres required",
        ));
    }
    let count = spheres.len() * (spheres.len() - 1) / 2;
    if count > controls.maximum_pairs || laws.len() != count {
        return Err(MixedSphereError::Invalid(
            "complete bounded pair laws required",
        ));
    }
    let mut ids: Vec<_> = spheres.iter().map(|sphere| sphere.id).collect();
    ids.sort_unstable();
    if ids.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(MixedSphereError::Invalid("duplicate sphere identity"));
    }
    if spheres
        .iter()
        .any(|sphere| !sphere.radius_m.is_finite() || sphere.radius_m <= 0.0)
    {
        return Err(MixedSphereError::Invalid(
            "positive finite sphere radius required",
        ));
    }
    let mut ordered = laws.to_vec();
    ordered.sort_unstable_by_key(|law| (law.first, law.second));
    let mut pairs = Vec::with_capacity(count);
    for (id, law) in ordered.iter().enumerate() {
        if law.first >= law.second
            || !law.friction_coefficient.is_finite()
            || law.friction_coefficient < 0.0
            || !law.restitution.is_finite()
            || !(0.0..=1.0).contains(&law.restitution)
            || (id > 0
                && (ordered[id - 1].first, ordered[id - 1].second) == (law.first, law.second))
        {
            return Err(MixedSphereError::Invalid(
                "invalid or duplicate sphere pair law",
            ));
        }
        let first = spheres
            .iter()
            .position(|sphere| sphere.id == law.first)
            .ok_or(MixedSphereError::Invalid(
                "pair law references unknown sphere",
            ))?;
        let second = spheres
            .iter()
            .position(|sphere| sphere.id == law.second)
            .ok_or(MixedSphereError::Invalid(
                "pair law references unknown sphere",
            ))?;
        if !(spheres[first].radius_m + spheres[second].radius_m).is_finite() {
            return Err(MixedSphereError::Invalid("unrepresentable radius sum"));
        }
        pairs.push(Pair {
            first,
            second,
            id: id as u64,
            law: *law,
        });
    }
    Ok(pairs)
}

fn base(state: &State) -> BaseState {
    match state {
        State::Free(state) => BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        State::Prescribed(state) => state.base,
    }
}
fn centre_velocity(state: &State) -> Result<Vec3, GaError> {
    let base = base(state);
    base.world_from_base
        .rotation()
        .rotate(base.twist_body.linear)
}
fn coordinates(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn path_input(
    pair: Pair,
    spheres: &[ContactSphere],
    centres: &[CentrePath],
    duration: f64,
) -> Result<AcceleratedSphereInput, MixedSphereError> {
    Ok(AcceleratedSphereInput {
        motion: SweptSphereInput {
            relative_position_world: coordinates(
                centres[pair.first].position_world - centres[pair.second].position_world,
            ),
            relative_velocity_world: coordinates(
                centres[pair.first].velocity_world - centres[pair.second].velocity_world,
            ),
            radius_sum: spheres[pair.first].radius_m + spheres[pair.second].radius_m,
            duration_s: duration,
        },
        relative_acceleration_world: coordinates(
            centres[pair.first].acceleration_world - centres[pair.second].acceleration_world,
        ),
    })
}

fn support_path(
    state: &PrescribedBaseMotionState,
    command: CentreCommand,
    duration: f64,
) -> Result<PrescribedBasePath, PrescribedBaseMotionError> {
    let initial = state.base;
    let velocity = initial
        .world_from_base
        .rotation()
        .rotate(initial.twist_body.linear)?;
    let sample = |time: f64| -> Result<BaseState, PrescribedBaseMotionError> {
        let rotation = if command.angular_velocity_world == Vec3::default() {
            // A translating support retains the exact stored orientation;
            // composing an identity must not re-round a fixed guide frame.
            initial.world_from_base.rotation()
        } else {
            So3::exp(So3Tangent {
                angular: command.angular_velocity_world.scale(time),
            })?
            .compose(initial.world_from_base.rotation())?
        };
        let position = initial.world_from_base.translation()
            + velocity.scale(time)
            + command.acceleration_world.scale(time).scale(0.5 * time);
        let world_velocity = velocity + command.acceleration_world.scale(time);
        let omega = rotation.inverse().rotate(command.angular_velocity_world)?;
        let linear = rotation.inverse().rotate(world_velocity)?;
        let acceleration =
            rotation.inverse().rotate(command.acceleration_world)? - omega.cross(linear);
        Ok(BaseState::prescribed(
            Se3::from_parts(rotation, position)?,
            Twist::new(omega, linear),
            Twist::new(Vec3::default(), acceleration),
        ))
    };
    Ok(PrescribedBasePath {
        midpoint: sample(0.5 * duration)?,
        endpoint: sample(duration)?,
    })
}

struct Observation {
    contacts: Vec<ArticulatedFrictionalContact>,
    keys: Vec<MixedSphereSceneContactKey>,
    minimum_gap: f64,
    maximum_gap: f64,
}

#[allow(clippy::too_many_arguments)]
fn observe(
    systems: &[System<'_>],
    geometry: &PreparedGeometry,
    planes: &PreparedPlanes,
    before: &[CentrePath],
    after: &[State],
    duration: f64,
    elapsed: f64,
    event_pair: Option<u64>,
    controls: MixedSphereControls,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<Observation>, MixedSphereError> {
    let spheres = &geometry.spheres;
    let pairs = &geometry.pairs;
    let snapshot = geometry.snapshot(
        systems,
        after,
        event_pair.is_some(),
        elapsed + duration,
        cancelled,
    )?;
    let kinematics = &snapshot.kinematics;
    let mut observation = Observation {
        contacts: Vec::new(),
        keys: Vec::new(),
        minimum_gap: f64::INFINITY,
        maximum_gap: 0.0,
    };
    for &pair in pairs {
        poll(cancelled)?;
        work.query(controls)?;
        let path = sphere_path(path_input(pair, spheres, before, duration)?).map_err(|error| {
            MixedSphereError::Path {
                pair: pair.ids(),
                time: elapsed,
                error,
            }
        })?;
        let relative = snapshot.centres[pair.first].position_world
            - snapshot.centres[pair.second].position_world;
        // Admit the ACTUAL native endpoint through the same geometry owner
        // that starts the next interval. Its deterministic norm can differ by
        // an ULP from std::hypot, and the native centre subtraction need not
        // equal the analytic path's rounded endpoint. Retry an earlier time
        // rather than publish a state the next original query would penetrate.
        work.query(controls)?;
        let endpoint = match sphere_path(path_input(pair, spheres, &snapshot.centres, 0.0)?) {
            Ok(path) => path,
            Err(SweptSphereError::InitialPenetration) => return Ok(None),
            Err(error) => {
                return Err(MixedSphereError::Path {
                    pair: pair.ids(),
                    time: elapsed + duration,
                    error,
                });
            }
        };
        let gap = endpoint.minimum_gap;
        if !gap.is_finite() || !path.minimum_gap.is_finite() {
            return Err(MixedSphereError::Invalid(
                "nonfinite complete sphere path gap",
            ));
        }
        observation.minimum_gap = observation.minimum_gap.min(path.minimum_gap).min(gap);
        if observation.minimum_gap < 0.0 {
            return Ok(None);
        }
        if event_pair.is_some() && gap <= controls.event_gap_tolerance_m {
            work.query(controls)?;
            let normal = unit_normal(relative)
                .map_err(|_| MixedSphereError::Invalid("unresolvable sphere normal"))?;
            let first_centre = geometry.attachments[pair.first].centre;
            let second_centre = geometry.attachments[pair.second].centre;
            let first_pose =
                kinematics[first_centre.system_index].world_from_link[first_centre.link_index];
            let second_pose =
                kinematics[second_centre.system_index].world_from_link[second_centre.link_index];
            let first_arm = first_centre.point_body
                + first_pose
                    .rotation()
                    .inverse()
                    .rotate(normal.scale(-spheres[pair.first].radius_m))?;
            let second_arm = second_centre.point_body
                + second_pose
                    .rotation()
                    .inverse()
                    .rotate(normal.scale(spheres[pair.second].radius_m))?;
            let first_point = first_pose.transform_point(first_arm)?;
            let second_point = second_pose.transform_point(second_arm)?;
            let material_gap = normal.dot(first_point - second_point);
            if !material_gap.is_finite() {
                return Err(MixedSphereError::Invalid("nonfinite actual sphere witness"));
            }
            observation.minimum_gap = observation.minimum_gap.min(material_gap);
            if material_gap < 0.0 {
                return Ok(None);
            }
            if material_gap <= controls.event_gap_tolerance_m {
                let axis = if normal.x.abs() < 0.8 {
                    Vec3::new(1.0, 0.0, 0.0)
                } else {
                    Vec3::new(0.0, 1.0, 0.0)
                };
                let contact = ArticulatedFrictionalContact {
                    id: pair.id,
                    point: LinkContactPoint {
                        point_body: first_arm,
                        ..first_centre
                    },
                    other: ContactCounterpart::Body(LinkContactPoint {
                        point_body: second_arm,
                        ..second_centre
                    }),
                    normal_world: relative,
                    tangent_world: normal.cross(axis),
                    activation_tolerance: controls.event_gap_tolerance_m,
                    friction_coefficient: pair.law.friction_coefficient,
                    normal_law: NormalVelocityLaw::NewtonRestitution {
                        coefficient: pair.law.restitution,
                    },
                };
                // The same original material-point and normal-line admission
                // used by the native event must succeed before accepting this
                // floating time. Retry a representable earlier neighbour;
                // never project an arm or relax the mechanical owner's gate.
                match prepare_contact_frame(
                    contact.point,
                    contact.other,
                    contact.normal_world,
                    contact.tangent_world,
                    pair.id as usize,
                    kinematics,
                ) {
                    Ok(_) => {}
                    Err(ContactGeometryError::NonCollinearContact { .. }) => return Ok(None),
                    Err(error) => {
                        return Err(ArticulatedTrajectoryError::Contact {
                            event_index: 0,
                            time: elapsed + duration,
                            error: error.into(),
                        }
                        .into());
                    }
                }
                observation.contacts.push(contact);
                observation
                    .keys
                    .push(MixedSphereSceneContactKey::SpherePair(pair.ids()));
                observation.maximum_gap = observation.maximum_gap.max(gap).max(material_gap);
            }
        }
    }
    if !planes.observe(
        spheres,
        &geometry.attachments,
        before,
        &snapshot.centres,
        kinematics,
        duration,
        elapsed,
        event_pair.is_some(),
        controls,
        work,
        cancelled,
        &mut observation,
    )? {
        return Ok(None);
    }
    if event_pair.is_some_and(|id| !observation.contacts.iter().any(|contact| contact.id == id)) {
        return Ok(None);
    }
    poll(cancelled)?;
    Ok(Some(observation))
}

#[derive(Default)]
struct Totals {
    actual: Impulses,
    sources: Impulses,
    smooth_momentum_defect: Impulses,
    kinetic: Dd,
    actuator: Dd,
    external: Dd,
    support: Dd,
    boundary: Dd,
    loss: Dd,
    smooth_defect: Dd,
}
impl Totals {
    fn smooth(&mut self, step: &MixedArticulatedTrajectoryResult) {
        self.actual.add(
            step.linear_momentum_change_world,
            step.angular_momentum_change_about_origin_world,
            Vec3::default(),
        );
        self.sources.add(
            step.external_linear_impulse_world,
            step.external_angular_impulse_about_origin_world,
            Vec3::default(),
        );
        self.smooth_momentum_defect.add(
            step.smooth_linear_momentum_residual,
            step.smooth_angular_momentum_residual,
            Vec3::default(),
        );
        self.kinetic = self.kinetic + Dd::from_f64(step.kinetic_energy_change);
        self.actuator = self.actuator + Dd::from_f64(step.actuator_work);
        self.external = self.external + Dd::from_f64(step.smooth_external_work);
        self.support = self.support + Dd::from_f64(step.smooth_support_work);
        self.smooth_defect = self.smooth_defect + Dd::from_f64(step.smooth_kinetic_work_residual);
    }

    fn impact(&mut self, impact: &Impact, origin: Vec3) {
        for response in &impact.responses {
            match response {
                MixedArticulatedImpulseResponse::Free(response) => {
                    self.actual.add(
                        response.linear_momentum_change_world,
                        response.angular_momentum_change_world,
                        origin,
                    );
                    self.sources.add(
                        response.external_linear_impulse_world,
                        response.external_angular_impulse_world,
                        origin,
                    );
                    self.kinetic = self.kinetic + Dd::from_f64(response.kinetic_energy_change);
                }
                MixedArticulatedImpulseResponse::Prescribed(response) => {
                    self.actual.add(
                        response.linear_momentum_change_world,
                        response.angular_momentum_change_world,
                        origin,
                    );
                    self.sources.add(
                        response.external_impulse_world.force,
                        response.external_impulse_world.torque,
                        origin,
                    );
                    self.sources.add(
                        response.support_impulse_world.force,
                        response.support_impulse_world.torque,
                        origin,
                    );
                    self.kinetic = self.kinetic + Dd::from_f64(response.kinetic_energy_change);
                }
            }
        }
        self.boundary = self.boundary + Dd::from_f64(impact.boundary_work);
        self.loss = self.loss + Dd::from_f64(impact.contact.dissipated_energy);
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        self,
        before: ArticulatedTrajectoryMomentum,
        after: ArticulatedTrajectoryMomentum,
        states: Vec<State>,
        events: Vec<MixedArticulatedTrajectoryEventResult>,
        start: f64,
        end: f64,
        intervals: usize,
        origin: Vec3,
    ) -> Result<MixedArticulatedTrajectoryResult, MixedSphereError> {
        let external = self.external + self.support + self.boundary;
        let energy_residual = (self.kinetic - self.actuator - external + self.loss).to_f64();
        let (linear_residual, angular_residual) = self.actual.difference(&self.sources);
        let values = [
            self.kinetic,
            self.actuator,
            self.external,
            self.support,
            self.boundary,
            external,
            self.loss,
            self.smooth_defect,
        ]
        .map(Dd::to_f64);
        if !values.iter().all(|v| v.is_finite())
            || !energy_residual.is_finite()
            || ![
                self.actual.linear,
                self.actual.angular,
                self.sources.linear,
                self.sources.angular,
                self.smooth_momentum_defect.linear,
                self.smooth_momentum_defect.angular,
                linear_residual,
                angular_residual,
            ]
            .into_iter()
            .all(finite_vector)
        {
            return Err(MixedSphereError::Invalid(
                "nonfinite complete mixed sphere balance",
            ));
        }
        Ok(MixedArticulatedTrajectoryResult {
            time_before_s: start,
            time_after_s: end,
            states_after: states,
            events,
            smooth_intervals: intervals,
            angular_origin_world: origin,
            momentum_before: before,
            momentum_after: after,
            linear_momentum_change_world: self.actual.linear,
            angular_momentum_change_about_origin_world: self.actual.angular,
            kinetic_energy_change: values[0],
            external_linear_impulse_world: self.sources.linear,
            external_angular_impulse_about_origin_world: self.sources.angular,
            actuator_work: values[1],
            smooth_external_work: values[2],
            smooth_support_work: values[3],
            contact_boundary_work: values[4],
            external_work: values[5],
            contact_dissipated_energy: values[6],
            smooth_kinetic_work_residual: values[7],
            smooth_linear_momentum_residual: self.smooth_momentum_defect.linear,
            smooth_angular_momentum_residual: self.smooth_momentum_defect.angular,
            kinetic_work_residual: energy_residual,
            linear_momentum_residual: linear_residual,
            angular_momentum_residual: angular_residual,
        })
    }
}
