//! Automatic impacts of a spherical shell carrying co-located articulated rotors.
//!
//! All link origins and centres of mass coincide; only fixed/revolute joints
//! and torque loads are admitted. The shell centre is consequently the whole
//! tree COM and follows a quadratic path under uniform gravity, independently
//! of internal rotation. Every declared translating plane is queried on that
//! complete path. Rotation, actuator work and impact recoil retain the original
//! articulated owners. This is floating analytic sphere/plane geometry, not an
//! interval CCD certificate or a sustained-contact integrator.

/// Automatic inelastic landing and sustained support with native articulated recoil.
pub mod support;

use super::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectoryError, ArticulatedTrajectoryEventResult,
    ArticulatedTrajectoryResult, ArticulatedTrajectorySystem, Impulses, admission, finite_vector,
    momentum, poll, smooth_error,
};
use crate::articulated::motion::{FreeFloatingMotionState, step_free_floating};
use crate::articulated::{BaseState, JointKind, forward_kinematics};
use crate::articulated_contacts::friction::{
    ArticulatedFrictionalContact, resolve_articulated_frictional_contacts,
};
use crate::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, LinkContactPoint,
    unit_normal,
};
use fs_contact::ballistic::{
    BallisticContactError, BallisticContactInput, first_ballistic_contact,
};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_ga::{Twist, Vec3};

/// One infinite prescribed plane, with its position at the window start.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarrierPlane {
    /// Stable unique contact identity; discovery and simultaneous solve use this order.
    pub id: u64,
    /// World point on the plane at time zero of THIS call, in metres.
    pub point_world: Vec3,
    /// Constant prescribed translation, including tangential conveyor motion, in m/s.
    pub velocity_world: Vec3,
    /// Nonzero normal into the permitted half-space, normalized at admission.
    pub normal_world: Vec3,
    /// Nonparallel tangent hint for the circular Coulomb law.
    pub tangent_world: Vec3,
    /// Finite nonnegative ideal Coulomb coefficient.
    pub friction_coefficient: f64,
    /// Newton restitution in [0, 1].
    pub restitution: f64,
}

/// Bounded work and geometric resolution for one complete carrier window.
#[derive(Debug, Clone, Copy)]
pub struct CarrierPlaneConfig {
    /// Smooth accuracy, original impact controls, duration and momentum origin.
    /// `maximum_smooth_intervals` counts EVERY attempted native smooth step,
    /// including the bounded floating-time retries.
    pub trajectory: ArticulatedTrajectoryConfig,
    /// Positive spatial activation resolution, in metres. Nearby planes enter
    /// one simultaneous solve; no negative gap is admitted at any resolution.
    pub activation_tolerance: f64,
    /// At most 64 events. Zero permits only event-free windows.
    pub maximum_events: usize,
    /// Total plane discovery/path/witness observations, in 1..=1_000_000.
    pub maximum_geometry_queries: usize,
}

/// Native dynamics evidence plus complete analytic plane-path observations.
#[derive(Debug, Clone)]
pub struct CarrierPlaneResult {
    /// Actual original-tree states, impulses, work, impact loss and balance defects.
    pub trajectory: ArticulatedTrajectoryResult,
    /// All attempted native smooth steps, including rejected time neighbours.
    pub smooth_steps_attempted: usize,
    /// Number of charged plane observations, including rejected trials.
    pub geometry_queries: usize,
    /// Smallest evaluated sphere/plane gap over all accepted quadratic paths,
    /// actual native endpoints and event material witnesses, in metres.
    /// This floating analytic observation is not an interval enclosure.
    pub minimum_gap: f64,
}

/// Atomic refusal; no partially advanced state is returned.
#[derive(Debug)]
pub enum CarrierPlaneError {
    /// Model geometry, loads, controls or finite-value admission refused.
    Invalid(&'static str),
    /// Original smooth/contact/balance/cancellation refusal, unchanged.
    Trajectory(ArticulatedTrajectoryError),
    /// The scalar analytic geometry owner could not resolve a plane encounter.
    Ballistic {
        /// Plane identity.
        plane_id: u64,
        /// Current window time, in seconds.
        time: f64,
        /// Original scalar geometry refusal.
        error: BallisticContactError,
    },
    /// Complete path or actual material geometry could not be made admissible
    /// within 64 outward floating-time neighbours; no point is projected.
    UnresolvedGeometry {
        /// Attempted event/step time, in seconds.
        time: f64,
    },
    /// Contact requires sustained support, or its approach cannot be resolved
    /// at the original solver's velocity resolution. No compliant substitute.
    UnresolvedContact {
        /// Plane requiring a different contact lane or smaller resolution.
        plane_id: u64,
        /// Current time, in seconds.
        time: f64,
    },
    /// A second event at the same representable time would require resolving
    /// chatter or revisiting an already simultaneous contact set.
    NoEventProgress,
    /// The caller's attempted-smooth-step bound was exhausted.
    SmoothBudget,
    /// The caller's event bound was exhausted; includes bounded Zeno refusal.
    EventBudget,
    /// The caller's geometry-observation bound was exhausted.
    GeometryBudget,
}
impl core::fmt::Display for CarrierPlaneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spherical articulated carrier: {self:?}")
    }
}
impl std::error::Error for CarrierPlaneError {}
impl From<ArticulatedTrajectoryError> for CarrierPlaneError {
    fn from(value: ArticulatedTrajectoryError) -> Self {
        Self::Trajectory(value)
    }
}

/// Advance one spherical articulated carrier through automatically found impacts.
///
/// The declared sphere is attached to root link zero and encloses every internal
/// link. Every link COM and every zero-joint translation must be exactly zero;
/// fixed and revolute joints alone are admitted. The inertia second moments must
/// be consistent with the declared enclosure, but inertias cannot prove an
/// unknown shape is enclosed. The caller owns that geometric model declaration.
/// External link forces must be exactly zero. Finite body torques and internal
/// actuator efforts are allowed; uniform gravity supplies all COM acceleration.
///
/// Each step checks every plane's whole quadratic relative path, including an
/// interior minimum, plus the unchanged native endpoint. Actual root-link radial
/// material points receive one simultaneous original Coulomb solve at each
/// encounter. A root rounded into overlap is retried at up to 64 earlier floating
/// times; no pose, contact arm or velocity is projected. Positive activation
/// gaps permit early impulses only within the declared spatial resolution.
///
/// At most 21 planes, 4096 links and 4096 attempted smooth steps are admitted.
/// Cancellation, a late failure and any exhausted budget publish no state.
/// Returned plane positions for a subsequent window must be advanced by the
/// caller to that window's start. This operation covers this declared spherical
/// surface against these planes; it does not cover other bodies, protruding
/// internal geometry, rotating planes, sustained support, rolling constraints,
/// interval-certified CCD or exact rotational trajectory error.
pub fn advance_carrier_planes(
    system: ArticulatedTrajectorySystem<'_>,
    radius: f64,
    planes: &[CarrierPlane],
    config: CarrierPlaneConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CarrierPlaneResult, CarrierPlaneError> {
    poll(&mut cancelled)?;
    let planes = admit(system, radius, planes, config)?;
    let systems = [system];
    let mut state = system.state.clone();
    let origin = config.trajectory.angular_origin_world;
    let before = momentum(
        &systems,
        core::slice::from_ref(&state),
        origin,
        0.0,
        &mut cancelled,
    )?;
    let mut impulses = Impulses::default();
    let (mut actuator_work, mut external_work, mut loss, mut smooth_defect) = (0.0, 0.0, 0.0, 0.0);
    let mut events = Vec::new();
    let mut work = Work {
        smooth: 0,
        geometry: 0,
        accepted: 0,
        minimum_gap: f64::INFINITY,
    };
    let mut time = 0.0;
    while time < config.trajectory.duration {
        poll(&mut cancelled)?;
        let horizon =
            (config.trajectory.duration - time).min(config.trajectory.maximum_smooth_step);
        let velocity = world_velocity(&state)?;
        let mut encounter: Option<(f64, u64)> = None;
        for plane in &planes {
            work.query(config)?;
            let input = path_input(
                plane,
                &state,
                velocity,
                radius,
                time,
                horizon,
                system.loads.gravity_world,
            )?;
            let hit =
                first_ballistic_contact(input).map_err(|error| CarrierPlaneError::Ballistic {
                    plane_id: plane.id,
                    time,
                    error,
                })?;
            work.minimum_gap = work.minimum_gap.min(input.gap);
            if input.gap <= config.activation_tolerance
                && input.relative_normal_speed.abs() <= config.trajectory.contact.velocity_tolerance
                && input.relative_normal_acceleration < 0.0
            {
                return Err(CarrierPlaneError::UnresolvedContact {
                    plane_id: plane.id,
                    time,
                });
            }
            if let Some(hit) = hit
                && hit.relative_normal_speed < 0.0
                && encounter.is_none_or(|(earliest, _)| hit.time_s < earliest)
            {
                encounter = Some((hit.time_s, plane.id));
            }
        }
        let mut duration = encounter.map_or(horizon, |(at, _)| at);
        let mut accepted = None;
        for _ in 0..=64 {
            poll(&mut cancelled)?;
            let end = time + duration;
            if !end.is_finite() || (duration > 0.0 && end <= time) {
                return Err(CarrierPlaneError::UnresolvedGeometry { time: end });
            }
            let step = if duration > 0.0 {
                if work.smooth == config.trajectory.maximum_smooth_intervals {
                    return Err(CarrierPlaneError::SmoothBudget);
                }
                work.smooth += 1;
                Some(
                    step_free_floating(
                        system.model,
                        &state,
                        system.loads,
                        duration,
                        &mut cancelled,
                    )
                    .map_err(|error| smooth_error(0, end, error))?,
                )
            } else {
                None
            };
            let candidate = step.as_ref().map_or(&state, |step| &step.state_after);
            let (valid, minimum, contacts) = observe(
                system,
                &state,
                candidate,
                &planes,
                radius,
                time,
                duration,
                encounter.is_some(),
                config,
                &mut work,
                &mut cancelled,
            )?;
            if valid {
                accepted = Some((step, end, minimum, contacts));
                break;
            }
            if duration == 0.0 {
                break;
            }
            duration = duration.next_down();
        }
        let (step, end, minimum, contacts) =
            accepted.ok_or(CarrierPlaneError::UnresolvedGeometry {
                time: time + duration,
            })?;
        if let Some(step) = step {
            impulses.add(
                step.external_linear_impulse_world,
                step.external_angular_impulse_world,
                origin,
            );
            actuator_work += step.actuator_work;
            external_work += step.external_work;
            smooth_defect += step.kinetic_work_residual;
            state = step.state_after;
            work.accepted += 1;
        }
        time = end;
        work.minimum_gap = work.minimum_gap.min(minimum);
        if let Some((_, plane_id)) = encounter {
            if events.len() == config.maximum_events {
                return Err(CarrierPlaneError::EventBudget);
            }
            if events
                .last()
                .is_some_and(|event: &ArticulatedTrajectoryEventResult| event.time == time)
            {
                return Err(CarrierPlaneError::NoEventProgress);
            }
            if !contacts.iter().any(|contact| contact.id == plane_id) {
                return Err(CarrierPlaneError::UnresolvedGeometry { time });
            }
            let plane = planes
                .iter()
                .find(|plane| plane.id == plane_id)
                .expect("admitted plane");
            if plane_normal(plane)?.dot(world_velocity(&state)? - plane.velocity_world)
                >= -config.trajectory.contact.velocity_tolerance
            {
                return Err(CarrierPlaneError::UnresolvedContact { plane_id, time });
            }
            let impact = resolve_articulated_frictional_contacts(
                &[ArticulatedContactSystem {
                    model: system.model,
                    base: state.base,
                    joint_position: &state.joint_position,
                    joint_velocity: &state.joint_velocity,
                }],
                &contacts,
                config.trajectory.contact,
                &mut cancelled,
            )
            .map_err(|error| {
                if matches!(error, ArticulatedContactError::Cancelled) {
                    ArticulatedTrajectoryError::Cancelled
                } else {
                    ArticulatedTrajectoryError::Contact {
                        event_index: events.len(),
                        time,
                        error,
                    }
                }
            })?;
            let response = &impact.responses[0];
            impulses.add(
                response.external_linear_impulse_world,
                response.external_angular_impulse_world,
                origin,
            );
            state.base = response.base_after;
            state
                .joint_velocity
                .clone_from(&response.joint_velocity_after);
            external_work += impact.plane_work;
            loss += impact.contact.dissipated_energy;
            events.push(ArticulatedTrajectoryEventResult { time, impact });
        }
    }
    let after = momentum(
        &systems,
        core::slice::from_ref(&state),
        origin,
        time,
        &mut cancelled,
    )?;
    let energy_residual =
        after.kinetic_energy - before.kinetic_energy - actuator_work - external_work + loss;
    let linear_residual = after.linear_world - before.linear_world - impulses.linear;
    let angular_residual =
        after.angular_about_origin_world - before.angular_about_origin_world - impulses.angular;
    if ![
        actuator_work,
        external_work,
        loss,
        smooth_defect,
        energy_residual,
        work.minimum_gap,
    ]
    .iter()
    .all(|value| value.is_finite())
        || ![
            impulses.linear,
            impulses.angular,
            linear_residual,
            angular_residual,
        ]
        .into_iter()
        .all(finite_vector)
    {
        return Err(CarrierPlaneError::Invalid("nonfinite carrier balance"));
    }
    poll(&mut cancelled)?;
    Ok(CarrierPlaneResult {
        trajectory: ArticulatedTrajectoryResult {
            states_after: vec![state],
            events,
            smooth_intervals: work.accepted,
            angular_origin_world: origin,
            momentum_before: before,
            momentum_after: after,
            external_linear_impulse_world: impulses.linear,
            external_angular_impulse_about_origin_world: impulses.angular,
            actuator_work,
            external_work,
            contact_dissipated_energy: loss,
            smooth_kinetic_work_residual: smooth_defect,
            kinetic_work_residual: energy_residual,
            linear_momentum_residual: linear_residual,
            angular_momentum_residual: angular_residual,
        },
        smooth_steps_attempted: work.smooth,
        geometry_queries: work.geometry,
        minimum_gap: work.minimum_gap,
    })
}

fn admit(
    system: ArticulatedTrajectorySystem<'_>,
    radius: f64,
    planes: &[CarrierPlane],
    config: CarrierPlaneConfig,
) -> Result<Vec<CarrierPlane>, CarrierPlaneError> {
    admission(&[system], &[], config.trajectory)?;
    if !radius.is_finite()
        || radius <= 0.0
        || planes.is_empty()
        || planes.len() > 21
        || !config.activation_tolerance.is_finite()
        || config.activation_tolerance <= 0.0
        || config.maximum_events > 64
        || !(1..=1_000_000).contains(&config.maximum_geometry_queries)
        || !config.trajectory.contact.velocity_tolerance.is_finite()
        || config.trajectory.contact.velocity_tolerance <= 0.0
    {
        return Err(CarrierPlaneError::Invalid(
            "invalid carrier, plane count or work bounds",
        ));
    }
    for link in system.model.links() {
        if link.inertia().center_of_mass() != Vec3::default()
            || link.parent_from_child_zero().translation() != Vec3::default()
            || !matches!(link.joint().kind(), JointKind::Fixed | JointKind::Revolute)
        {
            return Err(CarrierPlaneError::Invalid(
                "carrier links must share the base COM with fixed/revolute joints",
            ));
        }
        let inertia = link.inertia();
        let matrix = inertia.inertia_com().m;
        // sqrt(I / m) avoids underflow in I / m before comparison with a
        // tiny radius. The normalized trace cannot overflow.
        let scale = matrix[0].max(matrix[4]).max(matrix[8]);
        let rms_radius = (scale.sqrt() / inertia.mass().sqrt())
            * (0.5 * (matrix[0] / scale + matrix[4] / scale + matrix[8] / scale)).sqrt();
        if !rms_radius.is_finite() || rms_radius > radius {
            return Err(CarrierPlaneError::Invalid(
                "link mass moments cannot fit the declared spherical enclosure",
            ));
        }
    }
    if system.loads.external_body_wrench.len() != system.model.link_count()
        || system
            .loads
            .external_body_wrench
            .iter()
            .any(|wrench| wrench.force != Vec3::default())
    {
        return Err(CarrierPlaneError::Invalid(
            "only gravity may supply carrier translation",
        ));
    }
    let mut ordered = planes.to_vec();
    ordered.sort_unstable_by_key(|plane| plane.id);
    if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(CarrierPlaneError::Invalid(
            "duplicate carrier plane identity",
        ));
    }
    for plane in &ordered {
        let normal = plane_normal(plane)?;
        let tangent = unit_normal(plane.tangent_world)
            .map_err(|_| CarrierPlaneError::Invalid("invalid plane tangent"))?;
        unit_normal(normal.cross(tangent))
            .map_err(|_| CarrierPlaneError::Invalid("parallel plane tangent"))?;
        if !finite_vector(plane.point_world)
            || !finite_vector(plane.velocity_world)
            || !plane.friction_coefficient.is_finite()
            || plane.friction_coefficient < 0.0
            || !plane.restitution.is_finite()
            || !(0.0..=1.0).contains(&plane.restitution)
        {
            return Err(CarrierPlaneError::Invalid(
                "invalid plane motion or friction/restitution law",
            ));
        }
    }
    Ok(ordered)
}

fn world_velocity(state: &FreeFloatingMotionState) -> Result<Vec3, CarrierPlaneError> {
    state
        .base
        .world_from_base
        .rotation()
        .rotate(state.base.twist_body.linear)
        .map_err(|_| CarrierPlaneError::Invalid("nonfinite carrier centre velocity"))
}
fn plane_normal(plane: &CarrierPlane) -> Result<Vec3, CarrierPlaneError> {
    // Normalize the original hint exactly as the native impact owner does;
    // renormalizing an already rounded unit vector could change a tiny gap.
    unit_normal(plane.normal_world).map_err(|_| CarrierPlaneError::Invalid("invalid plane normal"))
}
fn plane_position(plane: &CarrierPlane, time: f64) -> Result<Vec3, CarrierPlaneError> {
    let value = plane.point_world + plane.velocity_world.scale(time);
    if !finite_vector(value) {
        return Err(CarrierPlaneError::Invalid(
            "nonfinite translating plane position",
        ));
    }
    Ok(value)
}
fn path_input(
    plane: &CarrierPlane,
    state: &FreeFloatingMotionState,
    velocity: Vec3,
    radius: f64,
    time: f64,
    duration: f64,
    gravity: Vec3,
) -> Result<BallisticContactInput, CarrierPlaneError> {
    let normal = plane_normal(plane)?;
    Ok(BallisticContactInput {
        gap: normal.dot(state.base.world_from_base.translation() - plane_position(plane, time)?)
            - radius,
        relative_normal_speed: normal.dot(velocity - plane.velocity_world),
        relative_normal_acceleration: normal.dot(gravity),
        duration_s: duration,
    })
}

struct Work {
    smooth: usize,
    geometry: usize,
    accepted: usize,
    minimum_gap: f64,
}
impl Work {
    fn query(&mut self, config: CarrierPlaneConfig) -> Result<(), CarrierPlaneError> {
        if self.geometry == config.maximum_geometry_queries {
            return Err(CarrierPlaneError::GeometryBudget);
        }
        self.geometry += 1;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn observe(
    system: ArticulatedTrajectorySystem<'_>,
    before: &FreeFloatingMotionState,
    after: &FreeFloatingMotionState,
    planes: &[CarrierPlane],
    radius: f64,
    time: f64,
    duration: f64,
    event: bool,
    config: CarrierPlaneConfig,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(bool, f64, Vec<ArticulatedFrictionalContact>), CarrierPlaneError> {
    let velocity = world_velocity(before)?;
    let root_pose = if event {
        poll(cancelled)?;
        Some(
            forward_kinematics(
                system.model,
                BaseState::prescribed(
                    after.base.world_from_base,
                    after.base.twist_body,
                    Twist::zero(),
                ),
                &after.joint_position,
                &after.joint_velocity,
            )
            .map_err(|error| smooth_error(0, time + duration, error.into()))?
            .world_from_link[0],
        )
    } else {
        None
    };
    let mut contacts = Vec::new();
    let mut minimum = f64::INFINITY;
    for plane in planes {
        poll(cancelled)?;
        work.query(config)?;
        let input = path_input(
            plane,
            before,
            velocity,
            radius,
            time,
            duration,
            system.loads.gravity_world,
        )?;
        let gap_at = |t: f64| {
            (0.5 * input.relative_normal_acceleration)
                .mul_add(t, input.relative_normal_speed)
                .mul_add(t, input.gap)
        };
        let mut path_gap = input.gap.min(gap_at(duration));
        if input.relative_normal_acceleration > 0.0 {
            let extremum = -input.relative_normal_speed / input.relative_normal_acceleration;
            if extremum > 0.0 && extremum < duration {
                path_gap = path_gap.min(gap_at(extremum));
            }
        }
        let point_world = plane_position(plane, time + duration)?;
        let normal = plane_normal(plane)?;
        let gap = normal.dot(after.base.world_from_base.translation() - point_world) - radius;
        if !path_gap.is_finite() || !gap.is_finite() {
            return Err(CarrierPlaneError::Invalid(
                "nonfinite complete carrier path gap",
            ));
        }
        minimum = minimum.min(path_gap).min(gap);
        if minimum < 0.0 {
            return Ok((false, minimum, contacts));
        }
        if let Some(root_pose) = root_pose
            && gap <= config.activation_tolerance
        {
            let point_body = root_pose
                .rotation()
                .inverse()
                .rotate(normal.scale(-radius))
                .map_err(|_| CarrierPlaneError::Invalid("nonfinite carrier surface arm"))?;
            let material_point = root_pose
                .transform_point(point_body)
                .map_err(|_| CarrierPlaneError::Invalid("nonfinite carrier surface point"))?;
            let witness_gap = normal.dot(material_point - point_world);
            if !witness_gap.is_finite() {
                return Err(CarrierPlaneError::Invalid("nonfinite carrier material gap"));
            }
            minimum = minimum.min(witness_gap);
            if witness_gap < 0.0 {
                return Ok((false, minimum, contacts));
            }
            if witness_gap <= config.activation_tolerance {
                contacts.push(ArticulatedFrictionalContact {
                    id: plane.id,
                    point: LinkContactPoint {
                        system_index: 0,
                        link_index: 0,
                        point_body,
                    },
                    other: ContactCounterpart::PrescribedPlane {
                        point_world,
                        velocity_world: plane.velocity_world,
                    },
                    normal_world: plane.normal_world,
                    tangent_world: plane.tangent_world,
                    activation_tolerance: config.activation_tolerance,
                    friction_coefficient: plane.friction_coefficient,
                    normal_law: NormalVelocityLaw::NewtonRestitution {
                        coefficient: plane.restitution,
                    },
                });
            }
        }
    }
    poll(cancelled)?;
    Ok((true, minimum, contacts))
}
