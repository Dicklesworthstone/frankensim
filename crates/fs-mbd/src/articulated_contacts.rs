//! Simultaneous point impacts over free articulated systems.
//!
//! One contact-space inverse mass is assembled from the existing ABA owner.
//! All impulses are solved together and accumulated before any physical state
//! is evaluated. The caller supplies the event configuration and contact set;
//! this is not collision detection, position repair, or a contact time stepper.
//! The normal-only, [`friction`] and [`tangential`] entries share native mechanics.

/// Coupled normal and circular Coulomb friction on original articulated states.
pub mod friction;

/// Stateful, work-conjugate tangential contact on original articulated states.
pub mod tangential;

/// Continuous free trajectories joined by caller-located simultaneous impacts.
pub mod trajectory;

use core::fmt;
use fs_contact::frictional_multicontact::FrictionalContactError;
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_ga::{GaError, Twist, Vec3, Wrench};
use fs_tribo::TriboError;

use crate::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState,
    FreeFloatingImpulseResponse, Kinematics, apply_free_floating_impulses, forward_kinematics,
    free_floating_impulse_delta,
};

/// Borrowed state of one independently free articulated tree.
#[derive(Debug, Clone, Copy)]
pub struct ArticulatedContactSystem<'a> {
    /// Validated articulated model.
    pub model: &'a ArticulatedModel,
    /// Base pose and body-coordinate twist at the event.
    pub base: FreeFloatingBaseState,
    /// Joint positions, unchanged by impact.
    pub joint_position: &'a [f64],
    /// Joint velocities before impact.
    pub joint_velocity: &'a [f64],
}

/// A material point on one link of one supplied system.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinkContactPoint {
    /// Index in the supplied system slice.
    pub system_index: usize,
    /// Link index within that system.
    pub link_index: usize,
    /// Link-coordinate point relative to its origin, in metres.
    pub point_body: Vec3,
}

/// Opposite side of a point contact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContactCounterpart {
    /// Dynamic point receiving the opposite impulse, including within a tree.
    Body(LinkContactPoint),
    /// Prescribed plane; its actuator can supply or absorb impact work.
    PrescribedPlane {
        /// Any world-coordinate point on the plane, in metres.
        point_world: Vec3,
        /// Prescribed plane material velocity at contact, in m/s.
        velocity_world: Vec3,
    },
}

/// One normal row in the simultaneous impact problem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArticulatedNormalContact {
    /// Point receiving the positive normal impulse.
    pub point: LinkContactPoint,
    /// Opposite dynamic point or prescribed plane.
    pub other: ContactCounterpart,
    /// Nonzero outward normal from `other` to `point`, normalized here.
    pub normal_world: Vec3,
    /// Maximum nonnegative gap admitted as active, in metres.
    pub activation_tolerance: f64,
    /// Requested restitution in [0, 1]; never silently changed on refusal.
    pub restitution: f64,
}

/// Complete candidate states are published only after global admission.
#[derive(Debug, Clone)]
pub struct ArticulatedContactResult {
    /// Updated states and physical momenta in input-system order.
    pub responses: Vec<FreeFloatingImpulseResponse>,
    /// Simultaneous scalar impulses in input-contact order, in N s.
    pub normal_impulses: Vec<f64>,
    /// Relative speeds measured from the actual returned states, in m/s.
    pub relative_normal_speeds_after: Vec<f64>,
    /// Full row-major contact inverse mass J M^-1 J^T, in inverse kilograms.
    pub delassus: Vec<f64>,
    /// Contact-law result, including predicted speeds and solver residuals.
    pub normal_law: MultiImpactResult,
    /// Signed work from all prescribed planes, in joules.
    pub plane_work: f64,
    /// Globally admitted impact loss, in joules; no heat partition is implied.
    pub dissipated_energy: f64,
    /// Actual kinetic-energy change minus plane work plus admitted loss, in J.
    pub energy_balance_residual: f64,
    /// Total linear-momentum change minus prescribed-plane impulse, in N s.
    pub linear_momentum_residual: Vec3,
    /// Total angular-momentum change minus plane angular impulse, in N m s.
    pub angular_momentum_residual: Vec3,
}

/// Atomic refusal of the complete multi-system event.
#[derive(Debug)]
pub enum ArticulatedContactError {
    /// Invalid geometry, indexing, or a request above the bounded event size.
    Invalid(&'static str),
    /// Dynamics or physical endpoint admission failed for one system.
    Dynamics {
        /// Input system whose evaluation refused.
        system_index: usize,
        /// Original dynamics refusal.
        error: ArticulatedError,
    },
    /// Geometry operation refused.
    Geometry(GaError),
    /// Shared simultaneous-contact law or generic solver refused.
    ContactLaw(MultiImpactError),
    /// Shared coupled normal/Coulomb law refused.
    FrictionLaw(FrictionalContactError),
    /// Shared stateful tangential constitutive law refused.
    TangentialLaw(TriboError),
    /// Shared two-axis constitutive/body-compliance solve refused.
    TangentialMatrix(fs_tribo::tangential::matrix::MatrixTangentialError),
    /// Redundant impulses do not determine work between independently moving
    /// supports. A common support translation has determined total work.
    UnresolvedSupportWork,
    /// Dynamic witnesses do not share a normal line of action.
    NonCollinearContact {
        /// Contact whose paired forces would create an artificial couple.
        index: usize,
    },
    /// Actual endpoint failed the declared speed or global balance allowance.
    Unresolved(&'static str),
    /// Cancellation was observed before publishing any candidate state.
    Cancelled,
}

impl fmt::Display for ArticulatedContactError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(reason) | Self::Unresolved(reason) => f.write_str(reason),
            Self::Dynamics {
                system_index,
                error,
            } => write!(f, "impact system {system_index}: {error}"),
            Self::Geometry(error) => write!(f, "impact geometry: {error}"),
            Self::ContactLaw(error) => write!(f, "simultaneous contact: {error}"),
            Self::FrictionLaw(error) => write!(f, "simultaneous frictional contact: {error}"),
            Self::TangentialLaw(error) => write!(f, "stateful tangential contact: {error}"),
            Self::TangentialMatrix(error) => write!(f, "two-axis tangential contact: {error}"),
            Self::UnresolvedSupportWork => {
                f.write_str("redundant impulses do not determine moving-support work")
            }
            Self::NonCollinearContact { index } => {
                write!(f, "contact {index} has no common normal line of action")
            }
            Self::Cancelled => f.write_str("articulated contact event cancelled"),
        }
    }
}
impl std::error::Error for ArticulatedContactError {}
impl From<GaError> for ArticulatedContactError {
    fn from(error: GaError) -> Self {
        Self::Geometry(error)
    }
}
impl From<MultiImpactError> for ArticulatedContactError {
    fn from(error: MultiImpactError) -> Self {
        if matches!(error, MultiImpactError::Cancelled) {
            Self::Cancelled
        } else {
            Self::ContactLaw(error)
        }
    }
}
impl From<FrictionalContactError> for ArticulatedContactError {
    fn from(error: FrictionalContactError) -> Self {
        if matches!(error, FrictionalContactError::Cancelled) {
            Self::Cancelled
        } else {
            Self::FrictionLaw(error)
        }
    }
}
impl From<TriboError> for ArticulatedContactError {
    fn from(error: TriboError) -> Self {
        Self::TangentialLaw(error)
    }
}
impl From<fs_tribo::tangential::matrix::MatrixTangentialError> for ArticulatedContactError {
    fn from(error: fs_tribo::tangential::matrix::MatrixTangentialError) -> Self {
        if matches!(
            error,
            fs_tribo::tangential::matrix::MatrixTangentialError::Cancelled
        ) {
            Self::Cancelled
        } else {
            Self::TangentialMatrix(error)
        }
    }
}

#[derive(Clone, Copy)]
struct Term {
    point: LinkContactPoint,
    direction: Wrench,
}
struct Prepared {
    first: Term,
    second: Option<Term>,
    axis_world: Vec3,
    world_point: Vec3,
    plane_speed: f64,
    /// Only active normal rows carry a nonpenetration speed bound.
    minimum_speed: Option<f64>,
}
impl Prepared {
    fn terms(&self) -> impl Iterator<Item = Term> + '_ {
        core::iter::once(self.first).chain(self.second)
    }
    fn velocity(&self, kinematics: &[Kinematics]) -> f64 {
        self.terms()
            .map(|term| {
                term.direction
                    .pairing(kinematics[term.point.system_index].body_twist[term.point.link_index])
            })
            .sum::<f64>()
            - self.plane_speed
    }
}

/// Resolve up to 64 supplied normal contacts over at most 128 articulated trees.
///
/// One global contact-space solve precedes one accumulated physical impulse
/// application per system. Each inverse-mass probe visits only its one or two
/// involved systems. Body-pair impulses are equal and opposite on one normal
/// line, preserving global angular as well as linear momentum. Prescribed
/// planes contribute explicit external impulse and signed work.
///
/// All gaps must be nonnegative. Singular/redundant impulse systems, invalid
/// restitution candidates, unresolved endpoint balances, and cancellation
/// refuse atomically. No friction, contact detection, penetration repair,
/// continuous forcing, or timestep advancement is performed.
pub fn resolve_articulated_contacts(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedNormalContact],
    config: MultiImpactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedContactResult, ArticulatedContactError> {
    if contacts.len() > 64 || contacts.len() > config.solver.max_variables || systems.len() > 128 {
        return Err(ArticulatedContactError::Invalid(
            "articulated contact event exceeds its size cap",
        ));
    }
    let mut workspace = EventWorkspace::new(systems, &mut cancelled)?;
    let mut prepared = Vec::with_capacity(contacts.len());
    let mut laws = Vec::with_capacity(contacts.len());
    for (index, contact) in contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        let normal = unit_normal(contact.normal_world)?;
        let (first, world_point) = prepare_point(contact.point, normal, &workspace.before)?;
        let (second, other_point, plane_speed) = match contact.other {
            ContactCounterpart::Body(point) => {
                let (term, position) = prepare_point(point, normal.scale(-1.0), &workspace.before)?;
                (Some(term), position, 0.0)
            }
            ContactCounterpart::PrescribedPlane {
                point_world,
                velocity_world,
            } => {
                vector(point_world)?;
                vector(velocity_world)?;
                (None, point_world, finite(normal.dot(velocity_world))?)
            }
        };
        let separation = world_point - other_point;
        let gap = finite(normal.dot(separation))?;
        if second.is_some() {
            let tangent = separation - normal.scale(gap);
            let scale = magnitude(world_point)
                .max(magnitude(other_point))
                .max(magnitude(separation));
            if magnitude(tangent) > 256.0 * f64::EPSILON * scale {
                return Err(ArticulatedContactError::NonCollinearContact { index });
            }
        }
        let item = Prepared {
            first,
            second,
            axis_world: normal,
            world_point,
            plane_speed,
            minimum_speed: (gap <= contact.activation_tolerance).then_some(0.0),
        };
        laws.push(NormalContactImpact {
            gap,
            activation_tolerance: contact.activation_tolerance,
            relative_normal_speed: finite(item.velocity(&workspace.before))?,
            restitution: contact.restitution,
        });
        prepared.push(item);
    }
    let delassus = workspace.delassus(systems, &prepared, &mut cancelled)?;
    let velocity_tolerance = config.solver.gradient_tolerance;
    let energy_tolerance = config.energy_tolerance;
    let law = solve_normal_impacts(&laws, &delassus, config, &mut cancelled)?;
    let mechanics = workspace.apply(
        systems,
        &prepared,
        EventResponse {
            impulses: &law.normal_impulses,
            predicted_speeds: &law.relative_normal_speeds_after,
            stored_energy_change: 0.0,
            dissipated_energy: law.dissipated_energy,
            velocity_tolerance,
            energy_tolerance,
        },
        &mut cancelled,
    )?;
    Ok(ArticulatedContactResult {
        responses: mechanics.responses,
        normal_impulses: law.normal_impulses.clone(),
        relative_normal_speeds_after: mechanics.speeds_after,
        delassus,
        plane_work: mechanics.plane_work,
        dissipated_energy: law.dissipated_energy,
        normal_law: law,
        energy_balance_residual: mechanics.energy_balance_residual,
        linear_momentum_residual: mechanics.linear_momentum_residual,
        angular_momentum_residual: mechanics.angular_momentum_residual,
    })
}

// The normal and frictional laws share one original articulated inverse-mass
// action and one independent physical endpoint audit. A row can be normal or
// tangent; no surrogate link mass or separate friction impulse update exists.
struct EventWorkspace {
    before: Vec<Kinematics>,
    loads: Vec<Vec<Wrench>>,
    joint_impulses: Vec<Vec<f64>>,
}
struct EventResponse<'a> {
    impulses: &'a [f64],
    predicted_speeds: &'a [f64],
    // Zero for impact-only laws. Stateful shear can release this storage into
    // positive body kinetic work while total physical dissipation stays >= 0.
    stored_energy_change: f64,
    dissipated_energy: f64,
    velocity_tolerance: f64,
    energy_tolerance: f64,
}
struct EventResult {
    responses: Vec<FreeFloatingImpulseResponse>,
    speeds_after: Vec<f64>,
    plane_work: f64,
    kinetic_energy_change: f64,
    energy_balance_residual: f64,
    linear_momentum_residual: Vec3,
    angular_momentum_residual: Vec3,
    endpoint_velocity_residual: f64,
}
impl EventWorkspace {
    fn new(
        systems: &[ArticulatedContactSystem<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ArticulatedContactError> {
        let mut before = Vec::with_capacity(systems.len());
        let mut loads = Vec::with_capacity(systems.len());
        let mut joint_impulses = Vec::with_capacity(systems.len());
        for (index, system) in systems.iter().enumerate() {
            poll(cancelled)?;
            before.push(dynamics(
                index,
                forward_kinematics(
                    system.model,
                    BaseState::prescribed(
                        system.base.world_from_base,
                        system.base.twist_body,
                        Twist::zero(),
                    ),
                    system.joint_position,
                    system.joint_velocity,
                ),
            )?);
            loads.push(vec![Wrench::default(); system.model.link_count()]);
            joint_impulses.push(vec![0.0; system.model.dof_count()]);
        }
        Ok(Self {
            before,
            loads,
            joint_impulses,
        })
    }

    fn delassus(
        &mut self,
        systems: &[ArticulatedContactSystem<'_>],
        prepared: &[Prepared],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, ArticulatedContactError> {
        let count = prepared.len();
        let mut delassus = vec![0.0; count * count];
        for (column, contact) in prepared.iter().enumerate() {
            for (side, term) in contact.terms().enumerate() {
                let system_index = term.point.system_index;
                if side == 1 && system_index == contact.first.point.system_index {
                    continue;
                }
                poll(cancelled)?;
                for term in contact
                    .terms()
                    .filter(|term| term.point.system_index == system_index)
                {
                    add_load(
                        &mut self.loads[system_index][term.point.link_index],
                        term.direction,
                        1.0,
                    )?;
                }
                let system = &systems[system_index];
                let delta = dynamics(
                    system_index,
                    free_floating_impulse_delta(
                        system.model,
                        system.joint_position,
                        &self.joint_impulses[system_index],
                        &self.loads[system_index],
                    ),
                )?;
                for (row, other) in prepared.iter().enumerate() {
                    for term in other
                        .terms()
                        .filter(|term| term.point.system_index == system_index)
                    {
                        delassus[row * count + column] = finite(
                            delassus[row * count + column]
                                + term.direction.pairing(delta.body[term.point.link_index]),
                        )?;
                    }
                }
                for term in contact
                    .terms()
                    .filter(|term| term.point.system_index == system_index)
                {
                    self.loads[system_index][term.point.link_index] = Wrench::default();
                }
            }
        }
        Ok(delassus)
    }

    fn apply(
        mut self,
        systems: &[ArticulatedContactSystem<'_>],
        prepared: &[Prepared],
        law: EventResponse<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EventResult, ArticulatedContactError> {
        let mut plane_work = 0.0;
        let mut external_linear = Vec3::default();
        let mut external_angular = Vec3::default();
        let mut linear_impulse_scale = 0.0;
        let mut angular_impulse_scale = 0.0;
        for (contact, &impulse) in prepared.iter().zip(law.impulses) {
            for term in contact.terms() {
                add_load(
                    &mut self.loads[term.point.system_index][term.point.link_index],
                    term.direction,
                    impulse,
                )?;
                // Retain the magnitudes before cancellation: equal/opposite
                // contacts and spin/orbital terms can sum to exactly zero.
                let pose =
                    self.before[term.point.system_index].world_from_link[term.point.link_index];
                let force = pose.rotation().rotate(term.direction.force)?.scale(impulse);
                let spin = pose
                    .rotation()
                    .rotate(term.direction.torque)?
                    .scale(impulse);
                let orbital = pose.translation().cross(force);
                linear_impulse_scale = finite(linear_impulse_scale + magnitude(force))?;
                angular_impulse_scale =
                    finite(angular_impulse_scale + magnitude(spin) + magnitude(orbital))?;
            }
            if contact.second.is_none() {
                let linear = contact.axis_world.scale(impulse);
                external_linear = external_linear + linear;
                external_angular = external_angular + contact.world_point.cross(linear);
                plane_work = finite(plane_work + impulse * contact.plane_speed)?;
            }
        }
        vector(external_linear)?;
        vector(external_angular)?;
        let mut responses = Vec::with_capacity(systems.len());
        let mut after = Vec::with_capacity(systems.len());
        let mut energy_change = 0.0;
        let mut energy_scale =
            plane_work.abs() + law.dissipated_energy.abs() + law.stored_energy_change.abs();
        let mut linear_change = Vec3::default();
        let mut angular_change = Vec3::default();
        let mut linear_scale = linear_impulse_scale + magnitude(external_linear);
        let mut angular_scale = angular_impulse_scale + magnitude(external_angular);
        for (index, system) in systems.iter().enumerate() {
            poll(cancelled)?;
            let response = dynamics(
                index,
                apply_free_floating_impulses(
                    system.model,
                    system.base,
                    system.joint_position,
                    system.joint_velocity,
                    &self.joint_impulses[index],
                    &self.loads[index],
                ),
            )?;
            after.push(dynamics(
                index,
                forward_kinematics(
                    system.model,
                    BaseState::prescribed(
                        response.base_after.world_from_base,
                        response.base_after.twist_body,
                        Twist::zero(),
                    ),
                    system.joint_position,
                    &response.joint_velocity_after,
                ),
            )?);
            let initial = response.momentum_before;
            let final_state = response.momentum_after;
            energy_change =
                finite(energy_change + (final_state.kinetic_energy - initial.kinetic_energy))?;
            energy_scale = finite(
                energy_scale + initial.kinetic_energy.abs() + final_state.kinetic_energy.abs(),
            )?;
            linear_change =
                linear_change + (final_state.linear_momentum_world - initial.linear_momentum_world);
            angular_change = angular_change
                + (final_state.angular_momentum_world - initial.angular_momentum_world);
            linear_scale = finite(
                linear_scale
                    + magnitude(initial.linear_momentum_world)
                    + magnitude(final_state.linear_momentum_world),
            )?;
            angular_scale = finite(
                angular_scale
                    + magnitude(initial.angular_momentum_world)
                    + magnitude(final_state.angular_momentum_world),
            )?;
            // A free tree may have zero total momentum while individual links
            // carry large opposite linear/spin momenta. Their cancellation is
            // not a zero roundoff scale. Retain the original link contributions
            // at both endpoints, including the separate COM and world arms.
            for kinematics in [&self.before[index], &after[index]] {
                for (link_index, link) in system.model.links().iter().enumerate() {
                    if link_index % 64 == 0 {
                        poll(cancelled)?;
                    }
                    let pose = kinematics.world_from_link[link_index];
                    let velocity = kinematics.body_twist[link_index];
                    let body_momentum = dynamics(index, link.inertia().momentum(velocity))?;
                    let linear = pose.rotation().rotate(body_momentum.force)?;
                    let spin = pose
                        .rotation()
                        .rotate(link.inertia().inertia_com().apply(velocity.angular))?;
                    let com_arm = pose
                        .rotation()
                        .rotate(link.inertia().center_of_mass().cross(body_momentum.force))?;
                    let world_arm = pose.translation().cross(linear);
                    linear_scale = finite(linear_scale + magnitude(linear))?;
                    angular_scale = finite(
                        angular_scale + magnitude(spin) + magnitude(com_arm) + magnitude(world_arm),
                    )?;
                }
            }
            responses.push(response);
        }
        let roundoff = 4096.0 * f64::EPSILON;
        let mut speeds_after = Vec::with_capacity(prepared.len());
        let mut endpoint_velocity_residual = 0.0_f64;
        for (index, contact) in prepared.iter().enumerate() {
            let actual = finite(contact.velocity(&after))?;
            let predicted = law.predicted_speeds[index];
            let allowance = law.velocity_tolerance
                + roundoff
                    * contact
                        .velocity(&self.before)
                        .abs()
                        .max(predicted.abs())
                        .max(actual.abs());
            if (actual - predicted).abs() > allowance
                || contact
                    .minimum_speed
                    .is_some_and(|minimum| actual < minimum - allowance)
            {
                return Err(ArticulatedContactError::Unresolved(
                    "actual contact velocity failed admission",
                ));
            }
            endpoint_velocity_residual = endpoint_velocity_residual.max((actual - predicted).abs());
            speeds_after.push(actual);
        }
        let relative_work = finite(energy_change - plane_work)?;
        let stored_and_mechanical_work = finite(relative_work + law.stored_energy_change)?;
        let energy_balance_residual = finite(stored_and_mechanical_work + law.dissipated_energy)?;
        let energy_allowance = finite(law.energy_tolerance + roundoff * energy_scale)?;
        if stored_and_mechanical_work > energy_allowance
            || energy_balance_residual.abs() > energy_allowance
        {
            return Err(ArticulatedContactError::Unresolved(
                "actual global kinetic and contact energy failed admission",
            ));
        }
        let linear_momentum_residual = linear_change - external_linear;
        let angular_momentum_residual = angular_change - external_angular;
        vector(linear_momentum_residual)?;
        vector(angular_momentum_residual)?;
        if magnitude(linear_momentum_residual) > roundoff * linear_scale
            || magnitude(angular_momentum_residual) > roundoff * angular_scale
        {
            return Err(ArticulatedContactError::Unresolved(
                "actual global impact momentum failed admission",
            ));
        }
        poll(cancelled)?;
        Ok(EventResult {
            responses,
            speeds_after,
            plane_work,
            kinetic_energy_change: energy_change,
            energy_balance_residual,
            linear_momentum_residual,
            angular_momentum_residual,
            endpoint_velocity_residual,
        })
    }
}

fn prepare_point(
    point: LinkContactPoint,
    normal: Vec3,
    kinematics: &[Kinematics],
) -> Result<(Term, Vec3), ArticulatedContactError> {
    let system = kinematics
        .get(point.system_index)
        .ok_or(ArticulatedContactError::Invalid(
            "contact system index is out of range",
        ))?;
    let pose =
        system
            .world_from_link
            .get(point.link_index)
            .ok_or(ArticulatedContactError::Invalid(
                "contact link index is out of range",
            ))?;
    vector(point.point_body)?;
    let direction = pose.rotation().inverse().rotate(normal)?;
    Ok((
        Term {
            point,
            direction: Wrench::new(point.point_body.cross(direction), direction),
        },
        pose.transform_point(point.point_body)?,
    ))
}
fn add_load(
    load: &mut Wrench,
    direction: Wrench,
    magnitude: f64,
) -> Result<(), ArticulatedContactError> {
    load.torque = load.torque + direction.torque.scale(magnitude);
    load.force = load.force + direction.force.scale(magnitude);
    vector(load.torque)?;
    vector(load.force)
}
fn unit_normal(normal: Vec3) -> Result<Vec3, ArticulatedContactError> {
    vector(normal)?;
    let scale = magnitude(normal);
    if scale == 0.0 {
        return Err(ArticulatedContactError::Invalid(
            "contact normal must be nonzero",
        ));
    }
    let scaled = Vec3::new(normal.x / scale, normal.y / scale, normal.z / scale);
    Ok(scaled.scale(1.0 / scaled.norm()))
}
fn magnitude(value: Vec3) -> f64 {
    value.x.abs().max(value.y.abs()).max(value.z.abs())
}
fn vector(value: Vec3) -> Result<(), ArticulatedContactError> {
    finite(value.x)?;
    finite(value.y)?;
    finite(value.z)?;
    Ok(())
}
fn finite(value: f64) -> Result<f64, ArticulatedContactError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ArticulatedContactError::Invalid(
            "nonfinite impact arithmetic",
        ))
    }
}
fn dynamics<T>(
    system_index: usize,
    result: Result<T, ArticulatedError>,
) -> Result<T, ArticulatedContactError> {
    result.map_err(|error| ArticulatedContactError::Dynamics {
        system_index,
        error,
    })
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ArticulatedContactError> {
    if cancelled() {
        Err(ArticulatedContactError::Cancelled)
    } else {
        Ok(())
    }
}
