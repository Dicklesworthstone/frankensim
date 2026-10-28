//! One frictionless point impact between a free articulated tree and a plane.
//!
//! Geometry and point inertia come from the articulated owner; the scalar
//! restitution law comes from `fs-contact`. The caller supplies the impact
//! configuration. This module neither detects a collision time nor repairs
//! penetration, advances positions, or resolves simultaneous contacts.

use core::fmt;
use fs_contact::impact::{
    NormalImpactError, NormalImpactInput, NormalImpactResult, solve_normal_impact,
};
use fs_ga::{GaError, Twist, Vec3, Wrench};

use crate::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState,
    FreeFloatingImpulseResponse, apply_free_floating_impulses, forward_kinematics,
    free_floating_impulse_delta,
};

/// One link-fixed point and the prescribed plane it may strike.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArticulatedPointPlaneContact {
    /// Contacting link in model order.
    pub link_index: usize,
    /// Point relative to that link's origin, in link coordinates, in metres.
    pub point_body: Vec3,
    /// Any point on the plane, in world coordinates, in metres.
    pub plane_point_world: Vec3,
    /// Nonzero outward normal into the admissible half-space; normalized here.
    pub normal_world: Vec3,
    /// Prescribed material-point velocity of the plane at contact, in m/s.
    pub plane_velocity_world: Vec3,
    /// Nonnegative maximum admitted contact gap, in metres. Negative gaps
    /// always refuse; this tolerance never authorizes penetration repair.
    pub activation_tolerance: f64,
    /// Newton normal restitution coefficient in the closed interval [0, 1].
    pub restitution: f64,
}

/// Physical endpoint state and single-contact momentum/energy accounting.
#[derive(Debug, Clone, PartialEq)]
pub struct ArticulatedPlaneImpactResult {
    /// Atomic physical impulse response, including all base and joint recoil.
    pub response: FreeFloatingImpulseResponse,
    /// Result returned by the shared scalar normal-impact law.
    pub normal_law: NormalImpactResult,
    /// Signed initial point/plane gap, in metres; nonnegative on success.
    pub gap: f64,
    /// Normal point velocity change per unit impulse, in inverse kilograms.
    pub inverse_effective_mass: f64,
    /// Initial relative normal speed, in m/s; negative means approach.
    pub relative_normal_speed_before: f64,
    /// Relative normal speed measured from the actual returned endpoint.
    pub relative_normal_speed_after: f64,
    /// Unchanged contacting point in world coordinates, in metres.
    pub contact_point_world: Vec3,
    /// Work supplied by the prescribed moving plane, in joules. Positive work
    /// can increase body kinetic energy even with passive restitution.
    pub plane_work: f64,
    /// Kinetic-energy change minus plane work plus impact dissipation, in J.
    pub energy_balance_residual: f64,
}

/// Refusal leaves the borrowed articulated state unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum ArticulatedImpactError {
    /// Articulated state, inertia solve, or endpoint speed admission failed.
    Dynamics(ArticulatedError),
    /// The shared contact law rejected the gap or constitutive inputs.
    ContactLaw(NormalImpactError),
    /// The selected link is outside the model.
    InvalidLinkIndex {
        /// Requested link index.
        link_index: usize,
        /// Number of links in the model.
        link_count: usize,
    },
    /// A supplied or derived quantity is not finite.
    NonFinite {
        /// Quantity that failed admission.
        field: &'static str,
    },
    /// The plane normal has zero length.
    InvalidNormal,
}

impl fmt::Display for ArticulatedImpactError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dynamics(error) => write!(formatter, "articulated impact: {error}"),
            Self::ContactLaw(error) => write!(formatter, "normal impact: {error}"),
            Self::InvalidLinkIndex {
                link_index,
                link_count,
            } => write!(
                formatter,
                "contact link {link_index} is outside {link_count} links"
            ),
            Self::NonFinite { field } => write!(formatter, "{field} must be finite"),
            Self::InvalidNormal => formatter.write_str("contact normal must be nonzero"),
        }
    }
}

impl std::error::Error for ArticulatedImpactError {}

impl From<ArticulatedError> for ArticulatedImpactError {
    fn from(error: ArticulatedError) -> Self {
        Self::Dynamics(error)
    }
}

impl From<GaError> for ArticulatedImpactError {
    fn from(error: GaError) -> Self {
        Self::Dynamics(ArticulatedError::Geometry(error))
    }
}

impl From<NormalImpactError> for ArticulatedImpactError {
    fn from(error: NormalImpactError) -> Self {
        Self::ContactLaw(error)
    }
}

/// Resolve one frictionless impact in O(n) work, with no dense joint mass matrix.
///
/// The normal effective mass includes the six free-base velocities and every
/// articulated joint. A unit-response probe is an inverse-mass operation, not
/// a physical event: only the final impulse is checked against endpoint speed
/// limits. The existing impulse owner performs the atomic state update.
///
/// A separating point or a gap above the activation tolerance receives zero
/// impulse. Every negative gap refuses. There is no position projection,
/// friction, persistent-contact force, CCD, or multicontact impact claim.
pub fn resolve_articulated_plane_impact(
    model: &ArticulatedModel,
    base: FreeFloatingBaseState,
    joint_position: &[f64],
    joint_velocity: &[f64],
    contact: ArticulatedPointPlaneContact,
) -> Result<ArticulatedPlaneImpactResult, ArticulatedImpactError> {
    if contact.link_index >= model.link_count() {
        return Err(ArticulatedImpactError::InvalidLinkIndex {
            link_index: contact.link_index,
            link_count: model.link_count(),
        });
    }
    check_vector(contact.point_body, "point_body")?;
    check_vector(contact.plane_point_world, "plane_point_world")?;
    check_vector(contact.plane_velocity_world, "plane_velocity_world")?;
    check_vector(contact.normal_world, "normal_world")?;
    let scale = contact
        .normal_world
        .x
        .abs()
        .max(contact.normal_world.y.abs())
        .max(contact.normal_world.z.abs());
    if scale == 0.0 {
        return Err(ArticulatedImpactError::InvalidNormal);
    }
    let scaled = Vec3::new(
        contact.normal_world.x / scale,
        contact.normal_world.y / scale,
        contact.normal_world.z / scale,
    );
    let normal = scaled.scale(1.0 / scaled.norm());
    let before = forward_kinematics(
        model,
        BaseState::prescribed(base.world_from_base, base.twist_body, Twist::zero()),
        joint_position,
        joint_velocity,
    )?;
    let pose = before.world_from_link[contact.link_index];
    let contact_point_world = pose.transform_point(contact.point_body)?;
    let gap = finite(
        normal.dot(contact_point_world - contact.plane_point_world),
        "gap",
    )?;
    let normal_body = pose.rotation().inverse().rotate(normal)?;
    let plane_speed = finite(
        normal.dot(contact.plane_velocity_world),
        "plane_normal_speed",
    )?;
    let before_twist = before.body_twist[contact.link_index];
    let speed_before = finite(
        normal_body.dot(point_velocity(before_twist, contact.point_body)) - plane_speed,
        "relative_normal_speed_before",
    )?;

    let joint_impulses = vec![0.0; model.dof_count()];
    let mut body_impulses = vec![Wrench::default(); model.link_count()];
    let direction = Wrench::new(contact.point_body.cross(normal_body), normal_body);
    body_impulses[contact.link_index] = direction;
    let probe =
        free_floating_impulse_delta(model, joint_position, &joint_impulses, &body_impulses)?;
    let inverse_effective_mass = finite(
        normal_body.dot(point_velocity(
            probe.body[contact.link_index],
            contact.point_body,
        )),
        "inverse_effective_mass",
    )?;
    let normal_law = solve_normal_impact(NormalImpactInput {
        gap,
        activation_tolerance: contact.activation_tolerance,
        relative_normal_speed: speed_before,
        inverse_effective_mass,
        restitution: contact.restitution,
    })?;
    let magnitude = normal_law.normal_impulse;
    body_impulses[contact.link_index] = Wrench::new(
        direction.torque.scale(magnitude),
        direction.force.scale(magnitude),
    );
    let response = apply_free_floating_impulses(
        model,
        base,
        joint_position,
        joint_velocity,
        &joint_impulses,
        &body_impulses,
    )?;
    let after = forward_kinematics(
        model,
        BaseState::prescribed(
            response.base_after.world_from_base,
            response.base_after.twist_body,
            Twist::zero(),
        ),
        joint_position,
        &response.joint_velocity_after,
    )?;
    let speed_after = finite(
        normal_body.dot(point_velocity(
            after.body_twist[contact.link_index],
            contact.point_body,
        )) - plane_speed,
        "relative_normal_speed_after",
    )?;
    let plane_work = finite(magnitude * plane_speed, "plane_work")?;
    let energy_balance_residual = finite(
        response.momentum_after.kinetic_energy
            - response.momentum_before.kinetic_energy
            - plane_work
            + normal_law.dissipated_energy,
        "energy_balance_residual",
    )?;
    Ok(ArticulatedPlaneImpactResult {
        response,
        normal_law,
        gap,
        inverse_effective_mass,
        relative_normal_speed_before: speed_before,
        relative_normal_speed_after: speed_after,
        contact_point_world,
        plane_work,
        energy_balance_residual,
    })
}

fn point_velocity(twist: Twist, point: Vec3) -> Vec3 {
    twist.linear + twist.angular.cross(point)
}

fn finite(value: f64, field: &'static str) -> Result<f64, ArticulatedImpactError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ArticulatedImpactError::NonFinite { field })
    }
}

fn check_vector(value: Vec3, field: &'static str) -> Result<(), ArticulatedImpactError> {
    finite(value.x, field)?;
    finite(value.y, field)?;
    finite(value.z, field)?;
    Ok(())
}
