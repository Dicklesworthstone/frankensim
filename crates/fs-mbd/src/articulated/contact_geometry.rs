//! Shared material-point geometry for dry and coupled-inertia contact owners.
//!
//! This extracts the original articulated point/normal preparation. It does
//! not select contacts, repair penetration, supply inertia or solve an impact.

use super::Kinematics;
use fs_ga::{GaError, Vec3, Wrench};
use std::borrow::Borrow;

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

/// One material point's body-coordinate wrench for a unit normal impulse.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointDirection {
    /// Original material point and system/link indices.
    pub point: LinkContactPoint,
    /// Force direction and its moment about the link origin.
    pub direction: Wrench,
}

/// Admitted geometry and measured speed of one caller-supplied normal contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalContactGeometry {
    /// Positive side of the normal row.
    pub first: PointDirection,
    /// Equal/opposite dynamic side, absent for a prescribed plane.
    pub second: Option<PointDirection>,
    /// Unit separating axis in world coordinates.
    pub axis_world: Vec3,
    /// World point receiving the positive impulse.
    pub world_point: Vec3,
    /// Prescribed support speed along the normal, zero for a dynamic pair.
    pub plane_speed: f64,
    /// Signed normal separation; the contact law admits nonnegative gaps only.
    pub gap: f64,
    /// Actual relative speed of the supplied points/support.
    pub relative_normal_speed: f64,
}

/// Original material-point geometry in an orthonormal contact frame.
/// The first row is normal; the others are tangent and normal cross tangent.
/// Every row retains the same signed NORMAL gap, while its speed and wrench
/// are along that row's axis. No constitutive coefficient or impact law enters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactFrameGeometry {
    pub rows: [NormalContactGeometry; 3],
    /// Other material point, or the first point's normal projection onto a plane.
    pub counterpart_point_world: Vec3,
}

/// Geometry-only refusal, independent of a contact law or mass operator.
#[derive(Clone, Debug, PartialEq)]
pub enum ContactGeometryError {
    /// Invalid point index or nonfinite input/derived arithmetic.
    Invalid(&'static str),
    /// Dynamic witnesses would create an artificial couple.
    NonCollinearContact {
        /// Contact index supplied to the preparation call.
        index: usize,
    },
    /// Original Lie-group geometry refusal.
    Geometry(GaError),
}
impl core::fmt::Display for ContactGeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => f.write_str(reason),
            Self::NonCollinearContact { index } => {
                write!(f, "contact {index} has noncollinear dynamic witnesses")
            }
            Self::Geometry(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for ContactGeometryError {}
impl From<GaError> for ContactGeometryError {
    fn from(error: GaError) -> Self {
        Self::Geometry(error)
    }
}

/// Transform one original material-point normal row without prescribing inertia.
pub fn prepare_point(
    point: LinkContactPoint,
    normal_world: Vec3,
    kinematics: &[impl Borrow<Kinematics>],
) -> Result<(PointDirection, Vec3), ContactGeometryError> {
    let system = kinematics
        .get(point.system_index)
        .ok_or(ContactGeometryError::Invalid(
            "contact system index is out of range",
        ))?
        .borrow();
    let pose =
        system
            .world_from_link
            .get(point.link_index)
            .ok_or(ContactGeometryError::Invalid(
                "contact link index is out of range",
            ))?;
    if system.body_twist.len() != system.world_from_link.len() {
        return Err(ContactGeometryError::Invalid(
            "contact pose/twist counts differ",
        ));
    }
    vector(point.point_body)?;
    let direction = pose.rotation().inverse().rotate(normal_world)?;
    Ok((
        PointDirection {
            point,
            direction: Wrench::new(point.point_body.cross(direction), direction),
        },
        pose.transform_point(point.point_body)?,
    ))
}

/// Normalize a nonzero finite contact axis using the original scaled norm.
pub fn unit_normal(normal: Vec3) -> Result<Vec3, ContactGeometryError> {
    vector(normal)?;
    let scale = magnitude(normal);
    if scale == 0.0 {
        return Err(ContactGeometryError::Invalid(
            "contact normal must be nonzero",
        ));
    }
    let scaled = Vec3::new(normal.x / scale, normal.y / scale, normal.z / scale);
    Ok(scaled.scale(1.0 / scaled.norm()))
}

/// Prepare the original point/normal/gap and shared line-of-action checks.
/// Negative gaps are retained for the existing contact law to refuse; they
/// are never clamped to an activation tolerance or repaired here.
pub fn prepare_normal_contact(
    contact: &ArticulatedNormalContact,
    index: usize,
    kinematics: &[impl Borrow<Kinematics>],
) -> Result<NormalContactGeometry, ContactGeometryError> {
    prepare_normal_geometry(
        contact.point,
        contact.other,
        contact.normal_world,
        index,
        kinematics,
    )
}

fn prepare_normal_geometry(
    point: LinkContactPoint,
    other: ContactCounterpart,
    normal_world: Vec3,
    index: usize,
    kinematics: &[impl Borrow<Kinematics>],
) -> Result<NormalContactGeometry, ContactGeometryError> {
    let normal = unit_normal(normal_world)?;
    let (first, world_point) = prepare_point(point, normal, kinematics)?;
    let (second, other_point, plane_speed) = match other {
        ContactCounterpart::Body(point) => {
            let (term, position) = prepare_point(point, normal.scale(-1.0), kinematics)?;
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
            return Err(ContactGeometryError::NonCollinearContact { index });
        }
    }
    let speed = core::iter::once(first)
        .chain(second)
        .map(|term| {
            term.direction.pairing(
                kinematics[term.point.system_index].borrow().body_twist[term.point.link_index],
            )
        })
        .sum::<f64>()
        - plane_speed;
    Ok(NormalContactGeometry {
        first,
        second,
        axis_world: normal,
        world_point,
        plane_speed,
        gap,
        relative_normal_speed: finite(speed)?,
    })
}

/// Prepare normal and two tangent rows using the native point and frame rules.
/// Dynamic witnesses must be on distinct rigid links and share a normal line.
/// A finite-gap tangential impulse can still create an artificial couple; the
/// mechanical consumer must audit its actual spatial angular momentum.
pub fn prepare_contact_frame(
    point: LinkContactPoint,
    other: ContactCounterpart,
    normal_world: Vec3,
    tangent_world: Vec3,
    index: usize,
    kinematics: &[impl Borrow<Kinematics>],
) -> Result<ContactFrameGeometry, ContactGeometryError> {
    if matches!(other, ContactCounterpart::Body(p)
        if p.system_index == point.system_index && p.link_index == point.link_index)
    {
        return Err(ContactGeometryError::Invalid(
            "contact requires distinct rigid links",
        ));
    }
    let normal = prepare_normal_geometry(point, other, normal_world, index, kinematics)?;
    // Scaling precedes crossing; the two crosses discard the hint's normal
    // component without a cancellation-prone nearly parallel subtraction.
    let hint = unit_normal(tangent_world)?;
    let second_tangent = unit_normal(normal.axis_world.cross(hint))?;
    let first_tangent = unit_normal(second_tangent.cross(normal.axis_world))?;
    let mut rows = [normal; 3];
    for (i, axis) in [first_tangent, second_tangent].into_iter().enumerate() {
        let (first, world_point) = prepare_point(point, axis, kinematics)?;
        let (second, plane_speed) = match other {
            ContactCounterpart::Body(point) => {
                let (term, _) = prepare_point(point, axis.scale(-1.0), kinematics)?;
                (Some(term), 0.0)
            }
            ContactCounterpart::PrescribedPlane { velocity_world, .. } => {
                (None, finite(axis.dot(velocity_world))?)
            }
        };
        let speed = first.direction.pairing(
            kinematics[first.point.system_index].borrow().body_twist[first.point.link_index],
        ) + second.map_or(0.0, |term| {
            term.direction.pairing(
                kinematics[term.point.system_index].borrow().body_twist[term.point.link_index],
            )
        }) - plane_speed;
        rows[i + 1] = NormalContactGeometry {
            first,
            second,
            axis_world: axis,
            world_point,
            plane_speed,
            gap: normal.gap,
            relative_normal_speed: finite(speed)?,
        };
    }
    let counterpart_point_world = match other {
        ContactCounterpart::Body(point) => kinematics[point.system_index].borrow().world_from_link
            [point.link_index]
            .transform_point(point.point_body)?,
        ContactCounterpart::PrescribedPlane { .. } => {
            normal.world_point - normal.axis_world.scale(normal.gap)
        }
    };
    vector(counterpart_point_world)?;
    Ok(ContactFrameGeometry {
        rows,
        counterpart_point_world,
    })
}

fn magnitude(value: Vec3) -> f64 {
    value.x.abs().max(value.y.abs()).max(value.z.abs())
}
fn vector(value: Vec3) -> Result<(), ContactGeometryError> {
    finite(value.x)?;
    finite(value.y)?;
    finite(value.z)?;
    Ok(())
}
fn finite(value: f64) -> Result<f64, ContactGeometryError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ContactGeometryError::Invalid("nonfinite impact arithmetic"))
    }
}
