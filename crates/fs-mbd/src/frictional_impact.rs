//! Simultaneous normal and Coulomb impulses on native rigid-body momenta.
//!
//! Contact-space inertia comes from the original masses, full body inertias,
//! poses and point arms. Normal/tangent cross terms are retained. Geometry is
//! supplied by the collision owner; this frozen-configuration event does not
//! advance poses, find collisions or replace a contact time integrator.

/// Configuration-admitted force/contact/drift on the original rigid bodies.
pub mod step;

/// Automatic swept sphere geometry and event-resolved native body trajectories.
pub mod spheres;

/// Automatic sphere/plane scenes with sustained native contact and path admission.
pub mod sphere_scene;

/// Directional friction ellipses on the same simultaneous native-body owner.
pub mod elliptic;

use crate::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactError, FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;

/// Original COM-relative material point on a dynamic body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidImpactPoint {
    pub body: usize,
    pub arm_body: Vec3,
}

/// Opposite material point or prescribed translating support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RigidImpactCounterpart {
    Body(RigidImpactPoint),
    PrescribedPlane {
        point_world: Vec3,
        velocity_world: Vec3,
    },
}

/// One oriented geometric witness, with source-resolved contact coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidFrictionalContact<F = f64> {
    pub id: u64,
    pub point: RigidImpactPoint,
    pub other: RigidImpactCounterpart,
    /// Points from the counterpart toward the positive-impulse body.
    pub normal_world: Vec3,
    /// Tangent hint; its normal component is removed before normalization.
    pub tangent_world: Vec3,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: F,
    pub normal_law: NormalVelocityLaw,
}

/// Bounded assembly and independent native-physics acceptance tolerances.
#[derive(Clone, Copy, Debug)]
pub struct RigidImpactControls {
    pub maximum_bodies: usize,
    pub maximum_contacts: usize,
    /// Dynamic witnesses must lie on the same normal line within this bound.
    /// Their actual separation remains in the angular-momentum check.
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_mps: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
    pub contact: FrictionalContactConfig,
}

/// Actual body observables and the simultaneous constitutive solve receipt.
#[derive(Clone, Debug)]
pub struct RigidFrictionalImpactReport<R = FrictionalContactResult> {
    /// All contact-space vectors below follow this canonical order.
    pub contact_ids: Vec<u64>,
    pub contact: R,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Signed work supplied by explicitly moving prescribed supports.
    pub support_work_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub support_impulse_ns: Vec3,
    pub support_angular_impulse_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    pub endpoint_velocity_residual_mps: f64,
}

/// A failure publishes no body update, including failure after impulse staging.
#[derive(Debug)]
pub enum RigidFrictionalImpactError<R = FrictionalContactResult, E = FrictionalContactError> {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Contact(E),
    /// A unique body velocity does not determine the work exchanged between
    /// independently moving supports when contact impulses remain indeterminate.
    /// Velocity-only admission permits one common support translation, whose
    /// work is that velocity dotted with the uniquely determined net impulse.
    UnresolvedSupportWork,
    Cancelled,
    Audit(Box<RigidFrictionalImpactReport<R>>),
}
impl<R, E> From<DynamicsError> for RigidFrictionalImpactError<R, E> {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<FrictionalContactError> for RigidFrictionalImpactError {
    fn from(error: FrictionalContactError) -> Self {
        Self::Contact(error)
    }
}
impl<R: core::fmt::Debug, E: core::fmt::Debug> core::fmt::Display
    for RigidFrictionalImpactError<R, E>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid frictional impact: {self:?}")
    }
}
impl<R: core::fmt::Debug, E: core::fmt::Debug> std::error::Error
    for RigidFrictionalImpactError<R, E>
{
}

/// Original preimpact point kinematics, reconstructed once from native bodies.
/// Query order follows contact IDs; `original_index` addresses the caller's
/// geometry, material binding and history slices without an ID search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidImpactKinematics {
    pub original_index: usize,
    pub contact_id: u64,
    /// Exact normalized world axes used in the shared mobility assembly:
    /// [normal, tangent 1, tangent 2], with tangent 2 = normal cross tangent 1.
    pub axes_world: [Vec3; 3],
    pub point_world: Vec3,
    /// Original dynamic witness, or the normal projection of `point_world`
    /// onto a prescribed plane. A plane's reference anchor is not its local
    /// contact point. This query metadata does not relocate either body arm.
    pub counterpart_point_world: Vec3,
    pub gap_m: f64,
    /// [normal, tangent 1, tangent 2], with positive normal speed separating.
    pub relative_velocity_m_s: [f64; 3],
}

/// Only physical observables enter the common body audit. Any numerical
/// witness remains in its original constitutive report and coordinate system.
struct PhysicalImpulseView<'a> {
    impulses: &'a [[f64; 3]],
    relative_velocities_after: &'a [[f64; 3]],
    relative_work: f64,
    dissipated_energy: f64,
    velocity_only: bool,
}
trait ImpactResponse {
    fn physical(&self) -> PhysicalImpulseView<'_>;
}
impl ImpactResponse for FrictionalContactResult {
    fn physical(&self) -> PhysicalImpulseView<'_> {
        PhysicalImpulseView {
            impulses: &self.impulses,
            relative_velocities_after: &self.relative_velocities_after,
            relative_work: self.relative_work,
            dissipated_energy: self.dissipated_energy,
            velocity_only: matches!(
                self.determinacy,
                FrictionalContactDeterminacy::UniqueVelocity(_)
            ),
        }
    }
}

struct Influence {
    body: usize,
    factor: [[f64; 6]; 3],
}
struct PreparedContact {
    axes: [Vec3; 3],
    positive_point_world: Vec3,
    influences: Vec<Influence>,
}

/// Apply one globally admitted normal/friction event to original body states.
///
/// The full contact matrix is assembled from shared native body influences.
/// Bodies absent from a contact do not enter its matrix products. Normal and
/// tangential impulses are solved simultaneously, never by subsequent clipping.
/// Actual post-impulse point velocities, energy and momentum are reconstructed
/// independently before publication. Dynamic contact witnesses are never moved
/// to an invented common point: a finite gap can therefore refuse the angular
/// balance even when a geometric activation tolerance admitted it.
///
/// Redundant contact rows may be admitted when the original contact operator
/// proves a unique body velocity. Their individual impulses then remain one
/// admissible representative, not uniquely determined local forces or heat.
/// This velocity-only route admits dynamic bodies and supports with one common
/// translation. Their total work follows the uniquely determined net impulse.
/// Independently moving supports require unique impulses to determine work.
/// Other unresolved manifolds retain a typed determinacy refusal. Impact loss
/// remains unallocated mechanical loss; a physical thermal/plastic/damage law
/// must determine its disposition. No finite-patch preload, spin friction or
/// pose integration is silently substituted for this point-impact law.
#[allow(clippy::too_many_lines)]
pub fn apply_rigid_frictional_impacts(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact],
    controls: RigidImpactControls,
    gate: &CancelGate,
) -> Result<RigidFrictionalImpactReport, RigidFrictionalImpactError> {
    apply_impacts_with(
        bodies,
        properties,
        contacts,
        controls,
        gate,
        |kinematics, delassus| {
            let laws: Vec<_> = kinematics
                .iter()
                .map(|query| {
                    let c = contacts[query.original_index];
                    FrictionalContact {
                        id: c.id,
                        gap: query.gap_m,
                        activation_tolerance: c.activation_tolerance_m,
                        relative_velocity: query.relative_velocity_m_s,
                        friction_coefficient: c.friction_coefficient,
                        normal_law: c.normal_law,
                    }
                })
                .collect();
            solve_frictional_contacts(&laws, delassus, controls.contact, || gate.is_requested())
        },
    )
}

#[allow(clippy::too_many_lines)]
fn apply_impacts_with<F: Copy, R: ImpactResponse, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact<F>],
    controls: RigidImpactControls,
    gate: &CancelGate,
    solve: impl FnOnce(&[RigidImpactKinematics], &[f64]) -> Result<R, E>,
) -> Result<RigidFrictionalImpactReport<R>, RigidFrictionalImpactError<R, E>> {
    let invalid = RigidFrictionalImpactError::Invalid;
    poll_with(gate)?;
    if properties.len() != bodies.len()
        || controls.maximum_bodies == 0
        || controls.maximum_contacts == 0
        || bodies.len() > controls.maximum_bodies
        || contacts.len() > controls.maximum_contacts
        || contacts.len() > controls.contact.maximum_contacts
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_mps,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid(
            "matching bounded bodies and positive finite controls required",
        ));
    }
    let dimension = contacts
        .len()
        .checked_mul(3)
        .ok_or(invalid("contact dimension overflow"))?;
    let entries = dimension
        .checked_mul(dimension)
        .ok_or(invalid("contact matrix overflow"))?;
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_by_key(|&i| contacts[i].id);
    let mut prepared = Vec::with_capacity(contacts.len());
    let mut kinematics = Vec::with_capacity(contacts.len());
    for (slot, &index) in order.iter().enumerate() {
        poll_with(gate)?;
        let c = contacts[index];
        if c.point.body >= bodies.len() || (slot > 0 && contacts[order[slot - 1]].id == c.id) {
            return Err(invalid(
                "valid body indices and distinct contact identities required",
            ));
        }
        let n = c.normal_world.normalized("contact.normal")?;
        let t = c
            .tangent_world
            .sub(n.scale(n.dot(c.tangent_world)))
            .normalized("contact.tangent")?;
        let axes = [n, t, n.cross(t)];
        let a =
            bodies[c.point.body].point_kinematics(properties[c.point.body], c.point.arm_body)?;
        let mut influences = vec![influence(bodies, properties, c.point, axes)];
        let (other_point, other_velocity) = match c.other {
            RigidImpactCounterpart::Body(point) => {
                if point.body >= bodies.len() || point.body == c.point.body {
                    return Err(invalid("counterpart must name another valid dynamic body"));
                }
                let b =
                    bodies[point.body].point_kinematics(properties[point.body], point.arm_body)?;
                let separation = a.point_world.sub(b.point_world);
                if magnitude(separation.sub(n.scale(n.dot(separation))))
                    > controls.geometry_tolerance_m
                {
                    return Err(invalid(
                        "dynamic contact witnesses do not share a normal line",
                    ));
                }
                let mut opposite = influence(bodies, properties, point, axes);
                for row in &mut opposite.factor {
                    for value in row {
                        *value = -*value;
                    }
                }
                influences.push(opposite);
                (b.point_world, b.point_velocity_world)
            }
            RigidImpactCounterpart::PrescribedPlane {
                point_world,
                velocity_world,
            } => {
                if !point_world.is_finite() || !velocity_world.is_finite() {
                    return Err(invalid(
                        "prescribed support geometry and velocity must be finite",
                    ));
                }
                (point_world, velocity_world)
            }
        };
        let relative = a.point_velocity_world.sub(other_velocity);
        let gap = n.dot(a.point_world.sub(other_point));
        let counterpart_point_world = match c.other {
            RigidImpactCounterpart::Body(_) => other_point,
            RigidImpactCounterpart::PrescribedPlane { .. } => a.point_world.sub(n.scale(gap)),
        };
        kinematics.push(RigidImpactKinematics {
            original_index: index,
            contact_id: c.id,
            axes_world: axes,
            point_world: a.point_world,
            counterpart_point_world,
            gap_m: gap,
            relative_velocity_m_s: axes.map(|axis| relative.dot(axis)),
        });
        prepared.push(PreparedContact {
            axes,
            positive_point_world: a.point_world,
            influences,
        });
    }
    let mut delassus = vec![0.0; entries];
    for i in 0..dimension {
        poll_with(gate)?;
        for j in i..dimension {
            let mut value = 0.0;
            for a in &prepared[i / 3].influences {
                for b in &prepared[j / 3].influences {
                    if a.body == b.body {
                        for k in 0..6 {
                            value += a.factor[i % 3][k] * b.factor[j % 3][k];
                        }
                    }
                }
            }
            delassus[i * dimension + j] = value;
            delassus[j * dimension + i] = value;
        }
    }
    let contact = solve(&kinematics, &delassus).map_err(RigidFrictionalImpactError::Contact)?;
    let physical = contact.physical();
    if physical.impulses.len() != contacts.len()
        || physical.relative_velocities_after.len() != contacts.len()
    {
        return Err(invalid(
            "constitutive response does not match contact dimensions",
        ));
    }
    if physical.velocity_only {
        // A null impulse changes no body's momentum. Consequently it changes
        // neither total support impulse nor U dot total_support_impulse for a
        // common prescribed translation U. Distinct support velocities need
        // a separate work-determinacy witness, which this owner does not infer.
        let mut common_support_velocity = None;
        for c in contacts {
            if let RigidImpactCounterpart::PrescribedPlane { velocity_world, .. } = c.other {
                if common_support_velocity.is_some_and(|u| u != velocity_world) {
                    return Err(RigidFrictionalImpactError::UnresolvedSupportWork);
                }
                common_support_velocity = Some(velocity_world);
            }
        }
    }
    let mut candidate = bodies.to_vec();
    let mut support_work = 0.0;
    let mut support_impulse = Vec3::ZERO;
    let mut support_angular = Vec3::ZERO;
    for (slot, &index) in order.iter().enumerate() {
        poll_with(gate)?;
        let c = contacts[index];
        let impulse = world_impulse(prepared[slot].axes, physical.impulses[slot]);
        candidate[c.point.body] = candidate[c.point.body]
            .apply_impulse_at_body_point(properties[c.point.body], c.point.arm_body, impulse)?
            .state_after;
        match c.other {
            RigidImpactCounterpart::Body(point) => {
                candidate[point.body] = candidate[point.body]
                    .apply_impulse_at_body_point(
                        properties[point.body],
                        point.arm_body,
                        impulse.scale(-1.0),
                    )?
                    .state_after;
            }
            RigidImpactCounterpart::PrescribedPlane { velocity_world, .. } => {
                support_work += impulse.dot(velocity_world);
                support_impulse = support_impulse.add(impulse);
                support_angular =
                    support_angular.add(prepared[slot].positive_point_world.cross(impulse));
            }
        }
    }
    let (before, linear_before, angular_before) = totals_with(bodies, properties, gate)?;
    let (after, linear_after, angular_after) = totals_with(&candidate, properties, gate)?;
    let mut relative_work = 0.0;
    let mut endpoint_error = 0.0_f64;
    for (slot, &index) in order.iter().enumerate() {
        poll_with(gate)?;
        let c = contacts[index];
        let a =
            candidate[c.point.body].point_kinematics(properties[c.point.body], c.point.arm_body)?;
        let other = match c.other {
            RigidImpactCounterpart::Body(point) => {
                candidate[point.body]
                    .point_kinematics(properties[point.body], point.arm_body)?
                    .point_velocity_world
            }
            RigidImpactCounterpart::PrescribedPlane { velocity_world, .. } => velocity_world,
        };
        let relative = a.point_velocity_world.sub(other);
        for axis in 0..3 {
            let velocity = relative.dot(prepared[slot].axes[axis]);
            endpoint_error = endpoint_error
                .max((velocity - physical.relative_velocities_after[slot][axis]).abs());
            relative_work += 0.5
                * (kinematics[slot].relative_velocity_m_s[axis] + velocity)
                * physical.impulses[slot][axis];
        }
    }
    let report = RigidFrictionalImpactReport {
        contact_ids: order.iter().map(|&i| contacts[i].id).collect(),
        kinetic_before_j: before,
        kinetic_after_j: after,
        support_work_j: support_work,
        relative_impulse_work_j: relative_work,
        work_energy_residual_j: after - before - support_work - relative_work,
        total_energy_residual_j: after - before - support_work + physical.dissipated_energy,
        support_impulse_ns: support_impulse,
        support_angular_impulse_nms: support_angular,
        linear_momentum_residual_ns: linear_after.sub(linear_before).sub(support_impulse),
        angular_momentum_residual_nms: angular_after.sub(angular_before).sub(support_angular),
        endpoint_velocity_residual_mps: endpoint_error,
        contact,
    };
    if [
        before,
        after,
        support_work,
        relative_work,
        report.work_energy_residual_j,
        report.total_energy_residual_j,
        endpoint_error,
        magnitude(report.linear_momentum_residual_ns),
        magnitude(report.angular_momentum_residual_nms),
    ]
    .iter()
    .any(|x| !x.is_finite())
        || endpoint_error > controls.velocity_tolerance_mps
        || magnitude(report.linear_momentum_residual_ns) > controls.linear_momentum_tolerance_ns
        || magnitude(report.angular_momentum_residual_nms) > controls.angular_momentum_tolerance_nms
        || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
        || report.total_energy_residual_j.abs() > controls.energy_tolerance_j
        || (relative_work - report.contact.physical().relative_work).abs()
            > controls.energy_tolerance_j
    {
        return Err(RigidFrictionalImpactError::Audit(Box::new(report)));
    }
    poll_with(gate)?;
    bodies.copy_from_slice(&candidate);
    Ok(report)
}

fn influence(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    point: RigidImpactPoint,
    axes: [Vec3; 3],
) -> Influence {
    let body = bodies[point.body];
    let property = properties[point.body];
    let inverse_mass_root = property.mass().sqrt().recip();
    let mut factor = [[0.0; 6]; 3];
    for row in 0..3 {
        let translation = components(axes[row].scale(inverse_mass_root));
        let direction = body.pose().orientation().rotate_world_to_body(axes[row]);
        let rotation =
            components(property.angular_mobility_factor(point.arm_body.cross(direction)));
        for axis in 0..3 {
            factor[row][axis] = translation[axis];
            factor[row][axis + 3] = rotation[axis];
        }
    }
    Influence {
        body: point.body,
        factor,
    }
}
fn world_impulse(axes: [Vec3; 3], local: [f64; 3]) -> Vec3 {
    axes[0]
        .scale(local[0])
        .add(axes[1].scale(local[1]))
        .add(axes[2].scale(local[2]))
}
fn totals(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    gate: &CancelGate,
) -> Result<(f64, Vec3, Vec3), RigidFrictionalImpactError> {
    totals_with(bodies, properties, gate)
}
fn totals_with<R, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    gate: &CancelGate,
) -> Result<(f64, Vec3, Vec3), RigidFrictionalImpactError<R, E>> {
    let mut energy = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    for (&body, &property) in bodies.iter().zip(properties) {
        poll_with(gate)?;
        energy += body.kinetic_energy(property)?;
        linear = linear.add(body.linear_momentum_world());
        angular = angular
            .add(
                body.pose()
                    .position_world()
                    .cross(body.linear_momentum_world()),
            )
            .add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(body.angular_momentum_body()),
            );
    }
    Ok((energy, linear, angular))
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn poll(gate: &CancelGate) -> Result<(), RigidFrictionalImpactError> {
    poll_with(gate)
}
fn poll_with<R, E>(gate: &CancelGate) -> Result<(), RigidFrictionalImpactError<R, E>> {
    if gate.is_requested() {
        Err(RigidFrictionalImpactError::Cancelled)
    } else {
        Ok(())
    }
}
