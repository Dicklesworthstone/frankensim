//! Simultaneous frictional impacts between one deforming solid and rigid bodies.
//!
//! The complete contact mobility is B_s M_s,ff^-1 B_s^T plus the original
//! rigid-body mass/full-inertia response. Cross-point blocks are retained: two
//! rigid bodies striking one solid are coupled even when they do not touch each
//! other. The existing fs-contact owner determines all normal and tangential
//! impulses together. Its determinacy and passivity refusals remain operative.
//!
//! Actual native velocities, spin/orbital momentum, stationary-support reaction
//! and impulse work are checked before returning either candidate. Geometry is
//! supplied to the core event; neither witness is relocated to hide a couple
//! from a finite gap. The `surface` owner detects current sphere/solid witnesses;
//! `step` composes bounded, path-admitted transient trajectories through the
//! original native integrators. No penetration repair, sustained-contact solve
//! or impact-loss heat partition is inferred. Material-point weights and rigid
//! arms remain the original ones.

/// Directional interface sources through the same coupled mechanical event.
pub mod elliptic;

/// Current finite-solid surface witnesses for native rigid contact.
pub mod surface;

/// Configuration-admitted coupled force/contact/trajectory candidates.
pub mod step;

use super::impulse::{FiniteSolidImpulse, FiniteSolidImpulseControls, FiniteSolidPoint};
use super::{Cx, FiniteSolidBody, FiniteSolidError};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResult,
    NormalVelocityLaw, solve_frictional_contacts,
};
use fs_mbd::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_solver::SolveReport;

/// A rigid witness and a P1 material witness on the same contact normal line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidRigidContact<F = f64> {
    pub id: u64,
    pub rigid_body: usize,
    /// Centre-of-mass-relative point in this rigid body's material frame [m].
    pub arm_body_m: Vec3,
    pub solid_point: FiniteSolidPoint,
    /// Points from the solid toward the rigid body; positive velocity separates.
    pub normal_world: Vec3,
    /// Its normal component is removed before forming the contact frame.
    pub tangent_world: Vec3,
    pub activation_tolerance_m: f64,
    /// Explicit coefficient resolved by the caller's interface law/source.
    pub friction_coefficient: F,
    pub normal_law: NormalVelocityLaw,
}

impl<F> FiniteSolidRigidContact<F> {
    /// Retain the original material witnesses while binding another friction law.
    #[must_use]
    pub fn with_friction<T>(self, friction_coefficient: T) -> FiniteSolidRigidContact<T> {
        FiniteSolidRigidContact {
            id: self.id,
            rigid_body: self.rigid_body,
            arm_body_m: self.arm_body_m,
            solid_point: self.solid_point,
            normal_world: self.normal_world,
            tangent_world: self.tangent_world,
            activation_tolerance_m: self.activation_tolerance_m,
            friction_coefficient,
            normal_law: self.normal_law,
        }
    }
}

/// Original native point velocities and frame supplied to material queries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidRigidContactKinematics {
    /// Original caller order, before canonical contact-ID sorting.
    pub original_index: usize,
    pub contact_id: u64,
    /// Normal, tangent, bitangent; the exact orthonormal frame used by the solve.
    pub axes_world: [Vec3; 3],
    pub rigid_point_world: Vec3,
    pub solid_point_world: Vec3,
    pub gap_m: f64,
    /// Rigid minus solid point velocity, projected onto `axes_world` [m/s].
    pub relative_velocity_m_s: [f64; 3],
}

/// Bounded original mass/contact solves and independent physical tolerances.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidRigidContactControls {
    pub maximum_rigid_bodies: usize,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
    pub solid_impulse: FiniteSolidImpulseControls,
    pub contact: FrictionalContactConfig,
}
impl Default for FiniteSolidRigidContactControls {
    fn default() -> Self {
        Self {
            maximum_rigid_bodies: 64,
            geometry_tolerance_m: 1e-10,
            velocity_tolerance_m_s: 1e-8,
            linear_tolerance_kg_m_s: 1e-9,
            angular_tolerance_kg_m2_s: 1e-9,
            energy_tolerance_j: 1e-8,
            solid_impulse: FiniteSolidImpulseControls::default(),
            contact: FrictionalContactConfig::default(),
        }
    }
}

/// All contact arrays follow stable contact-ID order, independent of input order.
#[derive(Clone, Debug)]
pub struct FiniteSolidRigidContactReport<R = FrictionalContactResult> {
    pub contact_ids: Vec<u64>,
    /// Full rigid-plus-solid mobility, row-major 3N by 3N [1/kg].
    pub contact_mobility_per_kg: Vec<f64>,
    pub solid_mass_solves: Vec<SolveReport>,
    pub contact: R,
    pub solid_impulse: FiniteSolidImpulse,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Evaluated from actual momentum increments, not subtraction of energies.
    pub kinetic_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub dissipation_residual_j: f64,
    /// Change minus the external reaction at the solid's stationary supports.
    pub linear_momentum_residual_kg_m_s: Vec3,
    pub angular_momentum_residual_kg_m2_s: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
}

#[derive(Debug)]
pub enum FiniteSolidRigidContactError<R = FrictionalContactResult, E = FrictionalContactError> {
    Invalid(&'static str),
    Solid(FiniteSolidError),
    Rigid(DynamicsError),
    Contact(E),
    Cancelled,
    /// Neither the solid nor any rigid state is published by a failed audit.
    Audit(Box<FiniteSolidRigidContactReport<R>>),
}
impl<R, E> From<FiniteSolidError> for FiniteSolidRigidContactError<R, E> {
    fn from(error: FiniteSolidError) -> Self {
        Self::Solid(error)
    }
}
impl<R, E> From<DynamicsError> for FiniteSolidRigidContactError<R, E> {
    fn from(error: DynamicsError) -> Self {
        Self::Rigid(error)
    }
}
impl From<FrictionalContactError> for FiniteSolidRigidContactError {
    fn from(error: FrictionalContactError) -> Self {
        Self::Contact(error)
    }
}
impl<R: core::fmt::Debug, E: core::fmt::Debug> core::fmt::Display
    for FiniteSolidRigidContactError<R, E>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite solid/rigid contact: {self:?}")
    }
}
impl<R: core::fmt::Debug, E: core::fmt::Debug> std::error::Error
    for FiniteSolidRigidContactError<R, E>
{
}

/// Physical coordinates only: transformed solver witnesses never enter bodies.
struct PhysicalImpulseView<'a> {
    impulses: &'a [[f64; 3]],
    relative_velocities_after: &'a [[f64; 3]],
    relative_work: f64,
    dissipated_energy: f64,
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
        }
    }
}

struct Prepared {
    original: usize,
    axes: [Vec3; 3],
    torque_body: [Vec3; 3],
    angular_response_body: [Vec3; 3],
}

impl<'a> FiniteSolidBody<'a> {
    /// Resolve one coupled impact with Coulomb friction and a declared normal law.
    ///
    /// Both candidates are returned together only after native endpoint and
    /// physical balance checks. Fixed solid supports supply external reaction
    /// impulses and zero work. A UniqueVelocity result retains the original
    /// solver's warning that individual redundant impulses are representatives.
    /// Sources, positions, time and accepted-step history are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn rigid_contact_candidate(
        &self,
        rigid_bodies: &[RigidBodyState],
        rigid_properties: &[MassProperties],
        contacts: &[FiniteSolidRigidContact],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, Vec<RigidBodyState>, FiniteSolidRigidContactReport),
        FiniteSolidRigidContactError,
    > {
        self.rigid_contact_candidate_with(
            rigid_bodies,
            rigid_properties,
            contacts,
            origin_m,
            controls,
            cx,
            |kinematics, matrix| {
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
                solve_frictional_contacts(&laws, matrix, controls.contact, || {
                    cx.checkpoint().is_err()
                })
            },
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn rigid_contact_candidate_with<F: Copy, R: ImpactResponse, CE>(
        &self,
        rigid_bodies: &[RigidBodyState],
        rigid_properties: &[MassProperties],
        contacts: &[FiniteSolidRigidContact<F>],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        cx: &Cx<'_>,
        solve: impl FnOnce(&[FiniteSolidRigidContactKinematics], &[f64]) -> Result<R, CE>,
    ) -> Result<
        (Self, Vec<RigidBodyState>, FiniteSolidRigidContactReport<R>),
        FiniteSolidRigidContactError<R, CE>,
    > {
        use FiniteSolidRigidContactError as E;
        poll::<R, CE>(cx)?;
        let origin = vector(origin_m);
        if !origin.is_finite()
            || rigid_bodies.len() != rigid_properties.len()
            || controls.maximum_rigid_bodies == 0
            || rigid_bodies.len() > controls.maximum_rigid_bodies
            || controls.contact.maximum_contacts == 0
            || controls.contact.maximum_contacts > 21
            || contacts.len() > controls.contact.maximum_contacts
            || [
                controls.geometry_tolerance_m,
                controls.velocity_tolerance_m_s,
                controls.linear_tolerance_kg_m_s,
                controls.angular_tolerance_kg_m2_s,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(E::Invalid(
                "bounded matching bodies and positive finite controls required",
            ));
        }
        let dimension = 3 * contacts.len();
        let mut order: Vec<_> = (0..contacts.len()).collect();
        order.sort_by_key(|&i| contacts[i].id);
        let mut prepared = Vec::with_capacity(contacts.len());
        let mut queries = Vec::with_capacity(contacts.len());
        let mut rows = Vec::with_capacity(dimension);
        for (slot, &index) in order.iter().enumerate() {
            poll::<R, CE>(cx)?;
            let c = contacts[index];
            if c.rigid_body >= rigid_bodies.len()
                || (slot > 0 && contacts[order[slot - 1]].id == c.id)
            {
                return Err(E::Invalid(
                    "valid rigid body and distinct contact IDs required",
                ));
            }
            let axes = contact_frame(c.normal_world, c.tangent_world).map_err(E::Invalid)?;
            let n = axes[0];
            let rigid = rigid_bodies[c.rigid_body];
            let property = rigid_properties[c.rigid_body];
            let kinematics = rigid.point_kinematics(property, c.arm_body_m)?;
            let point_s = vector(self.point_position_m(c.solid_point, cx)?);
            let velocity_s = vector(self.point_velocity_m_s(c.solid_point, cx)?);
            let separation = kinematics.point_world.sub(point_s);
            let gap = n.dot(separation);
            if magnitude(separation.sub(n.scale(gap))) > controls.geometry_tolerance_m {
                return Err(E::Invalid("contact witnesses do not share a normal line"));
            }
            let relative = kinematics.point_velocity_world.sub(velocity_s);
            let torque_body = axes.map(|axis| {
                c.arm_body_m
                    .cross(rigid.pose().orientation().rotate_world_to_body(axis))
            });
            let mut angular_response_body = [Vec3::ZERO; 3];
            for axis in 0..3 {
                angular_response_body[axis] =
                    property.angular_velocity_body_checked(torque_body[axis])?;
                rows.push((c.solid_point, array(axes[axis])));
            }
            queries.push(FiniteSolidRigidContactKinematics {
                original_index: index,
                contact_id: c.id,
                axes_world: axes,
                rigid_point_world: kinematics.point_world,
                solid_point_world: point_s,
                gap_m: gap,
                relative_velocity_m_s: axes.map(|axis| relative.dot(axis)),
            });
            prepared.push(Prepared {
                original: index,
                axes,
                torque_body,
                angular_response_body,
            });
        }
        let mobility = self.point_mobility_matrix_per_kg(
            &rows,
            3 * controls.contact.maximum_contacts,
            controls.solid_impulse.mass_solve,
            cx,
        )?;
        let mut matrix = mobility.matrix_per_kg;
        for i in 0..dimension {
            poll::<R, CE>(cx)?;
            let a = &prepared[i / 3];
            let body_a = contacts[a.original].rigid_body;
            for j in i..dimension {
                let b = &prepared[j / 3];
                if body_a == contacts[b.original].rigid_body {
                    let value = a.axes[i % 3].dot(b.axes[j % 3]) / rigid_properties[body_a].mass()
                        + a.torque_body[i % 3].dot(b.angular_response_body[j % 3]);
                    matrix[i * dimension + j] += value;
                    if i != j {
                        matrix[j * dimension + i] += value;
                    }
                }
            }
        }
        let contact = solve(&queries, &matrix).map_err(E::Contact)?;
        let physical = contact.physical();
        let mut rigid_next = rigid_bodies.to_vec();
        let mut solid_impulses = Vec::with_capacity(contacts.len());
        for (slot, p) in prepared.iter().enumerate() {
            poll::<R, CE>(cx)?;
            let c = contacts[p.original];
            let j = world(p.axes, physical.impulses[slot]);
            rigid_next[c.rigid_body] = rigid_next[c.rigid_body]
                .apply_impulse_at_body_point(rigid_properties[c.rigid_body], c.arm_body_m, j)?
                .state_after;
            solid_impulses.push((c.solid_point, array(j.scale(-1.0))));
        }
        // One combined solve retains all cross-contact work and support reactions.
        let (solid_next, solid_impulse) =
            self.point_impulses_candidate(&solid_impulses, origin_m, controls.solid_impulse, cx)?;
        let mut before = solid_impulse.kinetic_before_j;
        let mut after = solid_impulse.kinetic_after_j;
        let mut change = solid_impulse.kinetic_change_j;
        // Using the actual solid impulse identity avoids subtracting two large
        // total momenta. Its stationary support contribution cancels from the
        // externally balanced residual and remains available in the receipt.
        let mut linear_residual = vector(solid_impulse.applied.linear_kg_m_s)
            .add(vector(solid_impulse.linear_residual_kg_m_s));
        let mut angular_residual = vector(solid_impulse.applied.angular_kg_m2_s)
            .add(vector(solid_impulse.angular_residual_kg_m2_s));
        for ((&old, &new), &property) in rigid_bodies.iter().zip(&rigid_next).zip(rigid_properties)
        {
            poll::<R, CE>(cx)?;
            before += old.kinetic_energy(property)?;
            after += new.kinetic_energy(property)?;
            let dp = new.linear_momentum_world().sub(old.linear_momentum_world());
            let dl = new.angular_momentum_body().sub(old.angular_momentum_body());
            let v_average = old
                .center_of_mass_velocity_world(property)?
                .scale(0.5)
                .add(new.center_of_mass_velocity_world(property)?.scale(0.5));
            let w_average = property
                .angular_velocity_body_checked(old.angular_momentum_body())?
                .scale(0.5)
                .add(
                    property
                        .angular_velocity_body_checked(new.angular_momentum_body())?
                        .scale(0.5),
                );
            change += v_average.dot(dp) + w_average.dot(dl);
            linear_residual = linear_residual.add(dp);
            angular_residual = angular_residual
                .add(old.pose().position_world().sub(origin).cross(dp))
                .add(old.pose().orientation().rotate_body_to_world(dl));
        }
        let mut relative_work = 0.0;
        let mut endpoint_error = 0.0_f64;
        for (slot, p) in prepared.iter().enumerate() {
            poll::<R, CE>(cx)?;
            let c = contacts[p.original];
            let rigid_v = rigid_next[c.rigid_body]
                .point_kinematics(rigid_properties[c.rigid_body], c.arm_body_m)?
                .point_velocity_world;
            let solid_v = vector(solid_next.point_velocity_m_s(c.solid_point, cx)?);
            let relative = rigid_v.sub(solid_v);
            for axis in 0..3 {
                let velocity = relative.dot(p.axes[axis]);
                endpoint_error = endpoint_error
                    .max((velocity - physical.relative_velocities_after[slot][axis]).abs());
                relative_work += (0.5 * queries[slot].relative_velocity_m_s[axis] + 0.5 * velocity)
                    * physical.impulses[slot][axis];
            }
        }
        let contact_work = physical.relative_work;
        let report = FiniteSolidRigidContactReport {
            contact_ids: queries.iter().map(|c| c.contact_id).collect(),
            contact_mobility_per_kg: matrix,
            solid_mass_solves: mobility.solves,
            kinetic_before_j: before,
            kinetic_after_j: after,
            kinetic_change_j: change,
            relative_impulse_work_j: relative_work,
            work_energy_residual_j: change - relative_work,
            dissipation_residual_j: change + physical.dissipated_energy,
            linear_momentum_residual_kg_m_s: linear_residual,
            angular_momentum_residual_kg_m2_s: angular_residual,
            endpoint_velocity_residual_m_s: endpoint_error,
            contact,
            solid_impulse,
        };
        if [
            before,
            after,
            change,
            relative_work,
            endpoint_error,
            report.work_energy_residual_j,
            report.dissipation_residual_j,
            magnitude(linear_residual),
            magnitude(angular_residual),
        ]
        .iter()
        .any(|v| !v.is_finite())
            || endpoint_error > controls.velocity_tolerance_m_s
            || magnitude(linear_residual) > controls.linear_tolerance_kg_m_s
            || magnitude(angular_residual) > controls.angular_tolerance_kg_m2_s
            || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
            || report.dissipation_residual_j.abs() > controls.energy_tolerance_j
            || (relative_work - contact_work).abs() > controls.energy_tolerance_j
            || change > controls.energy_tolerance_j
        {
            return Err(E::Audit(Box::new(report)));
        }
        poll::<R, CE>(cx)?;
        Ok((solid_next, rigid_next, report))
    }
}

fn poll<R, E>(cx: &Cx<'_>) -> Result<(), FiniteSolidRigidContactError<R, E>> {
    cx.checkpoint()
        .map_err(|_| FiniteSolidRigidContactError::Cancelled)
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn unit(v: Vec3) -> Result<Vec3, &'static str> {
    let norm = magnitude(v);
    if !norm.is_finite() || norm <= 0.0 {
        return Err("finite independent contact axes required");
    }
    Ok(Vec3::new(v.x / norm, v.y / norm, v.z / norm))
}
fn contact_frame(normal: Vec3, tangent: Vec3) -> Result<[Vec3; 3], &'static str> {
    let n = unit(normal)?;
    let hint = unit(tangent)?;
    let cross = n.cross(hint);
    // Resolve the plane before normalizing it. A nearly parallel hint cannot
    // define a reliable Coulomb disk through projection-subtraction roundoff.
    if magnitude(cross) <= 128.0 * f64::EPSILON {
        return Err("unresolved contact tangent plane");
    }
    let b = unit(cross)?;
    let t = unit(b.cross(n))?;
    let axes = [n, t, b];
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            if (axes[i].dot(axes[j]) - expected).abs() > 128.0 * f64::EPSILON {
                return Err("unresolved orthonormal contact frame");
            }
        }
    }
    Ok(axes)
}
fn world(axes: [Vec3; 3], local: [f64; 3]) -> Vec3 {
    axes[0]
        .scale(local[0])
        .add(axes[1].scale(local[1]))
        .add(axes[2].scale(local[2]))
}
