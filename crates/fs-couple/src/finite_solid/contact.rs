//! Simultaneous frictional impacts among one deforming solid and rigid bodies.
//!
//! The complete contact mobility is B_s M_s,ff^-1 B_s^T plus the original
//! signed rigid-body mass/full-inertia response. Rigid/rigid and rigid/solid
//! contacts share one matrix, including every cross-point and shared-body block.
//! The existing fs-contact owner determines all normal and tangential
//! impulses together. Its determinacy and passivity refusals remain operative.
//!
//! Actual native velocities, spin/orbital momentum, stationary-support reaction
//! and impulse work are checked before returning either candidate. Geometry is
//! supplied to the core event; neither witness is relocated to hide a couple
//! from a finite gap. The `surface` owner detects current sphere/solid witnesses;
//! `scene` combines these with original rigid-pair witnesses. `step` composes
//! bounded transient trajectories and an explicit first-order force/contact/drift
//! lane for loaded support, with shared whole affine-path admission. No position
//! repair, general Moreau--Jean existence or impact-loss heat partition is
//! inferred. Material-point weights and rigid arms remain the original ones.

/// Directional interface sources through the same coupled mechanical event.
pub mod elliptic;

/// Mixed rigid/rigid and rigid/finite-solid events through one mechanical core.
pub mod mixed;

/// Automatic current witnesses and candidate paths for mixed sphere scenes.
pub mod scene;

/// Current finite-solid surface witnesses for native rigid contact.
pub mod surface;

/// Configuration-admitted coupled force/contact/trajectory candidates.
pub mod step;

/// Physical material directions transported by original rigid/solid geometry.
pub mod texture;

use super::impulse::{FiniteSolidImpulse, FiniteSolidImpulseControls, FiniteSolidPoint};
use super::{Cx, FiniteSolidBody, FiniteSolidError};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResult,
    NormalVelocityLaw, solve_frictional_contacts,
};
use fs_mbd::{DynamicsError, MassProperties, RigidBodyState, Vec3};
use fs_solver::SolveReport;
use mixed::{
    FiniteSolidContactCounterpart, FiniteSolidMixedContact, FiniteSolidMixedContactKinematics,
    RigidImpactPoint,
};

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
    /// Sum of native endpoint impulse errors; opposing body errors cannot cancel.
    pub absolute_linear_residual_kg_m_s: f64,
    /// Sum of native spin/orbital impulse errors, about the declared origin.
    pub absolute_angular_residual_kg_m2_s: f64,
    /// Sum of native actual-change minus commanded-impulse work magnitudes.
    pub absolute_work_residual_j: f64,
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
    rigid: [Option<RigidIncidence>; 2],
    solid_row: Option<usize>,
}

struct RigidIncidence {
    point: RigidImpactPoint,
    sign: f64,
    torque_body: [Vec3; 3],
    angular_response_body: [Vec3; 3],
}

impl RigidIncidence {
    fn new(
        point: RigidImpactPoint,
        sign: f64,
        axes: [Vec3; 3],
        body: RigidBodyState,
        property: MassProperties,
    ) -> Result<Self, DynamicsError> {
        let torque_body = axes.map(|axis| {
            point
                .arm_body
                .cross(body.pose().orientation().rotate_world_to_body(axis))
        });
        let mut angular_response_body = [Vec3::ZERO; 3];
        for axis in 0..3 {
            angular_response_body[axis] =
                property.angular_velocity_body_checked(torque_body[axis])?;
        }
        Ok(Self {
            point,
            sign,
            torque_body,
            angular_response_body,
        })
    }
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

    /// Compatibility adapter; all mechanical assembly and admission lives below.
    #[allow(clippy::too_many_arguments)]
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
        poll::<R, CE>(cx)?;
        if contacts.len() > controls.contact.maximum_contacts
            || controls.contact.maximum_contacts > 21
        {
            return Err(FiniteSolidRigidContactError::Invalid(
                "bounded contact count required",
            ));
        }
        let mixed: Vec<_> = contacts
            .iter()
            .copied()
            .map(FiniteSolidMixedContact::from)
            .collect();
        self.mixed_contact_candidate_with(
            rigid_bodies,
            rigid_properties,
            &mixed,
            origin_m,
            controls,
            cx,
            |queries, matrix| {
                let original: Vec<_> = queries
                    .iter()
                    .map(|q| FiniteSolidRigidContactKinematics {
                        original_index: q.original_index,
                        contact_id: q.contact_id,
                        axes_world: q.axes_world,
                        rigid_point_world: q.first_point_world,
                        solid_point_world: q.second_point_world,
                        gap_m: q.gap_m,
                        relative_velocity_m_s: q.relative_velocity_m_s,
                    })
                    .collect();
                solve(&original, matrix)
            },
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn mixed_contact_candidate_with<F: Copy, R: ImpactResponse, CE>(
        &self,
        rigid_bodies: &[RigidBodyState],
        rigid_properties: &[MassProperties],
        contacts: &[FiniteSolidMixedContact<F>],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        cx: &Cx<'_>,
        solve: impl FnOnce(&[FiniteSolidMixedContactKinematics], &[f64]) -> Result<R, CE>,
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
            if c.point.body >= rigid_bodies.len()
                || (slot > 0 && contacts[order[slot - 1]].id == c.id)
            {
                return Err(E::Invalid(
                    "valid rigid body and distinct contact IDs required",
                ));
            }
            let axes = contact_frame(c.normal_world, c.tangent_world).map_err(E::Invalid)?;
            let n = axes[0];
            let rigid = rigid_bodies[c.point.body];
            let property = rigid_properties[c.point.body];
            let kinematics = rigid.point_kinematics(property, c.point.arm_body)?;
            let first = RigidIncidence::new(c.point, 1.0, axes, rigid, property)?;
            let (second_point, second_velocity, second, solid_row) = match c.counterpart {
                FiniteSolidContactCounterpart::Solid(point) => {
                    let offset = rows.len();
                    rows.extend(axes.map(|axis| (point, array(axis))));
                    (
                        vector(self.point_position_m(point, cx)?),
                        vector(self.point_velocity_m_s(point, cx)?),
                        None,
                        Some(offset),
                    )
                }
                FiniteSolidContactCounterpart::Rigid(point) => {
                    if point.body >= rigid_bodies.len() || point.body == c.point.body {
                        return Err(E::Invalid(
                            "distinct valid rigid contact endpoints required",
                        ));
                    }
                    let body = rigid_bodies[point.body];
                    let property = rigid_properties[point.body];
                    let kinematics = body.point_kinematics(property, point.arm_body)?;
                    (
                        kinematics.point_world,
                        kinematics.point_velocity_world,
                        Some(RigidIncidence::new(point, -1.0, axes, body, property)?),
                        None,
                    )
                }
            };
            let separation = kinematics.point_world.sub(second_point);
            let gap = n.dot(separation);
            if magnitude(separation.sub(n.scale(gap))) > controls.geometry_tolerance_m {
                return Err(E::Invalid("contact witnesses do not share a normal line"));
            }
            let relative = kinematics.point_velocity_world.sub(second_velocity);
            queries.push(FiniteSolidMixedContactKinematics {
                original_index: index,
                contact_id: c.id,
                axes_world: axes,
                first_point_world: kinematics.point_world,
                second_point_world: second_point,
                gap_m: gap,
                relative_velocity_m_s: axes.map(|axis| relative.dot(axis)),
            });
            prepared.push(Prepared {
                original: index,
                axes,
                rigid: [Some(first), second],
                solid_row,
            });
        }
        let mobility = self.point_mobility_matrix_per_kg(
            &rows,
            3 * controls.contact.maximum_contacts,
            controls.solid_impulse.mass_solve,
            cx,
        )?;
        // The original consistent solid inverse is evaluated only for the
        // material-point rows. Rigid-only rows never cause extra mass solves.
        let solid_dimension = rows.len();
        let mut matrix = vec![0.0; dimension * dimension];
        for (i, a) in prepared.iter().enumerate() {
            if let Some(ar) = a.solid_row {
                for (j, b) in prepared.iter().enumerate() {
                    if let Some(br) = b.solid_row {
                        for ai in 0..3 {
                            for bj in 0..3 {
                                // Retain both native entries: this composition
                                // never repairs an asymmetric solid mass solve.
                                matrix[(3 * i + ai) * dimension + 3 * j + bj] =
                                    mobility.matrix_per_kg[(ar + ai) * solid_dimension + br + bj];
                            }
                        }
                    }
                }
            }
        }
        for i in 0..dimension {
            poll::<R, CE>(cx)?;
            let a = &prepared[i / 3];
            for j in i..dimension {
                let b = &prepared[j / 3];
                let mut value = 0.0;
                for ia in a.rigid.iter().flatten() {
                    for ib in b.rigid.iter().flatten() {
                        if ia.point.body == ib.point.body {
                            value += ia.sign
                                * ib.sign
                                * (a.axes[i % 3].dot(b.axes[j % 3])
                                    / rigid_properties[ia.point.body].mass()
                                    + ia.torque_body[i % 3].dot(ib.angular_response_body[j % 3]));
                        }
                    }
                }
                matrix[i * dimension + j] += value;
                if i != j {
                    matrix[j * dimension + i] += value;
                }
            }
        }
        let contact = solve(&queries, &matrix).map_err(E::Contact)?;
        let physical = contact.physical();
        let mut rigid_next = Vec::with_capacity(rigid_bodies.len());
        let mut linear_impulses = vec![Vec3::ZERO; rigid_bodies.len()];
        let mut angular_impulses = vec![Vec3::ZERO; rigid_bodies.len()];
        let mut solid_impulses = Vec::with_capacity(contacts.len());
        for (slot, p) in prepared.iter().enumerate() {
            poll::<R, CE>(cx)?;
            let c = contacts[p.original];
            let j = world(p.axes, physical.impulses[slot]);
            for incidence in p.rigid.iter().flatten() {
                let point = incidence.point;
                let impulse = j.scale(incidence.sign);
                let impulse_body = rigid_bodies[point.body]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(impulse);
                linear_impulses[point.body] = linear_impulses[point.body].add(impulse);
                angular_impulses[point.body] =
                    angular_impulses[point.body].add(point.arm_body.cross(impulse_body));
            }
            if let FiniteSolidContactCounterpart::Solid(point) = c.counterpart {
                solid_impulses.push((point, array(j.scale(-1.0))));
            }
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
        let mut absolute_linear = magnitude(vector(solid_impulse.linear_residual_kg_m_s));
        let mut absolute_angular = magnitude(vector(solid_impulse.angular_residual_kg_m2_s));
        let mut absolute_work = solid_impulse.work_residual_j.abs();
        for (index, (&old, &property)) in rigid_bodies.iter().zip(rigid_properties).enumerate() {
            poll::<R, CE>(cx)?;
            // One native momentum update per body retains all point moments
            // without repeatedly rounding its accepted momentum between rows.
            let receipt = old.apply_momentum_impulse(
                property,
                linear_impulses[index],
                angular_impulses[index],
            )?;
            rigid_next.push(receipt.state_after);
            before += receipt.work.kinetic_energy_before_j;
            after += receipt.work.kinetic_energy_after_j;
            change += receipt.work.kinetic_energy_change_j;
            let dp = receipt.linear_increment_world;
            let dl = receipt.angular_increment_body;
            let arm = old.pose().position_world().sub(origin);
            linear_residual = linear_residual.add(dp);
            angular_residual = angular_residual
                .add(arm.cross(dp))
                .add(old.pose().orientation().rotate_body_to_world(dl));
            absolute_linear += magnitude(receipt.linear_residual_world);
            absolute_angular += magnitude(arm.cross(receipt.linear_residual_world))
                + magnitude(
                    old.pose()
                        .orientation()
                        .rotate_body_to_world(receipt.angular_residual_body),
                );
            absolute_work += receipt.work.work_energy_residual_j.abs();
        }
        let mut relative_work = 0.0;
        let mut endpoint_error = 0.0_f64;
        for (slot, p) in prepared.iter().enumerate() {
            poll::<R, CE>(cx)?;
            let c = contacts[p.original];
            let first_velocity = rigid_next[c.point.body]
                .point_kinematics(rigid_properties[c.point.body], c.point.arm_body)?
                .point_velocity_world;
            let second_velocity = match c.counterpart {
                FiniteSolidContactCounterpart::Solid(point) => {
                    vector(solid_next.point_velocity_m_s(point, cx)?)
                }
                FiniteSolidContactCounterpart::Rigid(point) => {
                    rigid_next[point.body]
                        .point_kinematics(rigid_properties[point.body], point.arm_body)?
                        .point_velocity_world
                }
            };
            let relative = first_velocity.sub(second_velocity);
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
            absolute_linear_residual_kg_m_s: absolute_linear,
            absolute_angular_residual_kg_m2_s: absolute_angular,
            absolute_work_residual_j: absolute_work,
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
            absolute_linear,
            absolute_angular,
            absolute_work,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || endpoint_error > controls.velocity_tolerance_m_s
            || magnitude(linear_residual) > controls.linear_tolerance_kg_m_s
            || magnitude(angular_residual) > controls.angular_tolerance_kg_m2_s
            || absolute_linear > controls.linear_tolerance_kg_m_s
            || absolute_angular > controls.angular_tolerance_kg_m2_s
            || absolute_work > controls.energy_tolerance_j
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
