//! Simultaneous retained scalar Mindlin histories on shared native rigid bodies.
//!
//! Each source owns its fixed Hertz load, temperature, material ordering and
//! scalar shear history. Actual masses, full inertias and arms determine every
//! cross-contact longitudinal response. The original simultaneous tangential
//! solver advances all histories from the same accepted checkpoint; a sweep is
//! never another physical time step. Held COM forces and body torques enter
//! that solve, and every body receives one summed contact/external resultant.
//!
//! This is an explicitly scalar, fixed-normal, frozen-geometry momentum phase.
//! The complete physical two-tangent mobility is retained for inspection. Only
//! its longitudinal subspace enters the scalar constitutive solve, with an
//! exactly inactive second bookkeeping coordinate. Actual native endpoints and
//! midpoint displacements must satisfy the omitted normal/lateral/torsional
//! motion gates; no physical recoil is projected away to obtain acceptance.
//! Normal prestress is supplied by the source and surrounding held-normal
//! owner, with no normal impulse or work in this shear phase. Shared-body patch
//! centers must be farther apart than the sum of their sourced Hertz radii;
//! this conservative footprint admission prevents counting one loaded area
//! twice. Geometry/pose transport and normal equilibration remain with the
//! caller; inter-patch bulk elastic fields and rolling-memory transport are
//! not inferred from rigid-body coupling.

use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3, Wrench};
use fs_contact::interface_binding::mindlin::{SourceMindlinContact, SourceMindlinState};
use fs_contact::normal_patch::NormalPatchError;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::multiple::{
    MultipleTangentialConfig, MultipleTangentialError, PreparedTangentialContacts,
    TangentialContactInput,
};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput};

/// Coupled normal dynamics and retained shear on shared native bodies.
#[path = "multiple/loading.rs"]
pub mod loading;

/// Original ordered source/history and attachments for one independent patch.
#[derive(Clone, Copy, Debug)]
pub struct RigidMindlinIslandContact<'a> {
    pub contact_id: &'a str,
    pub law: &'a SourceMindlinContact,
    pub history: &'a SourceMindlinState,
    /// Original whole-array A/B body indices; their material order is retained.
    pub body_indices: [usize; 2],
    /// COM-relative material points in each body's original frame.
    pub arms_body: [Vec3; 2],
    /// Positive A-minus-B longitudinal motion increases the retained shear.
    pub frame: TangentFrame,
}

/// Explicit original-law/global-solve budgets and native dimensional gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidMindlinIslandControls {
    pub solve: TangentialSolveConfig,
    pub multiple: MultipleTangentialConfig,
    /// Participating native bodies, at most 256; spectators are not copied.
    pub maximum_bodies: usize,
    pub maximum_step_s: f64,
    pub maximum_slip_speed_m_s: f64,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Actual native motion and work, in canonical original contact-ID order.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinIslandContactReport {
    /// Longitudinal, lateral, normal components of A-minus-B velocity.
    pub relative_velocity_before_m_s: [f64; 3],
    pub relative_velocity_after_m_s: [f64; 3],
    pub relative_torsion_before_rad_s: f64,
    pub relative_torsion_after_rad_s: f64,
    pub longitudinal_driving_acceleration_m_s2: f64,
    /// Actual midpoint displacement in the same three directions.
    pub relative_displacement_m: [f64; 3],
    pub impulse_on_b_world_ns: Vec3,
    pub midpoint_impulse_work_j: f64,
    /// Original storage change plus heat plus actual contact work on bodies.
    pub contact_energy_residual_j: f64,
    pub compatibility_residual_m: f64,
}

/// One accepted-source candidate; neither source nor history is refreshed.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinIslandContactCandidate {
    pub contact_id: String,
    pub contact: TangentialStep<SourceMindlinState>,
    pub report: RigidMindlinIslandContactReport,
}

/// Shared-body accounting counts each native body and external wrench once.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinIslandReport {
    pub sweeps: usize,
    pub law_evaluations: usize,
    pub maximum_compatibility_residual_m: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Actual native increments contracted with midpoint velocities.
    pub kinetic_change_j: f64,
    pub stored_change_j: f64,
    pub dissipated_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub external_work_j: f64,
    pub work_energy_residual_j: f64,
    /// Delta(K + retained shear storage) + physical heat - external work.
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub external_linear_impulse_ns: Vec3,
    pub external_angular_impulse_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// Publish all original body and source-bound history candidates together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinIslandCandidate {
    /// Sorted original global body indices, paired with `bodies_after`.
    pub body_indices: Vec<usize>,
    pub bodies_after: Vec<RigidBodyState>,
    pub contacts: Vec<RigidMindlinIslandContactCandidate>,
    pub report: RigidMindlinIslandReport,
}

#[derive(Debug)]
pub enum RigidMindlinIslandError {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Normal(NormalPatchError),
    Native(MultipleTangentialError),
    /// Original physical native motion is outside the scalar held-normal phase.
    UnsupportedMotion {
        contact_id: String,
        coordinate: &'static str,
    },
    Unresolved(&'static str),
    Cancelled,
}
impl From<DynamicsError> for RigidMindlinIslandError {
    fn from(value: DynamicsError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<NormalPatchError> for RigidMindlinIslandError {
    fn from(value: NormalPatchError) -> Self {
        Self::Normal(value)
    }
}
impl From<MultipleTangentialError> for RigidMindlinIslandError {
    fn from(value: MultipleTangentialError) -> Self {
        match value {
            MultipleTangentialError::Cancelled => Self::Cancelled,
            value => Self::Native(value),
        }
    }
}
impl core::fmt::Display for RigidMindlinIslandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native fixed-normal Mindlin island: {self:?}")
    }
}
impl std::error::Error for RigidMindlinIslandError {}

#[derive(Clone, Debug)]
struct Patch<'a> {
    id: String,
    law: &'a SourceMindlinContact,
    history: SourceMindlinState,
    bodies: [usize; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
    before: [PointKinematics; 2],
    rates: [f64; 3],
    torsion: f64,
    drive: f64,
    released: f64,
    stored: f64,
}

/// Captured native geometry, original sources/histories and held external drive.
#[derive(Clone, Debug)]
pub struct PreparedRigidMindlinIsland<'a> {
    patches: Vec<Patch<'a>>,
    body_indices: Vec<usize>,
    before: Vec<RigidBodyState>,
    properties: Vec<MassProperties>,
    external: Vec<Wrench>,
    mobility: Vec<f64>,
    native: PreparedTangentialContacts,
    dt_s: f64,
    controls: RigidMindlinIslandControls,
}

impl<'a> PreparedRigidMindlinIsland<'a> {
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidMindlinIslandContact<'a>],
        dt_s: f64,
        controls: RigidMindlinIslandControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RigidMindlinIslandError> {
        Self::prepare(
            bodies, properties, contacts, None, dt_s, controls, cancelled,
        )
    }

    /// External COM forces are world-frame vectors; torques use body frames.
    /// The slice follows the entire original body array. A nonzero wrench on a
    /// nonparticipating body is refused, never silently dropped from work.
    #[allow(clippy::too_many_arguments)]
    pub fn new_driven(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidMindlinIslandContact<'a>],
        external: &[Wrench],
        dt_s: f64,
        controls: RigidMindlinIslandControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RigidMindlinIslandError> {
        Self::prepare(
            bodies,
            properties,
            contacts,
            Some(external),
            dt_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidMindlinIslandContact<'a>],
        external: Option<&[Wrench]>,
        dt_s: f64,
        controls: RigidMindlinIslandControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RigidMindlinIslandError> {
        poll(cancelled)?;
        controls.validate(dt_s)?;
        if bodies.len() != properties.len()
            || contacts.is_empty()
            || contacts.len() > controls.multiple.maximum_contacts
            || 4 * contacts.len() * contacts.len() > controls.multiple.maximum_matrix_terms
            || external.is_some_and(|w| w.len() != bodies.len())
        {
            return Err(invalid("body/contact shape or physical matrix budget"));
        }
        let mut ordered = contacts.to_vec();
        ordered.sort_by_key(|contact| contact.contact_id);
        let mut body_indices = Vec::with_capacity(2 * contacts.len());
        for (i, contact) in ordered.iter().enumerate() {
            poll(cancelled)?;
            if contact.contact_id.trim().is_empty()
                || (i > 0 && contact.contact_id == ordered[i - 1].contact_id)
                || contact.body_indices[0] == contact.body_indices[1]
                || contact.body_indices.iter().any(|&b| b >= bodies.len())
            {
                return Err(invalid(
                    "unique original contact ID and distinct native bodies required",
                ));
            }
            body_indices.extend(contact.body_indices);
        }
        body_indices.sort_unstable();
        body_indices.dedup();
        if body_indices.len() > controls.maximum_bodies {
            return Err(invalid("participating native body budget"));
        }
        if let Some(external) = external {
            for (index, &wrench) in external.iter().enumerate() {
                poll(cancelled)?;
                wrench.validate()?;
                if wrench != Wrench::ZERO && body_indices.binary_search(&index).is_err() {
                    return Err(invalid(
                        "external wrench on a body outside this contact island",
                    ));
                }
            }
        }
        let before: Vec<_> = body_indices.iter().map(|&i| bodies[i]).collect();
        let properties: Vec<_> = body_indices.iter().map(|&i| properties[i]).collect();
        let external: Vec<_> = body_indices
            .iter()
            .map(|&i| external.map_or(Wrench::ZERO, |w| w[i]))
            .collect();
        let mut patches = Vec::with_capacity(ordered.len());
        for contact in ordered {
            poll(cancelled)?;
            let indices = contact
                .body_indices
                .map(|i| body_indices.binary_search(&i).unwrap());
            let points = [0, 1].map(|i| {
                before[indices[i]].point_kinematics(properties[indices[i]], contact.arms_body[i])
            });
            let [a, b] = points;
            let points = [a?, b?];
            if norm(points[0].point_world.sub(points[1].point_world))?
                > controls.geometry_tolerance_m
            {
                return Err(invalid("original contact points do not coincide"));
            }
            let axes = [
                contact.frame.longitudinal(),
                contact.frame.lateral(),
                contact.frame.normal(),
            ]
            .map(vector);
            let (rates, torsion) = motion(points, axes)?;
            admit_motion(
                contact.contact_id,
                contact.law,
                rates,
                torsion,
                dt_s,
                controls,
            )?;
            let stored = number(contact.law.stored_energy_j(contact.history))?;
            if stored < 0.0 {
                return Err(invalid("original source-bound shear history required"));
            }
            let mut acceleration = [Vec3::ZERO; 2];
            for i in 0..2 {
                let body = indices[i];
                let alpha = before[body].pose().orientation().rotate_body_to_world(
                    properties[body].angular_velocity_body_checked(external[body].torque_body)?,
                );
                acceleration[i] = external[body]
                    .force_world
                    .scale(properties[body].mass().recip())
                    .add(alpha.cross(points[i].arm_world));
                norm(acceleration[i])?;
            }
            let drive = number(acceleration[0].sub(acceleration[1]).dot(axes[0]))?;
            let released = number(dt_s * rates[0] + (0.5 * dt_s) * (dt_s * drive))?;
            patches.push(Patch {
                id: contact.contact_id.into(),
                law: contact.law,
                history: contact.history.clone(),
                bodies: indices,
                arms: contact.arms_body,
                axes,
                before: points,
                rates,
                torsion,
                drive,
                released,
                stored,
            });
        }
        for (i, patch) in patches.iter().enumerate() {
            poll(cancelled)?;
            for other in &patches[..i] {
                let clearance = number(
                    patch.law.normal_response().patch_radius_m
                        + other.law.normal_response().patch_radius_m
                        + controls.geometry_tolerance_m,
                )?;
                for a in 0..2 {
                    for b in 0..2 {
                        if patch.bodies[a] == other.bodies[b]
                            && norm(patch.arms[a].sub(other.arms[b]))? <= clearance
                        {
                            return Err(invalid(
                                "independent shared-body Hertz footprints are not separated",
                            ));
                        }
                    }
                }
            }
        }
        let dimension = 2 * patches.len();
        let mut factors = Vec::with_capacity(dimension);
        for patch in &patches {
            for &axis in &patch.axes[..2] {
                let mut row = [(Vec3::ZERO, Vec3::ZERO); 2];
                for i in 0..2 {
                    let body = patch.bodies[i];
                    let direction = axis.scale(if i == 0 { 1.0 } else { -1.0 });
                    row[i] = (
                        direction.scale(1.0 / properties[body].mass().sqrt()),
                        properties[body].angular_mobility_factor(
                            patch.arms[i].cross(
                                before[body]
                                    .pose()
                                    .orientation()
                                    .rotate_world_to_body(direction),
                            ),
                        ),
                    );
                    norm(row[i].0)?;
                    norm(row[i].1)?;
                }
                factors.push(row);
            }
        }
        let mut mobility = vec![0.0; dimension * dimension];
        let mut compliance = vec![0.0; dimension * dimension];
        for row in 0..dimension {
            poll(cancelled)?;
            for column in 0..=row {
                let mut value = 0.0;
                for i in 0..2 {
                    for j in 0..2 {
                        if patches[row / 2].bodies[i] == patches[column / 2].bodies[j] {
                            value += factors[row][i].0.dot(factors[column][j].0)
                                + factors[row][i].1.dot(factors[column][j].1);
                        }
                    }
                }
                let value = number(value)?;
                mobility[row * dimension + column] = value;
                mobility[column * dimension + row] = value;
                if row % 2 == 0 && column % 2 == 0 {
                    // Exact supported scalar subspace, not a replacement for
                    // the physical lateral rows retained above and audited below.
                    let c = number((0.5 * dt_s) * (dt_s * value))?;
                    compliance[row * dimension + column] = c;
                    compliance[column * dimension + row] = c;
                }
            }
        }
        let native = PreparedTangentialContacts::new(
            patches.len(),
            compliance,
            controls.multiple,
            &mut *cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            patches,
            body_indices,
            before,
            properties,
            external,
            mobility,
            native,
            dt_s,
            controls,
        })
    }

    /// Full unchanged physical 2N-by-2N mobility, longitudinal/lateral per site
    /// in canonical contact-ID order, including every off-diagonal body response.
    pub fn tangential_mobility_kg_inverse(&self) -> &[f64] {
        &self.mobility
    }

    pub fn body_indices(&self) -> &[usize] {
        &self.body_indices
    }

    /// Resolve one immutable trial from all captured original accepted histories.
    /// Retrying after any refusal or cancellation starts from those same states.
    pub fn advance(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RigidMindlinIslandCandidate, RigidMindlinIslandError> {
        poll(cancelled)?;
        let inputs: Vec<_> = self
            .patches
            .iter()
            .map(|p| TangentialContactInput {
                law: p.law,
                state: &p.history,
                input: TangentialStepInput {
                    normal_n: p.law.normal_response().normal_force_n,
                    slip_increment_m: [p.released, 0.0],
                    dt_s: self.dt_s,
                    contact_temperature_k: Some(p.law.temperature_kelvin()),
                    nominal_pressure_pa: None,
                    compliance_m_per_n: 0.0,
                    solve: self.controls.solve,
                },
            })
            .collect();
        let native = self.native.advance(&inputs, &mut *cancelled)?;
        let count = self.before.len();
        let mut linear = vec![Vec3::ZERO; count];
        let mut angular = vec![Vec3::ZERO; count];
        let mut impulses = Vec::with_capacity(self.patches.len());
        for (patch, step) in self.patches.iter().zip(&native.contacts) {
            poll(cancelled)?;
            if step.traction_n[1] != 0.0 || step.slip_increment_m[1] != 0.0 {
                return Err(invalid("scalar source produced a transverse response"));
            }
            let impulse = patch.axes[0].scale(number(self.dt_s * step.traction_n[0])?);
            norm(impulse)?;
            impulses.push(impulse);
            for i in 0..2 {
                let body = patch.bodies[i];
                let signed = impulse.scale(if i == 0 { -1.0 } else { 1.0 });
                linear[body] = linear[body].add(signed);
                angular[body] = angular[body].add(
                    patch.arms[i].cross(
                        self.before[body]
                            .pose()
                            .orientation()
                            .rotate_world_to_body(signed),
                    ),
                );
            }
        }
        let mut after = Vec::with_capacity(count);
        let (mut kinetic_before, mut kinetic_after, mut kinetic_change, mut external_work) =
            (0.0, 0.0, 0.0, 0.0);
        let (mut delta_linear, mut delta_angular) = (Vec3::ZERO, Vec3::ZERO);
        let (mut external_linear, mut external_angular) = (Vec3::ZERO, Vec3::ZERO);
        for body in 0..count {
            poll(cancelled)?;
            let state = self.before[body];
            let prop = self.properties[body];
            let wrench = self.external[body];
            let ext_linear = wrench.force_world.scale(self.dt_s);
            let ext_angular = wrench.torque_body.scale(self.dt_s);
            let receipt = state.apply_momentum_impulse(
                prop,
                linear[body].add(ext_linear),
                angular[body].add(ext_angular),
            )?;
            audit(
                norm(receipt.linear_residual_world)?,
                self.controls.linear_momentum_tolerance_ns,
                "native per-body resultant linear impulse",
            )?;
            audit(
                norm(receipt.angular_residual_body)?,
                self.controls.angular_momentum_tolerance_nms,
                "native per-body resultant angular impulse",
            )?;
            let mean_v = state.center_of_mass_velocity_world(prop)?.scale(0.5).add(
                receipt
                    .state_after
                    .center_of_mass_velocity_world(prop)?
                    .scale(0.5),
            );
            let mean_w = prop
                .angular_velocity_body_checked(state.angular_momentum_body())?
                .scale(0.5)
                .add(
                    prop.angular_velocity_body_checked(
                        receipt.state_after.angular_momentum_body(),
                    )?
                    .scale(0.5),
                );
            external_work =
                number(external_work + ext_linear.dot(mean_v) + ext_angular.dot(mean_w))?;
            kinetic_before = number(kinetic_before + receipt.work.kinetic_energy_before_j)?;
            kinetic_after = number(kinetic_after + receipt.work.kinetic_energy_after_j)?;
            kinetic_change = number(kinetic_change + receipt.work.kinetic_energy_change_j)?;
            delta_linear = delta_linear.add(receipt.linear_increment_world);
            delta_angular = delta_angular
                .add(
                    state
                        .pose()
                        .position_world()
                        .cross(receipt.linear_increment_world),
                )
                .add(
                    state
                        .pose()
                        .orientation()
                        .rotate_body_to_world(receipt.angular_increment_body),
                );
            external_linear = external_linear.add(ext_linear);
            external_angular = external_angular
                .add(state.pose().position_world().cross(ext_linear))
                .add(state.pose().orientation().rotate_body_to_world(ext_angular));
            after.push(receipt.state_after);
        }
        let mut contacts = Vec::with_capacity(self.patches.len());
        let (mut stored, mut heat, mut contact_work) = (0.0, 0.0, 0.0);
        let mut maximum_compatibility = native.maximum_compatibility_residual_m;
        for ((patch, step), impulse) in self.patches.iter().zip(native.contacts).zip(impulses) {
            poll(cancelled)?;
            let points = [0, 1].map(|i| {
                after[patch.bodies[i]]
                    .point_kinematics(self.properties[patch.bodies[i]], patch.arms[i])
            });
            let [a, b] = points;
            let points = [a?, b?];
            let (rates, torsion) = motion(points, patch.axes)?;
            admit_motion(
                &patch.id,
                patch.law,
                rates,
                torsion,
                self.dt_s,
                self.controls,
            )?;
            let travel =
                [0, 1, 2].map(|axis| self.dt_s * (0.5 * patch.rates[axis] + 0.5 * rates[axis]));
            for (axis, name) in [(1, "lateral displacement"), (2, "normal displacement")] {
                if number(travel[axis])?.abs()
                    > self.controls.multiple.local.displacement_tolerance_m
                {
                    return Err(RigidMindlinIslandError::UnsupportedMotion {
                        contact_id: patch.id.clone(),
                        coordinate: name,
                    });
                }
            }
            let compatibility = number(travel[0] - step.slip_increment_m[0])?;
            audit(
                compatibility,
                self.controls.multiple.local.displacement_tolerance_m,
                "actual native longitudinal displacement",
            )?;
            maximum_compatibility = maximum_compatibility.max(compatibility.abs());
            let mean_relative = patch.before[0]
                .point_velocity_world
                .sub(patch.before[1].point_velocity_world)
                .scale(0.5)
                .add(
                    points[0]
                        .point_velocity_world
                        .sub(points[1].point_velocity_world)
                        .scale(0.5),
                );
            let work = number(-impulse.dot(mean_relative))?;
            let change = number(patch.law.stored_energy_j(&step.state) - patch.stored)?;
            audit(
                number(change - step.stored_change_j)?,
                self.controls.multiple.local.energy_tolerance_j,
                "original retained shear storage",
            )?;
            let residual = number(change + step.dissipated_j + work)?;
            audit(
                residual,
                self.controls.multiple.local.energy_tolerance_j,
                "actual native contact work/storage/heat",
            )?;
            stored = number(stored + change)?;
            heat = number(heat + step.dissipated_j)?;
            contact_work = number(contact_work + work)?;
            contacts.push(RigidMindlinIslandContactCandidate {
                contact_id: patch.id.clone(),
                contact: step,
                report: RigidMindlinIslandContactReport {
                    relative_velocity_before_m_s: patch.rates,
                    relative_velocity_after_m_s: rates,
                    relative_torsion_before_rad_s: patch.torsion,
                    relative_torsion_after_rad_s: torsion,
                    longitudinal_driving_acceleration_m_s2: patch.drive,
                    relative_displacement_m: travel,
                    impulse_on_b_world_ns: impulse,
                    midpoint_impulse_work_j: work,
                    contact_energy_residual_j: residual,
                    compatibility_residual_m: compatibility,
                },
            });
        }
        let linear_residual = delta_linear.sub(external_linear);
        let angular_residual = delta_angular.sub(external_angular);
        audit(
            norm(linear_residual)?,
            self.controls.linear_momentum_tolerance_ns,
            "island linear momentum",
        )?;
        audit(
            norm(angular_residual)?,
            self.controls.angular_momentum_tolerance_nms,
            "island angular momentum",
        )?;
        let work_residual = number(kinetic_change - contact_work - external_work)?;
        let energy_residual = number(kinetic_change + stored + heat - external_work)?;
        audit(
            work_residual,
            self.controls.energy_tolerance_j,
            "actual full-body work/kinetic change",
        )?;
        audit(
            energy_residual,
            self.controls.energy_tolerance_j,
            "full-body storage/heat/external work",
        )?;
        poll(cancelled)?;
        Ok(RigidMindlinIslandCandidate {
            body_indices: self.body_indices.clone(),
            bodies_after: after,
            contacts,
            report: RigidMindlinIslandReport {
                sweeps: native.sweeps,
                law_evaluations: native.law_evaluations,
                maximum_compatibility_residual_m: maximum_compatibility,
                kinetic_before_j: kinetic_before,
                kinetic_after_j: kinetic_after,
                kinetic_change_j: kinetic_change,
                stored_change_j: stored,
                dissipated_j: heat,
                midpoint_impulse_work_j: contact_work,
                external_work_j: external_work,
                work_energy_residual_j: work_residual,
                total_energy_residual_j: energy_residual,
                linear_momentum_change_ns: delta_linear,
                angular_momentum_change_nms: delta_angular,
                external_linear_impulse_ns: external_linear,
                external_angular_impulse_nms: external_angular,
                linear_momentum_residual_ns: linear_residual,
                angular_momentum_residual_nms: angular_residual,
            },
        })
    }
}

impl RigidMindlinIslandControls {
    fn validate(self, dt: f64) -> Result<(), RigidMindlinIslandError> {
        self.multiple.validate()?;
        if !(1..=256).contains(&self.maximum_bodies)
            || !(1..=128).contains(&self.solve.max_iterations)
            || [
                dt,
                self.maximum_step_s,
                self.maximum_slip_speed_m_s,
                self.geometry_tolerance_m,
                self.velocity_tolerance_m_s,
                self.angular_velocity_tolerance_rad_s,
                self.linear_momentum_tolerance_ns,
                self.angular_momentum_tolerance_nms,
                self.energy_tolerance_j,
                self.solve.maximum_force_n,
                self.solve.force_absolute_tolerance_n,
                self.solve.force_relative_tolerance,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.0)
            || self.solve.force_relative_tolerance >= 1.0
            || dt > self.maximum_step_s
        {
            return Err(invalid(
                "finite positive native gates and original solve budgets required",
            ));
        }
        Ok(())
    }
}

fn motion(
    points: [PointKinematics; 2],
    axes: [Vec3; 3],
) -> Result<([f64; 3], f64), RigidMindlinIslandError> {
    let relative = points[0]
        .point_velocity_world
        .sub(points[1].point_velocity_world);
    let mut rates = [0.0; 3];
    for i in 0..3 {
        rates[i] = number(relative.dot(axes[i]))?;
    }
    let torsion = number(
        points[0]
            .angular_velocity_world
            .sub(points[1].angular_velocity_world)
            .dot(axes[2]),
    )?;
    Ok((rates, torsion))
}
fn admit_motion(
    id: &str,
    law: &SourceMindlinContact,
    rates: [f64; 3],
    torsion: f64,
    dt: f64,
    controls: RigidMindlinIslandControls,
) -> Result<(), RigidMindlinIslandError> {
    law.normal_source()
        .validate_kinematics(law.normal_response().approach_m, rates[2], dt)?;
    for (value, tolerance, coordinate) in [
        (
            rates[1],
            controls.velocity_tolerance_m_s,
            "lateral velocity",
        ),
        (rates[2], controls.velocity_tolerance_m_s, "normal velocity"),
        (
            torsion,
            controls.angular_velocity_tolerance_rad_s,
            "relative torsion",
        ),
    ] {
        if value.abs() > tolerance {
            return Err(RigidMindlinIslandError::UnsupportedMotion {
                contact_id: id.into(),
                coordinate,
            });
        }
    }
    if rates[0].abs() > controls.maximum_slip_speed_m_s {
        return Err(invalid(
            "actual native shear speed exceeds caller admission",
        ));
    }
    Ok(())
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn invalid(what: &'static str) -> RigidMindlinIslandError {
    RigidMindlinIslandError::Invalid(what)
}
fn number(value: f64) -> Result<f64, RigidMindlinIslandError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable native island arithmetic"))
    }
}
fn norm(value: Vec3) -> Result<f64, RigidMindlinIslandError> {
    Ok(value.stable_norm("native Mindlin island vector")?)
}
fn audit(value: f64, allowance: f64, what: &'static str) -> Result<(), RigidMindlinIslandError> {
    if number(value)?.abs() <= allowance {
        Ok(())
    } else {
        Err(RigidMindlinIslandError::Unresolved(what))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RigidMindlinIslandError> {
    if cancelled() {
        Err(RigidMindlinIslandError::Cancelled)
    } else {
        Ok(())
    }
}
