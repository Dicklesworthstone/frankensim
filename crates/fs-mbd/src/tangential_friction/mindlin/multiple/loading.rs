//! Shared-body Hertz normal dynamics with simultaneous retained Mindlin shear.
//!
//! The physical normal/longitudinal/lateral Gram response is assembled from
//! actual native masses, full inertias and arms. Normal and shear blocks must
//! be mechanically orthogonal to arithmetic precision; all cross-contact terms
//! within either supported block remain. The original Newton--Krylov owner
//! resolves the global normal load line using exact source Hertz averages and
//! derivatives. Original simultaneous tangential solves then consume those SAME
//! affine indentation paths, including unloading-induced microslip.
//!
//! Every body receives one summed normal/shear/external impulse. Actual native
//! motion and separate normal potential, shear storage, physical heat and drive
//! work must close together. Only continuously compressed, nonoverlapping
//! source patches are admitted. Geometry, orientations, material conditions and
//! drives are held during this momentum phase; separation/contact events, pose
//! integration, bulk patch interaction, rotating-memory transport and genuinely
//! coupled normal/shear mechanical blocks require their surrounding owners.

use super::{
    RigidMindlinIslandControls, RigidMindlinIslandError, audit, norm, number, poll, vector,
};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3, Wrench};
use fs_contact::interface_binding::mindlin::loading::{
    SourceHertzMindlinContact, SourceHertzMindlinState,
    path::{SourceMindlinPathError, SourceMindlinPathInput, SourceMindlinPathStep},
};
use fs_contact::normal_patch::{NormalPatchError, PointNormalPatchReceipt};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonReport};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::multiple::{
    MultipleTangentialError, PreparedTangentialContacts, TangentialContactInput,
};
use fs_tribo::{TangentialLaw, TangentialStep, TangentialStepInput, TriboError};
use std::cell::RefCell;

#[path = "loading/normal.rs"]
mod normal;

/// Original ordered changing-load source/history and physical attachments.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzMindlinIslandContact<'a> {
    pub contact_id: &'a str,
    pub law: &'a SourceHertzMindlinContact,
    pub history: &'a SourceHertzMindlinState,
    /// Original whole-array A/B body indices; positive approach is A minus B.
    pub body_indices: [usize; 2],
    pub arms_body: [Vec3; 2],
    /// Normal points A toward B; the longitudinal axis carries scalar memory.
    pub frame: TangentFrame,
}

/// Existing shear/native controls and explicit original Newton work ceilings.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzMindlinIslandControls {
    /// The matrix-entry ceiling also covers the full physical 3N-by-3N response.
    pub native: RigidMindlinIslandControls,
    /// Newton state and residual are final indentation and load-line error [m].
    /// Actual physical closure is additionally checked against native gates.
    pub normal_solver: NewtonKrylovConfig,
    /// Original Newton attempts, 1 through 128.
    pub maximum_normal_iterations: usize,
    /// Total original Hertz path evaluations, including rejected trials.
    /// Between 1 and 1,048,576; Jacobian products reuse exact cached derivatives.
    pub maximum_normal_evaluations: usize,
}

/// One original contact candidate and independently measured native work.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzMindlinIslandContactCandidate {
    pub contact_id: String,
    pub contact: SourceMindlinPathStep,
    pub normal_average_force_n: f64,
    /// Normal, longitudinal, lateral components of A-minus-B point velocity.
    pub relative_velocity_before_m_s: [f64; 3],
    pub relative_velocity_after_m_s: [f64; 3],
    pub relative_torsion_before_rad_s: f64,
    pub relative_torsion_after_rad_s: f64,
    /// Normal and longitudinal acceleration from original external wrenches.
    pub driving_acceleration_m_s2: [f64; 2],
    pub relative_displacement_m: [f64; 3],
    pub impulse_on_b_world_ns: Vec3,
    pub midpoint_impulse_work_j: f64,
    pub contact_energy_residual_j: f64,
    pub compatibility_residual_m: [f64; 2],
}

/// Separate physical energy and actual native world-momentum balances.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzMindlinIslandReport {
    pub normal_iterations: usize,
    pub normal_source_evaluations: usize,
    pub maximum_normal_residual_m: f64,
    pub shear_sweeps: usize,
    pub shear_law_evaluations: usize,
    pub maximum_shear_residual_m: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Actual momentum increments contracted with native midpoint velocities.
    pub kinetic_change_j: f64,
    pub normal_stored_change_j: f64,
    pub shear_stored_change_j: f64,
    pub dissipated_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub external_work_j: f64,
    pub work_energy_residual_j: f64,
    /// Delta(K + normal storage + shear storage) + heat - external work.
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub external_linear_impulse_ns: Vec3,
    pub external_angular_impulse_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// All native body and original contact histories are accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzMindlinIslandCandidate {
    pub body_indices: Vec<usize>,
    pub bodies_after: Vec<RigidBodyState>,
    /// Canonical original contact-ID order.
    pub contacts: Vec<RigidHertzMindlinIslandContactCandidate>,
    pub report: RigidHertzMindlinIslandReport,
}

#[derive(Debug)]
pub enum RigidHertzMindlinIslandError {
    Mechanics(RigidMindlinIslandError),
    Newton(NewtonError),
    NormalNotConverged {
        report: NewtonReport,
        last_source_refusal: Option<NormalPatchError>,
    },
    NormalEvaluationBudget {
        evaluations: usize,
    },
    SourcePath {
        contact_id: String,
        error: SourceMindlinPathError,
    },
    Cancelled,
}
impl From<RigidMindlinIslandError> for RigidHertzMindlinIslandError {
    fn from(error: RigidMindlinIslandError) -> Self {
        match error {
            RigidMindlinIslandError::Cancelled => Self::Cancelled,
            error => Self::Mechanics(error),
        }
    }
}
impl From<DynamicsError> for RigidHertzMindlinIslandError {
    fn from(error: DynamicsError) -> Self {
        RigidMindlinIslandError::from(error).into()
    }
}
impl From<NormalPatchError> for RigidHertzMindlinIslandError {
    fn from(error: NormalPatchError) -> Self {
        RigidMindlinIslandError::from(error).into()
    }
}
impl From<MultipleTangentialError> for RigidHertzMindlinIslandError {
    fn from(error: MultipleTangentialError) -> Self {
        RigidMindlinIslandError::from(error).into()
    }
}
impl From<NewtonError> for RigidHertzMindlinIslandError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl core::fmt::Display for RigidHertzMindlinIslandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native changing-normal Mindlin island: {self:?}")
    }
}
impl std::error::Error for RigidHertzMindlinIslandError {}

#[derive(Clone, Debug)]
struct Patch<'a> {
    id: String,
    law: &'a SourceHertzMindlinContact,
    history: SourceHertzMindlinState,
    bodies: [usize; 2],
    arms: [Vec3; 2],
    /// Normal, longitudinal, lateral world axes.
    axes: [Vec3; 3],
    before: [PointKinematics; 2],
    rates: [f64; 3],
    torsion: f64,
    drive: [f64; 2],
    released: [f64; 2],
    stored: f64,
}

/// Captured original sources, accepted histories, native geometry and drive.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzMindlinIsland<'a> {
    patches: Vec<Patch<'a>>,
    body_indices: Vec<usize>,
    before: Vec<RigidBodyState>,
    properties: Vec<MassProperties>,
    external: Vec<Wrench>,
    mobility: Vec<f64>,
    normal_compliance: Vec<f64>,
    shear: PreparedTangentialContacts,
    dt_s: f64,
    controls: RigidHertzMindlinIslandControls,
}

impl<'a> PreparedRigidHertzMindlinIsland<'a> {
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidHertzMindlinIslandContact<'a>],
        dt_s: f64,
        controls: RigidHertzMindlinIslandControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RigidHertzMindlinIslandError> {
        Self::prepare(
            bodies, properties, contacts, None, dt_s, controls, cancelled,
        )
    }

    /// Held world COM forces and body-frame torques, in original global order.
    /// Nonparticipating bodies must have zero external wrenches.
    #[allow(clippy::too_many_arguments)]
    pub fn new_driven(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[RigidHertzMindlinIslandContact<'a>],
        external: &[Wrench],
        dt_s: f64,
        controls: RigidHertzMindlinIslandControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RigidHertzMindlinIslandError> {
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
        contacts: &[RigidHertzMindlinIslandContact<'a>],
        external: Option<&[Wrench]>,
        dt_s: f64,
        controls: RigidHertzMindlinIslandControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RigidHertzMindlinIslandError> {
        poll(cancelled)?;
        controls.native.validate(dt_s)?;
        controls.normal_solver.validate()?;
        if !(1..=128).contains(&controls.maximum_normal_iterations)
            || !(1..=1_048_576).contains(&controls.maximum_normal_evaluations)
            || controls.normal_solver.linear_restart > 128
            || controls.normal_solver.max_linear_cycles > 128
            || bodies.len() != properties.len()
            || contacts.is_empty()
            || contacts.len() > controls.native.multiple.maximum_contacts
            || 9 * contacts.len() * contacts.len() > controls.native.multiple.maximum_matrix_terms
            || external.is_some_and(|w| w.len() != bodies.len())
        {
            return Err(invalid("native shape, matrix or original Newton budgets"));
        }
        let mut ordered = contacts.to_vec();
        ordered.sort_by_key(|contact| contact.contact_id);
        let mut body_indices = Vec::with_capacity(2 * ordered.len());
        for (i, contact) in ordered.iter().enumerate() {
            poll(cancelled)?;
            if contact.contact_id.trim().is_empty()
                || (i > 0 && contact.contact_id == ordered[i - 1].contact_id)
                || contact.body_indices[0] == contact.body_indices[1]
                || contact.body_indices.iter().any(|&b| b >= bodies.len())
            {
                return Err(invalid(
                    "unique original IDs and distinct native bodies required",
                ));
            }
            body_indices.extend(contact.body_indices);
        }
        body_indices.sort_unstable();
        body_indices.dedup();
        if body_indices.len() > controls.native.maximum_bodies {
            return Err(invalid("participating native body budget"));
        }
        if let Some(external) = external {
            for (i, &wrench) in external.iter().enumerate() {
                poll(cancelled)?;
                wrench.validate()?;
                if wrench != Wrench::ZERO && body_indices.binary_search(&i).is_err() {
                    return Err(invalid("external work on a nonparticipating body"));
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
            let [a, b] = [0, 1].map(|i| {
                before[indices[i]].point_kinematics(properties[indices[i]], contact.arms_body[i])
            });
            let points = [a?, b?];
            if norm(points[0].point_world.sub(points[1].point_world))?
                > controls.native.geometry_tolerance_m
            {
                return Err(invalid("original contact points do not coincide"));
            }
            let axes = [
                contact.frame.normal(),
                contact.frame.longitudinal(),
                contact.frame.lateral(),
            ]
            .map(vector);
            let (rates, torsion) = motion(points, axes)?;
            let stored = number(contact.law.stored_energy_j(contact.history))?;
            if stored < 0.0 || contact.history.indentation_m() <= 0.0 {
                return Err(invalid(
                    "original continuously compressed source history required",
                ));
            }
            admit_motion(
                contact.contact_id,
                contact.law,
                contact.history.indentation_m(),
                rates,
                torsion,
                dt_s,
                controls.native,
            )?;
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
            let relative_acceleration = acceleration[0].sub(acceleration[1]);
            let drive = [
                number(relative_acceleration.dot(axes[0]))?,
                number(relative_acceleration.dot(axes[1]))?,
            ];
            let released = [
                number(dt_s * rates[0] + (0.5 * dt_s) * (dt_s * drive[0]))?,
                number(dt_s * rates[1] + (0.5 * dt_s) * (dt_s * drive[1]))?,
            ];
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
        separated(
            &patches,
            &patches
                .iter()
                .map(|p| p.history.contact_radius_m())
                .collect::<Vec<_>>(),
            controls.native.geometry_tolerance_m,
            cancelled,
        )?;
        let count = patches.len();
        let dimension = 3 * count;
        let mut factors = Vec::with_capacity(dimension);
        for patch in &patches {
            for &axis in &patch.axes {
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
        for row in 0..dimension {
            poll(cancelled)?;
            for column in 0..=row {
                let mut value = 0.0;
                for i in 0..2 {
                    for j in 0..2 {
                        if patches[row / 3].bodies[i] == patches[column / 3].bodies[j] {
                            value += factors[row][i].0.dot(factors[column][j].0)
                                + factors[row][i].1.dot(factors[column][j].1);
                        }
                    }
                }
                let value = number(value)?;
                mobility[row * dimension + column] = value;
                mobility[column * dimension + row] = value;
            }
        }
        let mut normal_compliance = vec![0.0; count * count];
        let mut shear_compliance = vec![0.0; 4 * count * count];
        for i in 0..count {
            poll(cancelled)?;
            for j in 0..count {
                for tangent in [1, 2] {
                    let r = 3 * i;
                    let c = 3 * j + tangent;
                    let allowance = 256.0
                        * f64::EPSILON
                        * mobility[r * dimension + r].sqrt()
                        * mobility[c * dimension + c].sqrt();
                    if mobility[r * dimension + c].abs() > allowance {
                        return Err(invalid(
                            "physical normal/shear cross-block requires a coupled solver",
                        ));
                    }
                }
                normal_compliance[i * count + j] =
                    number((0.5 * dt_s) * (dt_s * mobility[(3 * i) * dimension + 3 * j]))?;
                shear_compliance[(2 * i) * (2 * count) + 2 * j] =
                    number((0.5 * dt_s) * (dt_s * mobility[(3 * i + 1) * dimension + 3 * j + 1]))?;
            }
        }
        let shear = PreparedTangentialContacts::new(
            count,
            shear_compliance,
            controls.native.multiple,
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
            normal_compliance,
            shear,
            dt_s,
            controls,
        })
    }

    /// Complete unchanged 3N-by-3N physical normal/longitudinal/lateral mobility.
    pub fn contact_mobility_kg_inverse(&self) -> &[f64] {
        &self.mobility
    }
    pub fn body_indices(&self) -> &[usize] {
        &self.body_indices
    }

    /// Resolve both global load lines from the same captured accepted histories.
    pub fn advance(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RigidHertzMindlinIslandCandidate, RigidHertzMindlinIslandError> {
        poll(cancelled)?;
        let normal = normal::solve(self, cancelled)?;
        let mut responses = Vec::with_capacity(self.patches.len());
        for (patch, &indentation) in self.patches.iter().zip(&normal.indentations) {
            poll(cancelled)?;
            if indentation <= 0.0 {
                return Err(invalid("separation requires the contact-event owner"));
            }
            let rate = (indentation - patch.history.indentation_m()) / self.dt_s;
            responses.push(patch.law.fixed_source().normal_source().evaluate(
                "native-changing-mindlin-island",
                indentation,
                rate,
                self.dt_s,
            )?);
        }
        // Affine normal paths have no intermediate radius larger than both
        // endpoints. Check the whole step's independent-footprint admission.
        let radii: Vec<_> = self
            .patches
            .iter()
            .zip(&responses)
            .map(|(p, r)| p.history.contact_radius_m().max(r.patch_radius_m))
            .collect();
        separated(
            &self.patches,
            &radii,
            self.controls.native.geometry_tolerance_m,
            cancelled,
        )?;
        let laws: Vec<_> = self
            .patches
            .iter()
            .zip(&normal.indentations)
            .zip(&responses)
            .map(|((patch, &indentation), response)| NormalPathShear {
                law: patch.law,
                indentation,
                response,
                dt: self.dt_s,
                failure: RefCell::new(None),
            })
            .collect();
        let inputs: Vec<_> = self
            .patches
            .iter()
            .zip(&laws)
            .map(|(patch, law)| TangentialContactInput {
                law,
                state: &patch.history,
                input: TangentialStepInput {
                    normal_n: law.response.normal_force_n,
                    slip_increment_m: [patch.released[1], 0.0],
                    dt_s: self.dt_s,
                    contact_temperature_k: Some(patch.law.temperature_kelvin()),
                    nominal_pressure_pa: None,
                    compliance_m_per_n: 0.0,
                    solve: self.controls.native.solve,
                },
            })
            .collect();
        let solved = match self.shear.advance(&inputs, &mut *cancelled) {
            Ok(solved) => solved,
            Err(MultipleTangentialError::Cancelled) => {
                return Err(RigidHertzMindlinIslandError::Cancelled);
            }
            Err(error) => {
                for (patch, law) in self.patches.iter().zip(&laws) {
                    if let Some(error) = law.failure.borrow_mut().take() {
                        return Err(RigidHertzMindlinIslandError::SourcePath {
                            contact_id: patch.id.clone(),
                            error,
                        });
                    }
                }
                return Err(error.into());
            }
        };
        let count = self.before.len();
        let mut linear = vec![Vec3::ZERO; count];
        let mut angular = vec![Vec3::ZERO; count];
        let mut impulses = Vec::with_capacity(self.patches.len());
        for ((patch, step), path) in self.patches.iter().zip(&solved.contacts).zip(&normal.paths) {
            poll(cancelled)?;
            if step.traction_n[1] != 0.0 || step.slip_increment_m[1] != 0.0 {
                return Err(invalid("scalar source produced a lateral response"));
            }
            let impulse = patch.axes[0]
                .scale(number(self.dt_s * path.mean_force_n)?)
                .add(patch.axes[1].scale(number(self.dt_s * step.traction_n[0])?));
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
        let gates = self.controls.native;
        for body in 0..count {
            poll(cancelled)?;
            let state = self.before[body];
            let prop = self.properties[body];
            let ext_linear = self.external[body].force_world.scale(self.dt_s);
            let ext_angular = self.external[body].torque_body.scale(self.dt_s);
            let receipt = state.apply_momentum_impulse(
                prop,
                linear[body].add(ext_linear),
                angular[body].add(ext_angular),
            )?;
            audit(
                norm(receipt.linear_residual_world)?,
                gates.linear_momentum_tolerance_ns,
                "native per-body linear impulse",
            )?;
            audit(
                norm(receipt.angular_residual_body)?,
                gates.angular_momentum_tolerance_nms,
                "native per-body angular impulse",
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
        let (mut normal_change, mut shear_change, mut heat, mut contact_work) =
            (0.0, 0.0, 0.0, 0.0);
        let mut normal_residual = normal.residual_m;
        let mut shear_residual = solved.maximum_compatibility_residual_m;
        for ((((patch, step), response), path), impulse) in self
            .patches
            .iter()
            .zip(solved.contacts)
            .zip(responses)
            .zip(&normal.paths)
            .zip(impulses)
        {
            poll(cancelled)?;
            let [a, b] = [0, 1].map(|i| {
                after[patch.bodies[i]]
                    .point_kinematics(self.properties[patch.bodies[i]], patch.arms[i])
            });
            let points = [a?, b?];
            let (rates, torsion) = motion(points, patch.axes)?;
            let indentation = step.state.indentation_m();
            admit_motion(
                &patch.id,
                patch.law,
                indentation,
                rates,
                torsion,
                self.dt_s,
                gates,
            )?;
            let travel = [0, 1, 2].map(|i| self.dt_s * (0.5 * patch.rates[i] + 0.5 * rates[i]));
            let residual = [
                number(travel[0] - (indentation - patch.history.indentation_m()))?,
                number(travel[1] - step.slip_increment_m[0])?,
            ];
            for value in residual {
                audit(
                    value,
                    gates.multiple.local.displacement_tolerance_m,
                    "actual native normal/shear displacement",
                )?;
            }
            audit(
                number(travel[2])?,
                gates.multiple.local.displacement_tolerance_m,
                "unsupported native lateral displacement",
            )?;
            normal_residual = normal_residual.max(residual[0].abs());
            shear_residual = shear_residual.max(residual[1].abs());
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
            let normal_stored =
                number(response.reversible_energy_j - patch.history.normal_energy_j())?;
            let shear_stored = number(patch.law.stored_energy_j(&step.state) - patch.stored)?;
            let balance = number(normal_stored + shear_stored + step.dissipated_j + work)?;
            for (value, what) in [
                (
                    path.mean_force_n * travel[0] - normal_stored,
                    "actual original normal work",
                ),
                (
                    path.work_j - normal_stored,
                    "original normal path potential",
                ),
                (
                    shear_stored - step.stored_change_j,
                    "original retained shear storage",
                ),
                (
                    shear_stored + step.dissipated_j - step.traction_n[0] * travel[1],
                    "actual original shear work",
                ),
                (balance, "native contact work/storage/heat"),
            ] {
                audit(
                    number(value)?,
                    gates.multiple.local.energy_tolerance_j,
                    what,
                )?;
            }
            normal_change = number(normal_change + normal_stored)?;
            shear_change = number(shear_change + shear_stored)?;
            heat = number(heat + step.dissipated_j)?;
            contact_work = number(contact_work + work)?;
            contacts.push(RigidHertzMindlinIslandContactCandidate {
                contact_id: patch.id.clone(),
                contact: SourceMindlinPathStep {
                    tangential: step,
                    normal_response: response,
                    normal_work_j: normal_stored,
                    normal_stored_change_j: normal_stored,
                },
                normal_average_force_n: path.mean_force_n,
                relative_velocity_before_m_s: patch.rates,
                relative_velocity_after_m_s: rates,
                relative_torsion_before_rad_s: patch.torsion,
                relative_torsion_after_rad_s: torsion,
                driving_acceleration_m_s2: patch.drive,
                relative_displacement_m: travel,
                impulse_on_b_world_ns: impulse,
                midpoint_impulse_work_j: work,
                contact_energy_residual_j: balance,
                compatibility_residual_m: residual,
            });
        }
        let linear_residual = delta_linear.sub(external_linear);
        let angular_residual = delta_angular.sub(external_angular);
        audit(
            norm(linear_residual)?,
            gates.linear_momentum_tolerance_ns,
            "island linear momentum",
        )?;
        audit(
            norm(angular_residual)?,
            gates.angular_momentum_tolerance_nms,
            "island angular momentum",
        )?;
        let work_residual = number(kinetic_change - contact_work - external_work)?;
        let energy_residual =
            number(kinetic_change + normal_change + shear_change + heat - external_work)?;
        audit(
            work_residual,
            gates.energy_tolerance_j,
            "complete native contact/drive work",
        )?;
        audit(
            energy_residual,
            gates.energy_tolerance_j,
            "complete normal/shear/heat/drive energy",
        )?;
        poll(cancelled)?;
        Ok(RigidHertzMindlinIslandCandidate {
            body_indices: self.body_indices.clone(),
            bodies_after: after,
            contacts,
            report: RigidHertzMindlinIslandReport {
                normal_iterations: normal.iterations,
                normal_source_evaluations: normal.evaluations,
                maximum_normal_residual_m: normal_residual,
                shear_sweeps: solved.sweeps,
                shear_law_evaluations: solved.law_evaluations,
                maximum_shear_residual_m: shear_residual,
                kinetic_before_j: kinetic_before,
                kinetic_after_j: kinetic_after,
                kinetic_change_j: kinetic_change,
                normal_stored_change_j: normal_change,
                shear_stored_change_j: shear_change,
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

// A true tangential constitutive adapter, with the normal path fixed by the
// separate native normal owner. Original source admission and shear storage/
// dissipation remain unchanged; normal displacement is never a tangent axis.
struct NormalPathShear<'a> {
    law: &'a SourceHertzMindlinContact,
    indentation: f64,
    response: &'a PointNormalPatchReceipt,
    dt: f64,
    failure: RefCell<Option<SourceMindlinPathError>>,
}
impl TangentialLaw for NormalPathShear<'_> {
    type State = SourceHertzMindlinState;
    fn model_id(&self) -> &str {
        "source-hertz-mindlin-native-normal-path-v1"
    }
    fn initial_state(&self) -> Self::State {
        self.law.initial_state()
    }
    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        self.law.stored_energy_j(state)
    }
    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        *self.failure.borrow_mut() = None;
        if input.normal_n != self.response.normal_force_n
            || input.dt_s != self.dt
            || input.slip_increment_m[1] != 0.0
            || input.nominal_pressure_pa.is_some()
            || input.contact_temperature_k != Some(self.law.temperature_kelvin())
        {
            return Err(TriboError::InvalidInput {
                field: "original normal path and scalar shear conditions",
            });
        }
        self.law
            .advance_coupled_path(
                state,
                &SourceMindlinPathInput {
                    indentation_m: self.indentation,
                    free_slip_increment_m: input.slip_increment_m[0],
                    dt_s: self.dt,
                    compliance_m_per_n: input.compliance_m_per_n,
                    solve: input.solve,
                },
                &mut || false,
            )
            .map(|path| path.tangential)
            .map_err(|error| {
                *self.failure.borrow_mut() = Some(error);
                TriboError::InvalidInput {
                    field: "original simultaneous source Mindlin path refused",
                }
            })
    }
}

fn separated(
    patches: &[Patch<'_>],
    radii: &[f64],
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), RigidHertzMindlinIslandError> {
    for (i, patch) in patches.iter().enumerate() {
        poll(cancelled)?;
        for (j, other) in patches[..i].iter().enumerate() {
            let clearance = number(radii[i] + radii[j] + tolerance)?;
            for a in 0..2 {
                for b in 0..2 {
                    if patch.bodies[a] == other.bodies[b]
                        && norm(patch.arms[a].sub(other.arms[b]))? <= clearance
                    {
                        return Err(invalid(
                            "independent changing Hertz footprints are not separated",
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}
fn motion(
    points: [PointKinematics; 2],
    axes: [Vec3; 3],
) -> Result<([f64; 3], f64), RigidHertzMindlinIslandError> {
    let relative = points[0]
        .point_velocity_world
        .sub(points[1].point_velocity_world);
    let rates = [
        number(relative.dot(axes[0]))?,
        number(relative.dot(axes[1]))?,
        number(relative.dot(axes[2]))?,
    ];
    let torsion = number(
        points[0]
            .angular_velocity_world
            .sub(points[1].angular_velocity_world)
            .dot(axes[0]),
    )?;
    Ok((rates, torsion))
}
fn admit_motion(
    id: &str,
    law: &SourceHertzMindlinContact,
    indentation: f64,
    rates: [f64; 3],
    torsion: f64,
    dt: f64,
    gates: RigidMindlinIslandControls,
) -> Result<(), RigidHertzMindlinIslandError> {
    law.fixed_source()
        .normal_source()
        .validate_kinematics(indentation, rates[0], dt)?;
    for (value, allowance, coordinate) in [
        (rates[2], gates.velocity_tolerance_m_s, "lateral velocity"),
        (
            torsion,
            gates.angular_velocity_tolerance_rad_s,
            "relative torsion",
        ),
    ] {
        if value.abs() > allowance {
            return Err(RigidMindlinIslandError::UnsupportedMotion {
                contact_id: id.into(),
                coordinate,
            }
            .into());
        }
    }
    if rates[1].abs() > gates.maximum_slip_speed_m_s {
        return Err(invalid("actual native shear endpoint exceeds source phase"));
    }
    Ok(())
}
fn invalid(message: &'static str) -> RigidHertzMindlinIslandError {
    super::invalid(message).into()
}
