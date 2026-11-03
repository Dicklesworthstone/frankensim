//! Native free rigid bodies striking a continuing modal/material network.
//!
//! One original simultaneous Coulomb solve uses the complete contact mobility
//! `B_modal B_modal^T + B_rigid M_rigid^-1 B_rigid^T`, including shared modes,
//! shared bodies, off-centre arms and the original full body inertia. The SAME
//! rows transpose its impulses into native body and modal momenta. No proxy
//! striker mass, sequential normal/friction clip, or restarted modal history
//! replaces either receiving mechanical owner.
//!
//! These are explicitly selected material-point witnesses and frozen world
//! frames at one event. Modal point geometry is the supplied linear trace in
//! the original retained basis. Spatial momentum conservation extends only to
//! that basis's represented symmetries; a supported mode is not a free body.
//! There is no collision search, frame transport, sustained contact integrator
//! or automatic heat allocation. Redundant loads retain the native solver's
//! determinacy and work-witness qualifications.

use super::{
    CoupledModalSystem, ModalAcousticTimeModel, ModalCouplingError, audit_modal_momentum, checked,
    finite, invalid, poll,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResiduals,
    FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3};

/// One physical material point in a component's complete original modal basis.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalImpactPoint {
    pub component: usize,
    /// World position at zero modal displacement [m].
    pub reference_position_world_m: Vec3,
    /// World displacement per mass-normalized coordinate [1/sqrt(kg)].
    /// Every original mode is present, including explicit zero traces.
    pub shapes_world: Vec<Vec3>,
}

/// One rigid witness receiving +p and one modal witness receiving -p.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidModalImpactContact {
    pub id: u64,
    pub rigid_body: usize,
    /// Original COM-relative material arm in the rigid body's frame [m].
    pub arm_body_m: Vec3,
    pub modal_point: ModalImpactPoint,
    /// Points from modal to rigid witness; positive relative speed separates.
    pub normal_world: Vec3,
    /// Its normal component is removed before forming the orthonormal frame.
    pub tangent_world: Vec3,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidModalImpactConfig {
    pub maximum_bodies: usize,
    /// Original witnesses must share a normal line within this distance [m].
    pub geometry_tolerance_m: f64,
    /// Each actual witness separation crossed with its impulse must satisfy
    /// this bound [N m s]. Finite-gap friction never silently relocates a point.
    pub angular_impulse_tolerance_n_m_s: f64,
    pub contact: FrictionalContactConfig,
}
impl Default for RigidModalImpactConfig {
    fn default() -> Self {
        Self {
            maximum_bodies: 64,
            geometry_tolerance_m: 1e-10,
            angular_impulse_tolerance_n_m_s: 1e-10,
            contact: FrictionalContactConfig::default(),
        }
    }
}

/// Cached original modal rows/Gram and bounded native projection scratch.
/// Rigid geometry, row moments and full inertia response are recomputed at
/// EVERY event from the supplied current bodies/properties. Reprepare if the
/// modal trace, selected material witnesses or physical contact frame changes.
/// The original contact solver and returned receipt allocate their own data.
pub struct PreparedRigidModalImpacts {
    contacts: Vec<RigidModalImpactContact>,
    axes: Vec<[Vec3; 3]>,
    columns: Vec<Vec<f64>>,
    modal_gram: Vec<Dd>,
    delassus: Vec<f64>,
    layout: Vec<usize>,
    config: RigidModalImpactConfig,
    sources: Vec<FrictionalContact>,
    separation: Vec<[Dd; 3]>,
    before: Vec<[Dd; 3]>,
    moments: Vec<Vec3>,
    angular_response: Vec<Vec3>,
    requested_modal: Vec<Dd>,
    modal_impulses: Vec<f64>,
    requested_rigid: Vec<[Dd; 6]>,
    candidates: Vec<RigidBodyState>,
}

/// The native constitutive result plus independent actual-state admission.
/// All contact arrays use canonical contact-ID order. Rigid receipts retain
/// the caller's body order, including zero impulses on uninvolved bodies.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidModalImpactReceipt {
    pub sample: u64,
    pub contact_ids: Vec<u64>,
    pub axes_world: Vec<[Vec3; 3]>,
    pub gaps_m: Vec<f64>,
    pub contact_mobility_kg_inv: Vec<f64>,
    pub contacts: FrictionalContactResult,
    pub rigid_impulses: Vec<MomentumImpulseReceipt>,
    pub relative_velocities_before_m_s: Vec<[f64; 3]>,
    pub relative_velocities_after_m_s: Vec<[f64; 3]>,
    /// Polarized actual increments, never subtraction of large total energies.
    pub modal_kinetic_energy_change_j: f64,
    pub rigid_kinetic_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub contact_work_j: f64,
    pub energy_residual_j: f64,
    pub contact_work_residual_j: f64,
    pub maximum_relative_momentum_residual: f64,
    pub maximum_velocity_residual_m_s: f64,
    pub maximum_geometric_angular_impulse_defect_n_m_s: f64,
    pub actual_residuals: FrictionalContactResiduals,
}

#[derive(Debug)]
pub enum RigidModalImpactError {
    Coupling(ModalCouplingError),
    Rigid(DynamicsError),
    Contact(FrictionalContactError),
    Momentum {
        body: usize,
        /// 0..3 are world linear; 3..6 are central body angular momentum.
        coordinate: usize,
        relative: f64,
        tolerance: f64,
    },
    Velocity {
        contact: usize,
        axis: usize,
        residual_m_s: f64,
        tolerance_m_s: f64,
    },
}
impl From<ModalCouplingError> for RigidModalImpactError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Coupling(error)
    }
}
impl From<DynamicsError> for RigidModalImpactError {
    fn from(error: DynamicsError) -> Self {
        Self::Rigid(error)
    }
}
impl From<FrictionalContactError> for RigidModalImpactError {
    fn from(error: FrictionalContactError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for RigidModalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid/modal impact refused: {self:?}")
    }
}
impl std::error::Error for RigidModalImpactError {}

impl PreparedRigidModalImpacts {
    pub fn new(
        system: &CoupledModalSystem,
        mut contacts: Vec<RigidModalImpactContact>,
        config: RigidModalImpactConfig,
        gate: &CancelGate,
    ) -> Result<Self, RigidModalImpactError> {
        poll(Some(gate))?;
        system.require_time_domain()?;
        let count = contacts.len();
        if count > 21
            || count > config.contact.maximum_contacts
            || config.maximum_bodies == 0
            || !config.geometry_tolerance_m.is_finite()
            || config.geometry_tolerance_m <= 0.0
            || !config.angular_impulse_tolerance_n_m_s.is_finite()
            || config.angular_impulse_tolerance_n_m_s <= 0.0
        {
            return Err(invalid(
                "bounded rigid/modal contacts and positive finite controls required",
            )
            .into());
        }
        let dimension = 3 * count;
        let work = dimension
            .checked_mul(dimension)
            .and_then(|entries| entries.checked_mul(system.mode_count().checked_add(6)?))
            .and_then(|entries| entries.checked_add(config.maximum_bodies.checked_mul(6)?))
            .ok_or_else(|| invalid("rigid/modal preparation extent overflow"))?;
        if work > system.config.max_setup_terms {
            return Err(invalid("rigid/modal setup and scratch exceed max_setup_terms").into());
        }
        let _ = solve_frictional_contacts(&[], &[], config.contact, || gate.is_requested())?;
        contacts.sort_by_key(|contact| contact.id);
        let mut axes = Vec::with_capacity(count);
        let mut columns = Vec::with_capacity(dimension);
        let mut sources = Vec::with_capacity(count);
        for (index, contact) in contacts.iter().enumerate() {
            poll(Some(gate))?;
            if (index > 0 && contacts[index - 1].id == contact.id)
                || contact.rigid_body >= config.maximum_bodies
                || !contact.arm_body_m.is_finite()
                || !contact.modal_point.reference_position_world_m.is_finite()
                || !contact.activation_tolerance_m.is_finite()
                || contact.activation_tolerance_m < 0.0
                || !contact.friction_coefficient.is_finite()
                || contact.friction_coefficient < 0.0
                || !match contact.normal_law {
                    NormalVelocityLaw::NewtonRestitution { coefficient } => {
                        coefficient.is_finite() && (0.0..=1.0).contains(&coefficient)
                    }
                    NormalVelocityLaw::Target { speed } => speed.is_finite() && speed >= 0.0,
                }
            {
                return Err(
                    invalid("rigid/modal identity, witness, friction or normal law").into(),
                );
            }
            let model = system
                .models
                .get(contact.modal_point.component)
                .ok_or_else(|| invalid("modal material point is outside the original network"))?;
            if contact.modal_point.shapes_world.len() != model.modes().len()
                || contact
                    .modal_point
                    .shapes_world
                    .iter()
                    .any(|shape| !shape.is_finite())
            {
                return Err(invalid(
                    "world point traces must cover the complete original modal basis",
                )
                .into());
            }
            let normal = unit(contact.normal_world)?;
            let tangent = unit(
                contact
                    .tangent_world
                    .sub(normal.scale(checked(dot(normal, contact.tangent_world))?)),
            )?;
            let frame = [normal, tangent, normal.cross(tangent)];
            for axis in frame {
                let mut column = vec![0.0; system.mode_count()];
                for (mode, &shape) in contact.modal_point.shapes_world.iter().enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    column[system.offsets[contact.modal_point.component] + mode] =
                        -checked(dot(axis, shape))?;
                }
                columns.push(column);
            }
            axes.push(frame);
            sources.push(FrictionalContact {
                id: contact.id,
                gap: 0.0,
                activation_tolerance: contact.activation_tolerance_m,
                relative_velocity: [0.0; 3],
                friction_coefficient: contact.friction_coefficient,
                normal_law: contact.normal_law,
            });
        }
        let mut modal_gram = vec![Dd::ZERO; dimension * dimension];
        for row in 0..dimension {
            poll(Some(gate))?;
            for col in row..dimension {
                let mut value = Dd::ZERO;
                for (mode, (&a, &b)) in columns[row].iter().zip(&columns[col]).enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    value = value + Dd::from_f64(a) * Dd::from_f64(b);
                }
                checked(value)?;
                modal_gram[row * dimension + col] = value;
                modal_gram[col * dimension + row] = value;
            }
        }
        poll(Some(gate))?;
        Ok(Self {
            contacts,
            axes,
            columns,
            delassus: vec![0.0; dimension * dimension],
            modal_gram,
            layout: system
                .models
                .iter()
                .map(|model| model.modes().len())
                .collect(),
            config,
            sources,
            separation: vec![[Dd::ZERO; 3]; count],
            before: vec![[Dd::ZERO; 3]; count],
            moments: vec![Vec3::ZERO; dimension],
            angular_response: vec![Vec3::ZERO; dimension],
            requested_modal: vec![Dd::ZERO; system.mode_count()],
            modal_impulses: vec![0.0; system.mode_count()],
            requested_rigid: Vec::with_capacity(config.maximum_bodies),
            candidates: Vec::with_capacity(config.maximum_bodies),
        })
    }

    pub fn contact_count(&self) -> usize {
        self.contacts.len()
    }

    /// Original physical modal incidence, in canonical contact-ID order.
    pub fn modal_components(&self) -> impl Iterator<Item = usize> + '_ {
        self.contacts
            .iter()
            .map(|contact| contact.modal_point.component)
    }

    fn modal_trace(
        &self,
        index: usize,
        models: &[ModalAcousticTimeModel],
        velocity: bool,
        gate: &CancelGate,
    ) -> Result<[Dd; 3], ModalCouplingError> {
        let mut result = [Dd::ZERO; 3];
        for (mode, state) in models.iter().flat_map(|model| model.states()).enumerate() {
            if mode % 64 == 0 {
                poll(Some(gate))?;
            }
            let value = Dd::from_f64(if velocity {
                state.velocity_m_sqrt_kg_per_s
            } else {
                state.displacement_m_sqrt_kg
            });
            for (axis, sum) in result.iter_mut().enumerate() {
                *sum = *sum + Dd::from_f64(self.columns[3 * index + axis][mode]) * value;
            }
        }
        Ok(result)
    }

    fn velocity(
        &self,
        index: usize,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        models: &[ModalAcousticTimeModel],
        gate: &CancelGate,
    ) -> Result<[Dd; 3], RigidModalImpactError> {
        let mut velocity = self.modal_trace(index, models, true, gate)?;
        let body = self.contacts[index].rigid_body;
        let v = bodies[body].center_of_mass_velocity_world(properties[body])?;
        let omega =
            properties[body].angular_velocity_body_checked(bodies[body].angular_momentum_body())?;
        for (axis, value) in velocity.iter_mut().enumerate() {
            *value = *value
                + dot(self.axes[index][axis], v)
                + dot(self.moments[3 * index + axis], omega);
            checked(*value)?;
        }
        Ok(velocity)
    }
}

impl CoupledModalSystem {
    /// Resolve simultaneous Coulomb impacts between current native rigid bodies
    /// and this SAME continuing network. Publish both receivers only after the
    /// original constitutive/determinacy solve, per-coordinate actual momentum,
    /// actual Coulomb graph, geometric couple, work and native budgets pass.
    ///
    /// All rigid poses, modal q, material/radiation histories, clock and preceding
    /// interval receipts remain unchanged. Every refusal/cancellation publishes
    /// nothing. The complete receiving state may be reused after a failed call.
    pub fn apply_rigid_frictional_impacts(
        &mut self,
        bodies: &mut [RigidBodyState],
        properties: &[MassProperties],
        impacts: &mut PreparedRigidModalImpacts,
        gate: &CancelGate,
    ) -> Result<RigidModalImpactReceipt, RigidModalImpactError> {
        self.apply_rigid_frictional_impacts_with_admission(
            bodies,
            properties,
            impacts,
            gate,
            |_| Ok(()),
        )
    }

    pub(crate) fn apply_rigid_frictional_impacts_with_admission(
        &mut self,
        bodies: &mut [RigidBodyState],
        properties: &[MassProperties],
        impacts: &mut PreparedRigidModalImpacts,
        gate: &CancelGate,
        admit: impl FnOnce(&[ModalAcousticTimeModel]) -> Result<(), ModalCouplingError>,
    ) -> Result<RigidModalImpactReceipt, RigidModalImpactError> {
        poll(Some(gate))?;
        self.require_time_domain()?;
        if bodies.len() != properties.len()
            || bodies.len() > impacts.config.maximum_bodies
            || self.models.len() != impacts.layout.len()
            || self
                .models
                .iter()
                .zip(&impacts.layout)
                .any(|(model, &count)| model.modes().len() != count)
            || impacts
                .contacts
                .iter()
                .any(|contact| contact.rigid_body >= bodies.len())
        {
            return Err(
                invalid("prepared rigid/modal contact no longer matches native receivers").into(),
            );
        }
        let dimension = 3 * impacts.contacts.len();
        for index in 0..impacts.contacts.len() {
            poll(Some(gate))?;
            let contact = &impacts.contacts[index];
            let body = bodies[contact.rigid_body];
            let property = properties[contact.rigid_body];
            let orientation = body.pose().orientation();
            let arm_world = orientation.rotate_body_to_world(contact.arm_body_m);
            let centre = components(body.pose().position_world());
            let reference = components(contact.modal_point.reference_position_world_m);
            let arm = components(arm_world);
            let mut separation = impacts.modal_trace(index, &self.models, false, gate)?;
            for axis in 0..3 {
                // Subtract world origins BEFORE adding the local arm/deflection.
                // A distant origin must not erase either material witness.
                for coordinate in 0..3 {
                    separation[axis] = separation[axis]
                        + Dd::from_f64(components(impacts.axes[index][axis])[coordinate])
                            * (Dd::from_f64(centre[coordinate])
                                - Dd::from_f64(reference[coordinate])
                                + Dd::from_f64(arm[coordinate]));
                }
                let moment = contact
                    .arm_body_m
                    .cross(orientation.rotate_world_to_body(impacts.axes[index][axis]));
                if !moment.is_finite() {
                    return Err(invalid("rigid contact row moment is not representable").into());
                }
                impacts.moments[3 * index + axis] = moment;
                impacts.angular_response[3 * index + axis] =
                    property.angular_velocity_body_checked(moment)?;
            }
            if checked(separation[1])?.hypot(checked(separation[2])?)
                > impacts.config.geometry_tolerance_m
            {
                return Err(
                    invalid("original rigid/modal witnesses do not share a normal line").into(),
                );
            }
            impacts.separation[index] = separation;
            impacts.sources[index].gap = checked(separation[0])?;
            let before = impacts.velocity(index, bodies, properties, &self.models, gate)?;
            impacts.before[index] = before;
            impacts.sources[index].relative_velocity = [
                checked(before[0])?,
                checked(before[1])?,
                checked(before[2])?,
            ];
        }
        for (value, &modal) in impacts.delassus.iter_mut().zip(&impacts.modal_gram) {
            *value = checked(modal)?;
        }
        for row in 0..dimension {
            poll(Some(gate))?;
            for col in row..dimension {
                let body = impacts.contacts[row / 3].rigid_body;
                if body == impacts.contacts[col / 3].rigid_body {
                    let rigid = dot(
                        impacts.axes[row / 3][row % 3],
                        impacts.axes[col / 3][col % 3],
                    ) / Dd::from_f64(properties[body].mass())
                        + dot(impacts.moments[row], impacts.angular_response[col]);
                    let value = checked(impacts.modal_gram[row * dimension + col] + rigid)?;
                    impacts.delassus[row * dimension + col] = value;
                    impacts.delassus[col * dimension + row] = value;
                }
            }
        }
        let contacts = solve_frictional_contacts(
            &impacts.sources,
            &impacts.delassus,
            impacts.config.contact,
            || gate.is_requested(),
        )?;
        impacts.requested_modal.fill(Dd::ZERO);
        impacts.requested_rigid.clear();
        impacts.requested_rigid.resize(bodies.len(), [Dd::ZERO; 6]);
        let mut maximum_geometric_angular_impulse_defect_n_m_s = 0.0_f64;
        for index in 0..impacts.contacts.len() {
            poll(Some(gate))?;
            let p = contacts.impulses[index];
            let separation = impacts.separation[index];
            let couple = [
                separation[1] * Dd::from_f64(p[2]) - separation[2] * Dd::from_f64(p[1]),
                separation[2] * Dd::from_f64(p[0]) - separation[0] * Dd::from_f64(p[2]),
                separation[0] * Dd::from_f64(p[1]) - separation[1] * Dd::from_f64(p[0]),
            ];
            let defect = checked(couple[0])?
                .hypot(checked(couple[1])?)
                .hypot(checked(couple[2])?);
            if defect > impacts.config.angular_impulse_tolerance_n_m_s {
                return Err(invalid(
                    "separated rigid/modal witnesses create an unresolved impulse couple",
                )
                .into());
            }
            maximum_geometric_angular_impulse_defect_n_m_s =
                maximum_geometric_angular_impulse_defect_n_m_s.max(defect);
            for axis in 0..3 {
                let linear = components(impacts.axes[index][axis]);
                let angular = components(impacts.moments[3 * index + axis]);
                let requested = &mut impacts.requested_rigid[impacts.contacts[index].rigid_body];
                for coordinate in 0..3 {
                    requested[coordinate] =
                        requested[coordinate] + product(linear[coordinate], p[axis])?;
                    requested[3 + coordinate] =
                        requested[3 + coordinate] + product(angular[coordinate], p[axis])?;
                }
                for (mode, requested) in impacts.requested_modal.iter_mut().enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    *requested =
                        *requested + product(impacts.columns[3 * index + axis][mode], p[axis])?;
                }
            }
        }
        for (value, &requested) in impacts
            .modal_impulses
            .iter_mut()
            .zip(&impacts.requested_modal)
        {
            *value = checked(requested)?;
        }
        let changes_modal = impacts.modal_impulses.iter().any(|&impulse| impulse != 0.0);
        if changes_modal {
            self.stage_modal_impulse(&impacts.modal_impulses, Some(gate))?;
        }
        let after_models = if changes_modal {
            &self.candidates
        } else {
            &self.models
        };
        let mut maximum_relative_momentum_residual =
            audit_modal_momentum(self, after_models, &impacts.requested_modal, gate)?;
        let mut rigid_impulses = Vec::with_capacity(bodies.len());
        let mut rigid_change = Dd::ZERO;
        impacts.candidates.clear();
        for (index, (&body, &property)) in bodies.iter().zip(properties).enumerate() {
            poll(Some(gate))?;
            let requested = impacts.requested_rigid[index];
            let receipt = body.apply_momentum_impulse(
                property,
                Vec3::new(
                    checked(requested[0])?,
                    checked(requested[1])?,
                    checked(requested[2])?,
                ),
                Vec3::new(
                    checked(requested[3])?,
                    checked(requested[4])?,
                    checked(requested[5])?,
                ),
            )?;
            let old = momentum(body);
            let new = momentum(receipt.state_after);
            let v0 = body.center_of_mass_velocity_world(property)?;
            let v1 = receipt
                .state_after
                .center_of_mass_velocity_world(property)?;
            let w0 = property.angular_velocity_body_checked(body.angular_momentum_body())?;
            let w1 = property
                .angular_velocity_body_checked(receipt.state_after.angular_momentum_body())?;
            let old_flow = six(v0, w0);
            let new_flow = six(v1, w1);
            for coordinate in 0..6 {
                let actual = Dd::from_f64(new[coordinate]) - Dd::from_f64(old[coordinate]);
                let scale = checked(actual)?
                    .abs()
                    .max(checked(requested[coordinate])?.abs());
                let residual = checked(actual - requested[coordinate])?.abs();
                let relative = if scale == 0.0 { 0.0 } else { residual / scale };
                if relative > self.config.solve_relative_tolerance {
                    return Err(RigidModalImpactError::Momentum {
                        body: index,
                        coordinate,
                        relative,
                        tolerance: self.config.solve_relative_tolerance,
                    });
                }
                maximum_relative_momentum_residual =
                    maximum_relative_momentum_residual.max(relative);
                rigid_change = rigid_change
                    + actual
                        * (Dd::from_f64(0.5 * old_flow[coordinate])
                            + Dd::from_f64(0.5 * new_flow[coordinate]));
            }
            impacts.candidates.push(receipt.state_after);
            rigid_impulses.push(receipt);
        }
        let mut modal_change = Dd::ZERO;
        for (mode, (old, new)) in self
            .models
            .iter()
            .flat_map(|model| model.states())
            .zip(after_models.iter().flat_map(|model| model.states()))
            .enumerate()
        {
            if mode % 64 == 0 {
                poll(Some(gate))?;
            }
            let old = old.velocity_m_sqrt_kg_per_s;
            let new = new.velocity_m_sqrt_kg_per_s;
            modal_change = modal_change
                + (Dd::from_f64(new) - Dd::from_f64(old))
                    * (Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new));
        }
        let mut after_velocities = Vec::with_capacity(impacts.contacts.len());
        let mut actual_residuals = FrictionalContactResiduals::default();
        let mut maximum_velocity_residual_m_s = 0.0_f64;
        let mut work = Dd::ZERO;
        for index in 0..impacts.contacts.len() {
            let after =
                impacts.velocity(index, &impacts.candidates, properties, after_models, gate)?;
            let mut rounded = [0.0; 3];
            for axis in 0..3 {
                rounded[axis] = checked(after[axis])?;
                let predicted = contacts.relative_velocities_after[index][axis];
                let residual = checked(after[axis] - Dd::from_f64(predicted))?;
                let tolerance = finite(
                    impacts.config.contact.velocity_tolerance
                        + self.config.solve_relative_tolerance
                            * rounded[axis]
                                .abs()
                                .max(predicted.abs())
                                .max(impacts.sources[index].relative_velocity[axis].abs()),
                )?;
                if residual.abs() > tolerance {
                    return Err(RigidModalImpactError::Velocity {
                        contact: index,
                        axis,
                        residual_m_s: residual,
                        tolerance_m_s: tolerance,
                    });
                }
                maximum_velocity_residual_m_s = maximum_velocity_residual_m_s.max(residual.abs());
                work = work
                    + Dd::from_f64(contacts.impulses[index][axis])
                        * (impacts.before[index][axis] * Dd::from_f64(0.5)
                            + after[axis] * Dd::from_f64(0.5));
            }
            after_velocities.push(rounded);
            if impacts.sources[index].gap <= impacts.sources[index].activation_tolerance {
                super::frictional::audit_contact(
                    impacts.sources[index].friction_coefficient,
                    contacts.normal_targets[index],
                    contacts.impulses[index],
                    after,
                    impacts.config.contact,
                    &mut actual_residuals,
                )?;
            }
        }
        let r = actual_residuals;
        let controls = impacts.config.contact;
        if r.normal_impulse_violation > controls.impulse_tolerance
            || r.normal_velocity_violation > controls.velocity_tolerance
            || r.normal_complementarity > controls.complementarity_tolerance
            || r.cone_violation > controls.impulse_tolerance
            || r.sliding_impulse_residual > controls.impulse_tolerance
            || r.tangential_duality_gap > controls.complementarity_tolerance
        {
            return Err(FrictionalContactError::ResidualNotSatisfied(r).into());
        }
        let contact_work_j = checked(work)?;
        if contact_work_j > controls.energy_tolerance {
            return Err(FrictionalContactError::EnergyIncrease {
                work: contact_work_j,
                allowed: controls.energy_tolerance,
            }
            .into());
        }
        let change = modal_change + rigid_change;
        let kinetic_energy_change_j = checked(change)?;
        let energy_residual_j = checked(change - work)?;
        let contact_work_residual_j = checked(work - Dd::from_f64(contacts.relative_work))?;
        for (residual_j, scale) in [
            (
                energy_residual_j,
                kinetic_energy_change_j.abs().max(contact_work_j.abs()),
            ),
            (
                contact_work_residual_j,
                contact_work_j.abs().max(contacts.relative_work.abs()),
            ),
        ] {
            let tolerance_j = finite(
                self.config.energy_absolute_tolerance_j
                    + self.config.energy_relative_tolerance * scale,
            )?;
            if residual_j.abs() > tolerance_j {
                return Err(ModalCouplingError::EnergyBalance {
                    residual_j,
                    tolerance_j,
                }
                .into());
            }
        }
        let receipt = RigidModalImpactReceipt {
            sample: self.samples_rendered(),
            contact_ids: impacts.sources.iter().map(|source| source.id).collect(),
            axes_world: impacts.axes.clone(),
            gaps_m: impacts.sources.iter().map(|source| source.gap).collect(),
            contact_mobility_kg_inv: impacts.delassus.clone(),
            contacts,
            rigid_impulses,
            relative_velocities_before_m_s: impacts
                .sources
                .iter()
                .map(|source| source.relative_velocity)
                .collect(),
            relative_velocities_after_m_s: after_velocities,
            modal_kinetic_energy_change_j: checked(modal_change)?,
            rigid_kinetic_energy_change_j: checked(rigid_change)?,
            kinetic_energy_change_j,
            contact_work_j,
            energy_residual_j,
            contact_work_residual_j,
            maximum_relative_momentum_residual,
            maximum_velocity_residual_m_s,
            maximum_geometric_angular_impulse_defect_n_m_s,
            actual_residuals,
        };
        admit(after_models)?;
        poll(Some(gate))?;
        if changes_modal {
            self.publish_modal_impulse();
        }
        bodies.copy_from_slice(&impacts.candidates);
        Ok(receipt)
    }
}

fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn six(linear: Vec3, angular: Vec3) -> [f64; 6] {
    [
        linear.x, linear.y, linear.z, angular.x, angular.y, angular.z,
    ]
}
fn momentum(body: RigidBodyState) -> [f64; 6] {
    six(body.linear_momentum_world(), body.angular_momentum_body())
}
fn dot(a: Vec3, b: Vec3) -> Dd {
    components(a)
        .into_iter()
        .zip(components(b))
        .fold(Dd::ZERO, |sum, (a, b)| {
            sum + Dd::from_f64(a) * Dd::from_f64(b)
        })
}
fn product(a: f64, b: f64) -> Result<Dd, ModalCouplingError> {
    let value = Dd::from_f64(a) * Dd::from_f64(b);
    checked(value)?;
    if a != 0.0 && b != 0.0 && value.hi == 0.0 && value.lo == 0.0 {
        return Err(invalid(
            "rigid/modal projection loses a nonzero momentum term",
        ));
    }
    Ok(value)
}
fn unit(value: Vec3) -> Result<Vec3, ModalCouplingError> {
    if !value.is_finite() {
        return Err(invalid("finite contact frame required"));
    }
    let scale = value.x.abs().max(value.y.abs()).max(value.z.abs());
    if scale == 0.0 {
        return Err(invalid(
            "nonzero independent contact frame directions required",
        ));
    }
    let value = Vec3::new(value.x / scale, value.y / scale, value.z / scale);
    Ok(value.scale(value.x.hypot(value.y).hypot(value.z).recip()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
    use crate::render::schedule::force::coupled::ModalCouplingConfig;
    use fs_math::c64::C64;
    use fs_mbd::{Pose, UnitQuaternion};

    #[test]
    fn late_admission_and_final_cancellation_publish_neither_receiver() {
        let gate = CancelGate::new_clock_free();
        let model = ModalAcousticTimeModel::try_new_with_free_coordinates(
            2000,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        let config = ModalCouplingConfig {
            max_modes: 4,
            max_connections: 4,
            max_setup_terms: 10_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 100.0,
            maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9,
        };
        let mut network = CoupledModalSystem::new(vec![model], vec![], config, &gate).unwrap();
        let original = RigidBodyState::new(
            Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap();
        let mut bodies = [original];
        let properties = [MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()];
        let mut impacts = PreparedRigidModalImpacts::new(
            &network,
            vec![RigidModalImpactContact {
                id: 1,
                rigid_body: 0,
                arm_body_m: Vec3::ZERO,
                modal_point: ModalImpactPoint {
                    component: 0,
                    reference_position_world_m: Vec3::ZERO,
                    shapes_world: vec![Vec3::new(1.0, 0.0, 0.0)],
                },
                normal_world: Vec3::new(1.0, 0.0, 0.0),
                tangent_world: Vec3::new(0.0, 1.0, 0.0),
                activation_tolerance_m: 0.0,
                friction_coefficient: 0.0,
                normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
            }],
            RigidModalImpactConfig::default(),
            &gate,
        )
        .unwrap();
        let refused = network.apply_rigid_frictional_impacts_with_admission(
            &mut bodies,
            &properties,
            &mut impacts,
            &gate,
            |candidate| {
                assert_eq!(candidate[0].states()[0].velocity_m_sqrt_kg_per_s, -0.5);
                Err(invalid("physical cavity flow admission refused"))
            },
        );
        assert!(matches!(
            refused,
            Err(RigidModalImpactError::Coupling(
                ModalCouplingError::Invalid { .. }
            ))
        ));
        assert_eq!(bodies, [original]);
        assert_eq!(
            network.components()[0].states()[0].velocity_m_sqrt_kg_per_s,
            0.0
        );
        let refused = network.apply_rigid_frictional_impacts_with_admission(
            &mut bodies,
            &properties,
            &mut impacts,
            &gate,
            |_| {
                gate.request();
                Ok(())
            },
        );
        assert!(matches!(
            refused,
            Err(RigidModalImpactError::Coupling(
                ModalCouplingError::Cancelled
            ))
        ));
        assert_eq!(bodies, [original]);
        assert_eq!(
            network.components()[0].states()[0].velocity_m_sqrt_kg_per_s,
            0.0
        );
        let fresh_gate = CancelGate::new_clock_free();
        network
            .apply_rigid_frictional_impacts(&mut bodies, &properties, &mut impacts, &fresh_gate)
            .unwrap();
        assert_eq!(bodies[0].linear_momentum_world().x, -0.5);
        assert_eq!(
            network.components()[0].states()[0].velocity_m_sqrt_kg_per_s,
            -0.5
        );
        assert_eq!(network.samples_rendered(), 0);
        assert!(network.last_frame().is_none());
    }
}
