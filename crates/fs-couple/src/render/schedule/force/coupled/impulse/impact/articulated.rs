//! Native articulated mechanisms exchanging impulses with continuing modal points.
//!
//! Every original free/prescribed ABA response contributes to the SAME full
//! contact mobility as the modal rows. One native Coulomb solve determines all
//! impulses, which transpose through those rows into original link wrenches
//! and modal momenta. Joint coupling, full tensors, offset COMs and prescribed
//! support reactions are retained. No rigid proxy or applied unit-event probe
//! substitutes for the original inverse-mass action.
//!
//! Selected witnesses and world contact frames are frozen at the event. This
//! is not collision discovery, a time stepper or a compliant-contact limit.
//! Modal spatial momentum extends only to the original basis's represented
//! symmetries. Supported/truncated modes do not acquire free-body symmetries.
//! Support work is an explicit source; impact loss has no automatic heat split.

use super::point::{ModalContactInput, ModalContactRows, components, product};
pub use super::rigid::ModalImpactPoint;
use super::{
    CoupledModalSystem, ModalAcousticTimeModel, ModalCouplingError, audit_modal_momentum, checked,
    finite, invalid, poll,
};
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactError, FrictionalContactResiduals,
    FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;
use fs_ga::{Mat3, Twist, Vec3 as SpatialVec3, Wrench};
use fs_math::dd::Dd;
use fs_mbd::Vec3;
use fs_mbd::articulated::contact_geometry::LinkContactPoint;
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, ImpulseVelocityChange,
    Kinematics, apply_free_floating_impulses, apply_prescribed_base_impulses, forward_kinematics,
    free_floating_impulse_delta, prescribed_base_impulse_delta,
};
use fs_mbd::articulated_contacts::mixed::MixedArticulatedImpulseResponse;
use fs_mbd::articulated_contacts::trajectory::mixed::MixedArticulatedTrajectoryState;

/// One original link material point receives +p; the modal point receives -p.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedModalImpactContact {
    pub id: u64,
    /// Original system/link indices and point relative to the link origin.
    pub point: LinkContactPoint,
    pub modal_point: ModalImpactPoint,
    pub normal_world: Vec3,
    pub tangent_world: Vec3,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedModalImpactConfig {
    pub maximum_systems: usize,
    pub maximum_total_links: usize,
    pub geometry_tolerance_m: f64,
    pub angular_impulse_tolerance_n_m_s: f64,
    /// Each mechanism's actual spatial increment minus applied/contact and
    /// prescribed-support impulses must pass independently, before reduction.
    pub linear_momentum_tolerance_n_s: f64,
    pub angular_momentum_tolerance_n_m_s: f64,
    pub contact: FrictionalContactConfig,
}
impl Default for ArticulatedModalImpactConfig {
    fn default() -> Self {
        Self {
            maximum_systems: 64,
            maximum_total_links: 1024,
            geometry_tolerance_m: 1e-10,
            angular_impulse_tolerance_n_m_s: 1e-10,
            linear_momentum_tolerance_n_s: 1e-9,
            angular_momentum_tolerance_n_m_s: 1e-9,
            contact: FrictionalContactConfig::default(),
        }
    }
}

/// Original modal rows are prepared once. Native link geometry and ABA
/// response are recomputed from the actual configuration at every event.
/// The original inverse-mass owner, contact solver and receipts allocate;
/// this API makes no allocation-free or cached-articulated-factor claim.
pub struct PreparedArticulatedModalImpacts {
    contacts: Vec<ArticulatedModalImpactContact>,
    rows: ModalContactRows,
    config: ArticulatedModalImpactConfig,
    delassus: Vec<f64>,
    directions: Vec<Wrench>,
    separation: Vec<[Dd; 3]>,
    before: Vec<[Dd; 3]>,
    requested_modal: Vec<Dd>,
    modal_impulses: Vec<f64>,
}

/// Native responses remain in caller system order; all contact arrays use
/// canonical contact-ID order. With UniqueVelocity, individual contact and
/// support loads/work shares remain representatives of the admitted solution.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedModalImpactReceipt {
    pub sample: u64,
    pub contact_ids: Vec<u64>,
    pub axes_world: Vec<[Vec3; 3]>,
    pub gaps_m: Vec<f64>,
    pub contact_mobility_kg_inv: Vec<f64>,
    pub contacts: FrictionalContactResult,
    pub articulated_impulses: Vec<MixedArticulatedImpulseResponse>,
    pub relative_velocities_before_m_s: Vec<[f64; 3]>,
    pub relative_velocities_after_m_s: Vec<[f64; 3]>,
    pub modal_kinetic_energy_change_j: f64,
    pub articulated_kinetic_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    /// Work of the passive contact pair; separate from driven base work.
    pub contact_work_j: f64,
    pub support_work_j: f64,
    /// Actual combined kinetic change minus contact and support work.
    pub energy_residual_j: f64,
    pub contact_work_residual_j: f64,
    /// Signed aggregate of original prescribed reactions, about world origin.
    pub support_impulse_world: Wrench,
    /// Actual mechanisms' spatial changes minus their physical contact loads
    /// and prescribed support reactions. These exclude modal spatial momentum.
    pub linear_momentum_residual_world: SpatialVec3,
    pub angular_momentum_residual_world: SpatialVec3,
    pub maximum_relative_modal_momentum_residual: f64,
    pub maximum_relative_articulated_velocity_increment_residual: f64,
    pub maximum_velocity_residual_m_s: f64,
    pub maximum_geometric_angular_impulse_defect_n_m_s: f64,
    pub actual_residuals: FrictionalContactResiduals,
}

#[derive(Debug)]
pub enum ArticulatedModalImpactError {
    Coupling(ModalCouplingError),
    Dynamics {
        system: usize,
        source: ArticulatedError,
    },
    Contact(FrictionalContactError),
    VelocityIncrement {
        system: usize,
        coordinate: usize,
        relative: f64,
        tolerance: f64,
    },
    Momentum {
        system: usize,
        coordinate: usize,
        residual: f64,
        tolerance: f64,
    },
    Velocity {
        contact: usize,
        axis: usize,
        residual_m_s: f64,
        tolerance_m_s: f64,
    },
}
impl From<ModalCouplingError> for ArticulatedModalImpactError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Coupling(error)
    }
}
impl From<FrictionalContactError> for ArticulatedModalImpactError {
    fn from(error: FrictionalContactError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for ArticulatedModalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native articulated/modal impact refused: {self:?}")
    }
}
impl std::error::Error for ArticulatedModalImpactError {}

impl PreparedArticulatedModalImpacts {
    pub fn new(
        system: &CoupledModalSystem,
        mut contacts: Vec<ArticulatedModalImpactContact>,
        config: ArticulatedModalImpactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ArticulatedModalImpactError> {
        poll(Some(gate))?;
        system.require_time_domain()?;
        if contacts.len() > 21
            || contacts.len() > config.contact.maximum_contacts
            || config.maximum_systems == 0
            || config.maximum_systems > 128
            || config.maximum_total_links == 0
            || config.maximum_total_links > 8192
            || [
                config.geometry_tolerance_m,
                config.angular_impulse_tolerance_n_m_s,
                config.linear_momentum_tolerance_n_s,
                config.angular_momentum_tolerance_n_m_s,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid(
                "bounded articulated/modal contacts and positive finite controls required",
            )
            .into());
        }
        let dimension = 3 * contacts.len();
        let work = dimension
            .checked_mul(dimension)
            .and_then(|n| n.checked_mul(system.mode_count()))
            .ok_or_else(|| invalid("articulated/modal modal preparation extent overflow"))?;
        if work > system.config.max_setup_terms {
            return Err(invalid("articulated/modal preparation exceeds max_setup_terms").into());
        }
        let _ = solve_frictional_contacts(&[], &[], config.contact, || gate.is_requested())?;
        contacts.sort_by_key(|contact| contact.id);
        for contact in &contacts {
            if contact.point.system_index >= config.maximum_systems
                || spatial_components(contact.point.point_body)
                    .iter()
                    .any(|value| !value.is_finite())
            {
                return Err(
                    invalid("articulated/modal system index or original link point").into(),
                );
            }
        }
        let rows = ModalContactRows::new(
            system,
            &contacts
                .iter()
                .map(|contact| ModalContactInput {
                    id: contact.id,
                    point: &contact.modal_point,
                    normal_world: contact.normal_world,
                    tangent_world: contact.tangent_world,
                    activation_tolerance_m: contact.activation_tolerance_m,
                    friction_coefficient: contact.friction_coefficient,
                    normal_law: contact.normal_law,
                })
                .collect::<Vec<_>>(),
            gate,
        )?;
        let count = contacts.len();
        Ok(Self {
            contacts,
            rows,
            config,
            delassus: vec![0.0; dimension * dimension],
            directions: vec![Wrench::default(); dimension],
            separation: vec![[Dd::ZERO; 3]; count],
            before: vec![[Dd::ZERO; 3]; count],
            requested_modal: vec![Dd::ZERO; system.mode_count()],
            modal_impulses: vec![0.0; system.mode_count()],
        })
    }

    pub fn contact_count(&self) -> usize {
        self.contacts.len()
    }

    pub fn modal_components(&self) -> impl Iterator<Item = usize> + '_ {
        self.contacts
            .iter()
            .map(|contact| contact.modal_point.component)
    }

    fn velocity(
        &self,
        index: usize,
        kinematics: &Kinematics,
        models: &[ModalAcousticTimeModel],
        gate: &CancelGate,
    ) -> Result<[Dd; 3], ModalCouplingError> {
        let mut result = self.rows.trace(index, models, true, gate)?;
        let twist = kinematics.body_twist[self.contacts[index].point.link_index];
        for axis in 0..3 {
            result[axis] = result[axis] + pairing(self.directions[3 * index + axis], twist);
            checked(result[axis])?;
        }
        Ok(result)
    }
}

struct Mechanism<'a> {
    model: &'a ArticulatedModel,
    base: BaseState,
    free: bool,
    position: &'a [f64],
    velocity: &'a [f64],
    before: Kinematics,
    zero: Vec<f64>,
    loads: Vec<Wrench>,
    requested: Vec<[Dd; 6]>,
    external_world: [Dd; 6],
}
impl Mechanism<'_> {
    fn inverse_mass(&self) -> Result<ImpulseVelocityChange, ArticulatedError> {
        if self.free {
            free_floating_impulse_delta(self.model, self.position, &self.zero, &self.loads)
        } else {
            prescribed_base_impulse_delta(self.model, self.position, &self.zero, &self.loads)
                .map(|delta| delta.velocity)
        }
    }

    fn apply(&self) -> Result<MixedArticulatedImpulseResponse, ArticulatedError> {
        if self.free {
            apply_free_floating_impulses(
                self.model,
                FreeFloatingBaseState::new(self.base.world_from_base, self.base.twist_body),
                self.position,
                self.velocity,
                &self.zero,
                &self.loads,
            )
            .map(MixedArticulatedImpulseResponse::Free)
        } else {
            apply_prescribed_base_impulses(
                self.model,
                self.base,
                self.position,
                self.velocity,
                &self.zero,
                &self.loads,
            )
            .map(MixedArticulatedImpulseResponse::Prescribed)
        }
    }
}

impl CoupledModalSystem {
    /// Resolve one simultaneous event against actual free and prescribed trees.
    /// Both native receivers publish together after all physical/law gates and
    /// a final cancellation check. Configuration, clocks, material/radiation
    /// memories and prior finite-interval receipts are never reset or advanced.
    /// The caller must select states at the same physical event; prescribed
    /// clocks retain their independent epochs rather than being resynchronized.
    ///
    /// Native link points are measured from link origins, not assumed COMs.
    /// ABA probes enforce configuration/inertia admission but no fictitious
    /// unit-event speed limit. The ACTUAL combined endpoint retains all native
    /// joint speed limits, support impulses/work and momentum admission.
    pub fn apply_articulated_frictional_impacts(
        &mut self,
        states: &mut [MixedArticulatedTrajectoryState],
        models: &[&ArticulatedModel],
        impacts: &mut PreparedArticulatedModalImpacts,
        gate: &CancelGate,
    ) -> Result<ArticulatedModalImpactReceipt, ArticulatedModalImpactError> {
        self.apply_articulated_frictional_impacts_with_admission(
            states,
            models,
            impacts,
            gate,
            |_| Ok(()),
        )
    }

    pub(crate) fn apply_articulated_frictional_impacts_with_admission(
        &mut self,
        states: &mut [MixedArticulatedTrajectoryState],
        models: &[&ArticulatedModel],
        impacts: &mut PreparedArticulatedModalImpacts,
        gate: &CancelGate,
        admit: impl FnOnce(&[ModalAcousticTimeModel]) -> Result<(), ModalCouplingError>,
    ) -> Result<ArticulatedModalImpactReceipt, ArticulatedModalImpactError> {
        poll(Some(gate))?;
        self.require_time_domain()?;
        if states.len() != models.len()
            || states.len() > impacts.config.maximum_systems
            || !impacts.rows.matches(self)
            || impacts.contacts.iter().any(|contact| {
                models
                    .get(contact.point.system_index)
                    .is_none_or(|model| contact.point.link_index >= model.link_count())
            })
        {
            return Err(invalid(
                "prepared articulated/modal contact no longer matches native receivers",
            )
            .into());
        }
        let dimension = 3 * impacts.contacts.len();
        let links = models
            .iter()
            .try_fold(0usize, |sum, model| sum.checked_add(model.link_count()))
            .ok_or_else(|| invalid("articulated/modal link extent overflow"))?;
        let work = links
            .checked_mul(dimension.checked_add(1).unwrap_or(usize::MAX))
            .and_then(|n| n.checked_mul(36))
            .and_then(|n| {
                dimension
                    .checked_mul(dimension)?
                    .checked_mul(self.mode_count())?
                    .checked_add(n)
            })
            .ok_or_else(|| invalid("articulated/modal response work overflow"))?;
        if links > impacts.config.maximum_total_links || work > self.config.max_setup_terms {
            return Err(
                invalid("articulated/modal original response exceeds link/setup bounds").into(),
            );
        }
        let mut mechanisms = Vec::with_capacity(states.len());
        for (index, (&model, state)) in models.iter().zip(states.iter()).enumerate() {
            poll(Some(gate))?;
            let (base, free, position, velocity) = match state {
                MixedArticulatedTrajectoryState::Free(state) => (
                    BaseState::prescribed(
                        state.base.world_from_base,
                        state.base.twist_body,
                        Twist::zero(),
                    ),
                    true,
                    state.joint_position.as_slice(),
                    state.joint_velocity.as_slice(),
                ),
                MixedArticulatedTrajectoryState::Prescribed(state) => {
                    if !state.time_s.is_finite() {
                        return Err(invalid("finite original prescribed clock required").into());
                    }
                    (
                        state.base,
                        false,
                        state.joint_position.as_slice(),
                        state.joint_velocity.as_slice(),
                    )
                }
            };
            let before = dynamics(index, forward_kinematics(model, base, position, velocity))?;
            mechanisms.push(Mechanism {
                model,
                base,
                free,
                position,
                velocity,
                before,
                zero: vec![0.0; model.dof_count()],
                loads: vec![Wrench::default(); model.link_count()],
                requested: vec![[Dd::ZERO; 6]; model.link_count()],
                external_world: [Dd::ZERO; 6],
            });
        }
        for index in 0..impacts.contacts.len() {
            poll(Some(gate))?;
            let contact = &impacts.contacts[index];
            let receiver = &mechanisms[contact.point.system_index];
            let pose = receiver.before.world_from_link[contact.point.link_index];
            let origin = spatial_components(pose.translation());
            let reference = components(contact.modal_point.reference_position_world_m);
            let arm = spatial_components(dynamics(
                contact.point.system_index,
                pose.rotation()
                    .rotate(contact.point.point_body)
                    .map_err(ArticulatedError::from),
            )?);
            let mut separation = impacts.rows.trace(index, &self.models, false, gate)?;
            for axis in 0..3 {
                let direction = dynamics(
                    contact.point.system_index,
                    pose.rotation()
                        .inverse()
                        .rotate(spatial(impacts.rows.axes[index][axis]))
                        .map_err(ArticulatedError::from),
                )?;
                impacts.directions[3 * index + axis] =
                    Wrench::new(contact.point.point_body.cross(direction), direction);
                for coordinate in 0..3 {
                    separation[axis] = separation[axis]
                        + Dd::from_f64(components(impacts.rows.axes[index][axis])[coordinate])
                            * (Dd::from_f64(origin[coordinate])
                                - Dd::from_f64(reference[coordinate])
                                + Dd::from_f64(arm[coordinate]));
                }
            }
            if checked(separation[1])?.hypot(checked(separation[2])?)
                > impacts.config.geometry_tolerance_m
            {
                return Err(invalid(
                    "original articulated/modal witnesses do not share a normal line",
                )
                .into());
            }
            impacts.separation[index] = separation;
            impacts.rows.sources[index].gap = checked(separation[0])?;
            let before = impacts.velocity(index, &receiver.before, &self.models, gate)?;
            impacts.before[index] = before;
            impacts.rows.sources[index].relative_velocity = checked3(before)?;
        }
        for (value, &modal) in impacts.delassus.iter_mut().zip(&impacts.rows.gram) {
            *value = checked(modal)?;
        }
        for column in 0..dimension {
            poll(Some(gate))?;
            let point = impacts.contacts[column / 3].point;
            let receiver = &mut mechanisms[point.system_index];
            receiver.loads[point.link_index] = impacts.directions[column];
            let delta = dynamics(point.system_index, receiver.inverse_mass())?;
            receiver.loads[point.link_index] = Wrench::default();
            for row in 0..dimension {
                let other = impacts.contacts[row / 3].point;
                if other.system_index == point.system_index {
                    impacts.delassus[row * dimension + column] = checked(
                        impacts.rows.gram[row * dimension + column]
                            + pairing(impacts.directions[row], delta.body[other.link_index]),
                    )?;
                }
            }
        }
        let contacts = solve_frictional_contacts(
            &impacts.rows.sources,
            &impacts.delassus,
            impacts.config.contact,
            || gate.is_requested(),
        )?;
        impacts.requested_modal.fill(Dd::ZERO);
        let mut maximum_geometric_angular_impulse_defect_n_m_s = 0.0_f64;
        for index in 0..impacts.contacts.len() {
            poll(Some(gate))?;
            let p = contacts.impulses[index];
            let couple = cross(impacts.separation[index], p.map(Dd::from_f64));
            let defect = checked(couple[0])?
                .hypot(checked(couple[1])?)
                .hypot(checked(couple[2])?);
            if defect > impacts.config.angular_impulse_tolerance_n_m_s {
                return Err(invalid(
                    "separated articulated/modal witnesses create an unresolved impulse couple",
                )
                .into());
            }
            maximum_geometric_angular_impulse_defect_n_m_s =
                maximum_geometric_angular_impulse_defect_n_m_s.max(defect);
            let point = impacts.contacts[index].point;
            let receiver = &mut mechanisms[point.system_index];
            let mut world_force = [Dd::ZERO; 3];
            for axis in 0..3 {
                for (value, direction) in receiver.requested[point.link_index]
                    .iter_mut()
                    .zip(impacts.directions[3 * index + axis].to_array())
                {
                    *value = *value + product(direction, p[axis])?;
                }
                for (value, direction) in world_force
                    .iter_mut()
                    .zip(components(impacts.rows.axes[index][axis]))
                {
                    *value = *value + product(direction, p[axis])?;
                }
                for (mode, value) in impacts.requested_modal.iter_mut().enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    *value =
                        *value + product(impacts.rows.columns[3 * index + axis][mode], p[axis])?;
                }
            }
            let pose = receiver.before.world_from_link[point.link_index];
            let arm_world = matvec(pose.rotation().matrix(), ddvec(point.point_body));
            add3(&mut receiver.external_world[3..], world_force);
            add3(
                &mut receiver.external_world[..3],
                cross(ddvec(pose.translation()), world_force),
            );
            add3(
                &mut receiver.external_world[..3],
                cross(arm_world, world_force),
            );
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
        let maximum_relative_modal_momentum_residual =
            audit_modal_momentum(self, after_models, &impacts.requested_modal, gate)?;
        let mut responses = Vec::with_capacity(mechanisms.len());
        let mut after_kinematics = Vec::with_capacity(mechanisms.len());
        let mut articulated_change = Dd::ZERO;
        let mut support_work = Dd::ZERO;
        let mut support_world = [Dd::ZERO; 6];
        let mut momentum_residual = [Dd::ZERO; 6];
        let mut maximum_increment_residual = 0.0_f64;
        for (index, receiver) in mechanisms.iter_mut().enumerate() {
            poll(Some(gate))?;
            for (load, requested) in receiver.loads.iter_mut().zip(&receiver.requested) {
                *load = wrench(*requested)?;
            }
            let delta = dynamics(index, receiver.inverse_mass())?;
            let response = dynamics(index, receiver.apply())?;
            let (base_after, joint_after, support_body) = match &response {
                MixedArticulatedImpulseResponse::Free(response) => (
                    BaseState::prescribed(
                        response.base_after.world_from_base,
                        response.base_after.twist_body,
                        Twist::zero(),
                    ),
                    response.joint_velocity_after.as_slice(),
                    Wrench::default(),
                ),
                MixedArticulatedImpulseResponse::Prescribed(response) => (
                    response.base_after,
                    response.joint_velocity_after.as_slice(),
                    response.support_impulse_body,
                ),
            };
            for (coordinate, ((old, new), expected)) in receiver
                .base
                .twist_body
                .to_array()
                .into_iter()
                .chain(receiver.velocity.iter().copied())
                .zip(
                    base_after
                        .twist_body
                        .to_array()
                        .into_iter()
                        .chain(joint_after.iter().copied()),
                )
                .zip(
                    delta
                        .base
                        .to_array()
                        .into_iter()
                        .chain(delta.joint.iter().copied()),
                )
                .enumerate()
            {
                let actual = Dd::from_f64(new) - Dd::from_f64(old);
                let scale = checked(actual)?.abs().max(expected.abs());
                let residual = checked(actual - Dd::from_f64(expected))?.abs();
                let relative = if scale == 0.0 { 0.0 } else { residual / scale };
                if relative > self.config.solve_relative_tolerance {
                    return Err(ArticulatedModalImpactError::VelocityIncrement {
                        system: index,
                        coordinate,
                        relative,
                        tolerance: self.config.solve_relative_tolerance,
                    });
                }
                maximum_increment_residual = maximum_increment_residual.max(relative);
            }
            let after = dynamics(
                index,
                forward_kinematics(receiver.model, base_after, receiver.position, joint_after),
            )?;
            let (change, actual) =
                mechanical_increment(receiver.model, &receiver.before, &after, gate)?;
            articulated_change = articulated_change + change;
            let support = world_wrench(
                support_body.to_array().map(Dd::from_f64),
                receiver.base.world_from_base,
            );
            let source_work = pairing(support_body, receiver.base.twist_body);
            support_work = support_work + source_work;
            for coordinate in 0..6 {
                support_world[coordinate] = support_world[coordinate] + support[coordinate];
                let residual =
                    actual[coordinate] - receiver.external_world[coordinate] - support[coordinate];
                let tolerance = if coordinate < 3 {
                    impacts.config.angular_momentum_tolerance_n_m_s
                } else {
                    impacts.config.linear_momentum_tolerance_n_s
                };
                if checked(residual)?.abs() > tolerance {
                    return Err(ArticulatedModalImpactError::Momentum {
                        system: index,
                        coordinate,
                        residual: checked(residual)?,
                        tolerance,
                    });
                }
                momentum_residual[coordinate] = momentum_residual[coordinate] + residual;
            }
            let mut applied_work = source_work;
            for (link, requested) in receiver.requested.iter().enumerate() {
                for ((old, new), p) in receiver.before.body_twist[link]
                    .to_array()
                    .into_iter()
                    .zip(after.body_twist[link].to_array())
                    .zip(requested)
                {
                    applied_work =
                        applied_work + *p * (Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new));
                }
            }
            audit_energy(self, change - applied_work, change, applied_work)?;
            after_kinematics.push(after);
            responses.push(response);
        }
        let mut modal_change = Dd::ZERO;
        for (old, new) in self
            .models
            .iter()
            .flat_map(|model| model.states())
            .zip(after_models.iter().flat_map(|model| model.states()))
        {
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
            let after = impacts.velocity(
                index,
                &after_kinematics[impacts.contacts[index].point.system_index],
                after_models,
                gate,
            )?;
            let rounded = checked3(after)?;
            for axis in 0..3 {
                let predicted = contacts.relative_velocities_after[index][axis];
                let residual = checked(after[axis] - Dd::from_f64(predicted))?;
                let tolerance = finite(
                    impacts.config.contact.velocity_tolerance
                        + self.config.solve_relative_tolerance
                            * rounded[axis]
                                .abs()
                                .max(predicted.abs())
                                .max(impacts.rows.sources[index].relative_velocity[axis].abs()),
                )?;
                if residual.abs() > tolerance {
                    return Err(ArticulatedModalImpactError::Velocity {
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
            if impacts.rows.sources[index].gap <= impacts.rows.sources[index].activation_tolerance {
                super::frictional::audit_contact(
                    impacts.rows.sources[index].friction_coefficient,
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
        let change = modal_change + articulated_change;
        let energy_residual = change - work - support_work;
        audit_energy(self, energy_residual, change, work + support_work)?;
        let contact_work_residual = work - Dd::from_f64(contacts.relative_work);
        audit_energy(
            self,
            contact_work_residual,
            work,
            Dd::from_f64(contacts.relative_work),
        )?;
        let receipt = ArticulatedModalImpactReceipt {
            sample: self.samples_rendered(),
            contact_ids: impacts
                .rows
                .sources
                .iter()
                .map(|source| source.id)
                .collect(),
            axes_world: impacts.rows.axes.clone(),
            gaps_m: impacts
                .rows
                .sources
                .iter()
                .map(|source| source.gap)
                .collect(),
            contact_mobility_kg_inv: impacts.delassus.clone(),
            contacts,
            articulated_impulses: responses,
            relative_velocities_before_m_s: impacts
                .rows
                .sources
                .iter()
                .map(|source| source.relative_velocity)
                .collect(),
            relative_velocities_after_m_s: after_velocities,
            modal_kinetic_energy_change_j: checked(modal_change)?,
            articulated_kinetic_energy_change_j: checked(articulated_change)?,
            kinetic_energy_change_j: checked(change)?,
            contact_work_j,
            support_work_j: checked(support_work)?,
            energy_residual_j: checked(energy_residual)?,
            contact_work_residual_j: checked(contact_work_residual)?,
            support_impulse_world: wrench(support_world)?,
            linear_momentum_residual_world: spatial_checked(
                momentum_residual[3..].try_into().unwrap(),
            )?,
            angular_momentum_residual_world: spatial_checked(
                momentum_residual[..3].try_into().unwrap(),
            )?,
            maximum_relative_modal_momentum_residual,
            maximum_relative_articulated_velocity_increment_residual: maximum_increment_residual,
            maximum_velocity_residual_m_s,
            maximum_geometric_angular_impulse_defect_n_m_s,
            actual_residuals,
        };
        admit(after_models)?;
        poll(Some(gate))?;
        // Native functions preserve shape and boundary type. Validate these
        // publication invariants before touching either accepted receiver.
        for (state, response) in states.iter().zip(&receipt.articulated_impulses) {
            let valid = match (state, response) {
                (
                    MixedArticulatedTrajectoryState::Free(state),
                    MixedArticulatedImpulseResponse::Free(response),
                ) => state.joint_velocity.len() == response.joint_velocity_after.len(),
                (
                    MixedArticulatedTrajectoryState::Prescribed(state),
                    MixedArticulatedImpulseResponse::Prescribed(response),
                ) => state.joint_velocity.len() == response.joint_velocity_after.len(),
                _ => false,
            };
            if !valid {
                return Err(invalid("native articulated impulse changed its state layout").into());
            }
        }
        drop(mechanisms);
        poll(Some(gate))?;
        if changes_modal {
            self.publish_modal_impulse();
        }
        for (state, response) in states.iter_mut().zip(&receipt.articulated_impulses) {
            match (state, response) {
                (
                    MixedArticulatedTrajectoryState::Free(state),
                    MixedArticulatedImpulseResponse::Free(response),
                ) => {
                    state.base = response.base_after;
                    state
                        .joint_velocity
                        .copy_from_slice(&response.joint_velocity_after);
                }
                (
                    MixedArticulatedTrajectoryState::Prescribed(state),
                    MixedArticulatedImpulseResponse::Prescribed(response),
                ) => {
                    state.base = response.base_after;
                    state
                        .joint_velocity
                        .copy_from_slice(&response.joint_velocity_after);
                }
                _ => unreachable!("preflight verified native boundary type"),
            }
        }
        Ok(receipt)
    }
}

fn dynamics<T>(
    system: usize,
    value: Result<T, ArticulatedError>,
) -> Result<T, ArticulatedModalImpactError> {
    value.map_err(|source| ArticulatedModalImpactError::Dynamics { system, source })
}
fn spatial(value: Vec3) -> SpatialVec3 {
    SpatialVec3::new(value.x, value.y, value.z)
}
fn spatial_components(value: SpatialVec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn ddvec(value: SpatialVec3) -> [Dd; 3] {
    spatial_components(value).map(Dd::from_f64)
}
fn checked3(value: [Dd; 3]) -> Result<[f64; 3], ModalCouplingError> {
    Ok([checked(value[0])?, checked(value[1])?, checked(value[2])?])
}
fn spatial_checked(value: [Dd; 3]) -> Result<SpatialVec3, ModalCouplingError> {
    let value = checked3(value)?;
    Ok(SpatialVec3::new(value[0], value[1], value[2]))
}
fn wrench(value: [Dd; 6]) -> Result<Wrench, ModalCouplingError> {
    Ok(Wrench::new(
        spatial_checked([value[0], value[1], value[2]])?,
        spatial_checked([value[3], value[4], value[5]])?,
    ))
}
fn pairing(wrench: Wrench, twist: Twist) -> Dd {
    wrench
        .to_array()
        .into_iter()
        .zip(twist.to_array())
        .fold(Dd::ZERO, |sum, (p, v)| {
            sum + Dd::from_f64(p) * Dd::from_f64(v)
        })
}
fn add3(sum: &mut [Dd], value: [Dd; 3]) {
    for (sum, value) in sum.iter_mut().zip(value) {
        *sum = *sum + value;
    }
}
fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn matvec(matrix: Mat3, value: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|row| {
        (0..3).fold(Dd::ZERO, |sum, col| {
            sum + Dd::from_f64(matrix.m[3 * row + col]) * value[col]
        })
    })
}
fn world_wrench(body: [Dd; 6], pose: fs_ga::Se3) -> [Dd; 6] {
    let force = matvec(pose.rotation().matrix(), [body[3], body[4], body[5]]);
    let torque = matvec(pose.rotation().matrix(), [body[0], body[1], body[2]]);
    let orbital = cross(ddvec(pose.translation()), force);
    [
        torque[0] + orbital[0],
        torque[1] + orbital[1],
        torque[2] + orbital[2],
        force[0],
        force[1],
        force[2],
    ]
}

// Polarized actual changes at link COMs. Keep local spin, COM arm and world
// orbital moments separate until the signed reduction; never subtract rounded
// whole-tree momenta or total energies, including under a distant world origin.
fn mechanical_increment(
    model: &ArticulatedModel,
    before: &Kinematics,
    after: &Kinematics,
    gate: &CancelGate,
) -> Result<(Dd, [Dd; 6]), ModalCouplingError> {
    let mut energy = Dd::ZERO;
    let mut momentum = [Dd::ZERO; 6];
    for (index, link) in model.links().iter().enumerate() {
        poll(Some(gate))?;
        let inertia = link.inertia();
        let old = before.body_twist[index];
        let new = after.body_twist[index];
        let old_arm = cross(ddvec(old.angular), ddvec(inertia.center_of_mass()));
        let new_arm = cross(ddvec(new.angular), ddvec(inertia.center_of_mass()));
        let old_com = core::array::from_fn::<_, 3, _>(|i| ddvec(old.linear)[i] + old_arm[i]);
        let new_com = core::array::from_fn::<_, 3, _>(|i| ddvec(new.linear)[i] + new_arm[i]);
        let delta_w = core::array::from_fn(|i| ddvec(new.angular)[i] - ddvec(old.angular)[i]);
        let linear =
            core::array::from_fn(|i| Dd::from_f64(inertia.mass()) * (new_com[i] - old_com[i]));
        let spin = matvec(inertia.inertia_com(), delta_w);
        let pose = before.world_from_link[index];
        let linear_world = matvec(pose.rotation().matrix(), linear);
        add3(&mut momentum[3..], linear_world);
        add3(&mut momentum[..3], matvec(pose.rotation().matrix(), spin));
        add3(
            &mut momentum[..3],
            matvec(
                pose.rotation().matrix(),
                cross(ddvec(inertia.center_of_mass()), linear),
            ),
        );
        add3(
            &mut momentum[..3],
            cross(ddvec(pose.translation()), linear_world),
        );
        for axis in 0..3 {
            energy = energy
                + linear[axis]
                    * (old_com[axis] * Dd::from_f64(0.5) + new_com[axis] * Dd::from_f64(0.5))
                + spin[axis]
                    * (Dd::from_f64(0.5 * spatial_components(old.angular)[axis])
                        + Dd::from_f64(0.5 * spatial_components(new.angular)[axis]));
        }
    }
    checked(energy)?;
    for value in momentum {
        checked(value)?;
    }
    Ok((energy, momentum))
}

fn audit_energy(
    system: &CoupledModalSystem,
    residual: Dd,
    first: Dd,
    second: Dd,
) -> Result<(), ModalCouplingError> {
    let residual_j = checked(residual)?;
    let scale = checked(first)?.abs().max(checked(second)?.abs());
    let tolerance_j = finite(
        system.config.energy_absolute_tolerance_j + system.config.energy_relative_tolerance * scale,
    )?;
    if residual_j.abs() > tolerance_j {
        return Err(ModalCouplingError::EnergyBalance {
            residual_j,
            tolerance_j,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
    use crate::render::schedule::force::coupled::ModalCouplingConfig;
    use fs_ga::Se3;
    use fs_math::c64::C64;
    use fs_mbd::articulated::motion::FreeFloatingMotionState;
    use fs_mbd::articulated::{JointModel, Link, SpatialInertia};

    #[test]
    fn late_admission_and_final_cancellation_publish_neither_receiver() {
        let gate = CancelGate::new_clock_free();
        let modal = ModalAcousticTimeModel::try_new_with_free_coordinates(
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
        let mut network = CoupledModalSystem::new(vec![modal], vec![], config, &gate).unwrap();
        let model = ArticulatedModel::new(vec![Link::new(
            "body",
            None,
            Se3::identity(),
            JointModel::FIXED,
            SpatialInertia::new(
                1.0,
                SpatialVec3::default(),
                Mat3 {
                    m: [1., 0., 0., 0., 1., 0., 0., 0., 1.],
                },
            )
            .unwrap(),
        )])
        .unwrap();
        let original = MixedArticulatedTrajectoryState::Free(FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::identity(),
                Twist::new(SpatialVec3::default(), SpatialVec3::new(-1., 0., 0.)),
            ),
            joint_position: vec![],
            joint_velocity: vec![],
        });
        let mut states = [original.clone()];
        let mut impacts = PreparedArticulatedModalImpacts::new(
            &network,
            vec![ArticulatedModalImpactContact {
                id: 1,
                point: LinkContactPoint {
                    system_index: 0,
                    link_index: 0,
                    point_body: SpatialVec3::default(),
                },
                modal_point: ModalImpactPoint {
                    component: 0,
                    reference_position_world_m: Vec3::ZERO,
                    shapes_world: vec![Vec3::new(1., 0., 0.)],
                },
                normal_world: Vec3::new(1., 0., 0.),
                tangent_world: Vec3::new(0., 1., 0.),
                activation_tolerance_m: 0.,
                friction_coefficient: 0.,
                normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0. },
            }],
            ArticulatedModalImpactConfig::default(),
            &gate,
        )
        .unwrap();
        let refused = network.apply_articulated_frictional_impacts_with_admission(
            &mut states,
            &[&model],
            &mut impacts,
            &gate,
            |candidate| {
                assert_eq!(candidate[0].states()[0].velocity_m_sqrt_kg_per_s, -0.5);
                Err(invalid("physical cavity flow admission refused"))
            },
        );
        assert!(matches!(
            refused,
            Err(ArticulatedModalImpactError::Coupling(
                ModalCouplingError::Invalid { .. }
            ))
        ));
        assert_eq!(states, [original.clone()]);
        assert_eq!(
            network.components()[0].states()[0].velocity_m_sqrt_kg_per_s,
            0.
        );
        let refused = network.apply_articulated_frictional_impacts_with_admission(
            &mut states,
            &[&model],
            &mut impacts,
            &gate,
            |_| {
                gate.request();
                Ok(())
            },
        );
        assert!(matches!(
            refused,
            Err(ArticulatedModalImpactError::Coupling(
                ModalCouplingError::Cancelled
            ))
        ));
        assert_eq!(states, [original]);
        assert_eq!(
            network.components()[0].states()[0].velocity_m_sqrt_kg_per_s,
            0.
        );
        network
            .apply_articulated_frictional_impacts(
                &mut states,
                &[&model],
                &mut impacts,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let MixedArticulatedTrajectoryState::Free(after) = &states[0] else {
            panic!("boundary changed")
        };
        assert_eq!(after.base.twist_body.linear.x, -0.5);
        assert_eq!(
            network.components()[0].states()[0].velocity_m_sqrt_kg_per_s,
            -0.5
        );
        assert_eq!(network.samples_rendered(), 0);
        assert!(network.last_frame().is_none());
    }
}
