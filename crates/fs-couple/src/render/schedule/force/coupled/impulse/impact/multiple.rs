//! Simultaneous frictionless impacts through the original reciprocal ports.
//!
//! Every relative row B_i includes both attachment endpoints. The complete
//! Delassus matrix is W_ij = B_i B_j^T in the native unit-mass coordinates;
//! cross-contact/shared-mode terms are retained. One original fs-contact
//! complementarity solve supplies all impulses, followed by ONE B^T lambda
//! jump. There is no sequential scalar-impact ordering or unit-kick simulation.
//!
//! The original solver's redundant-row result determines velocity, not unique
//! individual loads. Its representative is never promoted to a local heat or
//! pressure allocation. This is an explicitly supplied instantaneous contact
//! set, not collision detection, persistent contact, friction or an adjoint.

use super::{
    CoupledModalSystem, ModalAcousticTimeModel, ModalCouplingError, ModalImpactContact,
    ModalImpactError, PreparedModalImpact, checked, finite, invalid, poll,
};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Fixed physical rows, complete contact mobility and reusable projection
/// scratch. Geometry and velocities are sampled from the continuing state at
/// every application. The original bounded contact solver still allocates its
/// candidate and returned receipt; no allocation-free solver claim is made.
pub struct PreparedModalImpacts {
    rows: Vec<PreparedModalImpact>,
    delassus: Vec<f64>,
    config: MultiImpactConfig,
    sources: Vec<NormalContactImpact>,
    before: Vec<Dd>,
    requested: Vec<Dd>,
    impulses: Vec<f64>,
}

/// Original joint contact response and independent actual-state audits.
/// Dissipation in `contacts` remains an unresolved physical loss account.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalMultiImpactReceipt {
    pub sample: u64,
    pub contacts: MultiImpactResult,
    pub gaps_m: Vec<f64>,
    pub relative_normal_speeds_before_m_s: Vec<f64>,
    pub relative_normal_speeds_after_m_s: Vec<f64>,
    /// Stable increment of the represented mechanical kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Sum of original impulses times actual midpoint relative velocities [J].
    pub contact_work_j: f64,
    /// Actual kinetic increment minus actual contact work [J].
    pub energy_residual_j: f64,
    /// Actual contact work minus the original solver's predicted work [J].
    pub contact_work_residual_j: f64,
    pub maximum_velocity_residual_m_s: f64,
    pub maximum_complementarity_residual_m_s: f64,
    pub maximum_relative_momentum_residual: f64,
}

/// No refusal publishes any of the participating native states.
#[derive(Debug)]
pub enum ModalMultiImpactError {
    Coupling(ModalCouplingError),
    Contact(MultiImpactError),
    Attachment {
        index: usize,
        source: ModalImpactError,
    },
    /// Measured endpoint speed disagrees with the original Delassus response.
    Velocity {
        index: usize,
        residual_m_s: f64,
        tolerance_m_s: f64,
    },
}

impl From<ModalCouplingError> for ModalMultiImpactError {
    fn from(value: ModalCouplingError) -> Self {
        Self::Coupling(value)
    }
}
impl From<MultiImpactError> for ModalMultiImpactError {
    fn from(value: MultiImpactError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for ModalMultiImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Coupling(error) => core::fmt::Display::fmt(error, f),
            Self::Contact(error) => core::fmt::Display::fmt(error, f),
            Self::Attachment { index, source } => write!(f, "modal impact {index}: {source}"),
            Self::Velocity {
                index,
                residual_m_s,
                tolerance_m_s,
            } => write!(
                f,
                "modal impact {index} velocity residual {residual_m_s:e} m/s exceeds {tolerance_m_s:e} m/s"
            ),
        }
    }
}
impl std::error::Error for ModalMultiImpactError {}

impl PreparedModalImpacts {
    /// Prepare at most 64 rows, within the original contact solver's capacity.
    /// The network's existing `max_setup_terms` also bounds n_contacts^2 times
    /// n_modes before allocating/assembling the full Gram mobility. Empty sets
    /// are explicit no-ops. An individual zero relative row remains invalid.
    pub fn new(
        system: &CoupledModalSystem,
        contacts: Vec<ModalImpactContact>,
        config: MultiImpactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalMultiImpactError> {
        poll(Some(gate))?;
        let count = contacts.len();
        if count > 64 || count > config.solver.max_variables {
            return Err(
                MultiImpactError::Invalid("modal contact count exceeds solver capacity").into(),
            );
        }
        let work = count
            .checked_mul(count)
            .and_then(|square| square.checked_mul(system.mode_count()))
            .ok_or_else(|| invalid("modal contact setup extent overflow"))?;
        if work > system.config.max_setup_terms {
            return Err(invalid("modal contact Gram work exceeds max_setup_terms").into());
        }
        // The original numerical owner admits its own configuration, including
        // the empty case, without a fabricated physical validation problem.
        let _ = solve_normal_impacts(&[], &[], config, || gate.is_requested())?;
        let mut rows = Vec::with_capacity(count);
        for (index, contact) in contacts.into_iter().enumerate() {
            rows.push(
                PreparedModalImpact::new(system, contact, gate)
                    .map_err(|source| ModalMultiImpactError::Attachment { index, source })?,
            );
        }
        let mut delassus = vec![0.0; count * count];
        for i in 0..count {
            for j in 0..=i {
                poll(Some(gate))?;
                let mut entry = Dd::ZERO;
                for (mode, (&a, &b)) in rows[i].column.iter().zip(&rows[j].column).enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    entry = entry + Dd::from_f64(a) * Dd::from_f64(b);
                }
                let entry = checked(entry)?;
                delassus[i * count + j] = entry;
                delassus[j * count + i] = entry;
            }
        }
        let sources = rows
            .iter()
            .map(|row| NormalContactImpact {
                gap: 0.0,
                activation_tolerance: row.activation_tolerance_m,
                relative_normal_speed: 0.0,
                restitution: row.restitution,
            })
            .collect();
        poll(Some(gate))?;
        Ok(Self {
            rows,
            delassus,
            config,
            sources,
            before: vec![Dd::ZERO; count],
            requested: vec![Dd::ZERO; system.mode_count()],
            impulses: vec![0.0; system.mode_count()],
        })
    }

    pub fn contact_count(&self) -> usize {
        self.rows.len()
    }

    /// Full row-major original relative mobility [1/kg].
    pub fn delassus_kg_inv(&self) -> &[f64] {
        &self.delassus
    }

    /// Endpoint incidence in original contact order, for physical-surface
    /// admission. Auxiliary acoustic coordinates are not new contact bodies.
    pub fn component_pairs(&self) -> impl ExactSizeIterator<Item = [usize; 2]> + '_ {
        self.rows.iter().map(PreparedModalImpact::components)
    }
}

impl CoupledModalSystem {
    /// Resolve all supplied contacts together and publish one native reset.
    /// Initially separating geometry-active rows participate in the same solve;
    /// a globally energy-injecting restitution choice refuses rather than being
    /// silently weakened. Redundant rows retain the original velocity-only
    /// determinacy result and its unresolved local load allocation.
    ///
    /// Actual per-coordinate momentum, contact velocities/complementarity and
    /// independently accumulated contact work must pass before publication.
    /// All original state/pressure/connection-force/energy gates remain active.
    /// Displacements, material/acoustic memories, held inputs, preceding frame
    /// and clock stay unchanged. Preparation scratch may change on refusal;
    /// accepted physical state never does.
    pub fn apply_normal_impacts(
        &mut self,
        impacts: &mut PreparedModalImpacts,
        gate: &CancelGate,
    ) -> Result<ModalMultiImpactReceipt, ModalMultiImpactError> {
        self.apply_normal_impacts_with_admission(impacts, gate, |_| Ok(()))
    }

    pub(crate) fn apply_normal_impacts_with_admission(
        &mut self,
        impacts: &mut PreparedModalImpacts,
        gate: &CancelGate,
        admit: impl FnOnce(&[ModalAcousticTimeModel]) -> Result<(), ModalCouplingError>,
    ) -> Result<ModalMultiImpactReceipt, ModalMultiImpactError> {
        poll(Some(gate))?;
        self.require_time_domain()?;
        if impacts.impulses.len() != self.mode_count() {
            return Err(invalid("prepared impacts no longer match the native mode count").into());
        }
        for (index, row) in impacts.rows.iter().enumerate() {
            let (gap, before) = row
                .sample(&self.models, gate)
                .map_err(|source| ModalMultiImpactError::Attachment { index, source })?;
            impacts.sources[index].gap = gap;
            impacts.sources[index].relative_normal_speed = checked(before)?;
            impacts.before[index] = before;
        }
        let contacts =
            solve_normal_impacts(&impacts.sources, &impacts.delassus, impacts.config, || {
                gate.is_requested()
            })?;
        for mode in 0..self.mode_count() {
            if mode % 64 == 0 {
                poll(Some(gate))?;
            }
            let mut requested = Dd::ZERO;
            for (row, &impulse) in impacts.rows.iter().zip(&contacts.normal_impulses) {
                let shape = row.column[mode];
                let term = Dd::from_f64(shape) * Dd::from_f64(impulse);
                if shape != 0.0 && impulse != 0.0 && term.hi == 0.0 && term.lo == 0.0 {
                    return Err(invalid("contact projection loses nonzero modal momentum").into());
                }
                requested = requested + term;
            }
            impacts.requested[mode] = requested;
            impacts.impulses[mode] = checked(requested)?;
        }
        let changes_state = contacts
            .normal_impulses
            .iter()
            .any(|&impulse| impulse != 0.0);
        let jump = if changes_state {
            Some(self.stage_modal_impulse(&impacts.impulses, Some(gate))?)
        } else {
            None
        };
        let after_models = if changes_state {
            &self.candidates
        } else {
            &self.models
        };
        let mut maximum_relative_momentum_residual =
            jump.map_or(0.0, |jump| jump.maximum_relative_momentum_residual);
        // Audit the FULL compensated reciprocal projection, including its
        // rounding tail, against represented endpoint increments. Contact
        // cancellation cannot excuse a swallowed small net momentum change.
        for (mode, (before, after)) in self
            .models
            .iter()
            .flat_map(|model| model.states())
            .zip(after_models.iter().flat_map(|model| model.states()))
            .enumerate()
        {
            if mode % 64 == 0 {
                poll(Some(gate))?;
            }
            let actual = Dd::from_f64(after.velocity_m_sqrt_kg_per_s)
                - Dd::from_f64(before.velocity_m_sqrt_kg_per_s);
            let requested = impacts.requested[mode];
            let scale = checked(actual)?.abs().max(checked(requested)?.abs());
            let residual = checked(actual - requested)?.abs();
            let relative = if scale == 0.0 { 0.0 } else { residual / scale };
            if relative > self.config.solve_relative_tolerance {
                return Err(ModalCouplingError::ImpulseResidual {
                    mode,
                    relative,
                    tolerance: self.config.solve_relative_tolerance,
                }
                .into());
            }
            maximum_relative_momentum_residual = maximum_relative_momentum_residual.max(relative);
        }
        let mut speeds_after = Vec::with_capacity(impacts.rows.len());
        let mut work = Dd::ZERO;
        let mut maximum_velocity_residual_m_s = 0.0_f64;
        let mut maximum_complementarity_residual_m_s = 0.0_f64;
        for (index, row) in impacts.rows.iter().enumerate() {
            let (_, after) = row
                .sample(after_models, gate)
                .map_err(|source| ModalMultiImpactError::Attachment { index, source })?;
            let after_f64 = checked(after)?;
            speeds_after.push(after_f64);
            let predicted = contacts.relative_normal_speeds_after[index];
            let velocity_residual = checked(after - Dd::from_f64(predicted))?;
            let tolerance = finite(
                impacts.config.solver.gradient_tolerance
                    + self.config.solve_relative_tolerance
                        * predicted
                            .abs()
                            .max(after_f64.abs())
                            .max(impacts.sources[index].relative_normal_speed.abs()),
            )?;
            if velocity_residual.abs() > tolerance {
                return Err(ModalMultiImpactError::Velocity {
                    index,
                    residual_m_s: velocity_residual,
                    tolerance_m_s: tolerance,
                });
            }
            maximum_velocity_residual_m_s =
                maximum_velocity_residual_m_s.max(velocity_residual.abs());
            let impulse = contacts.normal_impulses[index];
            if impacts.sources[index].gap <= impacts.sources[index].activation_tolerance {
                let gradient = checked(after - Dd::from_f64(contacts.restitution_targets[index]))?;
                let residual = if impulse > 0.0 {
                    gradient.abs()
                } else {
                    (-gradient).max(0.0)
                }
                .max((-after_f64).max(0.0));
                if residual > impacts.config.solver.gradient_tolerance {
                    return Err(MultiImpactError::UnresolvedContact { index, residual }.into());
                }
                maximum_complementarity_residual_m_s =
                    maximum_complementarity_residual_m_s.max(residual);
            }
            let mean = impacts.before[index] * Dd::from_f64(0.5) + after * Dd::from_f64(0.5);
            work = work + Dd::from_f64(impulse) * mean;
        }
        let contact_work_j = checked(work)?;
        // DD reduction permits an allowance on the actual NET work, without
        // borrowing tolerance from cancelling individual contact works.
        if contact_work_j > impacts.config.energy_tolerance {
            return Err(MultiImpactError::EnergyIncrease {
                work: contact_work_j,
                allowed: impacts.config.energy_tolerance,
            }
            .into());
        }
        let kinetic_energy_change_j = jump.map_or(0.0, |jump| jump.kinetic_energy_change_j);
        let energy_residual_j = checked(Dd::from_f64(kinetic_energy_change_j) - work)?;
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
        let receipt = ModalMultiImpactReceipt {
            sample: self.samples_rendered(),
            gaps_m: impacts.sources.iter().map(|source| source.gap).collect(),
            relative_normal_speeds_before_m_s: impacts
                .sources
                .iter()
                .map(|source| source.relative_normal_speed)
                .collect(),
            relative_normal_speeds_after_m_s: speeds_after,
            contacts,
            kinetic_energy_change_j,
            contact_work_j,
            energy_residual_j,
            contact_work_residual_j,
            maximum_velocity_residual_m_s,
            maximum_complementarity_residual_m_s,
            maximum_relative_momentum_residual,
        };
        admit(after_models)?;
        poll(Some(gate))?;
        if changes_state {
            self.publish_modal_impulse();
        }
        Ok(receipt)
    }
}
