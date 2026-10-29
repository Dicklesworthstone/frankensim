//! Simultaneous retained shear on the original coupled modal network.
//!
//! Every contact has two work-conjugate rows in a fixed physical tangent frame.
//! Setup condenses the full response `C_ij = b_i^T D_network b_j` once. The
//! shared constitutive owner solves all contacts against that matrix, starting
//! every trial from the same accepted histories. One sum of contact efforts
//! advances the original network; actual body motion must then close every
//! local work interval and the complete system balance before publication.
//!
//! Normal loads and geometry are prescribed during a sample. This is neither
//! contact discovery nor a vector extension of scalar Hertz--Mindlin memory.
//! Nonunique constitutive roots retain the shared solver's disclosed selection.

use super::super::super::network_response;
use super::super::{
    CoupledModalSystem, ModalAttachment, ModalCouplingError, ModalFrictionConfig,
    ModalFrictionError, dot, finite, invalid, limit, poll, scalar,
};
use super::planar::PlanarTangentialModalDrive;
use fs_exec::CancelGate;
use fs_math::det;
use fs_tribo::tangential::multiple::{
    MultipleTangentialConfig, MultipleTangentialError, PreparedTangentialContacts,
    TangentialContactInput,
};
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

/// All retained contacts and their original finite thermal network accept together.
#[path = "multiple/heat.rs"]
pub mod heat;

/// One two-axis material contact and its unchanged physical body projections.
pub struct ModalTangentialContact<L> {
    /// Both axes name the same positive-side component and modal basis.
    pub left: [ModalAttachment; 2],
    /// A missing negative side is a prescribed support; drive work stays explicit.
    pub right: Option<[ModalAttachment; 2]>,
    /// Original law, including any source binding and constitutive applicability.
    pub law: L,
    /// Per-contact physical limits and original constitutive root budgets.
    pub limits: ModalFrictionConfig,
}

/// Explicit joint solve and setup budgets, in addition to the original network's.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiTangentialModalConfig {
    pub coupling: MultipleTangentialConfig,
    /// Bound `d * (d + 1) * n`, with `d = 2 * contacts` and `n` original modes.
    /// This bounds the added projection/storage work. Original response solves
    /// and matrix admission keep their existing owners and budgets.
    pub maximum_projection_terms: usize,
}

/// A contact's work measured against the same accepted body endpoints as all others.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalTangentialContactFrame {
    pub traction_n: [f64; 2],
    pub slip_distance_m: [f64; 2],
    pub compatibility_residual_m: [f64; 2],
    pub contact_energy_j: f64,
    pub stored_change_j: f64,
    pub dissipated_j: f64,
    pub drive_work_j: f64,
    pub work_residual_j: f64,
    pub regime: FrictionRegime,
}

/// One accepted sample for the entire contact set and original modal network.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiTangentialModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub contact_energy_j: f64,
    pub contact_stored_change_j: f64,
    pub network_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub external_work_j: f64,
    pub drive_work_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Original construction order; no entry describes a partially accepted trial.
    pub contacts: Vec<ModalTangentialContactFrame>,
    pub sweeps: usize,
    pub law_evaluations: usize,
}

/// Preserve the original network/constitutive distinction on failed transactions.
#[derive(Debug)]
pub enum MultiTangentialModalError {
    Mechanical(ModalFrictionError),
    Interface(MultipleTangentialError),
}

impl core::fmt::Display for MultiTangentialModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Mechanical(error) => write!(f, "simultaneous modal shear: {error}"),
            Self::Interface(error) => write!(f, "simultaneous modal shear: {error}"),
        }
    }
}
impl std::error::Error for MultiTangentialModalError {}
impl From<ModalFrictionError> for MultiTangentialModalError {
    fn from(error: ModalFrictionError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<ModalCouplingError> for MultiTangentialModalError {
    fn from(error: ModalCouplingError) -> Self {
        Self::Mechanical(error.into())
    }
}
impl From<MultipleTangentialError> for MultiTangentialModalError {
    fn from(error: MultipleTangentialError) -> Self {
        Self::Interface(error)
    }
}

/// Retained contact histories interacting through one complete modal response.
///
/// Laws may have different parameters and sources. A caller-owned enum can
/// implement `TangentialLaw` when different constitutive families are required.
/// Each history remains paired with its original law and fixed attachment rows.
/// Setup retains all `2N * 2N` compliance entries; stepping neither diagonalizes
/// the response nor advances the bodies once per contact. Candidate histories
/// allocate according to their original law; no allocation-free claim is made.
pub struct MultiTangentialModalSystem<L: TangentialLaw> {
    network: CoupledModalSystem,
    laws: Vec<L>,
    states: Vec<L::State>,
    limits: Vec<ModalFrictionConfig>,
    columns: Vec<[Vec<f64>; 2]>,
    response: PreparedTangentialContacts,
    forces: Vec<f64>,
    last: Option<MultiTangentialModalFrame>,
}

struct MultiTangentialModalCandidate<S> {
    frame: MultiTangentialModalFrame,
    states: Vec<S>,
}

impl<L: TangentialLaw> MultiTangentialModalSystem<L> {
    /// Admit the complete contact set before allocating its quadratic response.
    pub fn new(
        network: CoupledModalSystem,
        contacts: Vec<ModalTangentialContact<L>>,
        config: MultiTangentialModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, MultiTangentialModalError> {
        poll(Some(gate))?;
        config.coupling.validate()?;
        if contacts.is_empty() || contacts.len() > config.coupling.maximum_contacts {
            return Err(invalid("contact count exceeds the nonempty joint shear budget").into());
        }
        let dimensions = contacts
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("tangential response extent overflow"))?;
        let matrix_terms = dimensions
            .checked_mul(dimensions)
            .ok_or_else(|| invalid("tangential matrix extent overflow"))?;
        if matrix_terms > config.coupling.maximum_matrix_terms {
            return Err(invalid("tangential matrix exceeds its explicit storage budget").into());
        }
        let terms = dimensions
            .checked_add(1)
            .and_then(|n| n.checked_mul(dimensions))
            .and_then(|n| n.checked_mul(network.mode_count()))
            .ok_or_else(|| invalid("tangential projection work overflow"))?;
        if terms > config.maximum_projection_terms {
            return Err(invalid("tangential projections exceed the setup work budget").into());
        }
        let mut columns = Vec::with_capacity(contacts.len());
        let mut laws = Vec::with_capacity(contacts.len());
        let mut states = Vec::with_capacity(contacts.len());
        let mut limits = Vec::with_capacity(contacts.len());
        let mut contact_energy = 0.0;
        for contact in contacts {
            poll(Some(gate))?;
            if contact.left[0].component != contact.left[1].component
                || contact
                    .right
                    .as_ref()
                    .is_some_and(|r| r[0].component != r[1].component)
            {
                return Err(
                    invalid("both tangent axes must retain each original contact body").into(),
                );
            }
            columns.push([
                scalar::admit_attachments(
                    &network,
                    &contact.left[0],
                    contact.right.as_ref().map(|r| &r[0]),
                    contact.limits,
                )?,
                scalar::admit_attachments(
                    &network,
                    &contact.left[1],
                    contact.right.as_ref().map(|r| &r[1]),
                    contact.limits,
                )?,
            ]);
            let state = contact.law.initial_state();
            let energy = finite(contact.law.stored_energy_j(&state))?;
            if energy < 0.0 {
                return Err(invalid("reversible contact storage must be nonnegative").into());
            }
            contact_energy = finite(contact_energy + energy)?;
            states.push(state);
            laws.push(contact.law);
            limits.push(contact.limits);
        }
        limit(
            "initial modal network and all contact storage",
            finite(network.total_energy_j()? + contact_energy)?,
            network.config.maximum_total_energy_j,
        )?;
        let mut compliance = vec![0.0; matrix_terms];
        for j in 0..dimensions {
            poll(Some(gate))?;
            let response = network_response(&network, &columns[j / 2][j % 2], gate)?;
            for i in 0..dimensions {
                compliance[i * dimensions + j] = dot(&columns[i / 2][i % 2], &response)?;
            }
        }
        let response =
            PreparedTangentialContacts::new(laws.len(), compliance, config.coupling, || {
                gate.is_requested()
            })?;
        poll(Some(gate))?;
        Ok(Self {
            forces: vec![0.0; network.mode_count()],
            network,
            laws,
            states,
            limits,
            columns,
            response,
            last: None,
        })
    }

    #[must_use]
    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }
    #[must_use]
    pub fn laws(&self) -> &[L] {
        &self.laws
    }
    #[must_use]
    pub fn states(&self) -> &[L::State] {
        &self.states
    }
    /// Full row-major physical response, with two consecutive rows per contact.
    #[must_use]
    pub fn compliance_m_per_n(&self) -> &[f64] {
        self.response.compliance_m_per_n()
    }
    #[must_use]
    pub const fn last_frame(&self) -> Option<&MultiTangentialModalFrame> {
        self.last.as_ref()
    }
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> {
        self.laws
            .iter()
            .zip(&self.states)
            .try_fold(self.network.total_energy_j()?, |sum, (law, state)| {
                finite(sum + law.stored_energy_j(state))
            })
    }

    pub fn step(
        &mut self,
        external: &[f64],
        drives: &[PlanarTangentialModalDrive],
    ) -> Result<&MultiTangentialModalFrame, MultiTangentialModalError> {
        self.step_inner(external, drives, None)
    }

    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        drives: &[PlanarTangentialModalDrive],
        gate: &CancelGate,
    ) -> Result<&MultiTangentialModalFrame, MultiTangentialModalError> {
        self.step_inner(external, drives, Some(gate))
    }

    fn step_inner(
        &mut self,
        external: &[f64],
        drives: &[PlanarTangentialModalDrive],
        gate: Option<&CancelGate>,
    ) -> Result<&MultiTangentialModalFrame, MultiTangentialModalError> {
        let candidate = self.stage_candidate(external, drives, gate)?;
        Ok(self.publish_candidate(candidate))
    }

    fn stage_candidate(
        &mut self,
        external: &[f64],
        drives: &[PlanarTangentialModalDrive],
        gate: Option<&CancelGate>,
    ) -> Result<MultiTangentialModalCandidate<L::State>, MultiTangentialModalError> {
        poll(gate)?;
        if drives.len() != self.laws.len()
            || drives.iter().any(|drive| {
                drive.speed_m_s.iter().any(|v| !v.is_finite())
                    || !drive.normal_force_n.is_finite()
                    || drive.normal_force_n < 0.0
                    || drive
                        .contact_temperature_k
                        .is_some_and(|t| !t.is_finite() || t <= 0.0)
                    || drive
                        .nominal_pressure_pa
                        .is_some_and(|p| !p.is_finite() || p < 0.0)
            })
        {
            return Err(invalid("every contact requires one finite physical drive").into());
        }
        let network_before = self.network.total_energy_j()?;
        self.network.prepare_prediction(external, gate)?;
        let mut sites = Vec::with_capacity(self.laws.len());
        for (i, drive) in drives.iter().enumerate() {
            poll(gate)?;
            let limits = self.limits[i];
            let mut free_slip = [0.0; 2];
            for axis in 0..2 {
                free_slip[axis] = finite(
                    drive.speed_m_s[axis] * self.network.dt
                        - dot(&self.columns[i][axis], &self.network.free_delta)?,
                )?;
            }
            sites.push(TangentialContactInput {
                law: &self.laws[i],
                state: &self.states[i],
                input: TangentialStepInput {
                    normal_n: drive.normal_force_n,
                    slip_increment_m: free_slip,
                    dt_s: self.network.dt,
                    contact_temperature_k: drive.contact_temperature_k,
                    nominal_pressure_pa: drive.nominal_pressure_pa,
                    compliance_m_per_n: 0.0,
                    solve: TangentialSolveConfig {
                        max_iterations: limits.max_iterations,
                        maximum_force_n: limits.maximum_force_n,
                        force_absolute_tolerance_n: limits.force_absolute_tolerance_n,
                        force_relative_tolerance: limits.force_relative_tolerance,
                    },
                },
            });
        }
        let result = self
            .response
            .advance(&sites, || gate.is_some_and(CancelGate::is_requested))?;
        poll(gate)?;
        self.forces.copy_from_slice(external);
        for (i, step) in result.contacts.iter().enumerate() {
            poll(gate)?;
            limit(
                "simultaneous contact traction",
                finite(det::hypot(step.traction_n[0], step.traction_n[1]))?,
                self.limits[i].maximum_force_n,
            )?;
            for ((force, bx), by) in self
                .forces
                .iter_mut()
                .zip(&self.columns[i][0])
                .zip(&self.columns[i][1])
            {
                *force = finite(*force + bx * step.traction_n[0] + by * step.traction_n[1])?;
            }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let mut frames = Vec::with_capacity(result.contacts.len());
        let mut storage = 0.0;
        let mut storage_change = 0.0;
        let mut heat = 0.0;
        let mut drive_work = 0.0;
        for (i, step) in result.contacts.iter().enumerate() {
            poll(gate)?;
            let limits = self.limits[i];
            let mut slip = [0.0; 2];
            let mut mismatch = [0.0; 2];
            let mut advance = [0.0; 2];
            for axis in 0..2 {
                advance[axis] = finite(drives[i].speed_m_s[axis] * self.network.dt)?;
                slip[axis] =
                    finite(advance[axis] - dot(&self.columns[i][axis], &self.network.free_delta)?)?;
                mismatch[axis] = finite(slip[axis] - step.slip_increment_m[axis])?;
                for residual in [mismatch[axis], step.compatibility_residual_m[axis]] {
                    limit(
                        "simultaneous contact actual slip compatibility",
                        finite(residual)?.abs(),
                        limits.sticking_tolerance_m,
                    )?;
                }
            }
            limit(
                "simultaneous contact slip speed",
                finite(det::hypot(slip[0], slip[1]) / self.network.dt)?,
                limits.maximum_slip_speed_m_s,
            )?;
            let before = finite(self.laws[i].stored_energy_j(&self.states[i]))?;
            let after = finite(self.laws[i].stored_energy_j(&step.state))?;
            let loss = finite(step.dissipated_j)?;
            if before < 0.0 || after < 0.0 || loss < 0.0 {
                return Err(invalid("every contact requires nonnegative storage and heat").into());
            }
            let change = finite(after - before)?;
            let supplied = dot(&step.traction_n, &advance)?;
            let contact_work = dot(&step.traction_n, &slip)?;
            let residual = finite(change + loss - contact_work)?;
            let scale = before
                .max(after)
                .max(loss)
                .max(contact_work.abs())
                .max(supplied.abs());
            let tolerance = finite(
                self.network.config.energy_absolute_tolerance_j
                    + self.network.config.energy_relative_tolerance * scale,
            )?;
            for checked in [
                residual,
                finite(change - step.stored_change_j)?,
                finite(step.closure_residual_j)?,
            ] {
                if checked.abs() > tolerance {
                    return Err(ModalCouplingError::EnergyBalance {
                        residual_j: checked,
                        tolerance_j: tolerance,
                    }
                    .into());
                }
            }
            storage = finite(storage + after)?;
            storage_change = finite(storage_change + change)?;
            heat = finite(heat + loss)?;
            drive_work = finite(drive_work + supplied)?;
            frames.push(ModalTangentialContactFrame {
                traction_n: step.traction_n,
                slip_distance_m: slip,
                compatibility_residual_m: mismatch,
                contact_energy_j: after,
                stored_change_j: change,
                dissipated_j: loss,
                drive_work_j: supplied,
                work_residual_j: residual,
                regime: step.regime,
            });
        }
        let original = &self.network.staged_frame;
        let network_after = finite(original.total_energy_j())?;
        let network_loss = finite(original.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        let after = finite(network_after + storage)?;
        limit(
            "modal network and simultaneous contact storage",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        let residual = finite(
            (network_after - network_before) + storage_change + network_loss + heat
                - external_work
                - drive_work,
        )?;
        let scale = after
            .max(network_before)
            .max(storage_change.abs())
            .max(network_loss)
            .max(heat)
            .max(external_work.abs())
            .max(drive_work.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let accepted = MultiTangentialModalFrame {
            sample: original.sample,
            observer_pressure_pa: original.observer_pressure_pa,
            network_energy_j: network_after,
            contact_energy_j: storage,
            contact_stored_change_j: storage_change,
            network_dissipation_j: network_loss,
            friction_dissipation_j: heat,
            external_work_j: external_work,
            drive_work_j: drive_work,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            contacts: frames,
            sweeps: result.sweeps,
            law_evaluations: result.law_evaluations,
        };
        poll(gate)?;
        Ok(MultiTangentialModalCandidate {
            frame: accepted,
            states: result.contacts.into_iter().map(|step| step.state).collect(),
        })
    }

    fn publish_candidate(
        &mut self,
        candidate: MultiTangentialModalCandidate<L::State>,
    ) -> &MultiTangentialModalFrame {
        self.network.publish_staged();
        self.states = candidate.states;
        self.last.insert(candidate.frame)
    }
}
