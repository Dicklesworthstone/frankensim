//! Whole closed reacting intervals on the original state-entropy basis.
//!
//! Let F be the shared thermal/mixing state function and alpha_i its omitted
//! species constants. The original source entropy changes by
//! dS = dF + sum_i(alpha_i * dm_i / M_i). Chemical species changes are physical:
//! only dm_i - dm_i,chem is a numerical continuity defect. The admitted change
//! is therefore dF + sum_i(alpha_i * dm_i,chem / M_i), not dF alone.

mod exterior;
pub use exterior::ReactiveExteriorEntropyReport;

use super::advection::ExteriorEntropyAccount;
use super::{
    Dd, EntropyTracker, FlowError, MAX_SPECIES, MixtureEuler1d, MolecularEntropyControl,
    poll, positive, scalar,
};

/// Original-to-final state entropy of a complete closed reacting-flow step or interval.
///
/// The original chemical species source remains part of the admitted entropy.
/// Only its numerical species-continuity defect is reported separately. This
/// measures the complete split endpoint; summing local owner allowances or
/// chemistry-only entropy changes does not replace this admission.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReactiveEntropyReport {
    /// Original source entropy change with only the linear source-reference
    /// contribution of numerical species-continuity defects removed [J/K].
    /// Actual chemical conversion, thermal change and mixing remain included.
    pub entropy_change_j_k: f64,
    /// Original source entropy change including that continuity defect [J/K].
    pub source_entropy_change_j_k: f64,
    /// Sum_i alpha_i * (dm_i - dm_i,chem) / M_i [J/K]. Chemical species changes
    /// are actual represented before/after changes at the original cell owner.
    pub reference_entropy_residual_j_k: f64,
    /// Largest species continuity residual divided by its original mass plus
    /// the sum of absolute chemical mass changes over both halves and all cells.
    /// Zero original inventory and zero chemical throughput require zero drift.
    pub maximum_relative_species_residual: f64,
    /// Whole-step endpoint attempts rejected by entropy or species admission.
    /// An interval includes these counts from every accepted physical step.
    /// Internal chemistry, molecular and advective rejections retain their
    /// original owners' separate counters.
    pub rejected_entropy_steps: usize,
}

/// Actual accepted chemistry transfers, independent of the step's entropy basis.
/// Only the original composition owner extracts and merges these receipts after
/// its complete physical step succeeds. No rounded report sums enter the account.
#[derive(Clone, Copy)]
pub(in crate::multicomponent) struct ReactionInventory {
    species_count: usize,
    chemical_change: [Dd; MAX_SPECIES],
    chemical_throughput: [Dd; MAX_SPECIES],
    rejected_entropy_steps: usize,
}

/// The composition owner supplies actual chemical endpoint rows and the complete
/// original/final fields. This account advances no state and owns no integrator.
pub(in crate::multicomponent) struct ReactiveEntropyAccount {
    state: EntropyTracker,
    species_count: usize,
    chemical_change: [Dd; MAX_SPECIES],
    chemical_throughput: [Dd; MAX_SPECIES],
    rejected_entropy_steps: usize,
    // Optional fixed exterior efforts in the ORIGINAL pre-chemistry basis.
    // Length is zero or one; closed accounts incur no additional allocation.
    exterior: Vec<ExteriorEntropyAccount>,
}

impl ReactiveEntropyAccount {
    pub(in crate::multicomponent) fn new(
        initial: &MixtureEuler1d,
        controls: MolecularEntropyControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        Ok(Self {
            state: EntropyTracker::new(initial, controls, cancelled)?,
            species_count: initial.gas().species_count(),
            chemical_change: [Dd::ZERO; MAX_SPECIES],
            chemical_throughput: [Dd::ZERO; MAX_SPECIES],
            rejected_entropy_steps: 0,
            exterior: Vec::new(),
        })
    }

    /// Each complete retry returns to the same immutable original thermal basis.
    /// Rejected endpoint temperatures and chemical transfers cannot become anchors.
    pub(in crate::multicomponent) fn begin_attempt(&mut self) {
        self.reset_endpoint();
        self.chemical_change.fill(Dd::ZERO);
        self.chemical_throughput.fill(Dd::ZERO);
    }

    fn reset_endpoint(&mut self) {
        self.state.current_temperature.copy_from_slice(&self.state.initial_temperature);
        self.state.trial_entropy = Dd::ZERO;
        self.state.trial_mass.fill(Dd::ZERO);
        self.state.trial_relative = 0.0;
    }

    /// Copy only after native whole-step success. Rejected attempts have already
    /// been discarded by begin_attempt; their work counter remains diagnostic.
    pub(in crate::multicomponent) fn accepted_inventory(&self) -> ReactionInventory {
        ReactionInventory {
            species_count: self.species_count,
            chemical_change: self.chemical_change,
            chemical_throughput: self.chemical_throughput,
            rejected_entropy_steps: self.rejected_entropy_steps,
        }
    }

    /// Accumulate accepted chemistry into an interval's immutable original
    /// account. Source-reference contributions and endpoint entropy use that
    /// original basis at completion; per-step entropy allowances are never summed.
    /// Refusal or cancellation leaves all previously merged transfers intact.
    pub(in crate::multicomponent) fn merge_inventory(
        &mut self,
        inventory: &ReactionInventory,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        poll(cancelled)?;
        if inventory.species_count != self.species_count {
            return Err(FlowError::EntropyDomain("reacting entropy source species axis"));
        }
        let rejected_entropy_steps = self.rejected_entropy_steps
            .checked_add(inventory.rejected_entropy_steps)
            .ok_or(FlowError::EntropyDomain("reacting entropy rejection counter overflow"))?;
        let mut change = self.chemical_change;
        let mut throughput = self.chemical_throughput;
        for k in 0..self.species_count {
            poll(cancelled)?;
            change[k] = change[k] + inventory.chemical_change[k];
            throughput[k] = throughput[k] + inventory.chemical_throughput[k];
            scalar(change[k])?;
            scalar(throughput[k])?;
        }
        poll(cancelled)?;
        self.chemical_change = change;
        self.chemical_throughput = throughput;
        self.rejected_entropy_steps = rejected_entropy_steps;
        Ok(())
    }

    /// Observe the actual represented cell update after the original chemistry
    /// owner succeeds. Form differences before rounding their extensive sums.
    /// Absolute per-cell/per-half throughput cannot cancel between reactions.
    pub(in crate::multicomponent) fn observe_chemistry_row(
        &mut self,
        before: &[f64],
        after: &[f64],
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        if before.len() < self.species_count
            || after.len() < self.species_count
            || !positive(volume)
        {
            return Err(FlowError::EntropyDomain("reacting entropy source dimensions or volume"));
        }
        let mut change = self.chemical_change;
        let mut throughput = self.chemical_throughput;
        let volume = Dd::from_f64(volume);
        for k in 0..self.species_count {
            poll(cancelled)?;
            if !before[k].is_finite() || before[k] < 0.0
                || !after[k].is_finite() || after[k] < 0.0
            {
                return Err(FlowError::EntropyDomain("reacting entropy source species state"));
            }
            let increment = volume * (Dd::from_f64(after[k]) - Dd::from_f64(before[k]));
            let absolute = if after[k] >= before[k] {
                increment
            } else {
                Dd::ZERO - increment
            };
            change[k] = change[k] + increment;
            throughput[k] = throughput[k] + absolute;
            scalar(change[k])?;
            scalar(throughput[k])?;
        }
        self.chemical_change = change;
        self.chemical_throughput = throughput;
        Ok(())
    }

    /// Admit the complete original-to-final field after both chemistry halves
    /// and any molecular halves. A single actual final decode also checks the
    /// chemistry interval when the last molecular half could leave that interval.
    pub(in crate::multicomponent) fn check_candidate(
        &mut self,
        initial: &MixtureEuler1d,
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveEntropyReport, FlowError> {
        if candidate.gas() != initial.gas()
            || initial.mesh != candidate.mesh
            || initial.boundaries != candidate.boundaries
        {
            return Err(FlowError::EntropyDomain("reacting entropy candidate owner or dimensions"));
        }
        self.check_candidate_rows(&initial.states, candidate, temperature_bounds, cancelled)
    }

    /// Check an entire interval directly against its existing conserved rollback
    /// buffer. The caller retains the original fixed source, mesh and boundaries.
    /// Reset only endpoint scratch: accepted chemical transfers and their absolute
    /// throughput survive, while no earlier endpoint temperature becomes an anchor.
    pub(in crate::multicomponent) fn check_candidate_rows(
        &mut self,
        initial: &[f64],
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveEntropyReport, FlowError> {
        self.reset_endpoint();
        let check = self.check_candidate_inner(initial, candidate, temperature_bounds, cancelled);
        if matches!(&check, Err(
            FlowError::EntropyDecrease { .. } | FlowError::SpeciesInventoryDrift { .. }
        )) {
            self.rejected_entropy_steps += 1;
        }
        check
    }

    /// Shared original-to-final state reduction. The closed and open accounts
    /// use the identical native source inversion and temperature-domain check.
    fn collect_candidate_rows(
        &mut self,
        initial: &[f64],
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        let gas = candidate.gas();
        let n = self.species_count;
        let width = n + 2;
        if gas.species_count() != n
            || initial.len() != candidate.states.len()
            || candidate.cell_count() != self.state.initial_temperature.len()
        {
            return Err(FlowError::EntropyDomain("reacting entropy candidate owner or dimensions"));
        }
        if temperature_bounds.is_some_and(|[lower, upper]| {
            !positive(lower) || !positive(upper) || lower > upper
        }) {
            return Err(FlowError::EntropyDomain("reacting entropy temperature interval"));
        }
        for (i, (old, new)) in initial.chunks_exact(width)
            .zip(candidate.states.chunks_exact(width)).enumerate()
        {
            poll(cancelled)?;
            let primitive = gas.decode(new, Some(self.state.initial_temperature[i]), cancelled)?;
            if let Some([lower, upper]) = temperature_bounds
                && !(lower..=upper).contains(&primitive.temperature)
            {
                return Err(FlowError::EntropyDomain("final molecular state is outside chemistry interval"));
            }
            self.state.candidate_cell(
                gas, i, old, new, primitive.temperature, candidate.mesh.volumes()[i], cancelled,
            )?;
        }
        Ok(())
    }

    fn check_candidate_inner(
        &mut self,
        initial: &[f64],
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveEntropyReport, FlowError> {
        self.collect_candidate_rows(initial, candidate, temperature_bounds, cancelled)?;
        let gas = candidate.gas();
        let n = self.species_count;
        let relative = self.state.inventory_with_reaction(
            n, true, Some((&self.chemical_change, &self.chemical_throughput)),
        )?;
        let (reference, reference_residual, chemical_reference) =
            self.state.reference_changes(gas, Some(&self.chemical_change), cancelled)?;
        let change = scalar(self.state.trial_entropy + chemical_reference)?;
        if change < -self.state.controls.absolute_tolerance_j_k {
            return Err(FlowError::EntropyDecrease {
                change_j_k: change,
                tolerance_j_k: self.state.controls.absolute_tolerance_j_k,
            });
        }
        // Reduce every report component before the complete composition can
        // publish. Failed source/reference arithmetic retains the original anchor.
        Ok(ReactiveEntropyReport {
            entropy_change_j_k: change,
            source_entropy_change_j_k: scalar(self.state.trial_entropy + reference)?,
            reference_entropy_residual_j_k: scalar(reference_residual)?,
            maximum_relative_species_residual: relative,
            rejected_entropy_steps: self.rejected_entropy_steps,
        })
    }
}
