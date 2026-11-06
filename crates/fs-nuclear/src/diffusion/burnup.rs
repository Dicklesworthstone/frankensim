//! Finite fuel and product feedback on the spatial neutron history.
//!
//! Each bound macroscopic channel is assembled from the original microscopic
//! burnup data: Sigma = mean_parent_population * sigma(T) / cell_volume.
//! Neutron fluence and the original product propagator are iterated from the
//! SAME initial histories until their integrated reaction counts agree within
//! the caller's allowance. No trial iterate is an accepted physical substep.
//! This closes the discrete reaction accounts; it does not make the frozen
//! temperature/group data or within-step flux approximation exact.

use std::{collections::BTreeSet, sync::Arc};

use super::{
    DiffusionBudget, DiffusionError, DiffusionModel, DiffusionReactionKind, DiffusionState,
    DiffusionStep,
};
use crate::burnup::{BurnupBudget, BurnupError, BurnupInventory, BurnupRate, BurnupStep};
use crate::{Accumulator, FissionError, nonnegative, positive};

/// Bind one incident-energy neutron channel to its original nuclear reaction.
/// The diffusion channel supplies the physical cell and incident group.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReactionBinding {
    pub diffusion_reaction: usize,
    pub burnup_reaction: usize,
}

/// Component work/accounts and a separate nonlinear reaction-consistency limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionBurnupBudget {
    pub diffusion: DiffusionBudget,
    pub burnup: BurnupBudget,
    /// Bounded Picard iterations, each starting from the accepted physical state.
    pub max_iterations: usize,
    /// Sum of absolute neutron/burnup channel-count mismatches [events/step].
    pub reaction_tolerance: f64,
    /// Population-iterate relaxation in (0, 1]; one avoids unnecessary damping.
    pub relaxation: f64,
}

impl DiffusionBurnupBudget {
    pub(crate) fn validate(self) -> Result<(), DiffusionBurnupError> {
        if self.max_iterations == 0
            || self.max_iterations > 256
            || !nonnegative(self.reaction_tolerance)
            || !positive(self.relaxation)
            || self.relaxation > 1.0
        {
            return Err(DiffusionBurnupError::Invalid(
                "bounded iterations, nonnegative allowance and relaxation in (0,1] required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DiffusionBurnupError {
    Invalid(&'static str),
    Diffusion(DiffusionError),
    Burnup(BurnupError),
    Fission(FissionError),
    /// No candidate is published when the discrete channel accounts disagree.
    Nonconvergence {
        iterations: usize,
        reaction_defect: f64,
    },
}
impl From<DiffusionError> for DiffusionBurnupError {
    fn from(value: DiffusionError) -> Self {
        Self::Diffusion(value)
    }
}
impl From<BurnupError> for DiffusionBurnupError {
    fn from(value: BurnupError) -> Self {
        Self::Burnup(value)
    }
}
impl From<FissionError> for DiffusionBurnupError {
    fn from(value: FissionError) -> Self {
        Self::Fission(value)
    }
}
impl core::fmt::Display for DiffusionBurnupError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial neutron/burnup coupling refused: {self:?}")
    }
}
impl std::error::Error for DiffusionBurnupError {}

/// Complete spatial neutron, precursor, fuel and radioactive-product history.
#[derive(Clone, Debug)]
pub struct DiffusionBurnupState {
    neutrons: DiffusionState,
    inventories: Vec<BurnupInventory>,
    bindings: Arc<[ReactionBinding]>,
}

#[derive(Clone, Debug)]
pub struct DiffusionBurnupStep {
    pub neutrons: DiffusionStep,
    /// Original burnup reports in diffusion-cell order.
    pub burnup: Vec<BurnupStep>,
    pub iterations: usize,
    pub reaction_defect: f64,
    pub reaction_defect_by_cell: Vec<f64>,
}

impl DiffusionBurnupState {
    /// Every neutron-induced burnup channel must be bound exactly once per
    /// incident group. Unbound diffusion channels are explicitly fixed
    /// background interactions, whose target inventories are not tracked here.
    pub fn new(
        neutrons: DiffusionState,
        inventories: Vec<BurnupInventory>,
        bindings: Vec<ReactionBinding>,
    ) -> Result<Self, DiffusionBurnupError> {
        let state = Self {
            neutrons,
            inventories,
            bindings: bindings.into(),
        };
        state.validate_bindings()?;
        Ok(state)
    }

    #[must_use]
    pub fn neutrons(&self) -> &DiffusionState {
        &self.neutrons
    }
    #[must_use]
    pub fn inventories(&self) -> &[BurnupInventory] {
        &self.inventories
    }
    #[must_use]
    pub fn bindings(&self) -> &[ReactionBinding] {
        &self.bindings
    }

    fn validate_bindings(&self) -> Result<(), DiffusionBurnupError> {
        let model = self.neutrons.model();
        if self.inventories.len() != model.cell_count() {
            return Err(DiffusionBurnupError::Invalid(
                "one original burnup inventory per diffusion cell required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut coverage = BTreeSet::new();
        for binding in self.bindings.iter() {
            let neutron = model
                .reactions()
                .get(binding.diffusion_reaction)
                .ok_or(DiffusionBurnupError::Invalid("unknown diffusion reaction"))?;
            let nuclear = self.inventories[neutron.cell]
                .model()
                .reactions()
                .get(binding.burnup_reaction)
                .ok_or(DiffusionBurnupError::Invalid("unknown burnup reaction"))?;
            let same_kind = matches!(
                (&neutron.kind, &nuclear.rate),
                (
                    DiffusionReactionKind::Capture,
                    BurnupRate::Capture { .. } | BurnupRate::TabulatedCapture { .. }
                ) | (
                    DiffusionReactionKind::Fission { .. },
                    BurnupRate::Fission { .. } | BurnupRate::TabulatedFission { .. }
                )
            );
            if !same_kind
                || !used.insert(binding.diffusion_reaction)
                || !coverage.insert((
                    neutron.cell,
                    binding.burnup_reaction,
                    neutron.incident_group,
                ))
            {
                return Err(DiffusionBurnupError::Invalid(
                    "distinct matching fission/capture channel bindings required",
                ));
            }
        }
        for (cell, inventory) in self.inventories.iter().enumerate() {
            if inventory.model().groups().as_ref() != model.groups().as_ref() {
                return Err(DiffusionBurnupError::Invalid(
                    "burnup and diffusion require the same named neutron grid",
                ));
            }
            for (reaction, channel) in inventory.model().reactions().iter().enumerate() {
                if matches!(channel.rate, BurnupRate::Decay { .. }) {
                    continue;
                }
                for group in 0..model.groups().len() {
                    if !coverage.contains(&(cell, reaction, group)) {
                        return Err(DiffusionBurnupError::Invalid(
                            "every burnup neutron channel and group must affect the neutron balance",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Actual temperatures belong to the original material cells, in model order.
    pub fn trial(
        &self,
        duration_s: f64,
        temperatures_k: &[f64],
        budget: DiffusionBurnupBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DiffusionBurnupStep), DiffusionBurnupError> {
        checkpoint()?;
        budget.validate()?;
        if !nonnegative(duration_s) || temperatures_k.len() != self.inventories.len() {
            return Err(DiffusionBurnupError::Invalid(
                "finite duration and one temperature per cell required",
            ));
        }
        if self.bindings.len() > budget.diffusion.limits.max_reactions {
            return Err(FissionError::Budget("coupled neutron reaction bindings").into());
        }
        // Admit current original budgets and all table domains even at zero flux.
        // Processed microscopic coefficients are evaluated once per trial, not
        // again during each nonlinear population iterate.
        for (cell, inventory) in self.inventories.iter().enumerate() {
            checkpoint()?;
            let zero_flux = crate::irradiation::GroupFlux::new(
                Arc::clone(inventory.model().groups()),
                vec![0.0; inventory.model().groups().len()],
            )?;
            inventory.model().rate_constants(
                &zero_flux,
                temperatures_k[cell],
                budget.burnup.limits,
                &mut checkpoint,
            )?;
        }
        let mut microscopic = Vec::with_capacity(self.bindings.len());
        for binding in self.bindings.iter() {
            checkpoint()?;
            let neutron = &self.neutrons.model().reactions()[binding.diffusion_reaction];
            let nuclear =
                &self.inventories[neutron.cell].model().reactions()[binding.burnup_reaction];
            let sigma = match &nuclear.rate {
                BurnupRate::Fission { cross_sections_m2 }
                | BurnupRate::Capture { cross_sections_m2 } => {
                    cross_sections_m2[neutron.incident_group]
                }
                BurnupRate::TabulatedFission { cross_sections }
                | BurnupRate::TabulatedCapture { cross_sections } => cross_sections
                    .cross_section_at_group(temperatures_k[neutron.cell], neutron.incident_group)?,
                BurnupRate::Decay { .. } => {
                    return Err(DiffusionBurnupError::Invalid(
                        "spontaneous decay is not incident neutron removal",
                    ));
                }
            };
            microscopic.push(sigma);
        }
        let mut means: Vec<Vec<f64>> = self
            .inventories
            .iter()
            .map(|inventory| inventory.populations().to_vec())
            .collect();
        let mut last_defect = 0.0;
        for iteration in 1..=budget.max_iterations {
            checkpoint()?;
            let mut data = self.neutrons.model().data().clone();
            for (binding, &sigma) in self.bindings.iter().zip(&microscopic) {
                checkpoint()?;
                let reaction = &mut data.reactions[binding.diffusion_reaction];
                let parent = self.inventories[reaction.cell].model().reactions()
                    [binding.burnup_reaction]
                    .parent;
                reaction.macroscopic_cross_section_m_inv = macroscopic(
                    means[reaction.cell][parent],
                    sigma,
                    data.cell_volumes_m3[reaction.cell],
                )?;
            }
            let model = if &data == self.neutrons.model().data() {
                Arc::clone(self.neutrons.model())
            } else {
                Arc::new(DiffusionModel::new_with_checkpoint(
                    data,
                    budget.diffusion.limits,
                    &mut checkpoint,
                )?)
            };
            let initial = self.neutrons.rebind_model(model)?;
            let (neutrons, neutron_step) =
                initial.trial(duration_s, budget.diffusion, &mut checkpoint)?;
            let mut candidates = Vec::with_capacity(self.inventories.len());
            let mut reports = Vec::with_capacity(self.inventories.len());
            let mut next_means = Vec::with_capacity(self.inventories.len());
            for (cell, inventory) in self.inventories.iter().enumerate() {
                checkpoint()?;
                let flux = &neutron_step.mean_fluxes[cell];
                let (candidate, report) = inventory.trial(
                    duration_s,
                    flux,
                    temperatures_k[cell],
                    budget.burnup,
                    &mut checkpoint,
                )?;
                let rates = inventory.model().rate_constants(
                    flux,
                    temperatures_k[cell],
                    budget.burnup.limits,
                    &mut checkpoint,
                )?;
                let mut removal = vec![Accumulator::new(0.0); inventory.populations().len()];
                let mut hazards = vec![Accumulator::new(0.0); inventory.populations().len()];
                for ((reaction, &rate), &count) in inventory
                    .model()
                    .reactions()
                    .iter()
                    .zip(&rates)
                    .zip(&report.reaction_counts)
                {
                    checkpoint()?;
                    removal[reaction.parent].add(count)?;
                    hazards[reaction.parent].add(rate)?;
                }
                let mut mean = inventory.populations().to_vec();
                if duration_s > 0.0 {
                    for parent in 0..mean.len() {
                        if hazards[parent].value() > 0.0 {
                            mean[parent] = mean_population(
                                removal[parent].value(),
                                hazards[parent].value(),
                                duration_s,
                            )?;
                        }
                    }
                }
                next_means.push(mean);
                candidates.push(candidate);
                reports.push(report);
            }
            let mut by_reaction: Vec<Vec<Accumulator>> = self
                .inventories
                .iter()
                .map(|inventory| vec![Accumulator::new(0.0); inventory.model().reactions().len()])
                .collect();
            for binding in self.bindings.iter() {
                checkpoint()?;
                let cell = self.neutrons.model().reactions()[binding.diffusion_reaction].cell;
                by_reaction[cell][binding.burnup_reaction]
                    .add(neutron_step.reaction_counts[binding.diffusion_reaction])?;
            }
            let mut total = Accumulator::new(0.0);
            let mut defects = Vec::with_capacity(self.inventories.len());
            for (cell, inventory) in self.inventories.iter().enumerate() {
                let mut defect = Accumulator::new(0.0);
                for (index, reaction) in inventory.model().reactions().iter().enumerate() {
                    checkpoint()?;
                    if !matches!(reaction.rate, BurnupRate::Decay { .. }) {
                        defect.add(reaction_count_defect(
                            &by_reaction[cell][index],
                            reports[cell].reaction_counts[index],
                        )?)?;
                    }
                }
                let defect = upper_value(&defect)?;
                defects.push(defect);
                total.add(defect)?;
            }
            last_defect = upper_value(&total)?;
            if last_defect <= budget.reaction_tolerance {
                checkpoint()?;
                return Ok((
                    Self {
                        neutrons,
                        inventories: candidates,
                        bindings: Arc::clone(&self.bindings),
                    },
                    DiffusionBurnupStep {
                        neutrons: neutron_step,
                        burnup: reports,
                        iterations: iteration,
                        reaction_defect: last_defect,
                        reaction_defect_by_cell: defects,
                    },
                ));
            }
            for (old, new) in means.iter_mut().zip(next_means) {
                checkpoint()?;
                for (old, new) in old.iter_mut().zip(new) {
                    *old = if *old <= new {
                        (new - *old).mul_add(budget.relaxation, *old)
                    } else {
                        (*old - new).mul_add(1.0 - budget.relaxation, new)
                    };
                    if !nonnegative(*old) {
                        return Err(FissionError::Arithmetic.into());
                    }
                }
            }
        }
        Err(DiffusionBurnupError::Nonconvergence {
            iterations: budget.max_iterations,
            reaction_defect: last_defect,
        })
    }

    pub fn advance(
        &mut self,
        duration_s: f64,
        temperatures_k: &[f64],
        budget: DiffusionBurnupBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DiffusionBurnupStep, DiffusionBurnupError> {
        let (next, report) = self.trial(duration_s, temperatures_k, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

fn macroscopic(population: f64, sigma: f64, volume: f64) -> Result<f64, FissionError> {
    if population == 0.0 || sigma == 0.0 {
        return Ok(0.0);
    }
    for value in [
        population / volume * sigma,
        population * sigma / volume,
        sigma / volume * population,
    ] {
        if positive(value) {
            return Ok(value);
        }
    }
    Err(FissionError::Arithmetic)
}

fn reaction_count_defect(neutrons: &Accumulator, burnup: f64) -> Result<f64, FissionError> {
    // Keep the retained low part through the subtraction: a small group's
    // events must not vanish beside another group's large count.
    let value = crate::sum(&[neutrons.hi, neutrons.lo, -burnup]).abs();
    if nonnegative(value) {
        Ok(value)
    } else {
        Err(FissionError::Arithmetic)
    }
}

fn upper_value(account: &Accumulator) -> Result<f64, FissionError> {
    let rounded = account.value();
    let remainder = crate::sum(&[account.hi, -rounded, account.lo]);
    let value = if remainder > 0.0 {
        rounded.next_up()
    } else {
        rounded
    };
    if nonnegative(value) {
        Ok(value)
    } else {
        Err(FissionError::Arithmetic)
    }
}

fn mean_population(count: f64, rate: f64, duration: f64) -> Result<f64, FissionError> {
    if count == 0.0 {
        return Ok(0.0);
    }
    let exposure = rate * duration;
    if positive(exposure) {
        let value = count / exposure;
        if positive(value) {
            return Ok(value);
        }
    }
    for value in [count / rate / duration, count / duration / rate] {
        if positive(value) {
            return Ok(value);
        }
    }
    Err(FissionError::Arithmetic)
}

#[cfg(test)]
mod tests;
