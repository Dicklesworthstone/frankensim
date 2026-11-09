//! Finite fuel and product feedback on the spatial neutron history.
//!
//! Each bound macroscopic channel is assembled from the original microscopic
//! burnup data: Sigma = mean_parent_population * sigma(T) / cell_volume.
//! Neutron fluence and the original product propagator are iterated from the
//! SAME initial histories until their integrated reaction counts agree within
//! the caller's allowance. No trial iterate is an accepted physical substep.
//! This closes the discrete reaction accounts; it does not make the frozen
//! temperature/group data or within-step flux approximation exact.
//! External neutron/precursor rates and continuous nuclide supply participate
//! through their original affine propagators. The product of interval-mean
//! flux and inventory still approximates their time correlation; simultaneous
//! source startup therefore needs physical timestep refinement.
//! Optional finite nuclide transport uses the original joint depletion owner;
//! it changes the spatial targets and places event heat in the receiving cell.

use std::{collections::BTreeSet, sync::Arc};

use super::{
    DiffusionBudget, DiffusionError, DiffusionReactionKind, DiffusionState, DiffusionStep,
};
use crate::burnup::{BurnupBudget, BurnupError, BurnupInventory, BurnupRate, BurnupStep};
use crate::{Accumulator, FissionError, nonnegative, positive, sum};

mod transport;
pub use transport::DiffusionBurnupTransportStep;
use transport::InventoryOwner;

/// Radioactive emission from the same finite, neutron-reacting inventories.
pub mod decay;

#[cfg(test)]
use super::DiffusionModel;

/// Bind one incident-energy neutron channel to its original nuclear reaction.
/// The diffusion channel supplies the physical cell and incident group.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReactionBinding {
    pub diffusion_reaction: usize,
    pub burnup_reaction: usize,
}

/// Constant external sources over one physical interval. Neutrons are
/// cell-major then group; precursors are cell-major then delayed family.
/// Each nuclide row follows that cell's original burnup inventory order.
/// All rates are extensive expected particles/s or nuclei/s, not densities.
/// Nuclide feed does not imply feed enthalpy or a change in bulk thermal mass.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionBurnupSources {
    pub neutron_rates_s: Vec<f64>,
    pub precursor_rates_s: Vec<f64>,
    pub nuclide_rates_s: Vec<Vec<f64>>,
}

impl DiffusionBurnupSources {
    /// Validate against the original spatial/inventory layout without preparing
    /// any propagator or changing a history. Zero rates still require full shape.
    pub fn validate(&self, state: &DiffusionBurnupState) -> Result<(), DiffusionBurnupError> {
        self.validate_with_checkpoint(state, || Ok(()))
    }

    pub fn validate_with_checkpoint(
        &self,
        state: &DiffusionBurnupState,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), DiffusionBurnupError> {
        checkpoint()?;
        if self.neutron_rates_s.len() != state.neutrons.neutron_counts().len()
            || self.precursor_rates_s.len() != state.neutrons.precursor_counts().len()
            || self.nuclide_rates_s.len() != state.inventories().len()
        {
            return Err(DiffusionBurnupError::Invalid(
                "external source vectors must match original neutron, precursor and cell layouts",
            ));
        }
        for &rate in self.neutron_rates_s.iter().chain(&self.precursor_rates_s) {
            checkpoint()?;
            if !nonnegative(rate) {
                return Err(DiffusionBurnupError::Invalid(
                    "finite nonnegative external particle rates required",
                ));
            }
        }
        for (rates, inventory) in self.nuclide_rates_s.iter().zip(state.inventories()) {
            checkpoint()?;
            if rates.len() != inventory.populations().len() {
                return Err(DiffusionBurnupError::Invalid(
                    "one external rate per original cell nuclide required",
                ));
            }
            for &rate in rates {
                checkpoint()?;
                if !nonnegative(rate) {
                    return Err(DiffusionBurnupError::Invalid(
                        "finite nonnegative external nuclide rates required",
                    ));
                }
            }
        }
        checkpoint()?;
        Ok(())
    }
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
    inventories: InventoryOwner,
    bindings: Arc<[ReactionBinding]>,
    decay_neutrons: Option<Arc<[Arc<decay::BurnupDecayNeutronEmission>]>>,
}

#[derive(Clone, Debug)]
pub struct DiffusionBurnupStep {
    pub neutrons: DiffusionStep,
    /// Original burnup reports in diffusion-cell order.
    pub burnup: Vec<BurnupStep>,
    /// External supply integrated by the original production owner, in the
    /// original cell/nuclide order. Reaction-created products are not included.
    pub produced_nuclei: Vec<Vec<f64>>,
    /// Scalar exports of the supplied counts per cell. Their projection and
    /// summation defects are included in that cell's burnup population roundoff.
    pub produced_nuclei_by_cell: Vec<f64>,
    /// Joint finite-material movement, when explicitly bound to this state.
    /// Nuclear-event reports and external supply remain in the original fields.
    pub transport: Option<DiffusionBurnupTransportStep>,
    /// Internally emitted decay neutrons, cell-major then group. These births
    /// are included in the neutron owner's total injected counts, separately
    /// from externally prescribed sources. None means no emission binding.
    pub decay_neutron_counts: Option<Vec<Vec<f64>>>,
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
            inventories: InventoryOwner::Independent(inventories),
            bindings: bindings.into(),
            decay_neutrons: None,
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
        self.inventories.inventories()
    }
    #[must_use]
    pub fn bindings(&self) -> &[ReactionBinding] {
        &self.bindings
    }

    fn validate_bindings(&self) -> Result<(), DiffusionBurnupError> {
        let model = self.neutrons.model();
        if self.inventories().len() != model.cell_count() {
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
            let nuclear = self.inventories()[neutron.cell]
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
                    DiffusionReactionKind::NeutronProduction { .. },
                    BurnupRate::NeutronProduction { .. }
                        | BurnupRate::TabulatedNeutronProduction { .. }
                ) | (
                    DiffusionReactionKind::Fission { .. },
                    BurnupRate::Fission { .. } | BurnupRate::TabulatedFission { .. }
                )
            );
            if let (
                DiffusionReactionKind::NeutronProduction { neutron_yields },
                BurnupRate::TabulatedNeutronProduction { cross_sections },
            ) = (&neutron.kind, &nuclear.rate)
            {
                let expected = cross_sections
                    .reaction()
                    .emitted_neutrons_per_event()
                    .ok_or(DiffusionBurnupError::Invalid(
                        "neutron-production table must retain a supported multiplicity identity",
                    ))?;
                let expected = f64::from(expected);
                if (sum(neutron_yields) - expected).abs() > 64.0 * f64::EPSILON * expected {
                    return Err(DiffusionBurnupError::Invalid(
                        "authored outgoing yields disagree with the neutron-production table identity",
                    ));
                }
            }
            if !same_kind
                || !used.insert(binding.diffusion_reaction)
                || !coverage.insert((
                    neutron.cell,
                    binding.burnup_reaction,
                    neutron.incident_group,
                ))
            {
                return Err(DiffusionBurnupError::Invalid(
                    "distinct matching fission, capture or nonfission neutron-production bindings required",
                ));
            }
        }
        for (cell, inventory) in self.inventories().iter().enumerate() {
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
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DiffusionBurnupStep), DiffusionBurnupError> {
        self.trial_optional_sources(duration_s, temperatures_k, None, budget, checkpoint)
    }

    /// Simultaneous external particles and continuous nuclide feed. Every
    /// nonlinear iterate starts from the original neutron and nuclide histories,
    /// so repeated trials neither inject the source twice nor reset product ages.
    pub fn trial_with_sources(
        &self,
        duration_s: f64,
        temperatures_k: &[f64],
        sources: &DiffusionBurnupSources,
        budget: DiffusionBurnupBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DiffusionBurnupStep), DiffusionBurnupError> {
        self.trial_optional_sources(
            duration_s,
            temperatures_k,
            Some(sources),
            budget,
            checkpoint,
        )
    }

    fn trial_optional_sources(
        &self,
        duration_s: f64,
        temperatures_k: &[f64],
        sources: Option<&DiffusionBurnupSources>,
        budget: DiffusionBurnupBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DiffusionBurnupStep), DiffusionBurnupError> {
        checkpoint()?;
        budget.validate()?;
        let inventories = self.inventories();
        if !nonnegative(duration_s) || temperatures_k.len() != inventories.len() {
            return Err(DiffusionBurnupError::Invalid(
                "finite duration and one temperature per cell required",
            ));
        }
        self.neutrons
            .model()
            .check_limits(budget.diffusion.limits)?;
        self.admit_transport(budget.burnup, &mut checkpoint)?;
        if let Some(sources) = sources {
            sources.validate_with_checkpoint(self, &mut checkpoint)?;
        }
        if self.bindings.len() > budget.diffusion.limits.max_reactions {
            return Err(FissionError::Budget("coupled neutron reaction bindings").into());
        }
        // Admit current original budgets and all table domains even at zero flux.
        // Processed microscopic coefficients are evaluated once per trial, not
        // again during each nonlinear population iterate.
        let mut spontaneous_rates = Vec::new();
        let mut flux_guesses = Vec::new();
        for (cell, inventory) in inventories.iter().enumerate() {
            checkpoint()?;
            let zero_flux = crate::irradiation::GroupFlux::new(
                Arc::clone(inventory.model().groups()),
                vec![0.0; inventory.model().groups().len()],
            )?;
            let rates = inventory.model().rate_constants(
                &zero_flux,
                temperatures_k[cell],
                budget.burnup.limits,
                &mut checkpoint,
            )?;
            if self.decay_neutrons.is_some() {
                spontaneous_rates.push(rates);
                flux_guesses.push(zero_flux);
            }
        }
        let mut microscopic = Vec::with_capacity(self.bindings.len());
        for binding in self.bindings.iter() {
            checkpoint()?;
            let neutron = &self.neutrons.model().reactions()[binding.diffusion_reaction];
            let nuclear = &inventories[neutron.cell].model().reactions()[binding.burnup_reaction];
            let sigma = match &nuclear.rate {
                BurnupRate::Fission { cross_sections_m2 }
                | BurnupRate::Capture { cross_sections_m2 }
                | BurnupRate::NeutronProduction { cross_sections_m2 } => {
                    cross_sections_m2[neutron.incident_group]
                }
                BurnupRate::TabulatedFission { cross_sections }
                | BurnupRate::TabulatedCapture { cross_sections }
                | BurnupRate::TabulatedNeutronProduction { cross_sections } => cross_sections
                    .cross_section_at_group(temperatures_k[neutron.cell], neutron.incident_group)?,
                BurnupRate::Decay { .. } => {
                    return Err(DiffusionBurnupError::Invalid(
                        "spontaneous decay is not incident neutron removal",
                    ));
                }
            };
            microscopic.push(sigma);
        }
        let mut means = Vec::with_capacity(inventories.len());
        for (cell, inventory) in inventories.iter().enumerate() {
            checkpoint()?;
            let mut mean = inventory.populations().to_vec();
            if let Some(sources) = sources {
                for (parent, (&initial, &rate)) in inventory
                    .populations()
                    .iter()
                    .zip(&sources.nuclide_rates_s[cell])
                    .enumerate()
                {
                    checkpoint()?;
                    mean[parent] = supplied_mean(initial, supplied_count(rate, duration_s)?)?;
                }
            }
            means.push(mean);
        }
        let mut last_defect = 0.0;
        for iteration in 1..=budget.max_iterations {
            checkpoint()?;
            let mut data = self.neutrons.model().data().clone();
            for (binding, &sigma) in self.bindings.iter().zip(&microscopic) {
                checkpoint()?;
                let reaction = &mut data.reactions[binding.diffusion_reaction];
                let parent =
                    inventories[reaction.cell].model().reactions()[binding.burnup_reaction].parent;
                reaction.macroscopic_cross_section_m_inv = macroscopic(
                    means[reaction.cell][parent],
                    sigma,
                    data.cell_volumes_m3[reaction.cell],
                )?;
            }
            let model = if &data == self.neutrons.model().data() {
                Arc::clone(self.neutrons.model())
            } else {
                Arc::new(self.neutrons.model().rebuild_with_data(
                    data,
                    budget.diffusion.limits,
                    &mut checkpoint,
                )?)
            };
            let initial = self.neutrons.rebind_model(model)?;
            let (neutrons, neutron_step, trial, decay_neutron_counts) =
                if self.decay_neutrons.is_some() {
                    let (neutrons, report, trial, births) = self.trial_decay_neutrons(
                        &initial,
                        duration_s,
                        &flux_guesses,
                        &spontaneous_rates,
                        &microscopic,
                        sources,
                        budget,
                        &mut checkpoint,
                    )?;
                    (neutrons, report, trial, Some(births))
                } else {
                    let (neutrons, neutron_step) = if let Some(sources) = sources {
                        initial.trial_with_sources(
                            duration_s,
                            &sources.neutron_rates_s,
                            &sources.precursor_rates_s,
                            budget.diffusion,
                            &mut checkpoint,
                        )?
                    } else {
                        initial.trial(duration_s, budget.diffusion, &mut checkpoint)?
                    };
                    // Every iterate restarts the original complete history.
                    let trial = self.inventories.trial(
                        duration_s,
                        &neutron_step.mean_fluxes,
                        temperatures_k,
                        sources.map(|sources| sources.nuclide_rates_s.as_slice()),
                        budget.burnup,
                        &mut checkpoint,
                    )?;
                    (neutrons, neutron_step, trial, None)
                };
            let mut by_reaction: Vec<Vec<Accumulator>> = inventories
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
            let mut defects = Vec::with_capacity(inventories.len());
            // Both trial rates may initially vanish for a newly born target
            // irradiated by a newly born neutron. Compare with their ACTUAL
            // interval means as well, so zero/zero event agreement cannot
            // falsely accept that inactive initial guess.
            let mean_defects = if self.decay_neutrons.is_some() {
                decay::mean_reaction_defects(
                    self,
                    &neutron_step,
                    &trial.means,
                    &microscopic,
                    &mut checkpoint,
                )?
            } else {
                vec![0.0; inventories.len()]
            };
            for (cell, inventory) in inventories.iter().enumerate() {
                let mut defect = Accumulator::new(mean_defects[cell]);
                for (index, reaction) in inventory.model().reactions().iter().enumerate() {
                    checkpoint()?;
                    if !matches!(reaction.rate, BurnupRate::Decay { .. }) {
                        defect.add(reaction_count_defect(
                            &by_reaction[cell][index],
                            trial.reports[cell].reaction_counts[index],
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
                        inventories: trial.owner,
                        bindings: Arc::clone(&self.bindings),
                        decay_neutrons: self.decay_neutrons.clone(),
                    },
                    DiffusionBurnupStep {
                        neutrons: neutron_step,
                        burnup: trial.reports,
                        produced_nuclei: trial.produced,
                        produced_nuclei_by_cell: trial.produced_by_cell,
                        transport: trial.transport,
                        decay_neutron_counts,
                        iterations: iteration,
                        reaction_defect: last_defect,
                        reaction_defect_by_cell: defects,
                    },
                ));
            }
            if self.decay_neutrons.is_some() {
                for (old, new) in flux_guesses.iter_mut().zip(&neutron_step.mean_fluxes) {
                    let mut values = Vec::with_capacity(new.values_m2_s().len());
                    for (&old, &new) in old.values_m2_s().iter().zip(new.values_m2_s()) {
                        checkpoint()?;
                        let value = if old <= new {
                            (new - old).mul_add(budget.relaxation, old)
                        } else {
                            (old - new).mul_add(1.0 - budget.relaxation, new)
                        };
                        values.push(value);
                    }
                    *old = crate::irradiation::GroupFlux::new(Arc::clone(old.groups()), values)?;
                }
            }
            for (old, new) in means.iter_mut().zip(trial.means) {
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

    /// Atomically publish all spatial particles, supplied nuclei, products and
    /// energy accounts after both numerical owners and channel closure accept.
    pub fn advance_with_sources(
        &mut self,
        duration_s: f64,
        temperatures_k: &[f64],
        sources: &DiffusionBurnupSources,
        budget: DiffusionBurnupBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DiffusionBurnupStep, DiffusionBurnupError> {
        let (next, report) =
            self.trial_with_sources(duration_s, temperatures_k, sources, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

fn supplied_count(rate: f64, duration: f64) -> Result<f64, FissionError> {
    let count = rate * duration;
    if !nonnegative(count) || (rate > 0.0 && duration > 0.0 && count == 0.0) {
        return Err(FissionError::Arithmetic);
    }
    Ok(count)
}

fn supplied_mean(initial: f64, produced: f64) -> Result<f64, FissionError> {
    let mean = initial + produced * 0.5;
    if !nonnegative(mean) || (mean == 0.0 && (initial > 0.0 || produced > 0.0)) {
        return Err(FissionError::Arithmetic);
    }
    Ok(mean)
}

fn export_supplied_count(
    produced: &[f64],
    report: &mut BurnupStep,
    budget: BurnupBudget,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<f64, DiffusionBurnupError> {
    let mut total = Accumulator::new(0.0);
    let mut error = Accumulator::new(report.population_roundoff_nuclei);
    for &count in produced {
        checkpoint()?;
        error.add(total.add(count)?.abs())?;
    }
    let rounded = total.value();
    if !nonnegative(rounded) {
        return Err(FissionError::Arithmetic.into());
    }
    error.add(crate::sum(&[total.hi, total.lo, -rounded]).abs())?;
    report.population_roundoff_nuclei = upper_value(&error)?;
    if report.population_roundoff_nuclei > budget.nuclei_tolerance {
        return Err(BurnupError::Budget("coupled supplied-nuclei scalar export roundoff").into());
    }
    Ok(rounded)
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

#[cfg(test)]
mod source_tests;

#[cfg(test)]
mod transport_tests;
