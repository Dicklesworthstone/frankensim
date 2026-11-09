//! Radioactive neutrons from the same finite stocks that undergo irradiation.
//!
//! Each Picard iterate freezes the admitted neutron-induced hazards and spatial
//! coefficients, then propagates the complete neutron/nuclide block once. A
//! source edge uses the parent's spontaneous-channel fractions, never its total
//! removal hazard. Capture therefore consumes the emitter without producing a
//! fictitious radioactive neutron. Daughter birth, continuous feed and emission
//! retain their ages in the joint positive propagator.
//!
//! The original reaction-consistency iteration also checks the actual mean
//! target and integrated flux, including initially empty targets/fields. This
//! closes the frozen-mean equations, not the nonlinear time correlation; physical
//! timestep refinement remains necessary. Sources are independent stationary
//! cell inventories. Material movement must not be silently approximated here.

use std::sync::Arc;

use fs_la::decay::{
    MAX_DECAY_EDGES, MAX_DECAY_NODES,
    reaction::{ReactionAdvance, ReactionEdge, ReactionPlan, source::ReactionSourcePlan},
};

use super::super::{
    Accounts, DiffusionEvolution, SOURCE_NEUTRONS, account_add, checked, owner_call,
    source::{PreparedPlan, SourceRates, injected_counts},
};
use super::transport::{InventoryOwner, InventoryTrial};
use super::{
    DiffusionBurnupBudget, DiffusionBurnupError, DiffusionBurnupSources, DiffusionBurnupState,
    DiffusionState, DiffusionStep, export_supplied_count, mean_population, upper_value,
};
use crate::burnup::{BurnupError, BurnupModel, BurnupRate};
use crate::irradiation::GroupFlux;
use crate::{Accumulator, FissionError, nonnegative, sum};

/// Parent-wide expected radioactive neutron yields on an exact burnup model.
/// Neutron-induced reactions and spontaneous decays share that model's stock.
#[derive(Clone, Debug)]
pub struct BurnupDecayNeutronEmission {
    data_id: String,
    model: Arc<BurnupModel>,
    neutron_yields: Vec<Vec<f64>>,
}

impl BurnupDecayNeutronEmission {
    pub fn new(
        data_id: String,
        model: Arc<BurnupModel>,
        neutron_yields: Vec<Vec<f64>>,
    ) -> Result<Self, DiffusionBurnupError> {
        if data_id.trim().is_empty() || neutron_yields.len() != model.nuclides().len() {
            return Err(DiffusionBurnupError::Invalid(
                "identified decay-neutron yields in original burnup nuclide order required",
            ));
        }
        for (nuclide, row) in model.nuclides().iter().zip(&neutron_yields) {
            if row.len() != model.groups().len()
                || !row.iter().copied().all(nonnegative)
                || !sum(row).is_finite()
                || (nuclide.law.decay_constant().value() == 0.0
                    && row.iter().any(|&value| value != 0.0))
            {
                return Err(DiffusionBurnupError::Invalid(
                    "finite nonnegative group yields and zero stable-nuclide emission required",
                ));
            }
        }
        Ok(Self {
            data_id,
            model,
            neutron_yields,
        })
    }

    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }
    pub fn neutron_yields(&self) -> &[Vec<f64>] {
        &self.neutron_yields
    }
}

impl DiffusionBurnupState {
    /// Enable radioactive birth on the original stocks without resetting any
    /// inventory, resident neutron, event count or accumulated energy account.
    /// Every cell supplies its exact original model, including explicit zero
    /// rows. An inventory-transport state is refused instead of freezing motion.
    pub fn with_decay_neutrons(
        mut self,
        emissions: Vec<Arc<BurnupDecayNeutronEmission>>,
    ) -> Result<Self, DiffusionBurnupError> {
        if self.transport().is_some() || emissions.len() != self.inventories().len() {
            return Err(DiffusionBurnupError::Invalid(
                "decay-neutron coupling requires all original independent cell inventories",
            ));
        }
        let mut nodes = self.neutrons.model.nodes;
        for (source, inventory) in emissions.iter().zip(self.inventories()) {
            if !Arc::ptr_eq(source.model(), inventory.model()) {
                return Err(DiffusionBurnupError::Invalid(
                    "decay-neutron source must retain its exact original burnup model",
                ));
            }
            nodes = nodes
                .checked_add(inventory.populations().len())
                .ok_or(BurnupError::Budget("joint neutron/burnup nodes"))?;
            if nodes > MAX_DECAY_NODES {
                return Err(BurnupError::Budget("joint neutron/burnup nodes").into());
            }
        }
        self.decay_neutrons = Some(emissions.into());
        Ok(self)
    }

    pub fn decay_neutrons(&self) -> Option<&[Arc<BurnupDecayNeutronEmission>]> {
        self.decay_neutrons.as_deref()
    }

    pub(super) fn trial_decay_neutrons(
        &self,
        neutron_initial: &DiffusionState,
        duration_s: f64,
        fluxes: &[GroupFlux],
        spontaneous_rates: &[Vec<f64>],
        microscopic: &[f64],
        sources: Option<&DiffusionBurnupSources>,
        budget: DiffusionBurnupBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(DiffusionState, DiffusionStep, InventoryTrial, Vec<Vec<f64>>), DiffusionBurnupError>
    {
        // Keep this joint solver monomorphic across the public trial adapters.
        // The local closure also satisfies the original owners' sized callback
        // interfaces without duplicating their numerical implementations.
        let mut poll = || checkpoint();
        let checkpoint = &mut poll;
        checkpoint()?;
        budget.burnup.validate()?;
        if !nonnegative(budget.diffusion.population_tolerance)
            || budget.burnup.limits.propagation != budget.diffusion.limits.propagation
        {
            return Err(DiffusionBurnupError::Invalid(
                "finite neutron allowance and identical joint propagation limits required",
            ));
        }
        let emissions = self.decay_neutrons().ok_or(DiffusionBurnupError::Invalid(
            "joint radioactive propagation requires an emission binding",
        ))?;
        let inventories = self.inventories();
        let neutron_nodes = neutron_initial.model.nodes;
        let groups = neutron_initial.model.groups().len();
        let limits = budget.burnup.limits.propagation;
        let nodes = inventories
            .iter()
            .try_fold(neutron_nodes, |nodes, inventory| {
                nodes.checked_add(inventory.populations().len())
            })
            .ok_or(BurnupError::Budget("joint neutron/burnup nodes"))?;
        if nodes > limits.maximum_nodes || nodes > MAX_DECAY_NODES {
            return Err(BurnupError::Budget("joint neutron/burnup nodes").into());
        }

        // Microscopic tables were evaluated once at the actual temperature by
        // the original outer iteration. Only their flux contraction changes.
        let mut rate_accounts: Vec<Vec<Accumulator>> = spontaneous_rates
            .iter()
            .map(|row| row.iter().copied().map(Accumulator::new).collect())
            .collect();
        for (binding, &sigma) in self.bindings.iter().zip(microscopic) {
            checkpoint()?;
            let reaction = &neutron_initial.model.reactions()[binding.diffusion_reaction];
            let rate = product(
                sigma,
                fluxes[reaction.cell].values_m2_s()[reaction.incident_group],
            )?;
            rate_accounts[reaction.cell][binding.burnup_reaction].add(rate)?;
        }
        let mut rates = neutron_initial.model.removal_rates.clone();
        let mut edges = neutron_initial.model.edges.clone();
        let mut channels = Vec::with_capacity(inventories.len());
        for (cell, inventory) in inventories.iter().enumerate() {
            checkpoint()?;
            let offset = rates.len();
            let reaction_rates = rate_accounts[cell]
                .iter()
                .map(|rate| rate.value())
                .collect::<Vec<_>>();
            let compiled = inventory
                .model()
                .prepare_channels(&reaction_rates, &mut *checkpoint)?;
            rates.extend_from_slice(&compiled.removal_rates);
            for edge in &compiled.edges {
                push_edge(
                    &mut edges,
                    ReactionEdge {
                        parent: offset + edge.parent,
                        daughter: offset + edge.daughter,
                        yield_per_removal: edge.yield_per_removal,
                    },
                    limits.maximum_edges,
                )?;
            }
            let mut fractions = vec![Accumulator::new(0.0); inventory.populations().len()];
            for (reaction, &fraction) in inventory
                .model()
                .reactions()
                .iter()
                .zip(&compiled.fractions)
            {
                checkpoint()?;
                if matches!(reaction.rate, BurnupRate::Decay { .. }) {
                    fractions[reaction.parent].add(fraction)?;
                }
            }
            for (parent, row) in emissions[cell].neutron_yields().iter().enumerate() {
                checkpoint()?;
                for (group, &yield_) in row.iter().enumerate() {
                    let coefficient = product(fractions[parent].value(), yield_)?;
                    if coefficient > 0.0 {
                        push_edge(
                            &mut edges,
                            ReactionEdge {
                                parent: offset + parent,
                                daughter: cell * groups + group,
                                yield_per_removal: coefficient,
                            },
                            limits.maximum_edges,
                        )?;
                    }
                }
            }
            channels.push(compiled);
        }

        let mut accounts = Accounts::new(neutron_nodes)?;
        let mut injected = injected_counts(
            sources.map(|source| SourceRates {
                neutrons: &source.neutron_rates_s,
                precursors: &source.precursor_rates_s,
            }),
            neutron_initial.model.neutron_nodes,
            neutron_nodes,
            duration_s,
            &mut accounts,
            &mut *checkpoint,
        )?;
        let mut supplied = injected.clone();
        let mut initial = Vec::with_capacity(nodes);
        initial.extend_from_slice(neutron_initial.neutron_counts());
        initial.extend_from_slice(neutron_initial.precursor_counts());
        let mut produced = Vec::with_capacity(inventories.len());
        let mut input_errors = vec![Accumulator::new(0.0); inventories.len()];
        for (cell, inventory) in inventories.iter().enumerate() {
            let mut row = Vec::with_capacity(inventory.populations().len());
            for (parent, &population) in inventory.populations().iter().enumerate() {
                checkpoint()?;
                initial.push(population);
                let rate = sources.map_or(0.0, |source| source.nuclide_rates_s[cell][parent]);
                let count = product(rate, duration_s)?;
                input_errors[cell].add(rate.mul_add(duration_s, -count).abs())?;
                row.push(count);
                supplied.push(count);
            }
            produced.push(row);
        }
        let plan = if supplied.iter().any(|&count| count > 0.0) {
            PreparedPlan::Source(owner_call(checkpoint, |cancelled| {
                ReactionSourcePlan::prepare_with_population_integrals(
                    &rates, &edges, duration_s, limits, cancelled,
                )
            })?)
        } else {
            PreparedPlan::Closed(owner_call(checkpoint, |cancelled| {
                ReactionPlan::prepare_with_population_integrals(
                    &rates, &edges, duration_s, limits, cancelled,
                )
            })?)
        };
        let evolved = plan.apply(&initial, &supplied, checkpoint)?;
        let residence = plan.integrals(&initial, &supplied, checkpoint)?;
        let mut candidates = Vec::with_capacity(inventories.len());
        let mut reports = Vec::with_capacity(inventories.len());
        let mut means = Vec::with_capacity(inventories.len());
        let mut produced_by_cell = Vec::with_capacity(inventories.len());
        let mut births = Vec::with_capacity(inventories.len());
        let mut offset = neutron_nodes;
        let mut nuclei_error = Accumulator::new(0.0);
        let mut energy_error = Accumulator::new(0.0);
        for (cell, inventory) in inventories.iter().enumerate() {
            checkpoint()?;
            let n = inventory.populations().len();
            let (candidate, mut report) =
                inventory.account_evolution(
                    ReactionAdvance {
                        populations: evolved.populations[offset..offset + n].to_vec(),
                        removals: evolved.removals[offset..offset + n].to_vec(),
                        balance_residuals: evolved.balance_residuals[offset..offset + n].to_vec(),
                        population_balance_residual: sum(
                            &evolved.balance_residuals[offset..offset + n]
                        ),
                    },
                    duration_s,
                    &channels[cell].fractions,
                    Some(&produced[cell]),
                    upper_value(&input_errors[cell])?,
                    budget.burnup,
                    &mut *checkpoint,
                )?;
            produced_by_cell.push(export_supplied_count(
                &produced[cell],
                &mut report,
                budget.burnup,
                checkpoint,
            )?);
            nuclei_error.add(report.population_roundoff_nuclei)?;
            energy_error.add(report.energy_roundoff_j)?;
            let mut decays = vec![Accumulator::new(0.0); n];
            let mut decay_count_errors = vec![Accumulator::new(0.0); n];
            for (reaction, &count) in inventory
                .model()
                .reactions()
                .iter()
                .zip(&report.reaction_counts)
            {
                if matches!(reaction.rate, BurnupRate::Decay { .. }) {
                    decay_count_errors[reaction.parent]
                        .add(decays[reaction.parent].add(count)?.abs())?;
                }
            }
            let mut born = vec![Accumulator::new(0.0); groups];
            for ((&decays, mut count_error), row) in decays
                .iter()
                .zip(decay_count_errors)
                .zip(emissions[cell].neutron_yields())
            {
                checkpoint()?;
                let count = checked(decays.value())?;
                count_error.add(checked(sum(&[decays.hi, decays.lo, -count]))?.abs())?;
                let count_error = upper_value(&count_error)?;
                for (group, &yield_) in row.iter().enumerate() {
                    // Parent-count defects are nuclei; convert them with this
                    // group's yield before charging the neutron allowance.
                    let source_error = accounts.multiply(count_error, yield_)?;
                    accounts.error.add(source_error)?;
                    let value = accounts.multiply(count, yield_)?;
                    account_add(&mut born[group], value, &mut accounts.error)?;
                }
            }
            let mut row = Vec::with_capacity(groups);
            for (group, count) in born.into_iter().enumerate() {
                let value = checked(count.value())?;
                accounts
                    .error
                    .add(checked(sum(&[count.hi, count.lo, -value]))?.abs())?;
                let index = cell * groups + group;
                let mut total = Accumulator::new(injected[index]);
                account_add(&mut total, value, &mut accounts.error)?;
                injected[index] = checked(total.value())?;
                accounts
                    .error
                    .add(checked(sum(&[total.hi, total.lo, -injected[index]]))?.abs())?;
                accounts.total(SOURCE_NEUTRONS, value)?;
                row.push(value);
            }
            births.push(row);
            let mut mean = Vec::with_capacity(n);
            for (&population, &integral) in inventory
                .populations()
                .iter()
                .zip(&residence[offset..offset + n])
            {
                checkpoint()?;
                mean.push(if duration_s == 0.0 {
                    population
                } else {
                    mean_population(integral, 1.0, duration_s)?
                });
            }
            means.push(mean);
            candidates.push(candidate);
            reports.push(report);
            offset += n;
        }
        if upper_value(&nuclei_error)? > budget.burnup.nuclei_tolerance
            || upper_value(&energy_error)? > budget.burnup.energy_tolerance_j
        {
            return Err(BurnupError::Budget("aggregate joint burnup accounts").into());
        }
        let (neutrons, report) = neutron_initial.account_evolution(
            duration_s,
            DiffusionEvolution {
                initial: initial[..neutron_nodes].to_vec(),
                final_counts: evolved.populations[..neutron_nodes].to_vec(),
                removals: evolved.removals[..neutron_nodes].to_vec(),
                integrals: residence[..neutron_nodes].to_vec(),
                injected,
                accounts,
            },
            budget.diffusion,
            None,
            &mut *checkpoint,
        )?;
        checkpoint()?;
        Ok((
            neutrons,
            report,
            InventoryTrial {
                owner: InventoryOwner::Independent(candidates),
                reports,
                means,
                produced,
                produced_by_cell,
                transport: None,
            },
            births,
        ))
    }
}

// An active neutron channel must agree with the interval means generated by
// this same joint solve, even when both old guesses were zero. This residual
// supplements the original neutron-versus-nuclide event consistency check.
pub(super) fn mean_reaction_defects(
    state: &DiffusionBurnupState,
    neutrons: &DiffusionStep,
    means: &[Vec<f64>],
    microscopic: &[f64],
    checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
) -> Result<Vec<f64>, DiffusionBurnupError> {
    let groups = state.neutrons.model.groups().len();
    let mut errors = vec![Accumulator::new(0.0); means.len()];
    for (binding, &sigma) in state.bindings.iter().zip(microscopic) {
        checkpoint()?;
        let channel = &state.neutrons.model.reactions()[binding.diffusion_reaction];
        let parent =
            state.inventories()[channel.cell].model().reactions()[binding.burnup_reaction].parent;
        let amount = means[channel.cell][parent];
        let fluence =
            neutrons.integrated_fluence_m2[channel.cell * groups + channel.incident_group];
        let expected = product3(amount, sigma, fluence)?;
        errors[channel.cell].add(
            checked(sum(&[
                expected,
                -neutrons.reaction_counts[binding.diffusion_reaction],
            ]))?
            .abs(),
        )?;
    }
    errors
        .iter()
        .map(|error| upper_value(error).map_err(Into::into))
        .collect()
}

fn product(a: f64, b: f64) -> Result<f64, FissionError> {
    let value = a * b;
    if nonnegative(value) && (value > 0.0 || a == 0.0 || b == 0.0) {
        Ok(value)
    } else {
        Err(FissionError::Arithmetic)
    }
}
fn product3(a: f64, b: f64, c: f64) -> Result<f64, FissionError> {
    if a == 0.0 || b == 0.0 || c == 0.0 {
        return Ok(0.0);
    }
    for value in [(a * b) * c, (a * c) * b, (b * c) * a] {
        if nonnegative(value) && value > 0.0 {
            return Ok(value);
        }
    }
    Err(FissionError::Arithmetic)
}
fn push_edge(
    edges: &mut Vec<ReactionEdge>,
    edge: ReactionEdge,
    maximum: usize,
) -> Result<(), DiffusionBurnupError> {
    if edges.len() >= maximum || edges.len() >= MAX_DECAY_EDGES {
        return Err(BurnupError::Budget("joint neutron/burnup edges").into());
    }
    edges.push(edge);
    Ok(())
}

#[cfg(test)]
mod tests;
