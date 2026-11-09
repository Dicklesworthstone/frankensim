//! Radioactive neutrons and delayed families from the finite irradiated stocks.
//!
//! Each Picard iterate freezes the admitted neutron-induced hazards and spatial
//! coefficients, then propagates the complete neutron/nuclide block once. A
//! source edge uses the parent's spontaneous-channel fractions, never its total
//! removal hazard. Capture therefore consumes the emitter without producing a
//! fictitious radioactive neutron or precursor. Daughter birth, continuous
//! feed and emission retain their ages in the joint positive propagator. Native
//! delayed families decay and move through their existing neutron operator;
//! creating a precursor never creates its eventual neutron prematurely.
//!
//! The original reaction-consistency iteration also checks the actual mean
//! target and integrated flux, including initially empty targets/fields. This
//! closes the frozen-mean equations, not the nonlinear time correlation; physical
//! timestep refinement remains necessary. Independent inventories and the
//! original finite material-transport owner share this same block. A moved
//! nucleus emits in the compartment where it actually decays; exterior outflow
//! removes its future retained emission without creating a neutron or heat.

use std::sync::Arc;

use fs_la::decay::{
    MAX_DECAY_EDGES, MAX_DECAY_NODES,
    reaction::{ReactionAdvance, ReactionEdge, ReactionPlan, source::ReactionSourcePlan},
};

use super::super::{
    Accounts, DelayedGroup, DiffusionEvolution, SOURCE_NEUTRONS, SOURCE_PRECURSORS,
    account_add, checked, owner_call,
    source::{PreparedPlan, SourceRates, injected_counts},
};
use super::transport::{InventoryOwner, InventoryTrial};
use super::{
    DiffusionBurnupBudget, DiffusionBurnupError, DiffusionBurnupSources, DiffusionBurnupState,
    DiffusionState, DiffusionStep, export_supplied_count, mean_population, upper_value,
};
use crate::burnup::transport::{
    BurnupTransportBudget, BurnupTransportEvolution, BurnupTransportLimits, scale_yield,
};
use crate::burnup::{BurnupError, BurnupInventory, BurnupModel, BurnupRate, BurnupStep};
use crate::irradiation::GroupFlux;
use crate::{Accumulator, FissionError, nonnegative, positive, sum};

/// Parent-wide radioactive neutron and precursor yields on an exact burnup model.
/// Neutron-induced reactions and spontaneous decays share that model's stock.
#[derive(Clone, Debug)]
pub struct BurnupDecayNeutronEmission {
    data_id: String,
    model: Arc<BurnupModel>,
    neutron_yields: Vec<Vec<f64>>,
    delayed_groups: Vec<DelayedGroup>,
    precursor_yields: Vec<Vec<f64>>,
}

impl BurnupDecayNeutronEmission {
    pub fn new(
        data_id: String,
        model: Arc<BurnupModel>,
        neutron_yields: Vec<Vec<f64>>,
    ) -> Result<Self, DiffusionBurnupError> {
        let precursor_yields = vec![vec![]; model.nuclides().len()];
        Self::new_with_precursors(data_id, model, neutron_yields, vec![], precursor_yields)
    }

    /// Supply immediate group yields and births into explicitly identified
    /// native delayed families per spontaneous parent decay. Family spectra
    /// describe the eventual emitted neutron; they do not create it at birth.
    /// Every nuclide has one row in both arrays, including stable zero rows.
    /// All family metadata must match the bound diffusion model exactly.
    pub fn new_with_precursors(
        data_id: String,
        model: Arc<BurnupModel>,
        neutron_yields: Vec<Vec<f64>>,
        delayed_groups: Vec<DelayedGroup>,
        precursor_yields: Vec<Vec<f64>>,
    ) -> Result<Self, DiffusionBurnupError> {
        if data_id.trim().is_empty()
            || neutron_yields.len() != model.nuclides().len()
            || precursor_yields.len() != model.nuclides().len()
        {
            return Err(DiffusionBurnupError::Invalid(
                "identified neutron and precursor yields in original burnup nuclide order required",
            ));
        }
        if delayed_groups.len() > MAX_DECAY_NODES {
            return Err(BurnupError::Budget("radioactive delayed-family width").into());
        }
        for family in &delayed_groups {
            if family.data_id.trim().is_empty()
                || !positive(family.decay_constant_s_inv)
                || family.spectrum.len() != model.groups().len()
                || !family.spectrum.iter().copied().all(nonnegative)
                || !positive(sum(&family.spectrum))
                || (sum(&family.spectrum) - 1.0).abs() > 64.0 * f64::EPSILON
            {
                return Err(DiffusionBurnupError::Invalid(
                    "identified delayed families require positive decay and a normalized spectrum",
                ));
            }
        }
        for ((nuclide, row), precursors) in model
            .nuclides()
            .iter()
            .zip(&neutron_yields)
            .zip(&precursor_yields)
        {
            if row.len() != model.groups().len()
                || precursors.len() != delayed_groups.len()
                || !row.iter().chain(precursors).copied().all(nonnegative)
                || !sum(row).is_finite()
                || !sum(precursors).is_finite()
                || (nuclide.law.decay_constant().value() == 0.0
                    && row.iter().chain(precursors).any(|&value| value != 0.0))
            {
                return Err(DiffusionBurnupError::Invalid(
                    "finite nonnegative group/family yields and zero stable-nuclide emission required",
                ));
            }
        }
        Ok(Self {
            data_id,
            model,
            neutron_yields,
            delayed_groups,
            precursor_yields,
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
    pub fn delayed_groups(&self) -> &[DelayedGroup] {
        &self.delayed_groups
    }
    pub fn precursor_yields(&self) -> &[Vec<f64>] {
        &self.precursor_yields
    }
}

pub(super) struct DecayBirths {
    pub(super) neutrons: Vec<Vec<f64>>,
    pub(super) precursors: Vec<Vec<f64>>,
}

impl DiffusionBurnupState {
    /// Enable radioactive birth on the original stocks without resetting any
    /// inventory, resident neutron, event count or accumulated energy account.
    /// Every cell supplies its exact original model, including explicit zero
    /// rows. Transported stocks retain their original movement and require the
    /// same parent-wide yields in every compartment: movement cannot change a
    /// nuclide's decay spectrum or family yields on this common neutron grid.
    /// An immediate-only source supplies zero precursor births even when the
    /// neutron model already contains delayed families from other reactions.
    pub fn with_decay_neutrons(
        mut self,
        emissions: Vec<Arc<BurnupDecayNeutronEmission>>,
    ) -> Result<Self, DiffusionBurnupError> {
        if emissions.len() != self.inventories().len() {
            return Err(DiffusionBurnupError::Invalid(
                "decay-neutron coupling requires every original cell inventory",
            ));
        }
        if self.transport().is_some()
            && emissions
                .windows(2)
                .any(|pair| {
                    pair[0].neutron_yields() != pair[1].neutron_yields()
                        || pair[0].delayed_groups() != pair[1].delayed_groups()
                        || pair[0].precursor_yields() != pair[1].precursor_yields()
                })
        {
            return Err(DiffusionBurnupError::Invalid(
                "transported nuclides require unchanged radioactive neutron and family yields",
            ));
        }
        let mut nodes = self.neutrons.model.nodes;
        for (source, inventory) in emissions.iter().zip(self.inventories()) {
            if !Arc::ptr_eq(source.model(), inventory.model()) {
                return Err(DiffusionBurnupError::Invalid(
                    "decay-neutron source must retain its exact original burnup model",
                ));
            }
            if !source.delayed_groups().is_empty()
                && source.delayed_groups() != self.neutrons.model.delayed_groups()
            {
                return Err(DiffusionBurnupError::Invalid(
                    "radioactive delayed families must match the original neutron model exactly",
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
    ) -> Result<(DiffusionState, DiffusionStep, InventoryTrial, DecayBirths), DiffusionBurnupError>
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
        let families = neutron_initial.model.delayed_groups().len();
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
        let reaction_rates = rate_accounts
            .iter()
            .map(|row| row.iter().map(|rate| rate.value()).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let transported = self.transport().is_some();
        let mut channels = Vec::with_capacity(inventories.len());
        let (nuclear_rates, nuclear_edges) = match &self.inventories {
            InventoryOwner::Transported { state, limits } => state.compile_channels(
                &reaction_rates,
                state.transport(),
                BurnupTransportLimits {
                    burnup: budget.burnup.limits,
                    ..*limits
                },
                checkpoint,
            )?,
            InventoryOwner::Independent(_) => {
                let mut rates = Vec::with_capacity(nodes - neutron_nodes);
                let mut edges = Vec::new();
                for (inventory, reaction_rates) in inventories.iter().zip(&reaction_rates) {
                    checkpoint()?;
                    let offset = rates.len();
                    let compiled = inventory
                        .model()
                        .prepare_channels(reaction_rates, &mut *checkpoint)?;
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
                    channels.push(compiled);
                }
                (rates, edges)
            }
        };
        let mut rates = neutron_initial.model.removal_rates.clone();
        rates.extend_from_slice(&nuclear_rates);
        let mut edges = neutron_initial.model.edges.clone();
        for edge in nuclear_edges {
            checkpoint()?;
            push_edge(
                &mut edges,
                ReactionEdge {
                    parent: neutron_nodes + edge.parent,
                    daughter: neutron_nodes + edge.daughter,
                    yield_per_removal: edge.yield_per_removal,
                },
                limits.maximum_edges,
            )?;
        }
        let mut offset = neutron_nodes;
        for (cell, inventory) in inventories.iter().enumerate() {
            checkpoint()?;
            let mut fractions = vec![Accumulator::new(0.0); inventory.populations().len()];
            for (index, reaction) in inventory.model().reactions().iter().enumerate() {
                checkpoint()?;
                if matches!(reaction.rate, BurnupRate::Decay { .. }) {
                    fractions[reaction.parent].add(if transported {
                        reaction_rates[cell][index]
                    } else {
                        channels[cell].fractions[index]
                    })?;
                }
            }
            for (parent, (row, precursors)) in emissions[cell]
                .neutron_yields()
                .iter()
                .zip(emissions[cell].precursor_yields())
                .enumerate()
            {
                checkpoint()?;
                let spontaneous = fractions[parent].value();
                for (output, &yield_) in row.iter().chain(precursors).enumerate() {
                    checkpoint()?;
                    let coefficient = if spontaneous == 0.0 || yield_ == 0.0 {
                        0.0
                    } else if transported {
                        // The original transport owner includes departure in
                        // this hazard. Only spontaneous reactions emit; neither
                        // transfer nor exterior outflow is a radioactive event.
                        scale_yield(yield_, spontaneous, rates[offset + parent])?
                    } else {
                        product(spontaneous, yield_)?
                    };
                    if coefficient > 0.0 {
                        push_edge(
                            &mut edges,
                            ReactionEdge {
                                parent: offset + parent,
                                daughter: if output < groups {
                                    cell * groups + output
                                } else {
                                    neutron_initial.model.neutron_nodes
                                        + cell * families
                                        + output - groups
                                },
                                yield_per_removal: coefficient,
                            },
                            limits.maximum_edges,
                        )?;
                    }
                }
            }
            offset += inventory.populations().len();
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
        let trial = if let InventoryOwner::Transported { state, limits } = &self.inventories {
            let (candidate, report) = state.account_evolution(
                duration_s,
                &reaction_rates,
                state.transport(),
                BurnupTransportEvolution {
                    evolved: ReactionAdvance {
                        populations: evolved.populations[neutron_nodes..].to_vec(),
                        removals: evolved.removals[neutron_nodes..].to_vec(),
                        balance_residuals: evolved.balance_residuals[neutron_nodes..].to_vec(),
                        population_balance_residual: sum(
                            &evolved.balance_residuals[neutron_nodes..]
                        ),
                    },
                    residence: residence[neutron_nodes..].to_vec(),
                    produced_nuclei: produced,
                    input_roundoff_nuclei: input_errors
                        .iter()
                        .map(upper_value)
                        .collect::<Result<_, _>>()?,
                },
                BurnupTransportBudget {
                    limits: BurnupTransportLimits {
                        burnup: budget.burnup.limits,
                        ..*limits
                    },
                    nuclei_tolerance: budget.burnup.nuclei_tolerance,
                    energy_tolerance_j: budget.burnup.energy_tolerance_j,
                },
                checkpoint,
            )?;
            InventoryTrial::from_transport(state, candidate, *limits, report, checkpoint)?
        } else {
            let mut candidates = Vec::with_capacity(inventories.len());
            let mut reports = Vec::with_capacity(inventories.len());
            let mut means = Vec::with_capacity(inventories.len());
            let mut produced_by_cell = Vec::with_capacity(inventories.len());
            let mut offset = neutron_nodes;
            let mut nuclei_error = Accumulator::new(0.0);
            let mut energy_error = Accumulator::new(0.0);
            for (cell, inventory) in inventories.iter().enumerate() {
                checkpoint()?;
                let n = inventory.populations().len();
                let (candidate, mut report) = inventory.account_evolution(
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
            InventoryTrial {
                owner: InventoryOwner::Independent(candidates),
                reports,
                means,
                produced,
                produced_by_cell,
                transport: None,
            }
        };
        let births = account_births(
            inventories,
            emissions,
            &trial.reports,
            groups,
            families,
            &mut injected,
            &mut accounts,
            checkpoint,
        )?;
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
        Ok((neutrons, report, trial, births))
    }
}

// Reported local radioactive events fund stationary and moved neutron/family sources.
// This attaches an injection account to the accepted neutron trajectory; it
// never adds the emitted particles to the propagator a second time.
fn account_births(
    inventories: &[BurnupInventory],
    emissions: &[Arc<BurnupDecayNeutronEmission>],
    reports: &[BurnupStep],
    groups: usize,
    families: usize,
    injected: &mut [f64],
    accounts: &mut Accounts,
    checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
) -> Result<DecayBirths, DiffusionBurnupError> {
    let mut births = DecayBirths {
        neutrons: Vec::with_capacity(inventories.len()),
        precursors: Vec::with_capacity(inventories.len()),
    };
    let neutron_nodes = inventories.len() * groups;
    for (cell, (inventory, report)) in inventories.iter().zip(reports).enumerate() {
        checkpoint()?;
        let n = inventory.populations().len();
        let mut decays = vec![Accumulator::new(0.0); n];
        let mut decay_count_errors = vec![Accumulator::new(0.0); n];
        for (reaction, &count) in inventory
            .model()
            .reactions()
            .iter()
            .zip(&report.reaction_counts)
        {
            checkpoint()?;
            if matches!(reaction.rate, BurnupRate::Decay { .. }) {
                decay_count_errors[reaction.parent]
                    .add(decays[reaction.parent].add(count)?.abs())?;
            }
        }
        let mut born = vec![Accumulator::new(0.0); groups + families];
        for (((&decays, mut count_error), row), precursors) in decays
            .iter()
            .zip(decay_count_errors)
            .zip(emissions[cell].neutron_yields())
            .zip(emissions[cell].precursor_yields())
        {
            checkpoint()?;
            let count = checked(decays.value())?;
            count_error.add(checked(sum(&[decays.hi, decays.lo, -count]))?.abs())?;
            let count_error = upper_value(&count_error)?;
            for (output, &yield_) in row.iter().chain(precursors).enumerate() {
                checkpoint()?;
                // Convert parent events to emitted particles before charging
                // their common native neutron/precursor count allowance.
                let source_error = accounts.multiply(count_error, yield_)?;
                accounts.error.add(source_error)?;
                let value = accounts.multiply(count, yield_)?;
                account_add(&mut born[output], value, &mut accounts.error)?;
            }
        }
        let mut neutron_births = Vec::with_capacity(groups);
        let mut precursor_births = Vec::with_capacity(families);
        for (output, count) in born.into_iter().enumerate() {
            checkpoint()?;
            let value = checked(count.value())?;
            accounts
                .error
                .add(checked(sum(&[count.hi, count.lo, -value]))?.abs())?;
            let index = if output < groups {
                cell * groups + output
            } else {
                neutron_nodes + cell * families + output - groups
            };
            let mut total = Accumulator::new(injected[index]);
            account_add(&mut total, value, &mut accounts.error)?;
            injected[index] = checked(total.value())?;
            accounts
                .error
                .add(checked(sum(&[total.hi, total.lo, -injected[index]]))?.abs())?;
            if output < groups {
                accounts.total(SOURCE_NEUTRONS, value)?;
                neutron_births.push(value);
            } else {
                accounts.total(SOURCE_PRECURSORS, value)?;
                precursor_births.push(value);
            }
        }
        births.neutrons.push(neutron_births);
        births.precursors.push(precursor_births);
    }
    Ok(births)
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

#[cfg(test)]
mod transport_tests;
