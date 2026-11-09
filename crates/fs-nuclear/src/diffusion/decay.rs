//! Radioactive chains and their emitted neutrons on one positive evolution.
//!
//! For source inventories N and spatial neutron/precursor counts n, the frozen
//! operator is [A_neutron, Y diag(lambda); 0, A_decay]. The existing reaction
//! propagator evolves this block once, including daughter ingrowth and neutron
//! transport/reactions DURING emission. Neither endpoint activity nor uniform
//! injection of the interval's final decay count replaces this birth history.
//!
//! Parent-wide neutron multiplicities are authored or imported from evaluated
//! spectra. Branching is already included in these yields. Source nuclides keep
//! their original decay/energy accounts; neutron reactions keep their original
//! residence, fluence, interface and leakage accounts. This count-level owner
//! does not interpret the emitted kinetic energy as extra local heat.
//!
//! The conventional precursor source equations are described in
//! <https://arfc.github.io/moltres/getting_started/theory.html>. Here the source
//! is a complete original nuclide chain, rather than a fitted precursor group.
//! Source nuclides are stationary and are not depleted by incident neutrons;
//! targets in the diffusion model remain its prescribed macroscopic background.
//! The original diffusion approximation and 64-node joint bound remain in force.

use std::sync::Arc;

use fs_la::decay::reaction::{ReactionEdge, ReactionPlan};
use fs_la::decay::{DecayAdvance, DecayLimits, MAX_DECAY_EDGES, MAX_DECAY_NODES};

use super::{
    Accounts, DiffusionBudget, DiffusionError, DiffusionEvolution, DiffusionState, DiffusionStep,
    SOURCE_NEUTRONS, account_add, checked, owner_call, upward_account_bound,
};
use crate::decay::{DecayBudget, DecayError, DecayInventory, DecayModel, DecayStep};
use crate::irradiation::NeutronGroups;
use crate::{Accumulator, FissionError, nonnegative, sum};

/// Disjoint radioactive emission, transported neutron energy and local heat.
pub mod energy;

/// Immutable parent-wide neutron yields on an original decay model and grid.
/// These are neutrons per parent transformation, not branch probabilities.
#[derive(Clone, Debug)]
pub struct DecayNeutronEmission {
    data_id: String,
    model: Arc<DecayModel>,
    groups: Arc<NeutronGroups>,
    neutron_yields: Vec<Vec<f64>>,
}

impl DecayNeutronEmission {
    /// Supply every original nuclide/group coefficient explicitly. Stable rows
    /// must be zero. No missing spectrum, multiplicity or energy is inferred.
    pub fn new(
        data_id: String,
        model: Arc<DecayModel>,
        groups: Arc<NeutronGroups>,
        neutron_yields: Vec<Vec<f64>>,
    ) -> Result<Self, DecayNeutronError> {
        if data_id.trim().is_empty()
            || neutron_yields.len() != model.nuclides().len()
            || groups.len() > MAX_DECAY_NODES
        {
            return Err(DecayNeutronError::Invalid(
                "identified neutron yields in the original nuclide/group layout required",
            ));
        }
        for (node, row) in model.nuclides().iter().zip(&neutron_yields) {
            if row.len() != groups.len()
                || !row.iter().copied().all(nonnegative)
                || !sum(row).is_finite()
                || (node.law.decay_constant().value() == 0.0
                    && row.iter().any(|&value| value != 0.0))
            {
                return Err(DecayNeutronError::Invalid(
                    "finite nonnegative group yields and zero stable-nuclide emission required",
                ));
            }
        }
        Ok(Self {
            data_id,
            model,
            groups,
            neutron_yields,
        })
    }

    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
    }
    pub fn groups(&self) -> &Arc<NeutronGroups> {
        &self.groups
    }
    pub fn neutron_yields(&self) -> &[Vec<f64>] {
        &self.neutron_yields
    }
}

/// One independent finite source inventory at an original spatial cell.
#[derive(Clone, Debug)]
pub struct DecayNeutronSource {
    pub cell: usize,
    pub inventory: DecayInventory,
    pub emission: Arc<DecayNeutronEmission>,
}

/// Original extensive accounts and one common joint propagation allowance.
/// Both component propagation limits must be identical. The joint node/edge
/// and work bounds apply once to the whole coupled preparation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayNeutronBudget {
    pub diffusion: DiffusionBudget,
    /// Population and energy defects are summed across all source inventories.
    pub decay: DecayBudget,
}

#[derive(Clone, Debug)]
struct PreparedDecayNeutrons {
    plan: ReactionPlan,
    limits: DecayLimits,
}

/// Accepted neutron and original nuclide histories. Clone shares only immutable
/// source data and an optional prepared map; no observation redraws a source.
#[derive(Clone, Debug)]
pub struct DecayNeutronState {
    neutrons: DiffusionState,
    sources: Vec<DecayNeutronSource>,
    nodes: usize,
    prepared: Option<Arc<PreparedDecayNeutrons>>,
}

/// One jointly evolved interval; no source count is a post-step impulse.
#[derive(Clone, Debug)]
pub struct DecayNeutronStep {
    pub neutrons: DiffusionStep,
    /// Original source order, each with its original nuclide and branch order.
    pub decay: Vec<DecayStep>,
    /// Integrated neutrons born per source and group. These are the internal
    /// source-to-neutron transfer, also in neutrons.injected_neutron_counts.
    pub emitted_neutrons: Vec<Vec<f64>>,
    pub decay_population_roundoff: f64,
    pub decay_energy_roundoff_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DecayNeutronError {
    Invalid(&'static str),
    Budget(&'static str),
    Decay(DecayError),
    Diffusion(DiffusionError),
    Fission(FissionError),
}
impl From<DecayError> for DecayNeutronError {
    fn from(error: DecayError) -> Self {
        Self::Decay(error)
    }
}
impl From<DiffusionError> for DecayNeutronError {
    fn from(error: DiffusionError) -> Self {
        Self::Diffusion(error)
    }
}
impl From<FissionError> for DecayNeutronError {
    fn from(error: FissionError) -> Self {
        Self::Fission(error)
    }
}
impl core::fmt::Display for DecayNeutronError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radioactive neutron source refused: {self:?}")
    }
}
impl std::error::Error for DecayNeutronError {}

impl DecayNeutronState {
    /// Bind existing histories without replaying earlier decays or discarding
    /// neutrons already in flight. Several independent sources may share a cell.
    pub fn new(
        neutrons: DiffusionState,
        sources: Vec<DecayNeutronSource>,
    ) -> Result<Self, DecayNeutronError> {
        let mut nodes = neutrons.model.nodes;
        if sources.len() > MAX_DECAY_NODES {
            return Err(DecayNeutronError::Budget("radioactive neutron sources"));
        }
        for source in &sources {
            if source.cell >= neutrons.model().cell_count()
                || !Arc::ptr_eq(source.inventory.model(), source.emission.model())
                || source.emission.groups().as_ref() != neutrons.model().groups().as_ref()
            {
                return Err(DecayNeutronError::Invalid(
                    "source requires its exact original decay model and the spatial neutron grid",
                ));
            }
            nodes = nodes
                .checked_add(source.inventory.populations().len())
                .ok_or(DecayNeutronError::Budget("joint decay/neutron nodes"))?;
            if nodes > MAX_DECAY_NODES {
                return Err(DecayNeutronError::Budget("joint decay/neutron nodes"));
            }
        }
        Ok(Self {
            neutrons,
            sources,
            nodes,
            prepared: None,
        })
    }

    pub fn neutrons(&self) -> &DiffusionState {
        &self.neutrons
    }
    pub fn sources(&self) -> &[DecayNeutronSource] {
        &self.sources
    }

    /// Evolve the one unchanged-coefficient block and independently audit both
    /// physical owners. Failure at any checkpoint publishes neither history.
    pub fn trial(
        &self,
        duration_s: f64,
        budget: DecayNeutronBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DecayNeutronStep), DecayNeutronError> {
        checkpoint()?;
        if !nonnegative(duration_s)
            || !nonnegative(budget.diffusion.population_tolerance)
            || !nonnegative(budget.decay.nuclei_tolerance)
            || !nonnegative(budget.decay.energy_tolerance_j)
            || budget.diffusion.limits.propagation != budget.decay.propagation
        {
            return Err(DecayNeutronError::Invalid(
                "finite nonnegative interval/accounts and identical joint propagation limits required",
            ));
        }
        self.neutrons.model.check_limits(budget.diffusion.limits)?;
        let limits = budget.decay.propagation;
        if self.nodes > limits.maximum_nodes {
            return Err(DecayNeutronError::Budget("joint decay/neutron nodes"));
        }
        let prepared = if let Some(prepared) = self.prepared.as_ref().filter(|prepared| {
            prepared.limits == limits
                && prepared.plan.duration_s().to_bits() == duration_s.to_bits()
        }) {
            Arc::clone(prepared)
        } else {
            Arc::new(self.prepare(duration_s, limits, &mut checkpoint)?)
        };
        let neutron_nodes = self.neutrons.model.nodes;
        let groups = self.neutrons.model.groups().len();
        let mut initial = Vec::with_capacity(self.nodes);
        initial.extend_from_slice(self.neutrons.neutron_counts());
        initial.extend_from_slice(self.neutrons.precursor_counts());
        for source in &self.sources {
            initial.extend_from_slice(source.inventory.populations());
        }
        let evolved = owner_call(&mut checkpoint, |cancelled| {
            prepared.plan.apply_with_cancel(&initial, cancelled)
        })?;
        let integrals = owner_call(&mut checkpoint, |cancelled| {
            prepared
                .plan
                .integrate_populations_with_cancel(&initial, cancelled)
        })?;
        let mut accounts = Accounts::new(neutron_nodes)?;
        let mut injected = vec![Accumulator::new(0.0); neutron_nodes];
        let mut emitted_neutrons = Vec::with_capacity(self.sources.len());
        let mut decay = Vec::with_capacity(self.sources.len());
        let mut sources = self.sources.clone();
        let mut offset = neutron_nodes;
        let mut population_defect = Accumulator::new(0.0);
        let mut energy_defect = Accumulator::new(0.0);
        for source in &mut sources {
            checkpoint()?;
            let length = source.inventory.populations().len();
            let counts = &evolved.removals[offset..offset + length];
            let mut emitted = vec![Accumulator::new(0.0); groups];
            for (&decays, row) in counts.iter().zip(source.emission.neutron_yields()) {
                checkpoint()?;
                for (group, &yield_) in row.iter().enumerate() {
                    let born = accounts.multiply(decays, yield_)?;
                    account_add(&mut emitted[group], born, &mut accounts.error)?;
                    account_add(
                        &mut injected[source.cell * groups + group],
                        born,
                        &mut accounts.error,
                    )?;
                    accounts.total(SOURCE_NEUTRONS, born)?;
                }
            }
            let emitted = export_counts(emitted, &mut accounts)?;
            emitted_neutrons.push(emitted);
            let mut balance = Accumulator::new(0.0);
            for index in offset..offset + length {
                balance.add(evolved.populations[index])?;
                balance.add(-initial[index])?;
            }
            let mut checkpoint_error = None;
            let result = source.inventory.account_evolution(
                DecayAdvance {
                    populations: evolved.populations[offset..offset + length].to_vec(),
                    decays: counts.to_vec(),
                    population_balance_residual: checked(balance.value())?,
                },
                duration_s,
                None,
                0.0,
                budget.decay,
                &mut || match checkpoint() {
                    Ok(()) => Ok(()),
                    Err(error) => {
                        checkpoint_error = Some(error);
                        Err(DecayError::Cancelled)
                    }
                },
            );
            if let Some(error) = checkpoint_error {
                return Err(error.into());
            }
            let (inventory, report) = result?;
            population_defect.add(report.population_roundoff_nuclei)?;
            energy_defect.add(report.energy_roundoff_j)?;
            source.inventory = inventory;
            decay.push(report);
            offset += length;
        }
        let decay_population_roundoff = upward_account_bound(population_defect)?;
        let decay_energy_roundoff_j = upward_account_bound(energy_defect)?;
        if decay_population_roundoff > budget.decay.nuclei_tolerance
            || decay_energy_roundoff_j > budget.decay.energy_tolerance_j
        {
            return Err(DecayNeutronError::Budget(
                "aggregate radioactive source accounts",
            ));
        }
        let injected = export_counts(injected, &mut accounts)?;
        let (neutrons, neutron_step) = self.neutrons.account_evolution(
            duration_s,
            DiffusionEvolution {
                initial: initial[..neutron_nodes].to_vec(),
                final_counts: evolved.populations[..neutron_nodes].to_vec(),
                removals: evolved.removals[..neutron_nodes].to_vec(),
                integrals: integrals[..neutron_nodes].to_vec(),
                injected,
                accounts,
            },
            budget.diffusion,
            None,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok((
            Self {
                neutrons,
                sources,
                nodes: self.nodes,
                prepared: Some(prepared),
            },
            DecayNeutronStep {
                neutrons: neutron_step,
                decay,
                emitted_neutrons,
                decay_population_roundoff,
                decay_energy_roundoff_j,
            },
        ))
    }

    /// Publish only the complete accepted source/transport history.
    pub fn advance(
        &mut self,
        duration_s: f64,
        budget: DecayNeutronBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayNeutronStep, DecayNeutronError> {
        let (next, report) = self.trial(duration_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    fn prepare(
        &self,
        duration_s: f64,
        limits: DecayLimits,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<PreparedDecayNeutrons, DecayNeutronError> {
        let mut rates = self.neutrons.model.removal_rates.clone();
        let mut edges = self.neutrons.model.edges.clone();
        let groups = self.neutrons.model.groups().len();
        let mut push = |edge| -> Result<(), DecayNeutronError> {
            if edges.len() >= limits.maximum_edges || edges.len() >= MAX_DECAY_EDGES {
                return Err(DecayNeutronError::Budget("joint decay/neutron edges"));
            }
            edges.push(edge);
            Ok(())
        };
        for source in &self.sources {
            checkpoint()?;
            let offset = rates.len();
            rates.extend(
                source
                    .emission
                    .model()
                    .nuclides()
                    .iter()
                    .map(|nuclide| nuclide.law.decay_constant().value()),
            );
            for branch in source.emission.model().branches() {
                checkpoint()?;
                push(ReactionEdge {
                    parent: offset + branch.parent,
                    daughter: offset + branch.daughter,
                    yield_per_removal: branch.probability,
                })?;
            }
            for (parent, row) in source.emission.neutron_yields().iter().enumerate() {
                checkpoint()?;
                for (group, &yield_) in row.iter().enumerate() {
                    if yield_ > 0.0 {
                        push(ReactionEdge {
                            parent: offset + parent,
                            daughter: source.cell * groups + group,
                            yield_per_removal: yield_,
                        })?;
                    }
                }
            }
        }
        let plan = owner_call(checkpoint, |cancelled| {
            ReactionPlan::prepare_with_population_integrals(
                &rates, &edges, duration_s, limits, cancelled,
            )
        })?;
        Ok(PreparedDecayNeutrons { plan, limits })
    }
}

fn export_counts(
    counts: Vec<Accumulator>,
    accounts: &mut Accounts,
) -> Result<Vec<f64>, DiffusionError> {
    counts
        .into_iter()
        .map(|count| {
            let value = checked(count.value())?;
            accounts
                .error
                .add(checked(sum(&[count.hi, count.lo, -value]))?.abs())?;
            Ok(value)
        })
        .collect()
}

#[cfg(test)]
mod tests;
