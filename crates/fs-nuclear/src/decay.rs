//! Identified radioactive inventories on the shared positive decay propagator.
//!
//! Every node retains an admitted `DecayLawV1`; stable daughters are explicit.
//! A transformation selects one represented daughter. Branch probabilities and
//! emitted energy are supplied by an identified dataset, never inferred from
//! chemical species or an effective fission-energy reservoir. Counts are
//! continuous expectations, not stochastic integer events.
//!
//! The numerical owner `fs_la::decay` integrates both daughter production and
//! each parent's transformations. Energy uses those integrated counts, including
//! daughters born and decayed within the interval. Endpoint activity times the
//! interval is not a substitute. Preparation is reusable across inventories
//! bound to the same immutable model; application is quadratic in node count.
//!
//! The base advance is closed and acyclic; [`production`] adds explicitly
//! prescribed production throughout an interval. Neither operation infers
//! irradiation, fission yields, rate feedback or radiation transport. One
//! retained daughter per event does not assert conservation of nuclear
//! mass/charge without emitted-particle data.

use std::{collections::BTreeSet, sync::Arc};

use fs_la::decay::{DecayEdge, DecayLimits, DecayPlan};
use fs_thermochem::decay::DecayLawV1;

use crate::{Accumulator, EnergyChannels, channel_values, nonnegative, sum};

/// Decay heat on the original material-backed thermal timeline.
#[cfg(feature = "thermal")]
pub mod thermal;

/// Prescribed production during decay and source-bound production schedules.
pub mod production;

/// A dataset-owned nuclide identity and its immutable admitted half-life law.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayNuclide {
    pub id: String,
    pub law: DecayLawV1,
}

/// One alternative decay channel, retaining its own energy destinations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayBranch {
    pub parent: usize,
    pub daughter: usize,
    /// Positive probability. All channels of an unstable parent sum to one.
    pub probability: f64,
    /// Joules per event in THIS branch, not probability-weighted coefficients.
    pub energy_per_decay: EnergyChannels,
}

/// Immutable original-order nuclide data, topology and emission coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayModel {
    data_id: String,
    nuclides: Vec<DecayNuclide>,
    branches: Vec<DecayBranch>,
    rates: Vec<f64>,
    edges: Vec<DecayEdge>,
}

impl DecayModel {
    /// Validate identities, channels, stable endpoints and the complete DAG.
    /// The data label is provenance supplied by the caller, not calibration.
    /// No branch normalization, isotope table or equilibrium prehistory is added.
    pub fn new(
        data_id: String,
        nuclides: Vec<DecayNuclide>,
        branches: Vec<DecayBranch>,
        limits: DecayLimits,
    ) -> Result<Self, DecayError> {
        if data_id.trim().is_empty() {
            return Err(DecayError::Invalid("a decay dataset identity is required"));
        }
        if nuclides.len() > limits.maximum_nodes
            || nuclides.len() > fs_la::decay::MAX_DECAY_NODES
            || branches.len() > limits.maximum_edges
            || branches.len() > fs_la::decay::MAX_DECAY_EDGES
        {
            return Err(DecayError::Budget("decay model size"));
        }
        let mut identities = BTreeSet::new();
        for nuclide in &nuclides {
            if nuclide.id.trim().is_empty() || !identities.insert(nuclide.id.as_str()) {
                return Err(DecayError::Invalid(
                    "nuclide identities must be nonempty and unique",
                ));
            }
        }
        let rates: Vec<_> = nuclides
            .iter()
            .map(|node| node.law.decay_constant().value())
            .collect();
        let mut edges = Vec::with_capacity(branches.len());
        for branch in &branches {
            branch
                .energy_per_decay
                .validate()
                .map_err(|_| DecayError::Invalid("branch energy must be finite and nonnegative"))?;
            edges.push(DecayEdge {
                parent: branch.parent,
                daughter: branch.daughter,
                branch_fraction: branch.probability,
            });
        }
        // Identity preparation still performs the numerical owner's complete
        // rate/topology/branch admission; do not maintain a second graph checker.
        DecayPlan::prepare(&rates, &edges, 0.0, limits, || false)?;
        Ok(Self {
            data_id,
            nuclides,
            branches,
            rates,
            edges,
        })
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    #[must_use]
    pub fn nuclides(&self) -> &[DecayNuclide] {
        &self.nuclides
    }
    #[must_use]
    pub fn branches(&self) -> &[DecayBranch] {
        &self.branches
    }

    /// Prepare once, then reuse for any inventories of this exact model and
    /// duration. Matrix preparation and representability limits belong to fs-la.
    pub fn prepare(
        self: &Arc<Self>,
        duration_s: f64,
        limits: DecayLimits,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<PreparedDecay, DecayError> {
        let mut checkpoint_error = None;
        let plan =
            DecayPlan::prepare(
                &self.rates,
                &self.edges,
                duration_s,
                limits,
                || match checkpoint() {
                    Ok(()) => false,
                    Err(error) => {
                        checkpoint_error = Some(error);
                        true
                    }
                },
            );
        if let Some(error) = checkpoint_error {
            return Err(error);
        }
        Ok(PreparedDecay {
            model: Arc::clone(self),
            plan: plan?,
            limits,
        })
    }
}

/// Shared immutable interval map, bound to the exact admitted physical model.
#[derive(Clone, Debug)]
pub struct PreparedDecay {
    model: Arc<DecayModel>,
    plan: DecayPlan,
    limits: DecayLimits,
}

impl PreparedDecay {
    #[must_use]
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
    }
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.plan.duration_s()
    }
    /// The original numerical map and its work diagnostics, without a second solver.
    #[must_use]
    pub fn propagator(&self) -> &DecayPlan {
        &self.plan
    }
}

/// Explicit numerical resources and noncancelling extensive accounting limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayBudget {
    pub propagation: DecayLimits,
    /// Sum of absolute population/branch/counter defects [expected nuclei].
    pub nuclei_tolerance: f64,
    /// Sum of energy product, summation and cumulative-account defects [J].
    /// These limits do not certify nuclear data or the propagator's forward error.
    pub energy_tolerance_j: f64,
}

impl DecayBudget {
    fn validate(self) -> Result<(), DecayError> {
        if !nonnegative(self.nuclei_tolerance) || !nonnegative(self.energy_tolerance_j) {
            return Err(DecayError::Invalid(
                "decay account tolerances must be finite and nonnegative",
            ));
        }
        Ok(())
    }
}

/// Accepted original-order populations and compensated lifetime emission accounts.
/// Clone copies all mutable history and shares only immutable model data.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayInventory {
    model: Arc<DecayModel>,
    populations: Vec<f64>,
    decays: Vec<Accumulator>,
    produced: Vec<Accumulator>,
    released: [Accumulator; 3],
}

impl DecayInventory {
    /// Supply every initial expected population explicitly, including daughters.
    pub fn new(model: Arc<DecayModel>, populations: Vec<f64>) -> Result<Self, DecayError> {
        if populations.len() != model.nuclides.len()
            || !populations.iter().copied().all(nonnegative)
            || !sum(&populations).is_finite()
        {
            return Err(DecayError::Invalid(
                "one finite nonnegative inventory per nuclide is required",
            ));
        }
        let decays = vec![Accumulator::new(0.0); populations.len()];
        let produced = vec![Accumulator::new(0.0); populations.len()];
        Ok(Self {
            model,
            populations,
            decays,
            produced,
            released: [Accumulator::new(0.0); 3],
        })
    }

    #[must_use]
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
    }
    #[must_use]
    pub fn populations(&self) -> &[f64] {
        &self.populations
    }
    #[must_use]
    pub fn cumulative_decays(&self) -> Vec<f64> {
        self.decays.iter().map(|count| count.value()).collect()
    }
    /// Original-order expected nuclei produced by explicit sources since construction.
    #[must_use]
    pub fn cumulative_produced_nuclei(&self) -> Vec<f64> {
        self.produced.iter().map(|count| count.value()).collect()
    }
    #[must_use]
    pub fn cumulative_released(&self) -> EnergyChannels {
        channel_values(self.released)
    }

    /// Current activity of each original nuclide, lambda*N [events/s = Bq].
    /// An unrepresentable positive activity refuses; stable activity is exactly zero.
    pub fn activity_bq(&self) -> Result<Vec<f64>, DecayError> {
        self.model
            .rates
            .iter()
            .zip(&self.populations)
            .map(|(&rate, &population)| product(rate, population))
            .collect()
    }

    /// Instantaneous emission power in the three original energy destinations.
    /// This is a diagnostic; use integrated `released` for an interval's heat.
    pub fn decay_power_w(&self) -> Result<EnergyChannels, DecayError> {
        let activity = self.activity_bq()?;
        let mut channels = [Accumulator::new(0.0); 3];
        for branch in &self.model.branches {
            let rate = product(activity[branch.parent], branch.probability)?;
            for (total, coefficient) in channels
                .iter_mut()
                .zip(channel_array(branch.energy_per_decay))
            {
                add(total, product(rate, coefficient)?)?;
            }
        }
        let result = channel_values(channels);
        checked(result.total_j())?;
        Ok(result)
    }

    /// Prepare and evaluate a complete candidate without modifying this history.
    /// Repeated intervals can avoid preparation with `trial_prepared`.
    pub fn trial(
        &self,
        duration_s: f64,
        budget: DecayBudget,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, DecayStep), DecayError> {
        budget.validate()?;
        let prepared = self
            .model
            .prepare(duration_s, budget.propagation, &mut checkpoint)?;
        self.trial_prepared(&prepared, budget, checkpoint)
    }

    /// Apply a source-bound prepared map. Solver limits must match preparation;
    /// extensive account tolerances may differ for inventories of different size.
    /// Cancellation and all refusals leave both inventory and propagator intact.
    pub fn trial_prepared(
        &self,
        prepared: &PreparedDecay,
        budget: DecayBudget,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, DecayStep), DecayError> {
        checkpoint()?;
        budget.validate()?;
        if !Arc::ptr_eq(&self.model, &prepared.model) {
            return Err(DecayError::Invalid(
                "prepared decay belongs to a different physical model",
            ));
        }
        if prepared.limits != budget.propagation {
            return Err(DecayError::Invalid(
                "prepared decay requires the original propagation limits",
            ));
        }
        let mut checkpoint_error = None;
        let evolved = prepared
            .plan
            .apply_with_cancel(&self.populations, || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            });
        if let Some(error) = checkpoint_error {
            return Err(error);
        }
        let evolved = evolved?;
        self.account_evolution(
            evolved,
            prepared.duration_s(),
            None,
            0.0,
            budget,
            &mut checkpoint,
        )
    }

    // All numerical propagators share the same physical event/energy ledger.
    // Produced counts are a separately authored influx, not daughter transfers.
    fn account_evolution(
        &self,
        evolved: fs_la::decay::DecayAdvance,
        duration_s: f64,
        produced: Option<&[f64]>,
        production_roundoff_nuclei: f64,
        budget: DecayBudget,
        checkpoint: &mut impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, DecayStep), DecayError> {
        let mut next = self.clone();
        let mut incoming = vec![Accumulator::new(0.0); self.populations.len()];
        let mut released = [Accumulator::new(0.0); 3];
        let mut population_error = Accumulator::new(evolved.population_balance_residual.abs());
        add(&mut population_error, production_roundoff_nuclei)?;
        let mut energy_error = Accumulator::new(0.0);
        let mut branch_decays = Vec::with_capacity(self.model.branches.len());
        for branch in &self.model.branches {
            checkpoint()?;
            let parent_decays = evolved.decays[branch.parent];
            let count = product(parent_decays, branch.probability)?;
            let count_error = parent_decays.mul_add(branch.probability, -count).abs();
            let incoming_error = add(&mut incoming[branch.daughter], count)?.abs();
            add(
                &mut population_error,
                checked(count_error + incoming_error)?,
            )?;
            branch_decays.push(count);
            add(
                &mut energy_error,
                checked(count_error * branch.energy_per_decay.total_j())?,
            )?;
            for (total, coefficient) in released
                .iter_mut()
                .zip(channel_array(branch.energy_per_decay))
            {
                let energy = product(count, coefficient)?;
                add(&mut energy_error, count.mul_add(coefficient, -energy).abs())?;
                let error = add(total, energy)?.abs();
                add(&mut energy_error, error)?;
            }
        }
        for (index, &population) in evolved.populations.iter().enumerate() {
            checkpoint()?;
            let production = produced.map_or(0.0, |counts| counts[index]);
            // Per-node balance catches missed tiny debits even when the total
            // population residual rounds to zero or different nodes cancel.
            let residual = sum(&[
                population,
                -self.populations[index],
                evolved.decays[index],
                -incoming[index].hi,
                -incoming[index].lo,
                -production,
            ]);
            add(&mut population_error, checked(residual)?.abs())?;
            let error = add(&mut next.decays[index], evolved.decays[index])?.abs();
            add(&mut population_error, error)?;
            let error = add(&mut next.produced[index], production)?.abs();
            add(&mut population_error, error)?;
        }
        let released = channel_values(released);
        checked(released.total_j())?;
        for (total, energy) in next.released.iter_mut().zip(channel_array(released)) {
            let error = add(total, energy)?.abs();
            add(&mut energy_error, error)?;
        }
        checked(channel_values(next.released).total_j())?;
        let population_roundoff_nuclei = checked(population_error.value())?;
        let energy_roundoff_j = checked(energy_error.value())?;
        if population_roundoff_nuclei > budget.nuclei_tolerance {
            return Err(DecayError::Budget("population accounting roundoff"));
        }
        if energy_roundoff_j > budget.energy_tolerance_j {
            return Err(DecayError::Budget("emission energy accounting roundoff"));
        }
        next.populations = evolved.populations;
        checkpoint()?;
        Ok((
            next,
            DecayStep {
                duration_s,
                decays: evolved.decays,
                branch_decays,
                released,
                population_balance_residual: evolved.population_balance_residual,
                population_roundoff_nuclei,
                energy_roundoff_j,
            },
        ))
    }

    /// Atomically publish a complete candidate; all errors preserve `self`.
    pub fn advance(
        &mut self,
        duration_s: f64,
        budget: DecayBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<DecayStep, DecayError> {
        let (next, report) = self.trial(duration_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    /// Atomically publish an application of a reusable prepared interval.
    pub fn advance_prepared(
        &mut self,
        prepared: &PreparedDecay,
        budget: DecayBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<DecayStep, DecayError> {
        let (next, report) = self.trial_prepared(prepared, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

/// Integrated transformations and extensive energy of one accepted interval.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayStep {
    pub duration_s: f64,
    /// Original-order counts, including daughters created and decayed in this step.
    pub decays: Vec<f64>,
    /// Counts in the model's original branch order.
    pub branch_decays: Vec<f64>,
    pub released: EnergyChannels,
    /// Final minus initial minus explicitly produced expected nuclei.
    pub population_balance_residual: f64,
    pub population_roundoff_nuclei: f64,
    pub energy_roundoff_j: f64,
}

/// Refusal without accepting a partial physical history.
#[derive(Clone, Debug, PartialEq)]
pub enum DecayError {
    Invalid(&'static str),
    Budget(&'static str),
    Arithmetic,
    Propagation(fs_la::decay::DecayError),
    Cancelled,
}

impl From<fs_la::decay::DecayError> for DecayError {
    fn from(error: fs_la::decay::DecayError) -> Self {
        match error {
            fs_la::decay::DecayError::Cancelled => Self::Cancelled,
            other => Self::Propagation(other),
        }
    }
}
impl core::fmt::Display for DecayError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radioactive inventory refused: {self:?}")
    }
}
impl std::error::Error for DecayError {}

fn channel_array(value: EnergyChannels) -> [f64; 3] {
    [value.deposited_j, value.radiation_j, value.escaped_j]
}
fn checked(value: f64) -> Result<f64, DecayError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DecayError::Arithmetic)
    }
}
fn product(a: f64, b: f64) -> Result<f64, DecayError> {
    let value = checked(a * b)?;
    if value == 0.0 && a > 0.0 && b > 0.0 {
        return Err(DecayError::Arithmetic);
    }
    Ok(value)
}
fn add(target: &mut Accumulator, value: f64) -> Result<f64, DecayError> {
    target.add(value).map_err(|_| DecayError::Arithmetic)
}

#[cfg(test)]
mod tests;
