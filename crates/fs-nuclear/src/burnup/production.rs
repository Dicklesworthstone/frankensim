//! Continuous external nuclide supply with simultaneous burnup and decay.
//!
//! The original reaction owner evolves N'=A N+s over a frozen-rate interval,
//! retaining the ages, subsequent reactions and decay heat of supplied nuclei.
//! Cycles, sinks, self-return and fission multiplicity keep their existing
//! meanings. Production is an explicit expected-nuclei influx, separate from
//! reaction-created products. Its rates do not infer a neutron history, feed
//! enthalpy, emitted-particle balance or a change in thermal bulk material mass.

use std::sync::Arc;

use fs_la::decay::reaction::source::ReactionSourcePlan;

use super::{
    Accumulator, BurnupBudget, BurnupError, BurnupInventory, BurnupLimits, BurnupModel, BurnupStep,
    add, checked, nonnegative, nonnegative_bound, product, sum,
};
use crate::FissionError;
use crate::irradiation::GroupFlux;

/// Reusable production/removal convolution and original reaction attribution.
/// The recorded flux and temperature identify the first preparation point;
/// equivalent complete channel-rate vectors can reuse the same immutable map.
#[derive(Clone, Debug)]
pub struct PreparedBurnupProduction {
    model: Arc<BurnupModel>,
    plan: ReactionSourcePlan,
    flux: GroupFlux,
    temperature_k: f64,
    limits: BurnupLimits,
    rates: Vec<f64>,
    fractions: Vec<f64>,
}

impl PreparedBurnupProduction {
    #[must_use]
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.plan.homogeneous_plan().duration_s()
    }
    #[must_use]
    pub fn flux(&self) -> &GroupFlux {
        &self.flux
    }
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    #[must_use]
    pub fn reaction_rate_constants(&self) -> &[f64] {
        &self.rates
    }
    #[must_use]
    pub fn propagator(&self) -> &ReactionSourcePlan {
        &self.plan
    }
}

impl BurnupModel {
    /// Prepare for arbitrary initial inventories and constant production rates.
    /// The admitted model, interval, channel rates and limits remain bound.
    pub fn prepare_production(
        self: &Arc<Self>,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        limits: BurnupLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<PreparedBurnupProduction, BurnupError> {
        let rates = self.rate_constants(flux, temperature_k, limits, &mut checkpoint)?;
        self.prepare_production_rates(duration_s, flux, temperature_k, limits, rates, checkpoint)
    }

    // Cache misses arrive with already admitted current rates. Neither table
    // interpolation nor neutron-group collapse is repeated during preparation.
    fn prepare_production_rates(
        self: &Arc<Self>,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        limits: BurnupLimits,
        rates: Vec<f64>,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<PreparedBurnupProduction, BurnupError> {
        let channels = self.prepare_channels(&rates, &mut checkpoint)?;
        let mut checkpoint_error = None;
        let plan = ReactionSourcePlan::prepare(
            &channels.removal_rates,
            &channels.edges,
            duration_s,
            limits.propagation,
            || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            },
        );
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        Ok(PreparedBurnupProduction {
            model: Arc::clone(self),
            plan: plan?,
            flux: flux.clone(),
            temperature_k,
            limits,
            rates,
            fractions: channels.fractions,
        })
    }
}

/// Supplied nuclei and the original reaction/event-energy accounts for one step.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupProductionStep {
    /// Authored supply integrated over the interval, in original nuclide order.
    pub produced_nuclei: Vec<f64>,
    /// Actual reactions of old and supplied nuclei. Node balances include supply;
    /// the population defect includes source products and cumulative supply.
    pub burnup: BurnupStep,
}

impl BurnupInventory {
    /// Trial constant external production under the current flux and temperature.
    /// Source-vector changes reuse a map when its complete reaction operator is
    /// unchanged. Every call still admits the current environment and work limits.
    pub fn trial_with_production(
        &self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        rates_nuclei_s: &[f64],
        budget: BurnupBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupProductionStep), BurnupError> {
        checkpoint()?;
        budget.validate()?;
        admit_rates(&self.model, rates_nuclei_s)?;
        let rates =
            self.model
                .rate_constants(flux, temperature_k, budget.limits, &mut checkpoint)?;
        let prepared = if let Some(prepared) = self.production_prepared.as_ref().filter(|p| {
            Arc::ptr_eq(&self.model, p.model())
                && p.duration_s().to_bits() == duration_s.to_bits()
                && p.limits == budget.limits
                && p.reaction_rate_constants() == rates.as_slice()
        }) {
            Arc::clone(prepared)
        } else {
            Arc::new(self.model.prepare_production_rates(
                duration_s,
                flux,
                temperature_k,
                budget.limits,
                rates,
                &mut checkpoint,
            )?)
        };
        let (mut next, report) =
            self.trial_production_prepared(&prepared, rates_nuclei_s, budget, checkpoint)?;
        next.production_prepared = Some(prepared);
        Ok((next, report))
    }

    /// Apply a model-bound convolution without matrix preparation. Production
    /// rates may change; the original physical model and work limits may not.
    pub fn trial_production_prepared(
        &self,
        prepared: &PreparedBurnupProduction,
        rates_nuclei_s: &[f64],
        budget: BurnupBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupProductionStep), BurnupError> {
        checkpoint()?;
        budget.validate()?;
        admit_rates(&self.model, rates_nuclei_s)?;
        if !Arc::ptr_eq(&self.model, prepared.model()) || budget.limits != prepared.limits {
            return Err(BurnupError::Invalid(
                "prepared production requires the original burnup model and work limits",
            ));
        }
        let duration_s = prepared.duration_s();
        let mut produced_nuclei = Vec::with_capacity(rates_nuclei_s.len());
        let mut production_error = Accumulator::new(0.0);
        for &rate in rates_nuclei_s {
            checkpoint()?;
            let count = product(rate, duration_s)?;
            add(
                &mut production_error,
                rate.mul_add(duration_s, -count).abs(),
            )?;
            produced_nuclei.push(count);
        }
        checked(sum(&produced_nuclei))?;
        let mut checkpoint_error = None;
        let evolved = prepared
            .plan
            .apply_with_cancel(&self.populations, &produced_nuclei, || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (mut next, burnup) = self.account_evolution(
            evolved?,
            duration_s,
            &prepared.fractions,
            Some(&produced_nuclei),
            nonnegative_bound(production_error)?,
            budget,
            &mut checkpoint,
        )?;
        // Accepted mode changes retain one map family. A refused candidate has
        // not touched either cache or any compensated physical history.
        next.prepared = None;
        Ok((
            next,
            BurnupProductionStep {
                produced_nuclei,
                burnup,
            },
        ))
    }

    /// Publish the entire supplied/irradiated/decayed inventory after acceptance.
    pub fn advance_with_production(
        &mut self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        rates_nuclei_s: &[f64],
        budget: BurnupBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupProductionStep, BurnupError> {
        let (next, report) = self.trial_with_production(
            duration_s,
            flux,
            temperature_k,
            rates_nuclei_s,
            budget,
            checkpoint,
        )?;
        *self = next;
        Ok(report)
    }

    /// Atomically publish a reusable production-map application.
    pub fn advance_production_prepared(
        &mut self,
        prepared: &PreparedBurnupProduction,
        rates_nuclei_s: &[f64],
        budget: BurnupBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupProductionStep, BurnupError> {
        let (next, report) =
            self.trial_production_prepared(prepared, rates_nuclei_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

/// Right-continuous change of all original-order external production rates.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupProductionChange {
    pub time_s: f64,
    pub rates_nuclei_s: Vec<f64>,
}

/// Bounded absolute-time production history bound to one immutable burnup model.
/// A rate vector has units of expected nuclei/s, not neutron flux or mass flow.
#[derive(Clone, Debug, PartialEq)]
pub struct BurnupProductionSchedule {
    model: Arc<BurnupModel>,
    initial_rates_nuclei_s: Vec<f64>,
    changes: Vec<BurnupProductionChange>,
}

impl BurnupProductionSchedule {
    pub fn new(
        model: Arc<BurnupModel>,
        initial_rates_nuclei_s: Vec<f64>,
        changes: Vec<BurnupProductionChange>,
        maximum_changes: usize,
    ) -> Result<Self, BurnupError> {
        if changes.len() > maximum_changes {
            return Err(BurnupError::Budget("production schedule changes"));
        }
        admit_rates(&model, &initial_rates_nuclei_s)?;
        let mut previous = 0.0;
        for change in &changes {
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(BurnupError::Invalid(
                    "production changes need strictly increasing positive times",
                ));
            }
            admit_rates(&model, &change.rates_nuclei_s)?;
            previous = change.time_s;
        }
        Ok(Self {
            model,
            initial_rates_nuclei_s,
            changes,
        })
    }

    #[must_use]
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }
    #[must_use]
    pub fn changes(&self) -> &[BurnupProductionChange] {
        &self.changes
    }
    /// Active rates at a finite nonnegative absolute time.
    pub fn rates_at(&self, time_s: f64) -> Result<&[f64], BurnupError> {
        if !nonnegative(time_s) {
            return Err(BurnupError::Invalid(
                "production schedule time must be finite and nonnegative",
            ));
        }
        let after = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        Ok(if after == 0 {
            &self.initial_rates_nuclei_s
        } else {
            &self.changes[after - 1].rates_nuclei_s
        })
    }
    /// First transition strictly after an admitted time; invalid times return None.
    #[must_use]
    pub fn next_after(&self, time_s: f64) -> Option<f64> {
        if !nonnegative(time_s) {
            return None;
        }
        self.changes
            .get(
                self.changes
                    .partition_point(|change| change.time_s <= time_s),
            )
            .map(|change| change.time_s)
    }
}

fn admit_rates(model: &BurnupModel, rates: &[f64]) -> Result<(), BurnupError> {
    if rates.len() != model.nuclides().len() || !rates.iter().copied().all(nonnegative) {
        return Err(BurnupError::Invalid(
            "one finite nonnegative external production rate per nuclide required",
        ));
    }
    // A rate sum need not be representable when finite-duration counts are.
    Ok(())
}

#[cfg(test)]
mod tests;
