//! Constant nuclide production with simultaneous radioactive decay and ingrowth.
//!
//! For an authored rate vector s, the numerical owner integrates N'=A N+s and
//! every parent's transformations over the same interval. Newly produced nuclei
//! have their actual distribution of decay ages; integrated production is not
//! inserted as an instantaneous birth. This matters for short-lived daughters.
//!
//! Rates describe prescribed production of tracked nuclides [expected nuclei/s].
//! They do not infer activation, neutron flux, isotope yields, irradiation-parent
//! depletion or bulk mass/enthalpy transport. A time-varying or depleted source
//! must supply its correct history through an appropriate owner. In particular,
//! a frozen average of an exponentially depleted irradiation rate is not exact.

use std::sync::Arc;

use fs_la::decay::{DecayLimits, source::DecaySourcePlan};

use super::{
    Accumulator, DecayBudget, DecayError, DecayInventory, DecayModel, DecayStep, add, checked,
    nonnegative, product, sum,
};

/// Reusable homogeneous and uniform-production maps for one physical model.
#[derive(Clone, Debug)]
pub struct PreparedDecayProduction {
    model: Arc<DecayModel>,
    plan: DecaySourcePlan,
    limits: DecayLimits,
}

impl PreparedDecayProduction {
    #[must_use]
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
    }
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.plan.homogeneous_plan().duration_s()
    }
    #[must_use]
    pub fn propagator(&self) -> &DecaySourcePlan {
        &self.plan
    }
}

impl DecayModel {
    /// Prepare an interval independently of both initial populations and source
    /// rates. Reusing this object avoids matrix preparation across source changes.
    pub fn prepare_production(
        self: &Arc<Self>,
        duration_s: f64,
        limits: DecayLimits,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<PreparedDecayProduction, DecayError> {
        let mut checkpoint_error = None;
        let plan = DecaySourcePlan::prepare(&self.rates, &self.edges, duration_s, limits, || {
            match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            }
        });
        if let Some(error) = checkpoint_error {
            return Err(error);
        }
        Ok(PreparedDecayProduction {
            model: Arc::clone(self),
            plan: plan?,
            limits,
        })
    }
}

/// One accepted interval's production and the original shared decay accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayProductionStep {
    /// Counts produced at the authored rates, in the model's original node order.
    pub produced_nuclei: Vec<f64>,
    /// Integrated decays and energy from both old and newly produced inventories.
    /// Population balance is final minus initial minus produced, and the
    /// noncancelling defect includes source products and lifetime birth counters.
    pub decay: DecayStep,
}

impl DecayInventory {
    /// Trial constant production and radioactive evolution on one shared interval.
    pub fn trial_with_production(
        &self,
        duration_s: f64,
        rates_nuclei_s: &[f64],
        budget: DecayBudget,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, DecayProductionStep), DecayError> {
        budget.validate()?;
        admit_rates(&self.model, rates_nuclei_s)?;
        let prepared =
            self.model
                .prepare_production(duration_s, budget.propagation, &mut checkpoint)?;
        self.trial_production_prepared(&prepared, rates_nuclei_s, budget, checkpoint)
    }

    /// Apply a reusable source convolution. Rates can change between applications;
    /// physical model identity, duration and propagation limits remain bound.
    pub fn trial_production_prepared(
        &self,
        prepared: &PreparedDecayProduction,
        rates_nuclei_s: &[f64],
        budget: DecayBudget,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, DecayProductionStep), DecayError> {
        checkpoint()?;
        budget.validate()?;
        admit_rates(&self.model, rates_nuclei_s)?;
        if !Arc::ptr_eq(&self.model, &prepared.model) {
            return Err(DecayError::Invalid(
                "prepared production belongs to a different physical model",
            ));
        }
        if prepared.limits != budget.propagation {
            return Err(DecayError::Invalid(
                "prepared production requires the original propagation limits",
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
            return Err(error);
        }
        let (next, decay) = self.account_evolution(
            evolved?,
            duration_s,
            Some(&produced_nuclei),
            checked(production_error.value())?,
            budget,
            &mut checkpoint,
        )?;
        Ok((
            next,
            DecayProductionStep {
                produced_nuclei,
                decay,
            },
        ))
    }

    /// Publish only an entirely accepted production/decay candidate.
    pub fn advance_with_production(
        &mut self,
        duration_s: f64,
        rates_nuclei_s: &[f64],
        budget: DecayBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<DecayProductionStep, DecayError> {
        let (next, report) =
            self.trial_with_production(duration_s, rates_nuclei_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    /// Atomically publish an application of a reusable production convolution.
    pub fn advance_production_prepared(
        &mut self,
        prepared: &PreparedDecayProduction,
        rates_nuclei_s: &[f64],
        budget: DecayBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<DecayProductionStep, DecayError> {
        let (next, report) =
            self.trial_production_prepared(prepared, rates_nuclei_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

/// Right-continuous change of the complete original-order production vector.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayProductionChange {
    pub time_s: f64,
    pub rates_nuclei_s: Vec<f64>,
}

/// Immutable absolute-time production history, bound to one physical decay model.
/// Transition lookup is binary search; intervals must be clipped at transitions.
#[derive(Clone, Debug, PartialEq)]
pub struct DecayProductionSchedule {
    model: Arc<DecayModel>,
    initial_rates_nuclei_s: Vec<f64>,
    changes: Vec<DecayProductionChange>,
}

impl DecayProductionSchedule {
    /// Admit explicit nonnegative production and a bounded increasing history.
    /// The initial vector applies at time zero; every change takes effect exactly
    /// at its authored time. Supply zero rates to shut production down.
    pub fn new(
        model: Arc<DecayModel>,
        initial_rates_nuclei_s: Vec<f64>,
        changes: Vec<DecayProductionChange>,
        maximum_changes: usize,
    ) -> Result<Self, DecayError> {
        if changes.len() > maximum_changes {
            return Err(DecayError::Budget("production schedule changes"));
        }
        admit_rates(&model, &initial_rates_nuclei_s)?;
        let mut previous = 0.0;
        for change in &changes {
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(DecayError::Invalid(
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
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
    }
    #[must_use]
    pub fn changes(&self) -> &[DecayProductionChange] {
        &self.changes
    }

    /// Active rates at an absolute finite nonnegative time [s].
    pub fn rates_at(&self, time_s: f64) -> Result<&[f64], DecayError> {
        if !nonnegative(time_s) {
            return Err(DecayError::Invalid(
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

fn admit_rates(model: &DecayModel, rates: &[f64]) -> Result<(), DecayError> {
    if rates.len() != model.nuclides.len() || !rates.iter().copied().all(nonnegative) {
        return Err(DecayError::Invalid(
            "one finite nonnegative production rate per original nuclide is required",
        ));
    }
    // Individual rate dimensions are admitted here. A sum of rates need not be
    // representable when their finite-duration produced counts are representable.
    Ok(())
}

#[cfg(test)]
mod tests;
