//! Original multigroup capture and radioactive rates prepared at authored times.
//!
//! Each frame evaluates the same immutable BurnupModel once. Its table laws,
//! absolute flux, nuclide identities and temperature domain remain authoritative.
//! The history prescribes nuclear temperatures; it does not infer feedback from
//! a thermal receiver or interpolate forcing across a discontinuity.

use std::sync::Arc;

use crate::burnup::{BurnupError, BurnupLimits, BurnupModel, BurnupRate, add};
use crate::irradiation::GroupFlux;
use crate::{Accumulator, FissionError, nonnegative};

/// Replace the incident spectrum and declared nuclear temperature at this time.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticActivationChange {
    pub time_s: f64,
    pub flux: GroupFlux,
    pub temperature_k: f64,
}

/// Bounded, right-continuous forcing for one original capture/decay network.
///
/// Every channel has exactly one retained product per actual event. Expected
/// fractional or multiple-product yields and all fission channels are refused;
/// their mean yields do not define a joint stochastic product distribution.
/// Neutron-producing channels are also refused because this event owner has no
/// distinct neutron-production event and energy accounts. They are not capture.
/// The original integer population is consequently conserved by every event.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticActivationSchedule {
    model: Arc<BurnupModel>,
    limits: BurnupLimits,
    initial_flux: GroupFlux,
    initial_temperature_k: f64,
    initial_rates: Vec<f64>,
    changes: Vec<StochasticActivationChange>,
    changed_rates: Vec<Vec<f64>>,
    outgoing: Vec<Vec<usize>>,
}

impl StochasticActivationSchedule {
    /// Evaluate and retain original reaction rates once per admitted frame.
    ///
    /// Times must be finite, positive and strictly increasing. Every spectrum
    /// and temperature passes through the original model, including fixed and
    /// tabulated microscopic data. Cancellation returns its original reason.
    pub fn new(
        model: Arc<BurnupModel>,
        initial_flux: GroupFlux,
        initial_temperature_k: f64,
        changes: Vec<StochasticActivationChange>,
        limits: BurnupLimits,
        maximum_changes: usize,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, BurnupError> {
        checkpoint()?;
        if changes.len() > maximum_changes {
            return Err(BurnupError::Budget("stochastic activation forcing changes"));
        }
        changes
            .len()
            .checked_add(1)
            .and_then(|frames| frames.checked_mul(model.reactions().len()))
            .ok_or(BurnupError::Budget("stochastic activation rate storage"))?;
        model.admit_limits(limits)?;
        let mut outgoing = vec![Vec::new(); model.nuclides().len()];
        let mut branches = vec![Accumulator::new(0.0); model.nuclides().len()];
        for (index, reaction) in model.reactions().iter().enumerate() {
            checkpoint()?;
            match &reaction.rate {
                BurnupRate::Decay { branch_fraction } => {
                    add(&mut branches[reaction.parent], *branch_fraction)?;
                }
                BurnupRate::Capture { .. } | BurnupRate::TabulatedCapture { .. } => {}
                BurnupRate::NeutronProduction { .. }
                | BurnupRate::TabulatedNeutronProduction { .. } => {
                    return Err(BurnupError::Invalid(
                        "stochastic activation does not support neutron-producing reaction accounts",
                    ));
                }
                BurnupRate::Fission { .. } | BurnupRate::TabulatedFission { .. } => {
                    return Err(BurnupError::Invalid(
                        "stochastic fission requires explicit joint product outcomes",
                    ));
                }
            }
            if reaction.products.len() != 1 || reaction.products[0].yield_per_event != 1.0 {
                return Err(BurnupError::Invalid(
                    "stochastic activation requires one retained unit-yield product per channel",
                ));
            }
            outgoing[reaction.parent].push(index);
        }
        for (nuclide, branches) in model.nuclides().iter().zip(branches) {
            checkpoint()?;
            if nuclide.law.decay_constant().value() > 0.0 && branches.value() != 1.0 {
                return Err(BurnupError::Invalid(
                    "stochastic decay branches must sum to one without normalization",
                ));
            }
        }
        let initial_rates = model.rate_constants(
            &initial_flux,
            initial_temperature_k,
            limits,
            &mut checkpoint,
        )?;
        let mut changed_rates = Vec::with_capacity(changes.len());
        let mut previous = 0.0;
        for change in &changes {
            checkpoint()?;
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(BurnupError::Invalid(
                    "activation changes need strictly increasing positive times",
                ));
            }
            changed_rates.push(model.rate_constants(
                &change.flux,
                change.temperature_k,
                limits,
                &mut checkpoint,
            )?);
            previous = change.time_s;
        }
        checkpoint()?;
        Ok(Self {
            model,
            limits,
            initial_flux,
            initial_temperature_k,
            initial_rates,
            changes,
            changed_rates,
            outgoing,
        })
    }

    #[must_use]
    pub fn model(&self) -> &Arc<BurnupModel> {
        &self.model
    }

    #[must_use]
    pub fn limits(&self) -> BurnupLimits {
        self.limits
    }

    #[must_use]
    pub fn changes(&self) -> &[StochasticActivationChange] {
        &self.changes
    }

    pub fn flux_at(&self, time_s: f64) -> Result<&GroupFlux, BurnupError> {
        let position = self.position_at(time_s)?;
        Ok(if position == 0 {
            &self.initial_flux
        } else {
            &self.changes[position - 1].flux
        })
    }

    pub fn temperature_at(&self, time_s: f64) -> Result<f64, BurnupError> {
        let position = self.position_at(time_s)?;
        Ok(if position == 0 {
            self.initial_temperature_k
        } else {
            self.changes[position - 1].temperature_k
        })
    }

    /// Original reaction order; this borrows the prepared frame without collapse.
    pub fn rate_constants_at(&self, time_s: f64) -> Result<&[f64], BurnupError> {
        Ok(self.rates(self.position_at(time_s)?))
    }

    #[must_use]
    pub fn next_after(&self, time_s: f64) -> Option<f64> {
        self.changes
            .get(self.position_at(time_s).ok()?)
            .map(|change| change.time_s)
    }

    pub(super) fn position_at(&self, time_s: f64) -> Result<usize, BurnupError> {
        if !nonnegative(time_s) {
            return Err(BurnupError::Invalid(
                "activation time must be finite and nonnegative",
            ));
        }
        Ok(self
            .changes
            .partition_point(|change| change.time_s <= time_s))
    }

    pub(super) fn rates(&self, position: usize) -> &[f64] {
        if position == 0 {
            &self.initial_rates
        } else {
            &self.changed_rates[position - 1]
        }
    }

    pub(super) fn outgoing(&self, node: usize) -> &[usize] {
        &self.outgoing[node]
    }
}
