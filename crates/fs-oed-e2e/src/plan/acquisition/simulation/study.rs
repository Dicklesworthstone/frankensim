//! Bounded stochastic campaign execution through the real ask/acquire/update API.

use fs_blake3::ContentHash;
use fs_exec::{AdmittedBudget, Cx};

use super::GaussianMeasurementSimulator;
use crate::{MAX_CAMPAIGN_SENSORS, ObjectiveValue, OedError};
use crate::plan::acquisition::AcquisitionError;
use crate::plan::measurement::{
    MeasurementAlternative, MeasurementDecision, MeasurementError, MeasurementStop, MeasurementStudy,
};

/// A finite segment either applied readings or encountered an actual study stop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SimulationProgress {
    applied: usize,
    stop: Option<MeasurementStop>,
}
impl SimulationProgress {
    /// Number of newly committed measurements in this successful call.
    #[must_use]
    pub const fn applied(&self) -> usize { self.applied }

    /// `None` means the segment limit was reached without another recommendation.
    /// It does not mean the study's decision threshold was achieved.
    #[must_use]
    pub const fn stop(&self) -> Option<MeasurementStop> { self.stop }
}

/// A complete fixed-truth simulation around an otherwise outcome-free study.
///
/// The planner receives only priors and declared instrument noise. Every update
/// uses a generated observation through the same acquisition and assimilation
/// operations used with an external provider. Splitting execution into segments
/// retains the study and addresses the next normal by the complete new request.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulatedStudy {
    study: MeasurementStudy,
    provider: GaussianMeasurementSimulator,
}
impl SimulatedStudy {
    /// Require exactly one schema-compatible truth per alternative before any
    /// samples can be generated. Unlike a partial instrument provider, a complete
    /// simulation cannot discover a missing truth halfway through its campaign.
    pub fn new(
        alternatives: Vec<MeasurementAlternative>,
        threshold: ObjectiveValue,
        max_readings: usize,
        provider: GaussianMeasurementSimulator,
        cx: &Cx<'_>,
    ) -> Result<Self, MeasurementError> {
        if alternatives.len() != provider.truths.len() {
            return Err(MeasurementError::Invalid("simulation must cover the complete alternative menu"));
        }
        let study = MeasurementStudy::new(alternatives, threshold, max_readings, cx)?;
        let mut budget = AdmittedBudget::admit_ambient(cx, 2_048 * provider.truths.len() as u64)
            .map_err(OedError::BudgetRefused)?;
        for truth in &provider.truths {
            budget.checkpoint("simulation menu matching", cx).map_err(OedError::BudgetRefused)?;
            budget.charge_cost("simulation menu matching", 2_048).map_err(OedError::BudgetRefused)?;
            let prior = study.posterior_mean(&truth.alternative)
                .ok_or(MeasurementError::Invalid("simulation truth names another alternative"))?;
            if prior.spec() != truth.value.spec() {
                return Err(MeasurementError::Invalid("simulation truth and prior quantity schemas differ"));
            }
        }
        budget.checkpoint("simulation study publication", cx).map_err(OedError::BudgetRefused)?;
        Ok(Self { study, provider })
    }

    /// Accepted state and evidence; synthetic truth is not installed in it.
    #[must_use]
    pub const fn study(&self) -> &MeasurementStudy { &self.study }

    /// Synthetic fixture, explicit seed and sampler identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash { self.provider.identity() }

    /// Apply at most `max_new_readings` actual stochastic observations.
    ///
    /// A failure retains the accepted prefix. The concrete stateless simulator
    /// can regenerate an acquired-but-unapplied sample exactly on retry. This is
    /// NOT a retry policy for hardware providers: their acquisitions must retain
    /// `AcquisitionError::After` data instead of being repeated.
    ///
    /// Segment overhead has one cumulative accountant. Recommendation, synthetic
    /// acquisition and numerical update retain their own existing compute scopes;
    /// this is not a claim of one aggregate campaign-lifetime cost/poll quota.
    pub fn advance(&mut self, max_new_readings: usize, cx: &Cx<'_>) -> Result<SimulationProgress, MeasurementError> {
        if max_new_readings > MAX_CAMPAIGN_SENSORS {
            return Err(OedError::TooManySensors { count: max_new_readings, max: MAX_CAMPAIGN_SENSORS }.into());
        }
        let mut budget = AdmittedBudget::admit_ambient(cx, max_new_readings as u64 + 1)
            .map_err(OedError::BudgetRefused)?;
        budget.checkpoint("simulation segment", cx).map_err(OedError::BudgetRefused)?;
        let mut applied = 0;
        for _ in 0..max_new_readings {
            budget.checkpoint("simulation next request", cx).map_err(OedError::BudgetRefused)?;
            budget.charge_cost("simulation next request", 1).map_err(OedError::BudgetRefused)?;
            let ask = match self.study.recommend(cx)? {
                MeasurementDecision::Stop(stop) => return Ok(SimulationProgress { applied, stop: Some(stop) }),
                MeasurementDecision::Acquire(request) => request,
            };
            let acquired = self.study.acquire(&ask, &mut self.provider, cx).map_err(|error| match error {
                AcquisitionError::Before(source)
                | AcquisitionError::Provider { source, .. }
                | AcquisitionError::After { source, .. } => source,
            })?;
            acquired.apply(&mut self.study, cx)?;
            applied += 1;
        }
        Ok(SimulationProgress { applied, stop: None })
    }

    /// Synthetic realized opportunity loss of the current mean-based decision.
    /// This diagnostic reads truth AFTER selection; it is never used by planning,
    /// stopping, or updating and is not an estimator of real-world decision risk.
    pub fn realized_loss(&self) -> Result<ObjectiveValue, MeasurementError> {
        let selected = self.provider.truths.iter()
            .find(|truth| truth.alternative == self.study.best_alternative())
            .expect("complete immutable simulation catalog");
        let best = self.provider.truths.iter().map(|truth| truth.value.value())
            .min_by(f64::total_cmp).expect("nonempty simulation catalog");
        let loss = selected.value.value() - best;
        if !loss.is_finite() {
            return Err(MeasurementError::Invalid("synthetic opportunity loss exceeds finite range"));
        }
        Ok(ObjectiveValue::from_raw(loss, selected.value.spec().decision_spec()))
    }
}
