//! Replicated randomized QMC for uniform and empirical Gaussian-copula marginals.
//!
//! Shares marginal admission, the exact normal-CDF transform and physical
//! checkpoint binding with `GaussianCopulaExecution`. `QmcExecution` continues
//! to own Sobol addressing, latent PSD admission, replicate statistics, work
//! limits, interruption and durable replay. No separate RNG or inverse CDF.
//!
//! Dependence is within each physical input vector. Between-replicate standard
//! error remains descriptive: dependent points inside a net cannot be treated
//! as iid Bernoulli trials or used in the Monte Carlo stopping policy.

use core::fmt::Display;
use fs_blake3::ContentHash;
use crate::{EmpiricalMarginal, ParameterUncertainty, QmcConfig, QmcExecution, QmcReport, UqCheckpointError, UqPlan};
use crate::product_copula::{CopulaMarginals, EvaluationError, admit_marginals};

/// Resumable replicated Sobol sampling of bounded, dependent physical inputs.
#[derive(Debug, Clone)]
pub struct GaussianCopulaQmcExecution {
    pub(super) marginals: CopulaMarginals,
    execution: QmcExecution,
}

impl GaussianCopulaQmcExecution {
    /// Admit an explicit QMC marginal plan, latent-normal matrix and net layout.
    ///
    /// The marginal plan must declare `Independent`; the separately supplied
    /// copula supplies its complete joint law. Other joint laws are not erased.
    /// The existing QMC owner admits 1..=10 dimensions and complete replicates.
    /// Its finite-grid/normal-quantile limitations still apply; no physical work
    /// or sample is generated here. Invalid supports, PSD or budgets refuse.
    pub fn new(marginal_plan: &UqPlan, latent_correlation: &[Vec<f64>], config: QmcConfig)
        -> Result<Self, &'static str>
    {
        let (marginals, latent) = admit_marginals(marginal_plan, latent_correlation)?;
        let execution = QmcExecution::new(&latent, config)?;
        Ok(Self { marginals, execution })
    }

    // Read-only access for the shared physical mean-control implementation.
    // Never expose the latent sampler as a mutable physical execution.
    pub(crate) fn control_source(&self) -> &QmcExecution { &self.execution }

    /// Physical marginal supports and units in declaration order.
    #[must_use]
    pub fn marginals(&self) -> &[ParameterUncertainty] { &self.marginals.parameters }

    /// Empirical declarations in evaluator order, including source counts and
    /// observed ranges. Draws are confined to the observed values.
    pub fn empirical_marginals(&self) -> impl Iterator<Item = (&ParameterUncertainty, &EmpiricalMarginal)> {
        self.marginals.empirical()
    }

    /// Physical QoI observations, including an unfinished replicate's prefix.
    #[must_use]
    pub fn observations(&self) -> &[f64] { self.execution.observations() }

    /// Completed or terminally refused evaluations; interruptions do not count.
    #[must_use]
    pub fn evaluations_attempted(&self) -> usize { self.execution.evaluations_attempted() }

    /// Physical-QoI statistics from complete independent replicates only.
    #[must_use]
    pub fn report(&self) -> QmcReport { self.execution.report() }

    /// Perform at most `allowance` additional physical evaluations.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluator: F) -> QmcReport
    where F: FnMut(&[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |x| evaluator(x).map(Some))
    }

    /// `None` means unfinished work and retries the exact same physical vector
    /// after resumption. Errors and nonfinite QoIs are terminal, never redrawn.
    /// Completed observations stay paid work; zero allowance is a no-op.
    pub fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, cancelled: C, mut evaluator: F,
    ) -> QmcReport
    where F: FnMut(&[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        let marginals = &self.marginals;
        self.execution.advance_interruptible(allowance, cancelled, |latent| {
            let physical = marginals.transform(latent).map_err(EvaluationError::Transform)?;
            evaluator(&physical).map_err(EvaluationError::Producer)
        })
    }

    /// Persist observations with the latent matrix, complete net layout,
    /// physical supports/units and caller-supplied physical model identity.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        self.execution.checkpoint(self.marginals.bound_model(model))
    }

    /// Re-admit and restore without replaying any physical prefix.
    /// A changed law, layout, support, unit or model refuses. Hash integrity
    /// does not authenticate a producer or establish physical validity.
    pub fn restore(marginal_plan: &UqPlan, latent_correlation: &[Vec<f64>], config: QmcConfig,
        model: ContentHash, bytes: &[u8]) -> Result<Self, UqCheckpointError>
    {
        let (marginals, latent) = admit_marginals(marginal_plan, latent_correlation)
            .map_err(UqCheckpointError::InvalidPlan)?;
        let execution = QmcExecution::restore(&latent, config, marginals.bound_model(model), bytes)?;
        Ok(Self { marginals, execution })
    }
}

#[cfg(test)]
mod tests;
