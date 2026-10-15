//! Bounded globalization of the existing Newton equation, not another stepper.
use super::{StepWorkspace, PhsError, PreparedStepError, PortHamiltonian,
    dimensions, residual_into, dissipation::Dissipation};

impl StepWorkspace {
    /// Opt into residual-decreasing Newton updates. Zero retains the original
    /// undamped traversal, including its callback count and arithmetic order.
    /// Otherwise try alpha=1,1/2,...,2^-halvings (at most 17 residual evaluations
    /// per Newton update). Iteration limits and both convergence thresholds are
    /// unchanged. A failed search uses the original stagnation admission.
    ///
    /// Full analytic/finite-difference and condensed/dense solves all use the
    /// same complete Gonzalez residual. Scratch is reused; callbacks must remain
    /// read-only and allocation-free for the hot path to avoid allocation.
    /// Nonfinite values, invalid constitutive callbacks and cancellation still
    /// refuse; their errors are not swallowed to obtain a smaller trial.
    /// No temporal error estimate, unique-root or real-time claim follows.
    pub fn set_backtracking_limit(&mut self, halvings: u8) -> Result<(), PhsError> {
        if halvings > 16 { return Err(dimensions("Newton backtracking needs 0..16 halvings")); }
        self.backtracking_halvings = halvings;
        Ok(())
    }

    /// (Extra residual evaluations, reduced updates accepted) in the latest
    /// solver call, including refused attempts. Counters are not wall-clock time.
    pub fn backtracking_work(&self) -> (usize, usize) { self.backtracking_work }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn backtracked_update<F>(
        &mut self, sys: &PortHamiltonian, x0: &[f64], dt: f64,
        current_norm: f64, tolerance: f64, dissipation: Option<Dissipation<'_>>,
        poll: &mut F,
    ) -> Result<bool, PreparedStepError>
    where F: FnMut() -> Result<(), PreparedStepError> {
        let mut alpha = 1.0;
        for halving in 0..=self.backtracking_halvings {
            poll()?;
            for i in 0..self.n { self.trial[i] = self.x[i] - alpha*self.delta[i]; }
            self.backtracking_work.0 = self.backtracking_work.0.saturating_add(1);
            let trial_norm = residual_into(sys, x0, &self.trial, &self.forcing, dt,
                &self.flow, &mut self.midpoint, &mut self.effort, &mut self.plus,
                dissipation, &mut self.nonlinear_loss)?;
            // Infinity-norm Armijo decrease; same norm/units as Newton admission.
            // Already converged finite candidates need not prove more decrease.
            if trial_norm <= (1.0 - 1.0e-4*alpha)*current_norm || trial_norm <= tolerance {
                poll()?;
                self.x.copy_from_slice(&self.trial);
                if halving > 0 {
                    self.backtracking_work.1 = self.backtracking_work.1.saturating_add(1);
                }
                return Ok(true);
            }
            alpha *= 0.5;
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests;
