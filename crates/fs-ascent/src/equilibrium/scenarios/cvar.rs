//! Empirical CVaR through the existing scenario physics and SQP engine.
//!
//! For s equally weighted, explicitly supplied realizations, minimize
//! eta + sum(z_i)/(s*(1-alpha)), with J_i(x)-eta-z_i <= 0 and -z_i <= 0.
//! Every original physical constraint still applies to EVERY realization.
//! This is the finite Rockafellar--Uryasev formulation, not a differentiable
//! approximation to sorting or a claim about unseen events. No risk statistic,
//! sampler, contact solver or optimizer is reimplemented here.
use super::*;

#[derive(Clone, Copy)]
pub(super) struct EmpiricalTail {
    alpha: f64,
    coefficient: f64,
}

impl ScenarioProblem<'_> {
    /// Select equally weighted empirical CVaR instead of the default maximum.
    /// The caller explicitly assigns mass 1/s to every declared scenario; this
    /// does not infer probabilities from tolerance endpoints or source labels.
    ///
    /// The dense cap now includes n+1+s decisions, 2*n common box faces,
    /// 2*s tail rows and s*c physical rows: 3*n+1+s*(3+c). Nothing is dropped
    /// when the tail mass s*(1-alpha) is fractional or smaller than one.
    /// Configuration consumes the problem before any study can cache it.
    pub fn with_cvar(mut self, alpha: f64) -> Result<Self, ScenarioError> {
        if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 {
            return Err(ScenarioError::Invalid("empirical CVaR alpha must be finite and strictly between zero and one"));
        }
        let count = self.scenarios.len();
        let coefficient = 1.0 / ((count as f64) * (1.0 - alpha));
        if !coefficient.is_finite() || coefficient <= 0.0 {
            return Err(ScenarioError::Invalid("empirical tail weight is not representable"));
        }
        let dimension = count.checked_mul(3 + self.problem.constraints().len())
            .and_then(|rows| self.lower.len().checked_mul(3)
                .and_then(|n| n.checked_add(1)).and_then(|n| n.checked_add(rows)))
            .ok_or(ScenarioError::Invalid("CVaR KKT dimension overflow"))?;
        if dimension > self.maximum_kkt_dimension {
            return Err(ScenarioError::Invalid("all CVaR slack, tail, box and physical rows must fit the dense KKT cap"));
        }
        self.tail = Some(EmpiricalTail { alpha, coefficient });
        Ok(self)
    }

    /// None means the original worst-case objective. Some(alpha) explicitly
    /// selects the empirical, equally weighted upper tail of the supplied set.
    #[must_use]
    pub fn cvar_alpha(&self) -> Option<f64> { self.tail.map(|tail| tail.alpha) }

    /// Recompute the Rockafellar--Uryasev score at the supplied threshold from
    /// actual re-solved losses, without trusting the optimizer's excess slacks.
    /// This is an UPPER BOUND on empirical CVaR for the given nominal design;
    /// it need not be the minimizing threshold after a budget or stalled run.
    /// It is not the exact order-statistic risk report owned by fs-robust.
    pub fn cvar_upper_bound(&self, evaluation: &ScenarioEvaluation, threshold: f64)
        -> Result<f64, ScenarioError>
    {
        let tail = self.tail.ok_or(ScenarioError::Invalid("CVaR score requires explicit CVaR configuration"))?;
        if !threshold.is_finite() || evaluation.scenarios.len() != self.scenarios.len() {
            return Err(ScenarioError::Invalid("CVaR score needs a finite threshold and complete scenario losses"));
        }
        let mut score = threshold;
        for scenario in &evaluation.scenarios {
            let excess = scenario.value - threshold;
            if !scenario.value.is_finite() || !excess.is_finite() {
                return Err(ScenarioError::Invalid("CVaR loss or excess is not finite"));
            }
            score += tail.coefficient * excess.max(0.0);
            if !score.is_finite() { return Err(ScenarioError::Invalid("CVaR score is not finite")); }
        }
        Ok(score)
    }

    pub(super) fn cvar_sample(&self, point: &[f64], evaluation: &ScenarioEvaluation,
        tail: EmpiricalTail) -> SqpSample
    {
        let n = self.lower.len();
        let count = self.scenarios.len();
        let width = n + 1 + count;
        let mut gradient = vec![0.0; width];
        gradient[n] = 1.0;
        gradient[n+1..].fill(tail.coefficient);
        let mut objective = point[n];
        for &slack in &point[n+1..] { objective += tail.coefficient * slack; }
        let mut ci = Vec::new(); let mut ji = Vec::new();
        let mut ce = Vec::new(); let mut je = Vec::new();
        for i in 0..n {
            ci.push(self.lower[i]-point[i]); ci.push(point[i]-self.upper[i]);
            let mut row = vec![0.0; width]; row[i] = -1.0; ji.extend_from_slice(&row);
            row[i] = 1.0; ji.extend_from_slice(&row);
        }
        for (i, physical) in evaluation.scenarios.iter().enumerate() {
            let slack = n + 1 + i;
            ci.push(physical.value - point[n] - point[slack]);
            let mut row = vec![0.0; width];
            row[..n].copy_from_slice(&physical.gradient); row[n] = -1.0; row[slack] = -1.0;
            ji.extend_from_slice(&row);
            ci.push(-point[slack]);
            row.fill(0.0); row[slack] = -1.0; ji.extend_from_slice(&row);
            // Reuse the original sign/scaling authority. Tail selection NEVER
            // exempts a realization from an equality or physical safety limit.
            let original = physical_sample(self.problem, physical);
            ci.extend_from_slice(&original.ci[2*n..]);
            for source in original.ji[2*n*n..].chunks_exact(n) {
                row.fill(0.0); row[..n].copy_from_slice(source); ji.extend_from_slice(&row);
            }
            ce.extend_from_slice(&original.ce);
            for source in original.je.chunks_exact(n) {
                row.fill(0.0); row[..n].copy_from_slice(source); je.extend_from_slice(&row);
            }
        }
        // SqpSample's original finite/shape admission also checks objective
        // overflow and malformed trials before any accepted state changes.
        SqpSample { f: objective, gradient, ce, ci, je, ji }
    }
}

#[cfg(test)]
mod tests;
