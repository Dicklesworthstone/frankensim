//! Sample-scoped stress constraints on the existing projected load family.
//!
//! No stress adjoint or alternative design velocity is introduced. The baseline
//! must satisfy the declared limit; every accepted update preserves it. Tight
//! limits can exhaust the compliance-generated candidate family without a step.

use super::*;
use crate::robust_stress::sample_nodal_solution_controlled;

/// The same gate admits the baseline and every candidate. Missing measurements
/// cannot satisfy an installed constraint, and zero-weight cases remain present.
pub(super) fn require_feasible(
    limit: Option<SampledStressLimit>,
    stress: Option<&RobustSampledStressEvaluation>,
) -> Result<(), CutFemError> {
    let Some(limit) = limit else { return Ok(()) };
    let Some(stress) = stress else {
        return Err(invalid("sampled stress constraint has no complete load-family evaluation"));
    };
    if !(stress.worst_sampled_von_mises.is_finite()
        && stress.worst_sampled_von_mises >= 0.0
        && stress.worst_sampled_von_mises <= limit.admitted_max())
    {
        return Err(invalid(format!(
            "sampled stress limit exceeded in case {}: measured {} > admitted {}",
            stress.worst_stress_case, stress.worst_sampled_von_mises, limit.admitted_max(),
        )));
    }
    Ok(())
}

impl MultiLoadProjectedOptimizer {
    /// Install a fixed sampled plane-strain von Mises constraint on a newly
    /// created optimizer, before candidate solves or updates. Reuses every cached baseline
    /// displacement; no extra PDE solves are charged or hidden from the budget.
    ///
    /// The already area-feasible baseline must pass. This does NOT restore
    /// stress feasibility from an overstressed baseline. All declared load cases
    /// participate, including cases with zero compliance-objective weight.
    /// Constraint configuration cannot change after attempts or accepted updates.
    ///
    /// # Errors
    /// Refuses malformed/repeated/late limits, unavailable samples or a baseline
    /// above the admitted limit. The optimizer is returned only after full success.
    pub fn with_sampled_stress_limit(mut self, limit: SampledStressLimit) -> Result<Self, CutFemError> {
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        if self.stress_limit.is_some() || self.next_iteration != 0
            || self.solves_started != self.kernel.load_cases.len()
        {
            return Err(invalid("sampled stress limit must be installed once before candidate solves or updates"));
        }
        let stress = match self.sample_stress_controlled(&self.current, |_, _| {
            ControlFlow::<Infallible>::Continue(())
        })? {
            ControlFlow::Continue(stress) => stress,
            ControlFlow::Break(never) => match never {},
        };
        require_feasible(Some(limit), Some(&stress))?;
        self.stress_limit = Some(limit);
        self.baseline_stress = Some(stress.clone());
        self.current_stress = Some(stress);
        Ok(self)
    }

    /// Fixed constraint declaration; `None` means stress is NOT assessed.
    #[must_use]
    pub const fn stress_limit(&self) -> Option<SampledStressLimit> { self.stress_limit }

    /// Complete stress evidence for the feasible initial design.
    #[must_use]
    pub fn baseline_stress(&self) -> Option<&RobustSampledStressEvaluation> {
        self.baseline_stress.as_ref()
    }

    /// Complete stress evidence for the exact retained accepted geometry.
    #[must_use]
    pub fn current_stress(&self) -> Option<&RobustSampledStressEvaluation> {
        self.current_stress.as_ref()
    }

    pub(super) fn sample_stress_controlled<B>(
        &self,
        state: &MultiState,
        mut control: impl FnMut(usize, usize) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, RobustSampledStressEvaluation>, CutFemError> {
        let count = self.kernel.load_cases.len();
        if count == 0 || state.solutions.len() != count || state.compliances.len() != count {
            return Err(invalid("sampled stress requires the complete matching displacement family"));
        }
        let mut maxima = Vec::with_capacity(count);
        let mut locations = Vec::with_capacity(count);
        let mut counts = Vec::with_capacity(count);
        let mut worst_stress = f64::NEG_INFINITY;
        let mut worst_case = 0;
        let mut weighted_sum = 0.0;
        let mut worst_weighted = 0.0_f64;
        for (case, solution) in state.solutions.iter().enumerate() {
            let (maximum, location, samples) = match sample_nodal_solution_controlled(
                &self.kernel.grid, &state.phi, solution, self.kernel.lambda, self.kernel.mu,
                |cell| control(case, cell),
            )? {
                ControlFlow::Continue(samples) => samples,
                ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
            };
            let weighted = self.kernel.load_cases[case].weight() * state.compliances[case];
            weighted_sum += weighted;
            worst_weighted = worst_weighted.max(weighted);
            if maximum > worst_stress {
                worst_stress = maximum;
                worst_case = case;
            }
            maxima.push(maximum);
            locations.push(location);
            counts.push(samples);
        }
        if !(weighted_sum.is_finite() && worst_weighted.is_finite() && worst_stress.is_finite()) {
            return Err(invalid("sampled stress family produced a non-finite aggregate"));
        }
        Ok(ControlFlow::Continue(RobustSampledStressEvaluation {
            case_compliances: state.compliances.clone(),
            case_sampled_max_von_mises: maxima,
            case_max_locations: locations,
            case_sample_counts: counts,
            weighted_sum_compliance: weighted_sum,
            worst_weighted_compliance: worst_weighted,
            objective: state.objective,
            worst_sampled_von_mises: worst_stress,
            worst_stress_case: worst_case,
            volume: state.volume,
            snapshot: fnv(&state.phi),
        }))
    }
}

#[cfg(test)]
mod tests;
