//! Restore the shared load-space optimizer from accepted numerical state.
use super::*;
use fs_ascent::projected_al::ProjectedAlCheckpoint;

/// Accepted AL state, original gradient experiment and the least-material (or
/// least-compliance) feasible incumbent. No operator or displacement is trusted
/// across restoration. Consumers must bind the same physical source and policy.
#[derive(Debug, Clone)]
pub struct LoadSpaceDesignCheckpoint3 {
    pub optimizer: ProjectedAlCheckpoint,
    pub best_feasible_density: Option<Vec<f64>>,
    pub history: Vec<LoadSpaceDesignIteration3>,
    pub gradient_check: LoadSpaceGradientCheck3,
    /// Fresh endpoint evaluations already charged to optimizer work.
    pub restoration_evaluations: usize,
}
impl LoadSpaceDesignCheckpoint3 {
    /// One accepted endpoint and, only when distinct, one feasible incumbent.
    #[must_use]
    pub fn restoration_cost(&self) -> usize {
        1 + usize::from(self.best_feasible_density.as_ref()
            .is_some_and(|rho| rho != &self.optimizer.point))
    }
}

impl<'a, 'callback, O: Sdf3Elasticity> LoadSpaceDesignStudy3<'a, 'callback, O> {
    #[must_use]
    pub fn checkpoint(&self) -> LoadSpaceDesignCheckpoint3 {
        LoadSpaceDesignCheckpoint3 {
            optimizer: self.state.checkpoint(),
            best_feasible_density: self.best.as_ref().map(|p| p.rho.clone()),
            history: self.history.clone(), gradient_check: self.check.clone(),
            restoration_evaluations: self.restoration_evaluations,
        }
    }

    /// Restore the original volume-constrained formulation without replaying
    /// earlier search steps or the initial finite-difference experiment.
    pub fn restore(study: &'a mut CutDensityStudy3<O>, basis: &'a [&'a [f64]],
        checkpoint: LoadSpaceDesignCheckpoint3, policy: LoadSpace3,
        options: EnvelopeDesignOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        Self::restore_family(study, basis, checkpoint, policy, options, None, control)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn restore_family(study: &'a mut CutDensityStudy3<O>, basis: &'a [&'a [f64]],
        checkpoint: LoadSpaceDesignCheckpoint3, policy: LoadSpace3,
        options: EnvelopeDesignOptions3, cap: Option<limits::LoadSpaceCap3>,
        control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        admit(study, basis, &checkpoint.optimizer.point, &policy, options)?;
        let n = study.cells();
        let check = &checkpoint.gradient_check;
        let numerical = &check.numerical;
        let p = study.params();
        let q = numerical.params;
        if check.policy != policy || check.robustness != options.robustness
            || [p.e_min, p.penal, p.beta, p.eta] != [q.e_min, q.penal, q.beta, q.eta]
            || numerical.options.step != options.gradient_check.step
            || numerical.options.relative_tolerance != options.gradient_check.relative_tolerance
            || !numerical.passed() || numerical.probes.len() > 2
            || !numerical.baseline_compliance.is_finite() || numerical.baseline_compliance < 0.0
            || !(0.0..=1.0).contains(&numerical.baseline_volume_fraction)
            || numerical.probes.iter().any(|probe| probe.active_densities > n
                || [probe.compliance_analytic, probe.compliance_difference,
                    probe.volume_analytic, probe.volume_difference].iter().any(|v| !v.is_finite())
                || probe.compliance_relative_error < 0.0 || probe.volume_relative_error < 0.0)
            || (numerical.probes.len() == 2
                && numerical.probes[0].direction == numerical.probes[1].direction)
        {
            return Err(ProjectedAlError::Invalid("load-space checkpoint changed its gradient experiment"));
        }
        let cost = checkpoint.restoration_cost();
        let restoration_evaluations = checkpoint.restoration_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("load-space restoration counter overflow"))?;
        let mut optimizer = checkpoint.optimizer;
        let previous_evaluations = optimizer.work.evaluations;
        optimizer.work.evaluations = previous_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("load-space evaluation counter overflow"))?;
        let state = ProjectedAlState::try_restore(optimizer,
            &vec![options.density_floor; n], &vec![1.0; n], options.optimizer)?;
        let history = checkpoint.history;
        if checkpoint.restoration_evaluations >= previous_evaluations
            || history.len() != state.work().iterations + 1
            || history.iter().enumerate().any(|(i, r)| r.iteration != i
                || !r.maximum.is_finite() || r.maximum < 0.0
                || !r.smooth_maximum.is_finite() || r.smooth_maximum < r.maximum
                || !(0.0..=1.0).contains(&r.dilated_volume)
                || (r.feasible && if cap.is_some() {
                    r.maximum > 1.0 || r.smooth_maximum-1.0 > options.optimizer.tolerance
                } else { r.dilated_volume/options.volume_limit-1.0 > options.optimizer.tolerance }))
        {
            return Err(ProjectedAlError::Invalid("invalid load-space checkpoint history or work"));
        }
        let best_row = history.iter().filter(|r| r.feasible).fold(None,
            |best: Option<&LoadSpaceDesignIteration3>, r| {
                if best.is_none_or(|b| if cap.is_some() { r.dilated_volume < b.dilated_volume }
                    else { r.smooth_maximum < b.smooth_maximum }) { Some(r) } else { best }
            });
        if best_row.is_some() != checkpoint.best_feasible_density.is_some() {
            return Err(ProjectedAlError::Invalid("load-space checkpoint lost its feasible incumbent"));
        }
        if let Some(rho) = &checkpoint.best_feasible_density {
            admit(study, basis, rho, &policy, options)?;
        }
        // The existing evaluator restores incoming scales on success AND error.
        // Compare all objective/constraint gradients, not only a scalar response.
        control.checkpoint("sdf3-load-space-restore").map_err(ProjectedAlError::Evaluation)?;
        let mut accepted = None;
        let evaluated = sample(study, state.point(), basis, &policy, options, cap.is_some(),
            &mut accepted, control).map_err(ProjectedAlError::Evaluation)?;
        let accepted = accepted.ok_or(ProjectedAlError::Invalid("missing restored load-space endpoint"))?;
        if evaluated.as_ref() != Some(state.sample())
            || history.last() != Some(&row(&accepted, state.work().iterations, options, cap.as_ref())) {
            return Err(ProjectedAlError::Invalid("restored load-space point changed its physical response"));
        }
        let best = match checkpoint.best_feasible_density {
            None => None,
            Some(rho) if rho == accepted.rho => Some(accepted.clone()),
            Some(rho) => {
                let mut best = None;
                sample(study, &rho, basis, &policy, options, cap.is_some(), &mut best, control)
                    .map_err(ProjectedAlError::Evaluation)?;
                best
            }
        };
        if let (Some(best), Some(expected)) = (&best, best_row)
            && row(best, expected.iteration, options, cap.as_ref()) != *expected {
            return Err(ProjectedAlError::Invalid("restored load-space incumbent changed its physical response"));
        }
        control.checkpoint("sdf3-load-space-restored").map_err(ProjectedAlError::Evaluation)?;
        Ok(Self { study, basis, control, policy, options, cap, state, accepted, best,
            history, check: checkpoint.gradient_check, restoration_evaluations })
    }
}

#[cfg(test)]
mod tests;
