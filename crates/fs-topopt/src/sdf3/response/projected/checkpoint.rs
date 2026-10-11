//! Recover the existing response optimizer, never a different numerical method.
use super::*;
use fs_ascent::projected_al::ProjectedAlCheckpoint;

/// Complete accepted optimizer state and least-objective feasible endpoint.
///
/// Bind this value to the original geometry, support, independent external
/// loads, prescribed-motion/reaction laws, targets, SIMP/filter model and all
/// numerical policies before restoration. Cached derivatives are checked by
/// fresh physics; displacement fields and factors are deliberately NOT stored.
/// History is retained evidence, not independently re-solved old iterations.
/// An authenticated consumer must retain the original cumulative budgets too.
#[derive(Debug, Clone)]
pub struct ProjectedResponseCheckpoint3 {
    /// Existing AL checkpoint includes multipliers, spectral steps and work.
    pub optimizer: ProjectedAlCheckpoint,
    /// Only a density is needed to reconstruct the independent incumbent.
    pub best_feasible_density: Option<Vec<f64>>,
    /// Original accepted physical observations, including iteration zero.
    pub history: Vec<ProjectedResponseIteration3>,
    /// Endpoint re-evaluations already included in optimizer.work.evaluations.
    pub restoration_evaluations: usize,
}
impl ProjectedResponseCheckpoint3 {
    /// Re-solve the accepted endpoint; a distinct feasible incumbent costs one
    /// additional full experiment family. This never replays old search steps.
    #[must_use]
    pub fn restoration_cost(&self) -> usize {
        1 + usize::from(self.best_feasible_density.as_ref()
            .is_some_and(|rho| rho != &self.optimizer.point))
    }
}

impl<'a, 'callback, O: AdaptiveSdf3Elasticity> ProjectedResponseStudy3<'a, 'callback, O> {
    /// Snapshot only matching, completely accepted optimizer/physics state.
    /// Rejected trials do not become incumbents, but their work stays charged.
    #[must_use]
    pub fn checkpoint(&self) -> ProjectedResponseCheckpoint3 {
        ProjectedResponseCheckpoint3 {
            optimizer: self.state.checkpoint(),
            best_feasible_density: self.best_feasible.as_ref().map(|e| e.rho.clone()),
            history: self.history.clone(),
            restoration_evaluations: self.restoration_evaluations,
        }
    }

    /// Restore a displacement-response fit with the SAME reference experiment.
    /// Both primal and aggregate adjoint are recomputed, including the CURRENT
    /// density's Nitsche lifting and its material derivative. Incoming scales
    /// survive every returned failure, including a rejected incumbent check.
    ///
    /// Reserve all endpoint evaluations in the original optimizer allowance
    /// before physics. Actual solves, including failed restoration, spend the
    /// supplied SolveControl; the consumer must retain that consumption after
    /// errors and cannot replace it with a renewed lifetime allowance.
    pub fn restore(
        study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseCase3<'a>],
        checkpoint: ProjectedResponseCheckpoint3, options: ProjectedResponseOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<ResponseError3>> {
        Self::restore_impl(study, cases, None, checkpoint, options, control)
    }

    /// The same restoration for mixed displacement/reaction experiments.
    /// Reactions retain their direct material derivative and original force or
    /// moment mode; they are not reconstructed as fixed displacement weights.
    pub fn restore_with_reactions(
        study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseCase3<'a>],
        reactions: &'a [&'a [ReactionTarget3<'a>]],
        checkpoint: ProjectedResponseCheckpoint3, options: ProjectedResponseOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<ResponseError3>> {
        Self::restore_impl(study, cases, Some(reactions), checkpoint, options, control)
    }

    fn restore_impl(
        study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseCase3<'a>],
        reactions: Option<&'a [&'a [ReactionTarget3<'a>]]>,
        checkpoint: ProjectedResponseCheckpoint3, options: ProjectedResponseOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<ResponseError3>> {
        control.checkpoint("response-projected-restore")
            .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
        let n = study.cells();
        options.optimizer.validate(n)?;
        if !options.volume_cap.is_finite() || !(0.0 < options.volume_cap && options.volume_cap <= 1.0)
            || !options.density_floor.is_finite() || !(0.0 < options.density_floor && options.density_floor < 1.0)
            || !options.objective_scale.is_finite() || options.objective_scale <= 0.0 {
            return Err(ProjectedAlError::Invalid("invalid response restoration policy"));
        }
        let valid_density = |rho: &[f64]| rho.len() == n && rho.iter()
            .all(|r| r.is_finite() && (options.density_floor..=1.0).contains(r));
        if checkpoint.best_feasible_density.as_ref().is_some_and(|rho| !valid_density(rho)) {
            return Err(ProjectedAlError::Invalid("invalid response incumbent density"));
        }
        let cost = checkpoint.restoration_cost();
        let restoration_evaluations = checkpoint.restoration_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("response restoration counter overflow"))?;
        let mut optimizer = checkpoint.optimizer;
        let previous_evaluations = optimizer.work.evaluations;
        optimizer.work.evaluations = previous_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("response restoration evaluation overflow"))?;
        // Admission checks the original callback cap BEFORE any endpoint solve.
        let restored = ProjectedAlState::try_restore(
            optimizer, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
        )?;
        let history = checkpoint.history;
        if checkpoint.restoration_evaluations >= previous_evaluations
            || restored.work().iterations.checked_add(1) != Some(history.len())
            || history.iter().enumerate().any(|(i, r)| {
                r.iteration != i || !r.objective.is_finite() || r.objective < 0.0
                    || !r.volume_fraction.is_finite() || !(0.0..=1.0).contains(&r.volume_fraction)
                    || !r.constraint_violation.is_finite()
                    || r.constraint_violation != (r.volume_fraction - options.volume_cap).max(0.0)
            }) {
            return Err(ProjectedAlError::Invalid("invalid response checkpoint history or work"));
        }
        let best_row = history.iter()
            .filter(|r| r.constraint_violation <= options.optimizer.tolerance)
            .fold(None, |best: Option<&ProjectedResponseIteration3>, r| {
                if best.is_none_or(|b| r.objective < b.objective) { Some(r) } else { best }
            });
        if best_row.is_some() != checkpoint.best_feasible_density.is_some() {
            return Err(ProjectedAlError::Invalid("response checkpoint lost its feasible incumbent"));
        }
        let previous_scales = study.operator.scales().to_vec();
        let mut result = Self::new_impl(study, cases, reactions, restored.point(), options, control)?;
        let rebuilt = (|| {
            if result.state.sample() != restored.sample()
                || history.last() != Some(&row(&result.accepted, restored.work().iterations, options.volume_cap)) {
                return Err(ProjectedAlError::Invalid("restored response endpoint changed its physics"));
            }
            let best = match checkpoint.best_feasible_density {
                None => None,
                Some(rho) if rho == result.accepted.rho => Some(result.accepted.clone()),
                Some(rho) => {
                    let mut evaluated = None;
                    sample(result.study, cases, reactions, &rho, options, &mut evaluated, result.control)
                        .map_err(ProjectedAlError::Evaluation)?;
                    evaluated
                }
            };
            if let (Some(evaluated), Some(expected)) = (&best, best_row) {
                if row(evaluated, expected.iteration, options.volume_cap) != *expected {
                    return Err(ProjectedAlError::Invalid("restored response incumbent changed its physics"));
                }
            }
            result.control.checkpoint("response-projected-restore-publish")
                .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
            Ok(best)
        })();
        let best = match rebuilt {
            Ok(best) => best,
            Err(error) => {
                result.study.operator.set_scales(&previous_scales).expect("incoming scales remain admitted");
                return Err(error);
            }
        };
        result.state = restored;
        result.best_feasible = best;
        result.history = history;
        result.restoration_evaluations = restoration_evaluations;
        Ok(result)
    }
}
