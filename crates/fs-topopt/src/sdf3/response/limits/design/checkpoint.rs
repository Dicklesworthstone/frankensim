//! Restore the original constrained optimizer; do not replay its old searches.
use super::*;
use fs_ascent::projected_al::ProjectedAlCheckpoint;

type Error = ProjectedAlError<ResponseError3>;
fn invalid(message: &'static str) -> Error { ProjectedAlError::Invalid(message) }
fn physical(error: impl Into<ResponseError3>) -> Error {
    ProjectedAlError::Evaluation(error.into())
}

/// Ordered physical-bound signature. Callback identity cannot be inferred from
/// function pointers; bind the complete source/geometry/executable externally.
/// Even a currently inactive limit's kind, sense, bound and scale are retained.
#[derive(Debug, Clone, PartialEq)]
pub struct MinimumMaterialLimit3 {
    pub reaction: bool,
    pub sense: ResponseLimitSense3,
    pub bound: f64,
    pub scale: f64,
}
fn signatures(cases: &[ResponseLimitCase3<'_>]) -> Vec<Vec<MinimumMaterialLimit3>> {
    cases.iter().map(|case| case.limits.iter().map(|limit| MinimumMaterialLimit3 {
        reaction: matches!(limit.response, BoundedResponse3::Reaction(_)),
        sense: limit.sense, bound: limit.bound, scale: limit.scale,
    }).collect()).collect()
}

/// A complete endpoint's compact observation witness. No displacement field,
/// adjoint field, quadrature or preconditioner is trusted from the checkpoint.
/// This is O(cells + limits), not O(cases * displacement degrees of freedom).
#[derive(Debug, Clone, PartialEq)]
pub struct MinimumMaterialEndpoint3 {
    pub rho: Vec<f64>,
    pub projected_rho: Vec<f64>,
    pub scales: Vec<f64>,
    pub volume_fraction: f64,
    pub mass_kg: Option<f64>,
    pub objective: f64,
    pub constraint: f64,
    pub worst_residual: f64,
    pub measurements: Vec<Vec<f64>>,
    pub residuals: Vec<Vec<f64>>,
}
impl MinimumMaterialEndpoint3 {
    fn from_evaluation(e: &ResponseLimitEvaluation3, options: MinimumMaterialOptions3) -> Self {
        Self { rho: e.rho.clone(), projected_rho: e.projected_rho.clone(), scales: e.scales.clone(),
            volume_fraction: e.volume_fraction, mass_kg: e.mass_kg,
            objective: objective_value(e.volume_fraction, e.mass_kg, options).expect("admitted endpoint objective"),
            constraint: e.constraint, worst_residual: e.worst_residual,
            measurements: e.measurements.clone(), residuals: e.residuals.clone() }
    }
    fn row(&self, iteration: usize) -> MinimumMaterialIteration3 {
        MinimumMaterialIteration3 { iteration, volume_fraction: self.volume_fraction,
            mass_kg: self.mass_kg, objective: self.objective,
            constraint: self.constraint, worst_residual: self.worst_residual }
    }
    fn validate(&self, n: usize, cases: &[ResponseLimitCase3<'_>], options: MinimumMaterialOptions3) -> Result<(), Error> {
        if self.rho.len() != n || self.projected_rho.len() != n || self.scales.len() != n
            || self.rho.iter().any(|r| !r.is_finite() || !(options.density_floor..=1.0).contains(r))
            || self.projected_rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || self.scales.iter().any(|s| !s.is_finite() || *s <= 0.0 || *s > 1.0)
            || !self.volume_fraction.is_finite() || !(0.0..=1.0).contains(&self.volume_fraction)
            || objective_value(self.volume_fraction, self.mass_kg, options) != Some(self.objective)
            || !self.constraint.is_finite() || !self.worst_residual.is_finite()
            || self.measurements.len() != cases.len() || self.residuals.len() != cases.len() {
            return Err(invalid("invalid minimum-material endpoint shape or values"));
        }
        for ((case, values), residuals) in cases.iter().zip(&self.measurements).zip(&self.residuals) {
            if values.len() != case.limits.len() || residuals.len() != case.limits.len()
                || values.iter().chain(residuals).any(|v| !v.is_finite()) {
                return Err(invalid("invalid minimum-material endpoint observations"));
            }
        }
        let worst = self.residuals.iter().flatten().copied().fold(f64::NEG_INFINITY, f64::max);
        if worst != self.worst_residual { return Err(invalid("checkpoint lost an original limit residual")); }
        Ok(())
    }
}

/// Caller-owned restart data for the SAME fixed-grid material-design problem.
/// Source identity, executable version, geometry, field ordering, preconditioner,
/// pure load/motion/observation laws, SIMP/filter and protected regions must be
/// bound externally, alongside original cumulative physical/geometry budgets.
/// Endpoint checks cannot prove that arbitrary callbacks are unchanged globally.
/// History is retained evidence, not independently re-solved historical points.
#[derive(Debug, Clone)]
pub struct MinimumMaterialCheckpoint3 {
    pub optimizer: ProjectedAlCheckpoint,
    /// Complete prepared body map and native source attribution, when attached.
    pub material_body_load_fingerprint: Option<fs_blake3::ContentHash>,
    pub options: MinimumMaterialOptions3,
    pub limits: Vec<Vec<MinimumMaterialLimit3>>,
    pub prescribed: Vec<bool>,
    /// Ordered native reference inventory binds a mass objective to every cell,
    /// including changes invisible in the current projected design. None for
    /// the geometric-volume objective; full source identity remains external.
    pub reference_cell_masses_kg: Option<Vec<f64>>,
    pub accepted: MinimumMaterialEndpoint3,
    /// Distinct least-objective feasible accepted point; may genuinely be absent.
    pub best_feasible: Option<MinimumMaterialEndpoint3>,
    pub history: Vec<MinimumMaterialIteration3>,
    /// Reserved endpoint evaluations already included in optimizer.work.
    /// Failed/cancelled restoration reservations stay charged on this value.
    pub restoration_evaluations: usize,
}
impl MinimumMaterialCheckpoint3 {
    /// One current endpoint plus one distinct feasible incumbent, never the
    /// complete optimization trajectory. Equal raw points share one solve.
    #[must_use]
    pub fn restoration_cost(&self) -> usize {
        1 + usize::from(self.best_feasible.as_ref().is_some_and(|best| best.rho != self.optimizer.point))
    }
}
fn same_policy(a: MinimumMaterialOptions3, b: MinimumMaterialOptions3) -> bool {
    let counts = |o: MinimumMaterialOptions3| [o.limits.max_cases, o.limits.max_limits,
        o.optimizer.max_dimension, o.optimizer.max_evaluations,
        o.optimizer.max_multiplier_updates, o.optimizer.max_backtracks];
    let floats = |o: MinimumMaterialOptions3| [o.limits.sharpness, o.density_floor, o.volume_scale,
        o.optimizer.tolerance, o.optimizer.initial_penalty, o.optimizer.max_penalty,
        o.optimizer.inner_tolerance, o.optimizer.min_spectral_step, o.optimizer.max_spectral_step].map(f64::to_bits);
    counts(a) == counts(b) && floats(a) == floats(b) && a.objective == b.objective
}
fn feasible_row(r: &MinimumMaterialIteration3, tolerance: f64) -> bool {
    r.constraint <= tolerance && r.worst_residual <= tolerance
}

impl<'a, 'callback, O: AdaptiveSdf3Elasticity> MinimumMaterialStudy3<'a, 'callback, O> {
    /// Snapshot accepted numerical state and individual-bound witnesses only.
    #[must_use]
    pub fn checkpoint(&self) -> MinimumMaterialCheckpoint3 {
        MinimumMaterialCheckpoint3 {
            optimizer: self.state.checkpoint(), options: self.options,
            material_body_load_fingerprint: self.study.material_body_load_fingerprint(),
            limits: signatures(self.cases), prescribed: self.cases.iter().map(|c| c.prescribed.is_some()).collect(),
            reference_cell_masses_kg: match self.options.objective {
                MinimumMaterialObjective3::VolumeFraction => None,
                MinimumMaterialObjective3::MassKg { .. } => Some(self.study.operator.reference_cell_masses_kg()
                    .expect("admitted source material inventory").to_vec()),
            },
            accepted: MinimumMaterialEndpoint3::from_evaluation(&self.accepted, self.options),
            best_feasible: self.best_feasible.as_ref().map(|e| MinimumMaterialEndpoint3::from_evaluation(e, self.options)),
            history: self.history.clone(), restoration_evaluations: self.restoration_evaluations,
        }
    }

    /// Rebuild the current endpoint and distinct feasible incumbent, including
    /// independent equilibria, original observations and complete derivatives.
    /// Preserve the original multiplier, penalty, spectral step and dual schedule.
    /// Exact same-executable values are required; no tolerance silently changes
    /// the physical bounds or the derivative sample that drives continuation.
    ///
    /// `checkpoint` is mutable so all endpoint reservations remain charged even
    /// after failed/cancelled restoration. Persist that updated value on failure;
    /// retrying an old copy would discard spent work. The shared SolveControl
    /// ALSO retains actual physical work, which the caller must persist together
    /// with the original remaining budgets. Reservations happen after structural
    /// admission but before any endpoint physics. An early refusal spends none.
    ///
    /// The incoming operator remains unchanged on EVERY returned error. Only
    /// after all checks and the final cancellation checkpoint are complete are
    /// the accepted scales installed. No stale fields or factors are imported.
    pub fn restore(study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseLimitCase3<'a>],
        checkpoint: &mut MinimumMaterialCheckpoint3, options: MinimumMaterialOptions3,
        control: &'a mut SolveControl<'callback>) -> Result<Self, Error> {
        if checkpoint.material_body_load_fingerprint != study.material_body_load_fingerprint() {
            return Err(invalid("minimum-material restart changed its material body load"));
        }
        control.checkpoint("minimum-material-restore-start").map_err(physical)?;
        let n = study.cells();
        options.optimizer.validate(n)?;
        if !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
            || !options.volume_scale.is_finite() || options.volume_scale <= 0.0 || !(1.0/options.volume_scale).is_finite()
            || !same_policy(checkpoint.options, options) {
            return Err(invalid("minimum-material restart policy changed"));
        }
        admit_objective(study, options)?;
        let inventory = match options.objective {
            MinimumMaterialObjective3::VolumeFraction => None,
            MinimumMaterialObjective3::MassKg { .. } => study.operator.reference_cell_masses_kg(),
        };
        if checkpoint.reference_cell_masses_kg.as_deref() != inventory {
            return Err(invalid("minimum-material restart changed its reference cell masses"));
        }
        admit_limits(study, cases, options.limits).map_err(ProjectedAlError::Evaluation)?;
        if checkpoint.limits != signatures(cases)
            || checkpoint.prescribed != cases.iter().map(|c| c.prescribed.is_some()).collect::<Vec<_>>() {
            return Err(invalid("minimum-material restart changed its original limit family"));
        }
        checkpoint.accepted.validate(n, cases, options)?;
        if let Some(best) = &checkpoint.best_feasible { best.validate(n, cases, options)?; }
        let lower = vec![options.density_floor; n]; let upper = vec![1.0; n];
        // Check the original numerical state before reservation; a malformed
        // derivative or work counter is never repaired by adding restore cost.
        let original = ProjectedAlState::try_restore::<ResponseError3>(checkpoint.optimizer.clone(), &lower, &upper, options.optimizer)?;
        let work = original.work();
        let minimum = work.iterations.checked_add(work.rejected_trials)
            .and_then(|v| v.checked_add(checkpoint.restoration_evaluations)).and_then(|v| v.checked_add(1))
            .ok_or(invalid("minimum-material checkpoint evaluation count overflow"))?;
        if original.point() != checkpoint.accepted.rho
            || original.sample().objective != checkpoint.accepted.objective
            || original.sample().constraint != checkpoint.accepted.constraint
            || work.evaluations < minimum || work.iterations.checked_add(1) != Some(checkpoint.history.len())
            || checkpoint.history.iter().enumerate().any(|(i, r)| r.iteration != i
                || !r.volume_fraction.is_finite() || !(0.0..=1.0).contains(&r.volume_fraction)
                || objective_value(r.volume_fraction, r.mass_kg, options) != Some(r.objective)
                || !r.constraint.is_finite() || !r.worst_residual.is_finite())
            || checkpoint.history.last() != Some(&checkpoint.accepted.row(work.iterations)) {
            return Err(invalid("minimum-material checkpoint history/sample mismatch"));
        }
        let best_row = checkpoint.history.iter().filter(|r| feasible_row(r, options.optimizer.tolerance))
            .fold(None, |best: Option<&MinimumMaterialIteration3>, r| {
                if best.is_none_or(|b| objective_inventory(r.volume_fraction, r.mass_kg, options.objective)
                    < objective_inventory(b.volume_fraction, b.mass_kg, options.objective)) { Some(r) } else { best }
            });
        match (best_row, &checkpoint.best_feasible) {
            (None, None) => {},
            (Some(r), Some(best)) if best.row(r.iteration) == *r => {},
            _ => return Err(invalid("minimum-material checkpoint lost or changed its feasible incumbent")),
        }
        if let Some(best) = &checkpoint.best_feasible {
            if best.rho == checkpoint.accepted.rho && best != &checkpoint.accepted {
                return Err(invalid("one checkpoint density has inconsistent endpoint observations"));
            }
        }
        let cost = checkpoint.restoration_cost();
        let mut reserved = checkpoint.optimizer.clone();
        reserved.work.evaluations = reserved.work.evaluations.checked_add(cost)
            .ok_or(invalid("minimum-material restoration evaluation overflow"))?;
        let restoration_evaluations = checkpoint.restoration_evaluations.checked_add(cost)
            .ok_or(invalid("minimum-material restoration counter overflow"))?;
        let state = ProjectedAlState::try_restore(reserved.clone(), &lower, &upper, options.optimizer)?;
        // Charge the full bounded reservation before any callback can stop it.
        // No mutation of the physical study has occurred at this boundary.
        checkpoint.optimizer = reserved;
        checkpoint.restoration_evaluations = restoration_evaluations;
        control.checkpoint("minimum-material-restore-reserved").map_err(physical)?;
        let mut accepted = evaluate(study, cases, state.point(), options, control)
            .map_err(ProjectedAlError::Evaluation)?;
        if optimization_sample(&accepted, options).map_err(ProjectedAlError::Evaluation)? != *state.sample()
            || MinimumMaterialEndpoint3::from_evaluation(&accepted, options) != checkpoint.accepted {
            return Err(invalid("restored minimum-material endpoint changed its physical problem"));
        }
        let best_feasible = match &checkpoint.best_feasible {
            None => None,
            Some(best) if best.rho == accepted.rho => Some(accepted.clone()),
            Some(best) => {
                control.checkpoint("minimum-material-restore-incumbent").map_err(physical)?;
                let e = evaluate(study, cases, &best.rho, options, control)
                    .map_err(ProjectedAlError::Evaluation)?;
                optimization_sample(&e, options).map_err(ProjectedAlError::Evaluation)?;
                if !e.feasible(options.optimizer.tolerance) || MinimumMaterialEndpoint3::from_evaluation(&e, options) != *best {
                    return Err(invalid("restored minimum-material incumbent changed its physical limits"));
                }
                Some(e)
            },
        };
        control.checkpoint("minimum-material-restore-publish").map_err(physical)?;
        study.operator.set_scales(&accepted.scales).map_err(physical)?;
        accepted.work = control.work();
        Ok(Self { study, cases, control, options, state, accepted, best_feasible,
            history: checkpoint.history.clone(), restoration_evaluations })
    }
}
