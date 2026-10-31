//! Repeated mixed-response fitting driven by actual two-grid observations.
//! The enriched probe and the selectively rebuilt design are different grids.
use super::*;
use super::super::refinement::{ResponseRefinement3, ResponseRefinementOptions3};
use crate::sdf3::regions::inherit_physical_regions;
use crate::sdf3_goal::GoalRefinementError3;
use fs_cutfem::octree3::{Octree3, OctreeError3};
use fs_dwr::elasticity3::{GoalMarking3, GoalOptions3};

/// Why the geometry factory is being called. Both calls must use the same
/// reference geometry, material and support laws and a shared geometry budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveResponseBuild3 {
    /// Enrich every current active leaf for an observational error probe.
    Probe,
    /// Rebuild only marked leaves plus the octree's required balance closure.
    Refit,
}

/// Fixed controls for a mixed-response adaptive campaign. The material law,
/// filter radius, physical regions and experiment laws stay fixed across grids.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveResponseOptions3 {
    /// Existing projected optimizer and physical response policy.
    pub fit: ProjectedResponseOptions3,
    /// True-residual and two-grid identity checks, not stopping tolerances.
    pub numerical: GoalOptions3,
    /// Maximum accepted optimizer updates on EACH grid, including the first.
    pub steps_per_grid: usize,
    /// Cumulative number of selectively promoted grids, excluding the first.
    pub max_refinements: usize,
    /// Existing geometric transfer allowance, shared by each transfer operation.
    pub max_transfer_terms: usize,
    /// Requested fraction of the summed absolute local residual mass.
    pub theta: f64,
    /// Mark cap; partial coverage is reported, not claimed to meet `theta`.
    pub max_marks: usize,
    /// Maximum absolute loss change divided by `fit.objective_scale`.
    pub objective_tolerance: f64,
    /// Maximum absolute change of ANY observation divided by its target scale.
    /// Includes zero-weight rows; cancellation between rows cannot pass this.
    pub response_tolerance: f64,
    /// Maximum absolute change of numerical material-volume fraction.
    pub volume_tolerance: f64,
}
impl Default for AdaptiveResponseOptions3 {
    fn default() -> Self {
        Self { fit: ProjectedResponseOptions3::default(), numerical: GoalOptions3::default(),
            steps_per_grid: 40, max_refinements: 3, max_transfer_terms: 2_000_000,
            theta: 0.5, max_marks: 8, objective_tolerance: 1e-3,
            response_tolerance: 1e-3, volume_tolerance: 1e-4 }
    }
}
impl AdaptiveResponseOptions3 {
    fn validate(self) -> Result<(), FitError> {
        if !self.theta.is_finite() || self.theta <= 0.0 || self.theta > 1.0
            || self.max_marks == 0 || self.max_transfer_terms == 0
            || [self.objective_tolerance, self.response_tolerance, self.volume_tolerance]
                .iter().any(|v| !v.is_finite() || *v < 0.0)
            || [self.numerical.residual_tolerance, self.numerical.identity_tolerance]
                .iter().any(|v| !v.is_finite() || *v <= 0.0 || *v >= 1.0) {
            return Err(ProjectedAlError::Invalid("invalid adaptive response controls"));
        }
        Ok(())
    }
}

/// Measured changes at ONE retained design. These are neither continuum bounds
/// nor optimization-stationarity tests. Each observation is compared separately.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveResponseChange3 {
    /// Absolute loss difference / declared fixed optimizer objective scale.
    pub objective: f64,
    /// Largest absolute displacement-observation change / its target scale.
    pub displacement: f64,
    /// Largest absolute force/moment-observation change / its target scale.
    pub reaction: f64,
    /// Absolute numerical projected-volume-fraction change.
    pub volume: f64,
}
impl AdaptiveResponseChange3 {
    fn meets(self, o: AdaptiveResponseOptions3) -> bool {
        self.objective <= o.objective_tolerance && self.displacement <= o.response_tolerance
            && self.reaction <= o.response_tolerance && self.volume <= o.volume_tolerance
    }
}

/// An outer-loop stop, deliberately separate from the per-grid optimizer stop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveResponseStop3 {
    /// All declared numerical two-grid comparisons passed at the CURRENT grid.
    /// Does not imply KKT convergence, continuum accuracy or unique recovery.
    TwoGridToleranceMet,
    /// The campaign's cumulative selective-refinement allowance was reached.
    RefinementLimit,
    /// This call used its requested refits; the same campaign can continue.
    Paused,
    /// Comparisons failed but no local residual marks were available. This can
    /// occur when offset/measure transfer dominates; it is NOT an accuracy pass.
    NoRefinementSignal,
}

/// One completed, atomically promoted refinement. No cross-grid descent claim.
#[derive(Debug, Clone)]
pub struct AdaptiveResponseRound3 {
    /// Active cells before and after this selective rebuild.
    pub cells: [usize; 2],
    /// Error comparisons on the SOURCE grid, not on the new fitted design.
    pub source_change: AdaptiveResponseChange3,
    /// Actual marks, achieved mass fraction and mark-cap coverage status.
    pub marking: GoalMarking3,
    /// Source, freshly rebuilt baseline, and retained fine objective respectively.
    pub objectives: [f64; 3],
    /// Latest fine AL report for this round; its KKT values belong to that
    /// final AL iterate, not necessarily the selected feasible incumbent.
    pub optimizer: ProjectedAlReport,
    /// Whether the selected endpoint is also the latest accepted fine AL point.
    pub selected_latest: bool,
    /// Cumulative physical work at this promotion, including all failed trials.
    pub work: SolveWork,
}

/// Failures preserve the last promoted source. A post-promotion cancellation
/// instead leaves the complete new grid installed, with matching fields/history.
#[derive(Debug)]
pub enum AdaptiveResponseError3<E> {
    /// Caller-owned, budgeted geometry factory failed.
    Build(E),
    /// Initial or candidate projected fitting failed.
    Fit(FitError),
    /// The physical two-grid estimate was refused.
    Refinement(GoalRefinementError3),
    /// Octree level/leaf cap, invalid refinement or cancellation.
    Tree(OctreeError3),
    /// Shared solve control stopped the campaign.
    Evaluation(EvaluationStop),
    /// Invalid factory result or absence of an accepted feasible endpoint.
    Invalid(&'static str),
}
impl<E: std::fmt::Debug> std::fmt::Display for AdaptiveResponseError3<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "adaptive response fitting stopped: {self:?}")
    }
}
impl<E: std::fmt::Debug> std::error::Error for AdaptiveResponseError3<E> {}

fn selected(fit: &ReferenceResponseFit3) -> &ResponseEvaluation3 {
    fit.best_feasible.as_ref().unwrap_or(&fit.accepted)
}
fn tree_matches<O: AdaptiveSdf3Elasticity>(tree: &Octree3, op: &O) -> bool {
    !op.adaptive().leaves().is_empty()
        && op.adaptive().leaves().iter().all(|leaf| tree.leaves().contains(leaf))
}
fn change(a: f64, b: f64, scale: f64) -> Result<f64, &'static str> {
    let d = ((a - b) / scale).abs();
    if !a.is_finite() || !b.is_finite() || !scale.is_finite() || scale <= 0.0 || !d.is_finite() {
        Err("nonfinite normalized two-grid response change")
    } else { Ok(d) }
}
fn changes(report: &ResponseRefinement3, cases: &[ReferenceResponseCase3<'_>],
    reactions: Option<&[&[ReactionTarget3<'_>]]>, objective_scale: f64)
    -> Result<AdaptiveResponseChange3, &'static str> {
    if report.cases.len() != cases.len() || reactions.is_some_and(|r| r.len() != cases.len()) {
        return Err("two-grid response case count mismatch");
    }
    let mut result = AdaptiveResponseChange3 {
        objective: change(report.coarse_objective, report.fine_objective, objective_scale)?,
        volume: change(report.coarse_volume_fraction, report.fine_volume_fraction, 1.0)?,
        displacement: 0.0, reaction: 0.0,
    };
    for (i, (r, case)) in report.cases.iter().zip(cases).enumerate() {
        let rt = reactions.map_or(&[][..], |rows| rows[i]);
        if r.coarse_responses.len() != case.targets.len() || r.fine_responses.len() != case.targets.len()
            || r.coarse_reactions.len() != rt.len() || r.fine_reactions.len() != rt.len() {
            return Err("two-grid response observation count mismatch");
        }
        for ((a, b), t) in r.coarse_responses.iter().zip(&r.fine_responses).zip(case.targets) {
            result.displacement = result.displacement.max(change(*a, *b, t.scale)?);
        }
        for ((a, b), t) in r.coarse_reactions.iter().zip(&r.fine_reactions).zip(rt) {
            result.reaction = result.reaction.max(change(*a, *b, t.scale)?);
        }
    }
    Ok(result)
}

/// Reusable fit/estimate/mark/refit owner for fixed reference experiments.
/// `new` fits the first grid using the existing projected optimizer. After a
/// baseline exists, its late error remains inspectable in `fit().outcome` and
/// `run` refuses to refine past it. Subsequent candidate errors do not replace
/// the source; a completed but unpromoted candidate is retained separately.
///
/// Pure callbacks and the SAME cumulative SolveControl and geometry control
/// must be used throughout the campaign. `run` preserves its last probe across
/// pauses, avoiding duplicate solves. A new grid always starts a fresh AL state;
/// this is not same-grid optimizer checkpoint restoration. Builders return only
/// elasticity; this owner rebuilds the same filter/SIMP policy and inherits
/// protected physical regions through the checked geometric parent map.
pub struct AdaptiveResponseStudy3<'a, O: AdaptiveSdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    tree: Octree3,
    cases: &'a [ReferenceResponseCase3<'a>],
    reactions: Option<&'a [&'a [ReactionTarget3<'a>]]>,
    options: AdaptiveResponseOptions3,
    fit: ReferenceResponseFit3,
    rounds: Vec<AdaptiveResponseRound3>,
    probe: Option<ResponseRefinement3>,
    unpromoted: Option<(Octree3, ResponseRefit3<O>)>,
}
impl<'a, O: AdaptiveSdf3Elasticity> AdaptiveResponseStudy3<'a, O> {
    /// Fit a first grid. `None` selects displacement-only experiments; otherwise
    /// supply one reaction slice per case, including empty slices where needed.
    /// Invalid policy/tree/family inputs fail before reference-law callbacks.
    /// The selected feasible incumbent is installed without relabeling the
    /// latest AL report's stationarity or resetting any spent solve budget.
    #[allow(clippy::too_many_arguments)]
    pub fn new(study: &'a mut CutDensityStudy3<O>, tree: Octree3,
        cases: &'a [ReferenceResponseCase3<'a>], reactions: Option<&'a [&'a [ReactionTarget3<'a>]]>,
        rho: &[f64], options: AdaptiveResponseOptions3, control: &mut SolveControl<'_>)
        -> Result<Self, FitError> {
        study.require_material_body_recipe("response-adapt-requires-material-acceleration")
            .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
        options.validate()?;
        if !tree_matches(&tree, &study.operator) {
            return Err(ProjectedAlError::Invalid("source operator does not match adaptive background"));
        }
        let fit = study.fit_reference_response_family(cases, reactions, rho, options.fit,
            options.steps_per_grid, control)?;
        study.operator.set_scales(&selected(&fit).scales).map_err(physics)?;
        Ok(Self { study, tree, cases, reactions, options, fit, rounds: Vec::new(),
            probe: None, unpromoted: None })
    }
    /// Retained source operator, aligned with `accepted()`, never a failed trial.
    #[must_use] pub fn study(&self) -> &CutDensityStudy3<O> { self.study }
    /// Background of the retained source, including inactive cells.
    #[must_use] pub fn tree(&self) -> &Octree3 { &self.tree }
    /// Least-objective feasible accepted design, or latest point when none exists.
    #[must_use] pub fn accepted(&self) -> &ResponseEvaluation3 { selected(&self.fit) }
    /// Full current-grid fitting prefix, including its actual AL stop/error.
    #[must_use] pub fn fit(&self) -> &ReferenceResponseFit3 { &self.fit }
    /// Completed selective promotions only; the initial grid is excluded.
    #[must_use] pub fn rounds(&self) -> &[AdaptiveResponseRound3] { &self.rounds }
    /// Latest complete error probe of exactly the current retained design.
    /// Cleared at every promotion, so stale source evidence cannot pass a stop.
    #[must_use] pub fn probe(&self) -> Option<&ResponseRefinement3> { self.probe.as_ref() }
    /// Candidate retained after a late fitting/promotion failure, never installed.
    #[must_use] pub fn unpromoted(&self) -> Option<&(Octree3, ResponseRefit3<O>)> { self.unpromoted.as_ref() }
    /// Componentwise changes for the current probe, if one has completed.
    pub fn measured_change(&self) -> Result<Option<AdaptiveResponseChange3>, &'static str> {
        self.probe.as_ref().map(|r| changes(r, self.cases, self.reactions,
            self.options.fit.objective_scale)).transpose()
    }

    /// Perform up to `additional_refinements` selective promotions, within the
    /// cumulative campaign cap. Even a zero-refit call may solve a current-grid
    /// probe; the final retained design is always probed before a tolerance pass.
    /// Tree and geometry failures, partial marking, iteration limits and absent
    /// indicators are distinct from numerical tolerance and KKT convergence.
    ///
    /// The factory must build the requested tree under unchanged geometry,
    /// supports/material and its original bounded QuadratureControl. No filter,
    /// load vector, primal/adjoint or optimizer state is transferred nodally.
    pub fn run<E>(&mut self, additional_refinements: usize,
        mut build: impl FnMut(&Octree3, AdaptiveResponseBuild3) -> Result<O, E>,
        control: &mut SolveControl<'_>) -> Result<AdaptiveResponseStop3, AdaptiveResponseError3<E>> {
        use AdaptiveResponseError3 as Error;
        control.checkpoint("response-adapt-start").map_err(Error::Evaluation)?;
        self.fit.outcome.as_ref().map_err(|e| Error::Fit(e.clone()))?;
        if self.fit.best_feasible.is_none() {
            return Err(Error::Invalid("initial fit has no accepted feasible design"));
        }
        let mut remaining = additional_refinements;
        loop {
            control.checkpoint("response-adapt-probe").map_err(Error::Evaluation)?;
            if self.probe.is_none() {
                let leaves = self.study.operator.adaptive().leaves().to_vec();
                let tree = self.tree.refined(&leaves, || poll(control)).map_err(Error::Tree)?;
                let op = build(&tree, AdaptiveResponseBuild3::Probe).map_err(Error::Build)?;
                if !tree_matches(&tree, &op) { return Err(Error::Invalid("probe factory returned another background")); }
                let options = ResponseRefinementOptions3 { response: self.options.fit.response,
                    numerical: self.options.numerical, max_transfer_terms: self.options.max_transfer_terms };
                let accepted = selected(&self.fit);
                let report = match self.reactions {
                    Some(rows) => self.study.estimate_response_enrichment_with_reactions(op, accepted,
                        self.cases, rows, options, control),
                    None => self.study.estimate_response_enrichment(op, accepted, self.cases, options, control),
                }.map_err(Error::Refinement)?;
                self.probe = Some(report);
            }
            let change = self.measured_change().map_err(Error::Invalid)?.expect("probe completed");
            if change.meets(self.options) { return Ok(AdaptiveResponseStop3::TwoGridToleranceMet); }
            if self.rounds.len() >= self.options.max_refinements { return Ok(AdaptiveResponseStop3::RefinementLimit); }
            if remaining == 0 { return Ok(AdaptiveResponseStop3::Paused); }
            let marking = self.probe.as_ref().expect("probe completed")
                .mark(self.options.theta, self.options.max_marks, || poll(control))
                .map_err(|e| Error::Refinement(GoalRefinementError3::from(e)))?;
            if marking.marked.is_empty() { return Ok(AdaptiveResponseStop3::NoRefinementSignal); }
            let tree = self.tree.refined(&marking.marked, || poll(control)).map_err(Error::Tree)?;
            control.checkpoint("response-adapt-rebuild").map_err(Error::Evaluation)?;
            let op = build(&tree, AdaptiveResponseBuild3::Refit).map_err(Error::Build)?;
            if !tree_matches(&tree, &op) { return Err(Error::Invalid("refit factory returned another background")); }
            let mut candidate = CutDensityStudy3::new(op, self.study.filter_radius, self.study.params);
            self.study.inherit_material_acceleration(&mut candidate, control).map_err(Error::Evaluation)?;
            if let Some(regions) = self.study.physical_regions() {
                let labels = {
                    let transfer = AdaptiveTransfer3::new(self.study.operator.adaptive(), candidate.operator.adaptive(),
                        self.options.max_transfer_terms, || poll(control)).map_err(|e| Error::Fit(physics(e)))?;
                    inherit_physical_regions(regions, transfer.parents(), control).map_err(Error::Evaluation)?
                };
                candidate = candidate.with_physical_regions(labels).map_err(Error::Evaluation)?;
            }
            let mut next = self.study.refit_reference_response_family(candidate, selected(&self.fit),
                self.cases, self.reactions, self.options.fit, self.options.steps_per_grid,
                self.options.max_transfer_terms, control).map_err(Error::Fit)?;
            let admission = (|| {
                next.fit.outcome.as_ref().map_err(|e| Error::Fit(e.clone()))?;
                let best = next.fit.best_feasible.as_ref()
                    .ok_or(Error::Invalid("refitted grid has no accepted feasible design"))?;
                control.checkpoint("response-adapt-promote").map_err(Error::Evaluation)?;
                next.study.operator.set_scales(&best.scales).map_err(|e| Error::Fit(physics(e)))
            })();
            if let Err(error) = admission {
                self.unpromoted = Some((tree, next));
                return Err(error);
            }
            let row = AdaptiveResponseRound3 {
                cells: [self.study.cells(), next.study.cells()], source_change: change, marking,
                objectives: [self.accepted().objective, next.fit.initial.objective, selected(&next.fit).objective],
                optimizer: next.fit.outcome.as_ref().expect("admitted outcome").clone(),
                selected_latest: selected(&next.fit).rho == next.fit.accepted.rho, work: control.work(),
            };
            // No callback between these updates. The selected physical fields,
            // background and history are one logical commit, not the last trial.
            self.rounds.push(row);
            *self.study = next.study;
            self.tree = tree;
            self.fit = next.fit;
            self.probe = None;
            self.unpromoted = None;
            remaining -= 1;
            control.checkpoint("response-adapt-accepted").map_err(Error::Evaluation)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_normalized_changes_refuse_overflow_and_nonfinite_values() {
        assert_eq!(change(3.0, 1.0, 4.0).unwrap(), 0.5);
        for (a, b, s) in [(f64::NAN, 0.0, 1.0), (1.0, f64::INFINITY, 1.0),
            (f64::MAX, -f64::MAX, 1.0), (1.0, 0.0, 0.0), (1.0, 0.0, f64::NAN)] {
            assert!(change(a, b, s).is_err());
        }
    }
    #[test]
    fn g0_small_total_loss_does_not_override_observation_or_volume_changes() {
        let o = AdaptiveResponseOptions3::default();
        let zero = AdaptiveResponseChange3 { objective: 0.0, displacement: 0.0, reaction: 0.0, volume: 0.0 };
        assert!(zero.meets(o));
        assert!(!AdaptiveResponseChange3 { displacement: 1.0, ..zero }.meets(o));
        assert!(!AdaptiveResponseChange3 { reaction: 1.0, ..zero }.meets(o));
        assert!(!AdaptiveResponseChange3 { volume: 1.0, ..zero }.meets(o));
    }
}
