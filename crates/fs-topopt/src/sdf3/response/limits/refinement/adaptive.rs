//! Adaptive minimum-material design with a fresh final-grid feasibility check.
//! Reuse the per-limit estimator and reference refitter; no new physics or AL.
use super::*;
use super::super::design::MinimumMaterialOptions3;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport};
use fs_cutfem::octree3::{Octree3, OctreeError3};

type FitError = ProjectedAlError<ResponseError3>;

/// The factory receives distinct diagnostic and selectively marked trees.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MinimumMaterialBuild3 { Probe, Refit }

/// Per-grid optimization and whole-campaign refinement controls. Original limit
/// laws, smoothing, material, filter and physical region policy remain fixed.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveMinimumMaterialOptions3 {
    pub fit: MinimumMaterialOptions3,
    pub numerical: GoalOptions3,
    pub steps_per_grid: usize,
    /// Total promotions, excluding the initial grid; never reset by run().
    pub max_refinements: usize,
    pub max_transfer_terms: usize,
    pub theta: f64,
    pub max_marks: usize,
    /// Each physical response difference divided by that limit's fixed scale.
    pub response_tolerance: f64,
    /// Absolute difference in numerical material-volume fraction.
    pub volume_tolerance: f64,
    /// Additional mass-mode change tolerance, in units of the chosen
    /// `MassKg::scale_kg`. The geometric volume tolerance remains separate.
    pub mass_tolerance: f64,
}
impl Default for AdaptiveMinimumMaterialOptions3 {
    fn default() -> Self {
        Self { fit: MinimumMaterialOptions3::default(), numerical: GoalOptions3::default(),
            steps_per_grid: 80, max_refinements: 3, max_transfer_terms: 2_000_000,
            theta: 0.5, max_marks: 8, response_tolerance: 1e-3, volume_tolerance: 1e-4, mass_tolerance: 1e-4 }
    }
}
impl AdaptiveMinimumMaterialOptions3 {
    fn validate(self) -> Result<(), FitError> {
        if !self.theta.is_finite() || self.theta <= 0.0 || self.theta > 1.0
            || self.max_marks == 0 || self.max_transfer_terms == 0
            || [self.response_tolerance, self.volume_tolerance, self.mass_tolerance].iter().any(|x| !x.is_finite() || *x < 0.0)
            || [self.numerical.residual_tolerance, self.numerical.identity_tolerance]
                .iter().any(|x| !x.is_finite() || *x <= 0.0 || *x >= 1.0) {
            return Err(ProjectedAlError::Invalid("invalid minimum-material adaptive policy"));
        }
        Ok(())
    }
}

/// Current retained design assessed on one enriched discretization, not all
/// finer grids or the continuum. Feasibility is separate from response change.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MinimumMaterialAssessment3 {
    pub response_change: f64,
    pub volume_change: f64,
    /// Absolute change in the selected, normalized inventory objective.
    pub objective_change: f64,
    /// Actual absolute kg change, absent for volume-only design.
    pub mass_change_kg: Option<f64>,
    /// BOTH the smooth aggregate and every original fine residual pass the
    /// same tolerance used by the minimum-material optimizer's feasible archive.
    pub fine_feasible: bool,
}
impl MinimumMaterialAssessment3 {
    fn meets(self, o: AdaptiveMinimumMaterialOptions3) -> bool {
        self.fine_feasible && self.response_change <= o.response_tolerance
            && self.volume_change <= o.volume_tolerance
            && match o.fit.objective {
                MinimumMaterialObjective3::VolumeFraction => true,
                MinimumMaterialObjective3::MassKg { .. } => self.objective_change <= o.mass_tolerance,
            }
    }
}

/// Outer stops do not replace the per-grid optimizer's actual AL/KKT report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveMinimumMaterialStop3 {
    /// The CURRENT source and probe satisfy all numerical limits and declared
    /// change tolerances. Not continuum feasibility or optimizer stationarity.
    TwoGridLimitsSatisfied,
    RefinementLimit,
    Paused,
    /// Failed numerical checks without local response/inventory marks; never a pass.
    NoRefinementSignal,
    NoFeasibleSource,
    /// The completed candidate has no accepted feasible point. Source retained.
    NoFeasibleCandidate,
}

/// A completed selective promotion. Inventory need not decrease across grids.
#[derive(Debug, Clone)]
pub struct MinimumMaterialAdaptiveRound3 {
    pub cells: [usize; 2],
    /// Source, freshly solved baseline and selected feasible candidate volumes.
    pub volumes: [f64; 3],
    /// Source, fresh baseline and selected feasible candidate mass in kg.
    /// Absent for the volume-only objective.
    pub masses_kg: Option<[f64; 3]>,
    /// The same three endpoints in the chosen normalized objective units.
    pub objectives: [f64; 3],
    pub source_assessment: MinimumMaterialAssessment3,
    pub marking: GoalMarking3,
    /// KKT evidence belongs to the latest AL point, not an earlier incumbent.
    pub optimizer: ProjectedAlReport,
    pub selected_latest: bool,
    pub work: SolveWork,
}
#[derive(Debug)]
pub enum AdaptiveMinimumMaterialError3<E> {
    Build(E), Fit(FitError), Refinement(GoalRefinementError3), Tree(OctreeError3),
    Evaluation(EvaluationStop), Invalid(&'static str),
}
impl<E: std::fmt::Debug> std::fmt::Display for AdaptiveMinimumMaterialError3<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "adaptive minimum-material design stopped: {self:?}")
    }
}
impl<E: std::fmt::Debug> std::error::Error for AdaptiveMinimumMaterialError3<E> {}

fn selected(fit: &ReferenceMinimumMaterialFit3) -> &ResponseLimitEvaluation3 {
    fit.best_feasible.as_ref().unwrap_or(&fit.accepted)
}
fn on_tree<O: AdaptiveSdf3Elasticity>(tree: &Octree3, op: &O) -> bool {
    !op.adaptive().leaves().is_empty()
        && op.adaptive().leaves().iter().all(|leaf| tree.leaves().contains(leaf))
}
fn assessment(e: &ResponseLimitRefinement3, options: AdaptiveMinimumMaterialOptions3)
    -> Result<MinimumMaterialAssessment3, &'static str> {
    if e.objective != options.fit.objective || e.volume_scale.to_bits() != options.fit.volume_scale.to_bits() {
        return Err("refinement inventory policy changed during the campaign");
    }
    let mut response_change = 0.0_f64;
    for row in &e.cases { for limit in &row.limits {
        let difference = (limit.fine_residual-limit.coarse_residual).abs();
        if !difference.is_finite() { return Err("nonfinite individual limit change"); }
        response_change = response_change.max(difference);
    } }
    let volume_change = (e.fine_volume_fraction-e.coarse_volume_fraction).abs();
    if !volume_change.is_finite() { return Err("nonfinite material-volume change"); }
    let (mass_change_kg, objective_change) = match options.fit.objective {
        MinimumMaterialObjective3::VolumeFraction => (None,
            inventory::normalized(volume_change, options.fit.volume_scale)
                .map_err(|_| "unrepresentable material-volume change")?),
        MinimumMaterialObjective3::MassKg { scale_kg } => {
            let change = (e.fine_mass_kg.ok_or("missing enriched material mass")?
                - e.coarse_mass_kg.ok_or("missing source material mass")?).abs();
            (Some(change), inventory::normalized(change, scale_kg)
                .map_err(|_| "unrepresentable physical mass change")?)
        },
    };
    Ok(MinimumMaterialAssessment3 { response_change, volume_change, objective_change, mass_change_kg,
        fine_feasible: e.fine_feasible(options.fit.optimizer.tolerance) })
}
fn endpoint_objective(e: &ResponseLimitEvaluation3, options: MinimumMaterialOptions3)
    -> Result<f64, &'static str> {
    let (value, scale) = match options.objective {
        MinimumMaterialObjective3::VolumeFraction => (e.volume_fraction, options.volume_scale),
        MinimumMaterialObjective3::MassKg { scale_kg } =>
            (e.mass_kg.ok_or("refitted mass objective lost its physical inventory")?, scale_kg),
    };
    inventory::normalized(value, scale).map_err(|_| "unrepresentable refitted material objective")
}

/// Fit -> assess EVERY limit -> mark -> selectively refit, with fixed reference
/// laws. Only raw density controls transfer into optimization; the diagnostic
/// probe instead inherits physical stiffness. Nonlinear force windows have no
/// general monotone density-restoration rule, so candidate feasibility is earned
/// by its own solves, never by scaling volume or trusting the source feasibility.
///
/// All calls must reuse the same SolveControl; the factory must capture one
/// bounded geometry control and unchanged geometry/material/support laws. This
/// owner reconstructs the same filter/SIMP policy. The existing refitter carries
/// protected solid/void regions through the exact geometric parent map.
/// The selected volume or mass objective is retained on every probe and refit.
/// Mass mode checks native kg inventory change and adds its localized cell
/// corrections to marking. Responses and geometric volume remain independent
/// stopping criteria; all checks concern this one enriched discretization.
///
/// Complete probes are cached across pauses. A candidate is installed only after
/// fitting succeeds, a feasible incumbent exists, and the promotion checkpoint
/// passes. Late failures preserve the source and any complete candidate. After
/// promotion the cached source probe is invalidated; the final grid must earn a
/// NEW assessment even at the refinement cap. Pauses retain campaign state, not
/// a partially running candidate optimizer; retrying that candidate spends work.
pub struct AdaptiveMinimumMaterialStudy3<'a, O: AdaptiveSdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    tree: Octree3,
    cases: &'a [ReferenceResponseLimitCase3<'a>],
    options: AdaptiveMinimumMaterialOptions3,
    fit: ReferenceMinimumMaterialFit3,
    rounds: Vec<MinimumMaterialAdaptiveRound3>,
    probe: Option<ResponseLimitRefinement3>,
    unpromoted: Option<(Octree3, MinimumMaterialRefit3<O>)>,
}
impl<'a, O: AdaptiveSdf3Elasticity> AdaptiveMinimumMaterialStudy3<'a, O> {
    pub fn new(study: &'a mut CutDensityStudy3<O>, tree: Octree3,
        cases: &'a [ReferenceResponseLimitCase3<'a>], rho: &[f64], options: AdaptiveMinimumMaterialOptions3,
        control: &mut SolveControl<'_>) -> Result<Self, FitError> {
        study.require_material_body_recipe("minimum-material-adapt-requires-material-acceleration")
            .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
        options.validate()?;
        if !on_tree(&tree, &study.operator) { return Err(ProjectedAlError::Invalid("initial minimum-material background mismatch")); }
        let fit = study.fit_reference_minimum_material(cases, rho, options.fit, options.steps_per_grid, control)?;
        study.operator.set_scales(&selected(&fit).scales)
            .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
        Ok(Self { study, tree, cases, options, fit, rounds: Vec::new(), probe: None, unpromoted: None })
    }
    #[must_use] pub fn study(&self) -> &CutDensityStudy3<O> { self.study }
    #[must_use] pub fn tree(&self) -> &Octree3 { &self.tree }
    /// Best feasible accepted point, or latest point when none exists. Consult
    /// feasible() and fit().outcome; the accessor itself grants no acceptance.
    #[must_use] pub fn accepted(&self) -> &ResponseLimitEvaluation3 { selected(&self.fit) }
    #[must_use] pub fn fit(&self) -> &ReferenceMinimumMaterialFit3 { &self.fit }
    #[must_use] pub fn rounds(&self) -> &[MinimumMaterialAdaptiveRound3] { &self.rounds }
    /// Evidence only for exactly the current selected source design.
    #[must_use] pub fn probe(&self) -> Option<&ResponseLimitRefinement3> { self.probe.as_ref() }
    #[must_use] pub fn unpromoted(&self) -> Option<&(Octree3, MinimumMaterialRefit3<O>)> { self.unpromoted.as_ref() }
    pub fn assessment(&self) -> Result<Option<MinimumMaterialAssessment3>, &'static str> {
        self.probe.as_ref().map(|e| assessment(e, self.options)).transpose()
    }

    /// Up to additional_refinements promotions, bounded by the cumulative cap.
    /// Even run(0) assesses an unprobed source. A satisfied stopping test includes
    /// fine feasibility, not merely a small change or a coarse feasible optimum.
    pub fn run<E>(&mut self, additional_refinements: usize,
        mut build: impl FnMut(&Octree3, MinimumMaterialBuild3) -> Result<O, E>,
        control: &mut SolveControl<'_>) -> Result<AdaptiveMinimumMaterialStop3, AdaptiveMinimumMaterialError3<E>> {
        use AdaptiveMinimumMaterialError3 as Error;
        use AdaptiveMinimumMaterialStop3 as Stop;
        control.checkpoint("minimum-material-adapt-start").map_err(Error::Evaluation)?;
        self.fit.outcome.as_ref().map_err(|e| Error::Fit(e.clone()))?;
        if self.fit.best_feasible.is_none() { return Ok(Stop::NoFeasibleSource); }
        let mut remaining = additional_refinements;
        loop {
            control.checkpoint("minimum-material-adapt-assess").map_err(Error::Evaluation)?;
            if self.probe.is_none() {
                let tree = self.tree.refined(self.study.operator.adaptive().leaves(), || poll(control)).map_err(Error::Tree)?;
                let op = build(&tree, MinimumMaterialBuild3::Probe).map_err(Error::Build)?;
                if !on_tree(&tree, &op) { return Err(Error::Invalid("probe factory returned a different background")); }
                let options = ResponseLimitRefinementOptions3 { limits: self.options.fit.limits,
                    numerical: self.options.numerical, max_transfer_terms: self.options.max_transfer_terms,
                    objective: self.options.fit.objective, volume_scale: self.options.fit.volume_scale };
                let probe = self.study.estimate_response_limit_enrichment(op, selected(&self.fit),
                    self.cases, options, control).map_err(Error::Refinement)?;
                self.probe = Some(probe);
            }
            let assessed = self.assessment().map_err(Error::Invalid)?.expect("complete probe");
            if assessed.meets(self.options) { return Ok(Stop::TwoGridLimitsSatisfied); }
            if self.rounds.len() >= self.options.max_refinements { return Ok(Stop::RefinementLimit); }
            if remaining == 0 { return Ok(Stop::Paused); }
            let marking = self.probe.as_ref().expect("complete probe")
                .mark(self.options.theta, self.options.max_marks, || poll(control))
                .map_err(|e| Error::Refinement(e.into()))?;
            if marking.marked.is_empty() { return Ok(Stop::NoRefinementSignal); }
            let tree = self.tree.refined(&marking.marked, || poll(control)).map_err(Error::Tree)?;
            control.checkpoint("minimum-material-adapt-rebuild").map_err(Error::Evaluation)?;
            let op = build(&tree, MinimumMaterialBuild3::Refit).map_err(Error::Build)?;
            if !on_tree(&tree, &op) { return Err(Error::Invalid("refit factory returned a different background")); }
            let mut candidate = CutDensityStudy3::new(op, self.study.filter_radius, self.study.params);
            self.study.inherit_material_acceleration(&mut candidate, control).map_err(Error::Evaluation)?;
            let mut next = self.study.refit_reference_minimum_material(candidate, selected(&self.fit), self.cases,
                self.options.fit, self.options.steps_per_grid, self.options.max_transfer_terms, control).map_err(Error::Fit)?;
            if let Err(error) = &next.fit.outcome {
                let error = Error::Fit(error.clone());
                self.unpromoted = Some((tree, next)); return Err(error);
            }
            let Some(best) = next.fit.best_feasible.as_ref() else {
                self.unpromoted = Some((tree, next)); return Ok(Stop::NoFeasibleCandidate);
            };
            let masses_kg = match self.options.fit.objective {
                MinimumMaterialObjective3::VolumeFraction => None,
                MinimumMaterialObjective3::MassKg { .. } => Some([
                    self.accepted().mass_kg.ok_or(Error::Invalid("source mass objective lost its inventory"))?,
                    next.fit.initial.mass_kg.ok_or(Error::Invalid("refit baseline lost its inventory"))?,
                    best.mass_kg.ok_or(Error::Invalid("refit incumbent lost its inventory"))?,
                ]),
            };
            let objectives = [endpoint_objective(self.accepted(), self.options.fit).map_err(Error::Invalid)?,
                endpoint_objective(&next.fit.initial, self.options.fit).map_err(Error::Invalid)?,
                endpoint_objective(best, self.options.fit).map_err(Error::Invalid)?];
            let admission = control.checkpoint("minimum-material-adapt-promote").map_err(Error::Evaluation)
                .and_then(|()| next.study.operator.set_scales(&best.scales)
                    .map_err(|e| Error::Fit(ProjectedAlError::Evaluation(e.into()))));
            if let Err(error) = admission { self.unpromoted = Some((tree, next)); return Err(error); }
            let round = MinimumMaterialAdaptiveRound3 {
                cells: [self.study.cells(), next.study.cells()],
                volumes: [self.accepted().volume_fraction, next.fit.initial.volume_fraction, best.volume_fraction],
                masses_kg, objectives,
                source_assessment: assessed, marking,
                optimizer: next.fit.outcome.as_ref().expect("complete fit").clone(),
                selected_latest: best.rho == next.fit.accepted.rho, work: control.work(),
            };
            // No callback inside this promotion: operator, fields, tree and
            // history change together, never to a partially accepted trial.
            self.rounds.push(round); *self.study = next.study; self.tree = tree; self.fit = next.fit;
            self.probe = None; self.unpromoted = None; remaining -= 1;
            control.checkpoint("minimum-material-adapt-accepted").map_err(Error::Evaluation)?;
        }
    }
}
