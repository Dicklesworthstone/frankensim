//! Mesh-checked sampled-stress constraints inside the existing projected search.
//!
//! The same level set, load, material and stress limit are used at every level.
//! Fine fields are prolonged, not re-projected. These are observed grid/sample
//! checks, not continuous stress or continuum discretization certificates.
use super::*;
use crate::refinement::prolongate_level_set;
use crate::resolution::{MeshCheckStage, ResolutionPolicy, ResolutionStage};
use crate::{Cantilever, GridSdf, OptimizeSettings};

/// Actual independent mechanics and stress probes on one declared grid.
#[derive(Debug, Clone, PartialEq)]
pub struct StressResolutionRung {
    pub level: u32,
    pub evaluation: SampledStressEvaluation,
}

/// All declared grids, including the optimization grid. No interpolated solution.
#[derive(Debug, Clone, PartialEq)]
pub struct StressResolutionReport {
    pub rungs: Vec<StressResolutionRung>,
}

impl StressResolutionReport {
    /// Recheck complete measurements against the unchanged physical constraints.
    /// Used by both live candidate admission and retained-report readers.
    pub fn refusal(&self, policy: ResolutionPolicy, target: f64, limit: SampledStressLimit)
        -> Result<Option<String>, CutFemError>
    {
        let first = self.rungs.first().ok_or_else(|| refused("empty stress-resolution report"))?;
        policy.validate(first.level)?;
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        if !target.is_finite() || !(0.0..=1.0).contains(&target) || target == 0.0
            || self.rungs.len() != policy.extra_levels as usize + 1
        { return Err(refused("invalid stress-resolution target or missing grid levels")); }
        for (i, rung) in self.rungs.iter().enumerate() {
            let s = &rung.evaluation;
            if rung.level != first.level + i as u32
                || !s.compliance.is_finite() || s.compliance < 0.0
                || !s.volume.is_finite() || s.volume <= 0.0
                || !s.sampled_max_von_mises.is_finite() || s.sampled_max_von_mises < 0.0
                || s.sample_count == 0 || s.max_location.iter().any(|x| !x.is_finite() || !(0.0..=1.0).contains(x))
            { return Err(refused("invalid stress-resolution measurements")); }
        }
        // Check stress first so a finer-grid violation is never reported merely
        // as compliance underresolution, and never silently repairs the limit.
        for rung in &self.rungs {
            if rung.evaluation.sampled_max_von_mises > limit.admitted_max() {
                return Ok(Some(format!("level {} sampled von Mises stress {:.17e} exceeds unchanged admitted limit {:.17e}",
                    rung.level, rung.evaluation.sampled_max_von_mises, limit.admitted_max())));
            }
        }
        for pair in self.rungs.windows(2) {
            let (a, b) = (&pair[0].evaluation, &pair[1].evaluation);
            let allowance = policy.absolute_compliance_tolerance + policy.relative_compliance_tolerance * b.compliance;
            let dc = (a.compliance - b.compliance).abs();
            let da = (a.volume - b.volume).abs();
            if !allowance.is_finite() || !dc.is_finite() || !da.is_finite() {
                return Err(refused("stress-resolution comparison overflowed"));
            }
            if dc > allowance || da > policy.area_tolerance || (b.volume - target).abs() > policy.area_tolerance {
                return Ok(Some(format!("level {} compliance or material area is unresolved under the declared mesh policy", pair[1].level)));
            }
        }
        Ok(None)
    }

    /// A candidate must improve on every MATCHING grid, not just the coarse one.
    pub fn comparison_refusal(&self, baseline: &Self, policy: ResolutionPolicy,
        target: f64, limit: SampledStressLimit, decrease: f64) -> Result<Option<String>, CutFemError>
    {
        if !decrease.is_finite() || !(0.0..1.0).contains(&decrease) {
            return Err(refused("invalid mesh-checked compliance decrease"));
        }
        if let Some(reason) = baseline.refusal(policy, target, limit)? { return Ok(Some(reason)); }
        if let Some(reason) = self.refusal(policy, target, limit)? { return Ok(Some(reason)); }
        if baseline.rungs[0].level != self.rungs[0].level {
            return Err(refused("stress-resolution comparison uses different grid levels"));
        }
        for (old, new) in baseline.rungs.iter().zip(&self.rungs) {
            if (old.evaluation.volume - new.evaluation.volume).abs() > policy.area_tolerance {
                return Ok(Some(format!("level {} candidate changes material area beyond the allowance", new.level)));
            }
            if !(new.evaluation.compliance < old.evaluation.compliance * (1.0 - decrease)) {
                return Ok(Some(format!("level {} does not reproduce the requested compliance decrease", new.level)));
            }
        }
        Ok(None)
    }
}

/// Measure unchanged geometry on every declared grid using the original stress
/// evaluator. The one solve per rung supplies compliance, area AND stress.
/// Cancellation in CG or cell sampling returns no partial report.
pub fn assess_stress_resolution_controlled<B>(geometry: &GridSdf, fixture: Cantilever,
    settings: OptimizeSettings, policy: ResolutionPolicy, poll_iters: usize,
    mut control: impl FnMut(ResolutionStage) -> ControlFlow<B>)
    -> Result<ControlFlow<B, StressResolutionReport>, CutFemError>
{
    policy.validate(settings.level)?;
    if poll_iters == 0 || geometry.n() != 1usize << settings.level {
        return Err(refused("stress resolution requires positive polling and a level-matched field"));
    }
    if let ControlFlow::Break(reason) = control(ResolutionStage::Prepare) { return Ok(ControlFlow::Break(reason)); }
    let mut field = geometry.clone();
    let mut rungs = Vec::with_capacity(policy.extra_levels as usize + 1);
    for level in settings.level..=settings.level + policy.extra_levels {
        if level != settings.level {
            if let ControlFlow::Break(reason) = control(ResolutionStage::Refine(level)) { return Ok(ControlFlow::Break(reason)); }
            field = prolongate_level_set(&field, &[])?.geometry;
        }
        let evaluation = match evaluate_sampled_stress_controlled(&field, fixture,
            OptimizeSettings { level, ..settings }, poll_iters,
            |stage| control(ResolutionStage::Evaluate { level, stage }))?
        {
            ControlFlow::Continue(value) => value,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        rungs.push(StressResolutionRung { level, evaluation });
    }
    if let ControlFlow::Break(reason) = control(ResolutionStage::Publish) { return Ok(ControlFlow::Break(reason)); }
    Ok(ControlFlow::Continue(StressResolutionReport { rungs }))
}

/// One assessed candidate; refused fields never replace the accepted optimizer.
#[derive(Debug, Clone)]
pub struct StressCandidateResolution {
    pub index: usize,
    pub report: Option<StressResolutionReport>,
    pub refusal: Option<String>,
}

#[derive(Debug, Clone)]
pub enum MeshCheckedStressProgress {
    /// No geometry search occurred. The original endpoint remains available.
    UnresolvedBaseline { baseline: StressResolutionReport, reason: String },
    /// All candidate work uses the existing bounded projection/search owner.
    Searched { baseline: StressResolutionReport, candidates: Vec<StressCandidateResolution>,
        update: ProjectedStressUpdate },
    /// No new physics or callbacks after the declared update count is reached.
    IterationLimit,
}

impl ProjectedStressOptimizer {
    /// Require same-grid improvement, observed resolution and the ORIGINAL stress
    /// limit at every grid before publishing a projected update. Stress remains
    /// sample-scoped. Work is at most `(extra_levels+1)*(max_candidates+1)`
    /// assessment solves plus the existing proposal/projected-field solves.
    /// Assembly, prolongation and individual quadrature operations are indivisible.
    pub fn advance_one_resolution_controlled<B>(&mut self, policy: ResolutionPolicy,
        mut control: impl FnMut(MeshCheckStage) -> ControlFlow<B>)
        -> Result<ControlFlow<B, MeshCheckedStressProgress>, CutFemError>
    {
        let settings = self.checkpoint().settings();
        policy.validate(settings.level)?;
        if self.checkpoint().is_complete() { return Ok(ControlFlow::Continue(MeshCheckedStressProgress::IterationLimit)); }
        let fixture = self.checkpoint().fixture();
        let search = self.optimizer.controls();
        let limit = self.limit;
        let baseline = match assess_stress_resolution_controlled(self.checkpoint().geometry(), fixture,
            settings, policy, search.poll_iters, |stage| control(MeshCheckStage::Baseline(stage)))?
        {
            ControlFlow::Continue(report) => report,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        if baseline.rungs[0].evaluation != self.current {
            return Err(refused("mesh-checked baseline differs from retained stress/mechanics"));
        }
        if let Some(reason) = baseline.refusal(policy, settings.volfrac, limit)? {
            return Ok(ControlFlow::Continue(MeshCheckedStressProgress::UnresolvedBaseline { baseline, reason }));
        }
        let controller = std::cell::RefCell::new(&mut control);
        let mut candidates = Vec::new();
        let mut stress_checks = Vec::new();
        let mut next = self.current.clone();
        let progress = self.optimizer.advance_one_admitted_controlled(|index, geometry, state, _| {
            let checked = assess_stress_resolution_controlled(geometry, fixture, settings,
                policy, search.poll_iters,
                |stage| (controller.borrow_mut())(MeshCheckStage::Candidate { index, stage }));
            let (report, evaluation, refusal) = match checked {
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Ok(ControlFlow::Continue(report)) => {
                    let first = &report.rungs[0].evaluation;
                    let reason = if first.snapshot != state.snapshot
                        || first.compliance.to_bits() != state.compliance.to_bits()
                        || first.volume.to_bits() != state.volume.to_bits()
                    { Some("stress-resolution solve differs from candidate mechanics".into()) }
                    else { report.comparison_refusal(&baseline, policy, settings.volfrac, limit,
                        search.min_relative_improvement)? };
                    let evaluation = first.clone();
                    if reason.is_none() { next = evaluation.clone(); }
                    (Some(report), Some(evaluation), reason)
                }
                Err(error) => (None, None, Some(format!("stress-resolution solve: {error}"))),
            };
            stress_checks.push(ProjectedStressCheck { index, evaluation, refusal: refusal.clone() });
            candidates.push(StressCandidateResolution { index, report, refusal: refusal.clone() });
            Ok(ControlFlow::Continue(refusal))
        }, |stage| (controller.borrow_mut())(MeshCheckStage::Optimizer(stage)))?;
        match progress {
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(progress) => {
                if matches!(&progress, ProjectedProgress::Accepted(_)) { self.current = next; }
                Ok(ControlFlow::Continue(MeshCheckedStressProgress::Searched { baseline, candidates,
                    update: ProjectedStressUpdate { progress, stress_checks } }))
            }
        }
    }
}

#[cfg(test)]
mod tests;
