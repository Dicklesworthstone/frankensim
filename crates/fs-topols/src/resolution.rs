//! Observed mesh-resolution gates for the existing projected optimizer.
//!
//! Each geometry is bilinearly prolonged, never re-projected or redistanced,
//! then solved by the original CutFEM evaluator at the same declared levels.
//! This is a measured grid-sensitivity screen, NOT a continuum error bound.
//! Agreement can miss a common bias; no Richardson order, stress certificate,
//! optimum, or physical validation follows from passing these gates.
use std::ops::ControlFlow;

use fs_cutfem::CutFemError;
use crate::evaluated::{DesignEvaluationStage, evaluate_compliance_design_controlled};
use crate::projected::{ProjectedOptimizer, ProjectedProgress, ProjectedStage};
use crate::refinement::prolongate_level_set;
use crate::{Cantilever, GridSdf, OptimizeSettings};

/// Extra work is bounded by `(extra_levels + 1) * (max_candidates + 1)` solves
/// per update, in addition to the original proposal/projection solves. The
/// geometry level is capped at seven. CG uses the existing poll-iteration cap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResolutionPolicy {
    pub extra_levels: u32,
    pub absolute_compliance_tolerance: f64,
    pub relative_compliance_tolerance: f64,
    /// Both adjacent-grid area changes AND each finer area relative to the
    /// optimizer's target are checked. No material is added to pass this gate.
    pub area_tolerance: f64,
}
impl ResolutionPolicy {
    /// Refuse malformed policy and excessive refinement before any solve.
    pub fn validate(self, level: u32) -> Result<(), CutFemError> {
        if !(1..=2).contains(&self.extra_levels) || !(1..=6).contains(&level)
            || level.checked_add(self.extra_levels).is_none_or(|fine| fine > 7)
            || !self.absolute_compliance_tolerance.is_finite()
            || self.absolute_compliance_tolerance < 0.0
            || !self.relative_compliance_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_compliance_tolerance)
            || (self.absolute_compliance_tolerance == 0.0 && self.relative_compliance_tolerance == 0.0)
            || !self.area_tolerance.is_finite() || self.area_tolerance <= 0.0
        {
            return Err(invalid("mesh resolution needs 1..=2 extra levels through level 7, nonnegative finite compliance tolerances (at least one positive, relative <1), and positive finite area tolerance"));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResolutionRung {
    pub level: u32,
    pub compliance: f64,
    pub area: f64,
    /// Identifies the actual field evaluated at this level, not a binary mask.
    pub snapshot: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ResolutionReport {
    pub rungs: Vec<ResolutionRung>,
    /// Maximum of the measured differences between adjacent declared grids.
    /// This is not an upper bound on continuum discretization error.
    pub max_compliance_change: f64,
    pub max_area_change: f64,
    pub agrees: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResolutionStage {
    Prepare,
    Refine(u32),
    Evaluate { level: u32, stage: DesignEvaluationStage },
    Publish,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshCheckStage {
    Baseline(ResolutionStage),
    Candidate { index: usize, stage: ResolutionStage },
    Optimizer(ProjectedStage),
}

#[derive(Debug, Clone)]
pub struct CandidateResolution {
    pub index: usize,
    pub report: Option<ResolutionReport>,
    pub refusal: Option<String>,
}

/// An unresolved baseline is different from a completed candidate search.
#[derive(Debug, Clone)]
pub enum MeshCheckedProgress {
    UnresolvedBaseline { baseline: ResolutionReport, reason: String },
    Searched {
        baseline: ResolutionReport,
        candidates: Vec<CandidateResolution>,
        progress: ProjectedProgress,
    },
}

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn summarize(rungs: Vec<ResolutionRung>, policy: ResolutionPolicy) -> Result<ResolutionReport, CutFemError> {
    let mut report = ResolutionReport {
        rungs, max_compliance_change: 0.0, max_area_change: 0.0, agrees: true,
    };
    for pair in report.rungs.windows(2) {
        let dc = (pair[1].compliance - pair[0].compliance).abs();
        let da = (pair[1].area - pair[0].area).abs();
        let allowance = policy.absolute_compliance_tolerance
            + policy.relative_compliance_tolerance * pair[1].compliance.abs();
        if !dc.is_finite() || !da.is_finite() || !allowance.is_finite() {
            return Err(invalid("non-finite mesh-resolution difference or allowance"));
        }
        report.max_compliance_change = report.max_compliance_change.max(dc);
        report.max_area_change = report.max_area_change.max(da);
        report.agrees &= dc <= allowance && da <= policy.area_tolerance;
    }
    Ok(report)
}

/// Measure one unchanged design at the optimization level and one or two finer
/// levels. Every rung is an actual independent solve; interpolated displacement
/// is never reported as a solution. Coarse nodes retain their original bits;
/// interpolated nodes introduce only the prolongator's disclosed rounding.
/// A callback break returns no partial report and leaves the source untouched.
/// Assembly, prolongation and one quadrature call remain indivisible.
pub fn assess_compliance_resolution_controlled<B>(
    geometry: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
    policy: ResolutionPolicy,
    poll_iters: usize,
    mut control: impl FnMut(ResolutionStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, ResolutionReport>, CutFemError> {
    policy.validate(settings.level)?;
    if poll_iters == 0 || geometry.n() != 1usize << settings.level {
        return Err(invalid("mesh resolution requires positive CG polling and a level-matched field"));
    }
    if let ControlFlow::Break(reason) = control(ResolutionStage::Prepare) {
        return Ok(ControlFlow::Break(reason));
    }
    let mut field = geometry.clone();
    let mut rungs = Vec::with_capacity(policy.extra_levels as usize + 1);
    for level in settings.level..=settings.level + policy.extra_levels {
        if level != settings.level {
            if let ControlFlow::Break(reason) = control(ResolutionStage::Refine(level)) {
                return Ok(ControlFlow::Break(reason));
            }
            field = prolongate_level_set(&field, &[])?.geometry;
        }
        let state = match evaluate_compliance_design_controlled(
            &field, fixture, OptimizeSettings { level, ..settings }, poll_iters,
            |stage| control(ResolutionStage::Evaluate { level, stage }),
        )? {
            ControlFlow::Continue(state) => state,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        rungs.push(ResolutionRung { level, compliance: state.compliance,
            area: state.volume, snapshot: state.snapshot });
    }
    let report = summarize(rungs, policy)?;
    if let ControlFlow::Break(reason) = control(ResolutionStage::Publish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue(report))
}

fn area_refusal(report: &ResolutionReport, target: f64, tolerance: f64) -> Option<String> {
    report.rungs.iter().skip(1).find(|rung| (rung.area - target).abs() > tolerance)
        .map(|rung| format!("level {} area {:.17e} differs from target {:.17e} beyond {:.17e}",
            rung.level, rung.area, target, tolerance))
}

fn candidate_refusal(baseline: &ResolutionReport, candidate: &ResolutionReport,
    policy: ResolutionPolicy, target: f64, decrease: f64) -> Option<String>
{
    if !candidate.agrees {
        return Some("candidate mesh-resolution changes exceed declared tolerances".into());
    }
    if let Some(reason) = area_refusal(candidate, target, policy.area_tolerance) { return Some(reason); }
    if baseline.rungs.len() != candidate.rungs.len() {
        return Some("candidate and baseline do not cover identical grid levels".into());
    }
    for (old, new) in baseline.rungs.iter().zip(&candidate.rungs) {
        if old.level != new.level { return Some("candidate and baseline grid levels differ".into()); }
        if (old.area - new.area).abs() > policy.area_tolerance {
            return Some(format!("level {} comparison changes material area beyond the allowance", old.level));
        }
        if !(new.compliance < old.compliance * (1.0 - decrease)) {
            return Some(format!("level {} does not reproduce the requested compliance decrease", old.level));
        }
    }
    None
}

impl ProjectedOptimizer {
    /// Apply mesh-resolution and same-grid improvement gates INSIDE the existing
    /// candidate admission hook, before geometry, ordinal or multiplier changes.
    /// The baseline is evaluated once per call; only otherwise admissible coarse
    /// candidates require the additional solves. Failed checks contract the next
    /// proposal using the original bounded search. No refinement-induced area
    /// projection, alternative optimizer, or continuum authority is introduced.
    pub fn advance_one_resolution_controlled<B>(
        &mut self, policy: ResolutionPolicy,
        mut control: impl FnMut(MeshCheckStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, MeshCheckedProgress>, CutFemError> {
        let settings = self.checkpoint().settings();
        let fixture = self.checkpoint().fixture();
        let search = self.controls();
        let baseline = match assess_compliance_resolution_controlled(
            self.checkpoint().geometry(), fixture, settings, policy, search.poll_iters,
            |stage| control(MeshCheckStage::Baseline(stage)),
        )? {
            ControlFlow::Continue(report) => report,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        let first = baseline.rungs[0];
        let current = self.current();
        if first.snapshot != current.snapshot || first.compliance.to_bits() != current.compliance.to_bits()
            || first.area.to_bits() != current.volume.to_bits()
        {
            return Err(invalid("independent baseline resolution solve differs from accepted mechanics"));
        }
        let refusal = if !baseline.agrees {
            Some("accepted baseline mesh-resolution changes exceed declared tolerances".into())
        } else { area_refusal(&baseline, settings.volfrac, policy.area_tolerance) };
        if let Some(reason) = refusal {
            return Ok(ControlFlow::Continue(MeshCheckedProgress::UnresolvedBaseline { baseline, reason }));
        }
        // Share a single callback between the original optimizer and its
        // constraint hook. The borrows are sequential and never escape a call.
        let controller = std::cell::RefCell::new(&mut control);
        let mut candidates = Vec::new();
        let progress = self.advance_one_admitted_controlled(|index, geometry, state, _| {
            let checked = assess_compliance_resolution_controlled(
                geometry, fixture, settings, policy, search.poll_iters,
                |stage| (controller.borrow_mut())(MeshCheckStage::Candidate { index, stage }),
            );
            let (report, refusal) = match checked {
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Ok(ControlFlow::Continue(report)) => {
                    let first = report.rungs[0];
                    let reason = if first.snapshot != state.snapshot
                        || first.compliance.to_bits() != state.compliance.to_bits()
                        || first.area.to_bits() != state.volume.to_bits() {
                        Some("candidate resolution solve differs from the projected mechanics".into())
                    } else { candidate_refusal(&baseline, &report, policy, settings.volfrac,
                        search.min_relative_improvement) };
                    (Some(report), reason)
                }
                Err(error) => (None, Some(format!("candidate resolution solve: {error}"))),
            };
            candidates.push(CandidateResolution { index, report, refusal: refusal.clone() });
            Ok(ControlFlow::Continue(refusal))
        }, |stage| (controller.borrow_mut())(MeshCheckStage::Optimizer(stage)))?;
        match progress {
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(progress) => Ok(ControlFlow::Continue(
                MeshCheckedProgress::Searched { baseline, candidates, progress })),
        }
    }
}

#[cfg(test)]
mod tests;
