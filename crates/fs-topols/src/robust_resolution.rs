//! Finer-grid admission for independent-load projected structural design.
//!
//! Every load is resolved on every declared grid. Weights select the compliance
//! objective, never which loads must pass stress or resolution checks. Geometry,
//! material, boundary loads, area target and stress allowance are unchanged.
//! These are observed mesh/sample checks, not continuous stress certificates.

use fs_cutfem::CutFemError;
use crate::{RobustAggregate, RobustLoadCase, SampledStressEvaluation, SampledStressLimit};
use crate::projected_stress::resolution::{StressResolutionReport, StressResolutionRung};
use crate::resolution::ResolutionPolicy;
use crate::robust_descent::{MultiLoadProjectedProgress, MultiLoadProjectedStage};

/// Actual independently solved responses, in the original load declaration order.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiLoadResolutionRung {
    /// Background-grid level; the first rung is the optimization grid.
    pub level: u32,
    /// Same-geometry compliance, stress maximum/location and sample count per load.
    pub cases: Vec<SampledStressEvaluation>,
}

/// A complete grid-by-load assessment; partial families are never published.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiLoadResolutionReport {
    /// Consecutive coarse-to-fine grids. Fine geometry is prolonged, not reprojected.
    pub rungs: Vec<MultiLoadResolutionRung>,
}

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

impl MultiLoadResolutionReport {
    /// Reapply the original single-load resolution/stress gates to EVERY case,
    /// including zero-weight cases. A tight aggregate cannot mask an unresolved
    /// or overstressed constituent. Empty, partial or mismatched families refuse.
    pub fn refusal(&self, policy: ResolutionPolicy, target: f64,
        limit: SampledStressLimit, loads: &[RobustLoadCase]) -> Result<Option<String>, CutFemError>
    {
        let first = self.rungs.first().ok_or_else(|| invalid("empty multi-load mesh assessment"))?;
        policy.validate(first.level)?;
        if !(1..=64).contains(&loads.len()) || !loads.iter().any(|case| case.weight() > 0.0)
            || self.rungs.len() != policy.extra_levels as usize + 1
        { return Err(invalid("incomplete multi-load mesh assessment")); }
        for rung in &self.rungs {
            if rung.cases.len() != loads.len() {
                return Err(invalid("mesh assessment omitted an independent load case"));
            }
            if rung.cases.iter().any(|case| case.snapshot != rung.cases[0].snapshot
                || case.volume.to_bits() != rung.cases[0].volume.to_bits())
            { return Err(invalid("mesh load family mixes different geometries")); }
        }
        for case in 0..loads.len() {
            let column = StressResolutionReport { rungs: self.rungs.iter().map(|rung|
                StressResolutionRung { level: rung.level, evaluation: rung.cases[case].clone() }).collect() };
            if let Some(reason) = column.refusal(policy, target, limit)? {
                return Ok(Some(format!("load case {case}: {reason}")));
            }
        }
        Ok(None)
    }

    /// Strict decrease of the declared aggregate on EVERY matching grid. Loads
    /// may trade off within that objective; none may escape the unweighted stress
    /// or individual resolution gates. Worst-weighted active cases may change.
    #[allow(clippy::too_many_arguments)]
    pub fn comparison_refusal(&self, baseline: &Self, policy: ResolutionPolicy,
        target: f64, limit: SampledStressLimit, loads: &[RobustLoadCase],
        aggregate: RobustAggregate, decrease: f64) -> Result<Option<String>, CutFemError>
    {
        if !decrease.is_finite() || !(0.0..1.0).contains(&decrease) {
            return Err(invalid("invalid multi-load mesh compliance decrease"));
        }
        if let Some(reason) = baseline.refusal(policy, target, limit, loads)? { return Ok(Some(reason)); }
        if let Some(reason) = self.refusal(policy, target, limit, loads)? { return Ok(Some(reason)); }
        if baseline.rungs[0].level != self.rungs[0].level {
            return Err(invalid("multi-load mesh comparison changed grid levels"));
        }
        for (old, new) in baseline.rungs.iter().zip(&self.rungs) {
            if (old.cases[0].volume - new.cases[0].volume).abs() > policy.area_tolerance {
                return Ok(Some(format!("level {} changes material area beyond the allowance", new.level)));
            }
            let a = aggregate_compliance(&old.cases, loads, aggregate)?;
            let b = aggregate_compliance(&new.cases, loads, aggregate)?;
            if !(b < a * (1.0 - decrease)) {
                return Ok(Some(format!("level {} does not reproduce the requested independent-load aggregate decrease", new.level)));
            }
        }
        Ok(None)
    }
}

fn aggregate_compliance(states: &[SampledStressEvaluation], loads: &[RobustLoadCase],
    aggregate: RobustAggregate) -> Result<f64, CutFemError>
{
    let mut value = 0.0_f64;
    for (state, load) in states.iter().zip(loads) {
        let term = load.weight() * state.compliance;
        if !term.is_finite() { return Err(invalid("mesh load objective overflowed")); }
        value = match aggregate {
            RobustAggregate::WeightedSum => value + term,
            RobustAggregate::WorstWeightedCase => value.max(term),
        };
    }
    if !value.is_finite() { return Err(invalid("mesh load objective overflowed")); }
    Ok(value)
}

/// Cooperative boundaries of a finer-grid assessment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiLoadResolutionStage {
    /// Before prolongation and preparing this grid's operator.
    Refine(u32),
    /// Before/after an actual independent equilibrium solve.
    CaseSolve { level: u32, case: usize, complete: bool },
    /// Progress inside the original CG/correction solver.
    CaseIterations { level: u32, case: usize, iterations: usize },
    /// Before one cell's stress probes on a solved displacement field.
    StressCell { level: u32, case: usize, cell: usize },
    /// Before returning a complete assessment, never a partial grid/load family.
    Publish,
}

/// The original projected search interleaved with baseline/candidate mesh work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiLoadMeshStage {
    /// Assessment of the unchanged accepted endpoint before geometry search.
    Baseline(MultiLoadResolutionStage),
    /// Assessment of one otherwise admissible candidate.
    Candidate { index: usize, stage: MultiLoadResolutionStage },
    /// Existing proposal, projection, case, stress and final publication boundaries.
    Optimizer(MultiLoadProjectedStage),
}

/// Complete measurements or numerical refusal for one assessed candidate.
#[derive(Debug, Clone)]
pub struct MultiLoadMeshCandidate {
    /// Original bounded-search candidate index.
    pub index: usize,
    /// Absent on numerical refusal; a partly solved family is not a report.
    pub report: Option<MultiLoadResolutionReport>,
    /// `None` only when every mesh/load gate passed.
    pub refusal: Option<String>,
}

/// Mesh-checked progress with unchanged accepted-state and work-budget semantics.
#[derive(Debug, Clone)]
pub enum MultiLoadMeshProgress {
    /// No geometry update: the current design did not pass all finer-grid gates.
    UnresolvedBaseline { baseline: MultiLoadResolutionReport, reason: String },
    /// Existing bounded search, with additional pre-publication mesh admission.
    Searched { baseline: MultiLoadResolutionReport, candidates: Vec<MultiLoadMeshCandidate>,
        progress: MultiLoadProjectedProgress },
    /// No complete baseline assessment fits in the remaining solve allowance.
    SolveBudget,
    /// The existing accepted-update count is complete; no physics was run.
    IterationLimit,
}
