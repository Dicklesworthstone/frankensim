//! Compare design gains on ONE locally refined analysis grid per round.
//!
//! Both unchanged geometries are solved on the same background cells and the
//! same one-level enrichment. Marking includes both designs and every weighted
//! load. Independently adapted grids are never compared as matching evidence.
//! These are Estimated numerical screens, not continuum or stress certificates.
use super::*;

/// Refinement and numerical-gain admission, separate from the optimizer's search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveGainPolicy {
    /// Bounds and endpoint-resolution tolerances for the existing local refiner.
    pub refinement: AdaptiveDwrPolicy,
    /// Allowed change in the gain, as a fraction of its positive coarse value.
    pub max_error_fraction: f64,
    /// Each assessed area must match the original target within this allowance;
    /// the two designs' assessed areas must also agree within it, in square metres.
    pub area_tolerance: f64,
}
impl AdaptiveGainPolicy {
    /// Refuse invalid controls before any physics or refinement allocation.
    pub fn validate(self, level: u32, cases: usize) -> Result<(), CutFemError> {
        self.refinement.validate(level, cases)?;
        if !self.max_error_fraction.is_finite() || self.max_error_fraction <= 0.0
            || self.max_error_fraction >= 1.0 || !self.area_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.area_tolerance) {
            return Err(invalid("adaptive gain needs a fraction in (0,1) and finite area allowance in [0,1)"));
        }
        Ok(())
    }

    /// Complete worst-case comparison: two designs, two solves per case per round.
    pub fn max_solves(self, cases: usize) -> Option<usize> {
        self.refinement.max_solves(cases)?.checked_mul(2)
    }
}

/// Which geometry is being assessed. Neither input is ever modified.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComparedDesign { Baseline, Candidate }

/// Cooperative boundaries; one case's two solves/integration remain indivisible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveGainStage {
    BeforeCase { round: usize, design: ComparedDesign, case: usize },
    RoundComplete(usize),
    BeforeRefine(usize),
    BeforePublish,
}

/// A refusal is bounded search evidence, never an optimum or convergence proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveGainStop {
    ResolvedGain, NonImproving, RoundBudget, CellBudget, LevelLimit, NoIndicators,
}
impl AdaptiveGainStop {
    /// Stable numerical stop label.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ResolvedGain => "resolved-gain", Self::NonImproving => "non-improving",
            Self::RoundBudget => "round-budget", Self::CellBudget => "cell-budget",
            Self::LevelLimit => "level-limit", Self::NoIndicators => "no-indicators",
        }
    }
}

/// Numerical arithmetic for one complete, matching-grid pair of load families.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveGainRow {
    pub baseline: AdaptiveDwrRow,
    pub candidate: AdaptiveDwrRow,
    /// Areas on this round's coarse grid and its actual enrichment, respectively.
    pub baseline_areas: [f64; 2],
    pub candidate_areas: [f64; 2],
    pub coarse_gain: f64,
    pub enriched_gain: f64,
    /// Sum of weighted absolute CHANGES of signed case corrections, not a
    /// difference of unsigned indicator masses or a cancellable aggregate.
    pub correction_change: f64,
    /// Sum of weighted absolute changes in each case's coarse/enriched gain.
    pub observed_gain_change: f64,
    pub gain_allowance: f64,
    pub endpoints_resolved: bool,
    pub areas_match: bool,
    pub decreasing: bool,
    pub admissible: bool,
}

/// Initial binding and final complete pair. All final indicators use `grid`.
#[derive(Debug, Clone)]
pub struct AdaptiveGainReport {
    pub stop: AdaptiveGainStop,
    pub rows: Vec<AdaptiveGainRow>,
    pub grid: Quadtree,
    pub baseline_estimates: Vec<ElasticityDwrEstimate>,
    pub candidate_estimates: Vec<ElasticityDwrEstimate>,
    /// Baseline then candidate, in the original independent-case order.
    pub initial_case_compliances: [Vec<f64>; 2],
    pub snapshots: [u64; 2],
    pub solves: usize,
}

/// Interleave the original candidate search with actual adaptive pair solves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveDescentStage {
    Optimizer(crate::robust_descent::MultiLoadProjectedStage),
    Candidate { index: usize, stage: AdaptiveGainStage },
}

/// Complete comparison or numerical refusal for one coarse-admissible candidate.
#[derive(Debug, Clone)]
pub struct AdaptiveDescentCandidate {
    pub index: usize,
    pub report: Option<AdaptiveGainReport>,
    pub refusal: Option<String>,
}

/// Original optimizer outcome plus the comparisons made before publication.
#[derive(Debug, Clone)]
pub struct AdaptiveDescentUpdate {
    pub progress: crate::robust_descent::MultiLoadProjectedProgress,
    pub candidates: Vec<AdaptiveDescentCandidate>,
}

// Public estimate fields are data, not provenance. Only the actual producer
// below constructs rows used to admit an optimizer update.
fn changes(values: impl Iterator<Item = (f64, [f64; 3], [f64; 3])>)
    -> Result<(f64, f64), CutFemError> {
    let (mut correction, mut observed) = (0.0, 0.0);
    for (weight, old, new) in values {
        if !weight.is_finite() || weight < 0.0
            || old.iter().chain(&new).any(|v| !v.is_finite()) {
            return Err(invalid("invalid independent adaptive gain measurement"));
        }
        let delta = (old[2] - new[2]).abs();
        let change = ((old[1] - new[1]) - (old[0] - new[0])).abs();
        if !delta.is_finite() || !change.is_finite() {
            return Err(invalid("independent adaptive gain difference overflow"));
        }
        correction += weight * delta;
        observed += weight * change;
    }
    if !correction.is_finite() || !observed.is_finite() {
        return Err(invalid("independent adaptive gain accumulation overflow"));
    }
    Ok((correction, observed))
}

#[allow(clippy::too_many_arguments)]
fn row(baseline: AdaptiveDwrRow, candidate: AdaptiveDwrRow,
    baseline_areas: [f64; 2], candidate_areas: [f64; 2], correction_change: f64,
    observed_gain_change: f64, target: f64, decrease: f64, policy: AdaptiveGainPolicy)
    -> Result<AdaptiveGainRow, CutFemError> {
    let coarse_gain = baseline.coarse_compliance - candidate.coarse_compliance;
    let enriched_gain = baseline.enriched_compliance - candidate.enriched_compliance;
    let gain_allowance = policy.max_error_fraction * coarse_gain.max(0.0);
    if baseline.cells != candidate.cells || baseline.max_cell_level != candidate.max_cell_level
        || [coarse_gain, enriched_gain, gain_allowance, correction_change, observed_gain_change]
            .iter().any(|v| !v.is_finite())
        || baseline_areas.iter().chain(&candidate_areas).any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(invalid("adaptive comparison needs finite responses on identical grids"));
    }
    let endpoints_resolved = [&baseline, &candidate].iter().all(|r|
        r.correction_magnitude.max(r.observed_change) <= r.tolerance);
    let areas_match = baseline_areas.iter().zip(&candidate_areas).all(|(a, b)|
        (a - target).abs() <= policy.area_tolerance
            && (b - target).abs() <= policy.area_tolerance
            && (a - b).abs() <= policy.area_tolerance);
    let decreasing = candidate.coarse_compliance < baseline.coarse_compliance * (1.0 - decrease)
        && candidate.enriched_compliance < baseline.enriched_compliance * (1.0 - decrease);
    let admissible = endpoints_resolved && areas_match && decreasing
        && correction_change.max(observed_gain_change) <= gain_allowance;
    Ok(AdaptiveGainRow { baseline, candidate, baseline_areas, candidate_areas,
        coarse_gain, enriched_gain, correction_change, observed_gain_change,
        gain_allowance, endpoints_resolved, areas_match, decreasing, admissible })
}

/// Adapt the SAME analysis grid for two geometries until their gain is resolved
/// or explicit work/level/cell bounds stop the comparison. Each round solves
/// BOTH independent load families afresh, including zero-weight conditions.
/// Endpoint error screens and area checks on both grids accompany the gain test.
/// No fields, loads, weights, target areas or stress constraints are changed.
///
/// The caller must reserve `policy.max_solves(cases.len())` before entering and
/// charge two slots after every permitted `BeforeCase` callback, even on a later
/// failure. A break returns no partial report. Assembly and each two-solve case
/// are indivisible, not a hard wall-time guarantee.
#[allow(clippy::too_many_arguments)]
pub fn compare_adaptive_designs_controlled<B>(
    baseline: &GridSdf, candidate: &GridSdf, cases: &[RobustLoadCase], settings: OptimizeSettings,
    policy: AdaptiveGainPolicy, decrease: f64,
    mut control: impl FnMut(AdaptiveGainStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, AdaptiveGainReport>, CutFemError> {
    policy.validate(settings.level, cases.len())?;
    if !decrease.is_finite() || !(0.0..1.0).contains(&decrease)
        || !settings.volfrac.is_finite() || !(0.0..=1.0).contains(&settings.volfrac)
        || settings.volfrac == 0.0 {
        return Err(invalid("adaptive gain needs a relative decrease in [0,1) and target area in (0,1]"));
    }
    let material = prepare(baseline, cases, settings, policy.refinement)?;
    prepare(candidate, cases, settings, policy.refinement)?;
    let mut grid = Quadtree::with_room(settings.level, policy.refinement.max_level);
    let mut rows = Vec::new();
    rows.try_reserve_exact(policy.refinement.max_rounds)
        .map_err(|_| invalid("adaptive comparison history allocation refused"))?;
    let mut initial_case_compliances = [Vec::new(), Vec::new()];
    let mut solves = 0;
    for round in 0..policy.refinement.max_rounds {
        let mut estimates = [Vec::with_capacity(cases.len()), Vec::with_capacity(cases.len())];
        for (side, (design, geometry)) in [(ComparedDesign::Baseline, baseline),
            (ComparedDesign::Candidate, candidate)].into_iter().enumerate() {
            for (case_index, case) in cases.iter().enumerate() {
                if let ControlFlow::Break(reason) = control(AdaptiveGainStage::BeforeCase {
                    round, design, case: case_index,
                }) { return Ok(ControlFlow::Break(reason)); }
                solves += 2;
                estimates[side].push(estimate_case(&grid, geometry, &material, *case)?);
            }
        }
        if round == 0 {
            for side in 0..2 {
                initial_case_compliances[side] = estimates[side].iter().map(|e| e.j_primal).collect();
            }
        }
        let (old, mut indicators) = summarize(&grid, cases, &estimates[0], policy.refinement)?;
        let (new, other) = summarize(&grid, cases, &estimates[1], policy.refinement)?;
        for (cell, mass) in other {
            let combined = indicators.entry(cell).or_default();
            *combined += mass;
            if !combined.is_finite() { return Err(invalid("paired adaptive marking overflow")); }
        }
        let (correction, observed) = changes(cases.iter().zip(estimates[0].iter().zip(&estimates[1]))
            .map(|(case, (old, new))| (case.weight(), [old.j_primal, old.j_enriched, old.eta_signed],
                [new.j_primal, new.j_enriched, new.eta_signed])))?;
        let fine = grid.refined_once();
        let measured = row(old, new,
            [material_volume(&grid, baseline), material_volume(&fine, baseline)],
            [material_volume(&grid, candidate), material_volume(&fine, candidate)],
            correction, observed, settings.volfrac, decrease, policy)?;
        let resolved = measured.admissible;
        let non_improving = measured.endpoints_resolved && !measured.decreasing;
        rows.push(measured);
        if let ControlFlow::Break(reason) = control(AdaptiveGainStage::RoundComplete(round)) {
            return Ok(ControlFlow::Break(reason));
        }
        let stop = if resolved { Some(AdaptiveGainStop::ResolvedGain) }
            else if non_improving { Some(AdaptiveGainStop::NonImproving) }
            else if round + 1 == policy.refinement.max_rounds { Some(AdaptiveGainStop::RoundBudget) }
            else {
                let marked = mark(&grid, &indicators, policy.refinement.mark_fraction)?;
                if marked.is_empty() {
                    Some(if indicators.values().any(|v| *v > 0.0) {
                        AdaptiveGainStop::LevelLimit
                    } else { AdaptiveGainStop::NoIndicators })
                } else {
                    if let ControlFlow::Break(reason) = control(AdaptiveGainStage::BeforeRefine(round)) {
                        return Ok(ControlFlow::Break(reason));
                    }
                    match refine(&grid, &marked, policy.refinement.max_cells) {
                        Some(next) => { grid = next; None }
                        None => Some(AdaptiveGainStop::CellBudget),
                    }
                }
            };
        if let Some(stop) = stop {
            if let ControlFlow::Break(reason) = control(AdaptiveGainStage::BeforePublish) {
                return Ok(ControlFlow::Break(reason));
            }
            let [baseline_estimates, candidate_estimates] = estimates;
            return Ok(ControlFlow::Continue(AdaptiveGainReport { stop, rows, grid,
                baseline_estimates, candidate_estimates, initial_case_compliances,
                snapshots: [snapshot(baseline), snapshot(candidate)], solves }));
        }
    }
    Err(invalid("adaptive comparison ended without a complete pair"))
}

#[cfg(test)]
mod tests;
