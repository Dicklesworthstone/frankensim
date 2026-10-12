//! Goal-driven local refinement of the BACKGROUND grid, not the material design.
//!
//! Every round solves every load on a 2:1-balanced grid and its one-level
//! enrichment through fs-dwr's canonical vector estimator. The original bilinear
//! SDF, loads, material and boundary conditions never change. Signed estimates
//! are retained, but stopping and marking use non-cancelling weighted magnitudes.
//! This is an Estimated resolution study, not a continuum or stress certificate.

use std::collections::{BTreeMap, BTreeSet};
use std::ops::ControlFlow;

use crate::{GridSdf, OptimizeSettings, RobustLoadCase, material_volume};
use fs_cutfem::{BoundaryTraction, CellKey, CutElasticity, CutFemError,
    CutStabilizationScaling, EdgeBand, MAX_PLANE_STRAIN_STIFFNESS_RATIO, Quadtree};
use fs_dwr::ElasticityDwrEstimate;
use fs_material::IsotropicElastic;

/// Shared-grid local refinement of a proposed design's compliance gain.
pub mod gain;

/// Work/accuracy controls. Cell caps cover the entire background, including void.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveDwrPolicy {
    /// Highest COARSE cell level. Estimator enrichment uses one further level.
    pub max_level: u32,
    /// Maximum complete all-load, two-solve estimation rounds, in 1..=16.
    pub max_rounds: usize,
    /// Maximum coarse leaves, in 4..=4096; enrichment has exactly four times as many.
    pub max_cells: usize,
    /// Absolute tolerance in the declared, unnormalized weighted objective's J.
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    /// Fraction of total indicator mass to mark, in (0,1].
    pub mark_fraction: f64,
}

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

impl AdaptiveDwrPolicy {
    /// Admission happens before grid allocation or any physical solve.
    pub fn validate(self, initial_level: u32, cases: usize) -> Result<(), CutFemError> {
        if !(1..=5).contains(&initial_level) || !(initial_level..=6).contains(&self.max_level)
            || !(1..=16).contains(&cases) || !(1..=16).contains(&self.max_rounds)
            || !(4..=4096).contains(&self.max_cells)
            || self.max_cells < (1usize << (2 * initial_level))
            || !self.absolute_tolerance.is_finite() || self.absolute_tolerance < 0.0
            || !self.relative_tolerance.is_finite() || !(0.0..1.0).contains(&self.relative_tolerance)
            || (self.absolute_tolerance == 0.0 && self.relative_tolerance == 0.0)
            || !self.mark_fraction.is_finite() || self.mark_fraction <= 0.0 || self.mark_fraction > 1.0
        {
            return Err(invalid("adaptive DWR needs initial level 1..=5, max level through 6, 1..=16 cases/rounds, 4..=4096 cells covering the initial grid, finite nonnegative tolerances (one positive, relative <1) and marking fraction in (0,1]"));
        }
        Ok(())
    }
    /// Worst-case solve allowance. Validate first; this method is overflow-safe.
    pub fn max_solves(self, cases: usize) -> Option<usize> {
        self.max_rounds.checked_mul(cases)?.checked_mul(2)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveDwrStop { EstimatedTolerance, RoundBudget, CellBudget, LevelLimit, NoIndicators }
impl AdaptiveDwrStop {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::EstimatedTolerance => "estimated-tolerance-met", Self::RoundBudget => "round-budget",
            Self::CellBudget => "cell-budget", Self::LevelLimit => "level-limit", Self::NoIndicators => "no-indicators",
        }
    }
}

/// Cancellation boundaries. Each case's two solves and integration remain indivisible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveDwrStage {
    BeforeCase { round: usize, case: usize },
    RoundComplete(usize),
    BeforeRefine(usize),
    BeforePublish,
}

/// Actual measurements from one complete load family, before any next refinement.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveDwrRow {
    pub cells: usize,
    pub max_cell_level: u32,
    pub coarse_compliance: f64,
    pub enriched_compliance: f64,
    pub eta_signed: f64,
    /// Sum of weighted absolute signed corrections, NOT a guaranteed upper bound.
    pub correction_magnitude: f64,
    /// Sum of weighted absolute coarse/enriched compliance differences.
    pub observed_change: f64,
    pub tolerance: f64,
    pub coarse_dofs: usize,
    pub enriched_dofs: usize,
    /// Largest recomputed Euclidean residual among the independent solves.
    pub relative_residual: f64,
}

/// Last COMPLETE grid/load family and the complete numerical history.
#[derive(Debug, Clone)]
pub struct AdaptiveDwrReport {
    pub stop: AdaptiveDwrStop,
    pub rows: Vec<AdaptiveDwrRow>,
    pub grid: Quadtree,
    /// In declared case order. The signed indicators belong to `grid` only.
    pub estimates: Vec<ElasticityDwrEstimate>,
    /// Unweighted compliances on the ORIGINAL grid, in declared case order.
    /// Retained separately so a consumer can bind every original response even
    /// after refinement; equal weighted sums must not conceal changed cases.
    pub initial_case_compliances: Vec<f64>,
    /// Geometry identity and area on the ORIGINAL optimization grid.
    pub snapshot: u64,
    pub original_area: f64,
    /// Cut-quadrature area on the last adaptive grid. Never substituted into the optimizer.
    pub assessed_area: f64,
    pub solves: usize,
}

// Same deterministic signed-field fingerprint used by the existing evaluators.
fn snapshot(phi: &GridSdf) -> u64 {
    phi.nodes().iter().flat_map(|x| x.to_bits().to_le_bytes()).fold(
        0xcbf2_9ce4_8422_2325, |h, b| (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3))
}

// Shared unchanged physical admission for final assessment and paired gain work.
fn prepare(phi: &GridSdf, cases: &[RobustLoadCase], settings: OptimizeSettings,
    policy: AdaptiveDwrPolicy) -> Result<IsotropicElastic, CutFemError> {
    policy.validate(settings.level, cases.len())?;
    if phi.n() != 1usize << settings.level || phi.nodes().iter().any(|v| !v.is_finite())
        || !settings.youngs.is_finite() || settings.youngs <= 0.0
        || !settings.poisson.is_finite() || settings.poisson <= -1.0 || settings.poisson >= 0.5
        || !cases.iter().any(|c| c.weight() > 0.0)
    { return Err(invalid("adaptive DWR requires a finite level-matched SDF and a positive objective weight")); }
    for case in cases {
        let [a, b] = case.interval();
        RobustLoadCase::new(case.edge(), a, b, case.traction(), case.weight())?;
    }
    let material = IsotropicElastic::new(settings.youngs, settings.poisson, 1.0)
        .map_err(|e| invalid(format!("adaptive material: {e}")))?;
    let (lambda, mu) = material.lame();
    let ratio = (lambda + 2.0 * mu) / mu;
    if !lambda.is_finite() || !mu.is_finite() || mu <= 0.0
        || !ratio.is_finite() || ratio > MAX_PLANE_STRAIN_STIFFNESS_RATIO {
        return Err(invalid("adaptive DWR material is outside the canonical plane-strain regime"));
    }
    Ok(material)
}

// One canonical problem definition for both adaptive consumers. Every estimate
// still calls the original vector CutFEM/DWR implementation, never a new solver.
fn estimate_case(grid: &Quadtree, phi: &GridSdf, material: &IsotropicElastic,
    case: RobustLoadCase) -> Result<ElasticityDwrEstimate, CutFemError> {
    let clamp = |x: f64, _: f64| x < 1e-9;
    let problem = CutElasticity {
        grid, sdf: phi, material, nitsche_beta: 20.0, ghost_gamma: 0.5,
        stabilization_scaling: CutStabilizationScaling::LongitudinalModulus,
        quad_depth: 2, clamp: Some(&clamp), boundary_traction: None,
        traction_free_interface: true, solver_tol: 1e-12, solver_max_iters: 60_000,
    };
    let [a, b] = case.interval();
    let support = EdgeBand::new(case.edge(), a, b).map_err(|e| invalid(e.to_string()))?;
    let traction = |_: f64, _: f64| case.traction();
    fs_dwr::estimate_elasticity_compliance_with_boundary_traction(
        &problem, &|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0],
        BoundaryTraction::EdgeBand { support, value: &traction })
}

/// Locally refine until measured goal sensitivity fits the requested tolerance or
/// an explicit bound stops progress. No geometry projection, redistancing, load
/// aggregation before equilibrium, or interpolated displacement is performed.
///
/// A break returns no partial report and never changes the input. The caller
/// must charge interrupted work and retain its accepted design independently.
/// At most 2*cases*max_rounds solves are attempted, each capped at 60000 iterations.
/// Grids are capped BEFORE each split, including all 2:1 balance closure splits.
/// This is goal-driven background refinement, not optimization or a stress check.
pub fn assess_adaptive_dwr_controlled<B>(
    phi: &GridSdf, cases: &[RobustLoadCase], settings: OptimizeSettings,
    policy: AdaptiveDwrPolicy, mut control: impl FnMut(AdaptiveDwrStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, AdaptiveDwrReport>, CutFemError> {
    let material = prepare(phi, cases, settings, policy)?;
    let mut grid = Quadtree::with_room(settings.level, policy.max_level);
    let original_area = material_volume(&grid, phi);
    if !original_area.is_finite() || original_area <= 0.0 { return Err(invalid("invalid original area")); }
    let mut rows = Vec::new();
    rows.try_reserve_exact(policy.max_rounds).map_err(|_| invalid("adaptive history allocation refused"))?;
    let mut initial_case_compliances = Vec::new();
    let mut solves = 0usize;
    for round in 0..policy.max_rounds {
        let mut estimates = Vec::with_capacity(cases.len());
        for (case_index, case) in cases.iter().enumerate() {
            if let ControlFlow::Break(reason) = control(AdaptiveDwrStage::BeforeCase { round, case: case_index }) {
                return Ok(ControlFlow::Break(reason));
            }
            solves += 2;
            estimates.push(estimate_case(&grid, phi, &material, *case)?);
        }
        let (row, indicators) = summarize(&grid, cases, &estimates, policy)?;
        if round == 0 {
            initial_case_compliances = estimates.iter().map(|e| e.j_primal).collect();
        }
        let tolerance_met = row.correction_magnitude.max(row.observed_change) <= row.tolerance;
        rows.push(row);
        if let ControlFlow::Break(reason) = control(AdaptiveDwrStage::RoundComplete(round)) {
            return Ok(ControlFlow::Break(reason));
        }
        let stop = if tolerance_met { Some(AdaptiveDwrStop::EstimatedTolerance) }
            else if round + 1 == policy.max_rounds { Some(AdaptiveDwrStop::RoundBudget) }
            else {
                let marked = mark(&grid, &indicators, policy.mark_fraction)?;
                if marked.is_empty() {
                    Some(if indicators.values().any(|v| *v > 0.0) {
                        AdaptiveDwrStop::LevelLimit
                    } else { AdaptiveDwrStop::NoIndicators })
                } else {
                    if let ControlFlow::Break(reason) = control(AdaptiveDwrStage::BeforeRefine(round)) {
                        return Ok(ControlFlow::Break(reason));
                    }
                    match refine(&grid, &marked, policy.max_cells) {
                        Some(next) => { grid = next; None }
                        None => Some(AdaptiveDwrStop::CellBudget),
                    }
                }
            };
        if let Some(stop) = stop {
            let assessed_area = material_volume(&grid, phi);
            if !assessed_area.is_finite() || assessed_area <= 0.0 {
                return Err(invalid("invalid adaptive-grid numerical area"));
            }
            if let ControlFlow::Break(reason) = control(AdaptiveDwrStage::BeforePublish) {
                return Ok(ControlFlow::Break(reason));
            }
            return Ok(ControlFlow::Continue(AdaptiveDwrReport {
                stop, rows, grid, estimates, initial_case_compliances,
                snapshot: snapshot(phi), original_area, assessed_area, solves,
            }));
        }
    }
    Err(invalid("adaptive DWR exhausted without a terminal measurement"))
}

fn summarize(grid: &Quadtree, cases: &[RobustLoadCase], estimates: &[ElasticityDwrEstimate],
    policy: AdaptiveDwrPolicy) -> Result<(AdaptiveDwrRow, BTreeMap<CellKey, f64>), CutFemError> {
    let mut row = AdaptiveDwrRow { cells: grid.leaf_count(),
        max_cell_level: grid.leaves().map(|c| c.0).max().unwrap_or(0),
        coarse_compliance: 0.0, enriched_compliance: 0.0, eta_signed: 0.0,
        correction_magnitude: 0.0, observed_change: 0.0, tolerance: 0.0,
        coarse_dofs: 0, enriched_dofs: 0, relative_residual: 0.0 };
    let mut indicators = BTreeMap::<CellKey, f64>::new();
    for (case, e) in cases.iter().zip(estimates) {
        if [e.j_primal, e.j_enriched, e.eta_signed, e.primal_relative_residual, e.enriched_relative_residual]
            .iter().any(|v| !v.is_finite()) || e.j_primal < 0.0 || e.j_enriched < 0.0
            || e.primal_relative_residual < 0.0 || e.enriched_relative_residual < 0.0
            || e.primal_relative_residual >= 1e-12 || e.enriched_relative_residual >= 1e-12 {
            return Err(invalid("adaptive DWR received invalid values or unconverged independent solves"));
        }
        let w = case.weight();
        row.coarse_compliance += w * e.j_primal;
        row.enriched_compliance += w * e.j_enriched;
        row.eta_signed += w * e.eta_signed;
        row.correction_magnitude += w * e.eta_signed.abs();
        row.observed_change += w * (e.j_enriched - e.j_primal).abs();
        row.coarse_dofs += e.dofs;
        row.enriched_dofs += e.enriched_dofs;
        row.relative_residual = row.relative_residual.max(e.primal_relative_residual).max(e.enriched_relative_residual);
        for (&cell, &eta) in &e.indicators {
            if !grid.is_leaf(cell) || !eta.is_finite() { return Err(invalid("invalid adaptive cell indicator")); }
            *indicators.entry(cell).or_default() += w * eta.abs();
        }
    }
    row.tolerance = policy.absolute_tolerance + policy.relative_tolerance * row.enriched_compliance.abs();
    if [row.coarse_compliance, row.enriched_compliance, row.eta_signed,
        row.correction_magnitude, row.observed_change, row.tolerance].iter().any(|v| !v.is_finite())
        || indicators.values().any(|v| !v.is_finite()) {
        return Err(invalid("adaptive DWR aggregation overflow"));
    }
    Ok((row, indicators))
}

// Stable bulk marking. Saturated cells retain their mass but cannot split.
fn mark(grid: &Quadtree, indicators: &BTreeMap<CellKey, f64>, fraction: f64)
    -> Result<Vec<CellKey>, CutFemError> {
    let total: f64 = indicators.values().sum();
    if !total.is_finite() { return Err(invalid("adaptive marking mass overflow")); }
    let mut ranked: Vec<_> = indicators.iter().filter(|(c, v)| c.0 < grid.max_level() && **v > 0.0)
        .map(|(&c, &v)| (c, v)).collect();
    ranked.sort_by(|(a, x), (b, y)| y.total_cmp(x).then_with(|| a.cmp(b)));
    let mut sum = 0.0;
    let mut marked = Vec::new();
    for (cell, mass) in ranked {
        marked.push(cell); sum += mass;
        if sum >= fraction * total { break; }
    }
    Ok(marked)
}

// Same edge-balance rule as Quadtree::balance, with a pre-split cap. All work
// occurs in a disposable copy; an over-budget balance closure changes nothing.
fn refine(grid: &Quadtree, marked: &[CellKey], max_cells: usize) -> Option<Quadtree> {
    let mut next = grid.clone();
    let split = |grid: &mut Quadtree, cell: CellKey| -> Option<()> {
        if !grid.is_leaf(cell) || cell.0 >= grid.max_level() || grid.leaf_count().checked_add(3)? > max_cells {
            return None;
        }
        grid.split(cell); Some(())
    };
    for &cell in marked { split(&mut next, cell)?; }
    loop {
        let mut due = BTreeSet::new();
        for cell in next.leaves() {
            for dir in 0..4 {
                if let Some(neighbor) = next.covering_neighbor(cell, dir)
                    && neighbor.0 + 1 < cell.0 { due.insert(neighbor); }
            }
        }
        if due.is_empty() { return Some(next); }
        for cell in due { split(&mut next, cell)?; }
    }
}

#[cfg(test)]
mod tests;
