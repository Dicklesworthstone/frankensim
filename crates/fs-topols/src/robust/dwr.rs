//! Independent-load DWR on one unchanged final geometry. Linear objective
//! weights are applied AFTER each equilibrium and residual assessment.

use super::*;
use fs_dwr::{ElasticityDwrEstimate, ElasticityResidualTerms};
use std::convert::Infallible;
use std::ops::ControlFlow;

/// One independently solved operating condition, including zero-weight cases.
#[derive(Debug, Clone)]
pub struct LoadCaseDwrAssessment {
    pub load: RobustLoadCase,
    pub estimate: ElasticityDwrEstimate,
}

/// Estimated discretization error in `sum_i weight_i * compliance_i`.
/// This is not an estimate of sampled/continuous stress or a certified bound.
#[derive(Debug, Clone)]
pub struct WeightedComplianceDwrAssessment {
    pub snapshot: u64,
    pub volume: f64,
    pub level: u32,
    pub cases: Vec<LoadCaseDwrAssessment>,
    pub coarse_compliance: f64,
    pub enriched_compliance: f64,
    pub eta_signed: f64,
    /// Sum of weighted absolute indicator masses; no cross-case cancellation.
    /// Marking evidence only, not a continuum-error interval radius.
    pub absolute_indicator_sum: f64,
    pub terms: ElasticityResidualTerms,
}

/// Boundaries around each indivisible two-solve case and final publication.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WeightedComplianceDwrStage {
    BeforeCase { case: usize },
    BeforePublish,
}

/// Assess every load independently on the original and one-level enriched grid.
///
/// # Errors
/// Admits 1..=16 cases and levels 1..=5, valid geometry/materials and at least
/// one positive weight. Propagates any per-case solve or estimator refusal.
pub fn assess_weighted_compliance_dwr(
    phi: &GridSdf,
    cases: &[RobustLoadCase],
    settings: OptimizeSettings,
) -> Result<WeightedComplianceDwrAssessment, CutFemError> {
    match assess_weighted_compliance_dwr_controlled(phi, cases, settings, |_| {
        ControlFlow::<Infallible>::Continue(())
    })? {
        ControlFlow::Continue(assessment) => Ok(assessment),
        ControlFlow::Break(never) => match never {},
    }
}

/// Cancellation between independently assessed cases; a break publishes no
/// partial family. Each case has two solves capped at 60,000 CG iterations and
/// the canonical 1e-12 recomputed Euclidean residual gate. Assembly, those
/// solves and residual integration are indivisible within a case. The original
/// bilinear field is sampled unchanged on both grids, never reprojected.
///
/// # Errors
/// Same admission and numerical refusals as [`assess_weighted_compliance_dwr`].
pub fn assess_weighted_compliance_dwr_controlled<B>(
    phi: &GridSdf,
    cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    mut control: impl FnMut(WeightedComplianceDwrStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, WeightedComplianceDwrAssessment>, CutFemError> {
    if !(1..=5).contains(&settings.level) || !(1..=16).contains(&cases.len()) {
        return Err(invalid("weighted DWR requires levels 1..=5 and 1..=16 independent cases"));
    }
    validate_geometry(phi, settings)?;
    if !cases.iter().any(|case| case.weight > 0.0) {
        return Err(invalid("weighted DWR requires a positive objective weight"));
    }
    for case in cases {
        RobustLoadCase::new(case.edge, case.start, case.end, case.traction, case.weight)?;
    }
    let material = material(settings)?;
    let grid = Quadtree::uniform(settings.level);
    let clamp = |x: f64, _: f64| x < 1e-9;
    let problem = CutElasticity {
        grid: &grid, sdf: phi, material: &material,
        nitsche_beta: 20.0, ghost_gamma: 0.5,
        stabilization_scaling: CutStabilizationScaling::LongitudinalModulus,
        quad_depth: 2, clamp: Some(&clamp), boundary_traction: None,
        traction_free_interface: true, solver_tol: SOLVER_TOL,
        solver_max_iters: SOLVER_MAX_ITERS,
    };
    let mut result = WeightedComplianceDwrAssessment {
        snapshot: snapshot(phi), volume: material_volume(&grid, phi), level: settings.level,
        cases: Vec::with_capacity(cases.len()), coarse_compliance: 0.0,
        enriched_compliance: 0.0, eta_signed: 0.0, absolute_indicator_sum: 0.0,
        terms: ElasticityResidualTerms::default(),
    };
    if !(result.volume.is_finite() && result.volume > 0.0) {
        return Err(invalid("weighted DWR requires positive finite material area"));
    }
    for (index, &case) in cases.iter().enumerate() {
        if let ControlFlow::Break(reason) = control(WeightedComplianceDwrStage::BeforeCase { case: index }) {
            return Ok(ControlFlow::Break(reason));
        }
        let support = EdgeBand::new(case.edge, case.start, case.end)
            .map_err(|error| invalid(format!("DWR load segment refused: {error}")))?;
        let traction = |_: f64, _: f64| case.traction;
        let estimate = fs_dwr::estimate_elasticity_compliance_with_boundary_traction(
            &problem, &|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0],
            BoundaryTraction::EdgeBand { support, value: &traction },
        )?;
        let w = case.weight;
        result.coarse_compliance += w * estimate.j_primal;
        result.enriched_compliance += w * estimate.j_enriched;
        result.eta_signed += w * estimate.eta_signed;
        result.absolute_indicator_sum += w * estimate.eta_abs;
        result.terms.bulk += w * estimate.terms.bulk;
        result.terms.nitsche += w * estimate.terms.nitsche;
        result.terms.outer_traction += w * estimate.terms.outer_traction;
        result.terms.ghost += w * estimate.terms.ghost;
        result.cases.push(LoadCaseDwrAssessment { load: case, estimate });
    }
    if [result.coarse_compliance, result.enriched_compliance, result.eta_signed,
        result.absolute_indicator_sum, result.terms.bulk, result.terms.nitsche,
        result.terms.outer_traction, result.terms.ghost, result.terms.total()]
        .iter().any(|value| !value.is_finite()) {
        return Err(invalid("weighted DWR aggregate is not finite"));
    }
    if let ControlFlow::Break(reason) = control(WeightedComplianceDwrStage::BeforePublish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue(result))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case(traction: [f64; 2], weight: f64) -> RobustLoadCase {
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, traction, weight).unwrap()
    }

    #[test]
    fn independent_weighted_dwr_matches_actual_objective_and_loads_cannot_cancel() {
        let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        let before = phi.nodes().to_vec();
        let settings = OptimizeSettings { level: 3, ..OptimizeSettings::default() };
        let cases = [case([0.0, -1.0], 0.25), case([0.0, 2.0], 0.75), case([1.0, 0.0], 0.0)];
        let actual = crate::evaluate_robust_sampled_stress(&phi, &cases, settings,
            RobustAggregate::WeightedSum).unwrap();
        let dwr = assess_weighted_compliance_dwr(&phi, &cases, settings).unwrap();
        assert_eq!(dwr.snapshot, actual.snapshot);
        assert_eq!(dwr.volume.to_bits(), actual.volume.to_bits());
        assert_eq!(dwr.coarse_compliance.to_bits(), actual.objective.to_bits());
        for (case, compliance) in dwr.cases.iter().zip(actual.case_compliances) {
            assert_eq!(case.estimate.j_primal.to_bits(), compliance.to_bits());
            assert!(case.estimate.primal_relative_residual < SOLVER_TOL);
            assert!(case.estimate.enriched_relative_residual < SOLVER_TOL);
        }
        let first = &dwr.cases[0].estimate;
        let second = &dwr.cases[1].estimate;
        for (a, b) in [(first.j_primal, second.j_primal), (first.j_enriched, second.j_enriched),
            (first.eta_signed, second.eta_signed), (first.eta_abs, second.eta_abs)] {
            assert!((4.0 * a - b).abs() < 1e-8 * b.abs().max(1e-8));
        }
        assert_eq!(dwr.eta_signed.to_bits(), (0.25 * first.eta_signed + 0.75 * second.eta_signed).to_bits());
        assert!((dwr.eta_signed - dwr.terms.total()).abs() < 1e-10);
        assert!(dwr.absolute_indicator_sum >= dwr.eta_signed.abs() - 1e-12);
        assert!(dwr.cases[2].estimate.j_primal > 0.0, "zero weight retains independent evidence");
        assert_eq!(phi.nodes(), before);
    }

    #[test]
    fn phase_cancellation_returns_no_partial_family_and_retains_geometry() {
        let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        let before = phi.nodes().to_vec();
        let cases = [case([0.0, -1.0], 1.0), case([0.0, 1.0], 1.0)];
        for stop in [WeightedComplianceDwrStage::BeforeCase { case: 0 },
            WeightedComplianceDwrStage::BeforeCase { case: 1 }, WeightedComplianceDwrStage::BeforePublish] {
            let result = assess_weighted_compliance_dwr_controlled(&phi, &cases,
                OptimizeSettings { level: 3, ..OptimizeSettings::default() }, |stage| {
                    if stage == stop { ControlFlow::Break(stage) } else { ControlFlow::Continue(()) }
                }).unwrap();
            assert!(matches!(result, ControlFlow::Break(stage) if stage == stop));
            assert_eq!(phi.nodes(), before);
        }
    }
}
