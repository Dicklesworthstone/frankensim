//! Scenario-robust final evaluation for level-set topology designs.
//!
//! Every load case is solved independently against the exact same returned
//! geometry. Forces are never summed before equilibrium, so opposite scenarios
//! cannot cancel and hide structural demand. This module is an evaluation and
//! publication seam; it does not yet claim a simultaneous multi-load shape
//! gradient or global robust optimum.

use crate::{
    Cantilever, GridSdf, GuardedSettings, OptimizeReport, OptimizeSettings,
    material_volume, optimize_compliance_evaluated,
};
use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutFemError, CutStabilizationScaling,
    DesignBoxEdge, EdgeBand, MAX_PLANE_STRAIN_STIFFNESS_RATIO, Quadtree,
};
use fs_material::IsotropicElastic;

const MATERIAL_STRAIN_LIMIT: f64 = 1.0;
const SOLVER_TOL: f64 = 1e-12;
const SOLVER_MAX_ITERS: usize = 60_000;

/// One independent edge-traction scenario.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RobustLoadCase {
    edge: DesignBoxEdge,
    start: f64,
    end: f64,
    traction: [f64; 2],
    weight: f64,
}

impl RobustLoadCase {
    /// Declare one independent traction segment.
    ///
    /// The level-set optimizer's canonical support clamps the left edge, so a
    /// robust load on that same edge is refused instead of being silently
    /// interpreted as a reaction-only load.
    pub fn new(
        edge: DesignBoxEdge,
        start: f64,
        end: f64,
        traction: [f64; 2],
        weight: f64,
    ) -> Result<Self, CutFemError> {
        if edge == DesignBoxEdge::Left {
            return Err(invalid("robust traction cannot be applied on the clamped left edge"));
        }
        EdgeBand::new(edge, start, end)
            .map_err(|error| invalid(format!("invalid robust load segment: {error}")))?;
        if !traction.iter().all(|value| value.is_finite())
            || traction.iter().all(|value| *value == 0.0)
        {
            return Err(invalid("robust traction must be finite and nonzero"));
        }
        if !(weight.is_finite() && weight >= 0.0) {
            return Err(invalid("robust load weight must be finite and nonnegative"));
        }
        Ok(Self { edge, start, end, traction, weight })
    }

    /// Objective weight. Weights are not normalized implicitly.
    #[must_use]
    pub const fn weight(self) -> f64 { self.weight }

    /// Declared traction vector.
    #[must_use]
    pub const fn traction(self) -> [f64; 2] { self.traction }

    /// Supported design-box edge.
    #[must_use]
    pub const fn edge(self) -> DesignBoxEdge { self.edge }

    /// Inclusive support interval in the edge's natural coordinate.
    #[must_use]
    pub const fn interval(self) -> [f64; 2] { [self.start, self.end] }
}

/// How independent scenario compliances become one scalar selection objective.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobustAggregate {
    /// `sum_i weight_i * compliance_i`; weights are not normalized.
    WeightedSum,
    /// `max_i weight_i * compliance_i`.
    WorstWeightedCase,
}

/// Independently re-solved scenario metrics on one exact geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct RobustEvaluation {
    /// One unweighted compliance per input load case, in input order.
    pub case_compliances: Vec<f64>,
    /// `sum_i weight_i * compliance_i`.
    pub weighted_sum_compliance: f64,
    /// `max_i weight_i * compliance_i`.
    pub worst_weighted_compliance: f64,
    /// Selected scalar objective.
    pub objective: f64,
    /// Cut-quadrature material area.
    pub volume: f64,
    /// FNV-64 fingerprint of the exact nodal level-set bits.
    pub snapshot: u64,
}

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn snapshot(phi: &GridSdf) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for value in phi.nodes() {
        for byte in value.to_bits().to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

fn material(settings: OptimizeSettings) -> Result<IsotropicElastic, CutFemError> {
    if !(settings.youngs.is_finite() && settings.youngs > 0.0) {
        return Err(invalid("robust evaluation requires finite positive Young's modulus"));
    }
    if !(settings.poisson.is_finite() && settings.poisson > -1.0 && settings.poisson < 0.5) {
        return Err(invalid("robust evaluation requires Poisson ratio in (-1, 0.5)"));
    }
    let material = IsotropicElastic::new(settings.youngs, settings.poisson, MATERIAL_STRAIN_LIMIT)
        .map_err(|error| invalid(format!("robust material card refused: {error}")))?;
    let (lambda, mu) = material.lame();
    let ratio = (lambda + 2.0 * mu) / mu;
    if !(lambda.is_finite() && mu.is_finite() && mu > 0.0
        && ratio.is_finite() && ratio <= MAX_PLANE_STRAIN_STIFFNESS_RATIO)
    {
        return Err(invalid(format!(
            "robust plane-strain stiffness ratio {ratio} exceeds the admitted regime"
        )));
    }
    Ok(material)
}

fn validate_geometry(phi: &GridSdf, settings: OptimizeSettings) -> Result<(), CutFemError> {
    let expected = 1usize.checked_shl(settings.level)
        .ok_or_else(|| invalid("robust grid level exceeds platform address space"))?;
    if phi.n() != expected {
        return Err(invalid(format!(
            "robust SDF has {} cells per side but level {} requires {expected}",
            phi.n(), settings.level
        )));
    }
    if phi.nodes().iter().any(|value| !value.is_finite()) {
        return Err(invalid("robust evaluation requires finite level-set nodes"));
    }
    Ok(())
}

/// Independently solve every declared load case on one exact level-set geometry.
///
/// # Errors
/// Refuses malformed geometry/material/load declarations or any canonical
/// CutFEM elasticity failure.
pub fn evaluate_robust_design(
    phi: &GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    aggregate: RobustAggregate,
) -> Result<RobustEvaluation, CutFemError> {
    validate_geometry(phi, settings)?;
    if load_cases.is_empty() {
        return Err(invalid("robust evaluation requires at least one load case"));
    }
    if !load_cases.iter().any(|case| case.weight > 0.0) {
        return Err(invalid("robust evaluation requires at least one positive load weight"));
    }
    for case in load_cases {
        RobustLoadCase::new(case.edge, case.start, case.end, case.traction, case.weight)?;
    }

    let material = material(settings)?;
    let grid = Quadtree::uniform(settings.level);
    let clamp = |x: f64, _y: f64| x < 1e-9;
    let solver = CutElasticity {
        grid: &grid,
        sdf: phi,
        material: &material,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        stabilization_scaling: CutStabilizationScaling::LongitudinalModulus,
        quad_depth: 2,
        clamp: Some(&clamp),
        boundary_traction: None,
        traction_free_interface: true,
        solver_tol: SOLVER_TOL,
        solver_max_iters: SOLVER_MAX_ITERS,
    };

    let mut case_compliances = Vec::with_capacity(load_cases.len());
    let mut weighted_sum = 0.0;
    let mut worst = 0.0_f64;
    for case in load_cases {
        let support = EdgeBand::new(case.edge, case.start, case.end)
            .map_err(|error| invalid(format!("robust load segment refused: {error}")))?;
        let traction_value = case.traction;
        let traction = move |_: f64, _: f64| traction_value;
        let solution = solver.solve_with_boundary_traction(
            &|_, _| [0.0, 0.0],
            &|_, _| [0.0, 0.0],
            BoundaryTraction::EdgeBand { support, value: &traction },
        )?;
        let compliance = solution.compliance();
        if !(compliance.is_finite() && compliance >= 0.0) {
            return Err(invalid("robust load case produced invalid compliance"));
        }
        let weighted = case.weight * compliance;
        weighted_sum += weighted;
        worst = worst.max(weighted);
        case_compliances.push(compliance);
    }
    if !(weighted_sum.is_finite() && worst.is_finite()) {
        return Err(invalid("robust aggregate compliance is not finite"));
    }
    let volume = material_volume(&grid, phi);
    if !(volume.is_finite() && volume > 0.0) {
        return Err(invalid("robust design has invalid material area"));
    }
    let objective = match aggregate {
        RobustAggregate::WeightedSum => weighted_sum,
        RobustAggregate::WorstWeightedCase => worst,
    };
    Ok(RobustEvaluation {
        case_compliances,
        weighted_sum_compliance: weighted_sum,
        worst_weighted_compliance: worst,
        objective,
        volume,
        snapshot: snapshot(phi),
    })
}

/// One robustly re-evaluated candidate from the bounded move-size family.
#[derive(Debug, Clone)]
pub struct RobustCandidate {
    /// Zero-based deterministic candidate index.
    pub index: usize,
    /// Interface move limit used for the complete nominal-driver trajectory.
    pub move_cells: f64,
    /// Independent scenario evaluation, when both evolution and robust replay succeeded.
    pub evaluation: Option<RobustEvaluation>,
    /// Whether the cut-quadrature area met the declared upper limit.
    pub volume_feasible: bool,
    /// Whether the selected robust objective met the requested improvement gate.
    pub improvement_gate: bool,
    /// Typed evolution or robust-replay refusal retained for diagnostics.
    pub refusal: Option<String>,
}

/// Why robust bounded selection stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobustStop {
    /// A scenario-feasible, robustly improving candidate was published.
    Accepted,
    /// Every robustly evaluated candidate violated the area limit.
    NoFeasibleCandidate,
    /// Feasible candidates existed but none met the robust improvement gate.
    NoImprovingCandidate,
    /// Every candidate refused before robust metrics were available.
    AllCandidatesRefused,
}

/// Transactional robust-selection evidence.
#[derive(Debug, Clone)]
pub struct RobustOptimizeReport {
    /// Robust replay of the exact starting geometry.
    pub baseline: RobustEvaluation,
    /// Bounded candidate attempts in deterministic move-size order.
    pub candidates: Vec<RobustCandidate>,
    /// Search termination reason.
    pub stop: RobustStop,
    /// Robust replay of the published design, if one was accepted.
    pub accepted: Option<RobustEvaluation>,
    /// Nominal-driver evolution trajectory for the accepted geometry.
    pub trajectory: Option<OptimizeReport>,
    /// Move limit used by the accepted nominal-driver trajectory.
    pub accepted_move_cells: Option<f64>,
}

fn validate_guarded_controls(
    settings: OptimizeSettings,
    guarded: GuardedSettings,
) -> Result<(), CutFemError> {
    if settings.iterations == 0 {
        return Err(invalid("robust guarded optimization requires at least one evolution iteration"));
    }
    if !(settings.move_cells.is_finite() && settings.move_cells > 0.0) {
        return Err(invalid("robust guarded optimization requires finite positive move_cells"));
    }
    if !(settings.volfrac.is_finite() && settings.volfrac > 0.0 && settings.volfrac <= 1.0) {
        return Err(invalid("robust guarded optimization requires volfrac in (0, 1]"));
    }
    if guarded.max_candidates == 0 || guarded.max_candidates > 64 {
        return Err(invalid("robust guarded max_candidates must lie in [1, 64]"));
    }
    if !(guarded.contraction.is_finite()
        && guarded.contraction > 0.0
        && guarded.contraction < 1.0)
    {
        return Err(invalid("robust guarded contraction must lie strictly between zero and one"));
    }
    if !(guarded.volume_tolerance.is_finite()
        && guarded.volume_tolerance >= 0.0
        && guarded.volume_tolerance <= 1.0)
    {
        return Err(invalid("robust guarded volume_tolerance must lie in [0, 1]"));
    }
    if !(guarded.min_relative_improvement.is_finite()
        && guarded.min_relative_improvement >= 0.0
        && guarded.min_relative_improvement < 1.0)
    {
        return Err(invalid(
            "robust guarded min_relative_improvement must lie in [0, 1)",
        ));
    }
    Ok(())
}

/// Optimize with the existing nominal level-set descent, but publish only after
/// bounded candidates are independently replayed across every robust load case.
///
/// The optimizer itself is intentionally unchanged: `driver` supplies the one
/// nominal cantilever trajectory used to generate each candidate. Robust load
/// cases decide publication only after that trajectory finishes. This makes the
/// method a robust bounded selector, not a simultaneous multi-load shape-gradient
/// method.
///
/// Every candidate starts from the same exact input geometry and differs only by
/// `move_cells *= guarded.contraction`. The caller geometry is replaced only
/// when a robustly re-solved candidate satisfies the area and improvement gates.
///
/// # Errors
/// Refuses malformed controls, driver inputs, baseline robust evaluation, or
/// other admission failures. Individual candidate refusals are retained in the
/// returned report and do not mutate the caller geometry.
pub fn optimize_compliance_robust_guarded(
    phi: &mut GridSdf,
    driver: Cantilever,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    guarded: GuardedSettings,
    aggregate: RobustAggregate,
) -> Result<RobustOptimizeReport, CutFemError> {
    validate_guarded_controls(settings, guarded)?;
    let baseline = evaluate_robust_design(phi, load_cases, settings, aggregate)?;
    let limit = settings.volfrac + guarded.volume_tolerance;
    let baseline_feasible = baseline.volume <= limit;
    let improvement_limit =
        baseline.objective * (1.0 - guarded.min_relative_improvement);
    let origin = phi.clone();
    let mut candidates = Vec::with_capacity(guarded.max_candidates);
    let mut best: Option<(GridSdf, RobustEvaluation, OptimizeReport, f64)> = None;
    let mut move_cells = settings.move_cells;
    let mut any_final = false;
    let mut any_feasible = false;

    for index in 0..guarded.max_candidates {
        let mut trial_settings = settings;
        trial_settings.move_cells = move_cells;
        let mut trial = origin.clone();
        match optimize_compliance_evaluated(&mut trial, driver, trial_settings) {
            Ok(report) => match evaluate_robust_design(&trial, load_cases, settings, aggregate) {
                Ok(evaluation) => {
                    any_final = true;
                    let volume_feasible = evaluation.volume <= limit;
                    any_feasible |= volume_feasible;
                    let improvement_gate =
                        !baseline_feasible || evaluation.objective <= improvement_limit;
                    if volume_feasible && improvement_gate {
                        let replace = best.as_ref().is_none_or(|(_, current, _, _)| {
                            evaluation.objective < current.objective
                        });
                        if replace {
                            best = Some((
                                trial,
                                evaluation.clone(),
                                report.trajectory.clone(),
                                move_cells,
                            ));
                        }
                    }
                    candidates.push(RobustCandidate {
                        index,
                        move_cells,
                        evaluation: Some(evaluation),
                        volume_feasible,
                        improvement_gate,
                        refusal: None,
                    });
                }
                Err(error) => candidates.push(RobustCandidate {
                    index,
                    move_cells,
                    evaluation: None,
                    volume_feasible: false,
                    improvement_gate: false,
                    refusal: Some(format!("{error:?}")),
                }),
            },
            Err(error) => candidates.push(RobustCandidate {
                index,
                move_cells,
                evaluation: None,
                volume_feasible: false,
                improvement_gate: false,
                refusal: Some(format!("{error:?}")),
            }),
        }
        move_cells *= guarded.contraction;
    }

    if let Some((accepted_phi, accepted, trajectory, accepted_move_cells)) = best {
        *phi = accepted_phi;
        return Ok(RobustOptimizeReport {
            baseline,
            candidates,
            stop: RobustStop::Accepted,
            accepted: Some(accepted),
            trajectory: Some(trajectory),
            accepted_move_cells: Some(accepted_move_cells),
        });
    }

    let stop = if !any_final {
        RobustStop::AllCandidatesRefused
    } else if !any_feasible {
        RobustStop::NoFeasibleCandidate
    } else {
        RobustStop::NoImprovingCandidate
    };
    Ok(RobustOptimizeReport {
        baseline,
        candidates,
        stop,
        accepted: None,
        trajectory: None,
        accepted_move_cells: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn beam(level: u32) -> GridSdf {
        GridSdf::from_fn(1usize << level, &|_, y| (y - 0.5).abs() - 0.35)
    }

    fn settings(level: u32) -> OptimizeSettings {
        OptimizeSettings { level, iterations: 0, ..OptimizeSettings::default() }
    }

    #[test]
    fn opposite_independent_loads_do_not_cancel() {
        let phi = beam(3);
        let down = RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.5).unwrap();
        let up = RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 1.0], 0.5).unwrap();
        let one = evaluate_robust_design(&phi, &[down], settings(3), RobustAggregate::WeightedSum).unwrap();
        let both = evaluate_robust_design(&phi, &[down, up], settings(3), RobustAggregate::WeightedSum).unwrap();
        assert!(one.case_compliances[0] > 0.0);
        assert!((both.case_compliances[0] - both.case_compliances[1]).abs()
            < 1e-10 * both.case_compliances[0].max(1.0));
        assert!((both.weighted_sum_compliance - 2.0 * one.weighted_sum_compliance).abs()
            < 1e-10 * both.weighted_sum_compliance.max(1.0));
    }

    #[test]
    fn worst_case_and_weighted_sum_are_distinct_declared_aggregates() {
        let phi = beam(3);
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.45, 0.55, [0.0, -2.0], 0.25).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Top, 0.75, 1.0, [1.0, 0.0], 0.75).unwrap(),
        ];
        let sum = evaluate_robust_design(&phi, &cases, settings(3), RobustAggregate::WeightedSum).unwrap();
        let worst = evaluate_robust_design(&phi, &cases, settings(3), RobustAggregate::WorstWeightedCase).unwrap();
        assert_eq!(sum.case_compliances, worst.case_compliances);
        assert!(sum.objective >= worst.objective);
        assert_eq!(sum.weighted_sum_compliance, sum.objective);
        assert_eq!(worst.worst_weighted_compliance, worst.objective);
        assert_eq!(sum.snapshot, worst.snapshot);
        assert_eq!(sum.volume, worst.volume);
    }

    #[test]
    fn malformed_load_sets_refuse() {
        let phi = beam(3);
        assert!(evaluate_robust_design(&phi, &[], settings(3), RobustAggregate::WeightedSum).is_err());
        let zero_weight = RobustLoadCase::new(
            DesignBoxEdge::Right, 0.4, 0.6, [0.0, -1.0], 0.0
        ).unwrap();
        assert!(evaluate_robust_design(
            &phi, &[zero_weight], settings(3), RobustAggregate::WeightedSum
        ).is_err());
        assert!(RobustLoadCase::new(
            DesignBoxEdge::Left, 0.4, 0.6, [0.0, -1.0], 1.0
        ).is_err());
        assert!(RobustLoadCase::new(
            DesignBoxEdge::Right, 0.4, 0.6, [0.0, 0.0], 1.0
        ).is_err());
    }

    #[test]
    fn robust_guard_refusal_does_not_mutate_caller_geometry() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Top, 0.75, 1.0, [0.5, 0.0], 1.0).unwrap(),
        ];
        let settings = OptimizeSettings {
            level: 3,
            iterations: 1,
            move_cells: 1e-9,
            volfrac: 1.0,
            nucleation_period: 0,
            ..OptimizeSettings::default()
        };
        let guarded = GuardedSettings {
            max_candidates: 1,
            contraction: 0.5,
            volume_tolerance: 0.0,
            min_relative_improvement: 0.5,
        };
        let report = optimize_compliance_robust_guarded(
            &mut phi,
            Cantilever { load: 1.0, band: 0.125 },
            &cases,
            settings,
            guarded,
            RobustAggregate::WorstWeightedCase,
        ).unwrap();
        assert_ne!(report.stop, RobustStop::Accepted);
        assert!(report.accepted.is_none());
        assert_eq!(phi.nodes(), before.as_slice());
        assert_eq!(report.candidates.len(), 1);
    }

    #[test]
    fn invalid_robust_guard_controls_refuse_without_mutation() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.4, 0.6, [0.0, -1.0], 1.0).unwrap(),
        ];
        let settings = OptimizeSettings {
            level: 3,
            iterations: 1,
            ..OptimizeSettings::default()
        };
        let guarded = GuardedSettings {
            max_candidates: 0,
            ..GuardedSettings::default()
        };
        assert!(optimize_compliance_robust_guarded(
            &mut phi,
            Cantilever { load: 1.0, band: 0.125 },
            &cases,
            settings,
            guarded,
            RobustAggregate::WeightedSum,
        ).is_err());
        assert_eq!(phi.nodes(), before.as_slice());
    }
}
