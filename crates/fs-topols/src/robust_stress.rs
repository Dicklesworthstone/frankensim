//! Scenario-robust sampled-stress evaluation and publication.
//!
//! Every declared load case is solved independently on the exact same geometry.
//! Compliance aggregation follows [`RobustAggregate`], while stress admission
//! uses the worst deterministic sampled plane-strain von Mises value across all
//! cases, including zero-objective-weight cases. Candidate generation uses the
//! simultaneous multi-load descent and final publication independently replays
//! every scenario. This is sample-scoped stress admission, not a continuous
//! maximum-stress certificate or a stress-adjoint/KKT optimum.

use crate::{
    GridSdf, GuardedSettings, OptimizeReport, OptimizeSettings, RobustAggregate,
    RobustLoadCase, SampledStressLimit, optimize_compliance_multi_load,
};
use fs_cutfem::quad::cut_cell_rules;
use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutElasticitySolution, CutFemError, CutSdf,
    CutStabilizationScaling, EdgeBand, MAX_PLANE_STRAIN_STIFFNESS_RATIO, Quadtree,
};
use fs_material::IsotropicElastic;
use std::convert::Infallible;
use std::ops::ControlFlow;

const MATERIAL_STRAIN_LIMIT: f64 = 1.0;
const SOLVER_TOL: f64 = 1e-12;
const SOLVER_MAX_ITERS: usize = 60_000;
const GAUSS: f64 = 0.577_350_269_189_625_8;

/// Independent robust mechanics plus sampled-stress evidence on one geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct RobustSampledStressEvaluation {
    /// Unweighted compliance for every declared scenario, input order.
    pub case_compliances: Vec<f64>,
    /// Sampled maximum von Mises stress for every scenario, input order.
    pub case_sampled_max_von_mises: Vec<f64>,
    /// First deterministic location attaining each scenario maximum.
    pub case_max_locations: Vec<[f64; 2]>,
    /// Number of retained material stress probes in each scenario.
    pub case_sample_counts: Vec<usize>,
    /// Weighted-sum compliance across scenarios.
    pub weighted_sum_compliance: f64,
    /// Worst weighted compliance across scenarios.
    pub worst_weighted_compliance: f64,
    /// Compliance objective selected by the declared aggregate.
    pub objective: f64,
    /// Worst sampled von Mises stress across all declared scenarios.
    pub worst_sampled_von_mises: f64,
    /// Input index of the first scenario attaining the worst sampled stress.
    pub worst_stress_case: usize,
    /// Cut-quadrature material area.
    pub volume: f64,
    /// FNV-64 fingerprint of exact level-set node bits.
    pub snapshot: u64,
}

/// One robust stress-limited bounded-search candidate.
#[derive(Debug, Clone)]
pub struct RobustStressCandidate {
    pub index: usize,
    pub move_cells: f64,
    pub evaluation: Option<RobustSampledStressEvaluation>,
    pub volume_feasible: bool,
    pub stress_feasible: bool,
    pub improvement_gate: bool,
    pub refusal: Option<String>,
}

/// Why robust sampled-stress publication stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobustStressStop {
    Accepted,
    NoConstraintFeasibleCandidate,
    NoImprovingCandidate,
    AllCandidatesRefused,
}

/// Transactional bounded-search evidence for robust stress-limited design.
#[derive(Debug, Clone)]
pub struct RobustStressReport {
    pub baseline: RobustSampledStressEvaluation,
    pub limit: SampledStressLimit,
    pub candidates: Vec<RobustStressCandidate>,
    pub stop: RobustStressStop,
    pub accepted: Option<RobustSampledStressEvaluation>,
    pub trajectory: Option<OptimizeReport>,
    pub accepted_move_cells: Option<f64>,
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

fn material(settings: OptimizeSettings) -> Result<(IsotropicElastic, f64, f64), CutFemError> {
    if !(settings.youngs.is_finite() && settings.youngs > 0.0) {
        return Err(invalid("robust stress evaluation requires finite positive Young's modulus"));
    }
    if !(settings.poisson.is_finite() && settings.poisson > -1.0 && settings.poisson < 0.5) {
        return Err(invalid("robust stress evaluation requires Poisson ratio in (-1, 0.5)"));
    }
    let material = IsotropicElastic::new(settings.youngs, settings.poisson, MATERIAL_STRAIN_LIMIT)
        .map_err(|error| invalid(format!("robust stress material card refused: {error}")))?;
    let (lambda, mu) = material.lame();
    let ratio = (lambda + 2.0 * mu) / mu;
    if !(lambda.is_finite() && mu.is_finite() && mu > 0.0
        && ratio.is_finite() && ratio <= MAX_PLANE_STRAIN_STIFFNESS_RATIO)
    {
        return Err(invalid(format!(
            "robust stress plane-strain stiffness ratio {ratio} exceeds the admitted regime"
        )));
    }
    Ok((material, lambda, mu))
}

fn validate_geometry(
    phi: &GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
) -> Result<Vec<EdgeBand>, CutFemError> {
    let expected = 1usize.checked_shl(settings.level)
        .ok_or_else(|| invalid("robust stress grid level exceeds platform address space"))?;
    if phi.n() != expected {
        return Err(invalid(format!(
            "robust stress SDF has {} cells per side but level {} requires {expected}",
            phi.n(), settings.level
        )));
    }
    if phi.nodes().iter().any(|value| !value.is_finite()) {
        return Err(invalid("robust stress evaluation requires finite level-set nodes"));
    }
    if load_cases.is_empty() || !load_cases.iter().any(|case| case.weight() > 0.0) {
        return Err(invalid("robust stress evaluation requires a load set with a positive objective weight"));
    }
    load_cases.iter().map(|case| {
        let [start, end] = case.interval();
        RobustLoadCase::new(case.edge(), start, end, case.traction(), case.weight())?;
        EdgeBand::new(case.edge(), start, end)
    }).collect()
}

fn strain_at(
    grid: &Quadtree,
    solution: &CutElasticitySolution,
    cell: (u32, u32, u32),
    point: [f64; 2],
) -> Option<[f64; 3]> {
    // A Q1 displacement gradient is discontinuous across cell boundaries.
    // Evaluate the quadrature OWNER's trace, including interface points exactly
    // on lattice lines; floor(point / h) can select an unrepresented neighbour.
    let (lo, hi) = grid.rect(cell);
    let corners = grid.corner_nodes(cell);
    let nodal = solution.nodal();
    let mut values = [[0.0; 2]; 4];
    for (index, corner) in corners.iter().enumerate() {
        values[index] = *nodal.get(corner)?;
    }
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = ((point[0] - lo[0]) / hx).clamp(0.0, 1.0);
    let eta = ((point[1] - lo[1]) / hy).clamp(0.0, 1.0);
    let gradients = [
        [-(1.0 - eta) / hx, -(1.0 - xi) / hy],
        [(1.0 - eta) / hx, -xi / hy],
        [eta / hx, xi / hy],
        [-eta / hx, (1.0 - xi) / hy],
    ];
    let mut grad_u = [[0.0; 2]; 2];
    for a in 0..4 {
        for component in 0..2 {
            grad_u[component][0] += gradients[a][0] * values[a][component];
            grad_u[component][1] += gradients[a][1] * values[a][component];
        }
    }
    Some([
        grad_u[0][0],
        grad_u[1][1],
        f64::midpoint(grad_u[0][1], grad_u[1][0]),
    ])
}

fn von_mises(lambda: f64, mu: f64, strain: [f64; 3]) -> f64 {
    let sxx = (lambda + 2.0 * mu) * strain[0] + lambda * strain[1];
    let syy = lambda * strain[0] + (lambda + 2.0 * mu) * strain[1];
    let szz = lambda * (strain[0] + strain[1]);
    let sxy = 2.0 * mu * strain[2];
    fs_material::tensor::von_mises(&[sxx, syy, szz, sxy, 0.0, 0.0])
}

fn full_cell_points(lo: [f64; 2], hi: [f64; 2]) -> [[f64; 2]; 5] {
    let mid = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
    let half = [0.5 * (hi[0] - lo[0]), 0.5 * (hi[1] - lo[1])];
    [
        [mid[0] - GAUSS * half[0], mid[1] - GAUSS * half[1]],
        [mid[0] + GAUSS * half[0], mid[1] - GAUSS * half[1]],
        [mid[0] + GAUSS * half[0], mid[1] + GAUSS * half[1]],
        [mid[0] - GAUSS * half[0], mid[1] + GAUSS * half[1]],
        mid,
    ]
}

fn sample_solution(
    grid: &Quadtree,
    phi: &GridSdf,
    solution: &CutElasticitySolution,
    lambda: f64,
    mu: f64,
) -> Result<(f64, [f64; 2], usize), CutFemError> {
    match sample_solution_controlled(grid, phi, solution, lambda, mu, |_| {
        ControlFlow::<Infallible>::Continue(())
    })? {
        ControlFlow::Continue(samples) => Ok(samples),
        ControlFlow::Break(never) => match never {},
    }
}

/// Sample a crate-owned matching displacement/geometry pair without another PDE
/// solve. Cancellation is polled before each quadrature cell. No partial maximum
/// escapes; callers must retain the previous complete state on Break or Err.
pub(crate) fn sample_solution_controlled<B>(
    grid: &Quadtree,
    phi: &GridSdf,
    solution: &CutElasticitySolution,
    lambda: f64,
    mu: f64,
    mut control: impl FnMut(usize) -> ControlFlow<B>,
) -> Result<ControlFlow<B, (f64, [f64; 2], usize)>, CutFemError> {
    let mut maximum = 0.0_f64;
    let mut location = [0.0, 0.0];
    let mut count = 0usize;
    let mut observe = |cell, point: [f64; 2]| -> Result<(), CutFemError> {
        let strain = strain_at(grid, solution, cell, point)
            .ok_or_else(|| invalid("material stress probe is missing its owning-cell displacement"))?;
        let stress = von_mises(lambda, mu, strain);
        if !(stress.is_finite() && stress >= 0.0) {
            return Err(invalid("robust sampled von Mises stress is invalid"));
        }
        count = count.checked_add(1)
            .ok_or_else(|| invalid("robust stress sample count overflowed"))?;
        if count == 1 || stress > maximum {
            maximum = stress;
            location = point;
        }
        Ok(())
    };
    for (ordinal, cell) in grid.leaves().enumerate() {
        if let ControlFlow::Break(reason) = control(ordinal) {
            return Ok(ControlFlow::Break(reason));
        }
        let (lo, hi) = grid.rect(cell);
        let enclosure = phi.enclose(lo, hi);
        if enclosure.lo() > 0.0 {
            continue;
        }
        if enclosure.hi() < 0.0 {
            for point in full_cell_points(lo, hi) {
                observe(cell, point)?;
            }
        } else {
            let rules = cut_cell_rules(phi, lo, hi, 2);
            // An interface-only exterior neighbour has no material-side Q1
            // displacement. Its trace belongs to the positive-volume cell on
            // the other side, not to a fabricated or silently skipped sample.
            if !rules.bulk.iter().any(|&(_, weight)| weight > 0.0) {
                continue;
            }
            for &(point, weight) in &rules.bulk {
                if weight > 0.0 { observe(cell, point)?; }
            }
            for &(point, weight, _) in &rules.iface {
                if weight > 0.0 { observe(cell, point)?; }
            }
        }
    }
    if count == 0 {
        return Err(invalid("robust stress evaluation found no material stress samples"));
    }
    Ok(ControlFlow::Continue((maximum, location, count)))
}

/// Re-solve every scenario on one exact geometry and bind robust compliance and
/// sampled stress evidence to the same snapshot.
///
/// # Errors
/// Refuses malformed declarations, canonical CutFEM failures, or invalid stress data.
pub fn evaluate_robust_sampled_stress(
    phi: &GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    aggregate: RobustAggregate,
) -> Result<RobustSampledStressEvaluation, CutFemError> {
    let supports = validate_geometry(phi, load_cases, settings)?;
    let (material, lambda, mu) = material(settings)?;
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
    let mut case_stresses = Vec::with_capacity(load_cases.len());
    let mut case_locations = Vec::with_capacity(load_cases.len());
    let mut case_counts = Vec::with_capacity(load_cases.len());
    let mut weighted_sum = 0.0;
    let mut worst_weighted = 0.0_f64;
    let mut worst_stress = f64::NEG_INFINITY;
    let mut worst_stress_case = 0usize;
    for (index, (case, support)) in load_cases.iter().zip(&supports).enumerate() {
        let value = case.traction();
        let traction = move |_: f64, _: f64| value;
        let solution = solver.solve_with_boundary_traction(
            &|_, _| [0.0, 0.0],
            &|_, _| [0.0, 0.0],
            BoundaryTraction::EdgeBand { support: *support, value: &traction },
        )?;
        let compliance = solution.compliance();
        if !(compliance.is_finite() && compliance >= 0.0) {
            return Err(invalid("robust stress load case produced invalid compliance"));
        }
        let (stress, location, count) = sample_solution(&grid, phi, &solution, lambda, mu)?;
        let weighted = case.weight() * compliance;
        weighted_sum += weighted;
        worst_weighted = worst_weighted.max(weighted);
        if stress > worst_stress {
            worst_stress = stress;
            worst_stress_case = index;
        }
        case_compliances.push(compliance);
        case_stresses.push(stress);
        case_locations.push(location);
        case_counts.push(count);
    }
    if !(weighted_sum.is_finite() && worst_weighted.is_finite() && worst_stress.is_finite()) {
        return Err(invalid("robust stress aggregate is non-finite"));
    }
    let objective = match aggregate {
        RobustAggregate::WeightedSum => weighted_sum,
        RobustAggregate::WorstWeightedCase => worst_weighted,
    };
    let volume = crate::optimize::material_volume(&grid, phi);
    if !(volume.is_finite() && volume > 0.0) {
        return Err(invalid("robust stress design has invalid material area"));
    }
    Ok(RobustSampledStressEvaluation {
        case_compliances,
        case_sampled_max_von_mises: case_stresses,
        case_max_locations: case_locations,
        case_sample_counts: case_counts,
        weighted_sum_compliance: weighted_sum,
        worst_weighted_compliance: worst_weighted,
        objective,
        worst_sampled_von_mises: worst_stress,
        worst_stress_case,
        volume,
        snapshot: snapshot(phi),
    })
}

fn validate_guarded(
    settings: OptimizeSettings,
    guarded: GuardedSettings,
    limit: SampledStressLimit,
) -> Result<(), CutFemError> {
    SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
    if settings.iterations == 0 {
        return Err(invalid("robust stress optimization requires at least one iteration"));
    }
    if !(settings.move_cells.is_finite() && settings.move_cells > 0.0) {
        return Err(invalid("robust stress optimization requires finite positive move_cells"));
    }
    if !(settings.volfrac.is_finite() && settings.volfrac > 0.0 && settings.volfrac <= 1.0) {
        return Err(invalid("robust stress optimization requires volfrac in (0, 1]"));
    }
    if guarded.max_candidates == 0 || guarded.max_candidates > 64 {
        return Err(invalid("robust stress max_candidates must lie in [1, 64]"));
    }
    if !(guarded.contraction.is_finite() && guarded.contraction > 0.0 && guarded.contraction < 1.0) {
        return Err(invalid("robust stress contraction must lie strictly between zero and one"));
    }
    if !(guarded.volume_tolerance.is_finite() && guarded.volume_tolerance >= 0.0
        && guarded.volume_tolerance <= 1.0)
    {
        return Err(invalid("robust stress volume_tolerance must lie in [0, 1]"));
    }
    if !(guarded.min_relative_improvement.is_finite()
        && guarded.min_relative_improvement >= 0.0
        && guarded.min_relative_improvement < 1.0)
    {
        return Err(invalid("robust stress min_relative_improvement must lie in [0, 1)"));
    }
    Ok(())
}

/// Bounded transactional simultaneous multi-load design with worst-scenario
/// sampled-stress admission and independent final replay.
///
/// # Errors
/// Refuses malformed controls or baseline evaluation. Candidate refusals are retained.
pub fn optimize_compliance_multi_load_stress_guarded(
    phi: &mut GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    guarded: GuardedSettings,
    aggregate: RobustAggregate,
    limit: SampledStressLimit,
) -> Result<RobustStressReport, CutFemError> {
    validate_guarded(settings, guarded, limit)?;
    let baseline = evaluate_robust_sampled_stress(phi, load_cases, settings, aggregate)?;
    let volume_limit = settings.volfrac + guarded.volume_tolerance;
    let stress_limit = limit.admitted_max();
    let baseline_constraints = baseline.volume <= volume_limit
        && baseline.worst_sampled_von_mises <= stress_limit;
    let improvement_limit = baseline.objective * (1.0 - guarded.min_relative_improvement);
    let origin = phi.clone();
    let mut candidates = Vec::with_capacity(guarded.max_candidates);
    let mut best: Option<(GridSdf, RobustSampledStressEvaluation, OptimizeReport, f64)> = None;
    let mut move_cells = settings.move_cells;
    let mut any_final = false;
    let mut any_constraint_feasible = false;

    for index in 0..guarded.max_candidates {
        let mut trial_settings = settings;
        trial_settings.move_cells = move_cells;
        let mut trial = origin.clone();
        match optimize_compliance_multi_load(&mut trial, load_cases, trial_settings, aggregate) {
            Ok(descent) => match evaluate_robust_sampled_stress(&trial, load_cases, settings, aggregate) {
                Ok(evaluation) => {
                    any_final = true;
                    let volume_feasible = evaluation.volume <= volume_limit;
                    let stress_feasible = evaluation.worst_sampled_von_mises <= stress_limit;
                    let constraint_feasible = volume_feasible && stress_feasible;
                    any_constraint_feasible |= constraint_feasible;
                    let improvement_gate = !baseline_constraints
                        || evaluation.objective <= improvement_limit;
                    if constraint_feasible && improvement_gate {
                        let replace = best.as_ref().is_none_or(|(_, current, _, _)| {
                            evaluation.objective < current.objective
                        });
                        if replace {
                            best = Some((trial, evaluation.clone(), descent.trajectory.clone(), move_cells));
                        }
                    }
                    candidates.push(RobustStressCandidate {
                        index,
                        move_cells,
                        evaluation: Some(evaluation),
                        volume_feasible,
                        stress_feasible,
                        improvement_gate,
                        refusal: None,
                    });
                }
                Err(error) => candidates.push(RobustStressCandidate {
                    index,
                    move_cells,
                    evaluation: None,
                    volume_feasible: false,
                    stress_feasible: false,
                    improvement_gate: false,
                    refusal: Some(format!("{error:?}")),
                }),
            },
            Err(error) => candidates.push(RobustStressCandidate {
                index,
                move_cells,
                evaluation: None,
                volume_feasible: false,
                stress_feasible: false,
                improvement_gate: false,
                refusal: Some(format!("{error:?}")),
            }),
        }
        move_cells *= guarded.contraction;
    }

    if let Some((accepted_phi, accepted, trajectory, accepted_move_cells)) = best {
        *phi = accepted_phi;
        return Ok(RobustStressReport {
            baseline,
            limit,
            candidates,
            stop: RobustStressStop::Accepted,
            accepted: Some(accepted),
            trajectory: Some(trajectory),
            accepted_move_cells: Some(accepted_move_cells),
        });
    }
    let stop = if !any_final {
        RobustStressStop::AllCandidatesRefused
    } else if !any_constraint_feasible {
        RobustStressStop::NoConstraintFeasibleCandidate
    } else {
        RobustStressStop::NoImprovingCandidate
    };
    Ok(RobustStressReport {
        baseline,
        limit,
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
    use fs_cutfem::DesignBoxEdge;

    fn beam(level: u32) -> GridSdf {
        GridSdf::from_fn(1usize << level, &|_, y| (y - 0.5).abs() - 0.35)
    }

    fn settings(level: u32) -> OptimizeSettings {
        OptimizeSettings { level, iterations: 1, nucleation_period: 0, ..OptimizeSettings::default() }
    }

    #[test]
    fn aligned_interface_samples_cover_every_positive_volume_owner() {
        let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.25);
        let grid = Quadtree::uniform(3);
        let cases = [RobustLoadCase::new(
            DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0,
        ).unwrap()];
        let evaluation = evaluate_robust_sampled_stress(
            &phi, &cases, settings(3), RobustAggregate::WeightedSum,
        ).expect("material-side traces on an exactly lattice-aligned interface");
        let mut expected = 0;
        for cell in grid.leaves() {
            let (lo, hi) = grid.rect(cell);
            let enclosure = phi.enclose(lo, hi);
            if enclosure.lo() > 0.0 { continue; }
            if enclosure.hi() < 0.0 {
                expected += 5;
            } else {
                let rules = cut_cell_rules(&phi, lo, hi, 2);
                let bulk = rules.bulk.iter().filter(|(_, weight)| *weight > 0.0).count();
                if bulk > 0 {
                    expected += bulk + rules.iface.iter().filter(|(_, weight, _)| *weight > 0.0).count();
                }
            }
        }
        assert!(expected > 0);
        assert_eq!(evaluation.case_sample_counts, [expected]);
    }

    #[test]
    fn opposite_scenarios_have_identical_sampled_stress() {
        let phi = beam(3);
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.5).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 1.0], 0.5).unwrap(),
        ];
        let evaluation = evaluate_robust_sampled_stress(
            &phi, &cases, settings(3), RobustAggregate::WeightedSum
        ).expect("opposite independent scenarios");
        let scale = evaluation.case_sampled_max_von_mises[0].max(1.0);
        assert!((evaluation.case_sampled_max_von_mises[0]
            - evaluation.case_sampled_max_von_mises[1]).abs() < 1e-10 * scale);
        assert_eq!(evaluation.case_sample_counts[0], evaluation.case_sample_counts[1]);
    }

    #[test]
    fn robust_stress_refusal_is_transactional() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.5).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.25, 0.4, [0.5, 0.0], 0.5).unwrap(),
        ];
        let report = optimize_compliance_multi_load_stress_guarded(
            &mut phi,
            &cases,
            OptimizeSettings {
                level: 3, iterations: 1, volfrac: 1.0, move_cells: 0.02,
                nucleation_period: 0, ..OptimizeSettings::default()
            },
            GuardedSettings {
                max_candidates: 2, contraction: 0.5, volume_tolerance: 0.0,
                min_relative_improvement: 0.0,
            },
            RobustAggregate::WeightedSum,
            SampledStressLimit::new(1e-30, 0.0).unwrap(),
        ).expect("bounded robust stress refusal is data");
        assert_ne!(report.stop, RobustStressStop::Accepted);
        assert!(report.accepted.is_none());
        assert_eq!(phi.nodes(), before);
    }
}
