//! Simultaneous multi-load level-set compliance descent.
//!
//! Each load case gets its own canonical CutFEM equilibrium solve on the same
//! current geometry. Shape energies and topological derivatives are aggregated
//! only AFTER those independent solves. Weighted-sum descent combines every
//! case; worst-weighted descent uses the current active scenario (lowest input
//! index at an exact tie), which is one valid active-branch subgradient rather
//! than a smooth max derivative.

use crate::fim::redistance;
use crate::gridsdf::GridSdf;
use crate::guarded::GuardedSettings;
use crate::optimize::{OptimizeReport, OptimizeSettings, material_volume};
use crate::robust::{
    RobustAggregate, RobustCandidate, RobustLoadCase, RobustOptimizeReport, RobustStop,
    evaluate_robust_design,
};
use crate::topder::{nucleate, topological_derivative};
use crate::veloext::extend_velocity;
use crate::weno::{Velocity, advect, build_band};
use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutElasticitySolution, CutFemError,
    CutStabilizationScaling, DesignBoxEdge, EdgeBand, MAX_PLANE_STRAIN_STIFFNESS_RATIO,
    Quadtree,
};
use fs_material::IsotropicElastic;
use std::fmt::Write as _;

const MATERIAL_STRAIN_LIMIT: f64 = 1.0;
const SOLVER_TOL: f64 = 1e-12;
const SOLVER_MAX_ITERS: usize = 60_000;

/// Multi-load evolution evidence. The ordinary trajectory's `compliance` field
/// stores the selected aggregate objective at each iteration.
#[derive(Debug, Clone)]
pub struct RobustDescentReport {
    /// Existing geometric/evolution evidence with aggregate objective rows.
    pub trajectory: OptimizeReport,
    /// Unweighted per-case compliances for every iteration, in input order.
    pub case_compliances: Vec<Vec<f64>>,
    /// Active case for [`RobustAggregate::WorstWeightedCase`], otherwise `None`.
    pub active_case: Vec<Option<usize>>,
    /// Aggregate used by the descent.
    pub aggregate: RobustAggregate,
}

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn material(settings: OptimizeSettings) -> Result<(IsotropicElastic, f64, f64), CutFemError> {
    if !(settings.youngs.is_finite() && settings.youngs > 0.0) {
        return Err(invalid("multi-load optimizer requires finite positive Young's modulus"));
    }
    if !(settings.poisson.is_finite() && settings.poisson > -1.0 && settings.poisson < 0.5) {
        return Err(invalid("multi-load optimizer requires Poisson ratio in (-1, 0.5)"));
    }
    let material = IsotropicElastic::new(settings.youngs, settings.poisson, MATERIAL_STRAIN_LIMIT)
        .map_err(|error| invalid(format!("multi-load material card refused: {error}")))?;
    let (lambda, mu) = material.lame();
    let ratio = (lambda + 2.0 * mu) / mu;
    if !(lambda.is_finite() && mu.is_finite() && mu > 0.0
        && ratio.is_finite() && ratio <= MAX_PLANE_STRAIN_STIFFNESS_RATIO)
    {
        return Err(invalid(format!(
            "multi-load plane-strain stiffness ratio {ratio} exceeds the admitted regime"
        )));
    }
    Ok((material, lambda, mu))
}

fn validate(
    phi: &GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
) -> Result<Vec<EdgeBand>, CutFemError> {
    if load_cases.is_empty() || !load_cases.iter().any(|case| case.weight() > 0.0) {
        return Err(invalid("multi-load descent requires at least one case and one positive weight"));
    }
    if !(settings.volfrac.is_finite() && settings.volfrac > 0.0 && settings.volfrac <= 1.0) {
        return Err(invalid("multi-load descent requires volfrac in (0, 1]"));
    }
    if !(settings.move_cells.is_finite() && settings.move_cells > 0.0) {
        return Err(invalid("multi-load descent requires finite positive move_cells"));
    }
    let expected = 1usize.checked_shl(settings.level)
        .ok_or_else(|| invalid("multi-load grid level exceeds platform address space"))?;
    if phi.n() != expected {
        return Err(invalid(format!(
            "multi-load SDF has {} cells per side but level {} requires {expected}",
            phi.n(), settings.level
        )));
    }
    if phi.nodes().iter().any(|value| !value.is_finite()) {
        return Err(invalid("multi-load descent requires finite level-set nodes"));
    }
    load_cases.iter().map(|case| {
        let [start, end] = case.interval();
        EdgeBand::new(case.edge(), start, end)
    }).collect()
}

fn mass_stiffness(n: usize) -> (fs_sparse::Csr, fs_sparse::Csr) {
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / n as f64;
    let stride = n + 1;
    let nn = stride * stride;
    let mut mc = fs_sparse::Coo::new(nn, nn);
    let mut kc = fs_sparse::Coo::new(nn, nn);
    let ke = [
        [2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0, -1.0 / 3.0],
        [-1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0, -1.0 / 6.0],
        [-1.0 / 6.0, -1.0 / 3.0, -1.0 / 6.0, 2.0 / 3.0],
    ];
    for cj in 0..n {
        for ci in 0..n {
            let ids = [
                ci + cj * stride,
                ci + 1 + cj * stride,
                ci + 1 + (cj + 1) * stride,
                ci + (cj + 1) * stride,
            ];
            for (a, &ia) in ids.iter().enumerate() {
                mc.push(ia, ia, 0.25 * h * h);
                for (b, &ib) in ids.iter().enumerate() {
                    kc.push(ia, ib, ke[a][b]);
                }
            }
        }
    }
    (mc.assemble(), kc.assemble())
}

fn retain_load_pads(phi: &mut GridSdf, supports: &[EdgeBand]) -> usize {
    let h = phi.h();
    let mut changed = 0usize;
    for support in supports {
        let (x_min, x_max, y_min, y_max) = match support.edge() {
            DesignBoxEdge::Right => (1.0 - 2.0*h, 1.0 + h, support.start() - h, support.end() + h),
            DesignBoxEdge::Top => (support.start() - h, support.end() + h, 1.0 - 2.0*h, 1.0 + h),
            DesignBoxEdge::Bottom => (support.start() - h, support.end() + h, -h, 2.0*h),
            DesignBoxEdge::Left => continue,
        };
        for j in 0..=phi.n() {
            for i in 0..=phi.n() {
                let [x, y] = phi.pos(i, j);
                let pad = (x_min - x).max(x - x_max).max(y_min - y).max(y - y_max);
                if pad <= 0.0 && pad < phi.node(i, j) {
                    *phi.node_mut(i, j) = pad;
                    changed += 1;
                }
            }
        }
    }
    changed
}

fn fnv(phi: &GridSdf) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for value in phi.nodes() {
        for byte in value.to_bits().to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

fn strain_at(
    grid: &Quadtree,
    sol: &CutElasticitySolution,
    p: [f64; 2],
) -> ([f64; 3], bool) {
    let level = grid.max_level();
    let nf = f64::from(1u32 << level);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let ci = ((p[0] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let cj = ((p[1] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    let cell = (level, ci, cj);
    let (lo, hi) = grid.rect(cell);
    let corners = grid.corner_nodes(cell);
    let nodal = sol.nodal();
    let mut vals = [[0.0f64; 2]; 4];
    for (a, c) in corners.iter().enumerate() {
        match nodal.get(c) {
            Some(u) => vals[a] = *u,
            None => return ([0.0; 3], false),
        }
    }
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = ((p[0] - lo[0]) / hx).clamp(0.0, 1.0);
    let et = ((p[1] - lo[1]) / hy).clamp(0.0, 1.0);
    let g = [
        [-(1.0 - et) / hx, -(1.0 - xi) / hy],
        [(1.0 - et) / hx, -xi / hy],
        [et / hx, xi / hy],
        [-et / hx, (1.0 - xi) / hy],
    ];
    let mut gu = [[0.0f64; 2]; 2];
    for a in 0..4 {
        for c in 0..2 {
            gu[c][0] += g[a][0] * vals[a][c];
            gu[c][1] += g[a][1] * vals[a][c];
        }
    }
    ([gu[0][0], gu[1][1], f64::midpoint(gu[0][1], gu[1][0])], true)
}

fn stress_energy(lambda: f64, mu: f64, eps: [f64; 3]) -> ([f64; 3], f64) {
    let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
    let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
    let sxy = 2.0 * mu * eps[2];
    let energy = 0.5 * (sxx * eps[0] + syy * eps[1] + 2.0 * sxy * eps[2]);
    ([sxx, syy, sxy], energy)
}

fn active_case(compliances: &[f64], loads: &[RobustLoadCase]) -> usize {
    let mut active = 0usize;
    let mut best = f64::NEG_INFINITY;
    for (index, (&compliance, load)) in compliances.iter().zip(loads).enumerate() {
        let value = load.weight() * compliance;
        if value > best {
            best = value;
            active = index;
        }
    }
    active
}

/// Run simultaneous multi-load level-set compliance descent.
///
/// Weighted-sum mode accumulates every case's self-adjoint strain-energy shape
/// field and topological derivative. Worst-weighted mode uses the current active
/// weighted compliance case; at an exact tie the lowest case index supplies the
/// active-branch subgradient. Every traction support is retained as non-design
/// material after each geometry update.
///
/// # Errors
/// Propagates typed material, load-support, and canonical CutFEM solve refusals.
pub fn optimize_compliance_multi_load(
    phi: &mut GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    aggregate: RobustAggregate,
) -> Result<RobustDescentReport, CutFemError> {
    let supports = validate(phi, load_cases, settings)?;
    let (material, lambda, mu) = material(settings)?;
    let n = phi.n();
    let grid = Quadtree::uniform(settings.level);
    let h = phi.h();
    let stride = n + 1;
    let (mass, stiffness) = mass_stiffness(n);
    let clamp = |x: f64, _y: f64| x < 1e-9;
    let mut ell = settings.ell0;
    let mut trajectory = OptimizeReport::default();
    let mut case_history = Vec::with_capacity(settings.iterations);
    let mut active_history = Vec::with_capacity(settings.iterations);

    for iter in 0..settings.iterations {
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
        let mut solutions = Vec::with_capacity(load_cases.len());
        let mut compliances = Vec::with_capacity(load_cases.len());
        for (case, support) in load_cases.iter().zip(&supports) {
            let value = case.traction();
            let traction = move |_: f64, _: f64| value;
            let solution = solver.solve_with_boundary_traction(
                &|_, _| [0.0, 0.0],
                &|_, _| [0.0, 0.0],
                BoundaryTraction::EdgeBand { support: *support, value: &traction },
            )?;
            compliances.push(solution.compliance());
            solutions.push(solution);
        }
        let active = match aggregate {
            RobustAggregate::WeightedSum => None,
            RobustAggregate::WorstWeightedCase => Some(active_case(&compliances, load_cases)),
        };
        let objective = match aggregate {
            RobustAggregate::WeightedSum => compliances.iter().zip(load_cases)
                .map(|(&c, load)| load.weight() * c).sum::<f64>(),
            RobustAggregate::WorstWeightedCase => {
                let index = active.expect("worst-case descent selects an active case");
                load_cases[index].weight() * compliances[index]
            }
        };
        if !objective.is_finite() {
            return Err(invalid("multi-load aggregate compliance became nonfinite"));
        }

        let mut energy = vec![0.0f64; stride * stride];
        let mut seeded = vec![false; stride * stride];
        for j in 0..=n {
            for i in 0..=n {
                let k = i + j * stride;
                let p = phi.pos(i, j);
                let g = phi.gradient_at(p);
                let gn = g[0].hypot(g[1]).max(1e-12);
                let q = [
                    (p[0] - 0.75 * h * g[0] / gn).clamp(0.0, 1.0),
                    (p[1] - 0.75 * h * g[1] / gn).clamp(0.0, 1.0),
                ];
                if phi.value_at(q) > 0.0 {
                    continue;
                }
                let mut value = 0.0;
                match active {
                    Some(index) => {
                        let (eps, ok) = strain_at(&grid, &solutions[index], q);
                        if ok {
                            value = load_cases[index].weight() * stress_energy(lambda, mu, eps).1;
                        }
                    }
                    None => {
                        for (solution, load) in solutions.iter().zip(load_cases) {
                            let (eps, ok) = strain_at(&grid, solution, q);
                            if ok {
                                value += load.weight() * stress_energy(lambda, mu, eps).1;
                            }
                        }
                    }
                }
                energy[k] = value;
                seeded[k] = phi.node(i, j).abs() <= 2.0 * h;
            }
        }
        extend_velocity(phi, &mut energy, &seeded);
        let (smooth, _iters) = fs_adjoint::sobolev::sobolev_smooth(
            &mass,
            &stiffness,
            settings.sobolev_alpha * h * h,
            &energy,
            1e-10,
        );
        let vn: Vec<f64> = smooth.iter().map(|w| w - ell).collect();
        let band = build_band(phi, settings.band_cells);
        let vmax = vn.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-12);
        advect(
            phi,
            &band,
            &Velocity::Normal(&vn),
            settings.move_cells * h / vmax,
            0.45,
        );
        let mut load_pad_nodes = retain_load_pads(phi, &supports);
        let audit = redistance(phi, settings.band_cells);
        load_pad_nodes += retain_load_pads(phi, &supports);
        let volume = material_volume(&grid, phi);
        #[allow(clippy::cast_precision_loss)]
        let w_mean = smooth.iter().sum::<f64>() / smooth.len() as f64;
        ell = (ell + settings.mu_al * w_mean.abs().max(1e-30)
            * (volume - settings.volfrac) / settings.volfrac).max(0.0);

        if settings.nucleation_period > 0 && iter > 0 && iter % settings.nucleation_period == 0 {
            let mut dt_field = vec![f64::INFINITY; stride * stride];
            for j in 0..=n {
                for i in 0..=n {
                    let k = i + j * stride;
                    let p = phi.pos(i, j);
                    if phi.value_at(p) > -2.0 * h {
                        continue;
                    }
                    let mut value = 0.0;
                    let mut any = false;
                    match active {
                        Some(index) => {
                            let (eps, ok) = strain_at(&grid, &solutions[index], p);
                            if ok {
                                let (stress, _) = stress_energy(lambda, mu, eps);
                                value = load_cases[index].weight()
                                    * topological_derivative(lambda, mu, stress, eps);
                                any = true;
                            }
                        }
                        None => {
                            for (solution, load) in solutions.iter().zip(load_cases) {
                                let (eps, ok) = strain_at(&grid, solution, p);
                                if ok {
                                    let (stress, _) = stress_energy(lambda, mu, eps);
                                    value += load.weight()
                                        * topological_derivative(lambda, mu, stress, eps);
                                    any = true;
                                }
                            }
                        }
                    }
                    if any { dt_field[k] = value; }
                }
            }
            let events = nucleate(
                phi,
                &dt_field,
                ell,
                settings.hole_radius_cells * h,
                6.0 * settings.hole_radius_cells * h,
                2,
            );
            if !events.is_empty() {
                load_pad_nodes += retain_load_pads(phi, &supports);
                let _ = redistance(phi, settings.band_cells);
                load_pad_nodes += retain_load_pads(phi, &supports);
            }
            trajectory.events.extend(events);
        }

        let snap = fnv(phi);
        let cases = compliances.iter().map(|value| format!("{value:.6e}"))
            .collect::<Vec<_>>().join(",");
        let active_json = active.map_or_else(|| "null".to_string(), |index| index.to_string());
        let mut row = String::new();
        let _ = write!(row,
            "{{\"iter\":{iter},\"compliance\":{objective:.6e},\"case_compliances\":[{cases}],\"active_case\":{active_json},\"volume\":{volume:.4},\"ell\":{ell:.4e},\"drift_h\":{:.2e},\"load_pad_nodes\":{load_pad_nodes},\"snapshot\":\"{snap:#018x}\"}}",
            audit.interface_drift_h);
        trajectory.rows.push(row);
        trajectory.compliance.push(objective);
        trajectory.volume.push(volume);
        trajectory.ell.push(ell);
        trajectory.audits.push(audit);
        trajectory.snapshots.push(snap);
        trajectory.load_pad_nodes.push(load_pad_nodes);
        case_history.push(compliances);
        active_history.push(active);
    }

    Ok(RobustDescentReport {
        trajectory,
        case_compliances: case_history,
        active_case: active_history,
        aggregate,
    })
}

fn validate_guarded(settings: OptimizeSettings, guarded: GuardedSettings) -> Result<(), CutFemError> {
    if settings.iterations == 0 {
        return Err(invalid("multi-load guarded optimization requires at least one evolution iteration"));
    }
    if guarded.max_candidates == 0 || guarded.max_candidates > 64 {
        return Err(invalid("multi-load guarded max_candidates must lie in [1, 64]"));
    }
    if !(guarded.contraction.is_finite() && guarded.contraction > 0.0 && guarded.contraction < 1.0) {
        return Err(invalid("multi-load guarded contraction must lie strictly between zero and one"));
    }
    if !(guarded.volume_tolerance.is_finite()
        && guarded.volume_tolerance >= 0.0
        && guarded.volume_tolerance <= 1.0)
    {
        return Err(invalid("multi-load guarded volume_tolerance must lie in [0, 1]"));
    }
    if !(guarded.min_relative_improvement.is_finite()
        && guarded.min_relative_improvement >= 0.0
        && guarded.min_relative_improvement < 1.0)
    {
        return Err(invalid("multi-load guarded min_relative_improvement must lie in [0, 1)"));
    }
    Ok(())
}

/// Bounded transactional candidate search whose CANDIDATE GENERATION already
/// follows the simultaneous multi-load shape field.
///
/// This differs from [`crate::robust::optimize_compliance_robust_guarded`],
/// which deliberately retains the older single nominal driver and uses robust
/// scenarios only for publication. Here every candidate trajectory is driven by
/// `optimize_compliance_multi_load`, then independently replayed once more under
/// all scenarios before publication. The second replay is the authoritative
/// objective/area binding for the returned geometry.
///
/// # Errors
/// Refuses malformed settings or baseline robust evaluation. Individual
/// candidate evolution/replay refusals are retained and never mutate `phi`.
pub fn optimize_compliance_multi_load_guarded(
    phi: &mut GridSdf,
    load_cases: &[RobustLoadCase],
    settings: OptimizeSettings,
    guarded: GuardedSettings,
    aggregate: RobustAggregate,
) -> Result<RobustOptimizeReport, CutFemError> {
    validate_guarded(settings, guarded)?;
    let baseline = evaluate_robust_design(phi, load_cases, settings, aggregate)?;
    let limit = settings.volfrac + guarded.volume_tolerance;
    let baseline_feasible = baseline.volume <= limit;
    let improvement_limit = baseline.objective * (1.0 - guarded.min_relative_improvement);
    let origin = phi.clone();
    let mut candidates = Vec::with_capacity(guarded.max_candidates);
    let mut best: Option<(GridSdf, crate::robust::RobustEvaluation, OptimizeReport, f64)> = None;
    let mut move_cells = settings.move_cells;
    let mut any_final = false;
    let mut any_feasible = false;

    for index in 0..guarded.max_candidates {
        let mut trial_settings = settings;
        trial_settings.move_cells = move_cells;
        let mut trial = origin.clone();
        match optimize_compliance_multi_load(&mut trial, load_cases, trial_settings, aggregate) {
            Ok(descent) => match evaluate_robust_design(&trial, load_cases, settings, aggregate) {
                Ok(evaluation) => {
                    any_final = true;
                    let volume_feasible = evaluation.volume <= limit;
                    any_feasible |= volume_feasible;
                    let improvement_gate = !baseline_feasible || evaluation.objective <= improvement_limit;
                    if volume_feasible && improvement_gate {
                        let replace = best.as_ref().is_none_or(|(_, current, _, _)| {
                            evaluation.objective < current.objective
                        });
                        if replace {
                            best = Some((trial, evaluation.clone(), descent.trajectory.clone(), move_cells));
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

    fn settings() -> OptimizeSettings {
        OptimizeSettings {
            level: 3,
            iterations: 1,
            nucleation_period: 0,
            move_cells: 0.1,
            ..OptimizeSettings::default()
        }
    }

    #[test]
    fn one_case_matches_the_existing_nominal_descent() {
        let mut nominal = beam(3);
        let mut multi = nominal.clone();
        let cfg = settings();
        crate::optimize::optimize_compliance(
            &mut nominal, crate::optimize::Cantilever { load: 1.0, band: 0.125 }, cfg
        ).unwrap();
        let case = RobustLoadCase::new(
            DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0
        ).unwrap();
        optimize_compliance_multi_load(
            &mut multi, &[case], cfg, RobustAggregate::WeightedSum
        ).unwrap();
        assert_eq!(nominal.nodes(), multi.nodes());
    }

    #[test]
    fn weighted_opposite_cases_have_nonzero_objective_and_one_history_row() {
        let mut phi = beam(3);
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.5).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 1.0], 0.5).unwrap(),
        ];
        let report = optimize_compliance_multi_load(
            &mut phi, &cases, settings(), RobustAggregate::WeightedSum
        ).unwrap();
        assert_eq!(report.trajectory.compliance.len(), 1);
        assert_eq!(report.case_compliances[0].len(), 2);
        assert!(report.trajectory.compliance[0] > 0.0);
        assert!((report.case_compliances[0][0] - report.case_compliances[0][1]).abs()
            < 1e-9 * report.case_compliances[0][0].max(1.0));
        assert_eq!(report.active_case, vec![None]);
    }

    #[test]
    fn worst_weighted_case_records_the_active_scenario() {
        let mut phi = beam(3);
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.4, 0.6, [0.0, -1.0], 0.25).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.20, 0.35, [1.0, 0.0], 1.0).unwrap(),
        ];
        let report = optimize_compliance_multi_load(
            &mut phi, &cases, settings(), RobustAggregate::WorstWeightedCase
        ).unwrap();
        let active = report.active_case[0].unwrap();
        assert!(active < cases.len());
        assert_eq!(report.trajectory.compliance[0], cases[active].weight()*report.case_compliances[0][active]);
        for (index, case) in cases.iter().enumerate() {
            assert!(report.trajectory.compliance[0] >= case.weight()*report.case_compliances[0][index]);
        }
    }

    #[test]
    fn top_and_bottom_load_pads_are_retained_as_strict_material() {
        let mut phi = beam(3);
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Top, 0.4, 0.6, [1.0, 0.0], 1.0).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Bottom, 0.4, 0.6, [-1.0, 0.0], 1.0).unwrap(),
        ];
        let supports = validate(&phi, &cases, settings()).unwrap();
        assert!(retain_load_pads(&mut phi, &supports) > 0);
        assert!(phi.value_at([0.5, 1.0]) < 0.0);
        assert!(phi.value_at([0.5, 0.0]) < 0.0);
    }

    #[test]
    fn multi_load_guarded_impossible_improvement_preserves_geometry() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.5).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 1.0], 0.5).unwrap(),
        ];
        let mut cfg = settings();
        cfg.volfrac = 1.0;
        cfg.move_cells = 1e-9;
        let guarded = GuardedSettings {
            max_candidates: 1, contraction: 0.5, volume_tolerance: 0.0,
            min_relative_improvement: 0.5,
        };
        let report = optimize_compliance_multi_load_guarded(
            &mut phi, &cases, cfg, guarded, RobustAggregate::WeightedSum
        ).unwrap();
        assert_ne!(report.stop, RobustStop::Accepted);
        assert_eq!(phi.nodes(), before.as_slice());
    }
}
