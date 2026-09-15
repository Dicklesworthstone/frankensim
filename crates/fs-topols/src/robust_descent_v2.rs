//! Evaluated-state simultaneous multi-load level-set compliance descent.
//!
//! Every load case is solved independently on one exact geometry. Shape and
//! hole sensitivities are formed only from those matching solutions, geometry
//! evolution happens on a trial copy, and the complete post-evolution geometry
//! is re-solved under every load before any trajectory row or caller mutation is
//! published. Weighted-sum and active worst-weighted objectives therefore bind
//! objective, per-case compliance, area, and snapshot to the same design.

use crate::fim::{RedistanceAudit, redistance};
use crate::gridsdf::GridSdf;
use crate::guarded::GuardedSettings;
use crate::optimize::{OptimizeReport, OptimizeSettings, material_volume};
use crate::robust::{
    RobustAggregate, RobustCandidate, RobustEvaluation, RobustLoadCase, RobustOptimizeReport,
    RobustStop, evaluate_robust_design,
};
use crate::topder::{NucleationEvent, nucleate, topological_derivative};
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

/// Simultaneous multi-load evolution evidence. Every row and case-compliance
/// vector refers to the successfully solved post-evolution geometry for that
/// iteration.
#[derive(Debug, Clone)]
pub struct RobustDescentReport {
    /// Geometric/evolution evidence with the selected aggregate in `compliance`.
    pub trajectory: OptimizeReport,
    /// Unweighted per-case compliances for each returned iteration, input order.
    pub case_compliances: Vec<Vec<f64>>,
    /// Active case for [`RobustAggregate::WorstWeightedCase`], else `None`.
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
    if !(settings.band_cells.is_finite() && settings.band_cells > 0.0
        && settings.move_cells.is_finite() && settings.move_cells >= 0.0
        && settings.move_cells <= 0.5 * settings.band_cells)
    {
        return Err(invalid(
            "multi-load move must be finite, nonnegative, and at most half the positive band width",
        ));
    }
    if !(settings.ell0.is_finite() && settings.ell0 >= 0.0
        && settings.mu_al.is_finite() && settings.mu_al > 0.0
        && settings.sobolev_alpha.is_finite() && settings.sobolev_alpha >= 0.0)
    {
        return Err(invalid("multi-load multiplier, penalty, or smoothing controls are invalid"));
    }
    if settings.nucleation_period > 0
        && !(settings.hole_radius_cells.is_finite() && settings.hole_radius_cells > 0.0)
    {
        return Err(invalid("enabled multi-load nucleation requires a finite positive hole radius"));
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

fn strain_at(grid: &Quadtree, sol: &CutElasticitySolution, p: [f64; 2]) -> ([f64; 3], bool) {
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

mod engine;
pub use engine::optimize_compliance_multi_load;
mod guarded_impl;
pub use guarded_impl::optimize_compliance_multi_load_guarded;

#[cfg(test)]
mod tests {
    use super::*;

    fn beam(level: u32) -> GridSdf {
        GridSdf::from_fn(1usize << level, &|_, y| (y - 0.5).abs() - 0.35)
    }

    fn loads() -> [RobustLoadCase; 2] {
        [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.6).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.25, 0.4, [0.5, 0.0], 0.4).unwrap(),
        ]
    }

    #[test]
    fn trajectory_metrics_bind_to_returned_geometry() {
        let settings = OptimizeSettings {
            level: 3,
            iterations: 1,
            move_cells: 0.05,
            nucleation_period: 0,
            ..OptimizeSettings::default()
        };
        let mut phi = beam(settings.level);
        let report = optimize_compliance_multi_load(
            &mut phi,
            &loads(),
            settings,
            RobustAggregate::WeightedSum,
        ).expect("evaluated multi-load step");
        let replay = evaluate_robust_design(
            &phi,
            &loads(),
            settings,
            RobustAggregate::WeightedSum,
        ).expect("independent final replay");
        assert_eq!(report.trajectory.snapshots.last(), Some(&replay.snapshot));
        assert_eq!(report.trajectory.compliance.last().unwrap().to_bits(), replay.objective.to_bits());
        assert_eq!(report.trajectory.volume.last().unwrap().to_bits(), replay.volume.to_bits());
        assert_eq!(report.case_compliances.last().unwrap(), &replay.case_compliances);
    }

    #[test]
    fn malformed_controls_do_not_mutate_geometry() {
        let original = beam(3);
        let mut phi = original.clone();
        let settings = OptimizeSettings {
            level: 3,
            iterations: 1,
            move_cells: f64::INFINITY,
            ..OptimizeSettings::default()
        };
        assert!(optimize_compliance_multi_load(
            &mut phi,
            &loads(),
            settings,
            RobustAggregate::WeightedSum,
        ).is_err());
        assert_eq!(phi.nodes(), original.nodes());
    }
}
