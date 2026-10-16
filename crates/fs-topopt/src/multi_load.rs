//! Bounded, independently loaded SIMP topology studies.
//!
//! Each load is an independent equilibrium, not a force summed with the
//! others. The OC denominator is the derivative of PROJECTED volume,
//! `F^T(H' * V / sum(V))`. Trial designs are accepted only after evaluating
//! their actual volume and compliance. History rows describe the same solved
//! design; the last row always describes the returned densities.
//!
//! This is a fixed-mesh density study, not free-boundary CutFEM topology,
//! a KKT certificate, or a mesh-converged/experimentally validated design.
//! Filter, load and adjoint CG solves poll at most every 32 iterations. Work
//! spent on failed or rejected trials counts against the shared linear budget.

use std::ops::ControlFlow;

use crate::control::{EvaluationStop, SolveBudget, SolveControl, SolveWork};
use crate::elasticity::DensityElasticity;
use crate::oc::assert_valid_oc_inputs;
use crate::pipeline::{DesignPipeline, LoadCase, assert_valid_load_cases};

/// Resource and stopping policy for one fixed-continuation study.
#[derive(Debug, Clone, Copy)]
pub struct MultiLoadOcOptions {
    /// Upper bound on projected material volume divided by domain volume.
    pub volume_fraction: f64,
    /// Maximum raw-density change in one accepted iteration.
    pub move_limit: f64,
    /// Maximum number of accepted updates (the initial solve is separate).
    pub max_iterations: usize,
    /// Stop after an accepted update changes every density by at most this.
    /// This is a design-change test, NOT a stationarity certificate.
    pub change_tolerance: f64,
    /// Absolute projected-volume feasibility tolerance.
    pub volume_tolerance: f64,
    /// Maximum number of halved steps after a rejected full trial.
    pub max_backtracks: usize,
}

impl Default for MultiLoadOcOptions {
    fn default() -> Self {
        Self {
            volume_fraction: 0.5,
            move_limit: 0.15,
            max_iterations: 30,
            change_tolerance: 1e-3,
            volume_tolerance: 1e-8,
            max_backtracks: 12,
        }
    }
}

/// Why the driver stopped. No variant asserts optimality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiLoadOcTermination {
    /// The requested number of updates has been used.
    IterationBudget,
    /// An accepted update met the caller's design-change threshold.
    DesignChange,
    /// The caller stopped at a checkpoint; completed work is retained.
    Cancelled,
    /// A per-solve or cumulative linear-work budget was consumed.
    LinearBudget,
    /// A trial evaluation failed numerically; the accepted prefix is retained.
    NumericalFailure,
    /// No feasible, non-increasing-compliance trial was found within budget.
    NoAcceptableStep,
}

/// One real equilibrium evaluation of an accepted design.
#[derive(Debug, Clone)]
pub struct MultiLoadOcIteration {
    /// Zero is the feasible starting design, one is the first update.
    pub iteration: usize,
    /// Weighted sum of independent load compliances.
    pub compliance: f64,
    /// Unweighted load compliances, in the original input order.
    pub case_compliances: Vec<f64>,
    /// Actual projected material volume fraction of this design.
    pub volume_fraction: f64,
    /// Maximum change from the preceding accepted raw design; zero initially.
    pub max_change: f64,
}

/// Last accepted design and its aligned physical iteration history.
#[derive(Debug, Clone)]
pub struct MultiLoadOcReport {
    /// Raw SIMP densities. Projection is deliberately not a binary threshold.
    pub rho: Vec<f64>,
    /// Projected physical densities of the last solved design. Empty when
    /// stopped before the first equilibrium; retained for solver-free export.
    pub projected_rho: Vec<f64>,
    /// Last accepted displacement fields, in load-case order. Empty before
    /// the first equilibrium. Cancellation never substitutes a rejected trial.
    pub displacements: Vec<Vec<f64>>,
    /// Includes the initial solve, unless stopped before any complete solve.
    pub history: Vec<MultiLoadOcIteration>,
    /// Explicit terminal reason.
    pub termination: MultiLoadOcTermination,
    /// Detailed interruption/failure; never a usable partial sensitivity.
    pub evaluation_stop: Option<EvaluationStop>,
    /// Work consumed by the shared control, including rejected/partial trials
    /// and any setup evaluations performed with that same control.
    pub work: SolveWork,
}

fn volume(
    pipeline: &DesignPipeline, rho: &[f64], cell_vol: &[f64], control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    Ok(volume_projection(pipeline, rho, cell_vol, control)?.0)
}

fn volume_projection(
    pipeline: &DesignPipeline, rho: &[f64], cell_vol: &[f64], control: &mut SolveControl<'_>,
) -> Result<(f64, Vec<f64>), EvaluationStop> {
    let (_, projected, _) = pipeline.try_forward(rho, control)?;
    let total: f64 = cell_vol.iter().sum();
    // Normalize before multiplication so representable fractions do not
    // overflow solely because the caller uses a large volume unit.
    let value: f64 = projected.iter().zip(cell_vol).map(|(r, v)| r * (v / total)).sum();
    if !value.is_finite() { return Err(EvaluationStop::Breakdown { stage: "volume" }); }
    Ok((value, projected))
}

/// Log-space OC update at a specified volume multiplier. Log space avoids a
/// fixed dimensional multiplier bracket (which fails when load units change).
fn trial_densities(rho: &[f64], log_ratio: &[f64], log_lambda: f64, step: f64) -> Vec<f64> {
    rho.iter().zip(log_ratio).map(|(&r, &ratio)| {
        let lower = (r - step).max(1e-3);
        let upper = (r + step).min(1.0);
        let exponent = 0.5 * (ratio - log_lambda);
        // Compare before exponentiation: all computed exponentials are then
        // bounded by admissible density ratios instead of overflowing.
        if exponent <= fs_math::det::ln(lower / r) {
            lower
        } else if exponent >= fs_math::det::ln(upper / r) {
            upper
        } else {
            (r * fs_math::det::exp(exponent)).clamp(lower, upper)
        }
    }).collect()
}

/// Optimize with the default component iteration caps and a cancellation hook.
/// The hook is now also polled inside filter, elasticity and adjoint solves.
/// A stopped result retains only fully evaluated accepted designs. Invalid
/// modeling inputs still panic; numerical stops are returned in the report.
#[allow(clippy::too_many_arguments)]
pub fn multi_load_optimality_criteria(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    cell_vol: &[f64],
    options: MultiLoadOcOptions,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> MultiLoadOcReport {
    let mut callback = |_| checkpoint();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    controlled_multi_load_optimality_criteria(pipeline, elasticity, loads, rho0, cell_vol,
        options, &mut control)
}

/// Optimize under one shared cumulative linear-work and cancellation budget.
/// The starting projected design must be feasible under the SAME volume cap
/// and loads. Continuation parameters are fixed throughout this call.
///
/// All stops restore the last accepted operator and fields. An empty history
/// means no complete baseline was solved. Linear-budget exhaustion is not
/// convergence; no trial gradient, including a partially solved load family,
/// is ever accepted. The current budget is observable through `control.work()`.
#[allow(clippy::too_many_arguments)]
pub fn controlled_multi_load_optimality_criteria(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    cell_vol: &[f64],
    options: MultiLoadOcOptions,
    control: &mut SolveControl<'_>,
) -> MultiLoadOcReport {
    assert_valid_load_cases(elasticity, loads);
    assert_valid_oc_inputs(elasticity, loads[0].force, rho0, cell_vol,
        options.volume_fraction, options.move_limit);
    pipeline.params.assert_valid();
    assert!(rho0.iter().all(|r| *r >= 1e-3), "OC starting densities must be at least 1e-3");
    assert!(options.change_tolerance.is_finite() && options.change_tolerance >= 0.0,
        "design-change tolerance must be finite and nonnegative");
    assert!(options.volume_tolerance.is_finite() && options.volume_tolerance >= 0.0
        && options.volume_tolerance < options.volume_fraction,
        "volume tolerance must be finite, nonnegative, and smaller than the cap");
    assert!(options.max_backtracks <= 64, "at most 64 OC backtracks are admitted");
    let mut report = MultiLoadOcReport {
        rho: rho0.to_vec(), projected_rho: Vec::new(), displacements: Vec::new(),
        history: Vec::new(), termination: MultiLoadOcTermination::IterationBudget,
        evaluation_stop: None, work: control.work(),
    };
    // This is also the rollback target for a stop before any equilibrium.
    let mut accepted_moduli = elasticity.moduli.clone();
    let outcome = (|| -> Result<(), EvaluationStop> {
        control.checkpoint("optimizer")?;
        let (initial_volume, initial_projection) = volume_projection(pipeline, rho0, cell_vol, control)?;
        assert!(initial_volume <= options.volume_fraction + options.volume_tolerance,
            "starting projected design exceeds the material budget");
        let mut current = pipeline.try_multi_load_compliance_and_gradient(elasticity, rho0, loads, control)?;
        accepted_moduli = elasticity.moduli.clone();
        report.projected_rho = initial_projection;
        report.displacements = std::mem::take(&mut current.displacements);
        report.history.push(MultiLoadOcIteration {
            iteration: 0, compliance: current.compliance,
            case_compliances: current.case_compliances.clone(),
            volume_fraction: initial_volume, max_change: 0.0,
        });
        for iteration in 0..options.max_iterations {
            control.checkpoint("optimizer")?;
            let (_, dv) = pipeline.try_volume_and_gradient(&report.rho, cell_vol, control)?;
            // OC requires a locally monotone material constraint. Do not hide
            // zero/negative slopes or non-finite sensitivities behind floors.
            if !dv.iter().all(|v| v.is_finite() && *v > 0.0) {
                return Err(EvaluationStop::Breakdown { stage: "volume-gradient" });
            }
            if !current.gradient.iter().all(|g| g.is_finite() && *g <= 0.0) {
                return Err(EvaluationStop::Breakdown { stage: "compliance-gradient" });
            }
            let log_ratio: Vec<f64> = current.gradient.iter().zip(&dv).map(|(&g, &v)| {
                if g < 0.0 { fs_math::det::ln(-g) - fs_math::det::ln(v) }
                else { f64::NEG_INFINITY }
            }).collect();
            let lower: Vec<f64> = report.rho.iter().map(|r| (r - options.move_limit).max(1e-3)).collect();
            let upper: Vec<f64> = report.rho.iter().zip(&log_ratio).map(|(&r, &ratio)| {
                if ratio.is_finite() { (r + options.move_limit).min(1.0) }
                else { (r - options.move_limit).max(1e-3) }
            }).collect();
            let mut candidate = upper;
            if volume(pipeline, &candidate, cell_vol, control)? > options.volume_fraction {
                if volume(pipeline, &lower, cell_vol, control)? > options.volume_fraction {
                    report.termination = MultiLoadOcTermination::NoAcceptableStep;
                    break;
                }
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for ((&r, &ratio), &floor) in report.rho.iter().zip(&log_ratio).zip(&lower) {
                    if ratio.is_finite() {
                        lo = lo.min(ratio - 2.0 * fs_math::det::ln((r + options.move_limit).min(1.0) / r));
                        hi = hi.max(ratio - 2.0 * fs_math::det::ln(floor / r));
                    }
                }
                if !lo.is_finite() || !hi.is_finite() {
                    return Err(EvaluationStop::Breakdown { stage: "multiplier" });
                }
                candidate = lower;
                for _ in 0..80 {
                    control.checkpoint("multiplier")?;
                    let mid = 0.5 * lo + 0.5 * hi;
                    let trial = trial_densities(&report.rho, &log_ratio, mid, options.move_limit);
                    if volume(pipeline, &trial, cell_vol, control)? > options.volume_fraction {
                        lo = mid;
                    } else {
                        hi = mid;
                        candidate = trial; // Retain the actually feasible side.
                    }
                }
            }
            let proposed = candidate;
            let mut accepted = false;
            let mut alpha = 1.0;
            for _ in 0..=options.max_backtracks {
                control.checkpoint("line-search")?;
                let trial: Vec<f64> = report.rho.iter().zip(&proposed)
                    .map(|(&r, &p)| r + alpha * (p - r)).collect();
                let (v, projected) = volume_projection(pipeline, &trial, cell_vol, control)?;
                if v <= options.volume_fraction + options.volume_tolerance {
                    let mut next = pipeline.try_multi_load_compliance_and_gradient(elasticity, &trial, loads, control)?;
                    if next.compliance <= current.compliance {
                        let change = report.rho.iter().zip(&trial)
                            .map(|(r, p)| (r - p).abs()).fold(0.0f64, f64::max);
                        report.rho = trial;
                        report.projected_rho = projected;
                        report.displacements = std::mem::take(&mut next.displacements);
                        report.history.push(MultiLoadOcIteration {
                            iteration: iteration + 1, compliance: next.compliance,
                            case_compliances: next.case_compliances.clone(),
                            volume_fraction: v, max_change: change,
                        });
                        accepted_moduli = elasticity.moduli.clone();
                        current = next;
                        accepted = true;
                        if change <= options.change_tolerance {
                            report.termination = MultiLoadOcTermination::DesignChange;
                        }
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if !accepted {
                report.termination = MultiLoadOcTermination::NoAcceptableStep;
                break;
            }
            if report.termination == MultiLoadOcTermination::DesignChange { break; }
        }
        Ok(())
    })();
    // One exit path covers cancellation and failure at every nested stage,
    // including an adjoint following already completed load solves.
    elasticity.moduli = accepted_moduli;
    report.work = control.work();
    if let Err(stop) = outcome {
        report.termination = match &stop {
            EvaluationStop::Cancelled => MultiLoadOcTermination::Cancelled,
            EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. } => MultiLoadOcTermination::LinearBudget,
            EvaluationStop::Breakdown { .. } => MultiLoadOcTermination::NumericalFailure,
        };
        report.evaluation_stop = Some(stop);
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_update_uses_projected_volume_derivative() {
        // At KKT, -dc_i/dv_i is constant, so every unconstrained density
        // remains fixed. Replacing dv by equal raw volumes fails this test.
        let rho = [0.25, 0.5, 0.75];
        let dv: [f64; 3] = [0.02, 0.1, 0.4];
        let dc: [f64; 3] = [-0.06, -0.3, -1.2];
        let ratio: Vec<f64> = dc.iter().zip(dv)
            .map(|(&c, v)| fs_math::det::ln(-c) - fs_math::det::ln(v)).collect();
        let actual = trial_densities(&rho, &ratio, fs_math::det::ln(3.0), 0.1);
        for (a, e) in actual.iter().zip(rho) { assert!((a - e).abs() < 1e-14); }
    }

    #[test]
    fn multiplier_is_scale_covariant_without_overflow() {
        let rho = [0.1, 0.5, 0.9];
        let ratio = [-2.0, 0.0, 2.0];
        let a = trial_densities(&rho, &ratio, 0.5, 0.15);
        for shift in [-1000.0, 1000.0] {
            let shifted: Vec<f64> = ratio.iter().map(|r| r + shift).collect();
            let b = trial_densities(&rho, &shifted, 0.5 + shift, 0.15);
            for (a, b) in a.iter().zip(b) { assert!((a - b).abs() < 1e-13); }
        }
    }

    #[test]
    fn zero_sensitivity_and_extreme_multiplier_respect_move_limits() {
        let rho = [0.001, 0.5, 1.0];
        let actual = trial_densities(&rho, &[f64::NEG_INFINITY; 3], -1000.0, 0.05);
        assert_eq!(actual, vec![0.001, 0.45, 0.95]);
        let actual = trial_densities(&rho, &[1000.0; 3], -1000.0, 0.05);
        for (a, e) in actual.iter().zip([0.051, 0.55, 1.0]) {
            assert!((a - e).abs() < 1e-14);
        }
    }
}
