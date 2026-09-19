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
//! Checkpoints occur between component evaluations and multiplier trials;
//! the existing filter and elasticity solves are NOT interruptible internally.

use std::ops::ControlFlow;

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
    /// Includes the initial solve, unless cancelled before any solve.
    pub history: Vec<MultiLoadOcIteration>,
    /// Explicit terminal reason.
    pub termination: MultiLoadOcTermination,
}

fn volume(pipeline: &DesignPipeline, rho: &[f64], cell_vol: &[f64]) -> f64 {
    let (_, projected, _) = pipeline.forward(rho);
    let total: f64 = cell_vol.iter().sum();
    // Normalize before multiplication so representable fractions do not
    // overflow solely because the caller uses a large volume unit.
    let value: f64 = projected.iter().zip(cell_vol).map(|(r, v)| r * (v / total)).sum();
    assert!(value.is_finite(), "projected material volume must be finite");
    value
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

/// Optimize one fixed-mesh design against independently applied load cases.
///
/// `checkpoint` returns `Break(())` to stop. A stopped result retains only
/// evaluated, feasible designs; cancellation before the initial solve has an
/// empty history and must not be presented as a solved result. Components
/// retain their existing panic-on-invalid-model/failed-solve contract.
///
/// The starting projected design must be feasible, so comparisons are against
/// a design under the SAME volume cap and loads, never full-material compliance.
/// Continuation must be explicit: changing projection parameters changes the
/// physical volume constraint and requires a newly admitted starting design.
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
        rho: rho0.to_vec(), history: Vec::new(),
        termination: MultiLoadOcTermination::IterationBudget,
    };
    if checkpoint().is_break() {
        report.termination = MultiLoadOcTermination::Cancelled;
        return report;
    }
    let initial_volume = volume(pipeline, rho0, cell_vol);
    assert!(initial_volume <= options.volume_fraction + options.volume_tolerance,
        "starting projected design exceeds the material budget");
    if checkpoint().is_break() {
        report.termination = MultiLoadOcTermination::Cancelled;
        return report;
    }
    let mut current = pipeline.multi_load_compliance_and_gradient(elasticity, rho0, loads);
    report.history.push(MultiLoadOcIteration {
        iteration: 0, compliance: current.compliance,
        case_compliances: current.case_compliances.clone(),
        volume_fraction: initial_volume, max_change: 0.0,
    });
    for iteration in 0..options.max_iterations {
        if checkpoint().is_break() {
            report.termination = MultiLoadOcTermination::Cancelled;
            break;
        }
        let (_, dv) = pipeline.volume_and_gradient(&report.rho, cell_vol);
        // Classical multiplicative OC is not valid for a locally nonmonotone
        // material constraint. Do not hide negative/zero slopes behind a floor.
        assert!(dv.iter().all(|v| v.is_finite() && *v > 0.0),
            "OC requires a strictly positive projected-volume gradient");
        let log_ratio: Vec<f64> = current.gradient.iter().zip(&dv).map(|(&g, &v)| {
            assert!(g.is_finite() && g <= 0.0,
                "OC requires finite nonpositive compliance sensitivities");
            if g < 0.0 { fs_math::det::ln(-g) - fs_math::det::ln(v) }
            else { f64::NEG_INFINITY }
        }).collect();
        let lower: Vec<f64> = report.rho.iter().map(|r| (r - options.move_limit).max(1e-3)).collect();
        let upper: Vec<f64> = report.rho.iter().zip(&log_ratio).map(|(&r, &ratio)| {
            if ratio.is_finite() { (r + options.move_limit).min(1.0) }
            else { (r - options.move_limit).max(1e-3) }
        }).collect();
        if checkpoint().is_break() {
            report.termination = MultiLoadOcTermination::Cancelled;
            break;
        }
        let mut candidate = upper;
        if volume(pipeline, &candidate, cell_vol) > options.volume_fraction {
            if volume(pipeline, &lower, cell_vol) > options.volume_fraction {
                report.termination = MultiLoadOcTermination::NoAcceptableStep;
                break;
            }
            // Endpoints correspond to every sensitive cell at its upper/lower
            // move limit. They scale with the actual objective, not its units.
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for ((&r, &ratio), &floor) in report.rho.iter().zip(&log_ratio).zip(&lower) {
                if ratio.is_finite() {
                    lo = lo.min(ratio - 2.0 * fs_math::det::ln((r + options.move_limit).min(1.0) / r));
                    hi = hi.max(ratio - 2.0 * fs_math::det::ln(floor / r));
                }
            }
            assert!(lo.is_finite() && hi.is_finite(), "OC multiplier bracket must be finite");
            candidate = lower;
            for _ in 0..80 {
                if checkpoint().is_break() {
                    report.termination = MultiLoadOcTermination::Cancelled;
                    return report;
                }
                let mid = 0.5 * lo + 0.5 * hi;
                let trial = trial_densities(&report.rho, &log_ratio, mid, options.move_limit);
                if volume(pipeline, &trial, cell_vol) > options.volume_fraction {
                    lo = mid;
                } else {
                    hi = mid;
                    candidate = trial; // Always retain the actually feasible side.
                }
            }
        }
        let proposed = candidate;
        let accepted_moduli = elasticity.moduli.clone();
        let mut accepted = false;
        let mut alpha = 1.0;
        for _ in 0..=options.max_backtracks {
            if checkpoint().is_break() {
                elasticity.moduli = accepted_moduli;
                report.termination = MultiLoadOcTermination::Cancelled;
                return report;
            }
            let trial: Vec<f64> = report.rho.iter().zip(&proposed)
                .map(|(&r, &p)| r + alpha * (p - r)).collect();
            let v = volume(pipeline, &trial, cell_vol);
            if v <= options.volume_fraction + options.volume_tolerance {
                if checkpoint().is_break() {
                    elasticity.moduli = accepted_moduli;
                    report.termination = MultiLoadOcTermination::Cancelled;
                    return report;
                }
                let next = pipeline.multi_load_compliance_and_gradient(elasticity, &trial, loads);
                if next.compliance <= current.compliance {
                    let change = report.rho.iter().zip(&trial)
                        .map(|(r, p)| (r - p).abs()).fold(0.0f64, f64::max);
                    report.rho = trial;
                    report.history.push(MultiLoadOcIteration {
                        iteration: iteration + 1, compliance: next.compliance,
                        case_compliances: next.case_compliances.clone(),
                        volume_fraction: v, max_change: change,
                    });
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
            elasticity.moduli = accepted_moduli;
            report.termination = MultiLoadOcTermination::NoAcceptableStep;
            break;
        }
        if report.termination == MultiLoadOcTermination::DesignChange { break; }
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
