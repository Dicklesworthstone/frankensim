//! Discretization and mesh convergence ladder evaluation with Richardson extrapolation.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.6`
//!
//! Provides automated multi-rung mesh refinement analysis, apparent convergence order
//! fitting via deterministic all-rung least squares, monotonicity checks, ASME V&V 20
//! Grid Convergence Index (GCI) calculation, and conservative evidence-color assignment.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;
use std::fmt::Write as _;

use crate::axis::{AnalysisPolicy, AxisVerdict};

/// All-rung leading-power least-squares diagnostics. The coefficient test
/// implements the constant E/h^p asymptotic requirement; it never averages
/// unequal refinement ratios. This is the leading-power classification seam
/// of Eça & Hoekstra (JCP 262, 2014), not their full multi-model uncertainty
/// procedure or an uncertainty certificate for scattered data.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerFit {
    /// Number of points included in the fit.
    pub points: usize,
    /// Positive fitted power, if identifiable in the declared search interval.
    pub observed_order: Option<f64>,
    /// Fitted zero-resolution intercept.
    pub continuum: Option<f64>,
    /// RMS fit residual divided by the full observable range.
    pub relative_rms: Option<f64>,
    /// Maximum relative deviation of adjacent error coefficients.
    pub coefficient_spread: Option<f64>,
    /// Classification and optional relative fine-grid GCI.
    pub verdict: AxisVerdict,
}

// At fixed p solve the two-parameter least-squares problem y = a + b x^p.
// Both coordinates and observations are scaled before entering this function.
fn linear_power_fit(data: &[(f64, f64)], p: f64) -> Option<(f64, f64, f64)> {
    let n = data.len() as f64;
    let mean_x = data
        .iter()
        .map(|(h, _)| fs_math::det::pow(*h, p))
        .sum::<f64>()
        / n;
    let mean_y = data.iter().map(|(_, y)| *y).sum::<f64>() / n;
    let mut variance = 0.0;
    let mut covariance = 0.0;
    for &(h, y) in data {
        let x = fs_math::det::pow(h, p) - mean_x;
        variance += x * x;
        covariance += x * (y - mean_y);
    }
    if !variance.is_finite() || variance <= f64::MIN_POSITIVE {
        return None;
    }
    let b = covariance / variance;
    let a = mean_y - b * mean_x;
    let sse = data
        .iter()
        .map(|(h, y)| {
            let e = y - a - b * fs_math::det::pow(*h, p);
            e * e
        })
        .sum::<f64>();
    (a.is_finite() && b.is_finite() && sse.is_finite()).then_some((a, b, sse))
}

/// Rungs within `ROUNDOFF_ULPS` ulps of their largest magnitude are flat:
/// `Unresolved`, with no observable order.
pub(crate) const ROUNDOFF_ULPS: f64 = 128.0;

/// Fit `f(h)=f0+C h^p` to ALL supplied points. Positive finite h values
/// must be ordered coarse to fine. Observations may have either sign.
///
/// GCI uses Fs=1.25 and the actual finest pair's ratio, only after monotonicity,
/// all-rung least-squares fit, formal-order agreement, and adjacent coefficient
/// consistency pass. Zero/roundoff-sized fine values have no relative GCI.
#[must_use]
pub fn fit_power_law(data: &[(f64, f64)], formal_order: f64, policy: &AnalysisPolicy) -> PowerFit {
    let mut result = PowerFit {
        points: data.len(),
        observed_order: None,
        continuum: None,
        relative_rms: None,
        coefficient_spread: None,
        verdict: AxisVerdict::InsufficientRungs,
    };
    if data.len() < 3 {
        return result;
    }
    if !policy.valid()
        || !formal_order.is_finite()
        || formal_order <= 0.0
        || formal_order > policy.maximum_order
        || data
            .iter()
            .any(|(h, y)| !h.is_finite() || *h <= 0.0 || !y.is_finite())
        || data
            .windows(2)
            .any(|w| w[0].0 / w[1].0 < policy.minimum_refinement_ratio)
    {
        result.verdict = AxisVerdict::Refused {
            reason: "invalid algebraic fit declaration or resolution sequence".into(),
        };
        return result;
    }
    let lo = data.iter().map(|(_, y)| *y).fold(f64::INFINITY, f64::min);
    let hi = data
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::NEG_INFINITY, f64::max);
    let scale = hi - lo;
    let roundoff = ROUNDOFF_ULPS * f64::EPSILON * lo.abs().max(hi.abs()).max(f64::MIN_POSITIVE);
    if !scale.is_finite() || scale <= roundoff {
        result.verdict = AxisVerdict::Unresolved;
        return result;
    }
    let scaled: Vec<_> = data
        .iter()
        .map(|(h, y)| (*h / data[0].0, (*y - lo) / scale))
        .collect();
    let differences: Vec<_> = scaled.windows(2).map(|w| w[0].1 - w[1].1).collect();
    let monotone = differences.iter().all(|d| *d > 0.0) || differences.iter().all(|d| *d < 0.0);
    // Deterministic bounded variable projection: coarse bracketing followed by
    // golden-section minimization. Every objective evaluation includes all rungs.
    let minimum = 0.05;
    let step = (policy.maximum_order - minimum) / 128.0;
    let mut best = 0_usize;
    let mut best_sse = f64::INFINITY;
    for i in 0..=128_usize {
        let p = minimum + i as f64 * step;
        let sse = linear_power_fit(&scaled, p).map_or(f64::INFINITY, |(_, _, s)| s);
        if sse < best_sse {
            best = i;
            best_sse = sse;
        }
    }
    let mut left = minimum + best.saturating_sub(1) as f64 * step;
    let mut right = minimum + (best + 1).min(128) as f64 * step;
    let golden = (fs_math::det::sqrt(5.0) - 1.0) * 0.5;
    for _ in 0..96 {
        let a = right - golden * (right - left);
        let b = left + golden * (right - left);
        let fa = linear_power_fit(&scaled, a).map_or(f64::INFINITY, |(_, _, s)| s);
        let fb = linear_power_fit(&scaled, b).map_or(f64::INFINITY, |(_, _, s)| s);
        if fa <= fb {
            right = b;
        } else {
            left = a;
        }
    }
    let p = 0.5 * (left + right);
    let Some((a, b, sse)) = linear_power_fit(&scaled, p) else {
        result.verdict = AxisVerdict::Scattered;
        return result;
    };
    let rms = fs_math::det::sqrt(sse / data.len() as f64);
    result.observed_order = Some(p);
    result.continuum = Some(lo + scale * a);
    result.relative_rms = Some(rms);
    let coefficient_spread = scaled
        .windows(2)
        .map(|w| {
            let c =
                (w[0].1 - w[1].1) / (fs_math::det::pow(w[0].0, p) - fs_math::det::pow(w[1].0, p));
            ((c - b) / b).abs()
        })
        .fold(0.0_f64, f64::max);
    result.coefficient_spread = Some(coefficient_spread);
    if !monotone {
        let alternating = differences.windows(2).all(|w| w[0] * w[1] < 0.0);
        result.verdict = if alternating {
            AxisVerdict::Oscillatory
        } else {
            AxisVerdict::Scattered
        };
    } else if rms > policy.relative_fit_tolerance {
        result.verdict = AxisVerdict::Scattered;
    } else if coefficient_spread > policy.asymptotic_tolerance
        || (p - formal_order).abs() > policy.order_tolerance
        || p <= minimum + step * 1e-6
        || p >= policy.maximum_order - step * 1e-6
        || data.windows(2).any(|w| (w[0].1 - w[1].1).abs() <= roundoff)
    {
        result.verdict = AxisVerdict::PreAsymptotic { observed_order: p };
    } else {
        let fine = data[data.len() - 1];
        let medium = data[data.len() - 2];
        let denominator = fs_math::det::pow(medium.0 / fine.0, p) - 1.0;
        let gci = (fine.1.abs() > roundoff && denominator > 0.0)
            .then(|| 1.25 * ((medium.1 - fine.1) / fine.1).abs() / denominator);
        result.verdict = AxisVerdict::Asymptotic {
            observed_order: p,
            gci: gci.filter(|g| g.is_finite()),
        };
    }
    result
}

/// Status of the discretization convergence evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvergenceStatus {
    /// Monotone asymptotic convergence observed with order near theoretical expectation.
    Asymptotic,
    /// Pre-asymptotic mesh regime (observed order significantly lower than theoretical).
    PreAsymptotic,
    /// Oscillatory or non-monotone QoI behavior across refinement rungs.
    Oscillatory,
    /// Fewer than three valid rungs completed; insufficient for order fitting or Richardson extrapolation.
    InsufficientRungs,
    /// Solver non-convergence or fatal failure on one or more rungs.
    FailedSolve,
}

impl ConvergenceStatus {
    /// Machine-readable status label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Asymptotic => "asymptotic",
            Self::PreAsymptotic => "pre-asymptotic",
            Self::Oscillatory => "oscillatory",
            Self::InsufficientRungs => "insufficient-rungs",
            Self::FailedSolve => "failed-solve",
        }
    }
}

/// A single completed mesh refinement rung in the convergence study.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshRung {
    /// Rung index (0 = coarsest mesh).
    pub ordinal: usize,
    /// Mesh identifier or asset name.
    pub mesh_id: String,
    /// Characteristic mesh element spacing.
    pub h: f64,
    /// Physical unit of mesh spacing (e.g. "m", "mm").
    pub h_unit: String,
    /// Total degrees of freedom in the solved system.
    pub dof: usize,
    /// Solver termination verdict ("converged", "exhausted", "failed").
    pub solver_status: String,
    /// Final solver residual norm.
    pub solver_residual: f64,
    /// Target Quantity of Interest name.
    pub qoi_name: String,
    /// Evaluated Quantity of Interest value on this mesh.
    pub qoi_value: f64,
    /// Physical unit of the Quantity of Interest.
    pub qoi_unit: String,
    /// Computation time consumed in seconds.
    pub budget_consumed_s: f64,
}

impl MeshRung {
    /// Create a new mesh rung.
    #[must_use]
    #[allow(clippy::too_many_arguments)] // Flat constructor binds all rung metadata fields directly.
    pub fn new(
        ordinal: usize,
        mesh_id: impl Into<String>,
        h: f64,
        h_unit: impl Into<String>,
        dof: usize,
        qoi_name: impl Into<String>,
        qoi_value: f64,
        qoi_unit: impl Into<String>,
    ) -> Self {
        Self {
            ordinal,
            mesh_id: mesh_id.into(),
            h,
            h_unit: h_unit.into(),
            dof,
            solver_status: "converged".to_string(),
            solver_residual: 1e-8,
            qoi_name: qoi_name.into(),
            qoi_value,
            qoi_unit: qoi_unit.into(),
            budget_consumed_s: 1.0,
        }
    }

    /// Whether this rung solved successfully with finite values.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.solver_status == "converged"
            && self.h.is_finite()
            && self.h > 0.0
            && self.qoi_value.is_finite()
    }
}

/// Complete plan and data for evaluating mesh convergence on a target QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvergencePlan {
    /// Name of the target physical quantity.
    pub target_qoi: String,
    /// Theoretical asymptotic convergence order (e.g. 2.0 for P1 elements in L2).
    pub theoretical_order: f64,
    /// Refinement ladder rungs ordered coarse to fine.
    pub rungs: Vec<MeshRung>,
}

impl ConvergencePlan {
    /// Create a new plan with declared theoretical convergence order.
    #[must_use]
    pub fn new(target_qoi: impl Into<String>, theoretical_order: f64) -> Self {
        Self {
            target_qoi: target_qoi.into(),
            theoretical_order,
            rungs: Vec::new(),
        }
    }

    /// Add a mesh rung to the plan.
    #[must_use]
    pub fn with_rung(mut self, rung: MeshRung) -> Self {
        self.rungs.push(rung);
        self
    }
}

/// The result of evaluating a discretization convergence ladder.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvergenceResult {
    /// Name of the evaluated Quantity of Interest.
    pub qoi_name: String,
    /// Overall convergence regime classification.
    pub status: ConvergenceStatus,
    /// Expected theoretical order.
    pub theoretical_order: f64,
    /// Numerically fitted apparent order of convergence.
    pub observed_order: Option<f64>,
    /// All-rung RMS fit residual divided by the observed data range.
    pub fit_residual: Option<f64>,
    /// Richardson extrapolated continuum value (h -> 0).
    pub richardson_extrapolated_qoi: Option<f64>,
    /// Grid Convergence Index (ASME V&V 20) discretization uncertainty estimate.
    pub discretization_error_gci: Option<f64>,
    /// Resulting evidence color classification.
    pub evidence_color: Color,
    /// Subset of rungs admitted into the fit.
    pub admitted_rungs: Vec<MeshRung>,
    /// Rejection reason if fitting failed or was refused.
    pub rejection_reason: Option<String>,
}

impl ConvergenceResult {
    /// Generate a deterministic BLAKE3 digest of the convergence evidence.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut text = String::with_capacity(1024);
        let _ = write!(
            text,
            "qoi={};status={};theor={:.4};obs={:?};extrap={:?};gci={:?};rungs={}",
            self.qoi_name,
            self.status.label(),
            self.theoretical_order,
            self.observed_order,
            self.richardson_extrapolated_qoi,
            self.discretization_error_gci,
            self.admitted_rungs.len()
        );
        for r in &self.admitted_rungs {
            let _ = write!(
                text,
                "|rung={}:{}:{:.6}:{:.6}",
                r.ordinal, r.mesh_id, r.h, r.qoi_value
            );
        }
        hash_domain("org.frankensim.convergence.result.v1", text.as_bytes())
    }
}

/// Evaluator for mesh discretization ladders.
pub struct ConvergenceEvaluator;

impl ConvergenceEvaluator {
    /// Evaluate the convergence plan across all supplied rungs.
    #[must_use]
    #[allow(clippy::too_many_lines)] // Single end-to-end Richardson/GCI evaluation sequence.
    pub fn evaluate(plan: &ConvergencePlan) -> ConvergenceResult {
        let valid_rungs: Vec<MeshRung> = plan
            .rungs
            .iter()
            .filter(|r| r.is_valid())
            .cloned()
            .collect();

        // Sort by h ascending (finest mesh first: h1 < h2 < h3)
        let mut sorted_rungs = valid_rungs;
        sorted_rungs.sort_by(|a, b| a.h.partial_cmp(&b.h).unwrap_or(std::cmp::Ordering::Equal));

        if sorted_rungs.len() < 3 {
            return ConvergenceResult {
                qoi_name: plan.target_qoi.clone(),
                status: ConvergenceStatus::InsufficientRungs,
                theoretical_order: plan.theoretical_order,
                observed_order: None,
                fit_residual: None,
                richardson_extrapolated_qoi: None,
                discretization_error_gci: None,
                evidence_color: Color::Estimated {
                    estimator: "insufficient-ladder".to_string(),
                    dispersion: 0.20,
                },
                admitted_rungs: sorted_rungs,
                rejection_reason: Some(
                    "at least three valid completed rungs are required to fit convergence order"
                        .to_string(),
                ),
            };
        }

        // Preserve the legacy result shape and finest-first admitted order,
        // while sharing the all-rung, nonuniform analysis with axis studies.
        let data: Vec<_> = sorted_rungs
            .iter()
            .rev()
            .map(|r| (r.h, r.qoi_value))
            .collect();
        let policy = AnalysisPolicy {
            maximum_order: 32.0,
            order_tolerance: 0.5,
            relative_fit_tolerance: 0.01,
            asymptotic_tolerance: 0.05,
            minimum_refinement_ratio: 1.1,
        };
        let fit = fit_power_law(&data, plan.theoretical_order, &policy);
        let (status, gci) = match fit.verdict {
            AxisVerdict::Asymptotic { gci, .. } => (ConvergenceStatus::Asymptotic, gci),
            AxisVerdict::Oscillatory | AxisVerdict::Scattered => {
                (ConvergenceStatus::Oscillatory, None)
            }
            _ => (ConvergenceStatus::PreAsymptotic, None),
        };
        let fine = sorted_rungs[0].qoi_value;
        let dispersion = data
            .iter()
            .map(|(_, y)| (*y - fine).abs())
            .fold(0.0, f64::max);
        let evidence_color = if let Some(gci) = gci {
            // Retain the legacy color contract. Axis reports intentionally make
            // no Verified claim from an empirical Richardson estimate.
            let half_width = gci * fine.abs();
            Color::Verified {
                lo: fine - half_width,
                hi: fine + half_width,
            }
        } else {
            Color::Estimated {
                estimator: "all-rung-convergence-no-gci".into(),
                dispersion,
            }
        };
        ConvergenceResult {
            qoi_name: plan.target_qoi.clone(),
            status,
            theoretical_order: plan.theoretical_order,
            observed_order: fit.observed_order,
            fit_residual: fit.relative_rms,
            richardson_extrapolated_qoi: if status == ConvergenceStatus::Asymptotic {
                fit.continuum
            } else {
                None
            },
            discretization_error_gci: gci,
            evidence_color,
            admitted_rungs: sorted_rungs,
            rejection_reason: (status != ConvergenceStatus::Asymptotic)
                .then(|| format!("{:?}", fit.verdict)),
        }
    }
}
