//! Discretization and mesh convergence ladder evaluation with Richardson extrapolation.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.6`
//!
//! Provides automated multi-rung mesh refinement analysis, apparent convergence order
//! fitting via deterministic log-space regression, monotonicity checks, ASME V&V 20
//! Grid Convergence Index (GCI) calculation, and conservative evidence-color assignment.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;
use std::fmt::Write as _;

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
    /// Absolute difference between observed and theoretical order.
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

        // Take the 3 finest rungs (h1 < h2 < h3)
        let r1 = &sorted_rungs[0];
        let r2 = &sorted_rungs[1];
        let r3 = &sorted_rungs[2];

        let h1 = r1.h;
        let h2 = r2.h;
        let h3 = r3.h;

        let f1 = r1.qoi_value;
        let f2 = r2.qoi_value;
        let f3 = r3.qoi_value;

        let e21 = f2 - f1;
        let e32 = f3 - f2;

        // Monotonicity check
        let is_monotone = (e21 * e32) > 0.0;
        if !is_monotone {
            let dispersion = (f1 - f2).abs().max((f2 - f3).abs());
            return ConvergenceResult {
                qoi_name: plan.target_qoi.clone(),
                status: ConvergenceStatus::Oscillatory,
                theoretical_order: plan.theoretical_order,
                observed_order: None,
                fit_residual: None,
                richardson_extrapolated_qoi: None,
                discretization_error_gci: Some(dispersion),
                evidence_color: Color::Estimated {
                    estimator: "oscillatory-ladder-envelope".to_string(),
                    dispersion,
                },
                admitted_rungs: sorted_rungs,
                rejection_reason: Some(
                    "non-monotone oscillatory behavior across mesh refinement ladder".to_string(),
                ),
            };
        }

        // Compute refinement ratio r21 and r32
        let r21 = h2 / h1;
        let r32 = h3 / h2;
        let avg_r = (r21 + r32) * 0.5;

        let ratio = (e32 / e21).abs();
        if ratio < 1e-12 || avg_r <= 1.0 + 1e-6 {
            return ConvergenceResult {
                qoi_name: plan.target_qoi.clone(),
                status: ConvergenceStatus::PreAsymptotic,
                theoretical_order: plan.theoretical_order,
                observed_order: None,
                fit_residual: None,
                richardson_extrapolated_qoi: None,
                discretization_error_gci: Some((f1 - f2).abs()),
                evidence_color: Color::Estimated {
                    estimator: "degenerate-refinement-ratio".to_string(),
                    dispersion: (f1 - f2).abs(),
                },
                admitted_rungs: sorted_rungs,
                rejection_reason: Some(
                    "refinement ratio or difference is near machine precision; cannot fit order"
                        .to_string(),
                ),
            };
        }

        // Apparent order of convergence p = ln(|e32/e21|) / ln(avg_r)
        let observed_p = fs_math::det::ln(ratio) / fs_math::det::ln(avg_r);
        if !observed_p.is_finite() || observed_p <= 0.0 {
            return ConvergenceResult {
                qoi_name: plan.target_qoi.clone(),
                status: ConvergenceStatus::PreAsymptotic,
                theoretical_order: plan.theoretical_order,
                observed_order: None,
                fit_residual: None,
                richardson_extrapolated_qoi: None,
                discretization_error_gci: Some((f1 - f2).abs()),
                evidence_color: Color::Estimated {
                    estimator: "non-finite-observed-order".to_string(),
                    dispersion: (f1 - f2).abs(),
                },
                admitted_rungs: sorted_rungs,
                rejection_reason: Some(
                    "observed convergence order is non-positive or non-finite".to_string(),
                ),
            };
        }

        // Richardson extrapolation: f_ext = f1 + (f1 - f2)/(r21^p - 1)
        let r_p = avg_r.powf(observed_p);
        let richardson_ext = if (r_p - 1.0).abs() > 1e-6 {
            Some(f1 + (f1 - f2) / (r_p - 1.0))
        } else {
            None
        };

        // ASME V&V 20 Grid Convergence Index: GCI = (1.25 * |(f1 - f2)/f1|) / (r21^p - 1)
        let rel_diff = if f1.abs() > 1e-12 {
            ((f1 - f2) / f1).abs()
        } else {
            (f1 - f2).abs()
        };
        let gci = if (r_p - 1.0).abs() > 1e-6 {
            Some((1.25 * rel_diff) / (r_p - 1.0).abs())
        } else {
            None
        };

        // Classify asymptotic vs pre-asymptotic
        let order_diff = (observed_p - plan.theoretical_order).abs();
        let status = if order_diff <= 0.5 && observed_p > 0.5 {
            ConvergenceStatus::Asymptotic
        } else {
            ConvergenceStatus::PreAsymptotic
        };

        let evidence_color = match status {
            ConvergenceStatus::Asymptotic => {
                let lo = f1 - gci.unwrap_or(0.01) * f1.abs();
                let hi = f1 + gci.unwrap_or(0.01) * f1.abs();
                Color::Verified {
                    lo: lo.min(hi),
                    hi: lo.max(hi),
                }
            }
            _ => Color::Estimated {
                estimator: "pre-asymptotic-gci".to_string(),
                dispersion: gci.unwrap_or(0.05),
            },
        };

        ConvergenceResult {
            qoi_name: plan.target_qoi.clone(),
            status,
            theoretical_order: plan.theoretical_order,
            observed_order: Some(observed_p),
            fit_residual: Some(order_diff),
            richardson_extrapolated_qoi: richardson_ext,
            discretization_error_gci: gci,
            evidence_color,
            admitted_rungs: sorted_rungs,
            rejection_reason: None,
        }
    }
}
