//! Matrix-free LSQR for rectangular, inconsistent, and rank-deficient systems.
//!
//! Golub--Kahan bidiagonalization with Givens updates (Paige and Saunders,
//! ACM TOMS 8, 1982). No normal matrix is formed. Each accepted iteration
//! measures BOTH `||b - A x||` and `||A^T (b - A x)||` using fresh applies;
//! least-squares stationarity is not confused with an exact fit. These are
//! floating-point diagnostics, not certified error or condition bounds.
//!
//! A clone is an in-memory checkpoint. Drivers cancel between calls to
//! [`LsqrState::step`]; a failed step leaves the checkpoint unchanged. Resume
//! with the SAME immutable operator (dimensions are checked, identity is the
//! caller's responsibility). Reproducibility additionally requires deterministic
//! operator applies. The initial iterate is zero, retaining the minimum-norm
//! Krylov construction for underdetermined systems in exact arithmetic.

use crate::{RectLinearOp, dot};
use std::fmt;

/// Explicit stopping tolerances for the two differently dimensioned residuals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LsqrConfig {
    /// Relative data-residual tolerance, scaled by `||b||`.
    pub relative_tolerance: f64,
    /// Absolute data-residual tolerance, in the units of `b`.
    pub absolute_tolerance: f64,
    /// Relative stationarity tolerance, scaled by the INITIAL `||A^T b||`.
    pub normal_relative_tolerance: f64,
    /// Absolute stationarity tolerance, in the units of `A^T b`.
    pub normal_absolute_tolerance: f64,
    /// Maximum accepted bidiagonalization steps; zero permits diagnosis only.
    pub max_iters: usize,
}

impl Default for LsqrConfig {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            absolute_tolerance: 0.0,
            normal_relative_tolerance: 1e-10,
            normal_absolute_tolerance: 0.0,
            max_iters: 1000,
        }
    }
}

/// Why an LSQR state stopped; exhaustion is never reported as convergence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LsqrStop {
    /// The true data residual met its tolerance.
    Compatible,
    /// The true normal residual met its tolerance, but the data residual did not.
    LeastSquares,
    /// The caller's iteration budget was exhausted.
    IterationLimit,
    /// Bidiagonalization exhausted its direction without meeting either target.
    Breakdown,
}

/// Freshly evaluated residuals at one accepted iterate (including iteration zero).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LsqrResidual {
    /// Euclidean norm of `b - A x`.
    pub residual_norm: f64,
    /// Euclidean norm of `A^T (b - A x)`.
    pub normal_residual_norm: f64,
}

/// Diagnostics for a paused or completed solve.
#[derive(Debug, Clone, PartialEq)]
pub struct LsqrReport {
    /// Number of accepted steps, not operator applications.
    pub iterations: usize,
    /// `None` means paused and still runnable.
    pub stop: Option<LsqrStop>,
    /// Initial residual followed by one entry per accepted step.
    pub history: Vec<LsqrResidual>,
}

impl LsqrReport {
    /// Whether either explicitly configured convergence target was met.
    #[must_use]
    pub fn converged(&self) -> bool {
        matches!(self.stop, Some(LsqrStop::Compatible | LsqrStop::LeastSquares))
    }
}

/// Solution plus the measured stopping evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct LsqrSolution {
    /// Final coefficient vector.
    pub x: Vec<f64>,
    /// Residual history and termination reason.
    pub report: LsqrReport,
}

/// Invalid input or a non-finite numerical operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LsqrError {
    /// Both dimensions must be nonzero.
    Empty,
    /// A vector or resumed operator did not have the declared dimension.
    Dimension {
        /// Mismatched dimension.
        name: &'static str,
        /// Required dimension.
        expected: usize,
        /// Supplied dimension.
        actual: usize,
    },
    /// Invalid nonnegative finite tolerance (relative tolerances must be below one).
    InvalidTolerance(&'static str),
    /// Input, operator output, or arithmetic was non-finite.
    NonFinite(&'static str),
}

impl fmt::Display for LsqrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("LSQR requires nonzero rows and columns"),
            Self::Dimension {
                name,
                expected,
                actual,
            } => write!(f, "LSQR {name}: expected {expected}, observed {actual}"),
            Self::InvalidTolerance(name) => write!(f, "invalid LSQR tolerance: {name}"),
            Self::NonFinite(stage) => write!(f, "non-finite LSQR value in {stage}"),
        }
    }
}

impl std::error::Error for LsqrError {}

/// Complete, cloneable LSQR state with no hidden operator or executor state.
#[derive(Debug, Clone, PartialEq)]
pub struct LsqrState {
    config: LsqrConfig,
    b: Vec<f64>,
    x: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    w: Vec<f64>,
    alpha: f64,
    rhobar: f64,
    phibar: f64,
    initial: LsqrResidual,
    report: LsqrReport,
}

impl LsqrState {
    /// Initialize from zero, validating inputs and measuring initial stationarity.
    pub fn new(
        a: &dyn RectLinearOp,
        b: &[f64],
        config: LsqrConfig,
    ) -> Result<Self, LsqrError> {
        if a.rows() == 0 || a.cols() == 0 {
            return Err(LsqrError::Empty);
        }
        dimension("right-hand side", a.rows(), b.len())?;
        for (name, value, relative) in [
            ("relative_tolerance", config.relative_tolerance, true),
            ("absolute_tolerance", config.absolute_tolerance, false),
            (
                "normal_relative_tolerance",
                config.normal_relative_tolerance,
                true,
            ),
            (
                "normal_absolute_tolerance",
                config.normal_absolute_tolerance,
                false,
            ),
        ] {
            if !value.is_finite() || value < 0.0 || (relative && value >= 1.0) {
                return Err(LsqrError::InvalidTolerance(name));
            }
        }
        let beta = checked_norm(b, "right-hand side")?;
        let mut normal = vec![0.0; a.cols()];
        a.apply_transpose(b, &mut normal);
        let initial = LsqrResidual {
            residual_norm: beta,
            normal_residual_norm: checked_norm(&normal, "initial transpose")?,
        };
        let mut u = b.to_vec();
        normalize(&mut u, beta);
        let mut v = vec![0.0; a.cols()];
        a.apply_transpose(&u, &mut v);
        let alpha = checked_norm(&v, "initial bidiagonalization")?;
        normalize(&mut v, alpha);
        let stop = stopping(initial, initial, config)
            .or_else(|| (config.max_iters == 0).then_some(LsqrStop::IterationLimit));
        Ok(Self {
            config,
            b: b.to_vec(),
            x: vec![0.0; a.cols()],
            u,
            w: v.clone(),
            v,
            alpha,
            rhobar: alpha,
            phibar: beta,
            initial,
            report: LsqrReport {
                iterations: 0,
                stop,
                history: vec![initial],
            },
        })
    }

    /// Current solution; updated only after a complete finite iteration.
    #[must_use]
    pub fn x(&self) -> &[f64] {
        &self.x
    }

    /// Borrow the diagnostics without copying the history.
    #[must_use]
    pub fn report(&self) -> &LsqrReport {
        &self.report
    }

    /// Advance one iteration. Return `true` only when an iterate was accepted.
    /// A stopped state is a no-op; an error leaves ALL state unchanged.
    pub fn step(&mut self, a: &dyn RectLinearOp) -> Result<bool, LsqrError> {
        dimension("operator rows", self.b.len(), a.rows())?;
        dimension("operator columns", self.x.len(), a.cols())?;
        if self.report.stop.is_some() {
            return Ok(false);
        }
        // Keep the old checkpoint intact until every numerical check passes.
        let mut u = vec![0.0; self.u.len()];
        a.apply(&self.v, &mut u);
        for (next, old) in u.iter_mut().zip(&self.u) {
            *next -= self.alpha * old;
        }
        let beta = checked_norm(&u, "bidiagonal forward apply")?;
        normalize(&mut u, beta);
        let mut v = vec![0.0; self.v.len()];
        if beta != 0.0 {
            a.apply_transpose(&u, &mut v);
            for (next, old) in v.iter_mut().zip(&self.v) {
                *next -= beta * old;
            }
        }
        let alpha = checked_norm(&v, "bidiagonal transpose apply")?;
        normalize(&mut v, alpha);
        let rho = checked_norm(&[self.rhobar, beta], "Givens rotation")?;
        if rho == 0.0 {
            self.report.stop = Some(LsqrStop::Breakdown);
            return Ok(false);
        }
        let c = self.rhobar / rho;
        let s = beta / rho;
        let theta = s * alpha;
        let rhobar = -c * alpha;
        let phi = c * self.phibar;
        let phibar = s * self.phibar;
        let mut x = self.x.clone();
        let mut w = v.clone();
        let step = phi / rho;
        let turn = theta / rho;
        finite(&[step, turn, rhobar, phibar], "Givens update")?;
        for ((xi, wi), old) in x.iter_mut().zip(&mut w).zip(&self.w) {
            *xi += step * old;
            *wi -= turn * old;
        }
        finite(&x, "solution update")?;
        finite(&w, "direction update")?;
        let residual = measure(a, &self.b, &x)?;
        let iterations = self.report.iterations + 1;
        let stop = stopping(residual, self.initial, self.config)
            .or_else(|| (iterations >= self.config.max_iters).then_some(LsqrStop::IterationLimit))
            .or_else(|| (alpha == 0.0 || beta == 0.0).then_some(LsqrStop::Breakdown));
        self.x = x;
        self.u = u;
        self.v = v;
        self.w = w;
        self.alpha = alpha;
        self.rhobar = rhobar;
        self.phibar = phibar;
        self.report.iterations = iterations;
        self.report.history.push(residual);
        self.report.stop = stop;
        Ok(true)
    }

    /// Run to a stopping reason, using the remaining iteration budget.
    pub fn run(&mut self, a: &dyn RectLinearOp) -> Result<(), LsqrError> {
        while self.step(a)? {}
        Ok(())
    }

    /// Consume a paused or completed state without another operator application.
    #[must_use]
    pub fn into_solution(self) -> LsqrSolution {
        LsqrSolution {
            x: self.x,
            report: self.report,
        }
    }
}

/// Solve a rectangular least-squares problem from zero.
pub fn lsqr(
    a: &dyn RectLinearOp,
    b: &[f64],
    config: LsqrConfig,
) -> Result<LsqrSolution, LsqrError> {
    let mut state = LsqrState::new(a, b, config)?;
    state.run(a)?;
    Ok(state.into_solution())
}

fn dimension(name: &'static str, expected: usize, actual: usize) -> Result<(), LsqrError> {
    if expected != actual {
        return Err(LsqrError::Dimension {
            name,
            expected,
            actual,
        });
    }
    Ok(())
}

fn finite(values: &[f64], stage: &'static str) -> Result<(), LsqrError> {
    if values.iter().any(|v| !v.is_finite()) {
        return Err(LsqrError::NonFinite(stage));
    }
    Ok(())
}

// Scale BEFORE squaring: finite vectors such as [1e200] and [1e-200]
// must not produce an infinite or zero norm merely because of the square.
fn checked_norm(values: &[f64], stage: &'static str) -> Result<f64, LsqrError> {
    finite(values, stage)?;
    let scale = values.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    if scale == 0.0 {
        return Ok(0.0);
    }
    let scaled: Vec<f64> = values.iter().map(|v| v / scale).collect();
    let norm = scale * fs_math::det::sqrt(dot(&scaled, &scaled));
    finite(&[norm], stage)?;
    Ok(norm)
}

fn normalize(values: &mut [f64], norm: f64) {
    if norm != 0.0 {
        for v in values {
            *v /= norm;
        }
    }
}

fn measure(a: &dyn RectLinearOp, b: &[f64], x: &[f64]) -> Result<LsqrResidual, LsqrError> {
    let mut r = vec![0.0; b.len()];
    a.apply(x, &mut r);
    for (ri, bi) in r.iter_mut().zip(b) {
        *ri = bi - *ri;
    }
    let residual_norm = checked_norm(&r, "true residual")?;
    let mut normal = vec![0.0; x.len()];
    a.apply_transpose(&r, &mut normal);
    Ok(LsqrResidual {
        residual_norm,
        normal_residual_norm: checked_norm(&normal, "true normal residual")?,
    })
}

// Difference form avoids overflow in absolute + relative * baseline.
fn within(value: f64, baseline: f64, relative: f64, absolute: f64) -> bool {
    value <= absolute || value - absolute <= relative * baseline
}

fn stopping(r: LsqrResidual, initial: LsqrResidual, c: LsqrConfig) -> Option<LsqrStop> {
    if within(
        r.residual_norm,
        initial.residual_norm,
        c.relative_tolerance,
        c.absolute_tolerance,
    ) {
        Some(LsqrStop::Compatible)
    } else if within(
        r.normal_residual_norm,
        initial.normal_residual_norm,
        c.normal_relative_tolerance,
        c.normal_absolute_tolerance,
    ) {
        Some(LsqrStop::LeastSquares)
    } else {
        None
    }
}
