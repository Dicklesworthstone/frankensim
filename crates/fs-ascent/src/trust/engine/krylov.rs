//! Inexact Newton work/accuracy policy, distinct from hard cumulative budgets.

/// Invalid accuracy or per-iteration Krylov work policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KrylovOptionsError {
    /// The relative residual tolerance must be finite and strictly in (0, 1).
    RelativeTolerance,
    /// At least one Krylov product must be allowed per attempted Newton step.
    ZeroProducts,
}
impl core::fmt::Display for KrylovOptionsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::RelativeTolerance => write!(f, "Krylov relative tolerance must be finite and in (0, 1)"),
            Self::ZeroProducts => write!(f, "Krylov per-iteration product cap must be positive"),
        }
    }
}
impl std::error::Error for KrylovOptionsError {}

/// Accuracy and work for one inexact Newton subproblem. Reaching this local
/// cap returns the current approximate step to the existing model/acceptance
/// checks, rather than discarding it as a cumulative-budget exhaustion does.
///
/// The cap counts only Krylov Hessian products. The final model check needs
/// ONE ADDITIONAL product from the run's hard cumulative Hessian allowance.
/// Negative curvature and original-coordinate trust-boundary truncation stay
/// active. The original gradient still controls
/// the outer stopping rule; an approximate inner solution is not convergence.
///
/// The numerical ceiling is also capped at twice the design dimension, as in
/// the original kernel. A preconditioned solve spends no inverse action just
/// to prepare an unused direction after the configured local product cap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KrylovOptions {
    relative_tolerance: f64,
    max_products: usize,
}
impl Default for KrylovOptions {
    fn default() -> Self { Self { relative_tolerance: 1e-8, max_products: usize::MAX } }
}
impl KrylovOptions {
    /// Admit a policy before solver work. The default is the existing `1e-8`
    /// residual tolerance and the kernel's dimension-based product ceiling.
    pub fn new(relative_tolerance: f64, max_products: usize) -> Result<Self, KrylovOptionsError> {
        if !relative_tolerance.is_finite() || relative_tolerance <= 0.0 || relative_tolerance >= 1.0 {
            return Err(KrylovOptionsError::RelativeTolerance);
        }
        if max_products == 0 { return Err(KrylovOptionsError::ZeroProducts); }
        Ok(Self { relative_tolerance, max_products })
    }
    /// Residual tolerance relative to the initial, unpreconditioned residual.
    pub fn relative_tolerance(self) -> f64 { self.relative_tolerance }
    /// Local Krylov-only cap; the final model product is additional.
    pub fn max_products(self) -> usize { self.max_products }
}
