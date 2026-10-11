//! Inverse-action admission for Euclidean-ball Steihaug PCG.
use super::SteihaugStop;
use std::panic::{AssertUnwindSafe, catch_unwind};

/// A refused inverse action. The accepted optimizer state remains resumable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PreconditionerError {
    /// The cumulative allowance has already been spent.
    Budget {
        /// Actual attempted applications, including discarded work.
        used: usize,
        /// Supplied cumulative ceiling.
        limit: usize,
    },
    /// The inverse action did not return one value per residual coordinate.
    Dimension {
        /// Required length.
        expected: usize,
        /// Returned length.
        actual: usize,
    },
    /// The returned action, before or after normalization, is not finite.
    NonFinite {
        /// First invalid coordinate.
        component: usize,
    },
    /// `r^T M^-1 r` was not finite and strictly positive for a nonzero residual.
    /// This is a local numerical refusal, not a global definiteness certificate.
    NonPositiveProduct,
    /// The inverse callback unwound. Its external side effects are not undone.
    CallbackPanicked,
}
impl core::fmt::Display for PreconditionerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Budget { used, limit } => write!(f, "preconditioner budget exhausted: {used}/{limit}"),
            Self::Dimension { expected, actual } => write!(f, "preconditioner needs {expected} values, received {actual}"),
            Self::NonFinite { component } => write!(f, "non-finite preconditioner coordinate {component}"),
            Self::NonPositiveProduct => write!(f, "preconditioner residual product is not finite and positive"),
            Self::CallbackPanicked => write!(f, "preconditioner callback panicked; accepted point retained"),
        }
    }
}
impl std::error::Error for PreconditionerError {}

pub(super) struct Action<'a> {
    callback: crate::FnHv<'a>,
    point: &'a [f64],
    used: &'a mut usize,
    limit: usize,
    normalization: Option<f64>,
}
impl<'a> Action<'a> {
    pub(super) fn new(callback: crate::FnHv<'a>, point: &'a [f64], used: &'a mut usize, limit: usize) -> Self {
        Self { callback, point, used, limit, normalization: None }
    }

    pub(super) fn apply(&mut self, residual: &[f64], hv_count: usize,
        paused: &mut dyn FnMut() -> bool) -> Result<Vec<f64>, SteihaugStop>
    {
        if paused() { return Err(SteihaugStop::Paused(hv_count)); }
        if *self.used >= self.limit {
            return Err(SteihaugStop::Preconditioner(PreconditionerError::Budget {
                used: *self.used, limit: self.limit,
            }));
        }
        *self.used += 1; // Precharge, including a callback that unwinds.
        let result = catch_unwind(AssertUnwindSafe(|| (self.callback)(self.point, residual)));
        let mut z = result.map_err(|_| SteihaugStop::Preconditioner(PreconditionerError::CallbackPanicked))?;
        if paused() { return Err(SteihaugStop::Paused(hv_count)); }
        if z.len() != residual.len() {
            return Err(SteihaugStop::Preconditioner(PreconditionerError::Dimension {
                expected: residual.len(), actual: z.len(),
            }));
        }
        finite(&z)?;
        // One constant multiplier for the ENTIRE Krylov solve: changing the
        // normalization at each residual would violate fixed-preconditioner CG.
        let scale = *self.normalization.get_or_insert_with(|| super::quadratic_scale(&z));
        if scale != 1.0 {
            for value in &mut z { *value /= scale; }
        }
        finite(&z)?;
        Ok(z)
    }
}
fn finite(z: &[f64]) -> Result<(), SteihaugStop> {
    match z.iter().position(|v| !v.is_finite()) {
        Some(component) => Err(SteihaugStop::Preconditioner(PreconditionerError::NonFinite { component })),
        None => Ok(()),
    }
}
pub(super) fn positive_product(r: &[f64], z: &[f64]) -> Result<f64, SteihaugStop> {
    let rz: f64 = r.iter().zip(z).map(|(r, z)| r * z).sum();
    if rz.is_finite() && rz > 0.0 { Ok(rz) }
    else { Err(SteihaugStop::Preconditioner(PreconditionerError::NonPositiveProduct)) }
}
