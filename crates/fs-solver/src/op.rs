//! The matrix-free operator abstraction. P6 doctrine: solvers see
//! APPLIES, never entries — assembled matrices enter only through the
//! [`CsrOp`] adapter (and the p-MG coarse level).

use fs_sparse::{Csr, SparseApplyError, ops};
use std::fmt;

pub mod multilevel;
pub mod two_level;

/// A controlled operator action refused before the solver commits its work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperatorApplyError {
    /// The cooperative cancellation predicate requested a pause.
    Interrupted,
    /// The operator's dense shape cannot be represented by `usize`.
    DimensionOverflow,
    /// An input or output buffer does not match the operator's shape.
    DimensionMismatch {
        /// `"input"` or `"output"`.
        buffer: &'static str,
        /// Required scalar length.
        expected: usize,
        /// Supplied scalar length.
        actual: usize,
    },
    /// The operator could not reserve its temporary output.
    Allocation,
}

impl fmt::Display for OperatorApplyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Interrupted => f.write_str("operator apply interrupted"),
            Self::DimensionOverflow => f.write_str("operator apply dimensions overflow usize"),
            Self::DimensionMismatch {
                buffer,
                expected,
                actual,
            } => {
                write!(
                    f,
                    "operator apply {buffer} length: expected {expected}, got {actual}"
                )
            }
            Self::Allocation => {
                f.write_str("operator apply could not reserve its temporary output")
            }
        }
    }
}

impl std::error::Error for OperatorApplyError {}

impl From<SparseApplyError<()>> for OperatorApplyError {
    fn from(error: SparseApplyError<()>) -> Self {
        match error {
            SparseApplyError::Interrupted(()) => Self::Interrupted,
            SparseApplyError::DimensionOverflow => Self::DimensionOverflow,
            SparseApplyError::DimensionMismatch {
                buffer,
                expected,
                actual,
            } => Self::DimensionMismatch {
                buffer,
                expected,
                actual,
            },
            SparseApplyError::Allocation => Self::Allocation,
        }
    }
}

/// A linear operator with its transpose — the adjoint hook is part of
/// the trait, not an afterthought, so every solver built on it is
/// adjoint-equipped by construction.
pub trait LinearOp {
    /// Dimension (square operators in v1).
    fn n(&self) -> usize;
    /// y = A·x.
    fn apply(&self, x: &[f64], y: &mut [f64]);
    /// Cooperative action for cancellable solvers. `true` requests a pause.
    ///
    /// The default checks before and after the existing action. Override to
    /// poll during long kernels. The caller must discard `y` on refusal:
    /// arbitrary operators need not make their partial output transactional.
    /// [`CsrOp`] does stage its output and checks within sparse rows.
    ///
    /// # Errors
    /// Returns interruption or a structured operator refusal. The predicate
    /// must be bounded and must not panic.
    fn apply_controlled(
        &self,
        x: &[f64],
        y: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), OperatorApplyError> {
        if cancelled() {
            return Err(OperatorApplyError::Interrupted);
        }
        self.apply(x, y);
        if cancelled() {
            return Err(OperatorApplyError::Interrupted);
        }
        Ok(())
    }
    /// y = Aᵀ·x. The default forwards to `apply` — CORRECT ONLY for
    /// symmetric operators; nonsymmetric implementations MUST override.
    ///
    /// NOT CHECKED. This crate has no symmetry probe: nothing at compile
    /// time or run time detects a nonsymmetric operator that inherited
    /// this default. The transposed-solve battery pins the crate's own
    /// adapters (`CsrOp::general`, which materializes a real transpose,
    /// and one nonsymmetric fs-opdsl fixture that does override); it
    /// cannot cover a client operator it never sees. An implementor who
    /// inherits the default on a nonsymmetric operator gets a transposed
    /// solve of A instead of Aᵀ that converges normally, and the
    /// reported residual is measured against the operator supplied —
    /// not against its transpose — so a small residual is no evidence
    /// that the adjoint system was the one solved. Verify with an
    /// explicit ⟨Av, w⟩ = ⟨v, Aᵀw⟩ probe before relying on it.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }

    /// Cooperative transpose action, with the contract of
    /// [`Self::apply_controlled`]. The default preserves an implementation's
    /// existing transpose action; it does not assume symmetry.
    ///
    /// # Errors
    /// Returns interruption or a structured operator refusal.
    fn apply_transpose_controlled(
        &self,
        x: &[f64],
        y: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), OperatorApplyError> {
        if cancelled() {
            return Err(OperatorApplyError::Interrupted);
        }
        self.apply_transpose(x, y);
        if cancelled() {
            return Err(OperatorApplyError::Interrupted);
        }
        Ok(())
    }
}

/// Assembled-matrix adapter. For nonsymmetric matrices the transpose
/// is materialized once at construction (deterministic fs-sparse
/// transpose).
pub struct CsrOp {
    a: Csr,
    at: Option<Csr>,
}

impl CsrOp {
    /// Wrap a symmetric matrix (transpose = self).
    #[must_use]
    pub fn symmetric(a: Csr) -> CsrOp {
        CsrOp { a, at: None }
    }

    /// Wrap a general matrix (transpose materialized).
    #[must_use]
    pub fn general(a: Csr) -> CsrOp {
        let at = ops::transpose(&a);
        CsrOp { a, at: Some(at) }
    }

    /// The wrapped matrix.
    #[must_use]
    pub fn matrix(&self) -> &Csr {
        &self.a
    }
}

impl LinearOp for CsrOp {
    fn n(&self) -> usize {
        self.a.nrows()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.a.spmv(x, y);
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        match &self.at {
            Some(at) => at.spmv(x, y),
            None => self.a.spmv(x, y),
        }
    }

    fn apply_controlled(
        &self,
        x: &[f64],
        y: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), OperatorApplyError> {
        self.a
            .spmv_with_checkpoint(x, y, || if cancelled() { Err(()) } else { Ok(()) })
            .map_err(OperatorApplyError::from)
    }

    fn apply_transpose_controlled(
        &self,
        x: &[f64],
        y: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), OperatorApplyError> {
        self.at
            .as_ref()
            .unwrap_or(&self.a)
            .spmv_with_checkpoint(x, y, || if cancelled() { Err(()) } else { Ok(()) })
            .map_err(OperatorApplyError::from)
    }
}
