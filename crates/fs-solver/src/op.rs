//! The matrix-free operator abstraction. P6 doctrine: solvers see
//! APPLIES, never entries — assembled matrices enter only through the
//! [`CsrOp`] adapter (and the p-MG coarse level).

use fs_sparse::{Csr, ops};

/// A linear operator with its transpose — the adjoint hook is part of
/// the trait, not an afterthought, so every solver built on it is
/// adjoint-equipped by construction.
pub trait LinearOp {
    /// Dimension (square operators in v1).
    fn n(&self) -> usize;
    /// y = A·x.
    fn apply(&self, x: &[f64], y: &mut [f64]);
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
}
