//! Reusable dense partial-pivot solve for small repeated systems.
//!
//! Unlike the blocked [`crate::factor::lu`] path, this kernel never packs GEMM panels.
//! Its unblocked FMA traversal is deterministic but is NOT bit-interchangeable
//! with the blocked factorization. The existing factor bit contract is unchanged.

/// Refusal from a prepared dense solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LuWorkspaceError {
    /// Storage lengths do not match the prepared square system.
    Dimension,
    /// The requested workspace could not be represented or allocated.
    Capacity,
    /// An input or intermediate is not finite.
    NonFinite,
    /// No successful factorization is available (including after a failed replacement).
    Unfactored,
    /// No nonzero pivot exists in the indicated column.
    Singular {
        /// Zero-based pivot column.
        index: usize,
    },
}

impl core::fmt::Display for LuWorkspaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prepared LU: {self:?}")
    }
}

impl std::error::Error for LuWorkspaceError {}

/// Reusable small-dense solve storage, allocated only at construction.
///
/// Row pivot ties choose the lowest row. [`Self::solve_into`] publishes no output
/// until both factorization and substitution succeed. Failed calls may overwrite
/// scratch, but never change the caller's matrix, right-hand side or output.
/// This is an unblocked algorithm; it does not promise a speed advantage at large n.
#[derive(Debug)]
pub struct LuWorkspace {
    n: usize,
    factors: Vec<f64>,
    rhs: Vec<f64>,
    pivots: Vec<usize>,
    factored: bool,
}

fn zeroed(len: usize) -> Result<Vec<f64>, LuWorkspaceError> {
    let mut values = Vec::new();
    values.try_reserve_exact(len).map_err(|_| LuWorkspaceError::Capacity)?;
    values.resize(len, 0.0);
    Ok(values)
}

impl LuWorkspace {
    /// Allocate for an n-by-n row-major matrix. An empty system is supported.
    ///
    /// # Errors
    /// Returns [`LuWorkspaceError::Capacity`] on extent overflow or allocation failure.
    pub fn new(n: usize) -> Result<Self, LuWorkspaceError> {
        let entries = n.checked_mul(n).ok_or(LuWorkspaceError::Capacity)?;
        let mut pivots = Vec::new();
        pivots.try_reserve_exact(n).map_err(|_| LuWorkspaceError::Capacity)?;
        pivots.resize(n, 0);
        Ok(Self { n, factors: zeroed(entries)?, rhs: zeroed(n)?, pivots, factored: false })
    }

    /// Fixed matrix dimension.
    #[must_use]
    pub const fn dimension(&self) -> usize { self.n }

    /// Factor a fresh matrix and retain its row permutation for repeated RHSs.
    /// No allocation occurs. A failed replacement invalidates the old factor;
    /// it cannot accidentally solve a new physical system with stale coefficients.
    ///
    /// # Errors
    /// Bad extent, nonfinite data/intermediates, or an exactly zero pivot.
    #[allow(clippy::needless_range_loop)]
    pub fn factor(&mut self, a: &[f64]) -> Result<(), LuWorkspaceError> {
        self.factored = false;
        let n = self.n;
        if a.len() != self.factors.len() { return Err(LuWorkspaceError::Dimension); }
        if a.iter().any(|v| !v.is_finite()) { return Err(LuWorkspaceError::NonFinite); }
        self.factors.copy_from_slice(a);
        let m = &mut self.factors;
        for k in 0..n {
            let mut pivot_row = k;
            let mut pivot_abs = m[k * n + k].abs();
            for row in k + 1..n {
                let candidate = m[row * n + k].abs();
                if candidate > pivot_abs { pivot_abs = candidate; pivot_row = row; }
            }
            if !pivot_abs.is_finite() { return Err(LuWorkspaceError::NonFinite); }
            if pivot_abs == 0.0 { return Err(LuWorkspaceError::Singular { index: k }); }
            self.pivots[k] = pivot_row;
            if pivot_row != k {
                for col in 0..n { m.swap(k * n + col, pivot_row * n + col); }
            }
            let pivot = m[k * n + k];
            for row in k + 1..n {
                let multiplier = m[row * n + k] / pivot;
                if !multiplier.is_finite() { return Err(LuWorkspaceError::NonFinite); }
                m[row * n + k] = multiplier;
                for col in k + 1..n {
                    m[row * n + col] = (-multiplier).mul_add(m[k * n + col], m[row * n + col]);
                }
            }
        }
        if m.iter().any(|v| !v.is_finite()) { return Err(LuWorkspaceError::NonFinite); }
        self.factored = true;
        Ok(())
    }

    /// Solve against the last successful factorization, without refactoring.
    /// All RHSs see the same pivots. Failed RHS admission/substitution leaves
    /// both the caller's output and the retained factor unchanged.
    ///
    /// # Errors
    /// Missing valid factor, bad extent, or nonfinite data/solution.
    #[allow(clippy::needless_range_loop)]
    pub fn solve_factored_into(&mut self, b: &[f64], out: &mut [f64])
        -> Result<(), LuWorkspaceError>
    {
        let n = self.n;
        if b.len() != n || out.len() != n { return Err(LuWorkspaceError::Dimension); }
        if !self.factored { return Err(LuWorkspaceError::Unfactored); }
        if b.iter().any(|v| !v.is_finite()) { return Err(LuWorkspaceError::NonFinite); }
        self.rhs.copy_from_slice(b);
        for k in 0..n { self.rhs.swap(k, self.pivots[k]); }
        let m = &self.factors;
        for row in 0..n {
            let mut value = self.rhs[row];
            for col in 0..row { value = (-m[row * n + col]).mul_add(self.rhs[col], value); }
            self.rhs[row] = value;
        }
        for row in (0..n).rev() {
            let mut value = self.rhs[row];
            for col in row + 1..n { value = (-m[row * n + col]).mul_add(self.rhs[col], value); }
            self.rhs[row] = value / m[row * n + row];
        }
        if self.rhs.iter().any(|v| !v.is_finite()) { return Err(LuWorkspaceError::NonFinite); }
        out.copy_from_slice(&self.rhs);
        Ok(())
    }

    /// Factor a fresh matrix and solve `a * x = b` into caller-owned storage.
    /// No allocation, deallocation, resizing, locks or I/O occurs. Arithmetic
    /// and lowest-row pivot tie breaking retain the original one-shot path.
    ///
    /// # Errors
    /// Refuses bad lengths, nonfinite data/intermediates and singular pivots.
    /// `out` remains bit-for-bit unchanged on every refusal.
    pub fn solve_into(&mut self, a: &[f64], b: &[f64], out: &mut [f64])
        -> Result<(), LuWorkspaceError>
    {
        if b.len() != self.n || out.len() != self.n {
            return Err(LuWorkspaceError::Dimension);
        }
        self.factor(a)?;
        self.solve_factored_into(b, out)
    }

}
