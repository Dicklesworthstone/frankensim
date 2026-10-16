//! Cooperative cancellation for the existing deterministic complex LU owner.
use super::{C64, LuComplex};

/// A cancelled operation is distinct from a numerically singular factorization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LuComplexError {
    /// The caller requested cancellation; no factor or solution is published.
    Cancelled,
    /// The magnitude pivot is zero at this elimination step.
    Singular {
        /// Zero-based elimination step, identical to `EigFailure::window_hi`.
        index: usize,
    },
}
impl core::fmt::Display for LuComplexError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => f.write_str("complex LU cancelled"),
            Self::Singular { index } => write!(f, "complex LU singular at step {index}"),
        }
    }
}
impl std::error::Error for LuComplexError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LuComplexError> {
    if cancelled() { Err(LuComplexError::Cancelled) } else { Ok(()) }
}

/// Factor a row-major complex matrix with cooperative cancellation.
///
/// Uses the same magnitude pivots, lowest-index tie breaks, and scalar operation
/// order as [`super::lu_complex`]. Polls before copying input, at every pivot,
/// before each elimination row, and before publishing the factor. Between
/// elimination polls at most one O(n) row is updated, not an O(n^2) pivot block.
/// Allocation/copy latency and OS scheduling are not wall-time bounded here.
///
/// The immutable input is unchanged on every outcome. No partial factor escapes.
/// The callback can borrow a caller's cancellation gate; no executor or thread
/// identity is introduced into the numerical result.
///
/// # Errors
/// [`LuComplexError::Cancelled`] or [`LuComplexError::Singular`].
/// # Panics
/// If n*n overflows or `a` does not contain exactly n*n entries.
pub fn lu_complex_with_cancel(
    a: &[C64], n: usize, mut cancelled: impl FnMut() -> bool,
) -> Result<LuComplex, LuComplexError> {
    let entries = n.checked_mul(n).expect("complex matrix extent overflow");
    assert_eq!(a.len(), entries, "a must be n*n = {entries}");
    poll(&mut cancelled)?;
    let mut m = a.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    for k in 0..n {
        poll(&mut cancelled)?;
        let (mut best, mut best_i) = (m[k * n + k].abs(), k);
        for i in k + 1..n {
            let v = m[i * n + k].abs();
            if v > best {
                best = v;
                best_i = i;
            }
        }
        if best == 0.0 {
            return Err(LuComplexError::Singular { index: k });
        }
        if best_i != k {
            for c in 0..n {
                m.swap(k * n + c, best_i * n + c);
            }
            perm.swap(k, best_i);
        }
        let piv = m[k * n + k];
        for i in k + 1..n {
            poll(&mut cancelled)?;
            let f = m[i * n + k] / piv;
            m[i * n + k] = f;
            for c in k + 1..n {
                let sub = f * m[k * n + c];
                m[i * n + c] = m[i * n + c] - sub;
            }
        }
    }
    poll(&mut cancelled)?;
    Ok(LuComplex { n, data: m, perm })
}

impl LuComplex {
    /// Solve A*x=b with cancellation and atomic publication into `b`.
    ///
    /// The existing O(n) permuted-RHS workspace stages the entire solve. Checks
    /// occur before allocation, at each forward/backward row and immediately
    /// before copying the completed solution back. Cancellation (or a callback
    /// panic) leaves `b` bitwise unchanged; the factor remains reusable.
    /// Neither reduction order nor floating-point operations depend on polling.
    ///
    /// # Errors
    /// [`LuComplexError::Cancelled`]; no partial solution is published.
    /// # Panics
    /// If `b` does not have the factor's dimension.
    pub fn solve_with_cancel(
        &self, b: &mut [C64], mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), LuComplexError> {
        let n = self.n;
        assert_eq!(b.len(), n, "b must have length {n}");
        poll(&mut cancelled)?;
        let mut x: Vec<C64> = self.perm.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut v = x[i];
            for (k, &bk) in x.iter().enumerate().take(i) {
                let sub = self.data[i * n + k] * bk;
                v = v - sub;
            }
            x[i] = v;
        }
        for i in (0..n).rev() {
            poll(&mut cancelled)?;
            let mut v = x[i];
            for (k, &bk) in x.iter().enumerate().take(n).skip(i + 1) {
                let sub = self.data[i * n + k] * bk;
                v = v - sub;
            }
            x[i] = v / self.data[i * n + i];
        }
        poll(&mut cancelled)?;
        b.copy_from_slice(&x);
        Ok(())
    }
}

#[cfg(test)]
mod tests;
