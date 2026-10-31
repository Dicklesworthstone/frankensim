//! Cooperative cancellation for the existing deterministic complex LU owner.
use super::{C64, LuComplex};

/// Registered capsule: the factor and solve bodies below recompiled with
/// native fused multiply-add on x86-64 CPUs that have it.
mod fma;
/// Blocked, thread-parallel elimination with the same operation graph.
mod blocked;

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

/// The complex product exactly as `C64`'s `Mul` forms it (two fused
/// multiply-adds, same operands and order). `C64`'s operator is not inlined
/// across crates, so the FMA capsule would otherwise keep a per-product call.
#[inline(always)]
#[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fns
fn product(a: C64, b: C64) -> C64 {
    C64::new(a.re.mul_add(b.re, -(a.im * b.im)), a.re.mul_add(b.im, a.im * b.re))
}

/// `C64`'s `Sub`, inlined for the same reason.
#[inline(always)]
#[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fns
fn difference(a: C64, b: C64) -> C64 {
    C64::new(a.re - b.re, a.im - b.im)
}

/// Elimination shared by the portable and the FMA-compiled paths: the original
/// magnitude pivots, lowest-index ties, operation order and poll points.
#[inline(always)]
#[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fns
fn factor_body(
    m: &mut [C64], perm: &mut [usize], n: usize, cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LuComplexError> {
    for k in 0..n {
        poll(cancelled)?;
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
            poll(cancelled)?;
            let f = m[i * n + k] / piv;
            m[i * n + k] = f;
            let (upper, lower) = m.split_at_mut(i * n);
            let pivot_row = &upper[k * n + k + 1..k * n + n];
            for (target, &source) in lower[k + 1..n].iter_mut().zip(pivot_row) {
                *target = difference(*target, product(f, source));
            }
        }
    }
    Ok(())
}

/// Forward and backward substitution shared by both paths, polling per row.
#[inline(always)]
#[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fns
fn solve_body(
    lu: &LuComplex, b: &[C64], cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, LuComplexError> {
    let n = lu.n;
    let mut x: Vec<C64> = lu.perm.iter().map(|&p| b[p]).collect();
    for i in 0..n {
        poll(cancelled)?;
        let mut v = x[i];
        for (k, &bk) in x.iter().enumerate().take(i) {
            v = difference(v, product(lu.data[i * n + k], bk));
        }
        x[i] = v;
    }
    for i in (0..n).rev() {
        poll(cancelled)?;
        let mut v = x[i];
        for (k, &bk) in x.iter().enumerate().take(n).skip(i + 1) {
            v = difference(v, product(lu.data[i * n + k], bk));
        }
        x[i] = v / lu.data[i * n + i];
    }
    Ok(x)
}

/// Factor a row-major complex matrix with cooperative cancellation.
///
/// Uses the same magnitude pivots, lowest-index tie breaks, and scalar operation
/// order as [`super::lu_complex`]. Polls before copying input, at every pivot,
/// before each elimination row, and before publishing the factor. Between
/// elimination polls at most one O(n) row is updated, not an O(n^2) pivot block.
/// Allocation/copy latency and OS scheduling are not wall-time bounded here.
/// On x86-64 CPUs with FMA the elimination runs with native fused instructions;
/// the result is bit-identical (each product is the same correctly rounded fma).
///
/// The immutable input is unchanged on every outcome. No partial factor escapes.
/// The callback can borrow a caller's cancellation gate; no executor or thread
/// identity is introduced into the numerical result.
///
/// Inside [`super::with_lu_workers`] with more than one declared worker, this
/// runs [`lu_complex_blocked_with_cancel`] on them instead: the same bits, with
/// that function's poll points. Undeclared, the poll points above hold.
///
/// # Errors
/// [`LuComplexError::Cancelled`] or [`LuComplexError::Singular`].
/// # Panics
/// If n*n overflows or `a` does not contain exactly n*n entries.
pub fn lu_complex_with_cancel(
    a: &[C64], n: usize, mut cancelled: impl FnMut() -> bool,
) -> Result<LuComplex, LuComplexError> {
    let workers = super::LU_WORKERS.with(core::cell::Cell::get);
    if workers > 1 {
        return lu_complex_blocked_with(a, n, blocked::BLOCK, workers, cancelled);
    }
    let entries = n.checked_mul(n).expect("complex matrix extent overflow");
    assert_eq!(a.len(), entries, "a must be n*n = {entries}");
    poll(&mut cancelled)?;
    let mut m = a.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    fma::factor(&mut m, &mut perm, n, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(LuComplex { n, data: m, perm })
}

/// Logical temporary slice-descriptor payload of the fixed-block complex LU.
///
/// This is additional to the input matrix, copied factor, pivot permutation
/// and caller vectors. The first U12 phase bounds every later phase: it has
/// BLOCK rows and the widest trailing region. Includes the per-worker inner
/// tuple arrays, outer Vec descriptors and split edges. Excludes allocator
/// bookkeeping, thread stacks and runtime metadata; it is not an RSS bound.
/// Returns None when the logical byte count cannot be represented.
pub fn lu_complex_blocked_workspace_bytes(n: usize, threads: usize) -> Option<usize> {
    let width = n.saturating_sub(blocked::BLOCK);
    if width == 0 {
        return Some(0);
    }
    let parts = threads.max(1).min(width);
    let tuples = blocked::BLOCK
        .checked_mul(parts)?
        .checked_mul(core::mem::size_of::<(&[C64], &mut [C64])>())?;
    let outer = parts.checked_mul(core::mem::size_of::<Vec<(&[C64], &mut [C64])>>())?;
    let edges = parts
        .checked_add(1)?
        .checked_mul(core::mem::size_of::<usize>())?;
    tuples.checked_add(outer)?.checked_add(edges)
}

/// Factor with the blocked, thread-parallel elimination.
///
/// The result is bit-identical to [`lu_complex_with_cancel`] for every thread
/// count: the same magnitude pivots, lowest-index ties and singular index, and
/// each element's same ascending-pivot sequence of fused updates. Only the
/// order in which different elements are updated changes. That order is
/// cache-blocked, and `threads` workers (at least one) share each block's
/// upper-block and trailing updates, so the result never depends on them.
///
/// Polls before copying input, at every pivot, before each block's upper and
/// trailing updates, and before publishing. Between polls at most one block's
/// update of the trailing matrix runs, about 64 (n - k)^2 complex updates
/// shared by the workers. The immutable input is unchanged on every outcome,
/// and no partial factor escapes.
///
/// # Errors
/// [`LuComplexError::Cancelled`] or [`LuComplexError::Singular`].
/// # Panics
/// If n*n overflows or `a` does not contain exactly n*n entries.
pub fn lu_complex_blocked_with_cancel(
    a: &[C64], n: usize, threads: usize, cancelled: impl FnMut() -> bool,
) -> Result<LuComplex, LuComplexError> {
    lu_complex_blocked_with(a, n, blocked::BLOCK, threads, cancelled)
}

/// [`lu_complex_blocked_with_cancel`] with an explicit block size, for tests.
fn lu_complex_blocked_with(
    a: &[C64], n: usize, block: usize, threads: usize, mut cancelled: impl FnMut() -> bool,
) -> Result<LuComplex, LuComplexError> {
    let entries = n.checked_mul(n).expect("complex matrix extent overflow");
    assert_eq!(a.len(), entries, "a must be n*n = {entries}");
    poll(&mut cancelled)?;
    let mut m = a.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    blocked::factor_blocked(&mut m, &mut perm, n, block, threads, &mut cancelled)?;
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
        assert_eq!(b.len(), self.n, "b must have length {}", self.n);
        poll(&mut cancelled)?;
        let x = fma::solve(self, b, &mut cancelled)?;
        poll(&mut cancelled)?;
        b.copy_from_slice(&x);
        Ok(())
    }
}

#[cfg(test)]
mod tests;
