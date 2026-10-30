//! Reusable adjacent-pivot LU for real tridiagonal systems.
//!
//! At each elimination step the larger of the current diagonal and the next
//! row's subdiagonal is selected. Equal magnitudes retain the current row.
//! Swapping adjacent rows introduces only a second superdiagonal, so both
//! factorization and substitution use O(n) work and storage. The elimination
//! structure is described by the LAPACK Users' Guide and DGTTRF documentation:
//! <https://netlib.org/lapack/explore-html/d6/d46/group__gttrf_ga8d1e46216e6c861c89bd4328b8be52a1.html>.
//! This is a native safe Rust implementation with no external numerical runtime.
//!
//! Storage is fallibly allocated once. Repeated factor/solve calls allocate
//! nothing, use fixed scalar order and fused multiply-add, and poll cancellation
//! every 128 rows. A solve stages its answer in retained scratch and checks
//! cancellation before publishing it. The final output copy is an uninterrupted
//! commit; a failed/cancelled solve never exposes a partially overwritten RHS.
//! Exactly zero pivots and nonfinite intermediates refuse. There is no arbitrary
//! small-pivot cutoff, condition estimate, scaling or forward-error certificate;
//! a physical/nonlinear caller must still check its original residual.

/// Admission, numerical or cancellation refusal from the structured solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TridiagonalError {
    /// Expected lengths are n-1, n, n-1 for lower/diagonal/upper, and n for RHSs.
    Dimension,
    /// A workspace extent or allocation could not be represented.
    Capacity,
    /// An input, factor or substituted value is not finite.
    NonFinite,
    /// No successful factorization is retained, including after a failed replacement.
    Unfactored,
    /// No nonzero pivot remains in this zero-based column.
    Singular { index: usize },
    /// Cancellation observed before factor or solution publication.
    Cancelled,
}
impl core::fmt::Display for TridiagonalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tridiagonal LU: {self:?}")
    }
}
impl std::error::Error for TridiagonalError {}

/// Fixed-size factors, adjacent pivots and transactional RHS scratch.
///
/// All matrix and RHS inputs remain unchanged. Failed replacement factorizations
/// invalidate the old factors, preventing accidental reuse for a new system.
/// Failed solves retain valid factors and leave the caller's output unchanged.
#[derive(Debug)]
pub struct TridiagonalWorkspace {
    diagonal: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    upper_two: Vec<f64>,
    swapped: Vec<bool>,
    rhs: Vec<f64>,
    factored: bool,
}

fn filled<T: Clone>(len: usize, value: T) -> Result<Vec<T>, TridiagonalError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| TridiagonalError::Capacity)?;
    values.resize(len, value);
    Ok(values)
}
fn finite(value: f64) -> Result<f64, TridiagonalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TridiagonalError::NonFinite)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TridiagonalError> {
    if cancelled() {
        Err(TridiagonalError::Cancelled)
    } else {
        Ok(())
    }
}

impl TridiagonalWorkspace {
    /// Fallibly allocate O(n) storage. Empty and scalar systems are supported.
    pub fn new(n: usize) -> Result<Self, TridiagonalError> {
        // A conservative aggregate byte bound also rejects huge extents before
        // any individually representable but pointless partial allocation.
        n.checked_mul(5 * size_of::<f64>() + size_of::<bool>())
            .filter(|&bytes| bytes <= isize::MAX as usize)
            .ok_or(TridiagonalError::Capacity)?;
        let adjacent = n.saturating_sub(1);
        Ok(Self {
            diagonal: filled(n, 0.0)?,
            lower: filled(adjacent, 0.0)?,
            upper: filled(adjacent, 0.0)?,
            upper_two: filled(n.saturating_sub(2), 0.0)?,
            swapped: filled(adjacent, false)?,
            rhs: filled(n, 0.0)?,
            factored: false,
        })
    }

    pub fn dimension(&self) -> usize {
        self.diagonal.len()
    }

    /// Retain a new partial-pivot factorization without allocating.
    /// `lower[i] = A[i+1,i]`, `diagonal[i] = A[i,i]`, `upper[i] = A[i,i+1]`.
    /// Input lengths and every value are admitted before elimination. Any
    /// refusal, including early cancellation, invalidates the previous factor.
    pub fn factor(
        &mut self,
        lower: &[f64],
        diagonal: &[f64],
        upper: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TridiagonalError> {
        self.factored = false;
        poll(&mut cancelled)?;
        let n = self.dimension();
        if diagonal.len() != n || lower.len() != self.lower.len() || upper.len() != self.upper.len()
        {
            return Err(TridiagonalError::Dimension);
        }
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            self.diagonal[i] = finite(diagonal[i])?;
            if i + 1 < n {
                self.lower[i] = finite(lower[i])?;
                self.upper[i] = finite(upper[i])?;
                self.swapped[i] = false;
            }
            if i + 2 < n {
                self.upper_two[i] = 0.0;
            }
        }
        for i in 0..n.saturating_sub(1) {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let diagonal = self.diagonal[i];
            let below = self.lower[i];
            if diagonal.abs() >= below.abs() {
                if diagonal == 0.0 {
                    return Err(TridiagonalError::Singular { index: i });
                }
                let multiplier = finite(below / diagonal)?;
                self.lower[i] = multiplier;
                self.diagonal[i + 1] =
                    finite((-multiplier).mul_add(self.upper[i], self.diagonal[i + 1]))?;
            } else {
                // The old next diagonal becomes U[i,i+1]. Its old upper
                // entry becomes U[i,i+2] and must remain in the back solve.
                let multiplier = finite(diagonal / below)?;
                let next_diagonal = self.diagonal[i + 1];
                self.diagonal[i] = below;
                self.lower[i] = multiplier;
                self.diagonal[i + 1] = finite((-multiplier).mul_add(next_diagonal, self.upper[i]))?;
                self.upper[i] = next_diagonal;
                if i + 2 < n {
                    self.upper_two[i] = self.upper[i + 1];
                    self.upper[i + 1] = finite(-multiplier * self.upper[i + 1])?;
                }
                self.swapped[i] = true;
            }
        }
        if n != 0 && self.diagonal[n - 1] == 0.0 {
            return Err(TridiagonalError::Singular { index: n - 1 });
        }
        poll(&mut cancelled)?;
        self.factored = true;
        Ok(())
    }

    fn substitute(
        &mut self,
        b: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), TridiagonalError> {
        poll(cancelled)?;
        let n = self.dimension();
        if b.len() != n {
            return Err(TridiagonalError::Dimension);
        }
        if !self.factored {
            return Err(TridiagonalError::Unfactored);
        }
        for (i, &value) in b.iter().enumerate() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            self.rhs[i] = finite(value)?;
        }
        // Row interchanges and their eliminations must be applied together.
        // Applying all pivots first would incorrectly permute earlier L factors.
        for i in 0..n.saturating_sub(1) {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            if self.swapped[i] {
                let previous = self.rhs[i];
                self.rhs[i] = self.rhs[i + 1];
                self.rhs[i + 1] = finite((-self.lower[i]).mul_add(self.rhs[i], previous))?;
            } else {
                self.rhs[i + 1] = finite((-self.lower[i]).mul_add(self.rhs[i], self.rhs[i + 1]))?;
            }
        }
        for i in (0..n).rev() {
            if (n - 1 - i) % 128 == 0 {
                poll(cancelled)?;
            }
            let mut value = self.rhs[i];
            if i + 1 < n {
                value = (-self.upper[i]).mul_add(self.rhs[i + 1], value);
            }
            if i + 2 < n {
                value = (-self.upper_two[i]).mul_add(self.rhs[i + 2], value);
            }
            self.rhs[i] = finite(value / self.diagonal[i])?;
        }
        poll(cancelled)
    }

    /// Solve using retained factors and publish into the caller's output only
    /// after all admission, finite arithmetic and cancellation checks succeed.
    pub fn solve_factored_into(
        &mut self,
        b: &[f64],
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TridiagonalError> {
        if out.len() != self.dimension() {
            return Err(TridiagonalError::Dimension);
        }
        self.substitute(b, &mut cancelled)?;
        out.copy_from_slice(&self.rhs);
        Ok(())
    }

    /// Replace one RHS with its solution transactionally, reusing the factors.
    /// No allocation occurs, and refusal leaves the RHS bit-for-bit unchanged.
    pub fn solve_factored_in_place(
        &mut self,
        rhs: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TridiagonalError> {
        self.substitute(rhs, &mut cancelled)?;
        rhs.copy_from_slice(&self.rhs);
        Ok(())
    }

    /// Factor a fresh tridiagonal matrix and solve into caller-owned storage.
    /// All work uses the previously allocated linear workspace. Failure leaves
    /// `out` unchanged; factorization validity follows [`Self::factor`].
    #[allow(clippy::too_many_arguments)]
    pub fn solve_into(
        &mut self,
        lower: &[f64],
        diagonal: &[f64],
        upper: &[f64],
        b: &[f64],
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TridiagonalError> {
        if b.len() != self.dimension() || out.len() != self.dimension() {
            return Err(TridiagonalError::Dimension);
        }
        self.factor(lower, diagonal, upper, &mut cancelled)?;
        self.solve_factored_into(b, out, cancelled)
    }
}

#[cfg(test)]
#[path = "tridiagonal/tests.rs"]
mod tests;
