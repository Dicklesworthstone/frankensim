//! Work-conjugate arithmetic on an original physical coordinate subspace.
//!
//! The factorization/domain owner supplies N and its independent constraint
//! rows. This module admits the actually represented physical lift, its scaled
//! column independence and A*N residual, then shares q=N*z, f_z=N^T*f and
//! N^T*L*N. It does not select rank, invent a nullspace, alter constraints or
//! replace an existing QR owner. Completeness of a nullspace remains with the
//! owner that admits independent rows and selects n-rank(A) columns.

use fs_math::dd::Dd;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SubspaceControls {
    /// Both dimensions must be at most this explicit ceiling, in 1..=1024.
    pub maximum_dimension: usize,
    /// Maximum |(A*S)*(S^-1*N)_j| / ||(A*S)_row||_1.
    /// Nonzero rows are independently normalized, so row units/scales cancel.
    pub maximum_constraint_residual: f64,
    /// Maximum absolute entry of (S^-1*N)^T*(S^-1*N)-I. Its product with the
    /// reduced dimension must be below one, retaining column independence.
    pub maximum_orthogonality_error: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SubspaceError {
    Invalid(&'static str),
    Constraint {
        row: usize,
        column: usize,
        residual: f64,
        maximum: f64,
    },
    Orthogonality {
        error: f64,
        maximum: f64,
    },
    Cancelled,
}
impl core::fmt::Display for SubspaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical subspace: {self:?}")
    }
}
impl std::error::Error for SubspaceError {}

/// Immutable row-major physical lift N. Reduced coordinates are dimensionless
/// when S carries the physical coordinate units. No massless reaction unknowns
/// enter any reduced operator.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearSubspace {
    physical: usize,
    reduced: usize,
    basis: Vec<f64>,
    maximum_constraint_residual: f64,
    maximum_orthogonality_error: f64,
}

impl LinearSubspace {
    /// Admit an already constructed physical basis against original finite
    /// homogeneous rows A and positive coordinate scales S. Empty reduced
    /// spaces, no rows and the completely empty physical space are supported.
    /// Rank decisions and affine-position admission belong to the caller.
    pub fn try_new(
        coordinate_scales: &[f64],
        reduced: usize,
        basis: Vec<f64>,
        rows: &[f64],
        controls: SubspaceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SubspaceError> {
        poll(&mut cancelled)?;
        let n = coordinate_scales.len();
        if !(1..=1024).contains(&controls.maximum_dimension)
            || n > controls.maximum_dimension
            || reduced > n
            || basis.len() != n * reduced
            || (n == 0 && !rows.is_empty())
            || (n != 0 && rows.len() % n != 0)
            || (n != 0 && rows.len() / n > controls.maximum_dimension)
            || coordinate_scales
                .iter()
                .any(|s| !s.is_finite() || *s <= 0.0)
            || basis.iter().chain(rows).any(|v| !v.is_finite())
            || [
                controls.maximum_constraint_residual,
                controls.maximum_orthogonality_error,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0 || *x >= 1.0)
            || reduced as f64 * controls.maximum_orthogonality_error >= 1.0
        {
            return Err(SubspaceError::Invalid(
                "bounded finite dimensions, positive scales and tolerances required",
            ));
        }
        let mut scaled = Vec::with_capacity(basis.len());
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..reduced {
                let value = basis[i * reduced + j] / coordinate_scales[i];
                if !value.is_finite() || (basis[i * reduced + j] != 0.0 && value == 0.0) {
                    return Err(SubspaceError::Invalid(
                        "scaled physical basis is not representable",
                    ));
                }
                scaled.push(value);
            }
        }
        let mut maximum_orthogonality_error = 0.0_f64;
        for i in 0..reduced {
            poll(&mut cancelled)?;
            for j in 0..=i {
                let mut sum = if i == j { -Dd::ONE } else { Dd::ZERO };
                for row in 0..n {
                    sum = sum + product(scaled[row * reduced + i], scaled[row * reduced + j])?;
                }
                let error = represented(sum)?.abs();
                maximum_orthogonality_error = maximum_orthogonality_error.max(error);
            }
        }
        if maximum_orthogonality_error > controls.maximum_orthogonality_error {
            return Err(SubspaceError::Orthogonality {
                error: maximum_orthogonality_error,
                maximum: controls.maximum_orthogonality_error,
            });
        }
        let mut maximum_constraint_residual = 0.0_f64;
        if n != 0 {
            for (row, original) in rows.chunks_exact(n).enumerate() {
                poll(&mut cancelled)?;
                let mut coefficients = Vec::with_capacity(n);
                let mut scale = 0.0_f64;
                for (&a, &s) in original.iter().zip(coordinate_scales) {
                    let value = represented(product(a, s)?)?;
                    scale = scale.max(value.abs());
                    coefficients.push(value);
                }
                if scale == 0.0 {
                    return Err(SubspaceError::Invalid(
                        "nonzero original constraint row required",
                    ));
                }
                let mut norm = 0.0;
                for value in &mut coefficients {
                    let next = *value / scale;
                    if *value != 0.0 && next == 0.0 {
                        return Err(SubspaceError::Invalid(
                            "normalized original row is not representable",
                        ));
                    }
                    *value = next;
                    norm += next.abs();
                }
                for column in 0..reduced {
                    let mut sum = Dd::ZERO;
                    for i in 0..n {
                        // Evaluate the ORIGINAL represented A and N, rather
                        // than reconstructing their product from rounded A*S
                        // and S^-1*N. Scaling precedes accumulation only.
                        sum = sum
                            + product(original[i], basis[i * reduced + column])?
                                / Dd::from_f64(scale);
                    }
                    let residual = represented(sum)?.abs() / norm;
                    if !residual.is_finite() || residual > controls.maximum_constraint_residual {
                        return Err(SubspaceError::Constraint {
                            row,
                            column,
                            residual,
                            maximum: controls.maximum_constraint_residual,
                        });
                    }
                    maximum_constraint_residual = maximum_constraint_residual.max(residual);
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            physical: n,
            reduced,
            basis,
            maximum_constraint_residual,
            maximum_orthogonality_error,
        })
    }

    pub const fn physical_dimension(&self) -> usize {
        self.physical
    }
    pub const fn dimension(&self) -> usize {
        self.reduced
    }
    pub fn basis(&self) -> &[f64] {
        &self.basis
    }
    pub const fn maximum_constraint_residual(&self) -> f64 {
        self.maximum_constraint_residual
    }
    pub const fn maximum_orthogonality_error(&self) -> f64 {
        self.maximum_orthogonality_error
    }

    /// q=N*z, in original physical coordinate order.
    pub fn lift(
        &self,
        reduced: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, SubspaceError> {
        poll(&mut cancelled)?;
        admit_vector(reduced, self.reduced)?;
        let mut physical = vec![0.0; self.physical];
        for (i, value) in physical.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            let mut sum = Dd::ZERO;
            for (j, &z) in reduced.iter().enumerate() {
                sum = sum + product(self.basis[i * self.reduced + j], z)?;
            }
            *value = represented(sum)?;
        }
        poll(&mut cancelled)?;
        Ok(physical)
    }

    /// f_z=N^T*f. This is the work-conjugate of `lift`, not coordinate truncation.
    pub fn project(
        &self,
        physical: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, SubspaceError> {
        poll(&mut cancelled)?;
        admit_vector(physical, self.physical)?;
        let mut reduced = vec![0.0; self.reduced];
        for (j, value) in reduced.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            let mut sum = Dd::ZERO;
            for (i, &f) in physical.iter().enumerate() {
                sum = sum + product(self.basis[i * self.reduced + j], f)?;
            }
            *value = represented(sum)?;
        }
        poll(&mut cancelled)?;
        Ok(reduced)
    }

    /// Full N^T*L*N of a row-major original operator. Signed and nonsymmetric
    /// entries are preserved; no modal truncation, averaging or repair occurs.
    pub fn congruence(
        &self,
        matrix: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, SubspaceError> {
        self.congruence_impl(matrix, false, &mut cancelled)
    }

    /// The same congruence for an exactly symmetric original operator such as
    /// native mass. Evaluate each bilinear entry once and mirror its identical
    /// mathematical counterpart. A nonsymmetric input is refused, never averaged.
    pub fn symmetric_congruence(
        &self,
        matrix: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, SubspaceError> {
        self.congruence_impl(matrix, true, &mut cancelled)
    }

    fn congruence_impl(
        &self,
        matrix: &[f64],
        symmetric: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, SubspaceError> {
        poll(cancelled)?;
        let (n, r) = (self.physical, self.reduced);
        admit_vector(matrix, n * n)?;
        if symmetric && (0..n).any(|i| (0..i).any(|j| matrix[i * n + j] != matrix[j * n + i])) {
            return Err(SubspaceError::Invalid(
                "symmetric congruence requires original exact symmetry",
            ));
        }
        let mut action = vec![Dd::ZERO; n * r];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..r {
                for k in 0..n {
                    action[i * r + j] =
                        action[i * r + j] + product(matrix[i * n + k], self.basis[k * r + j])?;
                }
                represented(action[i * r + j])?;
            }
        }
        let mut result = vec![0.0; r * r];
        for i in 0..r {
            poll(cancelled)?;
            for j in (if symmetric { i } else { 0 })..r {
                let mut sum = Dd::ZERO;
                for k in 0..n {
                    let term = Dd::from_f64(self.basis[k * r + i]) * action[k * r + j];
                    if self.basis[k * r + i] != 0.0
                        && action[k * r + j] != Dd::ZERO
                        && term == Dd::ZERO
                    {
                        return Err(SubspaceError::Invalid(
                            "physical subspace product is not representable",
                        ));
                    }
                    sum = sum + term;
                }
                let value = represented(sum)?;
                result[i * r + j] = value;
                if symmetric {
                    result[j * r + i] = value;
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SubspaceError> {
    if cancelled() {
        Err(SubspaceError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_vector(values: &[f64], n: usize) -> Result<(), SubspaceError> {
    if values.len() != n || values.iter().any(|v| !v.is_finite()) {
        Err(SubspaceError::Invalid(
            "complete finite original vector required",
        ))
    } else {
        Ok(())
    }
}
fn represented(value: Dd) -> Result<f64, SubspaceError> {
    let rounded = value.to_f64();
    if !rounded.is_finite() || (rounded == 0.0 && value.lo != 0.0) {
        Err(SubspaceError::Invalid(
            "physical subspace arithmetic is not representable",
        ))
    } else {
        Ok(rounded)
    }
}
fn product(a: f64, b: f64) -> Result<Dd, SubspaceError> {
    let value = Dd::from_f64(a) * Dd::from_f64(b);
    if !value.hi.is_finite()
        || !value.lo.is_finite()
        || (a != 0.0 && b != 0.0 && value.hi == 0.0 && value.lo == 0.0)
    {
        Err(SubspaceError::Invalid(
            "physical subspace product is not representable",
        ))
    } else {
        Ok(value)
    }
}

#[cfg(test)]
#[path = "subspace/tests.rs"]
mod tests;
