//! Reuse the exact active friction matrix across independent driving directions.
//!
//! Source admission, composition reduction, fresh right-hand sides and original
//! physical residual checks remain with the constitutive owner. Only an exactly
//! equal assembled matrix and dimension can reuse these numerical factors.

use super::{MaxwellStefanErrorV1, filled, poll, positive, sum};

#[derive(Clone, Debug)]
pub(super) struct Factorization {
    source: Vec<f64>,
    factors: Vec<f64>,
    pivots: Vec<usize>,
    dimension: usize,
}

impl Factorization {
    pub(super) fn new(maximum_dimension: usize) -> Result<Self, MaxwellStefanErrorV1> {
        Ok(Self {
            source: filled(maximum_dimension * maximum_dimension, 0.0)?,
            factors: filled(maximum_dimension * maximum_dimension, 0.0)?,
            pivots: filled(maximum_dimension, 0)?,
            dimension: 0,
        })
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        (self.source.capacity() + self.factors.capacity()) * core::mem::size_of::<f64>()
            + self.pivots.capacity() * core::mem::size_of::<usize>()
    }

    pub(super) fn solve(
        &mut self,
        matrix: &[f64],
        right: &mut [f64],
        scales: &mut [f64],
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MaxwellStefanErrorV1> {
        debug_assert_eq!(matrix.len(), n * n);
        debug_assert_eq!(right.len(), n);
        debug_assert_eq!(scales.len(), n);
        debug_assert!(n > 0 && n <= self.pivots.len());
        poll(cancelled)?;
        let reusable = self.dimension == n
            && matrix
                .iter()
                .zip(&self.source)
                .all(|(actual, saved)| actual.to_bits() == saved.to_bits());
        if !reusable {
            // A cancelled or refused refactorization must never leave a
            // partially overwritten matrix marked as usable.
            self.dimension = 0;
            self.factors[..n * n].copy_from_slice(matrix);
            let a = &mut self.factors[..n * n];
            for i in 0..n {
                scales[i] = a[i * n..(i + 1) * n]
                    .iter()
                    .map(|value| value.abs())
                    .fold(0.0, f64::max);
                if !positive(scales[i]) {
                    return Err(MaxwellStefanErrorV1::SingularSystem);
                }
            }
            for k in 0..n {
                poll(cancelled)?;
                let mut pivot = k;
                for i in k + 1..n {
                    if a[i * n + k].abs() / scales[i] > a[pivot * n + k].abs() / scales[pivot] {
                        pivot = i;
                    }
                }
                if a[pivot * n + k] == 0.0 || !a[pivot * n + k].is_finite() {
                    return Err(MaxwellStefanErrorV1::SingularSystem);
                }
                self.pivots[k] = pivot;
                if pivot != k {
                    // Earlier columns record the row operations as they
                    // occurred, not a conventionally permuted lower factor.
                    // Keeping them in place permits exact chronological RHS
                    // replay, including every later row interchange.
                    for j in k..n {
                        a.swap(k * n + j, pivot * n + j);
                    }
                    scales.swap(k, pivot);
                }
                for i in k + 1..n {
                    let factor = a[i * n + k] / a[k * n + k];
                    a[i * n + k] = factor;
                    for j in k + 1..n {
                        a[i * n + j] -= factor * a[k * n + j];
                    }
                }
            }
            self.source[..n * n].copy_from_slice(matrix);
            self.dimension = n;
        }

        // Apply the same swaps and subtractions, in the same order, that the
        // original scalar elimination applied to its driving vector. Factor
        // storage remains immutable even when this solve is cancelled.
        let a = &self.factors;
        for k in 0..n {
            poll(cancelled)?;
            right.swap(k, self.pivots[k]);
            for i in k + 1..n {
                right[i] -= a[i * n + k] * right[k];
            }
        }
        for i in (0..n).rev() {
            poll(cancelled)?;
            let rhs = right[i] - sum((i + 1..n).map(|j| a[i * n + j] * right[j]));
            right[i] = rhs / a[i * n + i];
            if !right[i].is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic("diffusion linear solve"));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
