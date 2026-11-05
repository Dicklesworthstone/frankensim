//! Repeated material action without local dense or global sparse matrices.
use super::{MAX_LOCAL, MaterialOperator, normalized, pair, zeros};
use crate::highorder::simplex::{LocalFn, duffy_quadrature};
use crate::weighted::{WeightedError, poll};
use fs_qty::Dims;

/// Prepared quadrature action of the original material-weighted H1 pairing.
///
/// The immutable plan caches oriented local basis descriptors, one shared
/// quadrature rule and the exact quadrature diagonal. Storage is linear in
/// the number of local degrees of freedom, rather than their square. Warm
/// application uses caller-owned output and scratch, with no allocation.
/// Coefficients and geometry retain the original operator's lifetime and units.
pub struct MaterialAction<'p, 'a, 'c> {
    operator: &'p MaterialOperator<'a, 'c>,
    dofs: Vec<(usize, LocalFn)>,
    quadrature: Vec<([f64; 4], f64)>,
    diagonal: Vec<f64>,
}

impl<'a, 'c> MaterialOperator<'a, 'c> {
    /// Prepare repeated density/material action and its Jacobi diagonal.
    ///
    /// `max_storage_bytes` bounds retained basis descriptors, quadrature and
    /// diagonal entries. Caller-owned input/output/scratch are not part of that
    /// bound. The assembly triplet cap does not constrain this matrix-free path.
    /// Existing order/cell/dof admission still applies. The original bounded
    /// basis-list helper makes temporary allocations of at most one element.
    pub fn prepare_action<'p>(
        &'p self,
        max_storage_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MaterialAction<'p, 'a, 'c>, WeightedError> {
        poll(cancelled)?;
        let basis_count = self
            .weights
            .len()
            .checked_mul(self.local_count)
            .ok_or(WeightedError::SizeOverflow)?;
        let points = (self.space.r + 2).pow(3);
        let bytes = basis_count
            .checked_mul(core::mem::size_of::<(usize, LocalFn)>())
            .and_then(|n| {
                n.checked_add(points.checked_mul(core::mem::size_of::<([f64; 4], f64)>())?)
            })
            .and_then(|n| n.checked_add(self.space.ndof.checked_mul(core::mem::size_of::<f64>())?))
            .ok_or(WeightedError::SizeOverflow)?;
        if bytes > max_storage_bytes {
            return Err(WeightedError::LimitExceeded {
                resource: "prepared material action bytes",
                requested: bytes,
                limit: max_storage_bytes,
            });
        }
        let mut dofs = Vec::new();
        dofs.try_reserve_exact(basis_count)
            .map_err(|_| WeightedError::Allocation)?;
        let quadrature = duffy_quadrature(self.space.r + 2);
        let mut diagonal = zeros(self.space.ndof)?;
        let mut values = [[0.0; 3]; MAX_LOCAL];
        for (cell, &weight) in self.weights.iter().enumerate() {
            poll(cancelled)?;
            let local_dofs = self.space.element_dofs(cell);
            let (scale, tensor, affine) = normalized(weight);
            let mut local_diagonal = [0.0; MAX_LOCAL];
            for &(lambda, w) in &quadrature {
                self.sample(cell, lambda, &local_dofs, &mut values, cancelled)?;
                let mut factor = w * (6.0 * self.geometry.vol_signed[cell].abs());
                if let Some(density) = affine {
                    factor *= density
                        .iter()
                        .zip(lambda)
                        .map(|(rho, l)| rho * l)
                        .sum::<f64>();
                }
                for i in 0..self.local_count {
                    local_diagonal[i] =
                        factor.mul_add(pair(tensor, values[i], values[i]), local_diagonal[i]);
                }
            }
            for (i, &(global, _)) in local_dofs.iter().enumerate() {
                let value = scale * local_diagonal[i];
                if !value.is_finite() || !(diagonal[global] + value).is_finite() {
                    return Err(WeightedError::NonFiniteEntry);
                }
                if value <= 0.0 {
                    return Err(WeightedError::InvalidGeometry {
                        cell,
                        reason: "positive material diagonal is not representable",
                    });
                }
                diagonal[global] += value;
            }
            dofs.extend(local_dofs);
        }
        poll(cancelled)?;
        Ok(MaterialAction {
            operator: self,
            dofs,
            quadrature,
            diagonal,
        })
    }
}

impl MaterialAction<'_, '_, '_> {
    /// Original global H1 coordinate count, including any isolated vertices.
    #[must_use]
    pub fn ndof(&self) -> usize {
        self.operator.space.ndof
    }

    /// Dimensions added by this operator to the primal field dimensions.
    #[must_use]
    pub const fn entry_dims(&self) -> Dims {
        self.operator.entry_dims
    }

    /// Dimensions of the dual result, retaining the original material units.
    pub fn output_dims(&self, input: Dims) -> Result<Dims, WeightedError> {
        self.entry_dims()
            .checked_plus(input)
            .ok_or(WeightedError::DimensionOverflow)
    }

    /// Diagonal of the same density/material quadrature operator.
    ///
    /// This can be used for Jacobi scaling after the caller's boundary/gauge
    /// treatment. Isolated degrees of freedom remain zero; no inverse or
    /// coercivity claim is made for an unconstrained diffusion problem.
    #[must_use]
    pub fn diagonal(&self) -> &[f64] {
        &self.diagonal
    }

    /// Compute `output = A * input` without allocating or forming a matrix.
    ///
    /// All three slices must have `ndof()` entries. Scratch may be overwritten
    /// on refusal, but output is unchanged unless the entire action succeeds.
    /// Each quadrature point contracts the trial field once, then projects it
    /// onto the test functions: work is linear, not quadratic, in local basis
    /// size. This is the original consistent operator, not density lumping.
    pub fn apply_into(
        &self,
        input: &[f64],
        output: &mut [f64],
        scratch: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WeightedError> {
        poll(cancelled)?;
        let op = self.operator;
        if input.len() != self.ndof() || output.len() != self.ndof() || scratch.len() != self.ndof()
        {
            return Err(WeightedError::InvalidGeometry {
                cell: 0,
                reason: "material action vectors must match the original H1 layout",
            });
        }
        for (i, (&value, entry)) in input.iter().zip(scratch.iter_mut()).enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            if !value.is_finite() {
                return Err(WeightedError::NonFiniteEntry);
            }
            *entry = 0.0;
        }
        let mut values = [[0.0; 3]; MAX_LOCAL];
        for (cell, &weight) in op.weights.iter().enumerate() {
            poll(cancelled)?;
            let dofs = &self.dofs[cell * op.local_count..(cell + 1) * op.local_count];
            let (scale, tensor, affine) = normalized(weight);
            let mut local = [0.0; MAX_LOCAL];
            for &(lambda, w) in &self.quadrature {
                op.sample(cell, lambda, dofs, &mut values, cancelled)?;
                let mut trial = [0.0; 3];
                for (i, &(global, _)) in dofs.iter().enumerate() {
                    for c in 0..3 {
                        trial[c] = input[global].mul_add(values[i][c], trial[c]);
                    }
                }
                if trial.iter().any(|v| !v.is_finite()) {
                    return Err(WeightedError::NonFiniteEntry);
                }
                let mut flux = [0.0; 3];
                for c in 0..3 {
                    flux[c] = tensor[c][0].mul_add(
                        trial[0],
                        tensor[c][1].mul_add(trial[1], tensor[c][2] * trial[2]),
                    );
                }
                let mut factor = w * (6.0 * op.geometry.vol_signed[cell].abs());
                if let Some(density) = affine {
                    factor *= density
                        .iter()
                        .zip(lambda)
                        .map(|(rho, l)| rho * l)
                        .sum::<f64>();
                }
                for i in 0..op.local_count {
                    let value = values[i][0].mul_add(
                        flux[0],
                        values[i][1].mul_add(flux[1], values[i][2] * flux[2]),
                    );
                    local[i] = factor.mul_add(value, local[i]);
                }
            }
            for (i, &(global, _)) in dofs.iter().enumerate() {
                let value = scale * local[i];
                let sum = scratch[global] + value;
                if !value.is_finite() || !sum.is_finite() {
                    return Err(WeightedError::NonFiniteEntry);
                }
                scratch[global] = sum;
            }
        }
        poll(cancelled)?;
        output.copy_from_slice(scratch);
        Ok(())
    }
}
