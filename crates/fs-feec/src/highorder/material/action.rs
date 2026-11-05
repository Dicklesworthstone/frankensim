//! Repeated material action without local dense or global sparse matrices.
use super::{MAX_LOCAL, MaterialOperator, normalized, pair, zeros};
use crate::highorder::simplex::{LocalFn, duffy_quadrature};
use crate::weighted::{WeightedError, poll};
use fs_qty::Dims;

/// Largest admitted Gauss-point count on each Duffy axis for a sampled field.
pub const MAX_FIELD_QUADRATURE: usize = 16;

struct SampledScalarField {
    values: Vec<f64>,
    scales: Vec<f64>,
}

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
    scalar_field: Option<SampledScalarField>,
    entry_dims: Dims,
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
        self.check_action_storage(self.space.r + 2, false, max_storage_bytes)?;
        self.prepare_cached_action(
            duffy_quadrature(self.space.r + 2),
            None,
            self.entry_dims,
            cancelled,
        )
    }

    /// Freeze a positive spatial scalar multiplier into the material action.
    ///
    /// The callback receives the original cell and its stored-order barycentric
    /// point exactly once per quadrature sample. For example, the cell weight
    /// can carry specific heat while this field supplies density, or a unit
    /// cell weight can be multiplied by `1/rho(x)` for a diffusion coefficient.
    /// Closed-body inertia uses original reference density on reference
    /// geometry; temperature changes must not independently create inventory.
    /// Query source errors propagate unchanged. Neither the callback nor its
    /// evolving source state is retained or resampled by warm applications.
    ///
    /// `quadrature_order` must be between `space.r + 2` and
    /// [`MAX_FIELD_QUADRATURE`]. The caller must choose sufficient order for
    /// the actual field and assess quadrature error. Positive samples establish
    /// the discrete weight; they do not certify an unsampled source domain or
    /// integrate an arbitrary constitutive law exactly. Source-history changes
    /// require a newly prepared action. `field_dims` multiplies the original
    /// material dimensions, rather than silently replacing their meaning.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_action_with_scalar_field<'p, E: From<WeightedError>>(
        &'p self,
        field_dims: Dims,
        quadrature_order: usize,
        max_storage_bytes: usize,
        field_at: &mut impl FnMut(usize, [f64; 4]) -> Result<f64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MaterialAction<'p, 'a, 'c>, E> {
        poll(cancelled)?;
        if quadrature_order < self.space.r + 2 || quadrature_order > MAX_FIELD_QUADRATURE {
            return Err(WeightedError::InvalidGeometry {
                cell: 0,
                reason: "sampled material quadrature order must be r+2..=16",
            }
            .into());
        }
        let entry_dims = self
            .entry_dims
            .checked_plus(field_dims)
            .ok_or(WeightedError::DimensionOverflow)?;
        self.check_action_storage(quadrature_order, true, max_storage_bytes)?;
        let quadrature = duffy_quadrature(quadrature_order);
        let count = self
            .weights
            .len()
            .checked_mul(quadrature.len())
            .ok_or(WeightedError::SizeOverflow)?;
        let mut field = SampledScalarField {
            values: Vec::new(),
            scales: Vec::new(),
        };
        field
            .values
            .try_reserve_exact(count)
            .map_err(|_| WeightedError::Allocation)?;
        field
            .scales
            .try_reserve_exact(self.weights.len())
            .map_err(|_| WeightedError::Allocation)?;
        for cell in 0..self.weights.len() {
            let start = field.values.len();
            let mut scale = 0.0_f64;
            for &(lambda, _) in &quadrature {
                poll(cancelled)?;
                let value = field_at(cell, lambda)?;
                if !value.is_finite() || value <= 0.0 {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "sampled scalar field must be finite and positive",
                    }
                    .into());
                }
                scale = scale.max(value);
                field.values.push(value);
            }
            for (i, value) in field.values[start..].iter_mut().enumerate() {
                if i % 256 == 0 {
                    poll(cancelled)?;
                }
                *value /= scale;
                if *value <= 0.0 {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "sampled scalar contrast is not representable",
                    }
                    .into());
                }
            }
            field.scales.push(scale);
        }
        Ok(self.prepare_cached_action(quadrature, Some(field), entry_dims, cancelled)?)
    }

    fn check_action_storage(
        &self,
        quadrature_order: usize,
        sampled: bool,
        max_storage_bytes: usize,
    ) -> Result<(), WeightedError> {
        let basis_count = self
            .weights
            .len()
            .checked_mul(self.local_count)
            .ok_or(WeightedError::SizeOverflow)?;
        let points = quadrature_order.pow(3);
        let samples = if sampled {
            self.weights
                .len()
                .checked_mul(points + 1)
                .ok_or(WeightedError::SizeOverflow)?
        } else {
            0
        };
        let bytes = basis_count
            .checked_mul(core::mem::size_of::<(usize, LocalFn)>())
            .and_then(|n| {
                n.checked_add(points.checked_mul(core::mem::size_of::<([f64; 4], f64)>())?)
            })
            .and_then(|n| n.checked_add(self.space.ndof.checked_mul(core::mem::size_of::<f64>())?))
            .and_then(|n| n.checked_add(samples.checked_mul(core::mem::size_of::<f64>())?))
            .ok_or(WeightedError::SizeOverflow)?;
        if bytes > max_storage_bytes {
            return Err(WeightedError::LimitExceeded {
                resource: "prepared material action bytes",
                requested: bytes,
                limit: max_storage_bytes,
            });
        }
        Ok(())
    }

    fn prepare_cached_action<'p>(
        &'p self,
        quadrature: Vec<([f64; 4], f64)>,
        scalar_field: Option<SampledScalarField>,
        entry_dims: Dims,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MaterialAction<'p, 'a, 'c>, WeightedError> {
        poll(cancelled)?;
        let basis_count = self
            .weights
            .len()
            .checked_mul(self.local_count)
            .ok_or(WeightedError::SizeOverflow)?;
        let mut dofs = Vec::new();
        dofs.try_reserve_exact(basis_count)
            .map_err(|_| WeightedError::Allocation)?;
        let mut diagonal = zeros(self.space.ndof)?;
        let mut values = [[0.0; 3]; MAX_LOCAL];
        for (cell, &weight) in self.weights.iter().enumerate() {
            poll(cancelled)?;
            let local_dofs = self.space.element_dofs(cell);
            let (scale, tensor, affine) = normalized(weight);
            let mut local_diagonal = [0.0; MAX_LOCAL];
            for (point, &(lambda, w)) in quadrature.iter().enumerate() {
                self.sample(cell, lambda, &local_dofs, &mut values, cancelled)?;
                let mut factor = w * (6.0 * self.geometry.vol_signed[cell].abs());
                if let Some(density) = affine {
                    factor *= density
                        .iter()
                        .zip(lambda)
                        .map(|(rho, l)| rho * l)
                        .sum::<f64>();
                }
                if let Some(field) = &scalar_field {
                    factor *= field.values[cell * quadrature.len() + point];
                }
                for i in 0..self.local_count {
                    local_diagonal[i] =
                        factor.mul_add(pair(tensor, values[i], values[i]), local_diagonal[i]);
                }
            }
            for (i, &(global, _)) in local_dofs.iter().enumerate() {
                let value = scaled_entry(
                    scale,
                    scalar_field.as_ref().map(|field| field.scales[cell]),
                    local_diagonal[i],
                );
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
            scalar_field,
            entry_dims,
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
        self.entry_dims
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
            for (point, &(lambda, w)) in self.quadrature.iter().enumerate() {
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
                if let Some(field) = &self.scalar_field {
                    factor *= field.values[cell * self.quadrature.len() + point];
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
                let value = scaled_entry(
                    scale,
                    self.scalar_field.as_ref().map(|field| field.scales[cell]),
                    local[i],
                );
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

/// Restore independent coefficient scales only after integration. Their raw
/// product may overflow when anticorrelated fields have a finite product.
/// Choose the multiplication order before rescaling a subnormal intermediate.
fn scaled_entry(base: f64, field: Option<f64>, value: f64) -> f64 {
    let Some(field) = field else {
        return base * value;
    };
    let (small, large) = (base.min(field), base.max(field));
    let (first, second) = if value.abs() <= 1.0 {
        (large, small)
    } else {
        (small, large)
    };
    let partial = value * first;
    if partial.is_subnormal() {
        // Exact powers of two keep the intermediate normal without changing
        // the represented product; there is no platform libm dependency.
        let up = f64::from_bits((1023_u64 + 512) << 52);
        let down = f64::from_bits((1023_u64 - 512) << 52);
        ((value * up) * first) * (second * down)
    } else {
        partial * second
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn factored_density_scaling_preserves_a_representable_subnormal() {
        let minimum = f64::from_bits(1);
        // 2 * 0.745 * 0.34 = 0.5066: the exact result rounds up to one
        // minimum subnormal. Rounding each unscaled product loses it to zero.
        assert_eq!(
            super::scaled_entry(0.745, Some(0.34), 2.0 * minimum),
            minimum
        );
        assert_eq!(
            super::scaled_entry(0.745, Some(0.34), -2.0 * minimum),
            -minimum
        );
    }
}
