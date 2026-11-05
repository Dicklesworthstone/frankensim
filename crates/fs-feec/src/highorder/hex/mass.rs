//! Spatially weighted consistent mass by tensor partial assembly.

use super::TensorSpace;
use crate::highorder::quad1d::{gauss_legendre, lobatto_shapes};
use crate::weighted::{WeightedError, poll};
use fs_qty::Dims;

/// Maximum number of Gauss points per axis in a prepared tensor mass.
pub const MAX_TENSOR_MASS_QUADRATURE: usize = 16;

/// Consistent Q_r mass with a frozen positive spatial material field.
///
/// The unit-cube coordinates are in metres. With a density coefficient in
/// kg/m³, this pairs acceleration with force and velocity with momentum.
/// Hierarchical Lobatto coefficients are retained: bubble coefficients are
/// not nodal samples, and the mass is not diagonal or row-sum lumped.
///
/// Six rectangular axis contractions apply Bᵀ diag(w rho) B without forming
/// element matrices. For q = O(r), work is O(r⁴) per cell and retained field
/// storage is O(q³) per cell. Warm application allocates nothing.
pub struct TensorMass<'a> {
    space: &'a TensorSpace,
    q: usize,
    basis: Vec<f64>,
    weights: Vec<f64>,
    scales: Vec<f64>,
    diagonal: Vec<f64>,
    entry_dims: Dims,
    workspace_len: usize,
}

impl TensorSpace {
    /// Freeze a spatial material coefficient into the consistent tensor mass.
    ///
    /// `coefficient_at(cell, point)` receives the original `[cx, cy, cz]` cell
    /// and a physical point in the unit cube exactly once per tensor Gauss
    /// sample. Cell-local discontinuities are allowed. Source errors propagate
    /// unchanged. Every sample must be finite and strictly positive. Neither
    /// the callback nor its mutable source is retained or queried during apply.
    ///
    /// For closed-body inertia, supply reference density on reference geometry;
    /// a temperature update alone must not replace the body's mass inventory.
    /// `coefficient_dims` may instead describe volumetric heat capacity or
    /// another positive scalar pairing. Matrix dimensions add metre³.
    ///
    /// `quadrature_order` is in `r + 1 ..= 16`. A polynomial coefficient of
    /// coordinate-wise degree k is integrated exactly when `2q-1 >= 2r+k`.
    /// General fields require the caller to assess quadrature error; positive
    /// samples do not certify the unsampled source domain. A changed material
    /// state requires a new plan. `max_storage_bytes` bounds retained arrays,
    /// excluding caller workspace and bounded per-cell construction scratch.
    pub fn prepare_weighted_mass<E: From<WeightedError>>(
        &self,
        coefficient_dims: Dims,
        quadrature_order: usize,
        max_storage_bytes: usize,
        coefficient_at: &mut impl FnMut([usize; 3], [f64; 3]) -> Result<f64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TensorMass<'_>, E> {
        poll(cancelled)?;
        let q = quadrature_order;
        if q < self.p || q > MAX_TENSOR_MASS_QUADRATURE {
            return Err(WeightedError::InvalidGeometry {
                cell: 0,
                reason: "tensor mass quadrature order must be r+1..=16",
            }
            .into());
        }
        let entry_dims = coefficient_dims
            .checked_plus(Dims([3, 0, 0, 0, 0, 0]))
            .ok_or(WeightedError::DimensionOverflow)?;
        let cells = self
            .m
            .checked_mul(self.m)
            .and_then(|n| n.checked_mul(self.m))
            .ok_or(WeightedError::SizeOverflow)?;
        let points = q * q * q; // q <= 16, admitted above.
        let samples = cells
            .checked_mul(points)
            .ok_or(WeightedError::SizeOverflow)?;
        let entries = samples
            .checked_add(cells)
            .and_then(|n| n.checked_add(q * self.p))
            .and_then(|n| n.checked_add(self.ndof))
            .ok_or(WeightedError::SizeOverflow)?;
        let bytes = entries
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or(WeightedError::SizeOverflow)?;
        if bytes > max_storage_bytes {
            return Err(WeightedError::LimitExceeded {
                resource: "prepared tensor mass bytes",
                requested: bytes,
                limit: max_storage_bytes,
            }
            .into());
        }
        let workspace_len = self
            .ndof
            .checked_add(2 * points)
            .ok_or(WeightedError::SizeOverflow)?;
        workspace_len
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or(WeightedError::SizeOverflow)?;
        let (nodes, gauss_weights) = gauss_legendre(q);
        let mut mass = TensorMass {
            space: self,
            q,
            basis: zeros(q * self.p)?,
            weights: zeros(samples)?,
            scales: zeros(cells)?,
            diagonal: zeros(self.ndof)?,
            entry_dims,
            workspace_len,
        };
        for (i, &x) in nodes.iter().enumerate() {
            mass.basis[i * self.p..(i + 1) * self.p].copy_from_slice(&lobatto_shapes(self.r, x).0);
        }
        let h = 1.0 / self.m as f64;
        let jac = (h / 2.0) * (h / 2.0) * (h / 2.0);
        let mut first = zeros(points)?;
        let mut second = zeros(points)?;
        for cx in 0..self.m {
            for cy in 0..self.m {
                for cz in 0..self.m {
                    let cell = (cx * self.m + cy) * self.m + cz;
                    let start = cell * points;
                    let mut scale = 0.0_f64;
                    for i in 0..q {
                        for j in 0..q {
                            for k in 0..q {
                                poll(cancelled)?;
                                let position = [cx, cy, cz].map(|c| c as f64 * h);
                                let point = [
                                    (nodes[i] + 1.0).mul_add(h / 2.0, position[0]),
                                    (nodes[j] + 1.0).mul_add(h / 2.0, position[1]),
                                    (nodes[k] + 1.0).mul_add(h / 2.0, position[2]),
                                ];
                                let value = coefficient_at([cx, cy, cz], point)?;
                                if !value.is_finite() || value <= 0.0 {
                                    return Err(WeightedError::InvalidCoefficient {
                                        cell,
                                        reason: "tensor mass samples must be finite and positive",
                                    }
                                    .into());
                                }
                                scale = scale.max(value);
                                mass.weights[start + (i * q + j) * q + k] = value;
                            }
                        }
                    }
                    mass.scales[cell] = scale;
                    for i in 0..q {
                        poll(cancelled)?;
                        for j in 0..q {
                            for k in 0..q {
                                let value = &mut mass.weights[start + (i * q + j) * q + k];
                                *value = (*value / scale)
                                    * (gauss_weights[i]
                                        * gauss_weights[j]
                                        * gauss_weights[k]
                                        * jac);
                                if !value.is_finite() || *value <= 0.0 {
                                    return Err(WeightedError::InvalidCoefficient {
                                        cell,
                                        reason: "tensor mass sample contrast is not representable",
                                    }
                                    .into());
                                }
                            }
                        }
                    }
                    // The same quadrature diagonal, via squared basis tables.
                    // This also uses O(r^4), rather than a basis/point double loop.
                    mass.project_z::<true>(&mass.weights[start..start + points], &mut first);
                    mass.project_y::<true>(&first, &mut second);
                    mass.project_x::<true>(&second, &mut first);
                    for x in 0..self.p {
                        let gx = self.lat1(cx, x);
                        for y in 0..self.p {
                            let gy = self.lat1(cy, y);
                            for z in 0..self.p {
                                let global = self.gid(gx, gy, self.lat1(cz, z));
                                let value = scale * first[(x * self.p + y) * self.p + z];
                                let sum = mass.diagonal[global] + value;
                                if !value.is_finite() || !sum.is_finite() {
                                    return Err(WeightedError::NonFiniteEntry.into());
                                }
                                if value <= 0.0 {
                                    return Err(WeightedError::InvalidCoefficient {
                                        cell,
                                        reason: "positive tensor mass diagonal is not representable",
                                    }
                                    .into());
                                }
                                mass.diagonal[global] = sum;
                            }
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        Ok(mass)
    }
}

impl TensorMass<'_> {
    /// Global hierarchical H1 coordinate count.
    #[must_use]
    pub const fn ndof(&self) -> usize {
        self.space.ndof
    }

    /// Gauss samples per axis; there are q³ samples per cell.
    #[must_use]
    pub const fn quadrature_order(&self) -> usize {
        self.q
    }

    /// Number of f64 entries needed for atomic, allocation-free application.
    /// This is `ndof + 2*q³`, independent of the number of cells in the grid.
    #[must_use]
    pub const fn workspace_len(&self) -> usize {
        self.workspace_len
    }

    /// Physical dimensions added to the input field by this mass operator.
    #[must_use]
    pub const fn entry_dims(&self) -> Dims {
        self.entry_dims
    }

    /// Physical dimensions of the resulting dual field.
    pub fn output_dims(&self, input: Dims) -> Result<Dims, WeightedError> {
        self.entry_dims
            .checked_plus(input)
            .ok_or(WeightedError::DimensionOverflow)
    }

    /// Consistent quadrature diagonal, suitable for Jacobi scaling.
    #[must_use]
    pub fn diagonal(&self) -> &[f64] {
        &self.diagonal
    }

    /// Apply the frozen consistent mass without allocation or source queries.
    ///
    /// Input and output have `ndof()` entries; workspace has `workspace_len()`
    /// entries. Workspace may change on error, but output changes only after
    /// every cell and the final cancellation checkpoint succeeds. Cell visits
    /// and contraction accumulation order are fixed for reproducible reuse.
    pub fn apply_into(
        &self,
        input: &[f64],
        output: &mut [f64],
        workspace: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WeightedError> {
        crate::highorder::fma::apply_tensor_mass_dispatch(self, input, output, workspace, cancelled)
    }

    // Inline into the existing FMA capsule, including the rectangular kernels,
    // so supported x86 CPUs use fused instructions instead of libm calls.
    #[allow(clippy::inline_always)]
    #[inline(always)]
    pub(crate) fn apply_body(
        &self,
        input: &[f64],
        output: &mut [f64],
        workspace: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WeightedError> {
        poll(cancelled)?;
        if input.len() != self.ndof()
            || output.len() != self.ndof()
            || workspace.len() != self.workspace_len
        {
            return Err(WeightedError::InvalidGeometry {
                cell: 0,
                reason: "tensor mass vectors must match the original H1 layout and workspace",
            });
        }
        let (result, local) = workspace.split_at_mut(self.ndof());
        let (first, second) = local.split_at_mut(self.q * self.q * self.q);
        for (i, (&value, entry)) in input.iter().zip(result.iter_mut()).enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            if !value.is_finite() {
                return Err(WeightedError::NonFiniteEntry);
            }
            *entry = 0.0;
        }
        let space = self.space;
        let p = space.p;
        let points = self.q * self.q * self.q;
        for cx in 0..space.m {
            for cy in 0..space.m {
                for cz in 0..space.m {
                    poll(cancelled)?;
                    let cell = (cx * space.m + cy) * space.m + cz;
                    for x in 0..p {
                        let gx = space.lat1(cx, x);
                        for y in 0..p {
                            let gy = space.lat1(cy, y);
                            for z in 0..p {
                                first[(x * p + y) * p + z] =
                                    input[space.gid(gx, gy, space.lat1(cz, z))];
                            }
                        }
                    }
                    self.evaluate_x(first, second);
                    self.evaluate_y(second, first);
                    self.evaluate_z(first, second);
                    for (value, &weight) in second
                        .iter_mut()
                        .zip(&self.weights[cell * points..(cell + 1) * points])
                    {
                        *value *= weight;
                    }
                    self.project_z::<false>(second, first);
                    self.project_y::<false>(first, second);
                    self.project_x::<false>(second, first);
                    for x in 0..p {
                        let gx = space.lat1(cx, x);
                        for y in 0..p {
                            let gy = space.lat1(cy, y);
                            for z in 0..p {
                                let global = space.gid(gx, gy, space.lat1(cz, z));
                                let value = self.scales[cell] * first[(x * p + y) * p + z];
                                let sum = result[global] + value;
                                if !value.is_finite() || !sum.is_finite() {
                                    return Err(WeightedError::NonFiniteEntry);
                                }
                                result[global] = sum;
                            }
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        output.copy_from_slice(result);
        Ok(())
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn evaluate_x(&self, input: &[f64], output: &mut [f64]) {
        let (p, q) = (self.space.p, self.q);
        for i in 0..q {
            for y in 0..p {
                for z in 0..p {
                    let mut sum = 0.0_f64;
                    for x in 0..p {
                        sum = self.basis[i * p + x].mul_add(input[(x * p + y) * p + z], sum);
                    }
                    output[(i * p + y) * p + z] = sum;
                }
            }
        }
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn evaluate_y(&self, input: &[f64], output: &mut [f64]) {
        let (p, q) = (self.space.p, self.q);
        for i in 0..q {
            for j in 0..q {
                for z in 0..p {
                    let mut sum = 0.0_f64;
                    for y in 0..p {
                        sum = self.basis[j * p + y].mul_add(input[(i * p + y) * p + z], sum);
                    }
                    output[(i * q + j) * p + z] = sum;
                }
            }
        }
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn evaluate_z(&self, input: &[f64], output: &mut [f64]) {
        let (p, q) = (self.space.p, self.q);
        for i in 0..q {
            for j in 0..q {
                for k in 0..q {
                    let mut sum = 0.0_f64;
                    for z in 0..p {
                        sum = self.basis[k * p + z].mul_add(input[(i * q + j) * p + z], sum);
                    }
                    output[(i * q + j) * q + k] = sum;
                }
            }
        }
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn project_z<const SQUARED: bool>(&self, input: &[f64], output: &mut [f64]) {
        let (p, q) = (self.space.p, self.q);
        for i in 0..q {
            for j in 0..q {
                for z in 0..p {
                    let mut sum = 0.0_f64;
                    for k in 0..q {
                        let b = factor::<SQUARED>(self.basis[k * p + z]);
                        sum = b.mul_add(input[(i * q + j) * q + k], sum);
                    }
                    output[(i * q + j) * p + z] = sum;
                }
            }
        }
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn project_y<const SQUARED: bool>(&self, input: &[f64], output: &mut [f64]) {
        let (p, q) = (self.space.p, self.q);
        for i in 0..q {
            for y in 0..p {
                for z in 0..p {
                    let mut sum = 0.0_f64;
                    for j in 0..q {
                        let b = factor::<SQUARED>(self.basis[j * p + y]);
                        sum = b.mul_add(input[(i * q + j) * p + z], sum);
                    }
                    output[(i * p + y) * p + z] = sum;
                }
            }
        }
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn project_x<const SQUARED: bool>(&self, input: &[f64], output: &mut [f64]) {
        let (p, q) = (self.space.p, self.q);
        for x in 0..p {
            for y in 0..p {
                for z in 0..p {
                    let mut sum = 0.0_f64;
                    for i in 0..q {
                        let b = factor::<SQUARED>(self.basis[i * p + x]);
                        sum = b.mul_add(input[(i * p + y) * p + z], sum);
                    }
                    output[(x * p + y) * p + z] = sum;
                }
            }
        }
    }
}

#[inline]
fn factor<const SQUARED: bool>(value: f64) -> f64 {
    if SQUARED { value * value } else { value }
}

fn zeros(len: usize) -> Result<Vec<f64>, WeightedError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| WeightedError::Allocation)?;
    values.resize(len, 0.0);
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispatched_mass_preserves_portable_contraction_bits() {
        let space = TensorSpace::new(2, 3);
        let mass = space
            .prepare_weighted_mass(
                Dims::NONE,
                5,
                1_000_000,
                &mut |_, p| Ok::<_, WeightedError>(1.0 + p[0] * p[1] + p[2]),
                &mut || false,
            )
            .unwrap();
        let input: Vec<_> = (0..space.ndof()).map(|i| (i % 19) as f64 - 7.5).collect();
        let mut portable = vec![0.0; space.ndof()];
        let mut dispatched = vec![0.0; space.ndof()];
        let mut workspace = vec![0.0; mass.workspace_len()];
        mass.apply_body(&input, &mut portable, &mut workspace, &mut || false)
            .unwrap();
        mass.apply_into(&input, &mut dispatched, &mut workspace, &mut || false)
            .unwrap();
        for (a, b) in portable.iter().zip(&dispatched) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}
