//! Material inertia in the original Nedelec and RT moment coordinates.

use super::{ElementBasis, Family, VecSpace, dot3, mpow};
use crate::highorder::material::scaled_entry;
use crate::highorder::simplex::duffy_quadrature;
use crate::weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, assemble_triplets, poll};
use fs_qty::Dims;
use fs_sparse::Csr;

mod admit;
mod field;

const MAX_LOCAL: usize = 84;
const MAX_MONOMIAL: usize = 35;

/// Prepared consistent vector mass for scalar, P1 scalar, or Cartesian SPD
/// material coefficients, optionally multiplied by a frozen positive spatial
/// scalar field, on the existing H(curl)/H(div) moment layout.
///
/// The physical pairing is `integral test . W trial dV`. Density on reference
/// geometry gives inertia and kinetic energy; other material coefficients
/// retain their declared dimensions. Cell-local weights may jump at interfaces.
/// No incidence, trace orientation, basis normalization or pressure/gauge is
/// changed. Material positivity alone does not certify a complete PDE model.
///
/// Warm action contracts local DOFs into three polynomials, samples each trial
/// vector once per quadrature point, then projects accumulated polynomial
/// moments back to the original DOFs. With n local DOFs, s scalar monomials and
/// Q points, work is O(n*s + Q*s). No dense element matrix or warm allocation
/// is required. Optional CSR assembly uses the same quadrature and coefficient.
pub struct VectorMaterialMass<'a, 'c> {
    space: &'a VecSpace<'c>,
    positions: &'a [[f64; 3]],
    weights: &'a [CellWeight],
    dofs: Vec<usize>,
    vol6: Vec<f64>,
    local_count: usize,
    quadrature: Vec<([f64; 4], f64)>,
    scalar_field: Option<field::SampledScalarField>,
    diagonal: Vec<f64>,
    coefficient_dims: Dims,
    limits: WeightedAssemblyLimits,
}

impl<'c> VecSpace<'c> {
    /// Prepare material mass using the positions that originally built this
    /// space. Reference density belongs on the original reference geometry.
    ///
    /// Scalar and affine coefficients must be positive; tensors must pass the
    /// existing interval SPD admission. The r+3 Duffy rule exactly integrates
    /// these coefficients against both native vector families at orders 1..=4.
    /// Geometry/layout and coefficient checks precede numerical preparation;
    /// the public basis tables must still describe their original space.
    /// This checks finite charts and layout, not an independent duality proof.
    ///
    /// `max_storage_bytes` bounds retained quadrature, local-to-global indices,
    /// volumes and diagonal. Small order-bounded construction scratch is extra.
    /// `max_triplets` applies only when `assemble` is requested. Input fields,
    /// geometry and materials stay borrowed immutably for the plan's lifetime.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_material_mass<'a>(
        &'a self,
        positions: &'a [[f64; 3]],
        weights: &'a [CellWeight],
        coefficient_dims: Dims,
        limits: WeightedAssemblyLimits,
        max_storage_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorMaterialMass<'a, 'c>, WeightedError> {
        poll(cancelled)?;
        coefficient_dims
            .checked_plus(Dims([3, 0, 0, 0, 0, 0]))
            .ok_or(WeightedError::DimensionOverflow)?;
        let admitted = admit::admit(
            self,
            positions,
            weights,
            limits,
            None,
            max_storage_bytes,
            cancelled,
        )?;
        let mut mass = VectorMaterialMass {
            space: self,
            positions,
            weights,
            dofs: admitted.dofs,
            vol6: admitted.vol6,
            local_count: admitted.local_count,
            quadrature: duffy_quadrature(admitted.quadrature_order),
            scalar_field: None,
            diagonal: zeros(self.ndof)?,
            coefficient_dims,
            limits,
        };
        mass.prepare_diagonal(cancelled)?;
        poll(cancelled)?;
        Ok(mass)
    }
}

impl VectorMaterialMass<'_, '_> {
    /// Original global moment coordinate count, including isolated entities.
    #[must_use]
    pub const fn ndof(&self) -> usize {
        self.space.ndof
    }

    /// Material dimensions, before geometric integration and moment scaling.
    #[must_use]
    pub const fn coefficient_dims(&self) -> Dims {
        self.coefficient_dims
    }

    /// Dimensions of coordinate `index` for a physical vector field.
    /// RT uses mean normal/cell moments (field units). Nedelec edge and face
    /// moments have field-times-length units; its interior moments have field
    /// units. This preserves the native mixed convention at orders >=3.
    pub fn coordinate_dims(&self, index: usize, field: Dims) -> Result<Dims, WeightedError> {
        field
            .checked_plus(Dims([self.moment_power(index)?, 0, 0, 0, 0, 0]))
            .ok_or(WeightedError::DimensionOverflow)
    }

    /// Dimensions of one mass entry, accounting for BOTH moment coordinates.
    /// A single uniform entry dimension would be wrong for high-order Nedelec.
    pub fn entry_dims(&self, row: usize, col: usize) -> Result<Dims, WeightedError> {
        self.coefficient_dims
            .checked_plus(Dims([
                3 - self.moment_power(row)? - self.moment_power(col)?,
                0,
                0,
                0,
                0,
                0,
            ]))
            .ok_or(WeightedError::DimensionOverflow)
    }

    /// Dimensions of one dual output for the declared physical input field.
    pub fn output_dims(&self, row: usize, field: Dims) -> Result<Dims, WeightedError> {
        let power = self.moment_power(row)?;
        self.coefficient_dims
            .checked_plus(field)
            .and_then(|dims| dims.checked_plus(Dims([3 - power, 0, 0, 0, 0, 0])))
            .ok_or(WeightedError::DimensionOverflow)
    }

    /// Dimensions of `u^T M u` for the declared physical vector field.
    pub fn energy_dims(&self, field: Dims) -> Result<Dims, WeightedError> {
        self.coefficient_dims
            .checked_plus(field)
            .and_then(|dims| dims.checked_plus(field))
            .and_then(|dims| dims.checked_plus(Dims([3, 0, 0, 0, 0, 0])))
            .ok_or(WeightedError::DimensionOverflow)
    }

    /// Diagonal of the same consistent quadrature operator, for Jacobi scaling.
    /// Unused global entities keep zero; no inverse/coercivity claim is implied.
    #[must_use]
    pub fn diagonal(&self) -> &[f64] {
        &self.diagonal
    }

    fn moment_power(&self, index: usize) -> Result<i8, WeightedError> {
        if index >= self.ndof() {
            return Err(WeightedError::InvalidGeometry {
                cell: 0,
                reason: "vector moment index outside the original layout",
            });
        }
        Ok(i8::from(
            self.space.family == Family::Nedelec && index < self.space.cell_off,
        ))
    }

    fn prepare_diagonal(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WeightedError> {
        let mut monomials = [0.0; MAX_MONOMIAL];
        let mut values = [[0.0; 3]; MAX_LOCAL];
        for (cell, &weight) in self.weights.iter().enumerate() {
            poll(cancelled)?;
            let (scale, tensor, affine) = weight.normalized_3d();
            let mut local = [0.0_f64; MAX_LOCAL];
            for (point, &(lambda, w)) in self.quadrature.iter().enumerate() {
                poll(cancelled)?;
                self.sample_monomials(cell, lambda, &mut monomials);
                basis_values(&self.space.elements[cell], &monomials, &mut values);
                let factor = self.point_factor(cell, point, w, affine, lambda);
                for i in 0..self.local_count {
                    local[i] =
                        factor.mul_add(dot3(values[i], tensor_action(tensor, values[i])), local[i]);
                }
            }
            for i in 0..self.local_count {
                let global = self.dofs[cell * self.local_count + i];
                let value = self.restore_scale(cell, scale, local[i]);
                let sum = self.diagonal[global] + value;
                if !value.is_finite() || !sum.is_finite() {
                    return Err(WeightedError::NonFiniteEntry);
                }
                if value <= 0.0 {
                    return Err(WeightedError::InvalidGeometry {
                        cell,
                        reason: "positive vector material diagonal is not representable",
                    });
                }
                self.diagonal[global] = sum;
            }
        }
        Ok(())
    }

    /// Assemble a numeric CSR in the original moment ordering. The plan's
    /// per-coordinate unit methods retain its physical interpretation.
    /// Triplet admission precedes allocation; cancellation publishes no matrix.
    pub fn assemble(&self, cancelled: &mut impl FnMut() -> bool) -> Result<Csr, WeightedError> {
        poll(cancelled)?;
        let count = self
            .weights
            .len()
            .checked_mul(self.local_count * self.local_count)
            .ok_or(WeightedError::SizeOverflow)?;
        self.limits.check(self.weights.len(), self.ndof(), count)?;
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(count)
            .map_err(|_| WeightedError::Allocation)?;
        let mut monomials = [0.0; MAX_MONOMIAL];
        let mut values = [[0.0; 3]; MAX_LOCAL];
        let mut weighted = [[0.0; 3]; MAX_LOCAL];
        for (cell, &weight) in self.weights.iter().enumerate() {
            poll(cancelled)?;
            let (scale, tensor, affine) = weight.normalized_3d();
            let mut local = [0.0_f64; MAX_LOCAL * MAX_LOCAL];
            for (point, &(lambda, w)) in self.quadrature.iter().enumerate() {
                poll(cancelled)?;
                self.sample_monomials(cell, lambda, &mut monomials);
                basis_values(&self.space.elements[cell], &monomials, &mut values);
                for i in 0..self.local_count {
                    weighted[i] = tensor_action(tensor, values[i]);
                }
                let factor = self.point_factor(cell, point, w, affine, lambda);
                for i in 0..self.local_count {
                    for j in i..self.local_count {
                        let entry = &mut local[i * self.local_count + j];
                        *entry = factor.mul_add(dot3(values[i], weighted[j]), *entry);
                    }
                }
            }
            let dofs = self.cell_dofs(cell);
            for (i, &row) in dofs.iter().enumerate() {
                poll(cancelled)?;
                for (j, &col) in dofs.iter().enumerate() {
                    let value = self.restore_scale(
                        cell,
                        scale,
                        local[i.min(j) * self.local_count + i.max(j)],
                    );
                    if !value.is_finite() {
                        return Err(WeightedError::NonFiniteEntry);
                    }
                    entries.push((row, col, value));
                }
            }
        }
        assemble_triplets(self.ndof(), entries, cancelled)
    }

    /// Apply the complete consistent material mass with no allocation.
    /// All slices have `ndof()` entries. Scratch may change on failure; output
    /// changes only after the full action and final cancellation check succeed.
    pub fn apply_into(
        &self,
        input: &[f64],
        output: &mut [f64],
        scratch: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WeightedError> {
        crate::highorder::fma::apply_vector_mass_dispatch(self, input, output, scratch, cancelled)
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    pub(crate) fn apply_body(
        &self,
        input: &[f64],
        output: &mut [f64],
        scratch: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WeightedError> {
        poll(cancelled)?;
        if input.len() != self.ndof() || output.len() != self.ndof() || scratch.len() != self.ndof()
        {
            return Err(WeightedError::InvalidGeometry {
                cell: 0,
                reason: "vector material action requires original-layout vectors",
            });
        }
        for (i, (&value, out)) in input.iter().zip(scratch.iter_mut()).enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            if !value.is_finite() {
                return Err(WeightedError::NonFiniteEntry);
            }
            *out = 0.0;
        }
        let mut monomials = [0.0; MAX_MONOMIAL];
        for (cell, &weight) in self.weights.iter().enumerate() {
            poll(cancelled)?;
            let element = &self.space.elements[cell];
            let s = element.monos.len();
            let dofs = self.cell_dofs(cell);
            let (scale, tensor, affine) = weight.normalized_3d();
            let mut trial = [[0.0_f64; MAX_MONOMIAL]; 3];
            for (function, &global) in element.funcs.iter().zip(dofs) {
                for c in 0..3 {
                    for m in 0..s {
                        trial[c][m] = input[global].mul_add(function.comp[c][m], trial[c][m]);
                    }
                }
            }
            let mut moments = [[0.0_f64; MAX_MONOMIAL]; 3];
            for (point, &(lambda, w)) in self.quadrature.iter().enumerate() {
                poll(cancelled)?;
                self.sample_monomials(cell, lambda, &mut monomials);
                let mut value = [0.0_f64; 3];
                for c in 0..3 {
                    for m in 0..s {
                        value[c] = trial[c][m].mul_add(monomials[m], value[c]);
                    }
                }
                let flux = tensor_action(tensor, value);
                let factor = self.point_factor(cell, point, w, affine, lambda);
                for c in 0..3 {
                    for m in 0..s {
                        moments[c][m] = (factor * monomials[m]).mul_add(flux[c], moments[c][m]);
                    }
                }
            }
            for (function, &global) in element.funcs.iter().zip(dofs) {
                let mut local = 0.0_f64;
                for c in 0..3 {
                    for m in 0..s {
                        local = function.comp[c][m].mul_add(moments[c][m], local);
                    }
                }
                let value = self.restore_scale(cell, scale, local);
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

    #[inline]
    fn cell_dofs(&self, cell: usize) -> &[usize] {
        &self.dofs[cell * self.local_count..(cell + 1) * self.local_count]
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn point_factor(
        &self,
        cell: usize,
        point: usize,
        weight: f64,
        affine: Option<[f64; 4]>,
        lambda: [f64; 4],
    ) -> f64 {
        let base = quadrature_factor(weight * self.vol6[cell], affine, lambda);
        self.scalar_field.as_ref().map_or(base, |field| {
            base * field.values[cell * self.quadrature.len() + point]
        })
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn restore_scale(&self, cell: usize, material: f64, value: f64) -> f64 {
        scaled_entry(
            material,
            self.scalar_field.as_ref().map(|field| field.scales[cell]),
            value,
        )
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn sample_monomials(&self, cell: usize, lambda: [f64; 4], output: &mut [f64; MAX_MONOMIAL]) {
        let mut point = [0.0_f64; 3];
        for (a, &vertex) in self.space.complex.tets[cell].iter().enumerate() {
            for c in 0..3 {
                point[c] = lambda[a].mul_add(self.positions[vertex as usize][c], point[c]);
            }
        }
        let element = &self.space.elements[cell];
        let xi = element.chart.local(point);
        for (m, &[a, b, c]) in element.monos.iter().enumerate() {
            output[m] = mpow(xi[0], a) * mpow(xi[1], b) * mpow(xi[2], c);
        }
    }
}

fn basis_values(
    element: &ElementBasis,
    monomials: &[f64; MAX_MONOMIAL],
    output: &mut [[f64; 3]; MAX_LOCAL],
) {
    for (i, function) in element.funcs.iter().enumerate() {
        output[i] = [0.0; 3];
        for m in 0..element.monos.len() {
            for c in 0..3 {
                output[i][c] = function.comp[c][m].mul_add(monomials[m], output[i][c]);
            }
        }
    }
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn quadrature_factor(base: f64, affine: Option<[f64; 4]>, lambda: [f64; 4]) -> f64 {
    affine.map_or(base, |rho| {
        base * rho.iter().zip(lambda).map(|(r, l)| r * l).sum::<f64>()
    })
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn tensor_action(tensor: [[f64; 3]; 3], value: [f64; 3]) -> [f64; 3] {
    tensor.map(|row| row[0].mul_add(value[0], row[1].mul_add(value[1], row[2] * value[2])))
}

fn zeros(n: usize) -> Result<Vec<f64>, WeightedError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| WeightedError::Allocation)?;
    values.resize(n, 0.0);
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_material_dispatch_preserves_portable_action_bits() {
        let (complex, positions) = crate::fixtures::single_tet();
        let weights = [CellWeight::AffineScalar([2.0, 3.0, 5.0, 7.0])];
        for family in [Family::Nedelec, Family::Rt] {
            let space = VecSpace::new(&complex, &positions, 3, family);
            let mass = space
                .prepare_material_mass(
                    &positions,
                    &weights,
                    Dims::NONE,
                    WeightedAssemblyLimits {
                        max_cells: 1,
                        max_dofs: 100,
                        max_triplets: 0,
                    },
                    100_000,
                    &mut || false,
                )
                .unwrap();
            let input: Vec<_> = (0..space.ndof).map(|i| (i % 13) as f64 - 5.5).collect();
            let mut portable = vec![0.0; space.ndof];
            let mut dispatched = vec![0.0; space.ndof];
            let mut scratch = vec![0.0; space.ndof];
            mass.apply_body(&input, &mut portable, &mut scratch, &mut || false)
                .unwrap();
            mass.apply_into(&input, &mut dispatched, &mut scratch, &mut || false)
                .unwrap();
            for (a, b) in portable.iter().zip(dispatched) {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }
    }
}
