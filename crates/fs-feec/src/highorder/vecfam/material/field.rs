//! Frozen spatial density on the original vector material operator.

use super::{VectorMaterialMass, admit, zeros};
use crate::highorder::material::MAX_FIELD_QUADRATURE;
use crate::highorder::simplex::duffy_quadrature;
use crate::highorder::vecfam::VecSpace;
use crate::weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, poll};
use fs_qty::Dims;

pub(super) struct SampledScalarField {
    pub(super) values: Vec<f64>,
    pub(super) scales: Vec<f64>,
}

impl<'c> VecSpace<'c> {
    /// Prepare vector mass with a frozen positive spatial scalar multiplier.
    ///
    /// The physical pairing is `integral test . W trial * field dV`. A unit
    /// cell weight and a reference-density field supply general graded inertia;
    /// an SPD tensor can retain anisotropy while the scalar field varies inside
    /// each cell. Original reference density belongs on reference geometry for
    /// a closed body's inventory. This is not evolving-density transport.
    ///
    /// `field_at` receives each original cell and its stored-corner-order
    /// barycentric point exactly once per quadrature sample. Its error type is
    /// preserved. The source and callback are not retained or queried by mass
    /// action, diagonal access, or CSR assembly. Changed source state requires
    /// preparing another operator. All three use the same frozen samples.
    ///
    /// `quadrature_order` is the Gauss count per Duffy axis, from `r+3` through
    /// [`MAX_FIELD_QUADRATURE`]. The caller owns quadrature convergence for its
    /// actual field. Positive samples establish a positive discrete weight;
    /// they do not certify the source away from these points, resolve an
    /// unaligned interface exactly, or prove a complete PDE formulation stable.
    ///
    /// Geometry, cell weights, dimensions and total retained storage are
    /// admitted before any source query. The storage bound includes one scalar
    /// per cell/point and one independent normalization scale per cell, in
    /// addition to the ordinary vector mass plan. Source and material scales
    /// stay separate until after integration. Warm action remains allocation
    /// free and O(n*s + Q*s) per cell in the original moment coordinates.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_material_mass_with_scalar_field<'a, E: From<WeightedError>>(
        &'a self,
        positions: &'a [[f64; 3]],
        weights: &'a [CellWeight],
        coefficient_dims: Dims,
        field_dims: Dims,
        quadrature_order: usize,
        limits: WeightedAssemblyLimits,
        max_storage_bytes: usize,
        field_at: &mut impl FnMut(usize, [f64; 4]) -> Result<f64, E>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorMaterialMass<'a, 'c>, E> {
        poll(cancelled)?;
        let coefficient_dims = coefficient_dims
            .checked_plus(field_dims)
            .ok_or(WeightedError::DimensionOverflow)?;
        coefficient_dims
            .checked_plus(Dims([3, 0, 0, 0, 0, 0]))
            .ok_or(WeightedError::DimensionOverflow)?;
        let admitted = admit::admit(
            self,
            positions,
            weights,
            limits,
            Some(quadrature_order),
            max_storage_bytes,
            cancelled,
        )?;
        let quadrature = duffy_quadrature(admitted.quadrature_order);
        let count = weights
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
            .try_reserve_exact(weights.len())
            .map_err(|_| WeightedError::Allocation)?;
        let diagonal = zeros(self.ndof)?;
        for cell in 0..weights.len() {
            let start = field.values.len();
            let mut scale = 0.0_f64;
            for &(lambda, _) in &quadrature {
                poll(cancelled)?;
                let value = field_at(cell, lambda)?;
                if !value.is_finite() || value <= 0.0 {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "sampled vector scalar field must be finite and positive",
                    }
                    .into());
                }
                scale = scale.max(value);
                field.values.push(value);
            }
            for (point, value) in field.values[start..].iter_mut().enumerate() {
                if point % 256 == 0 {
                    poll(cancelled)?;
                }
                *value /= scale;
                if *value <= 0.0 {
                    return Err(WeightedError::InvalidCoefficient {
                        cell,
                        reason: "sampled vector scalar contrast is not representable",
                    }
                    .into());
                }
            }
            field.scales.push(scale);
        }
        let mut mass = VectorMaterialMass {
            space: self,
            positions,
            weights,
            dofs: admitted.dofs,
            vol6: admitted.vol6,
            local_count: admitted.local_count,
            quadrature,
            scalar_field: Some(field),
            diagonal,
            coefficient_dims,
            limits,
        };
        mass.prepare_diagonal(cancelled)?;
        poll(cancelled)?;
        Ok(mass)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::highorder::vecfam::Family;

    #[test]
    fn sampled_vector_density_dispatch_preserves_portable_bits() {
        let (complex, positions) = crate::fixtures::single_tet();
        let weights = [CellWeight::AffineScalar([2.0, 3.0, 5.0, 7.0])];
        for family in [Family::Nedelec, Family::Rt] {
            let space = VecSpace::new(&complex, &positions, 3, family);
            let mass = space
                .prepare_material_mass_with_scalar_field(
                    &positions,
                    &weights,
                    Dims::NONE,
                    Dims::NONE,
                    7,
                    WeightedAssemblyLimits {
                        max_cells: 1,
                        max_dofs: 100,
                        max_triplets: 0,
                    },
                    100_000,
                    &mut |_, l| -> Result<f64, WeightedError> { Ok(1.0 + l[0] * l[2]) },
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
