//! Stress observations through the original conforming hanging-node map.
use super::*;
use crate::elastic3::stress::{AffineBulkStress3, BulkStressPoint3, BulkStressPullback3};

impl AdaptiveElasticity3 {
    /// Evaluate original bulk rules after exact constrained reconstruction.
    /// Cell ids use `leaves()` ordering. No nodal stress recovery or projection
    /// is inserted. Reconstruction is bounded by the admitted physical DOFs;
    /// stress evaluation then polls at each retained quadrature point.
    pub fn bulk_stress(
        &self,
        u: &[f64],
        max_points: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<BulkStressPoint3>, ElasticityError3> {
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        let full = self.physical_displacements(u)?;
        self.raw.bulk_stress(&full, max_points, checkpoint)
    }

    /// Apply the same conforming reconstruction to the baseline and every
    /// coordinate direction, then delegate whole-ball source-domain admission
    /// and stress sampling to the retained native owner. The total point-field
    /// allowance is checked before allocating any reconstructed field family.
    pub fn affine_bulk_stress(
        &self,
        offset: &[f64],
        basis: &[&[f64]],
        max_point_fields: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<AffineBulkStress3, ElasticityError3> {
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        self.raw
            .affine_stress_point_count(basis.len(), max_point_fields, &mut checkpoint)?;
        let full_offset = self.physical_displacements(offset)?;
        let mut full_basis = Vec::with_capacity(basis.len());
        for &field in basis {
            if checkpoint().is_break() {
                return Err(ElasticityError3::Cancelled);
            }
            full_basis.push(self.physical_displacements(field)?);
        }
        let borrowed: Vec<_> = full_basis.iter().map(Vec::as_slice).collect();
        self.raw
            .affine_bulk_stress(&full_offset, &borrowed, max_point_fields, checkpoint)
    }

    /// Exact physical-stress VJP in independent master coordinates. The direct
    /// material-scale derivative is unchanged by T; state cotangents use T^T.
    /// The supplied cotangents already include any desired volume weights.
    pub fn bulk_stress_pullback(
        &self,
        u: &[f64],
        derivatives: &[[f64; 6]],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<BulkStressPullback3, ElasticityError3> {
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        let full = self.physical_displacements(u)?;
        let mut result = self
            .raw
            .bulk_stress_pullback(&full, derivatives, &mut checkpoint)?;
        let mut displacement = vec![0.0; self.n()];
        self.restrict(&result.displacement, &mut displacement);
        if !displacement.iter().all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid(
                "reduced bulk stress pullback overflow",
            ));
        }
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        result.displacement = displacement;
        Ok(result)
    }

    /// Exact reference-material stress transpose, with hanging constraints and
    /// homogeneous master clamps applied. No current material-scale factor.
    pub fn reference_bulk_stress_pullback(
        &self,
        derivatives: &[[f64; 6]],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        let full = self
            .raw
            .reference_bulk_stress_pullback(derivatives, &mut checkpoint)?;
        let mut displacement = vec![0.0; self.n()];
        self.restrict(&full, &mut displacement);
        if !displacement.iter().all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid(
                "reduced reference stress pullback overflow",
            ));
        }
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        Ok(displacement)
    }
}
