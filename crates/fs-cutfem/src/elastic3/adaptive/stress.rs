//! Stress observations through the original conforming hanging-node map.
use super::*;
use crate::elastic3::stress::{BulkStressPoint3, BulkStressPullback3};

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
