//! Stress observations through the original conforming hanging-node map.
use super::*;
use crate::elastic3::stress::{AffineBulkStress3, BulkStressPoint3, BulkStressPullback3};

impl AdaptiveElasticity3 {
    /// Check one physical field against the original source-material domain,
    /// with exact hanging reconstruction and homogeneous master clamps.
    pub fn admit_material_domain(
        &self,
        u: &[f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_affine_material_domain(u, &[], checkpoint)
    }

    /// Admit the complete coefficient unit ball on the retained bulk points.
    /// Reconstruct only the current cell's 24 physical coefficients using the
    /// SAME row owner as `physical_displacements`; no expanded nodal fields or
    /// stress samples are allocated. The native owner checks finite inputs and
    /// source ceilings, with bounded cancellation and no operator mutation.
    pub fn admit_affine_material_domain(
        &self,
        offset: &[f64],
        basis: &[&[f64]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.raw.admit_affine_material_domain_mapped(
            Some(offset),
            basis,
            self.n(),
            |cell, field| {
                std::array::from_fn(|i| self.physical_coefficient(field, cell.nodes[i / 3], i % 3))
            },
            checkpoint,
        )
    }

    /// Admit the complete zero-offset coefficient unit ball with exact hanging
    /// reconstruction and master clamps. No global zero field or reconstructed
    /// family is allocated; source checks share the affine admission owner.
    pub fn admit_load_space_material_domain(
        &self,
        basis: &[&[f64]],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.raw.admit_affine_material_domain_mapped(
            None,
            basis,
            self.n(),
            |cell, field| {
                std::array::from_fn(|i| self.physical_coefficient(field, cell.nodes[i / 3], i % 3))
            },
            checkpoint,
        )
    }

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
