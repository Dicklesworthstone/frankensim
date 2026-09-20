//! The same Nitsche cell blocks and lifting under the existing Q1 constraints.
use super::*;
use crate::elastic3::dirichlet::EmbeddedDirichletOptions3;
use crate::quad3::surface::SurfaceOptions3;

impl AdaptiveElasticity3 {
    /// Homogeneous full-vector support on a sampled implicit-surface patch.
    /// Strong zero box clamps are optional. Surface geometry, bulk physics and
    /// ghost terms come from the same original builders; the new boundary terms
    /// enter their retained cell matrices BEFORE hanging-node elimination.
    ///
    /// `patch` must be pure and consistently selected across refinement levels.
    /// Method/penalty mismatch refuses transfer; patch/geometry equivalence is
    /// still a caller obligation. Resolve narrow patches and junctions explicitly.
    /// The positive penalty is not a coercivity or component-support certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn build_with_embedded_dirichlet(domain: HexCell, tree: &Octree3, sdf: &dyn CutSdf3,
        material: &IsotropicElastic, box_clamp: &dyn Fn([f64; 3]) -> bool,
        patch: &dyn Fn([f64; 3], [f64; 3]) -> bool,
        options: ElasticityOptions3, boundary: EmbeddedDirichletOptions3,
        surface: SurfaceOptions3, control: &mut QuadratureControl3<'_>) -> Result<Self, ElasticityError3> {
        boundary.validate()?;
        surface.validate()?;
        let mut op = Self::build_core(domain, tree, sdf, material, box_clamp, options, false, control)?;
        op.raw.integrate_surface(sdf, surface, control)?;
        op.raw.attach_embedded_dirichlet(material, patch, boundary, control)?;
        op.reference[3] = boundary.beta;
        Ok(op)
    }
    #[must_use]
    pub fn embedded_dirichlet_penalty(&self) -> Option<f64> { self.raw.embedded_dirichlet_penalty() }
    #[must_use]
    pub fn embedded_dirichlet_area(&self) -> Option<f64> { self.raw.embedded_dirichlet_area() }

    /// b_g reduced with the SAME T^T used by stiffness and other loads. Add it
    /// to external loads before solving; recompute after any density update.
    /// It is not an applied traction/reaction or a physical actuator-work value.
    pub fn prescribed_displacement_load(&self,
        prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        let raw = self.raw.prescribed_displacement_load(prescribed, &mut checkpoint)?;
        let mut rhs = vec![0.0; self.n()];
        self.restrict(&raw, &mut rhs);
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        if !rhs.iter().all(|x| x.is_finite()) { return Err(ElasticityError3::Invalid("reduced prescribed RHS overflow")); }
        Ok(rhs)
    }
    /// v^T db_g/dscale_c, using the physical reconstruction T v. Combine with
    /// the existing stiffness contraction for the chosen objective; a nonzero
    /// prescribed g is NOT a fixed independent load in the OC driver.
    pub fn prescribed_displacement_scale_work(&self,
        prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3], v: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        let physical = self.physical_displacements(v)?;
        self.raw.prescribed_displacement_scale_work(prescribed, &physical, checkpoint)
    }
}
