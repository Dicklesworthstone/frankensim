//! The same Nitsche cell blocks and lifting under the existing Q1 constraints.
use super::*;
use crate::elastic3::dirichlet::EmbeddedDirichletOptions3;
use crate::quad3::surface::SurfaceOptions3;
use crate::elastic3::surface::ReferenceLoad3;
use super::enrichment::CellResidual3;

/// Pure displacement law on the already selected reference support patch.
pub type PrescribedMotion3<'a> = dyn Fn([f64; 3], [f64; 3]) -> [f64; 3] + 'a;

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

    /// Exact discrete z^T (dK/dscale_c) u for independent primal/adjoint fields.
    /// Reconstruct both with the original Q1 constraint map. Nitsche terms are
    /// already in the cell block; each adjacent cell receives HALF the reference
    /// ghost contraction. Do not multiply by current scales: this is dK/dscale.
    ///
    /// Direct bilinear contraction avoids subtracting two large quadratic
    /// energies when primal and adjoint have very different magnitudes. Fields
    /// need not be solved; this operation by itself makes no residual claim.
    pub fn scale_bilinear_forms(&self, z: &[f64], u: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        let u = self.physical_displacements(u)?;
        let z = self.physical_displacements(z)?;
        let mut result = vec![0.0; self.cells()];
        for (id, cell) in self.raw.cells.iter().enumerate() {
            if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
            let local: [f64; 24] = std::array::from_fn(|j| u[3*cell.nodes[j/3]+j%3]);
            for i in 0..24 {
                let applied: f64 = cell.stiffness[i].iter().zip(&local).map(|(k,u)| k*u).sum();
                result[id] += z[3*cell.nodes[i/3]+i%3]*applied;
            }
        }
        for face in &self.raw.ghosts {
            if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
            let mut value = 0.0;
            for c in 0..3 {
                let ju: f64 = face.nodes.iter().zip(&face.jump).map(|(&n,&j)| j*u[3*n+c]).sum();
                let jz: f64 = face.nodes.iter().zip(&face.jump).map(|(&n,&j)| j*z[3*n+c]).sum();
                value += face.weight*jz*ju;
            }
            for &cell in &face.cells { result[cell] += 0.5*value; }
        }
        if !result.iter().all(|v| v.is_finite()) { return Err(ElasticityError3::Invalid("bilinear scale contraction overflow")); }
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        Ok(result)
    }

    /// Assemble external reference forces plus the CURRENT material's Nitsche
    /// lifting. None preserves the original reference-load path. This is a
    /// variational RHS, not a physical traction or actuator-work definition.
    pub fn reference_load_with_motion(&self, load: ReferenceLoad3<'_>,
        prescribed: Option<&PrescribedMotion3<'_>>,
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        if prescribed.is_some() && self.embedded_dirichlet_penalty().is_none() {
            return Err(ElasticityError3::Invalid("motion requires embedded Dirichlet support"));
        }
        let mut rhs = self.reference_load(load, &mut checkpoint)?;
        if let Some(g) = prescribed {
            let lifting = self.prescribed_displacement_load(g, &mut checkpoint)?;
            for (i, (r, b)) in rhs.iter_mut().zip(lifting).enumerate() {
                if i % 192 == 0 && checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
                *r += b;
                if !r.is_finite() { return Err(ElasticityError3::Invalid("motion RHS overflow")); }
            }
        }
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        Ok(rhs)
    }

    /// Localize f(w)+b_g(w)-a(u,w) with the SAME retained bulk, surface,
    /// Nitsche and ghost terms used in the solve. The lifting is linear in each
    /// cell scale: b_g(w) = sum_c scale_c * (w^T db_g/dscale_c). Its existing
    /// exact pullback therefore supplies localization without a second trace
    /// integrator. Adjoint equations use homogeneous motion (None).
    pub fn cell_residuals_with_motion(&self, u: &[f64], w: &[f64],
        load: ReferenceLoad3<'_>, prescribed: Option<&PrescribedMotion3<'_>>,
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<CellResidual3>, ElasticityError3> {
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        if prescribed.is_some() && self.embedded_dirichlet_penalty().is_none() {
            return Err(ElasticityError3::Invalid("motion requires embedded Dirichlet support"));
        }
        let mut residuals = self.cell_reference_residuals(u, w, load, &mut checkpoint)?;
        if let Some(g) = prescribed {
            let work = self.prescribed_displacement_scale_work(g, w, &mut checkpoint)?;
            for (i, ((r, b), &scale)) in residuals.iter_mut().zip(work).zip(self.scales()).enumerate() {
                if i % 64 == 0 && checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
                r.load += scale*b;
                if !r.load.is_finite() || !r.residual().is_finite() {
                    return Err(ElasticityError3::Invalid("motion residual overflow"));
                }
            }
        }
        if checkpoint().is_break() { return Err(ElasticityError3::Cancelled); }
        Ok(residuals)
    }
}
