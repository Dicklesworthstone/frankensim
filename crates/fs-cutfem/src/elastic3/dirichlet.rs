//! Symmetric Nitsche displacement conditions on a retained implicit surface.
//!
//! For basis v_i, t_i=sigma(v_i)n and gamma=beta*(lambda+2*mu)/h:
//! K_Gamma[i,j] = integral(-v_i.t_j - t_i.v_j + gamma*v_i.v_j).
//! b_g[i] = integral(-t_i.g + gamma*v_i.g). Each cell's reference
//! K_Gamma and b_g scale by its CURRENT material multiplier. Bulk and ghost
//! actions, hanging-node elimination and linear solvers are not replaced.
//!
//! The supplied patch is sampled at the retained numerical surface points.
//! This is not a certified patch measure, a guaranteed coercivity constant,
//! follower condition or shape derivative. Users must resolve patch junctions
//! and support every connected component. An insufficient penalty/support may
//! produce an indefinite/singular system; small residuals do not prove stability.
use super::*;
use crate::quad3::surface::{SurfaceOptions3, SurfacePoint3};

mod reaction;
pub use reaction::EmbeddedReaction3;

/// Explicit Nitsche penalty; it never uses the possibly tiny cut volume.
#[derive(Debug, Clone, Copy)]
pub struct EmbeddedDirichletOptions3 {
    /// Applied coefficient is beta*(lambda+2*mu)/min(background cell spans).
    /// A positive value is required, not a theorem that this value is sufficient.
    pub beta: f64,
}
impl Default for EmbeddedDirichletOptions3 {
    fn default() -> Self { Self { beta: 32.0 } }
}
impl EmbeddedDirichletOptions3 {
    pub(crate) fn validate(self) -> Result<(), ElasticityError3> {
        if !self.beta.is_finite() || self.beta <= 0.0 {
            return Err(ElasticityError3::Invalid("Nitsche beta must be finite and positive"));
        }
        Ok(())
    }
}

// Retain indices into existing surface quadrature, not copies of 24-vector
// trace stencils. This immutable selection is shared by every displacement law.
pub(super) struct DirichletData3 {
    points: Vec<Vec<usize>>,
    lame: [f64; 2],
    beta: f64,
    area: f64,
}
fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() { Err(ElasticityError3::Cancelled) } else { Ok(()) }
}
fn penalty(cell: &Cell3, lame: [f64; 2], beta: f64) -> Result<f64, ElasticityError3> {
    let h = (0..3).map(|a| cell.bounds.hi()[a]-cell.bounds.lo()[a]).fold(f64::INFINITY, f64::min);
    let gamma = beta * ((lame[0]+2.0*lame[1])/h);
    if !gamma.is_finite() || gamma <= 0.0 { return Err(ElasticityError3::Invalid("Nitsche penalty overflow")); }
    Ok(gamma)
}
// Stress of N_a e_c, contracted with the numerical outward unit normal.
fn traces(cell: &Cell3, point: &SurfacePoint3, lame: [f64; 2]) -> ([f64; 8], [[f64; 3]; 24]) {
    let (n, gradient) = q1(cell.bounds, point.position);
    let [lambda, mu] = lame;
    let traction = std::array::from_fn(|i| {
        let (a, c) = (i/3, i%3);
        let dn: f64 = gradient[a].iter().zip(point.normal).map(|(g, n)| g*n).sum();
        std::array::from_fn(|d| lambda*gradient[a][c]*point.normal[d]
            + mu*(point.normal[c]*gradient[a][d] + if c==d {dn} else {0.0}))
    });
    (n, traction)
}

impl CutElasticity3 {
    /// Build actual implicit-surface supports, with optional additional zero
    /// box clamps. `patch(p,n)` selects full-vector homogeneous displacement;
    /// all other embedded surface points remain natural. The box clamp may
    /// select nothing: support admission is deferred until this patch exists.
    ///
    /// A pure SDF supplies BOTH bulk and surface geometry. Only an unpublished
    /// operator is modified; cancellation/error returns no partially supported
    /// model. Existing surface/body load laws are unchanged: callers must return
    /// zero Neumann traction on the selected Dirichlet patch. Use
    /// `prescribed_displacement_load` for nonzero g and add it to external loads.
    #[allow(clippy::too_many_arguments)]
    pub fn build_with_embedded_dirichlet(
        domain: HexCell, counts: [usize; 3], sdf: &dyn CutSdf3,
        material: &IsotropicElastic, box_clamp: &dyn Fn([f64; 3]) -> bool,
        patch: &dyn Fn([f64; 3], [f64; 3]) -> bool,
        options: ElasticityOptions3, boundary: EmbeddedDirichletOptions3,
        surface: SurfaceOptions3, control: &mut QuadratureControl3<'_>,
    ) -> Result<Self, ElasticityError3> {
        boundary.validate()?;
        surface.validate()?;
        let mut op = Self::build_core(domain, counts, sdf, material, box_clamp, options, false, control)?;
        op.integrate_surface(sdf, surface, control)?;
        op.attach_embedded_dirichlet(material, patch, boundary, control)?;
        Ok(op)
    }

    // Only called by unpublished geometry builders, once. The same modified
    // cell blocks feed applies, exact diagonals, density contractions and MG.
    pub(super) fn attach_embedded_dirichlet(&mut self, material: &IsotropicElastic,
        patch: &dyn Fn([f64; 3], [f64; 3]) -> bool, options: EmbeddedDirichletOptions3,
        control: &mut QuadratureControl3<'_>) -> Result<(), ElasticityError3> {
        control.poll()?;
        if self.embedded.is_some() { return Err(ElasticityError3::Invalid("embedded supports already attached")); }
        let (lambda, mu) = material.lame();
        let mut data = DirichletData3 { points: Vec::with_capacity(self.cells.len()), lame: [lambda, mu], beta: options.beta, area: 0.0 };
        for cell in &mut self.cells {
            control.poll()?;
            let gamma = penalty(cell, data.lame, options.beta)?;
            let rule = cell.rules.surface().ok_or(ElasticityError3::Invalid("embedded support requires retained surface rules"))?;
            let mut selected = Vec::new();
            for (index, point) in rule.points().iter().enumerate() {
                control.poll()?;
                let active = patch(point.position, point.normal);
                control.poll()?;
                if !active { continue; }
                selected.push(index);
                data.area += point.weight;
                let (n, t) = traces(cell, point, data.lame);
                for i in 0..24 { for j in i..24 {
                    let value = -n[i/3]*t[j][i%3] - n[j/3]*t[i][j%3]
                        + if i%3==j%3 {gamma*n[i/3]*n[j/3]} else {0.0};
                    cell.stiffness[i][j] += point.weight*value;
                } }
            }
            for i in 0..24 { for j in i..24 {
                if !cell.stiffness[i][j].is_finite() { return Err(ElasticityError3::Invalid("Nitsche stiffness overflow")); }
                cell.stiffness[j][i] = cell.stiffness[i][j];
            } }
            data.points.push(selected);
        }
        if !data.area.is_finite() || data.area <= 0.0 {
            return Err(ElasticityError3::Invalid("embedded support has no positive numerical area"));
        }
        control.poll()?;
        self.embedded = Some(data);
        Ok(())
    }

    /// None for legacy natural-interface operators. This distinguishes the
    /// method/penalty; it does NOT certify equivalence of two selected patches.
    #[must_use]
    pub fn embedded_dirichlet_penalty(&self) -> Option<f64> { self.embedded.as_ref().map(|d| d.beta) }
    /// Retained numerical support area, not a coercivity or connectivity proof.
    #[must_use]
    pub fn embedded_dirichlet_area(&self) -> Option<f64> { self.embedded.as_ref().map(|d| d.area) }

    fn displacement_cells(&self, prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<Vec<[f64; 24]>, ElasticityError3> {
        poll(checkpoint)?;
        let data = self.embedded.as_ref().ok_or(ElasticityError3::Invalid("operator has no embedded Dirichlet support"))?;
        let mut loads = vec![[0.0;24]; self.cells.len()];
        for (id, cell) in self.cells.iter().enumerate() {
            poll(checkpoint)?;
            if data.points[id].is_empty() { continue; }
            let gamma = penalty(cell, data.lame, data.beta)?;
            let rule = cell.rules.surface().expect("builder retained surface");
            for &index in &data.points[id] {
                poll(checkpoint)?;
                let point = &rule.points()[index];
                let g = prescribed(point.position, point.normal);
                poll(checkpoint)?;
                if !g.iter().all(|x| x.is_finite()) { return Err(ElasticityError3::Invalid("nonfinite prescribed displacement")); }
                let (n, t) = traces(cell, point, data.lame);
                for i in 0..24 {
                    let tg: f64 = t[i].iter().zip(g).map(|(t,g)| t*g).sum();
                    loads[id][i] += point.weight*(-tg + gamma*n[i/3]*g[i%3]);
                }
            }
            if !loads[id].iter().all(|x| x.is_finite()) { return Err(ElasticityError3::Invalid("prescribed displacement load overflow")); }
        }
        poll(checkpoint)?;
        Ok(loads)
    }

    /// RHS lifting b_g for the CURRENT density. Add this to body/traction loads
    /// before solving. Recompute after EVERY scale change; g must be pure and
    /// compatible with any strong zero box clamps. The selected support is
    /// fixed, while g may differ between independent experiments.
    ///
    /// This is a variational lifting, not an applied physical traction or a
    /// reaction vector. A solver's (f+b_g)^T u is an augmented load functional,
    /// not external-load compliance f^T u or actuator work.
    pub fn prescribed_displacement_load(&self,
        prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        let cells = self.displacement_cells(prescribed, &mut checkpoint)?;
        let mut rhs = vec![0.0; self.n()];
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            for i in 0..24 { if !self.fixed[cell.nodes[i/3]] {
                rhs[3*cell.nodes[i/3]+i%3] += self.scales[id]*cells[id][i];
            } }
        }
        if !rhs.iter().all(|x| x.is_finite()) { return Err(ElasticityError3::Invalid("prescribed RHS overflow")); }
        poll(&mut checkpoint)?;
        Ok(rhs)
    }

    /// Exact discrete v^T db_g/dscale_c for every cell. No differentiation
    /// through the linear solver and no assumed density-independent lifting.
    /// For J=(f+b_g)^T u with fixed external f, dJ/dscale = 2*load_work(u)
    /// - scale_quadratic_forms(u). For an external-only observation q^T u,
    /// solve K z=q and use load_work(z)-z^T(dK/dscale)u instead.
    pub fn prescribed_displacement_scale_work(&self,
        prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3], v: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        if v.len()!=self.n() || !v.iter().all(|v| v.is_finite()) { return Err(ElasticityError3::Invalid("invalid displacement load pullback field")); }
        let cells = self.displacement_cells(prescribed, &mut checkpoint)?;
        let mut result = vec![0.0; self.cells.len()];
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            for i in 0..24 { if !self.fixed[cell.nodes[i/3]] {
                result[id] += cells[id][i]*v[3*cell.nodes[i/3]+i%3];
            } }
        }
        if !result.iter().all(|x| x.is_finite()) { return Err(ElasticityError3::Invalid("displacement load pullback overflow")); }
        poll(&mut checkpoint)?;
        Ok(result)
    }
}
