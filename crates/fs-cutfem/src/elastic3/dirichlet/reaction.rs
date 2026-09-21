//! Generalized numerical reaction on the retained embedded support.
use super::*;

/// A scalar support reaction and its two PARTIAL derivatives. No equilibrium
/// solve is performed and no residual, continuum flux, or actuator energy is
/// certified. The derivatives refer to the supplied field/current cell scales.
#[derive(Debug, Clone)]
pub struct EmbeddedReaction3 {
    /// Integral h dot [sigma_s(u)n - gamma_s(u-g)] on selected support points.
    /// Sign: traction exerted ON the solid. Constant h selects a force component;
    /// h=e_i cross (x-origin) selects a moment component about that origin.
    pub value: f64,
    /// Partial derivative in independent displacement coordinates at fixed s/g.
    /// Strong homogeneous clamp entries are zero; adaptive models return T^T q.
    pub displacement_gradient: Vec<f64>,
    /// Partial derivative per cell scale at fixed u/g; NOT the solved total
    /// derivative. Combine with the equilibrium adjoint and lifting derivative.
    pub scale_gradient: Vec<f64>,
}

impl CutElasticity3 {
    /// Evaluate the Nitsche-consistent embedded reaction, not raw sampled stress
    /// alone. With b_h the EXISTING prescribed-motion lifting for virtual mode h,
    /// R_h = -b_h(s)^T u + sum_c s_c integral gamma_c h dot g.
    /// Thus R_u=-b_h(s) and R_s[c]=integral gamma_c h dot g-b_h,c^T u.
    /// The same reference penalty, support selection, surface rules and trace
    /// assembly used by K and b_g define this functional; geometry is not rebuilt.
    ///
    /// `mode` and optional `prescribed` must be pure, finite reference laws,
    /// density-independent and compatible with any strong box clamps. None
    /// means g=0. Zeroing mode on part of the support selects a patch; narrow
    /// patches still require quadrature resolution. Only embedded-support
    /// reactions are included, never reactions at additional strong box clamps.
    ///
    /// For a solved scalar objective R, solve K z=R_u and use
    /// dR/ds_c=R_s[c]+z^T db_g/ds_c-z^T(dK/ds_c)u. This numerical reference
    /// reaction is not physical actuator work or a follower/shape derivative.
    pub fn embedded_reaction(&self, u: &[f64],
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        mode: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<EmbeddedReaction3, ElasticityError3> {
        poll(&mut checkpoint)?;
        if u.len() != self.n() || u.iter().any(|v| !v.is_finite()) {
            return Err(ElasticityError3::Invalid("invalid embedded reaction field"));
        }
        let local = self.displacement_cells(mode, &mut checkpoint)?;
        let data = self.embedded.as_ref().expect("lifting admitted embedded support");
        let mut result = EmbeddedReaction3 { value: 0.0, displacement_gradient: vec![0.0; self.n()],
            scale_gradient: vec![0.0; self.cells.len()] };
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let mut direct = 0.0;
            if let Some(g) = prescribed {
                if !data.points[id].is_empty() {
                    let gamma = penalty(cell, data.lame, data.beta)?;
                    let rule = cell.rules.surface().expect("builder retained surface");
                    for &index in &data.points[id] {
                        poll(&mut checkpoint)?;
                        let point = &rule.points()[index];
                        let h = mode(point.position, point.normal);
                        let value = g(point.position, point.normal);
                        poll(&mut checkpoint)?;
                        if h.iter().chain(&value).any(|v| !v.is_finite()) {
                            return Err(ElasticityError3::Invalid("nonfinite reaction mode or prescribed motion"));
                        }
                        direct += point.weight*gamma*h.iter().zip(value).map(|(a,b)| a*b).sum::<f64>();
                    }
                }
            }
            for i in 0..24 {
                let node = cell.nodes[i/3];
                if !self.fixed[node] {
                    let dof = 3*node+i%3;
                    direct -= local[id][i]*u[dof];
                    result.displacement_gradient[dof] -= self.scales[id]*local[id][i];
                }
            }
            result.scale_gradient[id] = direct;
            result.value += self.scales[id]*direct;
        }
        if !result.value.is_finite() || result.displacement_gradient.iter()
            .chain(&result.scale_gradient).any(|v| !v.is_finite()) {
            return Err(ElasticityError3::Invalid("embedded reaction arithmetic overflow"));
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }
}
