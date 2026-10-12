//! Complete matrix-free curvature through the existing density hygiene stack.
//!
//! Let `F` be the native Helmholtz cell filter and `m(t)` the native
//! projection/SIMP composition. For `t = F rho`,
//! `H_rho v = F^T [m' * H_m(m' * Fv) + m'' * Fv * g_m]`.
//! The second term is essential; pulling the stiffness Hessian through a
//! first-order Jacobian alone does not differentiate a nonlinear design map.
//!
//! This linearization deliberately refuses the stiffness clamp/floor region:
//! every projected density must lie strictly in `(1e-12, 1)`. No generalized
//! Hessian is invented at those kinks. First-order pipeline behavior is unchanged.

use crate::{DensityElasticity, DesignPipeline, SolveControl};
use crate::elasticity::second_order::{ComplianceError, ComplianceGoal, ComplianceSolveOptions};
use crate::filter::{heaviside_derivative, heaviside_second_derivative};
use crate::pipeline::LoadCase;

/// A retained smooth raw-density compliance and projected-volume linearization.
/// All filter, transpose-filter and physical tangent solves use the caller's
/// shared lifetime controller. No new linear solver or dense filter is built.
pub struct SimpComplianceGoal<'a> {
    pipeline: &'a DesignPipeline,
    physical: ComplianceGoal<'a>,
    projected: Vec<f64>,
    first: Vec<f64>,
    second: Vec<f64>,
    projection_first: Vec<f64>,
    projection_second: Vec<f64>,
    gradient: Vec<f64>,
}
impl std::fmt::Debug for SimpComplianceGoal<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SimpComplianceGoal").field("physical", &self.physical).finish_non_exhaustive()
    }
}
fn finite(value: f64, stage: &'static str) -> Result<f64, ComplianceError> {
    if value.is_finite() { Ok(value) } else { Err(ComplianceError::NonFinite(stage)) }
}

impl<'a> SimpComplianceGoal<'a> {
    pub fn new(pipeline: &'a DesignPipeline, elasticity: &'a DensityElasticity, rho: &[f64],
        loads: &[LoadCase<'_>], options: ComplianceSolveOptions, control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        control.checkpoint("simp-hessian-admission")?;
        let n = elasticity.cells(); let p = pipeline.params;
        if n == 0 || rho.len() != n || pipeline.filter.cells() != n {
            return Err(ComplianceError::InvalidInput("raw design, filter and elasticity shapes must agree"));
        }
        if !(p.e_min.is_finite() && p.e_min > 0.0 && p.e_min < 1.0
            && p.penal.is_finite() && p.penal >= 1.0 && p.beta.is_finite() && p.beta >= 0.0
            && p.eta.is_finite() && (0.0..=1.0).contains(&p.eta)) {
            return Err(ComplianceError::InvalidInput("invalid SIMP/projection parameters"));
        }
        if !(options.tolerance.is_finite() && options.tolerance > 0.0 && options.tolerance < 1.0) {
            return Err(ComplianceError::InvalidInput("relative tolerance must lie in (0, 1)"));
        }
        for (e, &r) in rho.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("simp-hessian-admission")?; }
            if !r.is_finite() || !(0.0..=1.0).contains(&r) {
                return Err(ComplianceError::InvalidInput("raw densities must lie in [0, 1]"));
            }
        }
        if loads.is_empty() || !loads.iter().any(|l| l.weight > 0.0) {
            return Err(ComplianceError::InvalidInput("at least one positive-weight load required"));
        }
        for load in loads {
            if load.force.len() != elasticity.n() || !(load.weight.is_finite() && load.weight >= 0.0) {
                return Err(ComplianceError::InvalidInput("invalid structural load shape or weight"));
            }
            for (d, &f) in load.force.iter().enumerate() {
                if d % 256 == 0 { control.checkpoint("simp-hessian-admission")?; }
                finite(f, "structural force")?;
            }
        }
        let (filtered, projected, moduli) = pipeline.try_forward(rho, control)?;
        let mut first = Vec::with_capacity(n); let mut second = Vec::with_capacity(n);
        let mut projection_first = Vec::with_capacity(n); let mut projection_second = Vec::with_capacity(n);
        for e in 0..n {
            if e % 256 == 0 { control.checkpoint("simp-curvature")?; }
            let rb = projected[e];
            if !(rb > 1e-12 && rb < 1.0) {
                return Err(ComplianceError::InvalidInput("SIMP Hessian requires projected density strictly in (1e-12, 1)"));
            }
            let dh = finite(heaviside_derivative(filtered[e], p.beta, p.eta), "projection slope")?;
            let ddh = finite(heaviside_second_derivative(filtered[e], p.beta, p.eta), "projection curvature")?;
            let dm = (1.0-p.e_min)*p.penal*fs_math::det::pow(rb, p.penal-1.0);
            let ddm = if p.penal == 1.0 { 0.0 } else {
                (1.0-p.e_min)*p.penal*(p.penal-1.0)*fs_math::det::pow(rb, p.penal-2.0)
            };
            first.push(finite(dm*dh, "SIMP first chain")?);
            second.push(finite(ddm*dh*dh+dm*ddh, "SIMP second chain")?);
            projection_first.push(dh); projection_second.push(ddh);
        }
        let physical = ComplianceGoal::new(elasticity, &moduli, loads, options, control)?;
        let mut local = Vec::with_capacity(n);
        for e in 0..n {
            if e % 256 == 0 { control.checkpoint("simp-gradient")?; }
            local.push(finite(first[e]*physical.gradient()[e], "SIMP gradient")?);
        }
        let gradient = pipeline.filter.try_apply_transpose(&local, control)?;
        control.checkpoint("simp-linearization-publish")?;
        Ok(Self { pipeline, physical, projected, first, second, projection_first, projection_second, gradient })
    }

    pub fn value(&self) -> f64 { self.physical.value() }
    pub fn gradient(&self) -> &[f64] { &self.gradient }
    pub fn physical(&self) -> &ComplianceGoal<'a> { &self.physical }
    pub fn projected_densities(&self) -> &[f64] { &self.projected }

    fn direction(&self, direction: &[f64], control: &mut SolveControl<'_>) -> Result<(), ComplianceError> {
        control.checkpoint("simp-direction")?;
        if direction.len() != self.first.len() {
            return Err(ComplianceError::InvalidInput("raw-density Hessian direction shape mismatch"));
        }
        for (e, &v) in direction.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("simp-direction")?; }
            finite(v, "raw-density direction")?;
        }
        Ok(())
    }

    /// One forward filter, one tangent per active load, one transpose filter.
    pub fn hessian_vector_product(&self, direction: &[f64], control: &mut SolveControl<'_>)
        -> Result<Vec<f64>, ComplianceError>
    {
        self.direction(direction, control)?;
        if direction.iter().all(|&v| v == 0.0) { return Ok(vec![0.0; direction.len()]); }
        let filtered = self.pipeline.filter.try_apply(direction, control)?;
        let mut tangent = Vec::with_capacity(direction.len());
        for (e, &v) in filtered.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("simp-modulus-tangent")?; }
            tangent.push(finite(self.first[e]*v, "SIMP modulus tangent")?);
        }
        let physical = self.physical.hessian_vector_product(&tangent, control)?;
        let mut local = Vec::with_capacity(direction.len());
        for e in 0..direction.len() {
            if e % 256 == 0 { control.checkpoint("simp-hessian-pullback")?; }
            local.push(finite(self.first[e]*physical[e]
                +self.second[e]*filtered[e]*self.physical.gradient()[e], "SIMP Hessian chain")?);
        }
        Ok(self.pipeline.filter.try_apply_transpose(&local, control)?)
    }

    fn weights(&self, volumes: &[f64], control: &mut SolveControl<'_>) -> Result<Vec<f64>, ComplianceError> {
        control.checkpoint("simp-volume-admission")?;
        if volumes.len() != self.first.len() {
            return Err(ComplianceError::InvalidInput("one physical volume per design cell required"));
        }
        let mut total = 0.0;
        for (e, &volume) in volumes.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("simp-volume-admission")?; }
            if !(volume.is_finite() && volume > 0.0) {
                return Err(ComplianceError::InvalidInput("cell volumes must be finite and positive"));
            }
            total = finite(total+volume, "total volume")?;
        }
        Ok(volumes.iter().map(|volume| volume/total).collect())
    }

    /// Projected physical volume, not mean stiffness or the raw volume vector.
    pub fn volume_value_gradient(&self, volumes: &[f64], control: &mut SolveControl<'_>)
        -> Result<(f64, Vec<f64>), ComplianceError>
    {
        let weights = self.weights(volumes, control)?;
        let mut value = 0.0; let mut local = Vec::with_capacity(weights.len());
        for (e, &w) in weights.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("simp-volume-gradient")?; }
            value += w*self.projected[e];
            local.push(finite(w*self.projection_first[e], "projected volume gradient")?);
        }
        let value = finite(value, "projected volume")?;
        Ok((value, self.pipeline.filter.try_apply_transpose(&local, control)?))
    }

    /// F^T diag(V/sum(V) * H'') F v: there is NO SIMP term in volume.
    pub fn volume_hessian_vector_product(&self, volumes: &[f64], direction: &[f64], control: &mut SolveControl<'_>)
        -> Result<Vec<f64>, ComplianceError>
    {
        self.direction(direction, control)?;
        let weights = self.weights(volumes, control)?;
        if direction.iter().all(|&v| v == 0.0) || self.projection_second.iter().all(|&h| h == 0.0) {
            return Ok(vec![0.0; direction.len()]);
        }
        let filtered = self.pipeline.filter.try_apply(direction, control)?;
        let mut local = Vec::with_capacity(direction.len());
        for e in 0..direction.len() {
            if e % 256 == 0 { control.checkpoint("simp-volume-hessian")?; }
            local.push(finite(weights[e]*self.projection_second[e]*filtered[e], "projected volume Hessian")?);
        }
        Ok(self.pipeline.filter.try_apply_transpose(&local, control)?)
    }
}
