//! Complete matrix-free curvature through the existing density hygiene stack.
//!
//! For `t = F rho` and the projection/SIMP map `m(t)`,
//! `H_rho v = F^T [m' * H_m(m' * Fv) + m'' * Fv * g_m]`.
//! The ordinary constructor uses one projection for compliance and inventory.
//! The robust constructor uses eroded stiffness and dilated volume and mass,
//! sharing ONE filtered field. The two projections have distinct derivatives.
//!
//! Stiffness projected densities must lie strictly in `(1e-12, 1)`: no
//! generalized Hessian is invented at the existing SIMP clamp/floor kinks.

use crate::{DensityElasticity, DesignPipeline, SolveControl};
use crate::elasticity::second_order::{ComplianceError, ComplianceGoal, ComplianceSolveOptions};
use crate::filter::{DensityFilter, heaviside, heaviside_derivative, heaviside_second_derivative};
use crate::pipeline::{LoadCase, SimpParams};
use crate::robust::RobustPipeline;

/// Smooth compliance, volume and material-mass linearization at one raw design.
/// All filter, transpose-filter and tangent solves share the caller's work
/// controller. There is no dense filter or parameter Hessian.
pub struct SimpComplianceGoal<'a> {
    filter: &'a DensityFilter,
    physical: ComplianceGoal<'a>,
    reference_cell_masses_kg: Option<&'a [f64]>,
    projected: Vec<f64>,
    volume_projected: Vec<f64>,
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

/// Resolve the immutable eroded/dilated pair before any numerical work.
/// This is also used by the live-IR provider; it must not reinterpret delta.
pub(crate) fn robust_parameters(p: SimpParams, offset: f64)
    -> Result<(SimpParams, f64), ComplianceError>
{
    let eroded = p.eta + offset;
    let dilated = p.eta - offset;
    if !p.eta.is_finite() || !offset.is_finite() || offset < 0.0
        || !(0.0..=1.0).contains(&eroded) || !(0.0..=1.0).contains(&dilated) {
        return Err(ComplianceError::InvalidInput("robust thresholds must satisfy 0 <= eta-delta <= eta+delta <= 1"));
    }
    Ok((SimpParams { eta: eroded, ..p }, dilated))
}

impl<'a> SimpComplianceGoal<'a> {
    /// Ordinary SIMP: compliance and volume use the same projection.
    pub fn new(pipeline: &'a DesignPipeline, elasticity: &'a DensityElasticity, rho: &[f64],
        loads: &[LoadCase<'_>], options: ComplianceSolveOptions, control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        Self::prepare(&pipeline.filter, pipeline.params, pipeline.params.eta,
            elasticity, rho, loads, options, control)
    }

    /// Eroded compliance with dilated physical volume, using the existing
    /// robust thresholds `eta +/- eta_offset` and one native Helmholtz solve.
    /// The nominal volume is NOT silently substituted for the dilated budget.
    /// Independent loads retain their weights and homogeneous supports.
    ///
    /// This is a discrete three-field design model, not a certified minimum
    /// feature size, binary topology or manufacturing-yield guarantee.
    pub fn new_robust(pipeline: &'a RobustPipeline, elasticity: &'a DensityElasticity, rho: &[f64],
        loads: &[LoadCase<'_>], options: ComplianceSolveOptions, control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        control.checkpoint("robust-simp-admission")?;
        let (p, volume_eta) = robust_parameters(pipeline.params, pipeline.eta_offset)?;
        Self::prepare(&pipeline.filter, p, volume_eta, elasticity, rho, loads, options, control)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn prepare(
        filter: &'a DensityFilter,
        p: SimpParams,
        volume_eta: f64,
        elasticity: &'a DensityElasticity,
        rho: &[f64],
        loads: &[LoadCase<'_>],
        options: ComplianceSolveOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        control.checkpoint("simp-hessian-admission")?;
        let n = elasticity.cells();
        if n == 0 || rho.len() != n || filter.cells() != n {
            return Err(ComplianceError::InvalidInput("raw design, filter and elasticity shapes must agree"));
        }
        if !(p.e_min.is_finite() && p.e_min > 0.0 && p.e_min < 1.0
            && p.penal.is_finite() && p.penal >= 1.0 && p.beta.is_finite() && p.beta >= 0.0
            && p.eta.is_finite() && (0.0..=1.0).contains(&p.eta)
            && volume_eta.is_finite() && (0.0..=1.0).contains(&volume_eta)) {
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
        if loads.is_empty() {
            return Err(ComplianceError::InvalidInput("at least one structural load case required"));
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
        let filtered = filter.try_apply(rho, control)?;
        let mut projected = Vec::with_capacity(n); let mut volume_projected = Vec::with_capacity(n);
        let mut moduli = Vec::with_capacity(n);
        let mut first = Vec::with_capacity(n); let mut second = Vec::with_capacity(n);
        let mut projection_first = Vec::with_capacity(n); let mut projection_second = Vec::with_capacity(n);
        for (e, &rt) in filtered.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("simp-curvature")?; }
            let rb = finite(heaviside(rt, p.beta, p.eta), "stiffness projection")?;
            if !(rb > 1e-12 && rb < 1.0) {
                return Err(ComplianceError::InvalidInput("SIMP Hessian requires projected density strictly in (1e-12, 1)"));
            }
            projected.push(rb);
            // Same native projection and SIMP arithmetic as try_forward. Since
            // rb is strictly inside the smooth interval, clamp/floor are idle.
            moduli.push(finite(p.e_min + (1.0-p.e_min)*fs_math::det::pow(rb, p.penal), "SIMP modulus")?);
            let dh = finite(heaviside_derivative(rt, p.beta, p.eta), "projection slope")?;
            let ddh = finite(heaviside_second_derivative(rt, p.beta, p.eta), "projection curvature")?;
            let dm = (1.0-p.e_min)*p.penal*fs_math::det::pow(rb, p.penal-1.0);
            let ddm = if p.penal == 1.0 { 0.0 } else {
                (1.0-p.e_min)*p.penal*(p.penal-1.0)*fs_math::det::pow(rb, p.penal-2.0)
            };
            first.push(finite(dm*dh, "SIMP first chain")?);
            second.push(finite(ddm*dh*dh+dm*ddh, "SIMP second chain")?);
            let (vb, dv, ddv) = if volume_eta == p.eta { (rb, dh, ddh) } else {
                (finite(heaviside(rt, p.beta, volume_eta), "volume projection")?,
                 finite(heaviside_derivative(rt, p.beta, volume_eta), "volume projection slope")?,
                 finite(heaviside_second_derivative(rt, p.beta, volume_eta), "volume projection curvature")?)
            };
            if !(0.0..=1.0).contains(&vb) {
                return Err(ComplianceError::InvalidInput("volume projection must lie in [0,1]"));
            }
            volume_projected.push(vb);
            projection_first.push(dv); projection_second.push(ddv);
        }
        let physical = ComplianceGoal::new(elasticity, &moduli, loads, options, control)?;
        let mut local = Vec::with_capacity(n);
        for e in 0..n {
            if e % 256 == 0 { control.checkpoint("simp-gradient")?; }
            local.push(finite(first[e]*physical.gradient()[e], "SIMP gradient")?);
        }
        let gradient = filter.try_apply_transpose(&local, control)?;
        control.checkpoint("simp-linearization-publish")?;
        Ok(Self {
            filter,
            physical,
            reference_cell_masses_kg: elasticity.reference_cell_masses_kg(),
            projected,
            volume_projected,
            first,
            second,
            projection_first,
            projection_second,
            gradient,
        })
    }

    pub fn value(&self) -> f64 { self.physical.value() }
    pub fn gradient(&self) -> &[f64] { &self.gradient }
    pub fn physical(&self) -> &ComplianceGoal<'a> { &self.physical }
    /// Projection used by stiffness (eroded for a robust goal).
    pub fn projected_densities(&self) -> &[f64] { &self.projected }
    /// Projection used by the volume constraint (dilated for a robust goal).
    pub fn volume_projected_densities(&self) -> &[f64] { &self.volume_projected }
    /// Number of independent structural load cases retained by the native goal.
    pub fn load_case_count(&self) -> usize { self.physical.load_case_count() }

    /// Raw-density value/gradient for one original load case. This reuses the
    /// retained primal displacement and spends only the transpose filter needed
    /// to pull that case's modulus gradient back to raw density coordinates.
    pub fn load_case_value_gradient(&self, index: usize, control: &mut SolveControl<'_>)
        -> Result<(f64, Vec<f64>), ComplianceError>
    {
        let value = self.physical.load_case_value(index)?;
        let modulus_gradient = self.physical.load_case_gradient(index, control)?;
        let mut local = Vec::with_capacity(self.first.len());
        for (element, (&chain, &gradient)) in self.first.iter().zip(&modulus_gradient).enumerate() {
            if element % 256 == 0 { control.checkpoint("simp-load-case-gradient")?; }
            local.push(finite(chain * gradient, "SIMP load-case gradient")?);
        }
        Ok((value, self.filter.try_apply_transpose(&local, control)?))
    }

    /// Complete raw-density Hessian action for one original load case:
    /// F^T[m' H_case(m' Fv) + m'' Fv g_case]. One nonzero action spends one
    /// forward filter, one physical tangent solve, and one transpose filter.
    pub fn load_case_hessian_vector_product(&self, index: usize, direction: &[f64],
        control: &mut SolveControl<'_>) -> Result<Vec<f64>, ComplianceError>
    {
        self.direction(direction, control)?;
        if direction.iter().all(|&value| value == 0.0) {
            // Still validate the case index; zero work must not turn a bad API
            // request into an apparently valid product.
            self.physical.load_case_value(index)?;
            return Ok(vec![0.0; direction.len()]);
        }
        let filtered = self.filter.try_apply(direction, control)?;
        let mut tangent = Vec::with_capacity(direction.len());
        for (element, &value) in filtered.iter().enumerate() {
            if element % 256 == 0 { control.checkpoint("simp-load-case-tangent")?; }
            tangent.push(finite(self.first[element] * value, "SIMP load-case modulus tangent")?);
        }
        let physical = self.physical.load_case_hessian_vector_product(index, &tangent, control)?;
        let modulus_gradient = self.physical.load_case_gradient(index, control)?;
        let mut local = Vec::with_capacity(direction.len());
        for element in 0..direction.len() {
            if element % 256 == 0 { control.checkpoint("simp-load-case-hessian-pullback")?; }
            local.push(finite(self.first[element] * physical[element]
                + self.second[element] * filtered[element] * modulus_gradient[element],
                "SIMP load-case Hessian chain")?);
        }
        Ok(self.filter.try_apply_transpose(&local, control)?)
    }

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
        let filtered = self.filter.try_apply(direction, control)?;
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
        Ok(self.filter.try_apply_transpose(&local, control)?)
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

    /// Physical projected volume. For robust goals this is DILATED volume,
    /// not eroded stiffness, nominal material, or an unweighted cell mean.
    pub fn volume_value_gradient(
        &self,
        volumes: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let weights = self.weights(volumes, control)?;
        self.projected_quantity_value_gradient(&weights, "simp-volume-gradient", control)
    }

    /// Material inventory in kilograms, with its exact raw-occupancy gradient.
    /// Uses the stiffness owner's reference cell masses, without normalization
    /// or modal void-mass interpolation. Robust goals use DILATED occupancy.
    pub fn mass_value_gradient(
        &self,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let masses = self
            .reference_cell_masses_kg
            .ok_or(ComplianceError::InvalidInput(
                "physical mass requires an assigned material field",
            ))?;
        self.projected_quantity_value_gradient(masses, "simp-mass-gradient", control)
    }

    fn projected_quantity_value_gradient(
        &self,
        weights: &[f64],
        stage: &'static str,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let mut value = 0.0; let mut local = Vec::with_capacity(weights.len());
        for (e, &w) in weights.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint(stage)?; }
            value += w*self.volume_projected[e];
            local.push(finite(w*self.projection_first[e], stage)?);
        }
        let value = finite(value, stage)?;
        Ok((value, self.filter.try_apply_transpose(&local, control)?))
    }

    /// F^T diag(V/sum(V) * H_volume'') F v; no SIMP term belongs in volume.
    pub fn volume_hessian_vector_product(
        &self,
        volumes: &[f64],
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.direction(direction, control)?;
        let weights = self.weights(volumes, control)?;
        self.projected_quantity_hessian_vector_product(
            &weights,
            direction,
            "simp-volume-hessian",
            control,
        )
    }

    /// F^T diag(m_ref * H_inventory'') F v, in kilograms. Includes the same
    /// robust dilation as the mass value; no stiffness or modal mass derivative.
    pub fn mass_hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.direction(direction, control)?;
        let masses = self
            .reference_cell_masses_kg
            .ok_or(ComplianceError::InvalidInput(
                "physical mass requires an assigned material field",
            ))?;
        self.projected_quantity_hessian_vector_product(
            masses,
            direction,
            "simp-mass-hessian",
            control,
        )
    }

    fn projected_quantity_hessian_vector_product(
        &self,
        weights: &[f64],
        direction: &[f64],
        stage: &'static str,
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        if direction.iter().all(|&v| v == 0.0) || self.projection_second.iter().all(|&h| h == 0.0) {
            return Ok(vec![0.0; direction.len()]);
        }
        let filtered = self.filter.try_apply(direction, control)?;
        let mut local = Vec::with_capacity(direction.len());
        for e in 0..direction.len() {
            if e % 256 == 0 { control.checkpoint(stage)?; }
            local.push(finite(weights[e]*self.projection_second[e]*filtered[e], stage)?);
        }
        Ok(self.filter.try_apply_transpose(&local, control)?)
    }
}
