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

use crate::elasticity::second_order::{ComplianceError, ComplianceGoal, ComplianceSolveOptions};
use crate::filter::DensityFilter;
use crate::pipeline::{LoadCase, SimpParams};
use crate::robust::RobustPipeline;
use crate::{DensityElasticity, DesignPipeline, SolveControl};

mod design_map;
use design_map::DesignMap;
mod thermal;
pub use thermal::SimpThermalComplianceGoal;

/// Smooth compliance, volume and material-mass linearization at one raw design.
/// All filter, transpose-filter and tangent solves share the caller's work
/// controller. There is no dense filter or parameter Hessian.
pub struct SimpComplianceGoal<'a> {
    map: DesignMap<'a>,
    physical: ComplianceGoal<'a>,
    gradient: Vec<f64>,
}
impl std::fmt::Debug for SimpComplianceGoal<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SimpComplianceGoal")
            .field("physical", &self.physical)
            .finish_non_exhaustive()
    }
}
fn finite(value: f64, stage: &'static str) -> Result<f64, ComplianceError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ComplianceError::NonFinite(stage))
    }
}

/// Resolve the immutable eroded/dilated pair before any numerical work.
/// This is also used by the live-IR provider; it must not reinterpret delta.
pub(crate) fn robust_parameters(
    p: SimpParams,
    offset: f64,
) -> Result<(SimpParams, f64), ComplianceError> {
    let eroded = p.eta + offset;
    let dilated = p.eta - offset;
    if !p.eta.is_finite()
        || !offset.is_finite()
        || offset < 0.0
        || !(0.0..=1.0).contains(&eroded)
        || !(0.0..=1.0).contains(&dilated)
    {
        return Err(ComplianceError::InvalidInput(
            "robust thresholds must satisfy 0 <= eta-delta <= eta+delta <= 1",
        ));
    }
    Ok((SimpParams { eta: eroded, ..p }, dilated))
}

impl<'a> SimpComplianceGoal<'a> {
    /// Ordinary SIMP: compliance and volume use the same projection.
    pub fn new(
        pipeline: &'a DesignPipeline,
        elasticity: &'a DensityElasticity,
        rho: &[f64],
        loads: &[LoadCase<'_>],
        options: ComplianceSolveOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        Self::prepare(
            &pipeline.filter,
            pipeline.params,
            pipeline.params.eta,
            elasticity,
            rho,
            loads,
            options,
            control,
        )
    }

    /// Eroded compliance with dilated physical volume, using the existing
    /// robust thresholds `eta +/- eta_offset` and one native Helmholtz solve.
    /// The nominal volume is NOT silently substituted for the dilated budget.
    /// Independent loads retain their weights and homogeneous supports.
    ///
    /// This is a discrete three-field design model, not a certified minimum
    /// feature size, binary topology or manufacturing-yield guarantee.
    pub fn new_robust(
        pipeline: &'a RobustPipeline,
        elasticity: &'a DensityElasticity,
        rho: &[f64],
        loads: &[LoadCase<'_>],
        options: ComplianceSolveOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        control.checkpoint("robust-simp-admission")?;
        let (p, volume_eta) = robust_parameters(pipeline.params, pipeline.eta_offset)?;
        Self::prepare(
            &pipeline.filter,
            p,
            volume_eta,
            elasticity,
            rho,
            loads,
            options,
            control,
        )
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
        let (map, moduli) = DesignMap::prepare(
            filter, p, volume_eta, elasticity, rho, loads, options, control,
        )?;
        let physical = ComplianceGoal::new(elasticity, &moduli, loads, options, control)?;
        let gradient = map.pullback_gradient(
            physical.gradient(),
            "simp-gradient",
            "SIMP gradient",
            control,
        )?;
        control.checkpoint("simp-linearization-publish")?;
        Ok(Self {
            map,
            physical,
            gradient,
        })
    }

    pub fn value(&self) -> f64 {
        self.physical.value()
    }
    pub fn gradient(&self) -> &[f64] {
        &self.gradient
    }
    pub fn physical(&self) -> &ComplianceGoal<'a> {
        &self.physical
    }
    /// Projection used by stiffness (eroded for a robust goal).
    pub fn projected_densities(&self) -> &[f64] {
        &self.map.projected
    }
    /// Projection used by volume and mass (dilated for a robust goal).
    pub fn volume_projected_densities(&self) -> &[f64] {
        &self.map.volume_projected
    }
    pub fn load_case_count(&self) -> usize {
        self.physical.load_case_count()
    }

    /// Unweighted raw-density value/gradient of one original load case.
    /// Retains its primal and spends only one transpose-filter solve.
    pub fn load_case_value_gradient(
        &self,
        index: usize,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let value = self.physical.load_case_value(index)?;
        let modulus_gradient = self.physical.load_case_gradient(index, control)?;
        let gradient = self.map.pullback_gradient(
            &modulus_gradient,
            "simp-load-case-gradient",
            "SIMP load-case gradient",
            control,
        )?;
        Ok((value, gradient))
    }

    /// F^T[m' H_case(m' Fv) + m'' Fv g_case], with one native tangent solve.
    pub fn load_case_hessian_vector_product(
        &self,
        index: usize,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.map.direction(direction, control)?;
        if direction.iter().all(|&value| value == 0.0) {
            self.physical.load_case_value(index)?;
            return Ok(vec![0.0; direction.len()]);
        }
        let (filtered, tangent) = self.map.tangent(
            direction,
            "simp-load-case-tangent",
            "SIMP load-case modulus tangent",
            control,
        )?;
        let physical = self
            .physical
            .load_case_hessian_vector_product(index, &tangent, control)?;
        let modulus_gradient = self.physical.load_case_gradient(index, control)?;
        self.map.pullback_hessian(
            &filtered,
            &physical,
            &modulus_gradient,
            "simp-load-case-hessian-pullback",
            "SIMP load-case Hessian chain",
            control,
        )
    }

    /// One forward filter, one tangent per active load, one transpose filter.
    pub fn hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.map.direction(direction, control)?;
        if direction.iter().all(|&v| v == 0.0) {
            return Ok(vec![0.0; direction.len()]);
        }
        let (filtered, tangent) = self.map.tangent(
            direction,
            "simp-modulus-tangent",
            "SIMP modulus tangent",
            control,
        )?;
        let physical = self.physical.hessian_vector_product(&tangent, control)?;
        self.map.pullback_hessian(
            &filtered,
            &physical,
            self.physical.gradient(),
            "simp-hessian-pullback",
            "SIMP Hessian chain",
            control,
        )
    }

    /// Physical projected volume; robust goals use the dilated field.
    pub fn volume_value_gradient(
        &self,
        volumes: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        self.map.volume_value_gradient(volumes, control)
    }
    /// Source-bound material inventory in kilograms, using the inventory projection.
    pub fn mass_value_gradient(
        &self,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        self.map.mass_value_gradient(control)
    }
    /// Exact inventory curvature; no stiffness SIMP term belongs in volume.
    pub fn volume_hessian_vector_product(
        &self,
        volumes: &[f64],
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.map
            .volume_hessian_vector_product(volumes, direction, control)
    }
    /// Exact kilogram inventory curvature, including robust dilation.
    pub fn mass_hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.map.mass_hessian_vector_product(direction, control)
    }
}
