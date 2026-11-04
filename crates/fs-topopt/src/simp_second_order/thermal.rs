//! Exact raw-density curvature of prescribed-temperature equivalent-load work.
//! Thermal and stiffness loads use the same native design weight. Temperature,
//! material state and source-integrated free strain remain fixed during design.

use super::{DesignMap, robust_parameters};
use crate::elasticity::second_order::thermal::{ThermalComplianceGoal, ThermalLoadCase};
use crate::elasticity::second_order::{ComplianceError, ComplianceSolveOptions};
use crate::filter::DensityFilter;
use crate::pipeline::{LoadCase, SimpParams};
use crate::{DensityElasticity, DesignPipeline, RobustPipeline, SolveControl};

/// Thermal equivalent-load compliance, physical volume and source-bound mass
/// at one raw design. This is f(rho)^T u, not stored elastic strain energy.
/// The original thermal goal retains source states and physical stress/reactions.
/// Exact Hessian actions include load derivatives and the full design-map
/// curvature. Positive thermal sensitivities require a general constrained
/// optimizer; they must not enter the dead-load monotone OC update.
pub struct SimpThermalComplianceGoal<'a> {
    map: DesignMap<'a>,
    physical: ThermalComplianceGoal<'a>,
    gradient: Vec<f64>,
}
impl std::fmt::Debug for SimpThermalComplianceGoal<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SimpThermalComplianceGoal")
            .field("value", &self.value())
            .field("loads", &self.load_case_count())
            .finish_non_exhaustive()
    }
}

impl<'a> SimpThermalComplianceGoal<'a> {
    /// Ordinary projection, with source-bound thermal states retained per case.
    pub fn new(
        pipeline: &'a DesignPipeline,
        elasticity: &'a DensityElasticity,
        rho: &[f64],
        loads: &[ThermalLoadCase<'_>],
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

    /// Eroded stiffness AND equivalent thermal loads, with dilated material
    /// volume/mass. All fields share the same original Helmholtz filter.
    pub fn new_robust(
        pipeline: &'a RobustPipeline,
        elasticity: &'a DensityElasticity,
        rho: &[f64],
        loads: &[ThermalLoadCase<'_>],
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
        loads: &[ThermalLoadCase<'_>],
        options: ComplianceSolveOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        // Reuse identical force/weight admission before the filter. Thermal
        // identity, field extent and strain admission stay with the source owner.
        let forces: Vec<_> = loads
            .iter()
            .map(|load| LoadCase {
                force: load.force,
                weight: load.weight,
            })
            .collect();
        let (map, moduli) = DesignMap::prepare(
            filter, p, volume_eta, elasticity, rho, &forces, options, control,
        )?;
        let physical = ThermalComplianceGoal::new(elasticity, &moduli, loads, options, control)?;
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
    pub fn physical(&self) -> &ThermalComplianceGoal<'a> {
        &self.physical
    }
    pub fn projected_densities(&self) -> &[f64] {
        &self.map.projected
    }
    pub fn volume_projected_densities(&self) -> &[f64] {
        &self.map.volume_projected
    }
    pub fn load_case_count(&self) -> usize {
        self.physical.case_compliances().len()
    }

    fn load_case_value(&self, index: usize) -> Result<f64, ComplianceError> {
        self.physical
            .case_compliances()
            .get(index)
            .copied()
            .ok_or(ComplianceError::InvalidInput(
                "load-case index out of range",
            ))
    }

    /// Unweighted case derivative, including its original thermal load even
    /// when its aggregate objective weight is zero. No new primal solve runs.
    pub fn load_case_value_gradient(
        &self,
        index: usize,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let value = self.load_case_value(index)?;
        let modulus_gradient = self.physical.load_case_gradient(index, control)?;
        let gradient = self.map.pullback_gradient(
            &modulus_gradient,
            "simp-load-case-gradient",
            "SIMP load-case gradient",
            control,
        )?;
        Ok((value, gradient))
    }

    /// Exact raw case curvature through both f(rho) and K(rho).
    pub fn load_case_hessian_vector_product(
        &self,
        index: usize,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.map.direction(direction, control)?;
        if direction.iter().all(|&value| value == 0.0) {
            self.load_case_value(index)?;
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
        let gradient = self.physical.load_case_gradient(index, control)?;
        self.map.pullback_hessian(
            &filtered,
            &physical,
            &gradient,
            "simp-load-case-hessian-pullback",
            "SIMP load-case Hessian chain",
            control,
        )
    }

    /// One filter, one native thermal tangent per active case, one transpose.
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
