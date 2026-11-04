//! Native thermal goal selection inside the original cache and work owner.

use super::*;
use crate::elasticity::second_order::thermal::ThermalLoadCase;
use crate::simp_second_order::SimpThermalComplianceGoal;

pub(super) enum CachedGoal<'a> {
    Dead(SimpComplianceGoal<'a>),
    Thermal(SimpThermalComplianceGoal<'a>),
}
impl CachedGoal<'_> {
    pub(super) fn value(&self) -> f64 {
        match self {
            Self::Dead(g) => g.value(),
            Self::Thermal(g) => g.value(),
        }
    }
    pub(super) fn gradient(&self) -> &[f64] {
        match self {
            Self::Dead(g) => g.gradient(),
            Self::Thermal(g) => g.gradient(),
        }
    }
    pub(super) fn volume_value_gradient(
        &self,
        volumes: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        match self {
            Self::Dead(g) => g.volume_value_gradient(volumes, control),
            Self::Thermal(g) => g.volume_value_gradient(volumes, control),
        }
    }
    pub(super) fn mass_value_gradient(
        &self,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        match self {
            Self::Dead(g) => g.mass_value_gradient(control),
            Self::Thermal(g) => g.mass_value_gradient(control),
        }
    }
    pub(super) fn load_case_value_gradient(
        &self,
        index: usize,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        match self {
            Self::Dead(g) => g.load_case_value_gradient(index, control),
            Self::Thermal(g) => g.load_case_value_gradient(index, control),
        }
    }
    pub(super) fn hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        match self {
            Self::Dead(g) => g.hessian_vector_product(direction, control),
            Self::Thermal(g) => g.hessian_vector_product(direction, control),
        }
    }
    pub(super) fn volume_hessian_vector_product(
        &self,
        volumes: &[f64],
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        match self {
            Self::Dead(g) => g.volume_hessian_vector_product(volumes, direction, control),
            Self::Thermal(g) => g.volume_hessian_vector_product(volumes, direction, control),
        }
    }
    pub(super) fn mass_hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        match self {
            Self::Dead(g) => g.mass_hessian_vector_product(direction, control),
            Self::Thermal(g) => g.mass_hessian_vector_product(direction, control),
        }
    }
    pub(super) fn load_case_hessian_vector_product(
        &self,
        index: usize,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        match self {
            Self::Dead(g) => g.load_case_hessian_vector_product(index, direction, control),
            Self::Thermal(g) => g.load_case_hessian_vector_product(index, direction, control),
        }
    }
}

impl<'a> SimpNewtonModel<'a> {
    /// Prescribed-temperature equivalent-load compliance with original volume
    /// and mass views. Forces and thermal fields are copied once, preserving
    /// their source identities; no numerical solve occurs during construction.
    /// Heat conduction, temperatures and material parameters are not design variables.
    pub fn new_thermal(
        pipeline: &'a DesignPipeline,
        elasticity: &'a DensityElasticity,
        loads: &[ThermalLoadCase<'_>],
        volumes: &[f64],
        options: ComplianceSolveOptions,
        budget: ComplianceNewtonBudget,
    ) -> Result<Self, PhysicsError> {
        let forces: Vec<_> = loads
            .iter()
            .map(|load| LoadCase {
                force: load.force,
                weight: load.weight,
            })
            .collect();
        Self::new(pipeline, elasticity, &forces, volumes, options, budget)?
            .with_thermal_fields(loads)
    }

    /// Robust thermal design: eroded stiffness/equivalent loads and dilated
    /// volume AND mass, on the same filter and shared native resource account.
    pub fn new_thermal_robust(
        pipeline: &'a RobustPipeline,
        elasticity: &'a DensityElasticity,
        loads: &[ThermalLoadCase<'_>],
        volumes: &[f64],
        options: ComplianceSolveOptions,
        budget: ComplianceNewtonBudget,
    ) -> Result<Self, PhysicsError> {
        let forces: Vec<_> = loads
            .iter()
            .map(|load| LoadCase {
                force: load.force,
                weight: load.weight,
            })
            .collect();
        Self::new_robust(pipeline, elasticity, &forces, volumes, options, budget)?
            .with_thermal_fields(loads)
    }

    /// Ordinary thermal stage sharing the original lifetime work ceilings.
    pub fn new_thermal_with_resources(
        pipeline: &'a DesignPipeline,
        elasticity: &'a DensityElasticity,
        loads: &[ThermalLoadCase<'_>],
        volumes: &[f64],
        options: ComplianceSolveOptions,
        resources: SimpNewtonResources,
    ) -> Result<Self, PhysicsError> {
        let mut model = Self::new_thermal(
            pipeline,
            elasticity,
            loads,
            volumes,
            options,
            resources.budget(),
        )?;
        model.resources = resources;
        Ok(model)
    }

    /// Robust thermal stage; sharing counters never shares another stage's cache.
    pub fn new_thermal_robust_with_resources(
        pipeline: &'a RobustPipeline,
        elasticity: &'a DensityElasticity,
        loads: &[ThermalLoadCase<'_>],
        volumes: &[f64],
        options: ComplianceSolveOptions,
        resources: SimpNewtonResources,
    ) -> Result<Self, PhysicsError> {
        let mut model = Self::new_thermal_robust(
            pipeline,
            elasticity,
            loads,
            volumes,
            options,
            resources.budget(),
        )?;
        model.resources = resources;
        Ok(model)
    }

    fn with_thermal_fields(mut self, loads: &[ThermalLoadCase<'_>]) -> Result<Self, PhysicsError> {
        if self.elasticity.material_field().is_none() {
            return Err(PhysicsError::InvalidInput(
                "thermal design requires an assigned native elastic material field",
            ));
        }
        let fields = loads
            .iter()
            .map(|load| {
                OwnedThermalField::new(load.strains, self.elasticity.cells())
                    .map_err(|error| self.native_error(error))
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.thermal = Some(fields);
        Ok(self)
    }

    pub(super) fn prepare_goal(
        &self,
        point: &[f64],
        forces: &[LoadCase<'_>],
        control: &mut SolveControl<'_>,
    ) -> Result<CachedGoal<'a>, ComplianceError> {
        if let Some(fields) = &self.thermal {
            let loads: Vec<_> = forces
                .iter()
                .zip(fields)
                .map(|(load, field)| ThermalLoadCase {
                    force: load.force,
                    weight: load.weight,
                    strains: field.borrowed(),
                })
                .collect();
            Ok(CachedGoal::Thermal(SimpThermalComplianceGoal::prepare(
                self.filter,
                self.params,
                self.volume_eta,
                self.elasticity,
                point,
                &loads,
                self.options,
                control,
            )?))
        } else {
            Ok(CachedGoal::Dead(SimpComplianceGoal::prepare(
                self.filter,
                self.params,
                self.volume_eta,
                self.elasticity,
                point,
                forces,
                self.options,
                control,
            )?))
        }
    }
}
