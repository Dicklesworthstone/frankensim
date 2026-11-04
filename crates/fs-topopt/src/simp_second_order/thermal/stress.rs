//! Physical thermal stress through the original filter/projection/SIMP map.
//! Robust goals retain eroded stiffness and stress alongside dilated inventory.
//! The physical sampled stress moment includes every independent load case.

use super::SimpThermalComplianceGoal;
use crate::elasticity::second_order::ComplianceError;
use crate::elasticity::second_order::thermal::stress::{
    StressLinearization, ThermalStressOptions,
};
use crate::SolveControl;

// Detached only inside the existing Newton cache, which binds the original
// physical owner and exact raw-design bits. Public callers retain the borrow.
pub(crate) struct SimpStressLinearization {
    pub(crate) physical: StressLinearization,
    pub(crate) gradient: Vec<f64>,
}

/// Source-bound physical stress and exact raw-density derivatives at one design.
/// Its source owns the original thermal/material fields, equilibrium and map.
pub struct SimpThermalStressGoal<'g, 'a> {
    source: &'g SimpThermalComplianceGoal<'a>,
    linearization: SimpStressLinearization,
}
impl<'g, 'a> SimpThermalStressGoal<'g, 'a> {
    pub fn value(&self) -> f64 {
        self.linearization.physical.value
    }
    pub fn gradient(&self) -> &[f64] {
        &self.linearization.gradient
    }
    pub fn maximum_ratio(&self) -> f64 {
        self.linearization.physical.maximum_ratio
    }
    pub fn options(&self) -> ThermalStressOptions {
        self.linearization.physical.options
    }
    pub fn source(&self) -> &'g SimpThermalComplianceGoal<'a> {
        self.source
    }
    /// Full Hessian chain, including design-map curvature. Uses one forward
    /// filter, two native tangent solves per case, and one transpose filter.
    pub fn hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.source.stress_hessian(&self.linearization, direction, control)
    }
}

impl<'a> SimpThermalComplianceGoal<'a> {
    /// Prepare the physical moment with one original adjoint per case and one
    /// transpose filter. The original primal fields are reused without solving.
    pub fn stress_goal<'g>(
        &'g self,
        options: ThermalStressOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<SimpThermalStressGoal<'g, 'a>, ComplianceError> {
        let linearization = self.prepare_stress(options, control)?;
        Ok(SimpThermalStressGoal { source: self, linearization })
    }

    pub(crate) fn prepare_stress(
        &self,
        options: ThermalStressOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<SimpStressLinearization, ComplianceError> {
        let physical = self.physical.prepare_stress(options, control)?;
        let gradient = self.map.pullback_gradient(
            &physical.gradient,
            "simp-thermal-stress-gradient",
            "SIMP thermal stress gradient",
            control,
        )?;
        control.checkpoint("simp-thermal-stress-publish")?;
        Ok(SimpStressLinearization { physical, gradient })
    }

    pub(crate) fn stress_hessian(
        &self,
        point: &SimpStressLinearization,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.map.direction(direction, control)?;
        if direction.iter().all(|&v| v == 0.0) {
            return Ok(vec![0.0; direction.len()]);
        }
        let (filtered, tangent) = self.map.tangent(
            direction,
            "simp-thermal-stress-tangent",
            "SIMP thermal stress tangent",
            control,
        )?;
        let physical = self.physical.stress_hessian(&point.physical, &tangent, control)?;
        self.map.pullback_hessian(
            &filtered,
            &physical,
            &point.physical.gradient,
            "simp-thermal-stress-hessian",
            "SIMP thermal stress Hessian chain",
            control,
        )
    }
}
