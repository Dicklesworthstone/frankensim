//! Physical thermal stress views of the original shared Newton owner.
//! One optional stress entry per exact raw-design point retains all adjoints.
//! Switching a criterion replaces that entry only after successful preparation.

use super::*;
use crate::simp_second_order::thermal::SimpThermalComplianceGoal;

impl<'a> SimpNewtonModel<'a> {
    /// Dimensionless physical von Mises moment, for an objective or constraint
    /// J <= 1. Uses all original thermal cases, including zero compliance weights.
    /// The allowance is an explicit caller-supplied Pa criterion, not a strength
    /// inferred from an anisotropic elastic tensor. It is a finite-cell stress
    /// measure, not a continuum safety certificate or a microstructure recovery.
    ///
    /// The original equilibrium, source-domain admission, design map, mass and
    /// volume remain shared with the other views. Preparation adds one adjoint
    /// per case and one transpose filter. Repeated reads of the same criterion
    /// reuse the retained gradient/adjoints; an HVP adds two solves per case and
    /// two filters. Criterion changes cannot reset the lifetime work account.
    pub fn thermal_stress_moment(
        &self,
        options: ThermalStressOptions,
        study: impl Into<String>,
    ) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        if self.thermal.is_none() {
            return Err(PhysicsError::InvalidInput(
                "thermal stress requires the original source-bound thermal model",
            ));
        }
        options.admit().map_err(|error| self.native_error(error))?;
        self.loads.len().checked_mul(2).and_then(|n| n.checked_add(2))
            .ok_or(PhysicsError::InvalidInput("too many thermal stress loads"))?;
        self.view(study.into(), Output::ThermalStress(options))
    }

    pub(super) fn at_stress<R>(
        &self,
        cached: &Cached<'a>,
        options: ThermalStressOptions,
        cx: Option<&Cx<'_>>,
        action: impl FnOnce(
            &SimpThermalComplianceGoal<'a>,
            &SimpStressLinearization,
        ) -> Result<R, PhysicsError>,
    ) -> Result<R, PhysicsError> {
        let CachedGoal::Thermal(goal) = &cached.goal else {
            return Err(PhysicsError::InvalidInput("thermal stress requires a thermal goal"));
        };
        let mut stress = cached.stress.borrow_mut();
        if !stress.as_ref().is_some_and(|point| point.physical.options == options) {
            let prepared = self.controlled(self.loads.len() + 1, cx, |control| {
                goal.prepare_stress(options, control)
            })?;
            poll(cx)?;
            *stress = Some(prepared);
        }
        action(goal, stress.as_ref().expect("successful thermal stress preparation"))
    }
}
