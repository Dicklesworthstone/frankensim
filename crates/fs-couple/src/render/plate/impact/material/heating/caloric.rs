//! Original finite-chart deposit shared by mechanical time owners.
use super::*;

/// A read-only trial; accepting mechanics is still required before publication.
#[derive(Clone, Copy)]
pub(super) struct CaloricDeposit {
    pub(super) state: EquilibriumPhaseState,
    pub(super) received_heat_j: f64,
    pub(super) received_correction_j: f64,
    pub(super) energy_change_j: f64,
}

impl CaloricReservoir<'_> {
    /// Preserve compensated cumulative work and the original adiabatic datum.
    /// The chart admits the actual candidate enthalpy before any state changes.
    pub(super) fn stage_deposit(
        &self,
        heat_j: f64,
        tolerance_j: f64,
    ) -> Result<CaloricDeposit, ImpactError> {
        if !heat_j.is_finite() || heat_j < 0.0 || !positive(tolerance_j) {
            return Err(invalid(
                "caloric deposit needs finite nonnegative heat and a positive allowance",
            ));
        }
        if heat_j == 0.0 {
            return Ok(CaloricDeposit {
                state: self.state,
                received_heat_j: self.received_heat_j,
                received_correction_j: self.received_correction_j,
                energy_change_j: 0.0,
            });
        }
        let increment = heat_j - self.received_correction_j;
        let total = self.received_heat_j + increment;
        let correction = (total - self.received_heat_j) - increment;
        let enthalpy = self.adiabatic_base.specific_enthalpy_j_kg()
            + (total - self.adiabatic_heat_j) / self.mass_kg;
        if !total.is_finite() || !enthalpy.is_finite() {
            return Err(invalid("contact thermal inventory overflow"));
        }
        let state = self
            .curve
            .state_at_specific_enthalpy(enthalpy)
            .map_err(|error| ImpactError::Owner(error.to_string()))?;
        let actual_h = state.specific_enthalpy_j_kg();
        let change = self.mass_kg * (actual_h - self.state.specific_enthalpy_j_kg());
        let inventory = self.mass_kg * (actual_h - self.initial.specific_enthalpy_j_kg());
        if !change.is_finite()
            || !inventory.is_finite()
            || (inventory - total - self.boundary_heat_j).abs() > tolerance_j
        {
            return Err(invalid(
                "contact thermal inventory exceeds its joule allowance",
            ));
        }
        Ok(CaloricDeposit {
            state,
            received_heat_j: total,
            received_correction_j: correction,
            energy_change_j: change,
        })
    }

    /// Infallible publication after all physical owners have accepted the step.
    pub(super) fn accept_deposit(&mut self, deposit: CaloricDeposit) {
        self.state = deposit.state;
        self.received_heat_j = deposit.received_heat_j;
        self.received_correction_j = deposit.received_correction_j;
    }
}

#[cfg(test)]
mod tests;
