//! The shared transactional heat-owner seam for the two admitted source laws.

use super::{FissionCellSource, FissionThermalError};
use crate::{EnergyChannels, FissionBudget, FissionError, FissionStep};
use fs_thermal::enthalpy::{CellId, MaterialId};

pub(super) mod sealed {
    pub trait Sealed {}
}

/// A cell-bound nuclear source using the shared thermal transaction.
///
/// This sealed interface keeps prescribed-event and neutron-irradiation laws
/// on one binding, scheduling, heat-credit and rollback implementation. Public
/// source types carry their own physical inventories and admitted budgets.
pub trait NuclearHeatSource: sealed::Sealed + Clone + core::fmt::Debug {
    /// Explicit work and numerical-accounting budget for this source law.
    type Budget: Copy + core::fmt::Debug;
    /// The source law's accepted interval report.
    type Step: Clone + core::fmt::Debug;

    /// Original thermal cell receiving only locally deposited energy.
    fn cell_id(&self) -> CellId;
    /// Original caloric material identity bound to the cell.
    fn material_id(&self) -> MaterialId;
    /// Advance-wide allowance for noncancelling nuclear and coupling defects.
    fn energy_tolerance_j(budget: Self::Budget) -> f64;
    /// Next authored discontinuity strictly after the accepted clock.
    fn next_transition(&self, time_s: f64) -> Option<f64>;
    /// Admit immutable source data, the initial environment and work limits.
    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError>;
    /// Refuse temperatures outside the declared nuclear-data domain.
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError>;
    /// Prepare an independent source candidate over the final thermal interval.
    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError>;
    /// Actual source release over an accepted interval.
    fn released(step: &Self::Step) -> EnergyChannels;
    /// Noncancelling numerical energy-account defects over that interval.
    fn energy_roundoff_j(step: &Self::Step) -> f64;
    /// Accepted fission events over that interval.
    fn fissions(step: &Self::Step) -> f64;
    /// Accepted competing captures over that interval.
    fn captures(step: &Self::Step) -> f64;
    /// The source's lifetime release, including all explicit initial reservoirs.
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError>;
}

impl sealed::Sealed for FissionCellSource {}
impl NuclearHeatSource for FissionCellSource {
    type Budget = FissionBudget;
    type Step = FissionStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.energy_tolerance_j
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.rates.next_after(time_s)
    }
    fn admit(&self, time_s: f64, _: f64, budget: Self::Budget) -> Result<(), FissionThermalError> {
        self.state
            .trial(0.0, self.rates.rate_at(time_s)?, budget, || Ok(()))?;
        Ok(())
    }
    fn validate_temperature(&self, _: f64) -> Result<(), FissionThermalError> {
        Ok(())
    }
    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        _: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        let (state, report) =
            self.state
                .trial(duration_s, self.rates.rate_at(time_s)?, budget, checkpoint)?;
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                state,
                rates: self.rates.clone(),
            },
            report,
        ))
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.energy_roundoff_j
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.fissions
    }
    fn captures(_: &Self::Step) -> f64 {
        0.0
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.state.cumulative_released())
    }
}
