//! Apply the original dissolution calorics to a material solution's live stocks.
//!
//! The thermal model must bind the exact ordered solubility laws and declared
//! solvent molar volume already held by the material owner. Its shared bath
//! accounts for sensible heat and the original dissolution enthalpy once, then
//! publishes the same finite inventory at its solved temperature. This adapter
//! supplies no additional heat law, transport integrator or phase kinetics.

use super::{FiniteSolubilityError, FiniteSolubilitySolution, MAX_OSMOTIC_SOLUTES};
use fs_qty::{Amount, Temperature, Volume};
use fs_thermochem::solution::solubility::caloric::{
    SolubilityCaloricBalanceV1, SolubilityCaloricErrorV1, SolubilityCaloricModelV1,
    SolubilityCaloricSolveV1,
};

/// Original thermal/material refusal or an incompatible source binding.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteSolubilityHeatError {
    /// The caloric model does not represent the original material pair laws.
    SourceMismatch,
    /// Original first-law, source-domain or cancellation refusal.
    Caloric(SolubilityCaloricErrorV1),
    /// Original dissolved-state or finite phase refusal.
    Material(FiniteSolubilityError),
}
impl core::fmt::Display for FiniteSolubilityHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite solubility heat step: {self:?}")
    }
}
impl std::error::Error for FiniteSolubilityHeatError {}
impl From<SolubilityCaloricErrorV1> for FiniteSolubilityHeatError {
    fn from(error: SolubilityCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}
impl From<FiniteSolubilityError> for FiniteSolubilityHeatError {
    fn from(error: FiniteSolubilityError) -> Self {
        Self::Material(error)
    }
}

/// A completed material endpoint and its original thermochemical energy balance.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteSolubilityHeatStep {
    solution: FiniteSolubilitySolution,
    balance: SolubilityCaloricBalanceV1,
}
impl FiniteSolubilityHeatStep {
    pub const fn solution(&self) -> &FiniteSolubilitySolution {
        &self.solution
    }
    pub const fn balance(&self) -> &SolubilityCaloricBalanceV1 {
        &self.balance
    }
    /// Continue material evolution without copying the accepted phase stocks.
    pub fn into_parts(self) -> (FiniteSolubilitySolution, SolubilityCaloricBalanceV1) {
        (self.solution, self.balance)
    }
}

impl FiniteSolubilitySolution {
    /// Close one isobaric heat input [J] against all current dissolved and solid
    /// stocks. Positive heat enters the bath; negative heat leaves it. Every
    /// original law and caloric curve must cover the requested solve bracket.
    ///
    /// Solvent volume and total component stocks stay fixed, including amounts
    /// previously changed by membrane transport. The input is borrowed and is
    /// unchanged on refusal or cancellation. The returned dissolved material
    /// view is rebuilt at the solved temperature with the same original laws.
    /// No simultaneous heat/mass transport or variable-pressure work is inferred.
    pub fn add_heat(
        &self,
        model: &SolubilityCaloricModelV1,
        heat_input_j: f64,
        controls: SolubilityCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteSolubilityHeatStep, FiniteSolubilityHeatError> {
        if cancelled() {
            return Err(SolubilityCaloricErrorV1::Cancelled.into());
        }
        if model.solvent_heat_capacity().species() != self.solvent()
            || model.solvent_molar_volume_m3_mol() != self.solution().solvent_molar_volume_m3_mol()
            || model.components().len() != self.components().len()
            || model
                .components()
                .iter()
                .zip(self.components())
                .any(|(thermal, material)| thermal.law() != material.law())
        {
            return Err(FiniteSolubilityHeatError::SourceMismatch);
        }
        let count = self.components().len();
        let mut totals = [Amount::new(0.0); MAX_OSMOTIC_SOLUTES];
        let mut dissolved = [Amount::new(0.0); MAX_OSMOTIC_SOLUTES];
        for (index, phase) in self.phase_partitions().iter().enumerate() {
            totals[index] = phase.total_solute;
            dissolved[index] = phase.dissolved;
        }
        let initial = model.prepare_at_solvent_volume(
            Temperature::new(self.solution().temperature_k()),
            Volume::new(self.solution().solvent_volume_m3()),
            &totals[..count],
            &dissolved[..count],
            &mut cancelled,
        )?;
        let balance = model.add_heat(&initial, heat_input_j, controls, &mut cancelled)?;
        let solution = self.with_temperature(balance.final_state().temperature())?;
        if cancelled() {
            return Err(SolubilityCaloricErrorV1::Cancelled.into());
        }
        Ok(FiniteSolubilityHeatStep { solution, balance })
    }
}

#[cfg(test)]
mod tests;
