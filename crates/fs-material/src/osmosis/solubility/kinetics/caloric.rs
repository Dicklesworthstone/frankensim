//! Finite-rate phase exchange and heat using the original material inventories.
//!
//! The thermochemical owner supplies all kinetic and caloric equations. This
//! adapter binds its exact ordered laws, conductances, surface geometry, growth
//! policies and solvent convention to the current material, then rebuilds its
//! aqueous view from actual accepted phases at the solved temperature. It introduces no additional
//! energy law or instantaneous equilibrium projection.

use super::{KineticSolubilityError, KineticSolubilitySolution, MAX_OSMOTIC_SOLUTES};
use fs_qty::{Temperature, Time, Volume};
use fs_thermochem::solution::solubility::caloric::{
    SolubilityCaloricErrorV1,
    kinetics::{
        SolubilityThermalKineticControlsV1, SolubilityThermalKineticOutcomeV1,
        SolubilityThermalKineticsV1,
    },
};
use std::sync::Arc;

/// Original thermal/material refusal or a different physical source binding.
#[derive(Clone, Debug, PartialEq)]
pub enum KineticSolubilityHeatError {
    SourceMismatch,
    Caloric(SolubilityCaloricErrorV1),
    Material(KineticSolubilityError),
}
impl core::fmt::Display for KineticSolubilityHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-rate soluble material heat step: {self:?}")
    }
}
impl std::error::Error for KineticSolubilityHeatError {}
impl From<SolubilityCaloricErrorV1> for KineticSolubilityHeatError {
    fn from(error: SolubilityCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}
impl From<KineticSolubilityError> for KineticSolubilityHeatError {
    fn from(error: KineticSolubilityError) -> Self {
        Self::Material(error)
    }
}

/// One accepted finite-rate material endpoint and its complete heat balance.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticSolubilityHeatStep {
    solution: KineticSolubilitySolution,
    balance: SolubilityThermalKineticOutcomeV1,
}
impl KineticSolubilityHeatStep {
    pub const fn solution(&self) -> &KineticSolubilitySolution {
        &self.solution
    }
    pub const fn balance(&self) -> &SolubilityThermalKineticOutcomeV1 {
        &self.balance
    }
    /// Continue material evolution without copying the accepted writable stocks.
    pub fn into_parts(self) -> (KineticSolubilitySolution, SolubilityThermalKineticOutcomeV1) {
        (self.solution, self.balance)
    }
}

impl KineticSolubilitySolution {
    /// Evolve finite-rate phase exchange with a uniformly delivered heat input.
    /// Positive heat enters the bath. Solvent volume and all component totals
    /// remain fixed; temperature and actual dissolved/solid amounts evolve
    /// together under the original thermochemical adaptive integrator.
    ///
    /// Every original law, conductance, surface geometry, growth policy and
    /// solvent convention must match. Refusal or cancellation leaves this
    /// material untouched.
    /// Explicitly bound virial phases refuse: the ideal caloric owner supplies
    /// neither their excess heat nor temperature-dependent interactions.
    /// The returned aqueous view uses the accepted dissolved amounts at the
    /// solved temperature, without equilibrium projection. Its zero-time phase
    /// receipts refresh target diagnostics; the window's heat and temporal
    /// diagnostics belong to `balance()`.
    pub fn advance_with_heat(
        &self,
        model: &SolubilityThermalKineticsV1,
        duration: Time,
        heat_input_j: f64,
        controls: SolubilityThermalKineticControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<KineticSolubilityHeatStep, KineticSolubilityHeatError> {
        if cancelled() {
            return Err(SolubilityCaloricErrorV1::Cancelled.into());
        }
        if self.solution().nonideal_model().is_some() {
            return Err(KineticSolubilityHeatError::Material(super::invalid(
                "nonideal heat evolution requires original excess calorics and temperature-dependent interactions",
            )));
        }
        if model.kinetics() != self.models()
            || model.caloric().solvent_molar_volume_m3_mol()
                != self.solution().solvent_molar_volume_m3_mol()
            || model.caloric().solvent_heat_capacity().species() != self.solvent()
        {
            return Err(KineticSolubilityHeatError::SourceMismatch);
        }
        let count = self.last_steps().len();
        let mut phases = [*self.last_steps()[0].state(); MAX_OSMOTIC_SOLUTES];
        for (phase, step) in phases.iter_mut().zip(self.last_steps()) {
            *phase = *step.state();
        }
        let initial = model.prepare_at_solvent_volume(
            Temperature::new(self.solution().temperature_k()),
            Volume::new(self.solution().solvent_volume_m3()),
            &phases[..count],
            &mut cancelled,
        )?;
        let balance = model.advance(&initial, duration, heat_input_j, controls, &mut cancelled)?;
        let endpoint = balance.final_state();
        let solution = Self::prepare(
            endpoint.temperature(),
            self.solution().solvent_molar_volume_m3_mol(),
            endpoint.solvent_volume(),
            Arc::clone(&self.models),
            endpoint.phases().iter().copied(),
        )?;
        if cancelled() {
            return Err(SolubilityCaloricErrorV1::Cancelled.into());
        }
        Ok(KineticSolubilityHeatStep { solution, balance })
    }
}

#[cfg(test)]
mod tests;
