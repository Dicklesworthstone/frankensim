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
        SensibleHeatReservoirV1, SolubilityThermalInteractionV1,
        SolubilityThermalKineticControlsV1, SolubilityThermalKineticOutcomeV1,
        SolubilityThermalKineticStateV1, SolubilityThermalKineticsV1,
        SolubilityThermalReservoirOutcomeV1,
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

/// Accepted material and finite reservoir with their original closed heat account.
/// The reservoir lives only in the native outcome; it is never reconstructed
/// from its rounded temperature when continuing a later window.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticSolubilityReservoirStep {
    solution: KineticSolubilitySolution,
    balance: SolubilityThermalReservoirOutcomeV1,
}
impl KineticSolubilityReservoirStep {
    pub const fn solution(&self) -> &KineticSolubilitySolution {
        &self.solution
    }
    pub const fn balance(&self) -> &SolubilityThermalReservoirOutcomeV1 {
        &self.balance
    }
    /// Retain the native reservoir's original datum and compensated heat for
    /// another physical window, alongside the accepted material inventories.
    pub fn into_parts(
        self,
    ) -> (
        KineticSolubilitySolution,
        SolubilityThermalReservoirOutcomeV1,
    ) {
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
    /// Nonideal phases require their original explicitly qualified thermal
    /// family, including excess heat and the pure-solid reference correction.
    /// A fixed-T-only virial fit still refuses. The SAME family Arc is attached
    /// to the original core thermal owner; a different interaction or caloric
    /// source cannot silently replace it.
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<KineticSolubilityHeatStep, KineticSolubilityHeatError> {
        let (solution, balance) = self.advance_thermal(
            model,
            cancelled,
            |model, initial, cancelled| {
                model.advance(initial, duration, heat_input_j, controls, cancelled)
            },
            SolubilityThermalKineticOutcomeV1::final_state,
        )?;
        Ok(KineticSolubilityHeatStep { solution, balance })
    }

    /// Exchange heat with a reservoir held at its prescribed temperature.
    ///
    /// The actual boundary power is `G * (T_reservoir - T_solution)`, with
    /// `conductance_w_k = G` in W/K for the complete contact or jacket. It is
    /// an explicit apparatus input, not a molecular phase-transfer conductance
    /// or an intrinsic solubility property. The native thermal solver couples
    /// this temperature feedback to the actual finite-rate phase trajectories.
    /// Its accepted boundary heat is reported by `balance().heat_input_j()`;
    /// positive heat enters this solution.
    ///
    /// The reservoir may lie outside the material's temperature support, but
    /// every solution trial must stay inside the original supported bracket.
    /// Fixed volume, component conservation, source identity and nonideal
    /// excess calorics follow the same rules as [`Self::advance_with_heat`].
    /// The reservoir remains prescribed; it has no evolving finite heat store.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_heat_transfer(
        &self,
        model: &SolubilityThermalKineticsV1,
        duration: Time,
        reservoir_temperature: Temperature,
        conductance_w_k: f64,
        controls: SolubilityThermalKineticControlsV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<KineticSolubilityHeatStep, KineticSolubilityHeatError> {
        let (solution, balance) = self.advance_thermal(
            model,
            cancelled,
            |model, initial, cancelled| {
                model.advance_with_heat_transfer(
                    initial,
                    duration,
                    reservoir_temperature,
                    conductance_w_k,
                    controls,
                    cancelled,
                )
            },
            SolubilityThermalKineticOutcomeV1::final_state,
        )?;
        Ok(KineticSolubilityHeatStep { solution, balance })
    }

    /// Exchange heat with an isolated finite reservoir of constant heat capacity.
    ///
    /// Both actual temperatures drive the complete contact conductance [W/K].
    /// The SAME native thermal controller couples the reservoir's changing
    /// energy to sensible heat, dissolution/precipitation and any qualified
    /// nonideal excess calorics. No external heat is supplied. Reservoir heat
    /// is equal and opposite to solution heat, and the returned native outcome
    /// independently checks their total energy residual.
    ///
    /// Continue with `balance().reservoir()` from the accepted outcome so its
    /// original datum and compensated accumulated heat survive rounding of
    /// the displayed temperature. Neither input is changed on refusal or
    /// cancellation. Source, volume and component rules remain those of
    /// [`Self::advance_with_heat`]; reservoir phase changes, solvent/species
    /// transfer and changing contact conductance are outside this boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_heat_reservoir(
        &self,
        model: &SolubilityThermalKineticsV1,
        duration: Time,
        reservoir: &SensibleHeatReservoirV1,
        conductance_w_k: f64,
        controls: SolubilityThermalKineticControlsV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<KineticSolubilityReservoirStep, KineticSolubilityHeatError> {
        let (solution, balance) = self.advance_thermal(
            model,
            cancelled,
            |model, initial, cancelled| {
                model.advance_with_finite_heat_reservoir(
                    initial,
                    duration,
                    reservoir,
                    conductance_w_k,
                    controls,
                    cancelled,
                )
            },
            |outcome| outcome.balance().final_state(),
        )?;
        Ok(KineticSolubilityReservoirStep { solution, balance })
    }

    /// Bind and reconstruct exactly once around any native heat boundary.
    /// Every path retains the same cancellation points and material ownership.
    fn advance_thermal<C: FnMut() -> bool, O>(
        &self,
        model: &SolubilityThermalKineticsV1,
        mut cancelled: C,
        advance: impl FnOnce(
            &SolubilityThermalKineticsV1,
            &SolubilityThermalKineticStateV1,
            &mut C,
        ) -> Result<O, SolubilityCaloricErrorV1>,
        endpoint: impl FnOnce(&O) -> &SolubilityThermalKineticStateV1,
    ) -> Result<(Self, O), KineticSolubilityHeatError> {
        if cancelled() {
            return Err(SolubilityCaloricErrorV1::Cancelled.into());
        }
        if self.solution().nonideal_model().is_some() && self.thermal_nonideal_source().is_none() {
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
        let prepared;
        let model = if let Some(source) = self.thermal_nonideal_source() {
            if model.caloric() != source.caloric() {
                return Err(KineticSolubilityHeatError::SourceMismatch);
            }
            let interactions: Arc<dyn SolubilityThermalInteractionV1> = source.clone();
            if model
                .interaction_source()
                .is_some_and(|original| !Arc::ptr_eq(original, &interactions))
            {
                return Err(KineticSolubilityHeatError::SourceMismatch);
            }
            prepared = model.with_interactions(interactions)?;
            &prepared
        } else {
            if model.interaction_source().is_some() {
                return Err(KineticSolubilityHeatError::SourceMismatch);
            }
            model
        };
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
        let balance = advance(model, &initial, &mut cancelled)?;
        let endpoint = endpoint(&balance);
        let mut solution = Self::prepare(
            endpoint.temperature(),
            self.solution().solvent_molar_volume_m3_mol(),
            endpoint.solvent_volume(),
            Arc::clone(&self.models),
            endpoint.phases().iter().copied(),
        )?;
        if let Some(source) = self.thermal_nonideal_source() {
            solution = solution.with_thermal_nonideal_thermodynamics(Arc::clone(source))?;
        }
        if cancelled() {
            return Err(SolubilityCaloricErrorV1::Cancelled.into());
        }
        Ok((solution, balance))
    }
}

#[cfg(test)]
mod tests;
