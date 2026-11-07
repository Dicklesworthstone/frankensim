//! Local multicomponent equilibrium storage on the native heat owner.
//!
//! Phase amounts and compositions respond to heat at fixed total composition,
//! pressure and reference geometry. The original branch continuation and
//! unregularized caloric response retain redistribution energy and entropy.
//! No sampled enthalpy chart or independent apparent heat capacity is fitted.
//!
//! The complete temperature bracket must admit the same resolved phase count
//! and support. A disappearing phase or singular response refuses. General
//! finite discovery is deliberately not a continuous ThermalLaw: unresolved
//! first-order enthalpy gaps require a separate genuine coexistence relation.
//! This owner carries the original seeded local-branch qualifications, not a
//! global Gibbs-minimum or mathematical branch-identity certificate.
//!
//! Construction resolves the two endpoint material phase bundles. Subsequent
//! caloric queries retain the original pure Cp sources without refreshing
//! density, viscosity or phase conductivity. Resolve accepted phases through
//! the original material law before consuming those transport properties.

/// Original equilibrium calorics with explicitly declared bulk conductivity.
pub mod conduction;
/// Original species and enthalpy/entropy references on stationary mass flow.
pub mod stream;

use super::super::storage::minimum_heat_capacity_j_per_k;
use super::super::storage::stream::StreamBinding;
use super::{
    IsobaricHeatCapacityLaw, MulticomponentEquilibriumCaloricSolveV1,
    MulticomponentEquilibriumCaloricStateV1, MulticomponentLiquidCaloricError,
    MulticomponentLiquidCaloricLaw, ResolvedMulticomponentLiquidEquilibriumCaloricState, Result,
    SourceOwner, Temperature, poll,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use std::sync::Arc;

/// Immutable same-feed, same-phase-count equilibrium thermal storage.
/// Timeline cells may carry any positive mass of this fixed overall composition.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidEquilibriumThermalStorage {
    owner: Arc<SourceOwner>,
    endpoints: [MulticomponentEquilibriumCaloricStateV1; 2],
    mass_kg: f64,
    lower_specific_h_j_kg: f64,
    upper_specific_h_j_kg: f64,
    maximum_slope: f64,
    controls: MulticomponentEquilibriumCaloricSolveV1,
    stream: Option<StreamBinding>,
}

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Prepare the original resolved local branch for spatial heat transport.
    /// Both endpoints must admit the original material sources. Preparation
    /// uses two bounded original continuations: initial-to-cold, then cold-to-hot;
    /// each retains `controls.maximum_phase_solves` as its whole-call limit.
    /// The original heat inverse subsequently uses that same per-query limit.
    ///
    /// Closed storage is zero at the cold endpoint. The reusable dT/dh bound
    /// comes from the original pure-source Cp minima, rounded conservatively.
    /// On the admitted unshifted stable branch, equilibrium redistribution
    /// contributes nonnegative Cp, so this bound includes all phase response.
    pub fn equilibrium_thermal_storage(
        &self,
        initial: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidEquilibriumThermalStorage> {
        poll(&mut cancelled)?;
        self.admit_equilibrium_caloric(initial)?;
        let cold = self.model().equilibrium_at_temperature(
            initial.caloric_state(),
            controls.thermal.temperature_bracket[0],
            controls,
            &mut cancelled,
        )?;
        self.liquid
            .resolve_equilibrium(cold.candidate().clone(), &mut cancelled)?;
        let hot = self.model().equilibrium_at_temperature(
            &cold,
            controls.thermal.temperature_bracket[1],
            controls,
            &mut cancelled,
        )?;
        self.liquid
            .resolve_equilibrium(hot.candidate().clone(), &mut cancelled)?;
        let mass_kg = initial.equilibrium().total_mass_kg();
        let span = self.model().equilibrium_enthalpy_change_j(&cold, &hot)? / mass_kg;
        let minimum_capacity = minimum_heat_capacity_j_per_k(
            &self.owner,
            cold.component_moles(),
            controls.thermal.temperature_bracket,
            &mut cancelled,
        )?;
        let maximum_slope = (mass_kg / minimum_capacity).next_up();
        if !positive(span) || !positive(maximum_slope) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "unrepresentable equilibrium storage span or capacity bound",
            ));
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentLiquidEquilibriumThermalStorage {
            owner: Arc::clone(&self.owner),
            endpoints: [cold, hot],
            mass_kg,
            lower_specific_h_j_kg: 0.0,
            upper_specific_h_j_kg: span,
            maximum_slope,
            controls,
            stream: None,
        })
    }
}

impl MulticomponentLiquidEquilibriumThermalStorage {
    pub const fn temperature_domain(&self) -> [Temperature; 2] {
        self.controls.thermal.temperature_bracket
    }
    /// Inclusive specific enthalpy interval [J/kg]. Closed storage is zero at
    /// the cold endpoint; stream storage retains the common pure-source datum.
    pub const fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        [self.lower_specific_h_j_kg, self.upper_specific_h_j_kg]
    }
    /// Mass of the original normalization inventory, not a constraint on cell mass.
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub fn component_moles(&self) -> &[f64] {
        self.endpoints[0].component_moles()
    }
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw] {
        &self.owner.sources
    }

    /// Exact original local Cp, including changing phase amounts/compositions.
    pub fn specific_heat_cp_j_kg_k(
        &self,
        temperature: Temperature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        let state = self.caloric_state_at_temperature(temperature, cancelled)?;
        let cp = state.isobaric_heat_capacity_j_per_k() / self.mass_kg;
        if !positive(cp) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "unrepresentable equilibrium specific heat",
            ));
        }
        Ok(cp)
    }

    /// Direct sensible-plus-redistribution increment from the cold endpoint,
    /// plus the prepared datum. Original source constants cancel before the
    /// increment is evaluated; only stream storage restores the shared datum.
    pub fn specific_enthalpy_j_kg(
        &self,
        temperature: Temperature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        let state = self.caloric_state_at_temperature(temperature, cancelled)?;
        if temperature == self.temperature_domain()[0] {
            return Ok(self.lower_specific_h_j_kg);
        }
        if temperature == self.temperature_domain()[1] {
            return Ok(self.upper_specific_h_j_kg);
        }
        let increment = self
            .owner
            .model
            .equilibrium_enthalpy_change_j(&self.endpoints[0], &state)?
            / self.mass_kg;
        let h = self.lower_specific_h_j_kg + increment;
        if !positive(increment)
            || !h.is_finite()
            || h <= self.lower_specific_h_j_kg
            || h >= self.upper_specific_h_j_kg
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "interior equilibrium storage enthalpy is unrepresentable",
            ));
        }
        Ok(h)
    }

    /// Original bounded equilibrium heat inverse with its exact analytic Cp.
    /// No phase search or alternate branch is substituted after refusal.
    pub fn caloric_state_at_specific_enthalpy(
        &self,
        h: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricStateV1> {
        poll(&mut cancelled)?;
        self.admit_enthalpy(h)?;
        if h == self.lower_specific_h_j_kg {
            return Ok(self.endpoints[0].clone());
        }
        if h == self.upper_specific_h_j_kg {
            return Ok(self.endpoints[1].clone());
        }
        let heat = (h - self.lower_specific_h_j_kg) * self.mass_kg;
        if !positive(heat) {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "positive equilibrium storage heat is unrepresentable",
            ));
        }
        let balance = self.owner.model.add_heat_equilibrium(
            &self.endpoints[0],
            heat,
            self.controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(balance.final_state().clone())
    }

    fn caloric_state_at_temperature(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricStateV1> {
        poll(&mut cancelled)?;
        let [lo, hi] = self.temperature_domain();
        if !temperature.value().is_finite()
            || temperature.value() < lo.value()
            || temperature.value() > hi.value()
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "temperature outside equilibrium storage domain",
            ));
        }
        if temperature == lo {
            return Ok(self.endpoints[0].clone());
        }
        if temperature == hi {
            return Ok(self.endpoints[1].clone());
        }
        let state = self.owner.model.equilibrium_at_temperature(
            &self.endpoints[0],
            temperature,
            self.controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(state)
    }

    fn admit_enthalpy(&self, h: f64) -> Result<()> {
        if !h.is_finite() || h < self.lower_specific_h_j_kg || h > self.upper_specific_h_j_kg {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "enthalpy outside equilibrium storage domain",
            ));
        }
        Ok(())
    }
}

impl ThermalLaw for MulticomponentLiquidEquilibriumThermalStorage {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.stream.as_ref().map(StreamBinding::caloric_identity)
    }
    fn advected_composition(
        &self,
    ) -> Option<fs_thermal::enthalpy::advection::species::CaloricComposition<'_>> {
        self.stream.as_ref().map(StreamBinding::composition)
    }
    fn advected_entropy_anchor(
        &self,
    ) -> Option<fs_thermal::enthalpy::advection::species::CaloricEntropyAnchor<'_>> {
        self.stream.as_ref().map(StreamBinding::entropy_anchor)
    }
    fn temperature_k(&self, h: f64) -> core::result::Result<f64, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h, || false)
            .map(|state| state.temperature().value())
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> core::result::Result<Option<f64>, DepositionError> {
        self.admit_enthalpy(initial_h)
            .and_then(|()| self.admit_enthalpy(final_h))
            .map_err(|_| DepositionError::MaterialDomain)?;
        if initial_h == final_h {
            return Ok(Some(0.0));
        }
        let initial = self
            .caloric_state_at_specific_enthalpy(initial_h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let final_state = self
            .caloric_state_at_specific_enthalpy(final_h, || false)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let change = self
            .owner
            .model
            .equilibrium_entropy_change_j_per_k(&initial, &final_state)
            .map_err(|_| DepositionError::MaterialDomain)?
            / self.mass_kg;
        if !change.is_finite()
            || (final_h > initial_h && change <= 0.0)
            || (final_h < initial_h && change >= 0.0)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(change))
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
