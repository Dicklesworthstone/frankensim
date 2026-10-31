//! One consumer seam for the original ideal and real-gas storage owners.
//!
//! The parcel retains mass, volume, composition and its constitutive sources.
//! Coupled consumers exchange heat through that owner and receive its complete
//! material state; they do not reconstruct pressure or calorics themselves.

use super::{GasEnergyError, GasEnergySolve, GasHeatStep, IsochoricGasParcel};
use crate::MaterialError;
use crate::fluid_transport::ThermalFluidTransport;
use crate::thermochemical_gas::{
    ResolvedThermochemicalGasState, failure,
    nonideal::{ResolvedThermochemicalRealGasState, energy::parcel::IsochoricRealGasParcel},
};
use fs_qty::{Pressure, Temperature};

/// Closed, homogeneous, fixed-volume gas storage with frozen composition.
///
/// Implementations must conserve mass and volume during `add_heat`, update all
/// material views together and leave themselves unchanged on failure. The
/// source-law changes evaluate the original parcel's fixed-density path to the
/// requested temperature, canceling common energy/entropy references before
/// rounding. This is a storage interface, not a new equation of state or heat
/// integrator; source/phase admission remains with the implementing owner.
pub trait IsochoricGasStorage: Clone {
    /// The same source-bearing endpoint used for heat transfer and fluid flow.
    type State: ThermalFluidTransport;

    /// Complete accepted material endpoint.
    fn state(&self) -> &Self::State;
    /// Accepted absolute temperature [K].
    fn temperature_k(&self) -> f64;
    /// Conserved mass [kg], established from the original EOS density and volume.
    fn mass_kg(&self) -> f64;
    /// Add signed heat [J] through the original bounded energy closure.
    fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError>;
    /// Original-law internal-energy change at fixed density [J/kg].
    fn specific_internal_energy_change_to_temperature_j_kg(
        &self,
        temperature_k: f64,
    ) -> Result<f64, GasEnergyError>;
    /// Original-law entropy change at fixed density [J/(kg K)].
    fn specific_entropy_change_to_temperature_j_kg_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, MaterialError>;
}

impl IsochoricGasStorage for IsochoricGasParcel {
    type State = ResolvedThermochemicalGasState;

    fn state(&self) -> &Self::State {
        IsochoricGasParcel::state(self)
    }
    fn temperature_k(&self) -> f64 {
        self.state().state().temperature
    }
    fn mass_kg(&self) -> f64 {
        IsochoricGasParcel::mass_kg(self)
    }
    fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError> {
        IsochoricGasParcel::add_heat(self, heat_j, controls, cancelled)
    }
    fn specific_internal_energy_change_to_temperature_j_kg(
        &self,
        temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        self.model()
            .thermodynamics()
            .specific_internal_energy_change(
                Temperature::new(self.temperature_k()),
                Temperature::new(temperature_k),
            )
            .map(|value| value.value())
            .map_err(|error| failure("isochoric storage energy change", error).into())
    }
    fn specific_entropy_change_to_temperature_j_kg_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, MaterialError> {
        let initial = self.state().state();
        // For this ideal owner only, fixed density gives p1/p0 = T1/T0.
        // The original mixture law supplies caloric and pressure entropy;
        // neither a constant-Cv fit nor an absolute-entropy subtraction is used.
        let pressure_pa = initial.pressure * (temperature_k / initial.temperature);
        self.model()
            .thermodynamics()
            .specific_entropy_change(
                Temperature::new(initial.temperature),
                Pressure::new(initial.pressure),
                Temperature::new(temperature_k),
                Pressure::new(pressure_pa),
            )
            .map(|value| value.value())
            .map_err(|error| failure("isochoric storage entropy change", error))
    }
}

impl IsochoricGasStorage for IsochoricRealGasParcel {
    type State = ResolvedThermochemicalRealGasState;

    fn state(&self) -> &Self::State {
        IsochoricRealGasParcel::state(self)
    }
    fn temperature_k(&self) -> f64 {
        self.state().temperature_k()
    }
    fn mass_kg(&self) -> f64 {
        IsochoricRealGasParcel::mass_kg(self)
    }
    fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError> {
        IsochoricRealGasParcel::add_heat(self, heat_j, controls, cancelled)
    }
    fn specific_internal_energy_change_to_temperature_j_kg(
        &self,
        temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        self.model().specific_internal_energy_change_at_density(
            self.density_kg_m3(),
            self.temperature_k(),
            temperature_k,
        )
    }
    fn specific_entropy_change_to_temperature_j_kg_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, MaterialError> {
        self.model().specific_entropy_change_between_densities(
            self.density_kg_m3(),
            self.temperature_k(),
            self.density_kg_m3(),
            temperature_k,
        )
    }
}
