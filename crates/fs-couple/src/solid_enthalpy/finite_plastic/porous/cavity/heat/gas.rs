//! Sealed heat operations on the original ideal and real-gas parcel owners.
//!
//! This seam adds no equation of state, caloric approximation or thermal solve.
//! The existing adiabatic cavity interface still owns gas volume, direct
//! mechanical storage increments and source acoustic response.

use fs_material::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, GasHeatStep,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::nonideal::energy::parcel::ClosedRealGasParcel;
use fs_solid::hyper3::pressure::cavity::AdiabaticCavityGas;

mod sealed {
    pub trait Gas {}
}

/// Original physical gas storage for explicit cavity-wall heat exchange.
///
/// Implementations are sealed to the native ideal and Peng--Robinson parcels.
/// Heat is applied at the ACTUAL current volume and fixed original inventory.
/// The original caloric source supplies its temperature, error admission and
/// direct energy receipt; real-gas heating retains density-dependent departures
/// and its homogeneous-phase assumptions. No fixed-Cv or fixed-gamma substitute
/// enters either the heat stage or subsequent mechanical isentrope.
pub trait ThermalCavityGas: AdiabaticCavityGas + sealed::Gas {
    /// Same-state physical gas temperature [K].
    fn temperature_k(&self) -> f64;

    /// Original conserved closed inventory [kg].
    fn mass_kg(&self) -> f64;

    /// Original source-owned signed heat at fixed current physical volume.
    ///
    /// Native control validation, direct caloric increments, zero-heat behavior,
    /// cancellation and atomic refusal are retained. This operation alone is
    /// not a simultaneous heat/motion integrator or a gas-wall transfer law.
    fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError>;
}

impl sealed::Gas for ClosedGasParcel {}
impl ThermalCavityGas for ClosedGasParcel {
    fn temperature_k(&self) -> f64 {
        self.state().state().temperature
    }

    fn mass_kg(&self) -> f64 {
        ClosedGasParcel::mass_kg(self)
    }

    fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError> {
        ClosedGasParcel::add_heat(self, heat_j, controls, cancelled)
    }
}

impl sealed::Gas for ClosedRealGasParcel {}
impl ThermalCavityGas for ClosedRealGasParcel {
    fn temperature_k(&self) -> f64 {
        self.state().temperature_k()
    }

    fn mass_kg(&self) -> f64 {
        ClosedRealGasParcel::mass_kg(self)
    }

    fn add_heat(
        &mut self,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError> {
        ClosedRealGasParcel::add_heat(self, heat_j, controls, cancelled)
    }
}
