//! Original material storage at the solid cavity's admitted volume.
//!
//! This sealed seam changes only which source owner supplies thermodynamics.
//! Geometry, mechanical pullback, acoustic factors and equilibrium stay shared.

use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::entropy::{GasEntropyError, GasEntropySolve};
use fs_material::thermochemical_gas::energy::parcel::mechanics::GasMechanicalLoad;
use fs_material::thermochemical_gas::energy::parcel::{ClosedGasParcel, GasVolumeStep};
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_material::thermochemical_gas::nonideal::energy::parcel::ClosedRealGasParcel;

mod sealed {
    pub trait Gas {}
}

/// Conserved, source-resolved gas storage for a reversible adiabatic cavity.
///
/// Implementations are sealed to the original ideal and Peng–Robinson parcel
/// owners. Each supplies its complete material state and the same shared
/// mechanical load. In particular, real-gas storage includes density-dependent
/// internal energy and its stiffness uses the actual `rho*c²`, never `gamma*p`.
/// The material owner admits the homogeneous phase branch and thermodynamic
/// domain; this interface adds neither a phase split nor a new equation of state.
pub trait AdiabaticCavityGas: sealed::Gas + Clone + core::fmt::Debug {
    /// Original coherent thermodynamic, acoustic and transport state.
    type State;

    /// Current physical volume [m³] of this unchanged mass inventory.
    fn volume_m3(&self) -> f64;

    /// Immutable endpoint and original-law work on the source isentrope.
    fn at_volume_isentropic(
        &self,
        volume_m3: f64,
        controls: GasEntropySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, GasVolumeStep), GasEntropyError>;

    /// Same-state pressure and acoustic compressibility on the admitted chart.
    fn mechanical_load_with_limit<'a>(
        &'a self,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
        maximum_coordinates: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMechanicalLoad<'a, Self::State>, GasEntropyError>;

    /// Direct original-law stored-energy change [J] between admitted endpoints.
    /// Both endpoints must retain the same inventory, model and auxiliary axes.
    fn stored_energy_change_to(
        &self,
        after: &Self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, GasEnergyError>;
}

impl sealed::Gas for ClosedGasParcel {}
impl AdiabaticCavityGas for ClosedGasParcel {
    type State = ResolvedThermochemicalGasState;

    fn volume_m3(&self) -> f64 {
        Self::volume_m3(self)
    }

    fn at_volume_isentropic(
        &self,
        volume_m3: f64,
        controls: GasEntropySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, GasVolumeStep), GasEntropyError> {
        let mut endpoint = self.clone();
        let step = endpoint.set_volume_isentropic(volume_m3, controls, cancelled)?;
        Ok((endpoint, step))
    }

    fn mechanical_load_with_limit<'a>(
        &'a self,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
        maximum_coordinates: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMechanicalLoad<'a, Self::State>, GasEntropyError> {
        Self::mechanical_load_with_limit(
            self,
            volume_m3,
            gradient,
            exterior_pressure_pa,
            maximum_coordinates,
            cancelled,
        )
    }

    fn stored_energy_change_to(
        &self,
        after: &Self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, GasEnergyError> {
        Self::stored_energy_change_to(self, after, cancelled)
    }
}

impl sealed::Gas for ClosedRealGasParcel {}
impl AdiabaticCavityGas for ClosedRealGasParcel {
    type State = ResolvedThermochemicalRealGasState;

    fn volume_m3(&self) -> f64 {
        Self::volume_m3(self)
    }

    fn at_volume_isentropic(
        &self,
        volume_m3: f64,
        controls: GasEntropySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, GasVolumeStep), GasEntropyError> {
        Self::at_volume_isentropic(self, volume_m3, controls, cancelled)
    }

    fn mechanical_load_with_limit<'a>(
        &'a self,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
        maximum_coordinates: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMechanicalLoad<'a, Self::State>, GasEntropyError> {
        Self::mechanical_load_with_limit(
            self,
            volume_m3,
            gradient,
            exterior_pressure_pa,
            maximum_coordinates,
            cancelled,
        )
    }

    fn stored_energy_change_to(
        &self,
        after: &Self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, GasEnergyError> {
        Self::stored_energy_change_to(self, after, cancelled)
    }
}
