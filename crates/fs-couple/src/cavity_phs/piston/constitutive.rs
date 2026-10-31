//! Sealed source-material seams for the original gas/piston integrator.
//!
//! Only existing material owners provide a gas, its complete resolved state,
//! and the original conservative heat/volume candidate. These traits introduce
//! no equation of state, pressure approximation or additional mechanics.

use super::{ClosedGasParcel, GasEnergyError, GasEnergySolve, GasHeatVolumeCandidate};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;

mod sealed {
    pub trait Gas {}
    pub trait State {}
}

/// Complete source-resolved state consumed by native piston work and acoustics.
/// Implementations are sealed to the original ideal and nonideal gas owners.
/// Real-fluid values never come from an ideal-reference provenance view.
pub trait PistonGasState: sealed::State + Clone + core::fmt::Debug + PartialEq {
    /// Actual absolute mechanical pressure [Pa].
    fn pressure_pa(&self) -> f64;
    /// Actual mass density [kg/m3].
    fn density_kg_m3(&self) -> f64;
    /// Same-state adiabatic sound speed [m/s], including EOS departures.
    fn sound_speed_m_s(&self) -> f64;
    /// Actual absolute temperature [K] for the accepted heat-port account.
    fn temperature_k(&self) -> f64;
}

/// Original conserved gas storage used by [`super::GasPiston`].
///
/// The constitutive owner supplies the work-conjugate midpoint pressure and
/// its discrete volume tangent. Piston geometry, impulse, acceptance and port
/// accounting remain entirely within the existing mechanical integrator.
pub trait PistonGas: sealed::Gas + Clone + core::fmt::Debug + PartialEq {
    /// Complete resolved state; source identities remain on the material owner.
    type State: PistonGasState;
    /// Accepted gas volume [m3], also the sole piston-position coordinate.
    fn volume_m3(&self) -> f64;
    /// Accepted coherent thermodynamic, acoustic and transport state.
    fn state(&self) -> &Self::State;
    /// Immutable conservative heat/volume proposal through the original owner.
    fn midpoint_heat_volume(
        &self,
        volume_m3: f64,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasHeatVolumeCandidate<Self, Self::State>, GasEnergyError>;
}

impl sealed::State for ResolvedThermochemicalGasState {}
impl PistonGasState for ResolvedThermochemicalGasState {
    fn pressure_pa(&self) -> f64 {
        self.state().pressure
    }
    fn density_kg_m3(&self) -> f64 {
        self.state().density
    }
    fn sound_speed_m_s(&self) -> f64 {
        self.state().sound_speed
    }
    fn temperature_k(&self) -> f64 {
        self.state().temperature
    }
}

impl sealed::State for ResolvedThermochemicalRealGasState {}
impl PistonGasState for ResolvedThermochemicalRealGasState {
    fn pressure_pa(&self) -> f64 {
        Self::pressure_pa(self)
    }
    fn density_kg_m3(&self) -> f64 {
        Self::density_kg_m3(self)
    }
    fn sound_speed_m_s(&self) -> f64 {
        Self::sound_speed_m_s(self)
    }
    fn temperature_k(&self) -> f64 {
        Self::temperature_k(self)
    }
}

impl sealed::Gas for ClosedGasParcel {}
impl PistonGas for ClosedGasParcel {
    type State = ResolvedThermochemicalGasState;
    fn volume_m3(&self) -> f64 {
        Self::volume_m3(self)
    }
    fn state(&self) -> &Self::State {
        Self::state(self)
    }
    fn midpoint_heat_volume(
        &self,
        volume_m3: f64,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasHeatVolumeCandidate<Self, Self::State>, GasEnergyError> {
        Self::midpoint_heat_volume(self, volume_m3, heat_j, controls, cancelled)
    }
}

impl sealed::Gas for IsochoricRealGasParcel {}
impl PistonGas for IsochoricRealGasParcel {
    type State = ResolvedThermochemicalRealGasState;
    fn volume_m3(&self) -> f64 {
        Self::volume_m3(self)
    }
    fn state(&self) -> &Self::State {
        Self::state(self)
    }
    fn midpoint_heat_volume(
        &self,
        volume_m3: f64,
        heat_j: f64,
        controls: GasEnergySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasHeatVolumeCandidate<Self, Self::State>, GasEnergyError> {
        Self::midpoint_heat_volume(self, volume_m3, heat_j, controls, cancelled)
    }
}
