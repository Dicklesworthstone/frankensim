//! Original industrial water state supplies every shared fluid response.

use super::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_thermochem::water::transport::If97TransportState;

impl FluidTransport for If97TransportState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        let original = self.thermodynamics();
        FluidTransportState::try_new(
            original.density().value(),
            original.cp().value(),
            self.thermal_conductivity().value(),
            self.dynamic_viscosity().value(),
        )
    }
}
impl ThermalFluidTransport for If97TransportState {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        let original = self.thermodynamics();
        ThermalFluidTransportState::try_new(
            self.transport_state()?,
            original.temperature().value(),
            Some(original.thermal_expansion().value()),
        )
    }
}
impl FluidAcoustics for If97TransportState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.thermodynamics().acoustic_state()
    }
}
impl ThermoviscousFluid for If97TransportState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        ThermoviscousFluidState::try_new(
            self.acoustic_state()?,
            self.transport_state()?,
            self.thermodynamics().cv().value(),
        )
    }
}
