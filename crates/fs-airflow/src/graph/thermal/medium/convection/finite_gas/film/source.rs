//! Original gas owners supplying source-pinned property-temperature evaluations.
use fs_couple::cavity_phs::piston::PistonGasState;
use fs_material::MaterialError;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonPhaseV1, ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel,
};
use fs_material::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};
use fs_qty::Temperature;

mod sealed {
    pub trait Source {}
}

/// Source bank used by the existing film convection evaluator.
/// Only original ideal and Peng--Robinson material owners provide this seam.
pub trait FilmGasSource: sealed::Source {
    type State: PistonGasState + ThermalFluidTransport;
    fn evaluate_property_temperature(
        &self,
        bulk: &Self::State,
        temperature: Temperature,
    ) -> Result<Self::State, MaterialError>;
}

impl sealed::Source for ThermochemicalGasModel {}
impl FilmGasSource for ThermochemicalGasModel {
    type State = ResolvedThermochemicalGasState;
    fn evaluate_property_temperature(
        &self,
        bulk: &Self::State,
        temperature: Temperature,
    ) -> Result<Self::State, MaterialError> {
        self.evaluate_temperature_from_snapshot(bulk, temperature)
    }
}

/// Explicit homogeneous branch for real-gas film-property evaluation.
/// The branch must reproduce the accepted bulk root. Stable-root selection is
/// refused by the material owner; it cannot silently switch film phases.
#[derive(Clone, Copy)]
pub struct PengRobinsonFilmSource<'a> {
    pub model: &'a ThermochemicalPengRobinsonModel,
    pub phase: PengRobinsonPhaseV1,
}
impl sealed::Source for PengRobinsonFilmSource<'_> {}
impl FilmGasSource for PengRobinsonFilmSource<'_> {
    type State = ResolvedThermochemicalRealGasState;
    fn evaluate_property_temperature(
        &self,
        bulk: &Self::State,
        temperature: Temperature,
    ) -> Result<Self::State, MaterialError> {
        self.model
            .evaluate_temperature_from_snapshot(bulk, temperature, self.phase)
    }
}
