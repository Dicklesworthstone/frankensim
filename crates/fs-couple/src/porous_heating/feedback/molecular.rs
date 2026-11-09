//! Original molecular reservoirs on the accepted isobaric thermal timeline.
//!
//! Each changed native cell temperature reevaluates the same equilibrium gas,
//! composition and molecular law once. Every frequency and coherent component
//! then consumes that immutable snapshot. Pore heat diffusion and molecular
//! compression share the existing material equation, without additive loss
//! fits. Background pressure/composition, specimen mass and exterior acoustic
//! normalization retain the parent owner's explicit quasisteady assumptions.

use super::{IsobaricPorousGas, PorousHeatingError, PorousTemperatureFluid};
use crate::cavity_modal::interior::heterogeneous::porous::AcousticVolumeRegionState;
use fs_duct::porous::PorousLayer;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::porous::{PorousFluidLaw, PorousFluidResponse};
use fs_material::relaxation::{RelaxationAbsorption, RelaxationState, RelaxingMedium};
use std::sync::Arc;

/// One evaluated original gas and molecular reservoir set. Cloned bins and
/// candidate caches share this immutable source; there are no self-references.
#[derive(Clone, Debug)]
pub struct PorousRelaxingState {
    source: Arc<RelaxationState>,
}
impl PorousRelaxingState {
    /// Complete original equilibrium gas, frozen heat capacity and reservoirs.
    pub fn source(&self) -> &RelaxationState {
        &self.source
    }
}
impl FluidAcoustics for PorousRelaxingState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, fs_material::MaterialError> {
        self.source.gas().acoustic_state()
    }
}
impl FluidTransport for PorousRelaxingState {
    fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
        self.source.gas().transport_state()
    }
}
impl ThermoviscousFluid for PorousRelaxingState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, fs_material::MaterialError> {
        self.source.gas().thermoviscous_state()
    }
}

/// Original gas specification, composition and relaxation law at fixed pressure.
/// The reference medium must exactly match the gas law at its own temperature.
/// Temperature support is checked by both original owners, including at each
/// later accepted state; endpoint admission does not extend either domain.
#[derive(Clone, Copy, Debug)]
pub struct IsobaricRelaxingPorousGas<'a> {
    gas: IsobaricPorousGas<'a>,
    law: &'a RelaxationAbsorption,
    reference: &'a RelaxingMedium,
}
impl<'a> IsobaricRelaxingPorousGas<'a> {
    pub fn try_new(
        gas: IsobaricPorousGas<'a>,
        law: &'a RelaxationAbsorption,
        reference: &'a RelaxingMedium,
    ) -> Result<Self, PorousHeatingError> {
        if gas.evaluate(reference.state().temperature)? != *reference.state() {
            return Err(PorousHeatingError::Invalid(
                "molecular pore reference differs from the original isobaric gas",
            ));
        }
        law.evaluate(reference)?;
        let source = Self {
            gas,
            law,
            reference,
        };
        for temperature in gas.temperature_support_k() {
            source.evaluate(temperature)?;
        }
        Ok(source)
    }
    /// Original pressure, equilibrium calorics, transport and thermal support.
    pub const fn gas(self) -> IsobaricPorousGas<'a> {
        self.gas
    }
    /// Original molecular rate and reservoir law, unchanged during feedback.
    pub const fn relaxation_law(self) -> &'a RelaxationAbsorption {
        self.law
    }
    /// Original reference gas and explicitly declared species/mole fractions.
    pub const fn reference_medium(self) -> &'a RelaxingMedium {
        self.reference
    }
}
impl PorousTemperatureFluid for IsobaricRelaxingPorousGas<'_> {
    type State = PorousRelaxingState;
    const IS_GAS: bool = true;

    fn evaluate(&self, temperature_k: f64) -> Result<Self::State, PorousHeatingError> {
        let gas = self.gas.evaluate(temperature_k)?;
        let medium = RelaxingMedium::new(gas, self.reference.components().to_vec())?;
        Ok(PorousRelaxingState {
            source: Arc::new(self.law.evaluate(&medium)?),
        })
    }
    fn temperature_k(&self, state: &Self::State) -> f64 {
        state.source.gas().temperature
    }
    fn volume_response(
        &self,
        state: &Self::State,
        material: Option<&PorousFluidLaw>,
        omega: f64,
    ) -> Result<AcousticVolumeRegionState, MaterialError> {
        match material {
            Some(law) => law
                .effective_relaxing(state.source(), omega)
                .map(AcousticVolumeRegionState::Porous),
            None => PorousFluidResponse::relaxing_fluid(state.source(), omega)
                .map(AcousticVolumeRegionState::MolecularGas),
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn prepare_layer<'a>(
        &self,
        material: &'a PorousFluidLaw,
        state: &Self::State,
        omega: f64,
        thickness_m: f64,
        area_m2: f64,
        reference_impedance: f64,
        tangential_wavenumber_rad_m: f64,
    ) -> Result<PorousLayer<'a>, PorousHeatingError> {
        Ok(PorousLayer::try_new_relaxing_oblique_shared(
            material,
            Arc::clone(&state.source),
            omega,
            thickness_m,
            area_m2,
            reference_impedance,
            tangential_wavenumber_rad_m,
        )?)
    }
}
