//! Source-preserving fluid continuation for the existing acoustic heat owner.
//!
//! These sources describe prescribed isobaric backgrounds. Temperature changes
//! refresh constitutive response; they neither transport pore mass nor change
//! the authored thermal masses, geometry or fixed exterior acoustic reference.

use super::{IsobaricPorousGas, PorousHeatingError};
use fs_matdb::MaterialCard;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::GasState;
use fs_material::liquid::density::{LiquidDensityError, ResolvedLiquidDensityState};
use fs_material::state_point::MaterialStatePointError;
use std::sync::Arc;

/// One immutable fluid source evaluated at the accepted native temperature.
///
/// Implementations retain the original source, composition, background pressure
/// and domain policy. A returned state must describe exactly the requested
/// temperature and one coherent acoustic/caloric/transport response. Cloning a
/// source or cached state must preserve that physical binding, not reevaluate it.
/// The original gas source remains the default of all heating owners.
pub trait PorousTemperatureFluid: Clone {
    /// Original state, including its source evidence when the material has it.
    type State: ThermoviscousFluid + Clone;
    /// Whether a state evaluation also belongs to the legacy gas work counter.
    const IS_GAS: bool = false;
    /// Evaluate all fluid properties together; preserve original source errors.
    fn evaluate(&self, temperature_k: f64) -> Result<Self::State, PorousHeatingError>;
    /// Actual absolute temperature carried by the returned native state [K].
    fn temperature_k(&self, state: &Self::State) -> f64;
}

impl PorousTemperatureFluid for IsobaricPorousGas<'_> {
    type State = GasState;
    const IS_GAS: bool = true;
    fn evaluate(&self, temperature_k: f64) -> Result<GasState, PorousHeatingError> {
        IsobaricPorousGas::evaluate(*self, temperature_k)
    }
    fn temperature_k(&self, state: &GasState) -> f64 {
        state.temperature
    }
}

// Refuse a material adapter that silently substitutes a reference temperature.
// The original porous material owner admits the complete thermoviscous tuple
// when constructing its frequency response; no second acoustic law is added.
pub(super) fn evaluate<F: PorousTemperatureFluid>(
    fluid: &F,
    temperature_k: f64,
) -> Result<F::State, PorousHeatingError> {
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(PorousHeatingError::Invalid(
            "porous fluid requires a finite positive native temperature",
        ));
    }
    let state = fluid.evaluate(temperature_k)?;
    if fluid.temperature_k(&state).to_bits() != temperature_k.to_bits() {
        return Err(PorousHeatingError::Invalid(
            "porous fluid source returned a different temperature",
        ));
    }
    Ok(state)
}

/// Shared immutable evidence-bearing liquid evaluation. Spectral bins and
/// candidate thermal clones share its original property receipts, rather than
/// copying every source bundle or resampling a material for each frequency.
#[derive(Clone, Debug)]
pub struct PorousLiquidState {
    source: Arc<ResolvedLiquidDensityState>,
}
impl PorousLiquidState {
    /// Original density law, local thermodynamics and retained property receipts.
    pub fn source(&self) -> &ResolvedLiquidDensityState {
        &self.source
    }
}
impl FluidAcoustics for PorousLiquidState {
    fn acoustic_state(&self) -> Result<FluidAcousticState, fs_material::MaterialError> {
        self.source.acoustic_state()
    }
}
impl FluidTransport for PorousLiquidState {
    fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
        self.source.transport_state()
    }
}
impl ThermoviscousFluid for PorousLiquidState {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, fs_material::MaterialError> {
        self.source.thermoviscous_state()
    }
}

/// A liquid-saturated specimen at the original liquid's prescribed pressure.
///
/// The native density owner derives compressibility, sound speed and Cv from
/// the SAME rho(T,p) source and Cp while refreshing viscosity and conductivity.
/// Its original card, selected claims, typed coordinates and nonthermal context
/// are retained. There is no ideal-gas substitution or independently fitted
/// acoustic temperature factor. This is a locally stable single-phase response,
/// not a phase-equilibrium, fluid-motion or finite-pore-inventory solver.
#[derive(Clone, Copy, Debug)]
pub struct IsobaricPorousLiquid<'a> {
    card: &'a MaterialCard,
    reference: &'a ResolvedLiquidDensityState,
    temperature_support_k: [f64; 2],
}
impl<'a> IsobaricPorousLiquid<'a> {
    /// Admit explicit thermal support using the exact original liquid source at
    /// both endpoints. Every later temperature is also checked by that source;
    /// endpoint admission does not fill holes or extend derivative support.
    /// The reference's typed pressure, composition and other coordinates stay
    /// fixed. The native thermal law must justify this isobaric, thermally
    /// lumped fluid/skeleton storage; its mass is never replaced by rho(T)V.
    pub fn try_new(
        card: &'a MaterialCard,
        reference: &'a ResolvedLiquidDensityState,
        temperature_support_k: [f64; 2],
    ) -> Result<Self, PorousHeatingError> {
        if !temperature_support_k
            .into_iter()
            .all(|t| t.is_finite() && t > 0.0)
            || temperature_support_k[0] > temperature_support_k[1]
        {
            return Err(PorousHeatingError::Invalid(
                "isobaric pore liquid needs explicit ordered temperature support",
            ));
        }
        let source = Self {
            card,
            reference,
            temperature_support_k,
        };
        for temperature in temperature_support_k {
            source.evaluate(temperature)?;
        }
        Ok(source)
    }
    /// Original immutable material card; source substitution remains a refusal.
    pub const fn card(self) -> &'a MaterialCard {
        self.card
    }
    /// Original admitted reference and its prepared density-derivative source.
    pub const fn reference(self) -> &'a ResolvedLiquidDensityState {
        self.reference
    }
    /// Fixed background pressure in the original source convention [Pa].
    pub fn pressure_pa(self) -> f64 {
        self.reference
            .material()
            .query_point()
            .iter()
            .find(|(axis, _)| axis == self.reference.pressure_axis())
            .expect("native liquid retains its pressure coordinate")
            .1
    }
    /// Caller-declared inclusive thermal support, within original source support.
    pub const fn temperature_support_k(self) -> [f64; 2] {
        self.temperature_support_k
    }
    fn refusal(self, temperature_k: f64, source: LiquidDensityError) -> PorousHeatingError {
        PorousHeatingError::LiquidSource {
            card_identity: self.reference.material().card_identity(),
            reference_query: self.reference.material().query(),
            temperature_k,
            source,
        }
    }
}
impl PorousTemperatureFluid for IsobaricPorousLiquid<'_> {
    type State = PorousLiquidState;
    fn evaluate(&self, temperature_k: f64) -> Result<Self::State, PorousHeatingError> {
        if !temperature_k.is_finite()
            || !(self.temperature_support_k[0]..=self.temperature_support_k[1])
                .contains(&temperature_k)
        {
            return Err(self.refusal(
                temperature_k,
                MaterialStatePointError::InvalidDerived {
                    quantity: "thermal pore-liquid temperature left its declared support",
                }
                .into(),
            ));
        }
        let point = self
            .reference
            .material()
            .query()
            .with(self.reference.temperature_axis(), temperature_k)
            .map_err(|source| {
                self.refusal(
                    temperature_k,
                    MaterialStatePointError::Query {
                        property: self.reference.temperature_axis().to_owned(),
                        source,
                    }
                    .into(),
                )
            })?;
        let state = self
            .reference
            .at_query_point(self.card, &point)
            .map_err(|source| self.refusal(temperature_k, source))?;
        Ok(PorousLiquidState {
            source: Arc::new(state),
        })
    }
    fn temperature_k(&self, state: &PorousLiquidState) -> f64 {
        state.source.response().temperature_k()
    }
}
