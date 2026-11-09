//! Resolved liquid heat, vapor phase dynamics and causal compact sound.
//!
//! The native coupled owner advances the spherical liquid temperature field,
//! vapor phase inventory and radial motion on ONE adaptive clock. Conductive
//! heat comes from that field, including Stefan enthalpy transport and the
//! resolved share of the original liquid viscous dissipation. The bath is an
//! explicit outer boundary; no prescribed bubble heat-transfer coefficient is
//! added here. The source caloric approximation and spatial/temporal limits
//! remain those admitted by the original material and native thermal owner.
//!
//! Acoustic propagation uses the same original bulk-water density and sound
//! speed and the actual displaced-liquid source, 4*pi*R^2*U_liquid. The observer
//! adds no radiation damping to the incompressible dynamics. A failed observer
//! or final cancellation discards the complete radial AND thermal proposal.

use super::source_sample;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::TimedCompactSource;
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits, MonopoleSample};
#[cfg(feature = "caloric-cavitation")]
use fs_flux::cavitation::vapor::conducting::gas::GasVaporLimits;
use fs_flux::cavitation::vapor::conducting::{
    ConductingVaporBubble, ConductingVaporBubbleError, ConductingVaporBubbleNumerics,
    ConductingVaporBubbleObservation, ConductingVaporBubbleStep,
};
use fs_flux::cavitation::vapor::liquid::SphericalLiquidField;
use fs_flux::cavitation::vapor::{
    SaturatedVaporBubble, VaporBubbleError, VaporBubbleLimits, VaporBubbleState,
};
use fs_material::MaterialError;
use fs_material::cavitation::SaturatedWaterCavitation;
#[cfg(feature = "caloric-cavitation")]
use fs_material::cavitation::vapor::gas::ConservedBubbleGas;
use fs_material::fluid_acoustics::FluidAcoustics;

#[derive(Clone, Debug, PartialEq)]
pub enum MaterialConductingVaporBubbleError {
    Material(MaterialError),
    Dynamics(ConductingVaporBubbleError),
    Acoustic(AeroacError),
    ObserverMediumMismatch,
}
impl core::fmt::Display for MaterialConductingVaporBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material conducting vapor bubble: {self:?}")
    }
}
impl std::error::Error for MaterialConductingVaporBubbleError {}
impl From<MaterialError> for MaterialConductingVaporBubbleError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<ConductingVaporBubbleError> for MaterialConductingVaporBubbleError {
    fn from(error: ConductingVaporBubbleError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<VaporBubbleError> for MaterialConductingVaporBubbleError {
    fn from(error: VaporBubbleError) -> Self {
        Self::Dynamics(ConductingVaporBubbleError::Physics(error))
    }
}
impl From<AeroacError> for MaterialConductingVaporBubbleError {
    fn from(error: AeroacError) -> Self {
        Self::Acoustic(error)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObservedConductingVaporBubbleStep {
    pub bubble: ConductingVaporBubbleStep,
    pub acoustic: MonopoleSample,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactConductingVaporBubbleStep {
    pub bubble: ConductingVaporBubbleStep,
    pub source: TimedCompactSource,
}

/// One original water source and its complete radial/liquid-thermal history.
pub struct MaterialConductingVaporBubble<'a> {
    source: &'a SaturatedWaterCavitation,
    model: ConductingVaporBubble,
}

impl<'a> MaterialConductingVaporBubble<'a> {
    /// Bind a field compiled from the same source and its actual initial cell
    /// temperatures. The native owner checks source and geometry coherence.
    pub fn new(
        source: &'a SaturatedWaterCavitation,
        state: VaporBubbleState,
        limits: VaporBubbleLimits,
        liquid: SphericalLiquidField,
        temperatures_k: &[f64],
    ) -> Result<Self, MaterialConductingVaporBubbleError> {
        let bubble = SaturatedVaporBubble::new(*source, state, limits)?;
        let model = ConductingVaporBubble::new(bubble, liquid, temperatures_k)?;
        Ok(Self { source, model })
    }

    /// Initialize the same resolved-liquid trajectory with a conserved
    /// noncondensable inventory. Gas and water retain separate original
    /// caloric sources; their shared temperature, pressure work and latent
    /// exchange are advanced by the one native owner.
    ///
    /// The authored wet-mixture transport bounds must cover the actual
    /// temperature, total pressure and composition. They qualify the
    /// common-temperature and rapid-mixing approximation; the dry-gas
    /// transport source is not silently substituted for the wet mixture.
    #[cfg(feature = "caloric-cavitation")]
    pub fn with_noncondensable_gas(
        self,
        gas: ConservedBubbleGas,
        limits: GasVaporLimits,
    ) -> Result<Self, MaterialConductingVaporBubbleError> {
        Ok(Self {
            source: self.source,
            model: self.model.with_noncondensable_gas(gas, limits)?,
        })
    }

    /// Immutable original gas inventory, if this bubble contains one.
    #[cfg(feature = "caloric-cavitation")]
    pub fn gas(&self) -> Option<&ConservedBubbleGas> {
        self.model.gas()
    }

    pub const fn source(&self) -> &'a SaturatedWaterCavitation {
        self.source
    }
    pub const fn model(&self) -> &ConductingVaporBubble {
        &self.model
    }
    pub fn state(&self) -> VaporBubbleState {
        self.model.state()
    }
    pub fn liquid_temperatures_k(&self) -> &[f64] {
        self.model.liquid_temperatures_k()
    }
    pub fn observe(
        &self,
        absolute_pressure_pa: f64,
        bath_temperature_k: f64,
    ) -> Result<ConductingVaporBubbleObservation, MaterialConductingVaporBubbleError> {
        Ok(self
            .model
            .observe(absolute_pressure_pa, bath_temperature_k)?)
    }

    pub fn compact_source_sample(
        &self,
        absolute_pressure_pa: f64,
        bath_temperature_k: f64,
    ) -> Result<TimedCompactSource, MaterialConductingVaporBubbleError> {
        Ok(source_sample(
            &self
                .observe(absolute_pressure_pa, bath_temperature_k)?
                .bubble,
        )?)
    }

    pub fn monopole_observer(
        &self,
        distance_m: f64,
        maximum_frequency_hz: f64,
        limits: MonopoleLimits,
    ) -> Result<CompactMonopoleObserver, MaterialConductingVaporBubbleError> {
        Ok(CompactMonopoleObserver::new(
            self.source.acoustic_state()?,
            distance_m,
            maximum_frequency_hz,
            limits,
        )?)
    }

    /// Bath temperature and pressure are external boundary functions. They
    /// must be deterministic and resolved by the caller's interval controls.
    pub fn advance<Pressure, Bath, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &ConductingVaporBubbleNumerics,
        pressure: &Pressure,
        bath: &Bath,
        cancelled: &mut Cancel,
    ) -> Result<ConductingVaporBubbleStep, MaterialConductingVaporBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Bath: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        Ok(self
            .model
            .advance_controlled(end_time_s, numerics, pressure, bath, cancelled)?)
    }

    pub fn advance_compact_source<Pressure, Bath, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &ConductingVaporBubbleNumerics,
        pressure: &Pressure,
        bath: &Bath,
        cancelled: &mut Cancel,
    ) -> Result<CompactConductingVaporBubbleStep, MaterialConductingVaporBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Bath: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        let proposal = self
            .model
            .prepare_advance_controlled(end_time_s, numerics, pressure, bath, cancelled)?;
        let source = source_sample(&proposal.step().bubble.observation)?;
        if cancelled() {
            return Err(ConductingVaporBubbleError::Cancelled.into());
        }
        Ok(CompactConductingVaporBubbleStep {
            bubble: proposal.commit(),
            source,
        })
    }

    /// Admit the entire liquid field and phase trajectory with the observer.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_observed<Pressure, Bath, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &ConductingVaporBubbleNumerics,
        pressure: &Pressure,
        bath: &Bath,
        observer: &CompactMonopoleObserver,
        cancelled: &mut Cancel,
    ) -> Result<ObservedConductingVaporBubbleStep, MaterialConductingVaporBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Bath: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        if observer.medium() != self.source.acoustic_state()? {
            return Err(MaterialConductingVaporBubbleError::ObserverMediumMismatch);
        }
        let proposal = self
            .model
            .prepare_advance_controlled(end_time_s, numerics, pressure, bath, cancelled)?;
        let endpoint = &proposal.step().bubble.observation;
        let acoustic = observer.sample(
            endpoint.state.time_s(),
            endpoint.state.radius_m(),
            endpoint.liquid_volume_flow_rate_m3_s2,
        )?;
        if cancelled() {
            return Err(ConductingVaporBubbleError::Cancelled.into());
        }
        Ok(ObservedConductingVaporBubbleStep {
            bubble: proposal.commit(),
            acoustic,
        })
    }
}
