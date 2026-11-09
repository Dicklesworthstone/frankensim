//! Nonlinear water phase exchange and sound from the displaced liquid.
//!
//! The native saturated-vapor owner advances radius, LIQUID interface velocity
//! and temperature with one signed heat port. Original IF97 phase mass,
//! calorics and saturation pressure drive the same trajectory. During phase
//! exchange the liquid velocity is not the geometric interface velocity:
//! the acoustic monopole is `Q_l = 4 pi R^2 U_l`, not `d(4 pi R^3/3)/dt`.
//! Its derivative is evaluated analytically by the native momentum owner.
//!
//! The original compact observer propagates that source through the SAME
//! liquid density and sound speed. It adds no radiation damping, heat or
//! independent phase evolution to the incompressible radial approximation.
//! The source and observer retain their own explicit low-Mach, compactness,
//! phase-equilibrium and thermal-model applicability conditions. Observer
//! refusal or cancellation discards the entire uncommitted physical window.

/// Resolved moving liquid heat storage on the same phase/acoustic trajectory.
pub mod conducting;

use super::compact_source_from_flow;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::TimedCompactSource;
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits, MonopoleSample};
use fs_flux::cavitation::vapor::{
    SaturatedVaporBubble, VaporBubbleError, VaporBubbleLimits, VaporBubbleNumerics,
    VaporBubbleObservation, VaporBubbleState, VaporBubbleStep, VaporBubbleThermalState,
};
use fs_material::MaterialError;
use fs_material::cavitation::{CavitationState, SaturatedWaterCavitation};
use fs_material::fluid_acoustics::FluidAcoustics;

#[derive(Clone, Debug, PartialEq)]
pub enum MaterialVaporBubbleError {
    Material(MaterialError),
    Dynamics(VaporBubbleError),
    Acoustic(AeroacError),
    ObserverMediumMismatch,
}
impl core::fmt::Display for MaterialVaporBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material vapor bubble: {self:?}")
    }
}
impl std::error::Error for MaterialVaporBubbleError {}
impl From<MaterialError> for MaterialVaporBubbleError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<VaporBubbleError> for MaterialVaporBubbleError {
    fn from(value: VaporBubbleError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<AeroacError> for MaterialVaporBubbleError {
    fn from(value: AeroacError) -> Self {
        Self::Acoustic(value)
    }
}

/// One accepted phase/radial window and its causal microphone endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObservedVaporBubbleStep {
    pub bubble: VaporBubbleStep,
    pub acoustic: MonopoleSample,
}

/// One accepted phase/radial window and its displaced-liquid source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactVaporBubbleStep {
    pub bubble: VaporBubbleStep,
    pub source: TimedCompactSource,
}

/// Original water source, native phase trajectory and passive sound binding.
/// There is one radial/thermal clock and no separately updated vapor parcel.
pub struct MaterialVaporBubble<'a> {
    source: &'a SaturatedWaterCavitation,
    model: SaturatedVaporBubble,
}

impl<'a> MaterialVaporBubble<'a> {
    pub fn new(
        source: &'a SaturatedWaterCavitation,
        state: VaporBubbleState,
        limits: VaporBubbleLimits,
    ) -> Result<Self, MaterialVaporBubbleError> {
        Ok(Self {
            source,
            model: SaturatedVaporBubble::new(*source, state, limits)?,
        })
    }

    pub const fn source(&self) -> &'a SaturatedWaterCavitation {
        self.source
    }
    pub fn material(&self) -> CavitationState {
        self.source.state()
    }
    pub const fn model(&self) -> &SaturatedVaporBubble {
        &self.model
    }
    pub fn state(&self) -> VaporBubbleState {
        self.model.state()
    }
    pub fn observe(
        &self,
        absolute_pressure_pa: f64,
        heat_power_w: f64,
    ) -> Result<VaporBubbleObservation, MaterialVaporBubbleError> {
        Ok(self.model.observe(absolute_pressure_pa, heat_power_w)?)
    }

    /// Exact native liquid volume flux and analytic derivative, on source time.
    /// Positive evaporation changes geometric bubble volume and displaced
    /// liquid volume differently; this adapter preserves that distinction.
    pub fn compact_source_sample(
        &self,
        absolute_pressure_pa: f64,
        heat_power_w: f64,
    ) -> Result<TimedCompactSource, MaterialVaporBubbleError> {
        source_sample(&self.observe(absolute_pressure_pa, heat_power_w)?).map_err(Into::into)
    }

    pub fn monopole_observer(
        &self,
        distance_m: f64,
        maximum_frequency_hz: f64,
        limits: MonopoleLimits,
    ) -> Result<CompactMonopoleObserver, MaterialVaporBubbleError> {
        Ok(CompactMonopoleObserver::new(
            self.source.acoustic_state()?,
            distance_m,
            maximum_frequency_hz,
            limits,
        )?)
    }

    /// Pressure and heat callbacks must be deterministic: rejected integration
    /// stages may revisit a time. Heat is signed power into the vapor system.
    pub fn advance<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &VaporBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<VaporBubbleStep, MaterialVaporBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Heat: Fn(f64, &VaporBubbleThermalState) -> Result<f64, VaporBubbleError>,
        Cancel: FnMut() -> bool,
    {
        Ok(self
            .model
            .advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?)
    }

    /// Admit the compact physical source before committing phase or motion.
    pub fn advance_compact_source<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &VaporBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<CompactVaporBubbleStep, MaterialVaporBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Heat: Fn(f64, &VaporBubbleThermalState) -> Result<f64, VaporBubbleError>,
        Cancel: FnMut() -> bool,
    {
        let proposal = self
            .model
            .prepare_advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?;
        let source = source_sample(&proposal.step().observation)?;
        if cancelled() {
            return Err(VaporBubbleError::Cancelled.into());
        }
        Ok(CompactVaporBubbleStep {
            bubble: proposal.commit(),
            source,
        })
    }

    /// Admit propagation through the same original medium before publication.
    /// The incompressible radial model receives no acoustic back-reaction.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_observed<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &VaporBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        observer: &CompactMonopoleObserver,
        cancelled: &mut Cancel,
    ) -> Result<ObservedVaporBubbleStep, MaterialVaporBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Heat: Fn(f64, &VaporBubbleThermalState) -> Result<f64, VaporBubbleError>,
        Cancel: FnMut() -> bool,
    {
        if observer.medium() != self.source.acoustic_state()? {
            return Err(MaterialVaporBubbleError::ObserverMediumMismatch);
        }
        let proposal = self
            .model
            .prepare_advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?;
        let endpoint = &proposal.step().observation;
        let acoustic = observer.sample(
            endpoint.state.time_s(),
            endpoint.state.radius_m(),
            endpoint.liquid_volume_flow_rate_m3_s2,
        )?;
        if cancelled() {
            return Err(VaporBubbleError::Cancelled.into());
        }
        Ok(ObservedVaporBubbleStep {
            bubble: proposal.commit(),
            acoustic,
        })
    }
}

fn source_sample(observation: &VaporBubbleObservation) -> Result<TimedCompactSource, AeroacError> {
    compact_source_from_flow(
        observation.state.time_s(),
        observation.liquid_volume_flow_m3_s,
        observation.liquid_volume_flow_rate_m3_s2,
    )
}
