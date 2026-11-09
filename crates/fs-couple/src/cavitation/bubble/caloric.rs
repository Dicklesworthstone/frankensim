//! Original finite gas calorics drive a physical bubble and its compact sound.
//!
//! The native radial owner advances radius, wall velocity and noncondensable
//! gas energy together. Its pressure and temperature come from the retained
//! gas material, and signed heat is an explicit external power port. The
//! original liquid coefficients, interface temperature and saturation pressure
//! stay frozen. Vapor phase kinetics and latent heat are not inferred from
//! total geometric volume or from heating the noncondensable gas.
//!
//! Acoustic endpoints use the native analytic volume acceleration and the
//! original liquid density and sound speed. Observation adds no force, heat,
//! radiation loss or second evolution of the gas. An unsuccessful observer or
//! cancellation discards the native proposal before any state is published.

use super::compact_source_from_observation;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::TimedCompactSource;
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits, MonopoleSample};
use fs_flux::cavitation::caloric::{
    CaloricBubble, CaloricBubbleEquation, CaloricBubbleError, CaloricBubbleLimits,
    CaloricBubbleNumerics, CaloricBubbleObservation, CaloricBubbleStep, CaloricGasObservation,
};
use fs_flux::cavitation::{BubblePressureSample, BubbleState, CavitationError};
use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// A refusal from one original physical owner; accepted history is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialCaloricBubbleError {
    Material(MaterialError),
    Dynamics(CaloricBubbleError),
    Acoustic(AeroacError),
    ObserverMediumMismatch,
}

impl core::fmt::Display for MaterialCaloricBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "caloric bubble liquid: {error}"),
            Self::Dynamics(error) => write!(f, "caloric bubble dynamics: {error}"),
            Self::Acoustic(error) => write!(f, "caloric bubble acoustics: {error}"),
            Self::ObserverMediumMismatch => {
                write!(f, "caloric bubble observer must use the original liquid medium")
            }
        }
    }
}
impl std::error::Error for MaterialCaloricBubbleError {}
impl From<MaterialError> for MaterialCaloricBubbleError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<CaloricBubbleError> for MaterialCaloricBubbleError {
    fn from(value: CaloricBubbleError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<AeroacError> for MaterialCaloricBubbleError {
    fn from(value: AeroacError) -> Self {
        Self::Acoustic(value)
    }
}

/// Complete gas/radial endpoint and its causally delayed microphone sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObservedCaloricBubbleStep {
    pub bubble: CaloricBubbleStep,
    pub acoustic: MonopoleSample,
}

/// Complete gas/radial endpoint and volume source on the physical source clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactCaloricBubbleStep {
    pub bubble: CaloricBubbleStep,
    pub source: TimedCompactSource,
}

/// One source-retaining liquid binding around the original finite-gas bubble.
///
/// This adapter changes no gas law and owns no second solver or caloric state.
/// A native borrowing proposal permits acoustic admission before publishing
/// radius, velocity, gas temperature/pressure, and their common accepted time.
pub struct MaterialCaloricBubble<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    model: CaloricBubble,
}

impl<'a, F: CavitationProperties + ?Sized> MaterialCaloricBubble<'a, F> {
    /// Bind the original liquid and a declared finite noncondensable gas.
    /// Gas volume must match the initial sphere; no inventory is resized.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        source: &'a F,
        gas: ClosedGasParcel,
        equation: CaloricBubbleEquation,
        state: BubbleState,
        initial_pressure: BubblePressureSample,
        limits: CaloricBubbleLimits,
        solve: GasEnergySolve,
    ) -> Result<Self, MaterialCaloricBubbleError> {
        let liquid = source.cavitation_state()?;
        let model = CaloricBubble::new(
            liquid,
            gas,
            equation,
            state,
            initial_pressure,
            limits,
            solve,
        )?;
        Ok(Self { source, model })
    }

    /// Original immutable liquid data and query identity.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Frozen liquid snapshot used in radial and acoustic response.
    pub fn material(&self) -> CavitationState {
        self.model.liquid()
    }
    /// Original accepted radial/caloric owner, available for read-only queries.
    pub const fn model(&self) -> &CaloricBubble {
        &self.model
    }
    /// Accepted radius, wall speed and physical time.
    pub fn state(&self) -> BubbleState {
        self.model.state()
    }
    /// Original finite gas material, mass and accepted thermodynamic state.
    pub fn gas(&self) -> &ClosedGasParcel {
        self.model.gas()
    }

    /// Read an endpoint using its actual signed heat power [W into gas].
    /// The heat rate contributes to the Keller--Miksis pressure derivative.
    pub fn observe(
        &self,
        pressure: BubblePressureSample,
        heat_power_w: f64,
    ) -> Result<CaloricBubbleObservation, MaterialCaloricBubbleError> {
        Ok(self.model.observe(pressure, heat_power_w)?)
    }

    /// Original analytic `Q` and `Q'` for coherent receiver-clock resampling.
    /// Source compactness, bandwidth and any imposed translation still need
    /// their own observer admission. No bubble translation is solved here.
    pub fn compact_source_sample(
        &self,
        pressure: BubblePressureSample,
        heat_power_w: f64,
    ) -> Result<TimedCompactSource, MaterialCaloricBubbleError> {
        Ok(compact_source_from_observation(
            &self.observe(pressure, heat_power_w)?.bubble,
        )?)
    }

    /// Prepare physical propagation through exactly the retained liquid.
    pub fn monopole_observer(
        &self,
        distance_m: f64,
        maximum_frequency_hz: f64,
        limits: MonopoleLimits,
    ) -> Result<CompactMonopoleObserver, MaterialCaloricBubbleError> {
        Ok(CompactMonopoleObserver::new(
            self.material().acoustic_state()?,
            distance_m,
            maximum_frequency_hz,
            limits,
        )?)
    }

    /// Advance the original gas/radial state atomically with explicit heat.
    /// Pressure and heat callbacks must be deterministic and free of side
    /// effects: rejected numerical stages can revisit a physical time.
    pub fn advance<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &CaloricBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<CaloricBubbleStep, MaterialCaloricBubbleError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, CaloricBubbleError>,
        Cancel: FnMut() -> bool,
    {
        Ok(self
            .model
            .advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?)
    }

    /// Admit one original volume source before committing the native proposal.
    pub fn advance_compact_source<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &CaloricBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<CompactCaloricBubbleStep, MaterialCaloricBubbleError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, CaloricBubbleError>,
        Cancel: FnMut() -> bool,
    {
        let proposal = self
            .model
            .prepare_advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?;
        let source = compact_source_from_observation(&proposal.step().observation.bubble)?;
        if cancelled() {
            return Err(CaloricBubbleError::from(CavitationError::Cancelled).into());
        }
        Ok(CompactCaloricBubbleStep {
            bubble: proposal.commit(),
            source,
        })
    }

    /// Commit gas, radius and the causal observer sample only after all pass.
    /// Observation never subtracts radiation work; the chosen radial equation
    /// already declares whether acoustic back-reaction is present.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_observed<Pressure, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &CaloricBubbleNumerics,
        pressure: &Pressure,
        heat: &Heat,
        observer: &CompactMonopoleObserver,
        cancelled: &mut Cancel,
    ) -> Result<ObservedCaloricBubbleStep, MaterialCaloricBubbleError>
    where
        Pressure: Fn(f64) -> BubblePressureSample,
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, CaloricBubbleError>,
        Cancel: FnMut() -> bool,
    {
        if observer.medium() != self.material().acoustic_state()? {
            return Err(MaterialCaloricBubbleError::ObserverMediumMismatch);
        }
        let proposal = self
            .model
            .prepare_advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?;
        let endpoint = &proposal.step().observation.bubble;
        let acoustic = observer.sample(
            endpoint.state.time_s(),
            endpoint.state.radius_m(),
            endpoint.volume_acceleration_m3_s2,
        )?;
        if cancelled() {
            return Err(CaloricBubbleError::from(CavitationError::Cancelled).into());
        }
        Ok(ObservedCaloricBubbleStep {
            bubble: proposal.commit(),
            acoustic,
        })
    }
}
