//! One accepted bubble history feeding the original compact acoustic observer.
//!
//! The material, bubble equation and sound propagation share rho and c. A
//! proposal is committed only after the complete physical interval and acoustic
//! observation are admitted. The observer is one-way: no extra radiation loss
//! is removed from an incompressible Rayleigh--Plesset trajectory. Vapor pressure
//! and liquid temperature stay frozen; no latent heat is manufactured from the
//! total volume of a bubble that also contains noncondensable gas.

use fs_aeroac::AeroacError;
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits, MonopoleSample};
use fs_flux::cavitation::{
    BubbleGas, BubbleLimits, BubbleNumerics, BubbleObservation, BubbleState, BubbleStep,
    CavitationError, RayleighPlessetBubble,
};
use fs_material::MaterialError;
use fs_material::cavitation::{CavitationProperties, CavitationState};
use fs_material::fluid_acoustics::FluidAcoustics;

/// Failure of an original constituent owner. No accepted state is changed.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialBubbleError {
    Material(MaterialError),
    Dynamics(CavitationError),
    Acoustic(AeroacError),
    /// An observer compiled from another medium cannot silently replace rho/c.
    ObserverMediumMismatch,
}
impl core::fmt::Display for MaterialBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "bubble material: {error}"),
            Self::Dynamics(error) => write!(f, "bubble dynamics: {error}"),
            Self::Acoustic(error) => write!(f, "bubble acoustics: {error}"),
            Self::ObserverMediumMismatch => write!(
                f,
                "bubble observer must use the same liquid density and sound speed"
            ),
        }
    }
}
impl std::error::Error for MaterialBubbleError {}
impl From<MaterialError> for MaterialBubbleError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<CavitationError> for MaterialBubbleError {
    fn from(error: CavitationError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<AeroacError> for MaterialBubbleError {
    fn from(error: AeroacError) -> Self {
        Self::Acoustic(error)
    }
}

/// Physical endpoint and the corresponding causal pressure sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObservedBubbleStep {
    pub bubble: BubbleStep,
    pub acoustic: MonopoleSample,
}

/// Source-retaining spherical bubble with a single authoritative state clock.
pub struct MaterialBubble<'a, F: CavitationProperties + ?Sized> {
    source: &'a F,
    material: CavitationState,
    model: RayleighPlessetBubble,
    state: BubbleState,
}

impl<'a, F: CavitationProperties + ?Sized> MaterialBubble<'a, F> {
    /// Admit the initial state against the actual initial absolute pressure.
    pub fn new(
        source: &'a F,
        gas: BubbleGas,
        limits: BubbleLimits,
        state: BubbleState,
        initial_absolute_pressure_pa: f64,
    ) -> Result<Self, MaterialBubbleError> {
        let material = source.cavitation_state()?;
        let model = RayleighPlessetBubble::try_new(material, gas, limits)?;
        model.observe(&state, initial_absolute_pressure_pa)?;
        Ok(Self {
            source,
            material,
            model,
            state,
        })
    }

    /// Original material, including source/query provenance when present.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Shared frozen liquid/vapor coefficients for this physical trajectory.
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    /// The original equation owner; no duplicate bubble solver is introduced.
    pub const fn model(&self) -> &RayleighPlessetBubble {
        &self.model
    }
    /// Accepted radius, wall velocity and physical time.
    pub const fn state(&self) -> BubbleState {
        self.state
    }

    /// Prepare an observer from exactly the liquid used in the bubble equation.
    /// The declared bandwidth must include the emitted signal's harmonics.
    pub fn monopole_observer(
        &self,
        distance_m: f64,
        maximum_frequency_hz: f64,
        limits: MonopoleLimits,
    ) -> Result<CompactMonopoleObserver, MaterialBubbleError> {
        Ok(CompactMonopoleObserver::new(
            self.material.acoustic_state()?,
            distance_m,
            maximum_frequency_hz,
            limits,
        )?)
    }

    /// Read the accepted physical source at the supplied absolute pressure.
    pub fn observe(
        &self,
        absolute_pressure_pa: f64,
    ) -> Result<BubbleObservation, MaterialBubbleError> {
        Ok(self.model.observe(&self.state, absolute_pressure_pa)?)
    }

    /// Advance the native physical history atomically, without an observer.
    pub fn advance<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<BubbleStep, MaterialBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        Ok(self.model.advance_controlled(
            &mut self.state,
            end_time_s,
            numerics,
            pressure,
            cancelled,
        )?)
    }

    /// Commit dynamics and observation together after both owners accept.
    ///
    /// Cancellation or observer-domain failure leaves radius, velocity and
    /// accepted time unchanged. A pressure callback must be deterministic and
    /// bounded; split calls at known drive discontinuities. A `Cx` cancellation
    /// check can be forwarded without introducing an alternate scheduler.
    pub fn advance_observed<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
        observer: &CompactMonopoleObserver,
        cancelled: &mut Cancel,
    ) -> Result<ObservedBubbleStep, MaterialBubbleError>
    where
        Pressure: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        if observer.medium() != self.material.acoustic_state()? {
            return Err(MaterialBubbleError::ObserverMediumMismatch);
        }
        let mut candidate = self.state;
        let bubble = self.model.advance_controlled(
            &mut candidate,
            end_time_s,
            numerics,
            pressure,
            cancelled,
        )?;
        let acoustic = observer.sample(
            candidate.time_s(),
            candidate.radius_m(),
            bubble.observation.volume_acceleration_m3_s2,
        )?;
        if cancelled() {
            return Err(MaterialBubbleError::Dynamics(CavitationError::Cancelled));
        }
        self.state = candidate;
        Ok(ObservedBubbleStep { bubble, acoustic })
    }
}
