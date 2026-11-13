//! One accepted bubble history feeding the original compact acoustic observer.
//!
//! The material, bubble equation and sound propagation share rho and c. A
//! proposal is committed only after the complete physical interval and acoustic
//! observation are admitted. The selected equation owns its liquid response:
//! Rayleigh--Plesset is incompressible, while Keller--Miksis includes the
//! first-order compressibility correction. The passive observer removes no
//! additional energy from either trajectory. The polytropic adapters below
//! retain frozen vapor pressure and liquid temperature; no latent heat is
//! inferred from total bubble volume containing noncondensable gas. The
//! separate [`vapor`] adapter carries an explicit saturated-water phase model
//! and its distinct displaced-liquid acoustic source.

/// Phase-changing saturated water vapor and its displaced-liquid sound source.
pub mod vapor;

/// Original finite gas calorics on the same radial and acoustic trajectory.
#[cfg(feature = "caloric-cavitation")]
pub mod caloric;

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::{CompactSourceStrength, TimedCompactSource};
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits, MonopoleSample};
use fs_flux::cavitation::{
    BubbleGas, BubbleLimits, BubbleNumerics, BubbleObservation, BubblePressureSample, BubbleState,
    BubbleStep, CavitationError, KellerMiksisBubble, KellerMiksisLimits, RayleighPlessetBubble,
    SphericalBubbleModel,
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

/// Accepted dynamics and the same physical volume source, before propagation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompactBubbleStep {
    pub bubble: BubbleStep,
    pub source: TimedCompactSource,
}

/// Source-retaining spherical bubble with a single authoritative state clock.
///
/// The default equation is Rayleigh--Plesset. [`Self::new_compressible`] selects
/// Keller--Miksis and requires a pressure value AND its analytic time derivative
/// at every evaluation. Both use the same accepted-state and observation path.
pub struct MaterialBubble<
    'a,
    F: CavitationProperties + ?Sized,
    M: SphericalBubbleModel = RayleighPlessetBubble,
> {
    source: &'a F,
    material: CavitationState,
    model: M,
    state: BubbleState,
}

impl<'a, F: CavitationProperties + ?Sized> MaterialBubble<'a, F, RayleighPlessetBubble> {
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
}

impl<'a, F: CavitationProperties + ?Sized> MaterialBubble<'a, F, KellerMiksisBubble> {
    /// Admit a weakly compressible liquid trajectory with the original material.
    ///
    /// The equation itself includes first-order acoustic back-reaction. A
    /// subsequent monopole observation adds no second radiation damping term.
    /// The pressure derivative must come from the prescribed physical drive;
    /// it is never approximated by differencing an opaque pressure callback.
    pub fn new_compressible(
        source: &'a F,
        gas: BubbleGas,
        limits: KellerMiksisLimits,
        state: BubbleState,
        initial_pressure: BubblePressureSample,
    ) -> Result<Self, MaterialBubbleError> {
        let material = source.cavitation_state()?;
        let model = KellerMiksisBubble::try_new(material, gas, limits)?;
        model.observe(&state, initial_pressure)?;
        Ok(Self {
            source,
            material,
            model,
            state,
        })
    }
}

impl<'a, F: CavitationProperties + ?Sized, M: SphericalBubbleModel> MaterialBubble<'a, F, M> {
    /// Original material, including source/query provenance when present.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Shared frozen liquid/vapor coefficients for this physical trajectory.
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    /// The original equation owner; no duplicate bubble solver is introduced.
    pub const fn model(&self) -> &M {
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

    /// Read the accepted source with the pressure input required by its equation.
    pub fn observe(&self, pressure: M::Pressure) -> Result<BubbleObservation, MaterialBubbleError> {
        Ok(self.model.observe(&self.state, pressure)?)
    }

    /// Export the accepted volume source and its analytic time derivative.
    ///
    /// The native radial equation supplies `Q = 4 pi R² R'` and `Q' = V''`.
    /// No force source or numerical differentiation is inferred. This sample
    /// stays on the mechanics clock; propagation and receiver-clock sampling
    /// use the original liquid density and sound speed from [`Self::material`].
    /// Finite-source compactness, observation bandwidth, and any prescribed
    /// translation must be admitted separately. Translation does not change
    /// this radial model into a solver for bubble translational dynamics.
    pub fn compact_source_sample(
        &self,
        pressure: M::Pressure,
    ) -> Result<TimedCompactSource, MaterialBubbleError> {
        Ok(compact_source_from_observation(&self.observe(pressure)?)?)
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
        Pressure: Fn(f64) -> M::Pressure,
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

    /// Commit the complete physical interval and its source endpoint together.
    ///
    /// The endpoint retains both `Q` and `Q'` from the accepted native equation
    /// for coherent interpolation on a later receiver clock. A failed export
    /// or final cancellation leaves the original accepted state unchanged.
    /// The equation retains ownership of liquid response and losses; exporting
    /// a source adds no extra radiation damping or synthetic tail.
    pub fn advance_compact_source<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &BubbleNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<CompactBubbleStep, MaterialBubbleError>
    where
        Pressure: Fn(f64) -> M::Pressure,
        Cancel: FnMut() -> bool,
    {
        let mut candidate = self.state;
        let bubble = self.model.advance_controlled(
            &mut candidate,
            end_time_s,
            numerics,
            pressure,
            cancelled,
        )?;
        let source = compact_source_from_observation(&bubble.observation)?;
        if cancelled() {
            return Err(MaterialBubbleError::Dynamics(CavitationError::Cancelled));
        }
        self.state = candidate;
        Ok(CompactBubbleStep { bubble, source })
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
        Pressure: Fn(f64) -> M::Pressure,
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

fn compact_source_from_observation(
    observation: &BubbleObservation,
) -> Result<TimedCompactSource, AeroacError> {
    compact_source_from_flow(
        observation.state.time_s(),
        observation.volume_rate_m3_s,
        observation.volume_acceleration_m3_s2,
    )
}

fn compact_source_from_flow(
    source_time_s: f64,
    volume_flow_m3_s: f64,
    volume_flow_rate_m3_s2: f64,
) -> Result<TimedCompactSource, AeroacError> {
    if !source_time_s.is_finite()
        || !volume_flow_m3_s.is_finite()
        || !volume_flow_rate_m3_s2.is_finite()
    {
        return Err(AeroacError::NonFinite {
            what: "accepted bubble compact source",
        });
    }
    Ok(TimedCompactSource {
        source_time_s,
        strength: CompactSourceStrength {
            volume_flow_m3_s,
            volume_flow_rate_m3_s2,
            force_on_fluid_n: [0.0; 3],
            force_rate_on_fluid_n_s: [0.0; 3],
        },
    })
}
