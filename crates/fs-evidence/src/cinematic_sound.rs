//! Fail-closed configuration for Euler-disc soundtrack synthesis.
//!
//! This module freezes inputs and authority boundaries; it does not synthesize
//! samples. Physically informed output remains model-derived unless a separate
//! acoustic calibration receipt is supplied and admitted.

use core::fmt;
use fs_blake3::{ContentHash, hash_domain};

use crate::cinematic::{
    CinematicClock, CinematicClockDomain, DeclaredAcousticCalibrationReceipt, SoundAuthority,
};
use crate::cinematic_config::{CinematicComponentRef, CinematicComponentRole};

/// Exact schema version accepted by [`SoundSynthesisConfig`].
pub const SOUND_SYNTHESIS_SCHEMA_VERSION: u16 = 1;
/// Frozen reference-master audio sample rate.
pub const SOUND_MASTER_SAMPLE_RATE_HZ: u32 = 48_000;
/// Frozen reference-master video rate numerator (denominator is one).
pub const SOUND_MASTER_VIDEO_RATE_HZ: u32 = 24;
/// Resource ceiling for admitted structural modes.
pub const MAX_SOUND_MODES: usize = 256;

const IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.sound-config.v1";

/// Listener coordinates. The reference composition requires camera-relative geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ListenerFrame {
    /// Position, forward, and up are expressed in the animated camera frame.
    AnimatedCamera = 1,
    /// World-space coordinates; represented so accidental use is refused explicitly.
    World = 2,
}

impl ListenerFrame {
    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::AnimatedCamera => "animated-camera",
            Self::World => "world",
        }
    }
}

/// Listener/microphone pose in metres with a right-handed orthonormal view basis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ListenerPose {
    /// Coordinate frame for every vector below.
    pub frame: ListenerFrame,
    /// Listener position in metres.
    pub position_m: [f64; 3],
    /// Unit vector toward the listener's forward direction.
    pub forward: [f64; 3],
    /// Unit vector toward the listener's up direction.
    pub up: [f64; 3],
}

/// Audio channel layout for the master.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum SoundChannelLayout {
    /// One full-band channel; not accepted by the v1 reference master.
    Mono = 1,
    /// Left/right stereo master.
    Stereo = 2,
}

impl SoundChannelLayout {
    /// Number of interleaved channels.
    #[must_use]
    pub const fn channels(self) -> u8 {
        self as u8
    }

    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Mono => "mono",
            Self::Stereo => "stereo",
        }
    }
}

/// Simulation channels allowed to drive the modal source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum SoundExcitationChannel {
    /// Normal contact force in newtons.
    ContactNormalForce = 1,
    /// Tangential contact force in newtons.
    ContactTangentialForce = 2,
    /// Base reaction force in newtons.
    BaseReactionForce = 3,
    /// Disc angular speed in radians per second.
    DiscAngularSpeed = 4,
    /// Disc precession rate in radians per second.
    PrecessionRate = 5,
}

/// Dimensionally explicit mapping from one simulation channel into the
/// normalized scalar source coordinate that drives every admitted mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SoundExcitationControl {
    /// Typed input channel, including its documented SI unit.
    pub channel: SoundExcitationChannel,
    /// Normalized source-coordinate units per input-channel unit.
    pub source_scale: f64,
}

/// One deterministic, damped structural mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SoundMode {
    /// Stable caller-assigned mode identifier. Zero is invalid.
    pub mode_id: u32,
    /// Undamped natural frequency in hertz.
    pub frequency_hz: f64,
    /// Dimensionless damping ratio, in `(0, 1]`.
    pub damping_ratio: f64,
    /// Signed digital-full-scale gain per normalized source-coordinate unit.
    pub gain: f64,
    /// Identity of material parameters used for this mode.
    pub material_identity: ContentHash,
    /// Identity of base/support parameters used for this mode.
    pub base_identity: ContentHash,
}

/// Declared room/spatialization treatment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SoundRoomResponse {
    /// No reflected-room contribution.
    Dry,
    /// Apply the room component's admitted SI impulse response.
    DeclaredImpulseResponse {
        /// Wet contribution in `[0, 1]`.
        wet_mix: f64,
    },
}

/// Amplitude reference carried by the waveform and derived level metrics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SoundAmplitudeReference {
    /// Non-calibrated floating-point digital full scale.
    DigitalFullScale {
        /// Reserved peak headroom below 0 dBFS.
        headroom_db: f64,
    },
    /// Physical pressure scale backed by the supplied calibration receipt.
    CalibratedPressure {
        /// Pascals RMS represented by 0 dBFS before headroom is applied.
        pascal_rms_at_full_scale: f64,
        /// Reserved peak headroom below 0 dBFS.
        headroom_db: f64,
    },
}

/// Source trajectory outcome presented to the synthesizer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum SoundTrajectoryDisposition {
    /// A localized physical terminal event ends accepted state.
    PhysicalTerminal = 1,
    /// The trajectory ends at its declared integration horizon.
    HorizonCensored = 2,
    /// Numerical refusal means no physical waveform may be emitted.
    NumericalRefusal = 3,
}

impl SoundTrajectoryDisposition {
    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::PhysicalTerminal => "physical-terminal",
            Self::HorizonCensored => "horizon-censored",
            Self::NumericalRefusal => "numerical-refusal",
        }
    }
}

/// End-of-trajectory audio policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SoundTerminalPolicy {
    /// End at the last accepted sample with a bounded deterministic fade.
    FadeAtLastAccepted {
        /// Fade length in sample frames.
        fade_sample_frames: u32,
    },
    /// Emit no waveform. Required for numerical refusal.
    Silence,
}

/// Explicit model assumptions disclosed with every informed soundtrack.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum SoundModelAssumption {
    /// Modal responses combine linearly.
    LinearModalSuperposition = 1,
    /// Damping is time invariant over the rendered interval.
    TimeInvariantDamping = 2,
    /// The declared excitation channels adequately represent source coupling.
    DeclaredExcitationCompleteness = 3,
    /// Room response is dry or exactly the declared response artifact.
    DeclaredRoomResponse = 4,
}

/// Complete caller input. No field has an implicit default.
#[derive(Debug, Clone, PartialEq)]
pub struct SoundSynthesisInput {
    /// Exact supported schema version.
    pub schema_version: u16,
    /// Separate sound authority tier.
    pub authority: SoundAuthority,
    /// Animation-grade trajectory reference.
    pub trajectory: CinematicComponentRef,
    /// Simulation-to-excitation control artifact.
    pub excitation: CinematicComponentRef,
    /// Modal/synthesis model artifact.
    pub sound_model: CinematicComponentRef,
    /// Listener/microphone artifact.
    pub microphone: CinematicComponentRef,
    /// Room/spatialization artifact.
    pub room: CinematicComponentRef,
    /// Master timeline artifact.
    pub timeline: CinematicComponentRef,
    /// Exact video-frame clock.
    pub video_clock: CinematicClock,
    /// Exact audio sample-frame clock.
    pub audio_clock: CinematicClock,
    /// Frozen master layout.
    pub channel_layout: SoundChannelLayout,
    /// Listener relative to the animated camera.
    pub listener: ListenerPose,
    /// Ordered, dimensionally explicit excitation controls.
    pub excitation_controls: Vec<SoundExcitationControl>,
    /// Ordered modes by strictly increasing `mode_id`.
    pub modes: Vec<SoundMode>,
    /// Room treatment.
    pub room_response: SoundRoomResponse,
    /// Amplitude/reference convention.
    pub amplitude_reference: SoundAmplitudeReference,
    /// Source trajectory outcome.
    pub trajectory_disposition: SoundTrajectoryDisposition,
    /// End-of-trajectory policy.
    pub terminal_policy: SoundTerminalPolicy,
    /// Resampler implementation/configuration identity.
    pub resampler_identity: ContentHash,
    /// Nonzero resampler version.
    pub resampler_version: u32,
    /// Reconstruction/anti-alias filter identity.
    pub filter_identity: ContentHash,
    /// Nonzero filter version.
    pub filter_version: u32,
    /// Explicit ordered assumption set.
    pub assumptions: Vec<SoundModelAssumption>,
    /// Required only for calibrated authority.
    pub calibration: Option<DeclaredAcousticCalibrationReceipt>,
}

/// Admitted immutable sound configuration and deterministic receipt identity.
#[derive(Debug, Clone, PartialEq)]
pub struct SoundSynthesisConfig {
    input: SoundSynthesisInput,
    identity: ContentHash,
}

/// Compact receipt whose configuration identity transitively binds every
/// admitted synthesis field while surfacing the principal upstream identities.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SoundSynthesisReceipt {
    /// Exact sound schema version.
    pub schema_version: u16,
    /// Identity of the complete [`SoundSynthesisInput`].
    pub configuration_identity: ContentHash,
    /// Admitted authority applying to waveform and derived metrics.
    pub authority: SoundAuthority,
    /// Accepted trajectory identity.
    pub trajectory_identity: ContentHash,
    /// Excitation-control artifact identity.
    pub excitation_identity: ContentHash,
    /// Synthesis-model artifact identity.
    pub sound_model_identity: ContentHash,
    /// Shared timeline artifact identity.
    pub timeline_identity: ContentHash,
}

/// Outputs that must carry the soundtrack's authority and assumptions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SoundOutputKind {
    /// Interleaved audio waveform.
    Waveform,
    /// Frequency, bandwidth, or spectral-centroid metric.
    SpectralMetric,
    /// Digital or calibrated acoustic level metric.
    LevelMetric,
}

/// Non-owning authority declaration attached to a waveform or derived metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SoundOutputDeclaration<'a> {
    /// Kind of output being declared.
    pub kind: SoundOutputKind,
    /// Authority inherited from the admitted synthesis configuration.
    pub authority: SoundAuthority,
    /// Complete admitted assumption set.
    pub assumptions: &'a [SoundModelAssumption],
    /// True only for calibrated structural acoustics.
    pub calibrated_acoustic_prediction: bool,
}

impl SoundSynthesisConfig {
    /// Validate and bind every synthesis parameter into one content identity.
    pub fn try_admit(input: SoundSynthesisInput) -> Result<Self, SoundSynthesisError> {
        validate(&input)?;
        let identity = sound_identity(&input);
        Ok(Self { input, identity })
    }

    /// Identity of the complete admitted synthesis configuration.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Admitted source input.
    #[must_use]
    pub fn input(&self) -> &SoundSynthesisInput {
        &self.input
    }

    /// Exact synthesis authority; it applies equally to waveform and metrics.
    #[must_use]
    pub const fn authority(&self) -> SoundAuthority {
        self.input.authority
    }

    /// Produce the stable synthesis receipt. Every omitted detail is committed
    /// transitively by `configuration_identity`, not treated as an implicit default.
    #[must_use]
    pub const fn receipt(&self) -> SoundSynthesisReceipt {
        SoundSynthesisReceipt {
            schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
            configuration_identity: self.identity,
            authority: self.input.authority,
            trajectory_identity: self.input.trajectory.identity(),
            excitation_identity: self.input.excitation.identity(),
            sound_model_identity: self.input.sound_model.identity(),
            timeline_identity: self.input.timeline.identity(),
        }
    }

    /// Declare the authority and assumptions for one waveform or derived metric.
    #[must_use]
    pub fn declare_output(&self, kind: SoundOutputKind) -> SoundOutputDeclaration<'_> {
        SoundOutputDeclaration {
            kind,
            authority: self.input.authority,
            assumptions: &self.input.assumptions,
            calibrated_acoustic_prediction: self.input.authority == SoundAuthority::Calibrated,
        }
    }

    /// Human/machine metadata generated from the admitted state.
    #[must_use]
    pub fn to_manifest_json(&self) -> String {
        format!(
            "{{\"schema_version\":{},\"identity\":\"{}\",\"authority\":\"{}\",\"synthesis_class\":\"{}\",\"sample_rate_hz\":{},\"channels\":{},\"channel_layout\":\"{}\",\"listener_frame\":\"{}\",\"trajectory_disposition\":\"{}\",\"mode_count\":{},\"calibrated_acoustic_prediction\":{}}}",
            SOUND_SYNTHESIS_SCHEMA_VERSION,
            self.identity.to_hex(),
            self.input.authority.code(),
            synthesis_class_code(self.input.authority),
            SOUND_MASTER_SAMPLE_RATE_HZ,
            self.input.channel_layout.channels(),
            self.input.channel_layout.code(),
            self.input.listener.frame.code(),
            self.input.trajectory_disposition.code(),
            self.input.modes.len(),
            self.input.authority == SoundAuthority::Calibrated,
        )
    }
}

/// Precise fail-closed admission errors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SoundSynthesisError {
    /// Unknown schema version.
    UnsupportedSchemaVersion(u16),
    /// A component had the wrong semantic role.
    WrongComponentRole {
        /// Required role.
        expected: CinematicComponentRole,
        /// Supplied role.
        got: CinematicComponentRole,
    },
    /// Audio clock is not exact 48 kHz sample time.
    InvalidAudioClock,
    /// Video clock is not exact 24/1 frame time.
    InvalidVideoClock,
    /// Video and audio timeline origins or endpoints differ.
    AudioVideoTimelineMismatch,
    /// Reference master is not stereo.
    InvalidChannelLayout,
    /// Listener must be relative to the animated camera.
    InvalidListenerFrame,
    /// Listener pose contains non-finite or non-orthonormal geometry.
    InvalidListenerPose,
    /// Physically informed/calibrated synthesis has no excitation channels.
    MissingExcitationChannels,
    /// Excitation channels are duplicated or out of canonical order.
    NonCanonicalExcitationChannels,
    /// Physically informed/calibrated synthesis has no modes.
    MissingModes,
    /// Mode resource ceiling exceeded.
    TooManyModes,
    /// Mode values or order are invalid.
    InvalidMode,
    /// Room treatment is invalid.
    InvalidRoomResponse,
    /// Amplitude reference/headroom is invalid for the authority tier.
    InvalidAmplitudeReference,
    /// Numerical refusal or accepted-state termination has an incompatible policy.
    InvalidTerminalPolicy,
    /// Algorithm identity/version is absent.
    InvalidAlgorithmReference,
    /// Assumptions are duplicated, unordered, or incomplete.
    InvalidAssumptions,
    /// Calibrated authority lacks calibration evidence.
    MissingCalibration,
    /// A lower authority tier attempted to carry calibration as implicit promotion.
    UnexpectedCalibration,
}

impl SoundSynthesisError {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::UnsupportedSchemaVersion(_) => "sound-unsupported-schema-version",
            Self::WrongComponentRole { .. } => "sound-wrong-component-role",
            Self::InvalidAudioClock => "sound-invalid-audio-clock",
            Self::InvalidVideoClock => "sound-invalid-video-clock",
            Self::AudioVideoTimelineMismatch => "sound-av-timeline-mismatch",
            Self::InvalidChannelLayout => "sound-invalid-channel-layout",
            Self::InvalidListenerFrame => "sound-invalid-listener-frame",
            Self::InvalidListenerPose => "sound-invalid-listener-pose",
            Self::MissingExcitationChannels => "sound-missing-excitation-channels",
            Self::NonCanonicalExcitationChannels => "sound-noncanonical-excitation-channels",
            Self::MissingModes => "sound-missing-modes",
            Self::TooManyModes => "sound-too-many-modes",
            Self::InvalidMode => "sound-invalid-mode",
            Self::InvalidRoomResponse => "sound-invalid-room-response",
            Self::InvalidAmplitudeReference => "sound-invalid-amplitude-reference",
            Self::InvalidTerminalPolicy => "sound-invalid-terminal-policy",
            Self::InvalidAlgorithmReference => "sound-invalid-algorithm-reference",
            Self::InvalidAssumptions => "sound-invalid-assumptions",
            Self::MissingCalibration => "sound-missing-calibration",
            Self::UnexpectedCalibration => "sound-unexpected-calibration",
        }
    }
}

impl fmt::Display for SoundSynthesisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {self:?}", self.code())
    }
}

impl std::error::Error for SoundSynthesisError {}

/// Full names used by the cinematic sound specification.
#[must_use]
pub const fn synthesis_class_code(authority: SoundAuthority) -> &'static str {
    match authority {
        SoundAuthority::Artistic => "artistic-sonification",
        SoundAuthority::PhysicallyInformed => "physically-informed-modal-synthesis",
        SoundAuthority::Calibrated => "calibrated-structural-acoustics",
    }
}

fn validate(input: &SoundSynthesisInput) -> Result<(), SoundSynthesisError> {
    if input.schema_version != SOUND_SYNTHESIS_SCHEMA_VERSION {
        return Err(SoundSynthesisError::UnsupportedSchemaVersion(
            input.schema_version,
        ));
    }
    check_role(input.trajectory, CinematicComponentRole::Trajectory)?;
    check_role(input.excitation, CinematicComponentRole::AudioExcitation)?;
    check_role(input.sound_model, CinematicComponentRole::SoundModel)?;
    check_role(input.microphone, CinematicComponentRole::Microphone)?;
    check_role(input.room, CinematicComponentRole::Room)?;
    check_role(input.timeline, CinematicComponentRole::Timeline)?;

    if input.audio_clock.domain() != CinematicClockDomain::Audio
        || input.audio_clock.ticks_per_second_numerator() != SOUND_MASTER_SAMPLE_RATE_HZ
        || input.audio_clock.ticks_per_second_denominator() != 1
    {
        return Err(SoundSynthesisError::InvalidAudioClock);
    }
    if input.video_clock.domain() != CinematicClockDomain::Video
        || input.video_clock.ticks_per_second_numerator() != SOUND_MASTER_VIDEO_RATE_HZ
        || input.video_clock.ticks_per_second_denominator() != 1
    {
        return Err(SoundSynthesisError::InvalidVideoClock);
    }
    if !same_rational_instant(
        input.video_clock.start_tick(),
        input.video_clock.ticks_per_second_numerator(),
        input.video_clock.ticks_per_second_denominator(),
        input.audio_clock.start_tick(),
        input.audio_clock.ticks_per_second_numerator(),
        input.audio_clock.ticks_per_second_denominator(),
    ) || !same_rational_instant(
        input.video_clock.end_tick_exclusive(),
        input.video_clock.ticks_per_second_numerator(),
        input.video_clock.ticks_per_second_denominator(),
        input.audio_clock.end_tick_exclusive(),
        input.audio_clock.ticks_per_second_numerator(),
        input.audio_clock.ticks_per_second_denominator(),
    ) {
        return Err(SoundSynthesisError::AudioVideoTimelineMismatch);
    }
    if input.channel_layout != SoundChannelLayout::Stereo {
        return Err(SoundSynthesisError::InvalidChannelLayout);
    }
    validate_listener(input.listener)?;

    let informed = input.authority != SoundAuthority::Artistic;
    if informed && input.excitation_controls.is_empty() {
        return Err(SoundSynthesisError::MissingExcitationChannels);
    }
    if !input
        .excitation_controls
        .windows(2)
        .all(|pair| pair[0].channel < pair[1].channel)
        || input
            .excitation_controls
            .iter()
            .any(|control| !control.source_scale.is_finite() || control.source_scale == 0.0)
    {
        return Err(SoundSynthesisError::NonCanonicalExcitationChannels);
    }
    if informed && input.modes.is_empty() {
        return Err(SoundSynthesisError::MissingModes);
    }
    if input.modes.len() > MAX_SOUND_MODES {
        return Err(SoundSynthesisError::TooManyModes);
    }
    validate_modes(&input.modes)?;

    if let SoundRoomResponse::DeclaredImpulseResponse { wet_mix } = input.room_response
        && (!wet_mix.is_finite() || !(0.0..=1.0).contains(&wet_mix))
    {
        return Err(SoundSynthesisError::InvalidRoomResponse);
    }
    validate_amplitude(input.authority, input.amplitude_reference)?;
    validate_terminal(input.trajectory_disposition, input.terminal_policy)?;
    if is_zero(input.resampler_identity)
        || is_zero(input.filter_identity)
        || input.resampler_version == 0
        || input.filter_version == 0
    {
        return Err(SoundSynthesisError::InvalidAlgorithmReference);
    }
    if !strictly_increasing(&input.assumptions)
        || (informed
            && ![
                SoundModelAssumption::LinearModalSuperposition,
                SoundModelAssumption::TimeInvariantDamping,
                SoundModelAssumption::DeclaredExcitationCompleteness,
                SoundModelAssumption::DeclaredRoomResponse,
            ]
            .iter()
            .all(|required| input.assumptions.binary_search(required).is_ok()))
    {
        return Err(SoundSynthesisError::InvalidAssumptions);
    }
    match (input.authority, input.calibration) {
        (SoundAuthority::Calibrated, None) => Err(SoundSynthesisError::MissingCalibration),
        (SoundAuthority::Artistic | SoundAuthority::PhysicallyInformed, Some(_)) => {
            Err(SoundSynthesisError::UnexpectedCalibration)
        }
        _ => Ok(()),
    }
}

fn check_role(
    component: CinematicComponentRef,
    expected: CinematicComponentRole,
) -> Result<(), SoundSynthesisError> {
    if component.role() == expected {
        Ok(())
    } else {
        Err(SoundSynthesisError::WrongComponentRole {
            expected,
            got: component.role(),
        })
    }
}

fn same_rational_instant(
    a_tick: i64,
    a_num: u32,
    a_den: u32,
    b_tick: i64,
    b_num: u32,
    b_den: u32,
) -> bool {
    i128::from(a_tick) * i128::from(a_den) * i128::from(b_num)
        == i128::from(b_tick) * i128::from(b_den) * i128::from(a_num)
}

fn validate_listener(listener: ListenerPose) -> Result<(), SoundSynthesisError> {
    if listener.frame != ListenerFrame::AnimatedCamera {
        return Err(SoundSynthesisError::InvalidListenerFrame);
    }
    if listener
        .position_m
        .iter()
        .chain(listener.forward.iter())
        .chain(listener.up.iter())
        .any(|value| !value.is_finite())
    {
        return Err(SoundSynthesisError::InvalidListenerPose);
    }
    let forward_norm = dot(listener.forward, listener.forward);
    let up_norm = dot(listener.up, listener.up);
    let orthogonality = dot(listener.forward, listener.up);
    if (forward_norm - 1.0).abs() > 1.0e-12
        || (up_norm - 1.0).abs() > 1.0e-12
        || orthogonality.abs() > 1.0e-12
    {
        return Err(SoundSynthesisError::InvalidListenerPose);
    }
    Ok(())
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn strictly_increasing<T: Ord>(values: &[T]) -> bool {
    values.windows(2).all(|pair| pair[0] < pair[1])
}

fn validate_modes(modes: &[SoundMode]) -> Result<(), SoundSynthesisError> {
    let mut previous_id = None;
    for mode in modes {
        if mode.mode_id == 0
            || previous_id.is_some_and(|previous| mode.mode_id <= previous)
            || !mode.frequency_hz.is_finite()
            || mode.frequency_hz <= 0.0
            || mode.frequency_hz >= f64::from(SOUND_MASTER_SAMPLE_RATE_HZ) / 2.0
            || !mode.damping_ratio.is_finite()
            || mode.damping_ratio <= 0.0
            || mode.damping_ratio > 1.0
            || !mode.gain.is_finite()
            || is_zero(mode.material_identity)
            || is_zero(mode.base_identity)
        {
            return Err(SoundSynthesisError::InvalidMode);
        }
        previous_id = Some(mode.mode_id);
    }
    Ok(())
}

fn validate_amplitude(
    authority: SoundAuthority,
    amplitude: SoundAmplitudeReference,
) -> Result<(), SoundSynthesisError> {
    let (headroom_db, pressure) = match amplitude {
        SoundAmplitudeReference::DigitalFullScale { headroom_db } => (headroom_db, None),
        SoundAmplitudeReference::CalibratedPressure {
            pascal_rms_at_full_scale,
            headroom_db,
        } => (headroom_db, Some(pascal_rms_at_full_scale)),
    };
    if !headroom_db.is_finite() || !(0.0..=60.0).contains(&headroom_db) {
        return Err(SoundSynthesisError::InvalidAmplitudeReference);
    }
    match (authority, pressure) {
        (SoundAuthority::Calibrated, Some(value)) if value.is_finite() && value > 0.0 => Ok(()),
        (SoundAuthority::Calibrated, _) | (_, Some(_)) => {
            Err(SoundSynthesisError::InvalidAmplitudeReference)
        }
        (_, None) => Ok(()),
    }
}

fn validate_terminal(
    disposition: SoundTrajectoryDisposition,
    policy: SoundTerminalPolicy,
) -> Result<(), SoundSynthesisError> {
    match (disposition, policy) {
        (SoundTrajectoryDisposition::NumericalRefusal, SoundTerminalPolicy::Silence) => Ok(()),
        (
            SoundTrajectoryDisposition::PhysicalTerminal
            | SoundTrajectoryDisposition::HorizonCensored,
            SoundTerminalPolicy::FadeAtLastAccepted { fade_sample_frames },
        ) if (1..=SOUND_MASTER_SAMPLE_RATE_HZ).contains(&fade_sample_frames) => Ok(()),
        _ => Err(SoundSynthesisError::InvalidTerminalPolicy),
    }
}

fn is_zero(identity: ContentHash) -> bool {
    identity.as_bytes().iter().all(|byte| *byte == 0)
}

fn sound_identity(input: &SoundSynthesisInput) -> ContentHash {
    let mut bytes = Vec::with_capacity(1024 + input.modes.len() * 100);
    push_u16(&mut bytes, input.schema_version);
    bytes.push(match input.authority {
        SoundAuthority::Artistic => 1,
        SoundAuthority::PhysicallyInformed => 2,
        SoundAuthority::Calibrated => 3,
    });
    for component in [
        input.trajectory,
        input.excitation,
        input.sound_model,
        input.microphone,
        input.room,
        input.timeline,
    ] {
        bytes.push(component.role() as u8);
        bytes.extend_from_slice(component.identity().as_bytes());
        push_u32(&mut bytes, component.version());
    }
    push_clock(&mut bytes, input.video_clock);
    push_clock(&mut bytes, input.audio_clock);
    bytes.push(input.channel_layout as u8);
    bytes.push(input.listener.frame as u8);
    for vector in [
        input.listener.position_m,
        input.listener.forward,
        input.listener.up,
    ] {
        for value in vector {
            push_f64(&mut bytes, value);
        }
    }
    push_u32(&mut bytes, input.excitation_controls.len() as u32);
    for control in &input.excitation_controls {
        bytes.push(control.channel as u8);
        push_f64(&mut bytes, control.source_scale);
    }
    push_u32(&mut bytes, input.modes.len() as u32);
    for mode in &input.modes {
        push_u32(&mut bytes, mode.mode_id);
        push_f64(&mut bytes, mode.frequency_hz);
        push_f64(&mut bytes, mode.damping_ratio);
        push_f64(&mut bytes, mode.gain);
        bytes.extend_from_slice(mode.material_identity.as_bytes());
        bytes.extend_from_slice(mode.base_identity.as_bytes());
    }
    match input.room_response {
        SoundRoomResponse::Dry => bytes.push(1),
        SoundRoomResponse::DeclaredImpulseResponse { wet_mix } => {
            bytes.push(2);
            push_f64(&mut bytes, wet_mix);
        }
    }
    match input.amplitude_reference {
        SoundAmplitudeReference::DigitalFullScale { headroom_db } => {
            bytes.push(1);
            push_f64(&mut bytes, headroom_db);
        }
        SoundAmplitudeReference::CalibratedPressure {
            pascal_rms_at_full_scale,
            headroom_db,
        } => {
            bytes.push(2);
            push_f64(&mut bytes, pascal_rms_at_full_scale);
            push_f64(&mut bytes, headroom_db);
        }
    }
    bytes.push(input.trajectory_disposition as u8);
    match input.terminal_policy {
        SoundTerminalPolicy::FadeAtLastAccepted { fade_sample_frames } => {
            bytes.push(1);
            push_u32(&mut bytes, fade_sample_frames);
        }
        SoundTerminalPolicy::Silence => bytes.push(2),
    }
    bytes.extend_from_slice(input.resampler_identity.as_bytes());
    push_u32(&mut bytes, input.resampler_version);
    bytes.extend_from_slice(input.filter_identity.as_bytes());
    push_u32(&mut bytes, input.filter_version);
    push_u32(&mut bytes, input.assumptions.len() as u32);
    bytes.extend(input.assumptions.iter().map(|assumption| *assumption as u8));
    match input.calibration {
        None => bytes.push(0),
        Some(calibration) => {
            bytes.push(1);
            bytes.extend_from_slice(calibration.dataset_identity().as_bytes());
            bytes.extend_from_slice(calibration.method_identity().as_bytes());
            bytes.extend_from_slice(calibration.validity_identity().as_bytes());
            push_u32(&mut bytes, calibration.version());
        }
    }
    hash_domain(IDENTITY_DOMAIN, &bytes)
}

fn push_clock(bytes: &mut Vec<u8>, clock: CinematicClock) {
    bytes.push(match clock.domain() {
        CinematicClockDomain::Simulation => 1,
        CinematicClockDomain::Video => 2,
        CinematicClockDomain::Audio => 3,
        CinematicClockDomain::Composition => 4,
        CinematicClockDomain::Timeless => 5,
    });
    push_u32(bytes, clock.ticks_per_second_numerator());
    push_u32(bytes, clock.ticks_per_second_denominator());
    bytes.extend_from_slice(&clock.start_tick().to_le_bytes());
    bytes.extend_from_slice(&clock.end_tick_exclusive().to_le_bytes());
}

fn push_u16(bytes: &mut Vec<u8>, value: u16) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_f64(bytes: &mut Vec<u8>, value: f64) {
    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
}
