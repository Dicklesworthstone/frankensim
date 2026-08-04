//! Authority and disclosure contract for simulation-driven cinematic artifacts.
//!
//! The types in this module keep four evidence classes distinct: simulated
//! state, raw Monte Carlo render estimates, biased visualization derivatives,
//! and sound. They provide structural admission and deterministic identities;
//! they do not authenticate a producer or validate a physical model.

use core::fmt;
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};

/// Current binary/manifest schema version.
pub const CINEMATIC_AUTHORITY_SCHEMA_VERSION: u16 = 1;

/// Maximum bytes in a machine-readable transform or synthesis label.
pub const MAX_CINEMATIC_LABEL_BYTES: usize = 128;

/// Maximum number of explicit no-claim declarations on one artifact.
pub const MAX_CINEMATIC_NO_CLAIMS: usize = 32;

const MAGIC: &[u8; 8] = b"FSCINAU1";
const IDENTITY_DOMAIN: &str = "org.frankensim.cinematic-authority-record.v1";

/// Current canonical Euler cinematic deliverable contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CinematicDeliverableContract {
    width_pixels: u32,
    height_pixels: u32,
    frames_per_second_numerator: u32,
    frames_per_second_denominator: u32,
    minimum_frame_count: u32,
    maximum_frame_count: u32,
    image_master: CinematicImageMaster,
    display_preview: CinematicDisplayPreview,
    audio_sample_rate_hz: u32,
    audio_channels: u8,
    audio_master: CinematicAudioMaster,
    sequence_manifest_required: bool,
    muxed_derivative: CinematicMuxedDerivative,
}

impl CinematicDeliverableContract {
    /// The frozen first Euler deliverable: 4K UHD, 24 fps, 8–12 seconds,
    /// float EXR image masters, display previews, stereo 48 kHz float WAV,
    /// a required sequence manifest, and an optional non-authoritative mux.
    #[must_use]
    pub const fn euler_disc_v1() -> Self {
        Self {
            width_pixels: 3_840,
            height_pixels: 2_160,
            frames_per_second_numerator: 24,
            frames_per_second_denominator: 1,
            minimum_frame_count: 8 * 24,
            maximum_frame_count: 12 * 24,
            image_master: CinematicImageMaster::OpenExrFloat,
            display_preview: CinematicDisplayPreview::DisplayReferred,
            audio_sample_rate_hz: 48_000,
            audio_channels: 2,
            audio_master: CinematicAudioMaster::WaveFloat32,
            sequence_manifest_required: true,
            muxed_derivative: CinematicMuxedDerivative::OptionalNonAuthoritative,
        }
    }

    /// Width in pixels.
    #[must_use]
    pub const fn width_pixels(self) -> u32 {
        self.width_pixels
    }

    /// Height in pixels.
    #[must_use]
    pub const fn height_pixels(self) -> u32 {
        self.height_pixels
    }

    /// Exact frames-per-second numerator.
    #[must_use]
    pub const fn frames_per_second_numerator(self) -> u32 {
        self.frames_per_second_numerator
    }

    /// Exact frames-per-second denominator.
    #[must_use]
    pub const fn frames_per_second_denominator(self) -> u32 {
        self.frames_per_second_denominator
    }

    /// Inclusive minimum number of frames.
    #[must_use]
    pub const fn minimum_frame_count(self) -> u32 {
        self.minimum_frame_count
    }

    /// Inclusive maximum number of frames.
    #[must_use]
    pub const fn maximum_frame_count(self) -> u32 {
        self.maximum_frame_count
    }

    /// Image master format.
    #[must_use]
    pub const fn image_master(self) -> CinematicImageMaster {
        self.image_master
    }

    /// Display preview policy.
    #[must_use]
    pub const fn display_preview(self) -> CinematicDisplayPreview {
        self.display_preview
    }

    /// Audio sample rate in Hz.
    #[must_use]
    pub const fn audio_sample_rate_hz(self) -> u32 {
        self.audio_sample_rate_hz
    }

    /// Number of interleaved audio channels.
    #[must_use]
    pub const fn audio_channels(self) -> u8 {
        self.audio_channels
    }

    /// Audio master format.
    #[must_use]
    pub const fn audio_master(self) -> CinematicAudioMaster {
        self.audio_master
    }

    /// Whether the image/audio sequence manifest is mandatory.
    #[must_use]
    pub const fn sequence_manifest_required(self) -> bool {
        self.sequence_manifest_required
    }

    /// Muxed media policy.
    #[must_use]
    pub const fn muxed_derivative(self) -> CinematicMuxedDerivative {
        self.muxed_derivative
    }

    /// Validate an exact frame/sample timeline against the frozen duration
    /// and synchronization contract. `audio_sample_frames` counts stereo
    /// frames, not individual interleaved scalar samples.
    pub fn validate_timeline(
        self,
        video_frame_count: u32,
        audio_sample_frames: u64,
    ) -> Result<(), CinematicDeliverableError> {
        if !(self.minimum_frame_count..=self.maximum_frame_count).contains(&video_frame_count) {
            return Err(CinematicDeliverableError::FrameCountOutOfRange {
                got: video_frame_count,
                minimum: self.minimum_frame_count,
                maximum: self.maximum_frame_count,
            });
        }
        let expected_audio_sample_frames = u64::from(video_frame_count)
            * u64::from(self.audio_sample_rate_hz)
            * u64::from(self.frames_per_second_denominator)
            / u64::from(self.frames_per_second_numerator);
        if audio_sample_frames != expected_audio_sample_frames {
            return Err(CinematicDeliverableError::AudioVideoClockMismatch {
                video_frame_count,
                audio_sample_frames,
                expected_audio_sample_frames,
            });
        }
        Ok(())
    }

    /// Deterministic machine-readable declaration of the frozen envelope.
    #[must_use]
    pub fn to_manifest_json(self) -> String {
        format!(
            "{{\"contract\":\"euler-disc-cinematic-v1\",\"width_pixels\":{},\"height_pixels\":{},\"frames_per_second_numerator\":{},\"frames_per_second_denominator\":{},\"minimum_frame_count\":{},\"maximum_frame_count\":{},\"image_master\":\"{}\",\"display_preview\":\"{}\",\"audio_sample_rate_hz\":{},\"audio_channels\":{},\"audio_master\":\"{}\",\"sequence_manifest_required\":{},\"muxed_derivative\":\"{}\"}}",
            self.width_pixels,
            self.height_pixels,
            self.frames_per_second_numerator,
            self.frames_per_second_denominator,
            self.minimum_frame_count,
            self.maximum_frame_count,
            self.image_master.code(),
            self.display_preview.code(),
            self.audio_sample_rate_hz,
            self.audio_channels,
            self.audio_master.code(),
            self.sequence_manifest_required,
            self.muxed_derivative.code(),
        )
    }
}

/// Image-master representation required by a deliverable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CinematicImageMaster {
    /// Per-frame floating-point OpenEXR master.
    OpenExrFloat,
}

impl CinematicImageMaster {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::OpenExrFloat => "openexr-float",
        }
    }
}

/// Display-preview representation required by a deliverable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CinematicDisplayPreview {
    /// Display-referred preview derived from, and never replacing, the master.
    DisplayReferred,
}

impl CinematicDisplayPreview {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::DisplayReferred => "display-referred-preview",
        }
    }
}

/// Audio-master representation required by a deliverable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CinematicAudioMaster {
    /// IEEE 32-bit floating-point samples in a RIFF/WAVE container.
    WaveFloat32,
}

impl CinematicAudioMaster {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::WaveFloat32 => "wave-float32",
        }
    }
}

/// Policy for a user-convenience muxed movie.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CinematicMuxedDerivative {
    /// Optional and explicitly subordinate to the manifest-bound masters.
    OptionalNonAuthoritative,
}

impl CinematicMuxedDerivative {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::OptionalNonAuthoritative => "optional-non-authoritative",
        }
    }
}

/// Stable meanings for terms that must not drift across the pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicTerm {
    /// Numerical behavior of the raw Monte Carlo image estimator.
    RenderConvergence,
    /// Reproduction under the same declared inputs and determinism envelope.
    DeterministicReplay,
    /// Aesthetic acceptance without scientific promotion.
    VisualApproval,
    /// Parameter fitting plus a declared evidence and validity domain.
    ModelCalibration,
    /// Independent comparison to empirical observations in a declared domain.
    ExperimentalValidation,
    /// Container/codec transformation of admitted masters.
    MediaEncoding,
}

impl CinematicTerm {
    /// Every frozen term in canonical order.
    pub const ALL: [Self; 6] = [
        Self::RenderConvergence,
        Self::DeterministicReplay,
        Self::VisualApproval,
        Self::ModelCalibration,
        Self::ExperimentalValidation,
        Self::MediaEncoding,
    ];

    /// Stable machine code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::RenderConvergence => "render-convergence",
            Self::DeterministicReplay => "deterministic-replay",
            Self::VisualApproval => "visual-approval",
            Self::ModelCalibration => "model-calibration",
            Self::ExperimentalValidation => "experimental-validation",
            Self::MediaEncoding => "media-encoding",
        }
    }

    /// Normative human definition.
    #[must_use]
    pub const fn definition(self) -> &'static str {
        match self {
            Self::RenderConvergence => {
                "Estimator diagnostics meet a declared image-space tolerance; this says nothing about physical-model validity."
            }
            Self::DeterministicReplay => {
                "The same admitted inputs, configuration, and declared execution envelope reproduce the same artifact identity."
            }
            Self::VisualApproval => {
                "A reviewer accepts the aesthetic result; approval cannot promote scientific or calibration authority."
            }
            Self::ModelCalibration => {
                "A declared procedure fitted model parameters to named evidence inside an explicit validity domain."
            }
            Self::ExperimentalValidation => {
                "An independent comparison established a stated physical claim against empirical observations in a declared domain."
            }
            Self::MediaEncoding => {
                "A container or codec transformed admitted masters into a convenience derivative without increasing their authority."
            }
        }
    }
}

/// Refusals from the frozen A/V timeline contract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicDeliverableError {
    /// Frame count lies outside the inclusive 8–12 second envelope.
    FrameCountOutOfRange {
        /// Supplied video frames.
        got: u32,
        /// Inclusive lower bound.
        minimum: u32,
        /// Inclusive upper bound.
        maximum: u32,
    },
    /// Audio sample-frame count does not end at the same exact time as video.
    AudioVideoClockMismatch {
        /// Supplied video frames.
        video_frame_count: u32,
        /// Supplied audio sample frames.
        audio_sample_frames: u64,
        /// Required audio sample frames for exact synchronization.
        expected_audio_sample_frames: u64,
    },
}

impl CinematicDeliverableError {
    /// Stable reason code for CLI and finalization diagnostics.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::FrameCountOutOfRange { .. } => "cinematic-frame-count-out-of-range",
            Self::AudioVideoClockMismatch { .. } => "cinematic-audio-video-clock-mismatch",
        }
    }
}

impl fmt::Display for CinematicDeliverableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {self:?}", self.code())
    }
}

impl std::error::Error for CinematicDeliverableError {}

/// Sound authority is separate from image/render authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SoundAuthority {
    /// Sound designed for communication or aesthetics without a physics claim.
    Artistic,
    /// Sound driven by simulation channels and declared synthesis models.
    PhysicallyInformed,
    /// Sound tied to an explicit, separately supplied calibration receipt.
    Calibrated,
}

impl SoundAuthority {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Artistic => "artistic",
            Self::PhysicallyInformed => "physically-informed",
            Self::Calibrated => "calibrated",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Artistic => 1,
            Self::PhysicallyInformed => 2,
            Self::Calibrated => 3,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, CinematicAuthorityError> {
        match tag {
            1 => Ok(Self::Artistic),
            2 => Ok(Self::PhysicallyInformed),
            3 => Ok(Self::Calibrated),
            _ => Err(CinematicAuthorityError::UnknownSoundAuthorityTag(tag)),
        }
    }
}

/// Non-interchangeable authority carried by a cinematic artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicAuthorityClass {
    /// Accepted simulation/model state. This is not experimental validation.
    SimulatedState,
    /// A raw Monte Carlo estimator prior to biased finishing.
    MonteCarloRender,
    /// A denoised, tone-mapped, composited, or encoded visual derivative.
    VisualizationDerivative,
    /// An audio artifact with its own authority tier.
    Sound(SoundAuthority),
}

impl CinematicAuthorityClass {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::SimulatedState => "simulated-state",
            Self::MonteCarloRender => "monte-carlo-render",
            Self::VisualizationDerivative => "visualization-derivative",
            Self::Sound(_) => "sound",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::SimulatedState => 1,
            Self::MonteCarloRender => 2,
            Self::VisualizationDerivative => 3,
            Self::Sound(_) => 4,
        }
    }

    const fn sound_tag(self) -> u8 {
        match self {
            Self::Sound(sound) => sound.tag(),
            _ => 0,
        }
    }

    fn from_tags(class: u8, sound: u8) -> Result<Self, CinematicAuthorityError> {
        match class {
            1 if sound == 0 => Ok(Self::SimulatedState),
            2 if sound == 0 => Ok(Self::MonteCarloRender),
            3 if sound == 0 => Ok(Self::VisualizationDerivative),
            4 => Ok(Self::Sound(SoundAuthority::from_tag(sound)?)),
            1..=3 => Err(CinematicAuthorityError::UnexpectedSoundAuthorityTag(sound)),
            _ => Err(CinematicAuthorityError::UnknownAuthorityClassTag(class)),
        }
    }
}

/// Concrete payload family described by an authority envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicArtifactKind {
    /// Time-indexed or instantaneous mechanics state.
    SimulationState,
    /// Raw radiance/film/AOV estimate.
    RenderEstimate,
    /// Denoised, display-referred, composited, or overlaid image.
    Visualization,
    /// Digital audio samples or an audio stem.
    Audio,
    /// A manifest/log describing payloads of the declared authority class.
    ManifestOrLog,
    /// A playable visual media derivative; component audio retains its record.
    MediaDerivative,
}

impl CinematicArtifactKind {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::SimulationState => "simulation-state",
            Self::RenderEstimate => "render-estimate",
            Self::Visualization => "visualization",
            Self::Audio => "audio",
            Self::ManifestOrLog => "manifest-or-log",
            Self::MediaDerivative => "media-derivative",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::SimulationState => 1,
            Self::RenderEstimate => 2,
            Self::Visualization => 3,
            Self::Audio => 4,
            Self::ManifestOrLog => 5,
            Self::MediaDerivative => 6,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, CinematicAuthorityError> {
        match tag {
            1 => Ok(Self::SimulationState),
            2 => Ok(Self::RenderEstimate),
            3 => Ok(Self::Visualization),
            4 => Ok(Self::Audio),
            5 => Ok(Self::ManifestOrLog),
            6 => Ok(Self::MediaDerivative),
            _ => Err(CinematicAuthorityError::UnknownArtifactKindTag(tag)),
        }
    }
}

/// Explicit unit contract for the payload values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicUnitContract {
    /// SI mechanics with radians for angular quantities.
    SiMechanicsRadians,
    /// Spectral radiance/film values under the renderer's declared convention.
    SpectralRadianceSi,
    /// Display-referred normalized code values under a declared color transform.
    DisplayEncoded,
    /// Digital audio amplitude relative to full scale.
    DigitalAudioFullScale,
    /// A manifest/log whose children carry their own unit contracts.
    ChildDeclared,
}

impl CinematicUnitContract {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::SiMechanicsRadians => "si-mechanics-radians-v1",
            Self::SpectralRadianceSi => "spectral-radiance-si-v1",
            Self::DisplayEncoded => "display-encoded-v1",
            Self::DigitalAudioFullScale => "digital-audio-full-scale-v1",
            Self::ChildDeclared => "child-declared-v1",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::SiMechanicsRadians => 1,
            Self::SpectralRadianceSi => 2,
            Self::DisplayEncoded => 3,
            Self::DigitalAudioFullScale => 4,
            Self::ChildDeclared => 5,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, CinematicAuthorityError> {
        match tag {
            1 => Ok(Self::SiMechanicsRadians),
            2 => Ok(Self::SpectralRadianceSi),
            3 => Ok(Self::DisplayEncoded),
            4 => Ok(Self::DigitalAudioFullScale),
            5 => Ok(Self::ChildDeclared),
            _ => Err(CinematicAuthorityError::UnknownUnitContractTag(tag)),
        }
    }
}

/// Clock namespace for the exact tick range carried by an artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicClockDomain {
    /// Simulation time.
    Simulation,
    /// Video-frame timeline.
    Video,
    /// Digital-audio sample timeline.
    Audio,
    /// Composition timeline shared by cuts, frames, and samples.
    Composition,
    /// Explicitly timeless metadata.
    Timeless,
}

impl CinematicClockDomain {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Simulation => "simulation",
            Self::Video => "video",
            Self::Audio => "audio",
            Self::Composition => "composition",
            Self::Timeless => "timeless",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Simulation => 1,
            Self::Video => 2,
            Self::Audio => 3,
            Self::Composition => 4,
            Self::Timeless => 5,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, CinematicAuthorityError> {
        match tag {
            1 => Ok(Self::Simulation),
            2 => Ok(Self::Video),
            3 => Ok(Self::Audio),
            4 => Ok(Self::Composition),
            5 => Ok(Self::Timeless),
            _ => Err(CinematicAuthorityError::UnknownClockDomainTag(tag)),
        }
    }
}

/// Exact rational clock and half-open tick range.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CinematicClock {
    domain: CinematicClockDomain,
    ticks_per_second_numerator: u32,
    ticks_per_second_denominator: u32,
    start_tick: i64,
    end_tick_exclusive: i64,
}

impl CinematicClock {
    /// Construct a checked clock. Instants use equal start/end ticks.
    pub fn try_new(
        domain: CinematicClockDomain,
        ticks_per_second_numerator: u32,
        ticks_per_second_denominator: u32,
        start_tick: i64,
        end_tick_exclusive: i64,
    ) -> Result<Self, CinematicAuthorityError> {
        if ticks_per_second_numerator == 0 || ticks_per_second_denominator == 0 {
            return Err(CinematicAuthorityError::InvalidClockRate);
        }
        if end_tick_exclusive < start_tick {
            return Err(CinematicAuthorityError::InvalidClockRange);
        }
        if domain == CinematicClockDomain::Timeless
            && (ticks_per_second_numerator != 1
                || ticks_per_second_denominator != 1
                || start_tick != 0
                || end_tick_exclusive != 0)
        {
            return Err(CinematicAuthorityError::InvalidTimelessClock);
        }
        Ok(Self {
            domain,
            ticks_per_second_numerator,
            ticks_per_second_denominator,
            start_tick,
            end_tick_exclusive,
        })
    }

    /// Explicit timeless metadata clock.
    #[must_use]
    pub const fn timeless() -> Self {
        Self {
            domain: CinematicClockDomain::Timeless,
            ticks_per_second_numerator: 1,
            ticks_per_second_denominator: 1,
            start_tick: 0,
            end_tick_exclusive: 0,
        }
    }

    /// Clock domain.
    #[must_use]
    pub const fn domain(self) -> CinematicClockDomain {
        self.domain
    }

    /// Rational ticks-per-second numerator.
    #[must_use]
    pub const fn ticks_per_second_numerator(self) -> u32 {
        self.ticks_per_second_numerator
    }

    /// Rational ticks-per-second denominator.
    #[must_use]
    pub const fn ticks_per_second_denominator(self) -> u32 {
        self.ticks_per_second_denominator
    }

    /// Inclusive first tick.
    #[must_use]
    pub const fn start_tick(self) -> i64 {
        self.start_tick
    }

    /// Exclusive last tick, or the same tick for an instant.
    #[must_use]
    pub const fn end_tick_exclusive(self) -> i64 {
        self.end_tick_exclusive
    }
}

/// Transform disposition prevents a derived artifact from masquerading as raw.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum CinematicTransformDisposition {
    /// Accepted model state without a render estimator.
    ModelState,
    /// Raw Monte Carlo estimator.
    MonteCarloEstimator,
    /// Explicitly biased image/display operation with a machine label.
    BiasedVisualization(String),
    /// Explicit sound synthesis operation with a machine label.
    SoundSynthesis(String),
}

impl CinematicTransformDisposition {
    /// Stable manifest code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::ModelState => "model-state",
            Self::MonteCarloEstimator => "monte-carlo-estimator",
            Self::BiasedVisualization(_) => "biased-visualization",
            Self::SoundSynthesis(_) => "sound-synthesis",
        }
    }

    /// Optional algorithm label for biased/synthesized transforms.
    #[must_use]
    pub fn label(&self) -> Option<&str> {
        match self {
            Self::BiasedVisualization(label) | Self::SoundSynthesis(label) => Some(label),
            Self::ModelState | Self::MonteCarloEstimator => None,
        }
    }

    const fn tag(&self) -> u8 {
        match self {
            Self::ModelState => 1,
            Self::MonteCarloEstimator => 2,
            Self::BiasedVisualization(_) => 3,
            Self::SoundSynthesis(_) => 4,
        }
    }
}

/// Claims that the Euler cinematic is explicitly forbidden to imply.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicNoClaim {
    /// No experiment-backed physical stopping-time prediction.
    PhysicalStoppingTime,
    /// No validated ring/cone/disc ranking.
    RingConeRanking,
    /// No experimental validation of the contact/rolling law.
    ContactLawValidation,
    /// No physical prediction of a terminal finite-time singularity.
    TerminalSingularityPrediction,
    /// No validation against the motivating Steve Mould video.
    MouldVideoCorrespondence,
    /// No calibrated acoustic waveform or SPL prediction.
    CalibratedAcousticPrediction,
    /// A biased image cannot become raw numerical/render evidence.
    RawEvidenceFromBiasedArtifact,
}

impl CinematicNoClaim {
    /// Stable manifest reason code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::PhysicalStoppingTime => "no-physical-stopping-time",
            Self::RingConeRanking => "no-validated-ring-cone-ranking",
            Self::ContactLawValidation => "no-contact-law-validation",
            Self::TerminalSingularityPrediction => "no-terminal-singularity-prediction",
            Self::MouldVideoCorrespondence => "no-mould-video-validation",
            Self::CalibratedAcousticPrediction => "no-calibrated-acoustic-prediction",
            Self::RawEvidenceFromBiasedArtifact => "no-raw-evidence-from-biased-artifact",
        }
    }

    /// Human-facing disclosure generated from the same closed enum.
    #[must_use]
    pub const fn statement(self) -> &'static str {
        match self {
            Self::PhysicalStoppingTime => {
                "This artifact does not validate or predict physical stopping time."
            }
            Self::RingConeRanking => {
                "This artifact does not validate a ring, cone, or disc performance ranking."
            }
            Self::ContactLawValidation => {
                "This artifact does not experimentally validate the contact or rolling law."
            }
            Self::TerminalSingularityPrediction => {
                "This artifact does not predict a physical terminal finite-time singularity."
            }
            Self::MouldVideoCorrespondence => {
                "Visual or audible resemblance does not validate correspondence with the Steve Mould video."
            }
            Self::CalibratedAcousticPrediction => {
                "This sound is not a calibrated acoustic waveform or sound-pressure prediction."
            }
            Self::RawEvidenceFromBiasedArtifact => {
                "A biased visualization derivative is not raw numerical or Monte Carlo evidence."
            }
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::PhysicalStoppingTime => 1,
            Self::RingConeRanking => 2,
            Self::ContactLawValidation => 3,
            Self::TerminalSingularityPrediction => 4,
            Self::MouldVideoCorrespondence => 5,
            Self::CalibratedAcousticPrediction => 6,
            Self::RawEvidenceFromBiasedArtifact => 7,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, CinematicAuthorityError> {
        match tag {
            1 => Ok(Self::PhysicalStoppingTime),
            2 => Ok(Self::RingConeRanking),
            3 => Ok(Self::ContactLawValidation),
            4 => Ok(Self::TerminalSingularityPrediction),
            5 => Ok(Self::MouldVideoCorrespondence),
            6 => Ok(Self::CalibratedAcousticPrediction),
            7 => Ok(Self::RawEvidenceFromBiasedArtifact),
            _ => Err(CinematicAuthorityError::UnknownNoClaimTag(tag)),
        }
    }
}

/// Declared external evidence needed for a new calibrated-sound record.
///
/// This is structural binding only. Authentication and in-domain validation
/// remain the responsibility of the calibration/evidence verifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DeclaredAcousticCalibrationReceipt {
    dataset_identity: ContentHash,
    method_identity: ContentHash,
    validity_identity: ContentHash,
    version: u32,
}

impl DeclaredAcousticCalibrationReceipt {
    /// Construct a structurally complete calibration declaration.
    pub fn try_new(
        dataset_identity: ContentHash,
        method_identity: ContentHash,
        validity_identity: ContentHash,
        version: u32,
    ) -> Result<Self, CinematicAuthorityError> {
        check_identity("calibration-dataset", dataset_identity)?;
        check_identity("calibration-method", method_identity)?;
        check_identity("calibration-validity", validity_identity)?;
        if version == 0 {
            return Err(CinematicAuthorityError::InvalidCalibrationVersion);
        }
        Ok(Self {
            dataset_identity,
            method_identity,
            validity_identity,
            version,
        })
    }

    /// Dataset artifact identity.
    #[must_use]
    pub const fn dataset_identity(self) -> ContentHash {
        self.dataset_identity
    }

    /// Calibration method/receipt identity.
    #[must_use]
    pub const fn method_identity(self) -> ContentHash {
        self.method_identity
    }

    /// Validity-domain identity.
    #[must_use]
    pub const fn validity_identity(self) -> ContentHash {
        self.validity_identity
    }

    /// Calibration schema/method version.
    #[must_use]
    pub const fn version(self) -> u32 {
        self.version
    }
}

/// Complete caller input. No field has an implicit default.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicAuthorityInput {
    /// Exact supported schema version.
    pub schema_version: u16,
    /// Payload family.
    pub artifact_kind: CinematicArtifactKind,
    /// Scientific/render/audio authority class.
    pub authority_class: CinematicAuthorityClass,
    /// Hash of the artifact receiving this record.
    pub artifact_identity: ContentHash,
    /// Hash of the immediate source artifact.
    pub source_identity: ContentHash,
    /// Hash of the operation/transform receipt.
    pub transform_identity: ContentHash,
    /// Stable machine-readable transform name.
    pub transform_name: String,
    /// Hash of the complete cinematic configuration.
    pub configuration_identity: ContentHash,
    /// Nonzero configuration schema version.
    pub configuration_version: u32,
    /// Explicit unit interpretation.
    pub unit_contract: CinematicUnitContract,
    /// Exact time/sample/frame clock.
    pub clock: CinematicClock,
    /// Raw/biased/synthesized disposition.
    pub transform_disposition: CinematicTransformDisposition,
    /// Applicable disclosure set. Order is canonicalized; duplicates refuse.
    pub no_claims: Vec<CinematicNoClaim>,
    /// Required only for calibrated sound.
    pub acoustic_calibration: Option<DeclaredAcousticCalibrationReceipt>,
}

/// Opaque, validated cinematic authority envelope.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicAuthorityRecord {
    schema_version: u16,
    artifact_kind: CinematicArtifactKind,
    authority_class: CinematicAuthorityClass,
    artifact_identity: ContentHash,
    source_identity: ContentHash,
    transform_identity: ContentHash,
    transform_name: String,
    configuration_identity: ContentHash,
    configuration_version: u32,
    unit_contract: CinematicUnitContract,
    clock: CinematicClock,
    transform_disposition: CinematicTransformDisposition,
    no_claims: Vec<CinematicNoClaim>,
    acoustic_calibration: Option<DeclaredAcousticCalibrationReceipt>,
}

impl CinematicAuthorityRecord {
    /// Admit one complete declaration. This validates structure and authority
    /// compatibility but does not authenticate the referenced hashes.
    pub fn try_new(mut input: CinematicAuthorityInput) -> Result<Self, CinematicAuthorityError> {
        if input.schema_version != CINEMATIC_AUTHORITY_SCHEMA_VERSION {
            return Err(CinematicAuthorityError::UnsupportedSchemaVersion(
                input.schema_version,
            ));
        }
        check_identity("artifact", input.artifact_identity)?;
        check_identity("source", input.source_identity)?;
        check_identity("transform", input.transform_identity)?;
        check_identity("configuration", input.configuration_identity)?;
        if input.configuration_version == 0 {
            return Err(CinematicAuthorityError::InvalidConfigurationVersion);
        }
        validate_label(&input.transform_name)?;
        validate_artifact_class(input.artifact_kind, input.authority_class)?;
        validate_disposition(&input.transform_disposition, input.authority_class)?;

        match (input.authority_class, input.acoustic_calibration) {
            (CinematicAuthorityClass::Sound(SoundAuthority::Calibrated), Some(_)) => {}
            (CinematicAuthorityClass::Sound(SoundAuthority::Calibrated), None) => {
                return Err(CinematicAuthorityError::MissingCalibrationReceipt);
            }
            (_, Some(_)) => return Err(CinematicAuthorityError::UnexpectedCalibrationReceipt),
            (_, None) => {}
        }

        if input.no_claims.len() > MAX_CINEMATIC_NO_CLAIMS {
            return Err(CinematicAuthorityError::TooManyNoClaims {
                got: input.no_claims.len(),
                max: MAX_CINEMATIC_NO_CLAIMS,
            });
        }
        input.no_claims.sort_unstable();
        if let Some(pair) = input.no_claims.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(CinematicAuthorityError::DuplicateNoClaim(pair[0]));
        }
        for required in required_no_claims(input.authority_class) {
            if input.no_claims.binary_search(&required).is_err() {
                return Err(CinematicAuthorityError::MissingNoClaim(required));
            }
        }

        Ok(Self {
            schema_version: input.schema_version,
            artifact_kind: input.artifact_kind,
            authority_class: input.authority_class,
            artifact_identity: input.artifact_identity,
            source_identity: input.source_identity,
            transform_identity: input.transform_identity,
            transform_name: input.transform_name,
            configuration_identity: input.configuration_identity,
            configuration_version: input.configuration_version,
            unit_contract: input.unit_contract,
            clock: input.clock,
            transform_disposition: input.transform_disposition,
            no_claims: input.no_claims,
            acoustic_calibration: input.acoustic_calibration,
        })
    }

    /// Derive a new immutable record from an admitted parent. The new record
    /// must name the parent's artifact as its immediate source.
    pub fn derive(
        parent: &Self,
        input: CinematicAuthorityInput,
    ) -> Result<Self, CinematicAuthorityError> {
        if input.source_identity != parent.artifact_identity {
            return Err(CinematicAuthorityError::SourceDoesNotMatchParent);
        }
        if !transition_allowed(parent.authority_class, input.authority_class) {
            return Err(CinematicAuthorityError::IllegalAuthorityTransition {
                from: parent.authority_class,
                to: input.authority_class,
            });
        }
        Self::try_new(input)
    }

    /// Schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u16 {
        self.schema_version
    }

    /// Payload family.
    #[must_use]
    pub const fn artifact_kind(&self) -> CinematicArtifactKind {
        self.artifact_kind
    }

    /// Authority class.
    #[must_use]
    pub const fn authority_class(&self) -> CinematicAuthorityClass {
        self.authority_class
    }

    /// Artifact identity.
    #[must_use]
    pub const fn artifact_identity(&self) -> ContentHash {
        self.artifact_identity
    }

    /// Immediate source identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Transform receipt identity.
    #[must_use]
    pub const fn transform_identity(&self) -> ContentHash {
        self.transform_identity
    }

    /// Stable transform name.
    #[must_use]
    pub fn transform_name(&self) -> &str {
        &self.transform_name
    }

    /// Complete configuration identity.
    #[must_use]
    pub const fn configuration_identity(&self) -> ContentHash {
        self.configuration_identity
    }

    /// Configuration schema version.
    #[must_use]
    pub const fn configuration_version(&self) -> u32 {
        self.configuration_version
    }

    /// Unit contract.
    #[must_use]
    pub const fn unit_contract(&self) -> CinematicUnitContract {
        self.unit_contract
    }

    /// Clock contract.
    #[must_use]
    pub const fn clock(&self) -> CinematicClock {
        self.clock
    }

    /// Transform/bias disposition.
    #[must_use]
    pub const fn transform_disposition(&self) -> &CinematicTransformDisposition {
        &self.transform_disposition
    }

    /// Canonically sorted disclosure set.
    #[must_use]
    pub fn no_claims(&self) -> &[CinematicNoClaim] {
        &self.no_claims
    }

    /// Optional declared acoustic calibration receipt.
    #[must_use]
    pub const fn acoustic_calibration(&self) -> Option<DeclaredAcousticCalibrationReceipt> {
        self.acoustic_calibration
    }

    /// Deterministic binary encoding used by manifests and record identities.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let label = self.transform_disposition.label().unwrap_or("").as_bytes();
        let mut out = Vec::with_capacity(256 + self.transform_name.len() + label.len());
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&self.schema_version.to_le_bytes());
        out.push(self.artifact_kind.tag());
        out.push(self.authority_class.tag());
        out.push(self.authority_class.sound_tag());
        out.extend_from_slice(self.artifact_identity.as_bytes());
        out.extend_from_slice(self.source_identity.as_bytes());
        out.extend_from_slice(self.transform_identity.as_bytes());
        out.extend_from_slice(self.configuration_identity.as_bytes());
        out.extend_from_slice(&self.configuration_version.to_le_bytes());
        out.push(self.unit_contract.tag());
        out.push(self.clock.domain.tag());
        out.extend_from_slice(&self.clock.ticks_per_second_numerator.to_le_bytes());
        out.extend_from_slice(&self.clock.ticks_per_second_denominator.to_le_bytes());
        out.extend_from_slice(&self.clock.start_tick.to_le_bytes());
        out.extend_from_slice(&self.clock.end_tick_exclusive.to_le_bytes());
        push_string(&mut out, &self.transform_name);
        out.push(self.transform_disposition.tag());
        push_string(&mut out, self.transform_disposition.label().unwrap_or(""));
        match self.acoustic_calibration {
            Some(receipt) => {
                out.push(1);
                out.extend_from_slice(receipt.dataset_identity.as_bytes());
                out.extend_from_slice(receipt.method_identity.as_bytes());
                out.extend_from_slice(receipt.validity_identity.as_bytes());
                out.extend_from_slice(&receipt.version.to_le_bytes());
            }
            None => out.push(0),
        }
        let no_claim_count = u16::try_from(self.no_claims.len()).unwrap_or(u16::MAX);
        out.extend_from_slice(&no_claim_count.to_le_bytes());
        out.extend(self.no_claims.iter().map(|claim| claim.tag()));
        out
    }

    /// Decode and revalidate the exact supported subset.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, CinematicAuthorityError> {
        let mut decoder = Decoder::new(bytes);
        if decoder.take(MAGIC.len())? != MAGIC {
            return Err(CinematicAuthorityError::BadMagic);
        }
        let schema_version = decoder.u16()?;
        if schema_version != CINEMATIC_AUTHORITY_SCHEMA_VERSION {
            return Err(CinematicAuthorityError::UnsupportedSchemaVersion(
                schema_version,
            ));
        }
        let artifact_kind = CinematicArtifactKind::from_tag(decoder.u8()?)?;
        let authority_class = CinematicAuthorityClass::from_tags(decoder.u8()?, decoder.u8()?)?;
        let artifact_identity = decoder.hash()?;
        let source_identity = decoder.hash()?;
        let transform_identity = decoder.hash()?;
        let configuration_identity = decoder.hash()?;
        let configuration_version = decoder.u32()?;
        let unit_contract = CinematicUnitContract::from_tag(decoder.u8()?)?;
        let clock = CinematicClock::try_new(
            CinematicClockDomain::from_tag(decoder.u8()?)?,
            decoder.u32()?,
            decoder.u32()?,
            decoder.i64()?,
            decoder.i64()?,
        )?;
        let transform_name = decoder.string()?;
        let disposition_tag = decoder.u8()?;
        let disposition_label = decoder.string()?;
        let transform_disposition = match disposition_tag {
            1 if disposition_label.is_empty() => CinematicTransformDisposition::ModelState,
            2 if disposition_label.is_empty() => CinematicTransformDisposition::MonteCarloEstimator,
            3 => CinematicTransformDisposition::BiasedVisualization(disposition_label),
            4 => CinematicTransformDisposition::SoundSynthesis(disposition_label),
            1 | 2 => return Err(CinematicAuthorityError::UnexpectedDispositionLabel),
            _ => {
                return Err(CinematicAuthorityError::UnknownDispositionTag(
                    disposition_tag,
                ));
            }
        };
        let acoustic_calibration = match decoder.u8()? {
            0 => None,
            1 => Some(DeclaredAcousticCalibrationReceipt::try_new(
                decoder.hash()?,
                decoder.hash()?,
                decoder.hash()?,
                decoder.u32()?,
            )?),
            tag => return Err(CinematicAuthorityError::UnknownCalibrationPresenceTag(tag)),
        };
        let count = usize::from(decoder.u16()?);
        if count > MAX_CINEMATIC_NO_CLAIMS {
            return Err(CinematicAuthorityError::TooManyNoClaims {
                got: count,
                max: MAX_CINEMATIC_NO_CLAIMS,
            });
        }
        let mut no_claims = Vec::with_capacity(count);
        for _ in 0..count {
            no_claims.push(CinematicNoClaim::from_tag(decoder.u8()?)?);
        }
        if !decoder.is_finished() {
            return Err(CinematicAuthorityError::TrailingBytes);
        }
        Self::try_new(CinematicAuthorityInput {
            schema_version,
            artifact_kind,
            authority_class,
            artifact_identity,
            source_identity,
            transform_identity,
            transform_name,
            configuration_identity,
            configuration_version,
            unit_contract,
            clock,
            transform_disposition,
            no_claims,
            acoustic_calibration,
        })
    }

    /// Domain-separated identity of the complete admitted record.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        hash_domain(IDENTITY_DOMAIN, &self.canonical_bytes())
    }

    /// Deterministic machine manifest. Human disclosures below use the same
    /// no-claim enum, preventing machine/human limitation drift.
    #[must_use]
    pub fn to_manifest_json(&self) -> String {
        let mut out = String::with_capacity(1024);
        let sound = match self.authority_class {
            CinematicAuthorityClass::Sound(tier) => tier.code(),
            _ => "none",
        };
        let label = self.transform_disposition.label().unwrap_or("");
        let _ = write!(
            out,
            "{{\"schema_version\":{},\"artifact_kind\":\"{}\",\"authority_class\":\"{}\",\"sound_authority\":\"{}\",\"artifact_identity\":\"{}\",\"source_identity\":\"{}\",\"transform_identity\":\"{}\",\"transform_name\":\"{}\",\"configuration_identity\":\"{}\",\"configuration_version\":{},\"unit_contract\":\"{}\",\"clock\":{{\"domain\":\"{}\",\"ticks_per_second_numerator\":{},\"ticks_per_second_denominator\":{},\"start_tick\":{},\"end_tick_exclusive\":{}}},\"transform_disposition\":\"{}\",\"transform_label\":\"{}\",\"acoustic_calibration\":",
            self.schema_version,
            self.artifact_kind.code(),
            self.authority_class.code(),
            sound,
            self.artifact_identity,
            self.source_identity,
            self.transform_identity,
            self.transform_name,
            self.configuration_identity,
            self.configuration_version,
            self.unit_contract.code(),
            self.clock.domain.code(),
            self.clock.ticks_per_second_numerator,
            self.clock.ticks_per_second_denominator,
            self.clock.start_tick,
            self.clock.end_tick_exclusive,
            self.transform_disposition.code(),
            label,
        );
        match self.acoustic_calibration {
            Some(receipt) => {
                let _ = write!(
                    out,
                    "{{\"dataset_identity\":\"{}\",\"method_identity\":\"{}\",\"validity_identity\":\"{}\",\"version\":{}}}",
                    receipt.dataset_identity,
                    receipt.method_identity,
                    receipt.validity_identity,
                    receipt.version,
                );
            }
            None => out.push_str("null"),
        }
        out.push_str(",\"no_claims\":[");
        for (index, claim) in self.no_claims.iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"code\":\"{}\",\"statement\":\"{}\"}}",
                claim.code(),
                claim.statement()
            );
        }
        out.push_str("]}");
        out
    }

    /// Human-facing disclosure/end-credit text derived from the same record.
    #[must_use]
    pub fn human_disclosure(&self) -> String {
        let mut out = format!(
            "Authority: {}. Source: {}. Transform: {} ({}).",
            self.authority_class.code(),
            self.source_identity,
            self.transform_name,
            self.transform_identity
        );
        for claim in &self.no_claims {
            out.push('\n');
            out.push_str(claim.statement());
        }
        out
    }
}

/// Closed, stable admission failures for agent and E2E diagnostics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicAuthorityError {
    /// Unknown/incompatible schema.
    UnsupportedSchemaVersion(u16),
    /// Wrong binary framing magic.
    BadMagic,
    /// Input ended before the declared frame completed.
    Truncated,
    /// Bytes remained after a complete record.
    TrailingBytes,
    /// Unknown payload tag.
    UnknownArtifactKindTag(u8),
    /// Unknown authority tag.
    UnknownAuthorityClassTag(u8),
    /// Unknown sound tier.
    UnknownSoundAuthorityTag(u8),
    /// Non-sound class carried a sound tag.
    UnexpectedSoundAuthorityTag(u8),
    /// Unknown units tag.
    UnknownUnitContractTag(u8),
    /// Unknown clock-domain tag.
    UnknownClockDomainTag(u8),
    /// Unknown disposition tag.
    UnknownDispositionTag(u8),
    /// Raw disposition unexpectedly carried a label.
    UnexpectedDispositionLabel,
    /// Unknown calibration-option tag.
    UnknownCalibrationPresenceTag(u8),
    /// Unknown no-claim declaration.
    UnknownNoClaimTag(u8),
    /// A required content identity was all zero.
    MissingIdentity(&'static str),
    /// Configuration versions start at one.
    InvalidConfigurationVersion,
    /// Calibration versions start at one.
    InvalidCalibrationVersion,
    /// Clock numerator and denominator must both be nonzero.
    InvalidClockRate,
    /// Clock range was inverted.
    InvalidClockRange,
    /// Timeless clocks must use the canonical 1/1, 0..0 encoding.
    InvalidTimelessClock,
    /// Transform label was empty, oversized, or outside the identifier grammar.
    InvalidTransformLabel,
    /// Payload kind and authority class disagree.
    IncompatibleArtifactClass,
    /// Raw/biased/synthesized disposition disagrees with authority.
    IncompatibleTransformDisposition,
    /// Calibrated sound omitted external calibration evidence.
    MissingCalibrationReceipt,
    /// A non-calibrated record attempted to carry calibration evidence.
    UnexpectedCalibrationReceipt,
    /// Required disclosure was absent.
    MissingNoClaim(CinematicNoClaim),
    /// Duplicate disclosures refuse rather than being silently normalized.
    DuplicateNoClaim(CinematicNoClaim),
    /// No-claim count exceeded the bounded schema.
    TooManyNoClaims {
        /// Number supplied by the caller or decoder.
        got: usize,
        /// Maximum admitted by this schema.
        max: usize,
    },
    /// A derived record did not name its parent artifact as source.
    SourceDoesNotMatchParent,
    /// Derivation attempted an illegal authority promotion or domain reversal.
    IllegalAuthorityTransition {
        /// Parent class.
        from: CinematicAuthorityClass,
        /// Requested child class.
        to: CinematicAuthorityClass,
    },
}

impl CinematicAuthorityError {
    /// Stable reason code used by structured diagnostics and negative twins.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::UnsupportedSchemaVersion(_) => "cinematic-authority-unsupported-schema",
            Self::BadMagic => "cinematic-authority-bad-magic",
            Self::Truncated => "cinematic-authority-truncated",
            Self::TrailingBytes => "cinematic-authority-trailing-bytes",
            Self::UnknownArtifactKindTag(_) => "cinematic-authority-unknown-artifact-kind",
            Self::UnknownAuthorityClassTag(_) => "cinematic-authority-unknown-class",
            Self::UnknownSoundAuthorityTag(_) => "cinematic-authority-unknown-sound-class",
            Self::UnexpectedSoundAuthorityTag(_) => "cinematic-authority-unexpected-sound-class",
            Self::UnknownUnitContractTag(_) => "cinematic-authority-unknown-unit-contract",
            Self::UnknownClockDomainTag(_) => "cinematic-authority-unknown-clock-domain",
            Self::UnknownDispositionTag(_) => "cinematic-authority-unknown-disposition",
            Self::UnexpectedDispositionLabel => "cinematic-authority-unexpected-disposition-label",
            Self::UnknownCalibrationPresenceTag(_) => {
                "cinematic-authority-unknown-calibration-presence"
            }
            Self::UnknownNoClaimTag(_) => "cinematic-authority-unknown-no-claim",
            Self::MissingIdentity(_) => "cinematic-authority-missing-identity",
            Self::InvalidConfigurationVersion => "cinematic-authority-invalid-config-version",
            Self::InvalidCalibrationVersion => "cinematic-authority-invalid-calibration-version",
            Self::InvalidClockRate => "cinematic-authority-invalid-clock-rate",
            Self::InvalidClockRange => "cinematic-authority-invalid-clock-range",
            Self::InvalidTimelessClock => "cinematic-authority-invalid-timeless-clock",
            Self::InvalidTransformLabel => "cinematic-authority-invalid-transform-label",
            Self::IncompatibleArtifactClass => "cinematic-authority-incompatible-artifact-class",
            Self::IncompatibleTransformDisposition => {
                "cinematic-authority-incompatible-transform-disposition"
            }
            Self::MissingCalibrationReceipt => "cinematic-authority-missing-calibration",
            Self::UnexpectedCalibrationReceipt => "cinematic-authority-unexpected-calibration",
            Self::MissingNoClaim(_) => "cinematic-authority-missing-no-claim",
            Self::DuplicateNoClaim(_) => "cinematic-authority-duplicate-no-claim",
            Self::TooManyNoClaims { .. } => "cinematic-authority-too-many-no-claims",
            Self::SourceDoesNotMatchParent => "cinematic-authority-source-parent-mismatch",
            Self::IllegalAuthorityTransition { .. } => "cinematic-authority-illegal-promotion",
        }
    }
}

impl fmt::Display for CinematicAuthorityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {self:?}", self.code())
    }
}

impl std::error::Error for CinematicAuthorityError {}

/// Required disclosures for a class. Callers may add stricter declarations.
#[must_use]
pub fn required_no_claims(class: CinematicAuthorityClass) -> Vec<CinematicNoClaim> {
    let mut claims = vec![
        CinematicNoClaim::PhysicalStoppingTime,
        CinematicNoClaim::RingConeRanking,
        CinematicNoClaim::ContactLawValidation,
        CinematicNoClaim::TerminalSingularityPrediction,
        CinematicNoClaim::MouldVideoCorrespondence,
    ];
    match class {
        CinematicAuthorityClass::VisualizationDerivative => {
            claims.push(CinematicNoClaim::RawEvidenceFromBiasedArtifact);
        }
        CinematicAuthorityClass::Sound(
            SoundAuthority::Artistic | SoundAuthority::PhysicallyInformed,
        ) => {
            claims.push(CinematicNoClaim::CalibratedAcousticPrediction);
        }
        CinematicAuthorityClass::SimulatedState
        | CinematicAuthorityClass::MonteCarloRender
        | CinematicAuthorityClass::Sound(SoundAuthority::Calibrated) => {}
    }
    claims.sort_unstable();
    claims
}

fn transition_allowed(from: CinematicAuthorityClass, to: CinematicAuthorityClass) -> bool {
    if from == to {
        return true;
    }
    matches!(
        (from, to),
        (
            CinematicAuthorityClass::SimulatedState,
            CinematicAuthorityClass::MonteCarloRender
                | CinematicAuthorityClass::VisualizationDerivative
                | CinematicAuthorityClass::Sound(
                    SoundAuthority::Artistic | SoundAuthority::PhysicallyInformed
                )
        ) | (
            CinematicAuthorityClass::MonteCarloRender,
            CinematicAuthorityClass::VisualizationDerivative
        ) | (
            CinematicAuthorityClass::Sound(SoundAuthority::PhysicallyInformed),
            CinematicAuthorityClass::Sound(SoundAuthority::Artistic | SoundAuthority::Calibrated)
        ) | (
            CinematicAuthorityClass::Sound(SoundAuthority::Calibrated),
            CinematicAuthorityClass::Sound(
                SoundAuthority::Artistic | SoundAuthority::PhysicallyInformed
            )
        )
    )
}

fn validate_artifact_class(
    artifact: CinematicArtifactKind,
    class: CinematicAuthorityClass,
) -> Result<(), CinematicAuthorityError> {
    let compatible = artifact == CinematicArtifactKind::ManifestOrLog
        || matches!(
            (artifact, class),
            (
                CinematicArtifactKind::SimulationState,
                CinematicAuthorityClass::SimulatedState
            ) | (
                CinematicArtifactKind::RenderEstimate,
                CinematicAuthorityClass::MonteCarloRender
            ) | (
                CinematicArtifactKind::Visualization | CinematicArtifactKind::MediaDerivative,
                CinematicAuthorityClass::VisualizationDerivative
            ) | (
                CinematicArtifactKind::Audio,
                CinematicAuthorityClass::Sound(_)
            )
        );
    if compatible {
        Ok(())
    } else {
        Err(CinematicAuthorityError::IncompatibleArtifactClass)
    }
}

fn validate_disposition(
    disposition: &CinematicTransformDisposition,
    class: CinematicAuthorityClass,
) -> Result<(), CinematicAuthorityError> {
    let compatible = matches!(
        (disposition, class),
        (
            CinematicTransformDisposition::ModelState,
            CinematicAuthorityClass::SimulatedState
        ) | (
            CinematicTransformDisposition::MonteCarloEstimator,
            CinematicAuthorityClass::MonteCarloRender
        ) | (
            CinematicTransformDisposition::BiasedVisualization(_),
            CinematicAuthorityClass::VisualizationDerivative
        ) | (
            CinematicTransformDisposition::SoundSynthesis(_),
            CinematicAuthorityClass::Sound(_)
        )
    );
    if !compatible {
        return Err(CinematicAuthorityError::IncompatibleTransformDisposition);
    }
    if let Some(label) = disposition.label() {
        validate_label(label)?;
    }
    Ok(())
}

fn validate_label(label: &str) -> Result<(), CinematicAuthorityError> {
    let bytes = label.as_bytes();
    if bytes.is_empty()
        || bytes.len() > MAX_CINEMATIC_LABEL_BYTES
        || !bytes[0].is_ascii_lowercase()
        || !bytes.iter().all(|byte| {
            byte.is_ascii_lowercase()
                || byte.is_ascii_digit()
                || matches!(byte, b'-' | b'_' | b'.' | b'/')
        })
    {
        return Err(CinematicAuthorityError::InvalidTransformLabel);
    }
    Ok(())
}

fn check_identity(
    field: &'static str,
    identity: ContentHash,
) -> Result<(), CinematicAuthorityError> {
    if identity.as_bytes().iter().all(|byte| *byte == 0) {
        Err(CinematicAuthorityError::MissingIdentity(field))
    } else {
        Ok(())
    }
}

fn push_string(out: &mut Vec<u8>, value: &str) {
    let len = u16::try_from(value.len()).expect("validated cinematic label fits u16");
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(value.as_bytes());
}

struct Decoder<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Decoder<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], CinematicAuthorityError> {
        let end = self
            .offset
            .checked_add(count)
            .ok_or(CinematicAuthorityError::Truncated)?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or(CinematicAuthorityError::Truncated)?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, CinematicAuthorityError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, CinematicAuthorityError> {
        Ok(u16::from_le_bytes(
            self.take(2)?
                .try_into()
                .map_err(|_| CinematicAuthorityError::Truncated)?,
        ))
    }

    fn u32(&mut self) -> Result<u32, CinematicAuthorityError> {
        Ok(u32::from_le_bytes(
            self.take(4)?
                .try_into()
                .map_err(|_| CinematicAuthorityError::Truncated)?,
        ))
    }

    fn i64(&mut self) -> Result<i64, CinematicAuthorityError> {
        Ok(i64::from_le_bytes(
            self.take(8)?
                .try_into()
                .map_err(|_| CinematicAuthorityError::Truncated)?,
        ))
    }

    fn hash(&mut self) -> Result<ContentHash, CinematicAuthorityError> {
        ContentHash::from_slice(self.take(32)?).ok_or(CinematicAuthorityError::Truncated)
    }

    fn string(&mut self) -> Result<String, CinematicAuthorityError> {
        let len = usize::from(self.u16()?);
        if len > MAX_CINEMATIC_LABEL_BYTES {
            return Err(CinematicAuthorityError::InvalidTransformLabel);
        }
        let bytes = self.take(len)?;
        let value = core::str::from_utf8(bytes)
            .map_err(|_| CinematicAuthorityError::InvalidTransformLabel)?;
        Ok(value.to_owned())
    }

    const fn is_finished(&self) -> bool {
        self.offset == self.bytes.len()
    }
}
