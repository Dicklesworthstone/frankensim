//! Deterministic offline stereo spatialization for Euler-disc modal stems.
//!
//! Renders one or more point-like mono sources relative to a listener. It
//! applies sampled propagation delay, inverse-distance attenuation, a frozen
//! equal-power pan law, listener-microphone directivity, and an optional bounded
//! stereo room impulse response.  Inputs may be generic mono samples or one
//! component of [`ModalStemFrame`].  A separate dry-bypass transaction copies
//! already-stereo frames bit-for-bit.
//!
//! This is a parameterized production sound transform, not calibrated
//! acoustics.  It does not claim BEM, HRTF, head shadowing, diffraction,
//! occlusion, radiated pressure, absolute SPL, or a measured room response.
//! The opt-in retarded policies reconstruct the Doppler source clock for tracks
//! subsonic relative to the declared medium, including uniform ambient flow;
//! the two deposition policies do not.
//! `PhysicallyParameterized` identifies parameter
//! provenance; it does not promote any of those no-claims.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_math::det;

use crate::{audio_artifact::StereoSample, modal_synthesis::ModalStemFrame};

mod retarded;

/// Version of the geometry, delay, pan, directivity, and convolution semantics.
pub const SPATIAL_AUDIO_ALGORITHM_VERSION: u32 = 3;
/// Maximum frames between execution-scope cancellation polls.
pub const SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES: usize = 256;
/// Hard source-count ceiling above caller-controlled budgets.
pub const MAX_SPATIAL_AUDIO_SOURCES: usize = 32;
/// Hard sample-rate ceiling [Hz].
pub const MAX_SPATIAL_AUDIO_SAMPLE_RATE_HZ: u32 = 384_000;
/// Hard room-response length ceiling per channel [samples].
pub const MAX_SPATIAL_AUDIO_ROOM_IR_TAPS: usize = 65_536;

const CONFIG_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.spatial-audio-config.v1";
const INPUT_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.spatial-audio-input.v1";
const OUTPUT_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.spatial-audio-output.v1";
const ROOM_IR_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.stereo-room-ir.v1";
const DRY_BYPASS_CONFIG_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-cinematic.spatial-audio-dry-bypass-config.v1";
const POSE_UNIT_TOLERANCE: f64 = 1.0e-9;

/// Provenance tier of spatial parameters and source material.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpatialAudioAuthority {
    /// Values came from declared physical dimensions or measurements.
    ///
    /// This is parameter provenance only, not a calibrated-acoustics claim.
    PhysicallyParameterized,
    /// At least one value was selected for presentation rather than measured.
    Artistic,
}

impl SpatialAudioAuthority {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::PhysicallyParameterized => "physically-parameterized",
            Self::Artistic => "artistic",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::PhysicallyParameterized => 1,
            Self::Artistic => 2,
        }
    }

    const fn combine(self, other: Self) -> Self {
        if matches!(self, Self::Artistic) || matches!(other, Self::Artistic) {
            Self::Artistic
        } else {
            Self::PhysicallyParameterized
        }
    }
}

/// Sampled propagation-delay policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpatialDelayPolicy {
    /// Place each input sample at `input_index + ceil(distance * fs / c)`.
    ///
    /// This never advances an arrival relative to the continuous delay, at the
    /// cost of at most one sample of additional delay.
    IntegerCeiling,
    /// Split each input sample between the floor and ceiling arrival samples.
    ///
    /// Both deposits are at or after the emission sample.  This is a causal
    /// sampled linear interpolator, not a band-limited fractional-delay filter.
    LinearFloorCeil,
    /// Sample the dry source at its causal emission time on the receiver clock.
    ///
    /// Source and receiver positions are piecewise linear on the physical sample
    /// clock and must remain strictly subsonic in the declared resting medium.
    /// Listener axes are read at reception. A normalized 32-tap Lanczos kernel
    /// interpolates within the declared emission interval; samples outside that
    /// interval and incomplete endpoint stencils use explicit zeros. No
    /// sample-density or Doppler amplitude multiplier is introduced.
    /// This is a finite-kernel approximation, not an ideal antialias filter or
    /// the moving-boundary acoustic pressure law.
    RetardedSinc {
        /// Declared dry-source bandwidth [Hz], at most one quarter of the sample
        /// rate. This sampling margin keeps the finite kernel in its tested band.
        /// The producer must establish this band; it is not inferred from samples.
        source_bandwidth_hz: f64,
        /// Declared received bandwidth [Hz], strictly below output Nyquist.
        /// Must cover admitted Doppler expansion and geometric modulation.
        receiver_bandwidth_hz: f64,
    },
    /// Receiver-clock reconstruction in a homogeneous, uniformly moving medium.
    ///
    /// Source and listener tracks remain laboratory positions. The causal ray
    /// satisfies `|x_receiver(t) - x_source(tau) - U * (t-tau)| = c * (t-tau)`;
    /// both endpoint velocities must be subsonic relative to `U`. Panning,
    /// directivity, and the declared inverse-distance gain use this medium-frame
    /// ray, not the laboratory displacement. The same finite-kernel and source
    /// bandwidth requirements as [`Self::RetardedSinc`] apply. This adds no
    /// convective pressure coefficient, shear/refraction, or moving-room model.
    RetardedSincUniformFlow {
        /// Declared dry-source bandwidth [Hz], at most one quarter sample rate.
        source_bandwidth_hz: f64,
        /// Declared received bandwidth [Hz], strictly below output Nyquist.
        /// Must cover admitted Doppler expansion and geometric modulation.
        receiver_bandwidth_hz: f64,
        /// Constant medium velocity in the laboratory coordinate system [m/s].
        medium_velocity_m_per_s: [f64; 3],
    },
}

impl SpatialDelayPolicy {
    const fn tag(self) -> u8 {
        match self {
            Self::IntegerCeiling => 1,
            Self::LinearFloorCeil => 2,
            Self::RetardedSinc { .. } => 3,
            Self::RetardedSincUniformFlow { .. } => 4,
        }
    }

    const fn is_retarded(self) -> bool {
        self.retarded_bands().is_some()
    }

    const fn retarded_bands(self) -> Option<(f64, f64)> {
        match self {
            Self::RetardedSinc {
                source_bandwidth_hz,
                receiver_bandwidth_hz,
            }
            | Self::RetardedSincUniformFlow {
                source_bandwidth_hz,
                receiver_bandwidth_hz,
                ..
            } => Some((source_bandwidth_hz, receiver_bandwidth_hz)),
            Self::IntegerCeiling | Self::LinearFloorCeil => None,
        }
    }

    const fn medium_velocity_m_per_s(self) -> [f64; 3] {
        match self {
            Self::RetardedSincUniformFlow {
                medium_velocity_m_per_s,
                ..
            } => medium_velocity_m_per_s,
            _ => [0.0; 3],
        }
    }
}

/// Explicit publication horizon for propagation and room-response tails.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpatialOutputHorizon {
    /// Publish the complete propagation and room-response tail.
    PreserveTail,
    /// Publish exactly the common input-frame horizon.
    ///
    /// Later deposits are deterministically discarded. Diagnostics retain the
    /// natural horizon and discarded count, and this policy is identity-bound.
    ClampToInputFrames,
}

impl SpatialOutputHorizon {
    const fn tag(self) -> u8 {
        match self {
            Self::PreserveTail => 1,
            Self::ClampToInputFrames => 2,
        }
    }
}

/// Listener-microphone polar response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MicrophoneDirectivity {
    /// Unit gain in every direction.
    Omnidirectional,
    /// First-order cardioid with an explicit rear-axis amplitude floor.
    Cardioid {
        /// Amplitude gain on the exact rear axis, in `[0, 1]`.
        rear_floor_gain: f64,
    },
}

impl MicrophoneDirectivity {
    fn hash_into(self, hasher: &mut DomainHasher) {
        match self {
            Self::Omnidirectional => hasher.update(&[1]),
            Self::Cardioid { rear_floor_gain } => {
                hasher.update(&[2]);
                hash_f64(hasher, rear_floor_gain);
            }
        }
    }

    fn gain(self, forward_cosine: f64) -> f64 {
        match self {
            Self::Omnidirectional => 1.0,
            Self::Cardioid { rear_floor_gain } => {
                rear_floor_gain + (1.0 - rear_floor_gain) * 0.5 * (1.0 + forward_cosine)
            }
        }
    }
}

/// One dry modal component selected as a mono source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpatialStemComponent {
    /// Disc-body stem.
    Disc,
    /// Glass-plate stem.
    GlassPlate,
    /// Base-assembly stem.
    BaseAssembly,
}

impl SpatialStemComponent {
    const fn tag(self) -> u8 {
        match self {
            Self::Disc => 1,
            Self::GlassPlate => 2,
            Self::BaseAssembly => 3,
        }
    }

    fn sample(self, frame: ModalStemFrame) -> f64 {
        match self {
            Self::Disc => frame.disc_fs,
            Self::GlassPlate => frame.glass_plate_fs,
            Self::BaseAssembly => frame.base_assembly_fs,
        }
    }
}

/// Borrowed mono source signal.
#[derive(Debug, Clone, Copy)]
pub enum SpatialMonoSignal<'a> {
    /// Generic mono samples in digital-full-scale coordinates.
    Samples(&'a [f64]),
    /// One component selected without copying from modal synthesis output.
    ModalStemFrames {
        /// Dry component frames.
        frames: &'a [ModalStemFrame],
        /// Component to read from every frame.
        component: SpatialStemComponent,
    },
}

impl SpatialMonoSignal<'_> {
    fn len(self) -> usize {
        match self {
            Self::Samples(samples) => samples.len(),
            Self::ModalStemFrames { frames, .. } => frames.len(),
        }
    }

    fn sample(self, index: usize) -> f64 {
        match self {
            Self::Samples(samples) => samples[index],
            Self::ModalStemFrames { frames, component } => component.sample(frames[index]),
        }
    }

    fn hash_header(self, hasher: &mut DomainHasher) {
        match self {
            Self::Samples(_) => hasher.update(&[1]),
            Self::ModalStemFrames { component, .. } => {
                hasher.update(&[2, component.tag()]);
            }
        }
    }
}

/// Static or sample-synchronous point-source positions [m].
#[derive(Debug, Clone, Copy)]
pub enum SourcePositionTrack<'a> {
    /// One position used for the whole transaction.
    Static([f64; 3]),
    /// One position at every source-emission frame.
    PerFrame(&'a [[f64; 3]]),
}

impl SourcePositionTrack<'_> {
    fn position(self, index: usize) -> [f64; 3] {
        match self {
            Self::Static(position) => position,
            Self::PerFrame(positions) => positions[index],
        }
    }

    fn validate_len(self, expected: usize) -> bool {
        match self {
            Self::Static(_) => true,
            Self::PerFrame(positions) => positions.len() == expected,
        }
    }
}

/// Listener position and orthonormal forward/right axes in world coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ListenerPose {
    /// Listener/microphone position [m].
    pub position_m: [f64; 3],
    /// Unit vector toward the listener's front.
    pub forward_unit: [f64; 3],
    /// Unit vector toward the listener's right channel.
    pub right_unit: [f64; 3],
}

/// Static or sample-synchronous listener poses.
#[derive(Debug, Clone, Copy)]
pub enum ListenerPoseTrack<'a> {
    /// One listener pose used for the whole transaction.
    Static(ListenerPose),
    /// One listener pose at every physical sample time.
    ///
    /// Deposition policies read emission-time poses and require the source
    /// frame count. Retarded policies read reception-time poses; include enough
    /// poses to receive the final source sample, even when publication is clamped.
    PerFrame(&'a [ListenerPose]),
}

impl ListenerPoseTrack<'_> {
    fn pose(self, index: usize) -> ListenerPose {
        match self {
            Self::Static(pose) => pose,
            Self::PerFrame(poses) => poses[index],
        }
    }

    fn validate_len(self, expected: usize) -> bool {
        match self {
            Self::Static(_) => true,
            Self::PerFrame(poses) => poses.len() == expected,
        }
    }
}

/// One identified mono source and its position track.
#[derive(Debug, Clone, Copy)]
pub struct SpatialAudioSource<'a> {
    /// Nonzero upstream signal identity.
    pub source_identity: ContentHash,
    /// Mono samples or modal stem frames.
    pub signal: SpatialMonoSignal<'a>,
    /// World-space source position at each emission frame.
    pub positions: SourcePositionTrack<'a>,
    /// Explicit finite linear amplitude gain applied before spatialization.
    ///
    /// This is identity-bound and is never inferred or normalized.
    pub gain_linear: f64,
    /// Authority of this source's samples and positions.
    pub authority: SpatialAudioAuthority,
}

/// Explicit resource and amplitude ceilings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialAudioBudget {
    /// Maximum admitted point sources.
    pub maximum_sources: usize,
    /// Maximum sum of mono input frames across all sources.
    pub maximum_total_input_frames: u64,
    /// Maximum stereo frames after delay and room-response tail.
    pub maximum_output_frames: u64,
    /// Maximum room-response taps per channel.
    pub maximum_room_ir_taps: usize,
    /// Maximum deterministic scalar work estimate.
    pub maximum_work_units: u64,
    /// Maximum aggregate bytes for renderer-owned sample buffers.
    pub maximum_owned_sample_bytes: u64,
    /// Maximum allowed absolute output sample; exceeding it refuses, never clips.
    pub maximum_abs_output_fs: f64,
}

impl SpatialAudioBudget {
    /// Practical 12-second 48-kHz preview budget with a short room response.
    pub const PREVIEW: Self = Self {
        maximum_sources: 8,
        maximum_total_input_frames: 8 * 576_000,
        maximum_output_frames: 700_000,
        maximum_room_ir_taps: 8_192,
        maximum_work_units: 200_000_000,
        maximum_owned_sample_bytes: 96 * 1024 * 1024,
        maximum_abs_output_fs: 1.0,
    };
}

/// Immutable point-source spatialization configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialAudioConfig {
    /// Explicit output and input sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Declared propagation speed relative to the medium [m/s].
    pub speed_of_sound_m_per_s: f64,
    /// Clamp distance for the inverse-distance amplitude law [m].
    ///
    /// Gain is `minimum_distance / max(actual_distance, minimum_distance)`.
    pub minimum_distance_m: f64,
    /// Sampled propagation-delay policy.
    pub delay_policy: SpatialDelayPolicy,
    /// Whether propagation and room-response tails are retained or explicitly
    /// truncated to the common source horizon.
    pub output_horizon: SpatialOutputHorizon,
    /// Listener-microphone polar response.
    pub microphone_directivity: MicrophoneDirectivity,
    /// Authority of the propagation, attenuation, pan, and directivity values.
    pub authority: SpatialAudioAuthority,
    /// Explicit resource and output-amplitude ceilings.
    pub budget: SpatialAudioBudget,
}

/// Content-identified bounded stereo room impulse response.
#[derive(Debug, Clone, PartialEq)]
pub struct StereoRoomImpulseResponse {
    identity: ContentHash,
    sample_rate_hz: u32,
    left_taps: Vec<f64>,
    right_taps: Vec<f64>,
    authority: SpatialAudioAuthority,
}

impl StereoRoomImpulseResponse {
    /// Validate, own, and identify a same-length stereo response.
    ///
    /// Convolution is channel-wise and pure wet: tap zero must contain any
    /// desired direct path.  No normalization or hidden dry signal is added.
    pub fn try_new(
        sample_rate_hz: u32,
        left_taps: Vec<f64>,
        right_taps: Vec<f64>,
        authority: SpatialAudioAuthority,
        cx: &Cx<'_>,
    ) -> Result<Self, SpatialAudioError> {
        validate_sample_rate(sample_rate_hz)?;
        if left_taps.is_empty() {
            return Err(invalid("nonempty room impulse response"));
        }
        if left_taps.len() != right_taps.len() {
            return Err(SpatialAudioError::LengthMismatch {
                field: "room impulse response channels",
                expected: left_taps.len(),
                actual: right_taps.len(),
            });
        }
        if left_taps.len() > MAX_SPATIAL_AUDIO_ROOM_IR_TAPS {
            return Err(limit(
                "room impulse response taps",
                left_taps.len() as u64,
                MAX_SPATIAL_AUDIO_ROOM_IR_TAPS as u64,
            ));
        }
        let mut hasher = DomainHasher::new(ROOM_IR_IDENTITY_DOMAIN);
        hasher.update(&SPATIAL_AUDIO_ALGORITHM_VERSION.to_le_bytes());
        hasher.update(&sample_rate_hz.to_le_bytes());
        hasher.update(&[authority.tag()]);
        hasher.update(&(left_taps.len() as u64).to_le_bytes());
        for (index, (&left, &right)) in left_taps.iter().zip(&right_taps).enumerate() {
            if index % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint(cx)?;
            }
            validate_finite(left, "left room impulse response", index)?;
            validate_finite(right, "right room impulse response", index)?;
            hash_f64(&mut hasher, left);
            hash_f64(&mut hasher, right);
        }
        checkpoint(cx)?;
        let identity = hasher.finalize();
        if is_zero_hash(identity) {
            return Err(SpatialAudioError::InvalidIdentity("room impulse response"));
        }
        Ok(Self {
            identity,
            sample_rate_hz,
            left_taps,
            right_taps,
            authority,
        })
    }

    /// Content identity of the rate, authority, and exact stereo taps.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact sample rate [Hz].
    #[must_use]
    pub const fn sample_rate_hz(&self) -> u32 {
        self.sample_rate_hz
    }

    /// Number of taps in each channel.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.left_taps.len()
    }

    /// Whether the response contains no taps.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.left_taps.is_empty()
    }

    /// Parameter-provenance tier of this response.
    #[must_use]
    pub const fn authority(&self) -> SpatialAudioAuthority {
        self.authority
    }
}

/// Borrowed inputs to one atomic offline spatialization transaction.
#[derive(Debug, Clone, Copy)]
pub struct SpatialAudioRenderInput<'a> {
    /// Ordered point sources. Source order is part of input identity and sum order.
    pub sources: &'a [SpatialAudioSource<'a>],
    /// Listener poses on the physical sample clock; see [`ListenerPoseTrack`].
    pub listener: ListenerPoseTrack<'a>,
    /// Optional exact stereo room response.
    pub room_ir: Option<&'a StereoRoomImpulseResponse>,
}

/// Diagnostics over one complete spatialization result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialAudioDiagnostics {
    /// Common number of emission frames per source.
    pub input_frames_per_source: u64,
    /// Direct-path frames including propagation-delay tail.
    pub direct_output_frames: u64,
    /// Final frames including the optional room-response tail.
    pub final_output_frames: u64,
    /// Final frame count before the explicit output-horizon policy.
    pub natural_final_output_frames: u64,
    /// Natural tail frames intentionally omitted from publication.
    pub discarded_tail_frames: u64,
    /// Largest propagation distance encountered [m]. In uniform flow this is
    /// the acoustic path length in the medium frame, not laboratory separation.
    pub maximum_distance_m: f64,
    /// Largest sampled propagation delay encountered [frames].
    pub maximum_delay_frames: f64,
    /// Evaluated source/receiver pairs whose attenuation used the minimum clamp.
    pub minimum_distance_clamp_count: u64,
    /// Largest absolute final sample.
    pub sample_peak_fs: f64,
}

/// Atomically published stereo samples and binding identities.
#[derive(Debug, Clone, PartialEq)]
pub struct SpatialAudioOutput {
    config_identity: ContentHash,
    input_identity: ContentHash,
    output_identity: ContentHash,
    sample_rate_hz: u32,
    authority: SpatialAudioAuthority,
    room_ir_identity: Option<ContentHash>,
    samples: Vec<StereoSample>,
    diagnostics: SpatialAudioDiagnostics,
}

impl SpatialAudioOutput {
    /// Complete point-render or dry-bypass configuration identity.
    #[must_use]
    pub const fn config_identity(&self) -> ContentHash {
        self.config_identity
    }

    /// Identity of the exact signals, tracks, listener poses, and room response.
    #[must_use]
    pub const fn input_identity(&self) -> ContentHash {
        self.input_identity
    }

    /// Identity of configuration, inputs, and exact output sample bits.
    #[must_use]
    pub const fn output_identity(&self) -> ContentHash {
        self.output_identity
    }

    /// Exact stereo sample rate [Hz].
    #[must_use]
    pub const fn sample_rate_hz(&self) -> u32 {
        self.sample_rate_hz
    }

    /// Combined parameter/source authority; any artistic input makes it artistic.
    #[must_use]
    pub const fn authority(&self) -> SpatialAudioAuthority {
        self.authority
    }

    /// Optional room-response content identity.
    #[must_use]
    pub const fn room_ir_identity(&self) -> Option<ContentHash> {
        self.room_ir_identity
    }

    /// Complete final stereo samples.
    #[must_use]
    pub fn samples(&self) -> &[StereoSample] {
        &self.samples
    }

    /// Complete transaction diagnostics.
    #[must_use]
    pub const fn diagnostics(&self) -> SpatialAudioDiagnostics {
        self.diagnostics
    }
}

/// Typed refusal from configuration, spatialization, convolution, or bypass.
#[derive(Debug, Clone, PartialEq)]
pub enum SpatialAudioError {
    /// A scalar or enum combination violates the v1 contract.
    InvalidConfig(&'static str),
    /// An upstream identity is the all-zero sentinel.
    InvalidIdentity(&'static str),
    /// A signal, position, or pose track has the wrong length.
    LengthMismatch {
        /// Named input.
        field: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        actual: usize,
    },
    /// A scalar input or derived value is NaN or infinite.
    NonFinite {
        /// Named value.
        field: &'static str,
        /// Sample or tap index.
        index: usize,
    },
    /// Listener axes are not unit length and mutually orthogonal.
    InvalidListenerPose {
        /// Emission-frame index.
        frame: usize,
        /// Failed invariant.
        reason: &'static str,
    },
    /// A caller-controlled or hard resource ceiling was exceeded.
    ResourceLimit {
        /// Named resource.
        resource: &'static str,
        /// Requested amount.
        requested: u64,
        /// Admitted maximum.
        limit: u64,
    },
    /// A renderer-owned buffer could not be reserved.
    AllocationFailed(&'static str),
    /// A final sample exceeds the explicit amplitude ceiling.
    OutputAmplitudeExceeded {
        /// Output frame.
        frame: usize,
        /// `left` or `right`.
        channel: &'static str,
        /// Absolute sample magnitude.
        magnitude_fs: f64,
        /// Configured ceiling.
        limit_fs: f64,
    },
    /// The execution scope requested cancellation before atomic publication.
    Cancelled,
}

impl fmt::Display for SpatialAudioError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig(field) => write!(formatter, "invalid spatial-audio {field}"),
            Self::InvalidIdentity(field) => write!(formatter, "invalid {field} identity"),
            Self::LengthMismatch {
                field,
                expected,
                actual,
            } => write!(
                formatter,
                "{field} length mismatch: expected {expected}, got {actual}"
            ),
            Self::NonFinite { field, index } => {
                write!(formatter, "non-finite {field} at index {index}")
            }
            Self::InvalidListenerPose { frame, reason } => {
                write!(
                    formatter,
                    "invalid listener pose at frame {frame}: {reason}"
                )
            }
            Self::ResourceLimit {
                resource,
                requested,
                limit,
            } => write!(
                formatter,
                "spatial-audio {resource} request {requested} exceeds limit {limit}"
            ),
            Self::AllocationFailed(resource) => {
                write!(formatter, "could not reserve spatial-audio {resource}")
            }
            Self::OutputAmplitudeExceeded {
                frame,
                channel,
                magnitude_fs,
                limit_fs,
            } => write!(
                formatter,
                "spatial-audio {channel} sample {frame} magnitude {magnitude_fs} exceeds {limit_fs} FS"
            ),
            Self::Cancelled => formatter.write_str("spatial-audio operation cancelled"),
        }
    }
}

impl std::error::Error for SpatialAudioError {}

/// Validated immutable point-source spatializer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OfflineSpatializer {
    config: SpatialAudioConfig,
    identity: ContentHash,
}

impl OfflineSpatializer {
    /// Validate and content-identify one immutable configuration.
    pub fn try_new(config: SpatialAudioConfig, cx: &Cx<'_>) -> Result<Self, SpatialAudioError> {
        checkpoint(cx)?;
        validate_config(config)?;
        let identity = config_identity(config);
        if is_zero_hash(identity) {
            return Err(SpatialAudioError::InvalidIdentity(
                "spatial-audio configuration",
            ));
        }
        checkpoint(cx)?;
        Ok(Self { config, identity })
    }

    /// Complete configuration identity, including every budget ceiling.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }

    /// Exact admitted configuration.
    #[must_use]
    pub const fn config(self) -> SpatialAudioConfig {
        self.config
    }

    /// Render all sources, optional room response, and tails atomically.
    pub fn spatialize(
        &self,
        input: SpatialAudioRenderInput<'_>,
        cx: &Cx<'_>,
    ) -> Result<SpatialAudioOutput, SpatialAudioError> {
        spatialize_with_checkpoint(*self, input, &mut || checkpoint(cx))
    }
}

/// Copy already-stereo dry frames exactly while retaining typed identities,
/// validation, cancellation, and the same amplitude/memory ceilings.
///
/// No pan, gain, delay, directivity, room response, normalization, or sample
/// conversion is applied.  Every output `f64` bit pattern therefore equals the
/// corresponding input bit pattern, including signed zero.
pub fn bypass_dry_stereo(
    frames: &[StereoSample],
    source_identity: ContentHash,
    sample_rate_hz: u32,
    authority: SpatialAudioAuthority,
    budget: SpatialAudioBudget,
    cx: &Cx<'_>,
) -> Result<SpatialAudioOutput, SpatialAudioError> {
    bypass_dry_stereo_with_checkpoint(
        frames,
        source_identity,
        sample_rate_hz,
        authority,
        budget,
        &mut || checkpoint(cx),
    )
}

fn spatialize_with_checkpoint(
    spatializer: OfflineSpatializer,
    input: SpatialAudioRenderInput<'_>,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<SpatialAudioOutput, SpatialAudioError> {
    checkpoint_fn()?;
    let mut preflight = preflight(spatializer, input, checkpoint_fn)?;
    let direct_len = usize::try_from(preflight.direct_output_frames)
        .map_err(|_| invalid("addressable direct output frame count"))?;
    let final_len = usize::try_from(preflight.final_output_frames)
        .map_err(|_| invalid("addressable final output frame count"))?;

    let mut direct = Vec::<StereoAccumulator>::new();
    direct
        .try_reserve_exact(direct_len)
        .map_err(|_| SpatialAudioError::AllocationFailed("direct accumulators"))?;
    direct.resize(direct_len, StereoAccumulator::ZERO);

    if spatializer.config.delay_policy.is_retarded() {
        retarded::render(
            spatializer.config,
            input,
            &mut direct,
            &mut preflight,
            checkpoint_fn,
        )?;
    } else {
        for (source_index, source) in input.sources.iter().copied().enumerate() {
            for frame in 0..preflight.frame_count {
                if frame % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
                    checkpoint_fn()?;
                }
                let sample = source.signal.sample(frame);
                let source_position = source.positions.position(frame);
                let listener = input.listener.pose(frame);
                let geometry = geometry(
                    source_position,
                    listener,
                    spatializer.config.minimum_distance_m,
                    spatializer.config.speed_of_sound_m_per_s,
                    spatializer.config.sample_rate_hz,
                    spatializer.config.microphone_directivity,
                    frame,
                )?;
                let left = sample
                    * source.gain_linear
                    * geometry.attenuation_gain
                    * geometry.microphone_gain
                    * det::sqrt((1.0 - geometry.pan) * 0.5);
                let right = sample
                    * source.gain_linear
                    * geometry.attenuation_gain
                    * geometry.microphone_gain
                    * det::sqrt((1.0 + geometry.pan) * 0.5);
                deposit(
                    &mut direct,
                    frame,
                    geometry.delay_frames,
                    left,
                    right,
                    spatializer.config.delay_policy,
                    spatializer.config.output_horizon,
                    source_index,
                )?;
            }
        }
    }
    checkpoint_fn()?;

    let mut dry = Vec::<StereoSample>::new();
    dry.try_reserve_exact(direct_len)
        .map_err(|_| SpatialAudioError::AllocationFailed("direct stereo samples"))?;
    for (index, accumulator) in direct.into_iter().enumerate() {
        if index % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        dry.push(accumulator.finish());
    }

    let mut final_samples = if let Some(room_ir) = input.room_ir {
        convolve_room(&dry, final_len, room_ir, checkpoint_fn)?
    } else {
        dry
    };
    checkpoint_fn()?;
    let sample_peak_fs = validate_final_samples(
        &final_samples,
        spatializer.config.budget.maximum_abs_output_fs,
        checkpoint_fn,
    )?;
    let output_identity = output_identity(
        spatializer.identity,
        preflight.input_identity,
        spatializer.config.sample_rate_hz,
        preflight.authority,
        input.room_ir.map(StereoRoomImpulseResponse::identity),
        &final_samples,
        checkpoint_fn,
    )?;
    checkpoint_fn()?;

    Ok(SpatialAudioOutput {
        config_identity: spatializer.identity,
        input_identity: preflight.input_identity,
        output_identity,
        sample_rate_hz: spatializer.config.sample_rate_hz,
        authority: preflight.authority,
        room_ir_identity: input.room_ir.map(StereoRoomImpulseResponse::identity),
        samples: core::mem::take(&mut final_samples),
        diagnostics: SpatialAudioDiagnostics {
            input_frames_per_source: preflight.frame_count as u64,
            direct_output_frames: preflight.direct_output_frames,
            final_output_frames: preflight.final_output_frames,
            natural_final_output_frames: preflight.natural_final_output_frames,
            discarded_tail_frames: preflight.discarded_tail_frames,
            maximum_distance_m: preflight.maximum_distance_m,
            maximum_delay_frames: preflight.maximum_delay_frames,
            minimum_distance_clamp_count: preflight.minimum_distance_clamp_count,
            sample_peak_fs,
        },
    })
}

fn bypass_dry_stereo_with_checkpoint(
    frames: &[StereoSample],
    source_identity: ContentHash,
    sample_rate_hz: u32,
    authority: SpatialAudioAuthority,
    budget: SpatialAudioBudget,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<SpatialAudioOutput, SpatialAudioError> {
    checkpoint_fn()?;
    validate_sample_rate(sample_rate_hz)?;
    validate_budget(budget)?;
    if is_zero_hash(source_identity) {
        return Err(SpatialAudioError::InvalidIdentity("dry stereo source"));
    }
    if frames.is_empty() {
        return Err(invalid("nonempty dry stereo frames"));
    }
    check_limit(
        "dry input frames",
        frames.len() as u64,
        budget.maximum_total_input_frames,
    )?;
    check_limit(
        "dry output frames",
        frames.len() as u64,
        budget.maximum_output_frames,
    )?;
    let bytes = (frames.len() as u64)
        .checked_mul(core::mem::size_of::<StereoSample>() as u64)
        .ok_or_else(|| invalid("dry output byte count"))?;
    check_limit(
        "dry owned sample bytes",
        bytes,
        budget.maximum_owned_sample_bytes,
    )?;
    let work = (frames.len() as u64)
        .checked_mul(4)
        .ok_or_else(|| invalid("dry bypass work"))?;
    check_limit("dry bypass work", work, budget.maximum_work_units)?;

    let config_identity = dry_bypass_config_identity(sample_rate_hz, budget);
    let mut input_hasher = DomainHasher::new(INPUT_IDENTITY_DOMAIN);
    input_hasher.update(&SPATIAL_AUDIO_ALGORITHM_VERSION.to_le_bytes());
    input_hasher.update(config_identity.as_bytes());
    input_hasher.update(source_identity.as_bytes());
    input_hasher.update(&[authority.tag()]);
    input_hasher.update(&(frames.len() as u64).to_le_bytes());
    let mut samples = Vec::new();
    samples
        .try_reserve_exact(frames.len())
        .map_err(|_| SpatialAudioError::AllocationFailed("dry bypass samples"))?;
    let mut peak = 0.0_f64;
    for (index, frame) in frames.iter().copied().enumerate() {
        if index % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        validate_finite(frame.left_fs, "dry left sample", index)?;
        validate_finite(frame.right_fs, "dry right sample", index)?;
        check_amplitude(frame.left_fs, index, "left", budget.maximum_abs_output_fs)?;
        check_amplitude(frame.right_fs, index, "right", budget.maximum_abs_output_fs)?;
        peak = peak.max(frame.left_fs.abs()).max(frame.right_fs.abs());
        hash_f64(&mut input_hasher, frame.left_fs);
        hash_f64(&mut input_hasher, frame.right_fs);
        samples.push(frame);
    }
    checkpoint_fn()?;
    let input_identity = input_hasher.finalize();
    let output_identity = output_identity(
        config_identity,
        input_identity,
        sample_rate_hz,
        authority,
        None,
        &samples,
        checkpoint_fn,
    )?;
    checkpoint_fn()?;
    Ok(SpatialAudioOutput {
        config_identity,
        input_identity,
        output_identity,
        sample_rate_hz,
        authority,
        room_ir_identity: None,
        diagnostics: SpatialAudioDiagnostics {
            input_frames_per_source: frames.len() as u64,
            direct_output_frames: frames.len() as u64,
            final_output_frames: frames.len() as u64,
            natural_final_output_frames: frames.len() as u64,
            discarded_tail_frames: 0,
            maximum_distance_m: 0.0,
            maximum_delay_frames: 0.0,
            minimum_distance_clamp_count: 0,
            sample_peak_fs: peak,
        },
        samples,
    })
}

#[derive(Debug, Clone, Copy)]
struct Preflight {
    frame_count: usize,
    direct_output_frames: u64,
    final_output_frames: u64,
    natural_final_output_frames: u64,
    discarded_tail_frames: u64,
    maximum_distance_m: f64,
    maximum_delay_frames: f64,
    minimum_distance_clamp_count: u64,
    authority: SpatialAudioAuthority,
    input_identity: ContentHash,
}

fn preflight(
    spatializer: OfflineSpatializer,
    input: SpatialAudioRenderInput<'_>,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<Preflight, SpatialAudioError> {
    if input.sources.is_empty() {
        return Err(invalid("nonempty spatial source set"));
    }
    check_limit(
        "source count",
        input.sources.len() as u64,
        spatializer.config.budget.maximum_sources as u64,
    )?;
    let frame_count = input.sources[0].signal.len();
    if frame_count == 0 {
        return Err(invalid("nonempty spatial source signals"));
    }
    let retarded = spatializer.config.delay_policy.is_retarded();
    let listener_frames = match input.listener {
        ListenerPoseTrack::Static(_) => frame_count,
        ListenerPoseTrack::PerFrame(poses) => poses.len(),
    };
    if (retarded && listener_frames < frame_count)
        || (!retarded && !input.listener.validate_len(frame_count))
    {
        let actual = match input.listener {
            ListenerPoseTrack::Static(_) => frame_count,
            ListenerPoseTrack::PerFrame(poses) => poses.len(),
        };
        return Err(SpatialAudioError::LengthMismatch {
            field: "listener pose track",
            expected: frame_count,
            actual,
        });
    }
    let total_input_frames = (frame_count as u64)
        .checked_mul(input.sources.len() as u64)
        .ok_or_else(|| invalid("total spatial source frames"))?;
    check_limit(
        "total input frames",
        total_input_frames,
        spatializer.config.budget.maximum_total_input_frames,
    )?;

    if let Some(room_ir) = input.room_ir {
        if room_ir.sample_rate_hz != spatializer.config.sample_rate_hz {
            return Err(invalid("room impulse response sample rate"));
        }
        check_limit(
            "room impulse response taps",
            room_ir.len() as u64,
            spatializer.config.budget.maximum_room_ir_taps as u64,
        )?;
    }

    let mut hasher = DomainHasher::new(INPUT_IDENTITY_DOMAIN);
    hasher.update(&SPATIAL_AUDIO_ALGORITHM_VERSION.to_le_bytes());
    hasher.update(spatializer.identity.as_bytes());
    hasher.update(&(input.sources.len() as u64).to_le_bytes());
    hasher.update(&(frame_count as u64).to_le_bytes());
    let mut maximum_distance_m = 0.0_f64;
    let mut maximum_delay_frames = 0.0_f64;
    let mut maximum_arrival_offset = 0_u64;
    let mut minimum_distance_clamp_count = 0_u64;
    let mut authority = spatializer.config.authority;

    for (source_index, source) in input.sources.iter().copied().enumerate() {
        if is_zero_hash(source.source_identity) {
            return Err(SpatialAudioError::InvalidIdentity("spatial source"));
        }
        if source.signal.len() != frame_count {
            return Err(SpatialAudioError::LengthMismatch {
                field: "spatial source signal",
                expected: frame_count,
                actual: source.signal.len(),
            });
        }
        if !source.positions.validate_len(frame_count) {
            let actual = match source.positions {
                SourcePositionTrack::Static(_) => frame_count,
                SourcePositionTrack::PerFrame(positions) => positions.len(),
            };
            return Err(SpatialAudioError::LengthMismatch {
                field: "source position track",
                expected: frame_count,
                actual,
            });
        }
        authority = authority.combine(source.authority);
        if !source.gain_linear.is_finite() || source.gain_linear < 0.0 {
            return Err(invalid("finite nonnegative source gain"));
        }
        hasher.update(&(source_index as u64).to_le_bytes());
        hasher.update(source.source_identity.as_bytes());
        hasher.update(&[source.authority.tag()]);
        hash_f64(&mut hasher, source.gain_linear);
        source.signal.hash_header(&mut hasher);
        if retarded {
            hasher.update(&[match source.positions {
                SourcePositionTrack::Static(_) => 1,
                SourcePositionTrack::PerFrame(_) => 2,
            }]);
        }
        for frame in 0..frame_count {
            if frame % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let sample = source.signal.sample(frame);
            validate_finite(sample, "spatial source sample", frame)?;
            hash_f64(&mut hasher, sample);
            let source_position = source.positions.position(frame);
            hash_vec3(&mut hasher, source_position);
            validate_vec3(source_position, "source position", frame)?;
            let listener = input.listener.pose(frame);
            validate_listener(listener, frame)?;
            if source_index == 0 {
                hash_listener(&mut hasher, listener);
            }
            if retarded {
                continue;
            }
            let geometry = geometry(
                source_position,
                listener,
                spatializer.config.minimum_distance_m,
                spatializer.config.speed_of_sound_m_per_s,
                spatializer.config.sample_rate_hz,
                spatializer.config.microphone_directivity,
                frame,
            )?;
            maximum_distance_m = maximum_distance_m.max(geometry.distance_m);
            maximum_delay_frames = maximum_delay_frames.max(geometry.delay_frames);
            if geometry.distance_m < spatializer.config.minimum_distance_m {
                minimum_distance_clamp_count = minimum_distance_clamp_count
                    .checked_add(1)
                    .ok_or_else(|| invalid("minimum-distance clamp count"))?;
            }
            let arrival_offset =
                delay_tail_offset(geometry.delay_frames, spatializer.config.delay_policy)?;
            maximum_arrival_offset = maximum_arrival_offset.max(arrival_offset);
        }
    }
    if retarded {
        // The finite receiver tail is part of the actual physical input, not a
        // held endpoint. Charge its work before traversing and identity-bind it.
        check_limit(
            "receiver pose frames",
            listener_frames as u64,
            spatializer.config.budget.maximum_total_input_frames,
        )?;
        hasher.update(&(listener_frames as u64).to_le_bytes());
        hasher.update(&[match input.listener {
            ListenerPoseTrack::Static(_) => 1,
            ListenerPoseTrack::PerFrame(_) => 2,
        }]);
        for frame in frame_count..listener_frames {
            if frame % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let listener = input.listener.pose(frame);
            validate_listener(listener, frame)?;
            hash_listener(&mut hasher, listener);
        }
    }
    if let Some(room_ir) = input.room_ir {
        authority = authority.combine(room_ir.authority);
        hasher.update(&[1]);
        hasher.update(room_ir.identity.as_bytes());
    } else {
        hasher.update(&[0]);
    }
    checkpoint_fn()?;
    let input_identity = hasher.finalize();

    let (natural_direct_output_frames, refinement_work) = if retarded {
        retarded::admit(spatializer.config, input, frame_count, checkpoint_fn)?
    } else {
        (
            (frame_count as u64)
                .checked_add(maximum_arrival_offset)
                .ok_or_else(|| invalid("direct output frame count"))?,
            0,
        )
    };
    let room_tail = input
        .room_ir
        .map_or(0_u64, |room_ir| room_ir.len() as u64 - 1);
    let natural_final_output_frames = natural_direct_output_frames
        .checked_add(room_tail)
        .ok_or_else(|| invalid("final output frame count"))?;
    let (direct_output_frames, final_output_frames) = match spatializer.config.output_horizon {
        SpatialOutputHorizon::PreserveTail => {
            (natural_direct_output_frames, natural_final_output_frames)
        }
        SpatialOutputHorizon::ClampToInputFrames => {
            let horizon = frame_count as u64;
            (horizon, horizon)
        }
    };
    let discarded_tail_frames = natural_final_output_frames.saturating_sub(final_output_frames);
    check_limit(
        "output frames",
        final_output_frames,
        spatializer.config.budget.maximum_output_frames,
    )?;
    validate_memory_budget(
        direct_output_frames,
        final_output_frames,
        input.room_ir.is_some(),
        spatializer.config.budget,
    )?;
    validate_work_budget(
        total_input_frames,
        direct_output_frames,
        final_output_frames,
        input.room_ir.map_or(0, StereoRoomImpulseResponse::len),
        if retarded {
            retarded_work_units(
                input.sources.len() as u64,
                direct_output_frames,
                listener_frames as u64,
                refinement_work,
            )?
        } else {
            0
        },
        spatializer.config.budget,
    )?;
    checkpoint_fn()?;

    Ok(Preflight {
        frame_count,
        direct_output_frames,
        final_output_frames,
        natural_final_output_frames,
        discarded_tail_frames,
        maximum_distance_m,
        maximum_delay_frames,
        minimum_distance_clamp_count,
        authority,
        input_identity,
    })
}

#[derive(Debug, Clone, Copy)]
struct SpatialGeometry {
    distance_m: f64,
    delay_frames: f64,
    attenuation_gain: f64,
    microphone_gain: f64,
    pan: f64,
}

fn geometry(
    source_position: [f64; 3],
    listener: ListenerPose,
    minimum_distance_m: f64,
    speed_of_sound_m_per_s: f64,
    sample_rate_hz: u32,
    directivity: MicrophoneDirectivity,
    frame: usize,
) -> Result<SpatialGeometry, SpatialAudioError> {
    let delta = [
        source_position[0] - listener.position_m[0],
        source_position[1] - listener.position_m[1],
        source_position[2] - listener.position_m[2],
    ];
    validate_vec3(delta, "source-listener displacement", frame)?;
    let distance_squared = dot(delta, delta);
    validate_finite(distance_squared, "source-listener squared distance", frame)?;
    let distance_m = det::sqrt(distance_squared);
    validate_finite(distance_m, "source-listener distance", frame)?;
    let delay_frames = distance_m * sample_rate_hz as f64 / speed_of_sound_m_per_s;
    validate_finite(delay_frames, "propagation delay", frame)?;
    if delay_frames < 0.0 || delay_frames > u64::MAX as f64 {
        return Err(invalid("addressable propagation delay"));
    }
    let attenuation_gain = minimum_distance_m / distance_m.max(minimum_distance_m);
    let (pan, forward_cosine) = if distance_m == 0.0 {
        // Direction is undefined at coincidence. V1 freezes a centred, on-axis
        // response while the attenuation law uses its finite distance clamp.
        (0.0, 1.0)
    } else {
        let inverse_distance = 1.0 / distance_m;
        let direction = [
            delta[0] * inverse_distance,
            delta[1] * inverse_distance,
            delta[2] * inverse_distance,
        ];
        (
            dot(direction, listener.right_unit).clamp(-1.0, 1.0),
            dot(direction, listener.forward_unit).clamp(-1.0, 1.0),
        )
    };
    let microphone_gain = directivity.gain(forward_cosine);
    validate_finite(microphone_gain, "microphone gain", frame)?;
    Ok(SpatialGeometry {
        distance_m,
        delay_frames,
        attenuation_gain,
        microphone_gain,
        pan,
    })
}

fn delay_tail_offset(
    delay_frames: f64,
    policy: SpatialDelayPolicy,
) -> Result<u64, SpatialAudioError> {
    let offset = match policy {
        SpatialDelayPolicy::IntegerCeiling => delay_frames.ceil(),
        SpatialDelayPolicy::LinearFloorCeil => {
            let floor = delay_frames.floor();
            if delay_frames.to_bits() == floor.to_bits() {
                floor
            } else {
                floor + 1.0
            }
        }
        SpatialDelayPolicy::RetardedSinc { .. }
        | SpatialDelayPolicy::RetardedSincUniformFlow { .. } => {
            return Err(invalid("retarded delays use receiver-clock reconstruction"));
        }
    };
    if !offset.is_finite() || offset < 0.0 || offset > u64::MAX as f64 {
        return Err(invalid("sampled delay offset"));
    }
    Ok(offset as u64)
}

fn deposit(
    output: &mut [StereoAccumulator],
    emission_frame: usize,
    delay_frames: f64,
    left: f64,
    right: f64,
    policy: SpatialDelayPolicy,
    output_horizon: SpatialOutputHorizon,
    source_index: usize,
) -> Result<(), SpatialAudioError> {
    validate_finite(left, "spatialized left contribution", emission_frame)?;
    validate_finite(right, "spatialized right contribution", emission_frame)?;
    match policy {
        SpatialDelayPolicy::IntegerCeiling => {
            let offset = delay_tail_offset(delay_frames, policy)?;
            let index = emission_frame
                .checked_add(
                    usize::try_from(offset)
                        .map_err(|_| invalid("addressable integer propagation-delay offset"))?,
                )
                .ok_or_else(|| invalid("integer propagation-delay output index"))?;
            let Some(accumulator) = output.get_mut(index) else {
                if matches!(output_horizon, SpatialOutputHorizon::ClampToInputFrames) {
                    return Ok(());
                }
                return Err(SpatialAudioError::ResourceLimit {
                    resource: "integer propagation-delay output index",
                    requested: index as u64,
                    limit: output.len().saturating_sub(1) as u64,
                });
            };
            accumulator.add(left, right);
        }
        SpatialDelayPolicy::LinearFloorCeil => {
            let lower_offset = delay_frames.floor();
            if lower_offset > usize::MAX as f64 {
                return Err(invalid("addressable fractional-delay offset"));
            }
            let lower = emission_frame
                .checked_add(lower_offset as usize)
                .ok_or_else(|| invalid("fractional-delay lower output index"))?;
            let fraction = delay_frames - lower_offset;
            let lower_weight = 1.0 - fraction;
            let Some(lower_accumulator) = output.get_mut(lower) else {
                if matches!(output_horizon, SpatialOutputHorizon::ClampToInputFrames) {
                    return Ok(());
                }
                return Err(SpatialAudioError::ResourceLimit {
                    resource: "fractional-delay lower output index",
                    requested: lower as u64,
                    limit: output.len().saturating_sub(1) as u64,
                });
            };
            lower_accumulator.add(left * lower_weight, right * lower_weight);
            if fraction != 0.0 {
                let upper = lower
                    .checked_add(1)
                    .ok_or_else(|| invalid("fractional-delay upper output index"))?;
                let Some(upper_accumulator) = output.get_mut(upper) else {
                    if matches!(output_horizon, SpatialOutputHorizon::ClampToInputFrames) {
                        return Ok(());
                    }
                    return Err(SpatialAudioError::ResourceLimit {
                        resource: "fractional-delay upper output index",
                        requested: upper as u64,
                        limit: output.len().saturating_sub(1) as u64,
                    });
                };
                upper_accumulator.add(left * fraction, right * fraction);
            }
        }
        SpatialDelayPolicy::RetardedSinc { .. }
        | SpatialDelayPolicy::RetardedSincUniformFlow { .. } => {
            return Err(invalid("retarded signals cannot be deposited"));
        }
    }
    let _ = source_index;
    Ok(())
}

fn convolve_room(
    direct: &[StereoSample],
    final_len: usize,
    room_ir: &StereoRoomImpulseResponse,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<Vec<StereoSample>, SpatialAudioError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(final_len)
        .map_err(|_| SpatialAudioError::AllocationFailed("room-response output"))?;
    for output_index in 0..final_len {
        if output_index % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        let first_tap = output_index.saturating_sub(direct.len().saturating_sub(1));
        let last_tap = output_index.min(room_ir.len() - 1);
        let mut left = CompensatedSum::ZERO;
        let mut right = CompensatedSum::ZERO;
        if first_tap <= last_tap {
            for tap in first_tap..=last_tap {
                if tap % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
                    checkpoint_fn()?;
                }
                let source_index = output_index - tap;
                left.add(direct[source_index].left_fs * room_ir.left_taps[tap]);
                right.add(direct[source_index].right_fs * room_ir.right_taps[tap]);
            }
        }
        output.push(StereoSample {
            left_fs: left.finish(),
            right_fs: right.finish(),
        });
    }
    checkpoint_fn()?;
    Ok(output)
}

fn validate_final_samples(
    samples: &[StereoSample],
    limit_fs: f64,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<f64, SpatialAudioError> {
    let mut peak = 0.0_f64;
    for (index, sample) in samples.iter().copied().enumerate() {
        if index % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        validate_finite(sample.left_fs, "final left sample", index)?;
        validate_finite(sample.right_fs, "final right sample", index)?;
        check_amplitude(sample.left_fs, index, "left", limit_fs)?;
        check_amplitude(sample.right_fs, index, "right", limit_fs)?;
        peak = peak.max(sample.left_fs.abs()).max(sample.right_fs.abs());
    }
    checkpoint_fn()?;
    Ok(peak)
}

fn output_identity(
    config_identity: ContentHash,
    input_identity: ContentHash,
    sample_rate_hz: u32,
    authority: SpatialAudioAuthority,
    room_ir_identity: Option<ContentHash>,
    samples: &[StereoSample],
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<ContentHash, SpatialAudioError> {
    let mut hasher = DomainHasher::new(OUTPUT_IDENTITY_DOMAIN);
    hasher.update(&SPATIAL_AUDIO_ALGORITHM_VERSION.to_le_bytes());
    hasher.update(config_identity.as_bytes());
    hasher.update(input_identity.as_bytes());
    hasher.update(&sample_rate_hz.to_le_bytes());
    hasher.update(&[authority.tag()]);
    match room_ir_identity {
        Some(identity) => {
            hasher.update(&[1]);
            hasher.update(identity.as_bytes());
        }
        None => hasher.update(&[0]),
    }
    hasher.update(&(samples.len() as u64).to_le_bytes());
    for (index, sample) in samples.iter().copied().enumerate() {
        if index % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        hash_f64(&mut hasher, sample.left_fs);
        hash_f64(&mut hasher, sample.right_fs);
    }
    checkpoint_fn()?;
    let identity = hasher.finalize();
    if is_zero_hash(identity) {
        return Err(SpatialAudioError::InvalidIdentity("spatial-audio output"));
    }
    Ok(identity)
}

fn validate_config(config: SpatialAudioConfig) -> Result<(), SpatialAudioError> {
    validate_sample_rate(config.sample_rate_hz)?;
    if !config.speed_of_sound_m_per_s.is_finite() || config.speed_of_sound_m_per_s <= 0.0 {
        return Err(invalid("positive finite speed of sound [m/s]"));
    }
    if !config.minimum_distance_m.is_finite() || config.minimum_distance_m <= 0.0 {
        return Err(invalid("positive finite minimum distance [m]"));
    }
    if let Some((source_bandwidth_hz, receiver_bandwidth_hz)) = config.delay_policy.retarded_bands()
    {
        let nyquist = 0.5 * f64::from(config.sample_rate_hz);
        if !source_bandwidth_hz.is_finite()
            || source_bandwidth_hz < 0.0
            || source_bandwidth_hz > 0.5 * nyquist
            || !receiver_bandwidth_hz.is_finite()
            || receiver_bandwidth_hz < 0.0
            || receiver_bandwidth_hz >= nyquist
        {
            return Err(invalid(
                "retarded source band at most fs/4 and receiver band below Nyquist",
            ));
        }
        validate_vec3(
            config.delay_policy.medium_velocity_m_per_s(),
            "uniform medium velocity [m/s]",
            0,
        )?;
    }
    if let MicrophoneDirectivity::Cardioid { rear_floor_gain } = config.microphone_directivity
        && (!rear_floor_gain.is_finite() || !(0.0..=1.0).contains(&rear_floor_gain))
    {
        return Err(invalid("cardioid rear-floor gain in [0,1]"));
    }
    validate_budget(config.budget)
}

fn validate_sample_rate(sample_rate_hz: u32) -> Result<(), SpatialAudioError> {
    if sample_rate_hz == 0 || sample_rate_hz > MAX_SPATIAL_AUDIO_SAMPLE_RATE_HZ {
        return Err(invalid("sample rate in (0,384000] Hz"));
    }
    Ok(())
}

fn validate_budget(budget: SpatialAudioBudget) -> Result<(), SpatialAudioError> {
    if budget.maximum_sources == 0 || budget.maximum_sources > MAX_SPATIAL_AUDIO_SOURCES {
        return Err(invalid("source budget in [1,32]"));
    }
    if budget.maximum_total_input_frames == 0 {
        return Err(invalid("nonzero total-input-frame budget"));
    }
    if budget.maximum_output_frames == 0 {
        return Err(invalid("nonzero output-frame budget"));
    }
    if budget.maximum_room_ir_taps == 0
        || budget.maximum_room_ir_taps > MAX_SPATIAL_AUDIO_ROOM_IR_TAPS
    {
        return Err(invalid("room-response tap budget in [1,65536]"));
    }
    if budget.maximum_work_units == 0 {
        return Err(invalid("nonzero work budget"));
    }
    if budget.maximum_owned_sample_bytes == 0 {
        return Err(invalid("nonzero owned-sample byte budget"));
    }
    if !budget.maximum_abs_output_fs.is_finite() || budget.maximum_abs_output_fs <= 0.0 {
        return Err(invalid("positive finite output-amplitude ceiling"));
    }
    Ok(())
}

fn validate_listener(listener: ListenerPose, frame: usize) -> Result<(), SpatialAudioError> {
    validate_vec3(listener.position_m, "listener position", frame)?;
    validate_vec3(listener.forward_unit, "listener forward axis", frame)?;
    validate_vec3(listener.right_unit, "listener right axis", frame)?;
    let forward_norm_squared = dot(listener.forward_unit, listener.forward_unit);
    let right_norm_squared = dot(listener.right_unit, listener.right_unit);
    if (forward_norm_squared - 1.0).abs() > POSE_UNIT_TOLERANCE {
        return Err(SpatialAudioError::InvalidListenerPose {
            frame,
            reason: "forward axis is not unit length",
        });
    }
    if (right_norm_squared - 1.0).abs() > POSE_UNIT_TOLERANCE {
        return Err(SpatialAudioError::InvalidListenerPose {
            frame,
            reason: "right axis is not unit length",
        });
    }
    if dot(listener.forward_unit, listener.right_unit).abs() > POSE_UNIT_TOLERANCE {
        return Err(SpatialAudioError::InvalidListenerPose {
            frame,
            reason: "forward and right axes are not orthogonal",
        });
    }
    Ok(())
}

fn validate_vec3(
    vector: [f64; 3],
    field: &'static str,
    index: usize,
) -> Result<(), SpatialAudioError> {
    for component in vector {
        validate_finite(component, field, index)?;
    }
    Ok(())
}

fn validate_finite(value: f64, field: &'static str, index: usize) -> Result<(), SpatialAudioError> {
    if !value.is_finite() {
        return Err(SpatialAudioError::NonFinite { field, index });
    }
    Ok(())
}

fn validate_memory_budget(
    direct_frames: u64,
    final_frames: u64,
    has_room: bool,
    budget: SpatialAudioBudget,
) -> Result<(), SpatialAudioError> {
    let direct_accumulators = direct_frames
        .checked_mul(core::mem::size_of::<StereoAccumulator>() as u64)
        .ok_or_else(|| invalid("direct accumulator byte count"))?;
    let direct_samples = direct_frames
        .checked_mul(core::mem::size_of::<StereoSample>() as u64)
        .ok_or_else(|| invalid("direct sample byte count"))?;
    let final_samples = if has_room {
        final_frames
            .checked_mul(core::mem::size_of::<StereoSample>() as u64)
            .ok_or_else(|| invalid("room output byte count"))?
    } else {
        0
    };
    let requested = direct_accumulators
        .checked_add(direct_samples)
        .and_then(|value| value.checked_add(final_samples))
        .ok_or_else(|| invalid("aggregate owned sample byte count"))?;
    check_limit(
        "owned sample bytes",
        requested,
        budget.maximum_owned_sample_bytes,
    )
}

fn retarded_work_units(
    sources: u64,
    direct_frames: u64,
    listener_frames: u64,
    refinement_work: u64,
) -> Result<u64, SpatialAudioError> {
    sources
        .checked_mul(direct_frames)
        .and_then(|n| n.checked_mul(96))
        .and_then(|n| n.checked_add(listener_frames))
        .and_then(|n| n.checked_add(refinement_work))
        .ok_or_else(|| invalid("retarded reconstruction and admission work"))
}

fn validate_work_budget(
    total_input_frames: u64,
    direct_frames: u64,
    final_frames: u64,
    ir_taps: usize,
    additional_work: u64,
    budget: SpatialAudioBudget,
) -> Result<(), SpatialAudioError> {
    let preflight_and_render = total_input_frames
        .checked_mul(3)
        .ok_or_else(|| invalid("spatial source work"))?;
    let sample_finalize_and_hash = direct_frames
        .checked_add(
            final_frames
                .checked_mul(2)
                .ok_or_else(|| invalid("final sample work"))?,
        )
        .ok_or_else(|| invalid("sample finalization work"))?;
    let convolution = if ir_taps == 0 {
        0
    } else {
        direct_frames
            .checked_mul(ir_taps as u64)
            .and_then(|value| value.checked_mul(2))
            .ok_or_else(|| invalid("room convolution work"))?
    };
    let requested = preflight_and_render
        .checked_add(sample_finalize_and_hash)
        .and_then(|value| value.checked_add(convolution))
        .and_then(|value| value.checked_add(additional_work))
        .ok_or_else(|| invalid("aggregate spatial-audio work"))?;
    check_limit("work units", requested, budget.maximum_work_units)
}

fn config_identity(config: SpatialAudioConfig) -> ContentHash {
    let mut hasher = DomainHasher::new(CONFIG_IDENTITY_DOMAIN);
    hasher.update(&SPATIAL_AUDIO_ALGORITHM_VERSION.to_le_bytes());
    hasher.update(&config.sample_rate_hz.to_le_bytes());
    hash_f64(&mut hasher, config.speed_of_sound_m_per_s);
    hash_f64(&mut hasher, config.minimum_distance_m);
    hasher.update(&[config.delay_policy.tag()]);
    if let Some((source_bandwidth_hz, receiver_bandwidth_hz)) = config.delay_policy.retarded_bands()
    {
        hash_f64(&mut hasher, source_bandwidth_hz);
        hash_f64(&mut hasher, receiver_bandwidth_hz);
        hash_vec3(&mut hasher, config.delay_policy.medium_velocity_m_per_s());
    }
    hasher.update(&[config.output_horizon.tag()]);
    config.microphone_directivity.hash_into(&mut hasher);
    hasher.update(&[config.authority.tag()]);
    hash_budget(&mut hasher, config.budget);
    hasher.finalize()
}

fn dry_bypass_config_identity(sample_rate_hz: u32, budget: SpatialAudioBudget) -> ContentHash {
    let mut hasher = DomainHasher::new(DRY_BYPASS_CONFIG_IDENTITY_DOMAIN);
    hasher.update(&SPATIAL_AUDIO_ALGORITHM_VERSION.to_le_bytes());
    hasher.update(&sample_rate_hz.to_le_bytes());
    hash_budget(&mut hasher, budget);
    hasher.finalize()
}

fn hash_budget(hasher: &mut DomainHasher, budget: SpatialAudioBudget) {
    hasher.update(&(budget.maximum_sources as u64).to_le_bytes());
    hasher.update(&budget.maximum_total_input_frames.to_le_bytes());
    hasher.update(&budget.maximum_output_frames.to_le_bytes());
    hasher.update(&(budget.maximum_room_ir_taps as u64).to_le_bytes());
    hasher.update(&budget.maximum_work_units.to_le_bytes());
    hasher.update(&budget.maximum_owned_sample_bytes.to_le_bytes());
    hash_f64(hasher, budget.maximum_abs_output_fs);
}

fn hash_listener(hasher: &mut DomainHasher, listener: ListenerPose) {
    hash_vec3(hasher, listener.position_m);
    hash_vec3(hasher, listener.forward_unit);
    hash_vec3(hasher, listener.right_unit);
}

fn hash_vec3(hasher: &mut DomainHasher, vector: [f64; 3]) {
    for value in vector {
        hash_f64(hasher, value);
    }
}

fn hash_f64(hasher: &mut DomainHasher, value: f64) {
    hasher.update(&value.to_bits().to_le_bytes());
}

fn check_amplitude(
    sample: f64,
    frame: usize,
    channel: &'static str,
    limit_fs: f64,
) -> Result<(), SpatialAudioError> {
    if sample.abs() > limit_fs {
        return Err(SpatialAudioError::OutputAmplitudeExceeded {
            frame,
            channel,
            magnitude_fs: sample.abs(),
            limit_fs,
        });
    }
    Ok(())
}

fn check_limit(
    resource: &'static str,
    requested: u64,
    admitted: u64,
) -> Result<(), SpatialAudioError> {
    if requested > admitted {
        return Err(limit(resource, requested, admitted));
    }
    Ok(())
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0].mul_add(right[0], left[1].mul_add(right[1], left[2] * right[2]))
}

fn is_zero_hash(identity: ContentHash) -> bool {
    identity.as_bytes().iter().all(|byte| *byte == 0)
}

fn invalid(field: &'static str) -> SpatialAudioError {
    SpatialAudioError::InvalidConfig(field)
}

fn limit(resource: &'static str, requested: u64, admitted: u64) -> SpatialAudioError {
    SpatialAudioError::ResourceLimit {
        resource,
        requested,
        limit: admitted,
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SpatialAudioError> {
    cx.checkpoint().map_err(|_| SpatialAudioError::Cancelled)
}

#[derive(Debug, Clone, Copy)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    const ZERO: Self = Self {
        sum: 0.0,
        correction: 0.0,
    };

    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn finish(self) -> f64 {
        self.sum + self.correction
    }
}

#[derive(Debug, Clone, Copy)]
struct StereoAccumulator {
    left: CompensatedSum,
    right: CompensatedSum,
}

impl StereoAccumulator {
    const ZERO: Self = Self {
        left: CompensatedSum::ZERO,
        right: CompensatedSum::ZERO,
    };

    fn add(&mut self, left: f64, right: f64) {
        self.left.add(left);
        self.right.add(right);
    }

    fn finish(self) -> StereoSample {
        StereoSample {
            left_fs: self.left.finish(),
            right_fs: self.right.finish(),
        }
    }
}

#[cfg(test)]
mod tests {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_blake3::hash_domain;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    use super::*;

    fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x5350_4154_4941_4c41,
                    kernel_id: 0x4555_4c45,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn identity(label: &str) -> ContentHash {
        hash_domain("org.frankensim.test.spatial-audio.v1", label.as_bytes())
    }

    fn budget() -> SpatialAudioBudget {
        SpatialAudioBudget {
            maximum_sources: 8,
            maximum_total_input_frames: 100_000,
            maximum_output_frames: 100_000,
            maximum_room_ir_taps: 1_024,
            maximum_work_units: 10_000_000,
            maximum_owned_sample_bytes: 16 * 1024 * 1024,
            maximum_abs_output_fs: 10.0,
        }
    }

    fn config(
        sample_rate_hz: u32,
        speed_of_sound_m_per_s: f64,
        policy: SpatialDelayPolicy,
    ) -> SpatialAudioConfig {
        SpatialAudioConfig {
            sample_rate_hz,
            speed_of_sound_m_per_s,
            minimum_distance_m: 1.0,
            delay_policy: policy,
            output_horizon: SpatialOutputHorizon::PreserveTail,
            microphone_directivity: MicrophoneDirectivity::Omnidirectional,
            authority: SpatialAudioAuthority::PhysicallyParameterized,
            budget: budget(),
        }
    }

    fn listener_at(position_m: [f64; 3]) -> ListenerPose {
        ListenerPose {
            position_m,
            forward_unit: [0.0, 0.0, 1.0],
            right_unit: [1.0, 0.0, 0.0],
        }
    }

    fn render_single(
        config: SpatialAudioConfig,
        samples: &[f64],
        positions: SourcePositionTrack<'_>,
        listener: ListenerPoseTrack<'_>,
        room_ir: Option<&StereoRoomImpulseResponse>,
    ) -> Result<SpatialAudioOutput, SpatialAudioError> {
        render_single_with_gain(config, samples, positions, listener, room_ir, 1.0)
    }

    fn render_single_with_gain(
        config: SpatialAudioConfig,
        samples: &[f64],
        positions: SourcePositionTrack<'_>,
        listener: ListenerPoseTrack<'_>,
        room_ir: Option<&StereoRoomImpulseResponse>,
        gain_linear: f64,
    ) -> Result<SpatialAudioOutput, SpatialAudioError> {
        with_cx(false, |cx| {
            let spatializer = OfflineSpatializer::try_new(config, cx)?;
            let source = SpatialAudioSource {
                source_identity: identity("source"),
                signal: SpatialMonoSignal::Samples(samples),
                positions,
                gain_linear,
                authority: SpatialAudioAuthority::PhysicallyParameterized,
            };
            spatializer.spatialize(
                SpatialAudioRenderInput {
                    sources: &[source],
                    listener,
                    room_ir,
                },
                cx,
            )
        })
    }

    fn first_nonzero(samples: &[StereoSample]) -> usize {
        samples
            .iter()
            .position(|sample| sample.left_fs != 0.0 || sample.right_fs != 0.0)
            .expect("impulse must arrive")
    }

    #[test]
    fn static_left_right_and_center_use_frozen_equal_power_pan() {
        let cfg = config(8, 8.0, SpatialDelayPolicy::IntegerCeiling);
        let listener = ListenerPoseTrack::Static(listener_at([0.0; 3]));
        let left = render_single(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([-1.0, 0.0, 0.0]),
            listener,
            None,
        )
        .unwrap();
        let right = render_single(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([1.0, 0.0, 0.0]),
            listener,
            None,
        )
        .unwrap();
        let center = render_single(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([0.0, 0.0, 1.0]),
            listener,
            None,
        )
        .unwrap();
        assert_eq!(left.samples()[1].right_fs, 0.0);
        assert!(left.samples()[1].left_fs > 0.99);
        assert_eq!(right.samples()[1].left_fs, 0.0);
        assert!(right.samples()[1].right_fs > 0.99);
        assert_eq!(center.samples()[1].left_fs, center.samples()[1].right_fs);
    }

    #[test]
    fn near_far_and_speed_of_sound_scale_gain_and_delay() {
        let near = render_single(
            config(8, 8.0, SpatialDelayPolicy::IntegerCeiling),
            &[1.0],
            SourcePositionTrack::Static([0.0, 0.0, 1.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap();
        let far = render_single(
            config(8, 8.0, SpatialDelayPolicy::IntegerCeiling),
            &[1.0],
            SourcePositionTrack::Static([0.0, 0.0, 2.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap();
        let faster = render_single(
            config(8, 16.0, SpatialDelayPolicy::IntegerCeiling),
            &[1.0],
            SourcePositionTrack::Static([0.0, 0.0, 2.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap();
        assert_eq!(first_nonzero(near.samples()), 1);
        assert_eq!(first_nonzero(far.samples()), 2);
        assert_eq!(first_nonzero(faster.samples()), 1);
        assert!(near.diagnostics().sample_peak_fs > far.diagnostics().sample_peak_fs);
    }

    #[test]
    fn source_gain_is_applied_once_and_bound_into_input_identity() {
        let cfg = config(8, 8.0, SpatialDelayPolicy::IntegerCeiling);
        let unity = render_single_with_gain(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([-1.0, 0.0, 0.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
            1.0,
        )
        .unwrap();
        let quarter = render_single_with_gain(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([-1.0, 0.0, 0.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
            0.25,
        )
        .unwrap();
        assert_eq!(
            quarter.samples()[1].left_fs,
            unity.samples()[1].left_fs * 0.25
        );
        assert_ne!(quarter.input_identity(), unity.input_identity());
        assert_ne!(quarter.output_identity(), unity.output_identity());
    }

    #[test]
    fn clamped_horizon_is_identity_bound_and_reports_discarded_delay_and_ir_tail() {
        with_cx(false, |cx| {
            let room = StereoRoomImpulseResponse::try_new(
                8,
                vec![1.0, 0.5, 0.25],
                vec![1.0, 0.5, 0.25],
                SpatialAudioAuthority::Artistic,
                cx,
            )
            .unwrap();
            let preserve_cfg = config(8, 8.0, SpatialDelayPolicy::IntegerCeiling);
            let preserved = render_single(
                preserve_cfg,
                &[1.0],
                SourcePositionTrack::Static([0.0, 0.0, 1.0]),
                ListenerPoseTrack::Static(listener_at([0.0; 3])),
                Some(&room),
            )
            .unwrap();
            let mut clamp_cfg = preserve_cfg;
            clamp_cfg.output_horizon = SpatialOutputHorizon::ClampToInputFrames;
            let clamped = render_single(
                clamp_cfg,
                &[1.0],
                SourcePositionTrack::Static([0.0, 0.0, 1.0]),
                ListenerPoseTrack::Static(listener_at([0.0; 3])),
                Some(&room),
            )
            .unwrap();
            assert_eq!(preserved.samples().len(), 4);
            assert_eq!(clamped.samples().len(), 1);
            assert_eq!(clamped.diagnostics().natural_final_output_frames, 4);
            assert_eq!(clamped.diagnostics().discarded_tail_frames, 3);
            assert_eq!(clamped.diagnostics().final_output_frames, 1);
            assert_ne!(clamped.config_identity(), preserved.config_identity());
            assert_ne!(clamped.output_identity(), preserved.output_identity());
        });
    }

    #[test]
    fn listener_motion_changes_pan_at_sample_synchronous_emission_poses() {
        let poses = [
            listener_at([-0.6, 0.0, -0.8]),
            listener_at([0.6, 0.0, -0.8]),
        ];
        let output = render_single(
            config(8, 8.0, SpatialDelayPolicy::IntegerCeiling),
            &[1.0, 1.0],
            SourcePositionTrack::Static([0.0; 3]),
            ListenerPoseTrack::PerFrame(&poses),
            None,
        )
        .unwrap();
        assert!(output.samples()[1].right_fs > output.samples()[1].left_fs);
        assert!(output.samples()[2].left_fs > output.samples()[2].right_fs);
    }

    #[test]
    fn coincidence_is_finite_centered_and_uses_minimum_distance_clamp() {
        let output = render_single(
            config(8, 8.0, SpatialDelayPolicy::LinearFloorCeil),
            &[1.0],
            SourcePositionTrack::Static([0.0; 3]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap();
        assert_eq!(output.samples().len(), 1);
        assert_eq!(output.samples()[0].left_fs, output.samples()[0].right_fs);
        assert_eq!(output.diagnostics().minimum_distance_clamp_count, 1);
        assert_eq!(output.diagnostics().maximum_delay_frames, 0.0);
    }

    #[test]
    fn cardioid_rejects_rear_axis_without_claiming_hrtf() {
        let mut cfg = config(8, 8.0, SpatialDelayPolicy::IntegerCeiling);
        cfg.microphone_directivity = MicrophoneDirectivity::Cardioid {
            rear_floor_gain: 0.0,
        };
        let front = render_single(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([0.0, 0.0, 1.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap();
        let rear = render_single(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([0.0, 0.0, -1.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap();
        assert!(front.diagnostics().sample_peak_fs > 0.7);
        assert_eq!(rear.diagnostics().sample_peak_fs, 0.0);
    }

    #[test]
    fn fractional_delay_and_room_ir_retain_impulse_and_declared_tail() {
        with_cx(false, |cx| {
            let room = StereoRoomImpulseResponse::try_new(
                8,
                vec![1.0, 0.5],
                vec![1.0, -0.5],
                SpatialAudioAuthority::Artistic,
                cx,
            )
            .unwrap();
            let output = render_single(
                config(8, 16.0, SpatialDelayPolicy::LinearFloorCeil),
                &[1.0],
                SourcePositionTrack::Static([0.0, 0.0, 1.0]),
                ListenerPoseTrack::Static(listener_at([0.0; 3])),
                Some(&room),
            )
            .unwrap();
            assert_eq!(output.samples().len(), 3);
            assert_eq!(output.authority(), SpatialAudioAuthority::Artistic);
            assert_eq!(output.room_ir_identity(), Some(room.identity()));
            assert!(output.samples()[0].left_fs > 0.0);
            assert!(output.samples()[1].left_fs > 0.0);
            assert!(output.samples()[2].left_fs > 0.0);
            assert!(output.samples()[2].right_fs < 0.0);
        });
    }

    #[test]
    fn room_ir_rate_must_match_and_identity_changes_with_taps() {
        with_cx(false, |cx| {
            let first = StereoRoomImpulseResponse::try_new(
                16,
                vec![1.0],
                vec![1.0],
                SpatialAudioAuthority::PhysicallyParameterized,
                cx,
            )
            .unwrap();
            let changed = StereoRoomImpulseResponse::try_new(
                16,
                vec![0.5],
                vec![1.0],
                SpatialAudioAuthority::PhysicallyParameterized,
                cx,
            )
            .unwrap();
            assert_ne!(first.identity(), changed.identity());
            let error = render_single(
                config(8, 8.0, SpatialDelayPolicy::IntegerCeiling),
                &[1.0],
                SourcePositionTrack::Static([0.0, 0.0, 1.0]),
                ListenerPoseTrack::Static(listener_at([0.0; 3])),
                Some(&first),
            )
            .unwrap_err();
            assert_eq!(
                error,
                SpatialAudioError::InvalidConfig("room impulse response sample rate")
            );
        });
    }

    #[test]
    fn modal_stem_selection_is_direct_and_replay_is_bit_stable() {
        let frames = [
            ModalStemFrame {
                disc_fs: 0.25,
                glass_plate_fs: 0.5,
                base_assembly_fs: -0.75,
            },
            ModalStemFrame {
                disc_fs: -0.125,
                glass_plate_fs: 0.75,
                base_assembly_fs: 0.5,
            },
        ];
        with_cx(false, |cx| {
            let spatializer =
                OfflineSpatializer::try_new(config(8, 8.0, SpatialDelayPolicy::IntegerCeiling), cx)
                    .unwrap();
            let source = SpatialAudioSource {
                source_identity: identity("modal-stem"),
                signal: SpatialMonoSignal::ModalStemFrames {
                    frames: &frames,
                    component: SpatialStemComponent::GlassPlate,
                },
                positions: SourcePositionTrack::Static([0.0, 0.0, 1.0]),
                gain_linear: 1.0,
                authority: SpatialAudioAuthority::Artistic,
            };
            let input = SpatialAudioRenderInput {
                sources: &[source],
                listener: ListenerPoseTrack::Static(listener_at([0.0; 3])),
                room_ir: None,
            };
            let first = spatializer.spatialize(input, cx).unwrap();
            let replay = spatializer.spatialize(input, cx).unwrap();
            assert_eq!(first, replay);
            assert_eq!(first.authority(), SpatialAudioAuthority::Artistic);
        });
    }

    #[test]
    fn dry_bypass_is_bit_exact_including_signed_zero() {
        let frames = [
            StereoSample {
                left_fs: -0.0,
                right_fs: 0.25,
            },
            StereoSample {
                left_fs: -0.5,
                right_fs: 0.0,
            },
        ];
        let output = with_cx(false, |cx| {
            bypass_dry_stereo(
                &frames,
                identity("dry"),
                48_000,
                SpatialAudioAuthority::PhysicallyParameterized,
                budget(),
                cx,
            )
        })
        .unwrap();
        for (actual, expected) in output.samples().iter().zip(frames) {
            assert_eq!(actual.left_fs.to_bits(), expected.left_fs.to_bits());
            assert_eq!(actual.right_fs.to_bits(), expected.right_fs.to_bits());
        }
    }

    #[test]
    fn nonfinite_and_output_clipping_refuse_the_whole_transaction() {
        let nonfinite = render_single(
            config(8, 8.0, SpatialDelayPolicy::IntegerCeiling),
            &[f64::NAN],
            SourcePositionTrack::Static([0.0, 0.0, 1.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap_err();
        assert_eq!(
            nonfinite,
            SpatialAudioError::NonFinite {
                field: "spatial source sample",
                index: 0,
            }
        );
        assert_eq!(
            render_single_with_gain(
                config(8, 8.0, SpatialDelayPolicy::IntegerCeiling),
                &[1.0],
                SourcePositionTrack::Static([0.0, 0.0, 1.0]),
                ListenerPoseTrack::Static(listener_at([0.0; 3])),
                None,
                f64::INFINITY,
            )
            .unwrap_err(),
            SpatialAudioError::InvalidConfig("finite nonnegative source gain")
        );

        let mut cfg = config(8, 8.0, SpatialDelayPolicy::IntegerCeiling);
        cfg.budget.maximum_abs_output_fs = 0.5;
        let clipped = render_single(
            cfg,
            &[1.0],
            SourcePositionTrack::Static([-1.0, 0.0, 0.0]),
            ListenerPoseTrack::Static(listener_at([0.0; 3])),
            None,
        )
        .unwrap_err();
        assert!(matches!(
            clipped,
            SpatialAudioError::OutputAmplitudeExceeded {
                frame: 1,
                channel: "left",
                ..
            }
        ));
    }

    fn retarded_config() -> SpatialAudioConfig {
        let mut result = config(
            8_192,
            256.0,
            SpatialDelayPolicy::RetardedSinc {
                source_bandwidth_hz: 256.0,
                receiver_bandwidth_hz: 1_024.0,
            },
        );
        // Keep distance gain exactly one for the analytical clock/amplitude
        // checks; those tests are about reconstruction, not source radiation.
        result.minimum_distance_m = 1_024.0;
        result
    }

    /// G1/G3: the exact longitudinal two-motion law is the independent oracle.
    #[test]
    fn retarded_source_and_listener_motion_preserve_pitch_and_sample_amplitude() {
        let count = 2_048;
        let fs = 8_192.0;
        let c = 256.0;
        let frequency = 128.0;
        let separation = 64.0;
        let signal: Vec<f64> = (0..count)
            .map(|i| (core::f64::consts::TAU * frequency * i as f64 / fs).sin())
            .collect();
        for (source_velocity, receiver_velocity, flow_velocity) in [
            (64.0, -32.0, 0.0),
            (-64.0, 32.0, 0.0),
            (192.0, 192.0, 0.0),
            (0.0, -32.0, 0.0),
            (0.0, 0.0, 128.0),
            (0.0, 0.0, -128.0),
            (64.0, -32.0, 32.0),
            // Laboratory speeds exceed c, but both are subsonic in this air.
            (320.0, 352.0, 384.0),
        ] {
            let positions: Vec<[f64; 3]> = (0..count)
                .map(|i| [0.0, 0.0, source_velocity * i as f64 / fs])
                .collect();
            let listeners: Vec<ListenerPose> = (0..16_384)
                .map(|i| listener_at([0.0, 0.0, separation + receiver_velocity * i as f64 / fs]))
                .collect();
            let mut cfg = retarded_config();
            if flow_velocity != 0.0 {
                cfg.delay_policy = SpatialDelayPolicy::RetardedSincUniformFlow {
                    source_bandwidth_hz: 256.0,
                    receiver_bandwidth_hz: 1_024.0,
                    medium_velocity_m_per_s: [0.0, 0.0, flow_velocity],
                };
            }
            let output = render_single(
                cfg,
                &signal,
                SourcePositionTrack::PerFrame(&positions),
                ListenerPoseTrack::PerFrame(&listeners),
                None,
            )
            .unwrap();
            let mut checked = 0;
            for (frame, sample) in output.samples().iter().enumerate() {
                let reception = frame as f64 / fs;
                let downstream_speed = c + flow_velocity;
                let emission = (reception * (downstream_speed - receiver_velocity) - separation)
                    / (downstream_speed - source_velocity);
                if emission * fs < 32.0 || emission * fs > (count - 33) as f64 {
                    continue;
                }
                let expected =
                    (core::f64::consts::TAU * frequency * emission).sin() * det::sqrt(0.5);
                assert!(
                    (sample.left_fs - expected).abs() < 1.0e-4,
                    "vs={source_velocity}, vr={receiver_velocity}, U={flow_velocity}, frame={frame}, got={}, expected={expected}",
                    sample.left_fs
                );
                assert_eq!(sample.left_fs.to_bits(), sample.right_fs.to_bits());
                checked += 1;
            }
            assert!(checked > 1_000, "complete causal test window");
        }
    }

    /// G1: receding sources need their actual compressed received band, not
    /// the expansion a source with the same speed would have on approach.
    #[test]
    fn retarded_directional_admission_renders_receding_sources_and_checks_every_span() {
        let count = 512;
        let signal: Vec<f64> = (0..count)
            .map(|i| (core::f64::consts::TAU * 128.0 * i as f64 / 8_192.0).sin())
            .collect();
        let positions: Vec<[f64; 3]> = (0..count)
            .map(|i| [0.0, 0.0, -192.0 * i as f64 / 8_192.0])
            .collect();
        let cfg = SpatialAudioConfig {
            delay_policy: SpatialDelayPolicy::RetardedSinc {
                source_bandwidth_hz: 1_800.0,
                receiver_bandwidth_hz: 1_200.0,
            },
            ..retarded_config()
        };
        let listener = ListenerPoseTrack::Static(listener_at([0.0, 0.0, 64.0]));
        let output = render_single(
            cfg,
            &signal,
            SourcePositionTrack::PerFrame(&positions),
            listener,
            None,
        )
        .unwrap();
        let mut checked = 0;
        for (frame, sample) in output.samples().iter().enumerate() {
            let emission = (frame as f64 - 2_048.0) / 1.75;
            if !(32.0..=(count - 33) as f64).contains(&emission) {
                continue;
            }
            let expected =
                (core::f64::consts::TAU * 128.0 * emission / 8_192.0).sin() * det::sqrt(0.5);
            assert!((sample.left_fs - expected).abs() < 3.0e-4);
            assert_eq!(sample.left_fs.to_bits(), sample.right_fs.to_bits());
            checked += 1;
        }
        assert!(checked > 700);
        for turn_around in [false, true] {
            let approaching: Vec<[f64; 3]> = (0..count)
                .map(|i| {
                    let travel = if turn_around {
                        (i as f64 - 256.0).abs() - 256.0
                    } else {
                        i as f64
                    };
                    [0.0, 0.0, 192.0 * travel / 8_192.0]
                })
                .collect();
            assert!(
                render_single(
                    cfg,
                    &signal,
                    SourcePositionTrack::PerFrame(&approaching),
                    listener,
                    None,
                )
                .is_err(),
                "approach, including after a receding span, must retain aliasing refusal"
            );
        }
    }

    /// G1: receiver recession is a numerator change, distinct from source
    /// recession. Physical reception poses and the finite tail remain required.
    #[test]
    fn retarded_directional_admission_uses_the_moving_receiver_clock() {
        let count = 512;
        let signal: Vec<f64> = (0..count)
            .map(|i| (core::f64::consts::TAU * 128.0 * i as f64 / 8_192.0).sin())
            .collect();
        let poses: Vec<ListenerPose> = (0..8_192)
            .map(|i| listener_at([0.0, 0.0, 16.0 + 192.0 * f64::from(i) / 8_192.0]))
            .collect();
        let cfg = SpatialAudioConfig {
            delay_policy: SpatialDelayPolicy::RetardedSinc {
                source_bandwidth_hz: 1_800.0,
                receiver_bandwidth_hz: 600.0,
            },
            ..retarded_config()
        };
        let output = render_single(
            cfg,
            &signal,
            SourcePositionTrack::Static([0.0; 3]),
            ListenerPoseTrack::PerFrame(&poses),
            None,
        )
        .unwrap();
        let mut checked = 0;
        for (frame, sample) in output.samples().iter().enumerate() {
            let emission = 0.25 * frame as f64 - 512.0;
            if !(32.0..=(count - 33) as f64).contains(&emission) {
                continue;
            }
            let expected =
                (core::f64::consts::TAU * 128.0 * emission / 8_192.0).sin() * det::sqrt(0.5);
            assert!((sample.left_fs - expected).abs() < 3.0e-4);
            checked += 1;
        }
        assert!(checked > 1_700);
        let approaching: Vec<ListenerPose> = (0..8_192)
            .map(|i| listener_at([0.0, 0.0, 16.0 - 192.0 * f64::from(i) / 8_192.0]))
            .collect();
        assert!(
            render_single(
                cfg,
                &signal,
                SourcePositionTrack::Static([0.0; 3]),
                ListenerPoseTrack::PerFrame(&approaching),
                None,
            )
            .is_err()
        );
    }

    #[test]
    fn retarded_directional_refinement_charges_each_source_scan_only_when_needed() {
        let signal = [0.25; 32];
        let poses: Vec<ListenerPose> = (0..4_096)
            .map(|i| listener_at([0.0, 0.0, 1.0 + 192.0 * f64::from(i) / 8_192.0]))
            .collect();
        let sources: [SpatialAudioSource<'_>; 3] = core::array::from_fn(|i| SpatialAudioSource {
            source_identity: identity("budgeted fixed source"),
            signal: SpatialMonoSignal::Samples(&signal),
            positions: SourcePositionTrack::Static([0.0, 0.0, -(i as f64)]),
            gain_linear: 1.0,
            authority: SpatialAudioAuthority::PhysicallyParameterized,
        });
        let input = SpatialAudioRenderInput {
            sources: &sources,
            listener: ListenerPoseTrack::PerFrame(&poses),
            room_ir: None,
        };
        let mut cfg = SpatialAudioConfig {
            delay_policy: SpatialDelayPolicy::RetardedSinc {
                source_bandwidth_hz: 1_800.0,
                receiver_bandwidth_hz: 3_500.0,
            },
            output_horizon: SpatialOutputHorizon::ClampToInputFrames,
            ..retarded_config()
        };
        with_cx(false, |cx| {
            // Ask the public admission result for the unchanged wide-band work.
            // The long receiver history remains input even for clamped output.
            cfg.budget.maximum_work_units = 1;
            let original_work = match OfflineSpatializer::try_new(cfg, cx)
                .unwrap()
                .spatialize(input, cx)
            {
                Err(SpatialAudioError::ResourceLimit {
                    resource: "work units",
                    requested,
                    ..
                }) => requested,
                other => panic!("expected original work-budget diagnostic, got {other:?}"),
            };
            cfg.budget.maximum_work_units = original_work;
            let wide = OfflineSpatializer::try_new(cfg, cx)
                .unwrap()
                .spatialize(input, cx)
                .unwrap();

            cfg.delay_policy = SpatialDelayPolicy::RetardedSinc {
                source_bandwidth_hz: 1_800.0,
                receiver_bandwidth_hz: 600.0,
            };
            let per_source_scan = (poses.len() - 1) as u64;
            assert!(matches!(
                OfflineSpatializer::try_new(cfg, cx).unwrap().spatialize(input, cx),
                Err(SpatialAudioError::ResourceLimit {
                    resource: "work units", requested, ..
                }) if requested == original_work + per_source_scan
            ));
            let full_work = original_work + sources.len() as u64 * per_source_scan;
            cfg.budget.maximum_work_units = full_work - 1;
            assert!(matches!(
                OfflineSpatializer::try_new(cfg, cx).unwrap().spatialize(input, cx),
                Err(SpatialAudioError::ResourceLimit {
                    resource: "work units", requested, ..
                }) if requested == full_work
            ));
            cfg.budget.maximum_work_units = full_work;
            let narrow = OfflineSpatializer::try_new(cfg, cx)
                .unwrap()
                .spatialize(input, cx)
                .unwrap();
            assert_eq!(narrow.samples(), wide.samples());
        });
    }

    /// G1: fixed laboratory endpoints retain their source frequency. Wind
    /// changes transit time and the acoustic wave normal, including crosswind.
    #[test]
    fn uniform_flow_static_delay_and_direction_follow_the_acoustic_ray() {
        let signal: Vec<f64> = (0..128).map(|i| (f64::from(i) * 0.1).sin()).collect();
        let mut previous_identity = None;
        for (wind, delay) in [(0.0, 768), (128.0, 512), (-128.0, 1_536)] {
            let mut cfg = retarded_config();
            cfg.delay_policy = SpatialDelayPolicy::RetardedSincUniformFlow {
                source_bandwidth_hz: 256.0,
                // Equal source and receiver bands admit stationary endpoints
                // even when the wind is large: their Doppler ratio is one.
                receiver_bandwidth_hz: 256.0,
                medium_velocity_m_per_s: [0.0, 0.0, wind],
            };
            let output = render_single(
                cfg,
                &signal,
                SourcePositionTrack::Static([0.0; 3]),
                ListenerPoseTrack::Static(listener_at([0.0, 0.0, 24.0])),
                None,
            )
            .unwrap();
            assert_eq!(output.samples().len(), signal.len() + delay);
            for (i, value) in signal.iter().copied().enumerate() {
                assert_eq!(
                    output.samples()[i + delay].left_fs.to_bits(),
                    (value * det::sqrt(0.5)).to_bits()
                );
            }
            assert_eq!(output.diagnostics().maximum_delay_frames, delay as f64);
            if let Some(previous) = previous_identity {
                assert_ne!(output.config_identity(), previous);
            }
            previous_identity = Some(output.config_identity());
        }

        let mut cfg = retarded_config();
        cfg.delay_policy = SpatialDelayPolicy::RetardedSincUniformFlow {
            source_bandwidth_hz: 0.0,
            receiver_bandwidth_hz: 0.0,
            medium_velocity_m_per_s: [192.0, 0.0, 0.0],
        };
        cfg.microphone_directivity = MicrophoneDirectivity::Cardioid {
            rear_floor_gain: 0.0,
        };
        let listener = ListenerPose {
            forward_unit: [0.0, 0.0, -1.0],
            ..listener_at([0.0, 0.0, 2.0])
        };
        let output = render_single(
            cfg,
            &[1.0; 256],
            SourcePositionTrack::Static([0.0; 3]),
            ListenerPoseTrack::Static(listener),
            None,
        )
        .unwrap();
        // n = (-M, 0, sqrt(1-M^2)), not the laboratory +z direction.
        let longitudinal = (1.0_f64 - 0.75_f64.powi(2)).sqrt();
        let delay = 64.0 / longitudinal;
        let sample = output.samples()[delay.ceil() as usize + 64];
        let cardioid = 0.5 * (1.0 + longitudinal);
        assert!((sample.left_fs - cardioid * 0.125_f64.sqrt()).abs() < 1.0e-12);
        assert!((sample.right_fs - cardioid * 0.875_f64.sqrt()).abs() < 1.0e-12);
        assert!((output.diagnostics().maximum_distance_m - 2.0 / longitudinal).abs() < 1.0e-12);
    }

    /// G3: a common Galilean boost changes neither clock nor observation ray.
    #[test]
    fn uniform_flow_common_boost_preserves_the_complete_waveform() {
        let signal: Vec<f64> = (0..2_048).map(|i| (f64::from(i) * 0.0625).sin()).collect();
        let source: Vec<[f64; 3]> = (0..signal.len())
            .map(|i| [-4.0, 1.0, 64.0 * i as f64 / 8_192.0])
            .collect();
        let receiver: Vec<ListenerPose> = (0..16_384)
            .map(|i| listener_at([1.0, 2.0, 64.0 - 32.0 * f64::from(i) / 8_192.0]))
            .collect();
        let original = render_single(
            retarded_config(),
            &signal,
            SourcePositionTrack::PerFrame(&source),
            ListenerPoseTrack::PerFrame(&receiver),
            None,
        )
        .unwrap();
        let boost = [512.0, -256.0, 128.0];
        let boosted_source: Vec<[f64; 3]> = source
            .iter()
            .enumerate()
            .map(|(i, p)| core::array::from_fn(|j| p[j] + boost[j] * i as f64 / 8_192.0))
            .collect();
        let boosted_receiver: Vec<ListenerPose> = receiver
            .iter()
            .enumerate()
            .map(|(i, p)| ListenerPose {
                position_m: core::array::from_fn(|j| {
                    p.position_m[j] + boost[j] * i as f64 / 8_192.0
                }),
                ..*p
            })
            .collect();
        let mut cfg = retarded_config();
        cfg.delay_policy = SpatialDelayPolicy::RetardedSincUniformFlow {
            source_bandwidth_hz: 256.0,
            receiver_bandwidth_hz: 1_024.0,
            medium_velocity_m_per_s: boost,
        };
        let boosted = render_single(
            cfg,
            &signal,
            SourcePositionTrack::PerFrame(&boosted_source),
            ListenerPoseTrack::PerFrame(&boosted_receiver),
            None,
        )
        .unwrap();
        assert_eq!(boosted.samples(), original.samples());
        assert_eq!(
            boosted.diagnostics().maximum_delay_frames,
            original.diagnostics().maximum_delay_frames
        );
    }

    #[test]
    fn uniform_flow_admission_uses_medium_relative_mach_and_bandwidth() {
        let signal = [0.0; 64];
        let mut cfg = retarded_config();
        let moving: Vec<[f64; 3]> = (0..signal.len())
            .map(|i| [0.0, 0.0, 192.0 * i as f64 / 8_192.0])
            .collect();
        let listener = ListenerPoseTrack::Static(listener_at([0.0, 0.0, 1.0]));
        for wind in [-128.0, -32.0, 256.0] {
            cfg.delay_policy = SpatialDelayPolicy::RetardedSincUniformFlow {
                source_bandwidth_hz: 256.0,
                receiver_bandwidth_hz: 1_024.0,
                medium_velocity_m_per_s: [0.0, 0.0, wind],
            };
            // -128: source is supersonic relative to air. -32: the admitted
            // receiver band cannot cover compression. +256: listener is sonic.
            assert!(
                render_single(
                    cfg,
                    &signal,
                    SourcePositionTrack::PerFrame(&moving),
                    listener,
                    None
                )
                .is_err()
            );
        }
        cfg.delay_policy = SpatialDelayPolicy::RetardedSincUniformFlow {
            source_bandwidth_hz: 256.0,
            receiver_bandwidth_hz: 1_024.0,
            medium_velocity_m_per_s: [f64::NAN, 0.0, 0.0],
        };
        with_cx(false, |cx| {
            assert!(OfflineSpatializer::try_new(cfg, cx).is_err());
        });
    }

    #[test]
    fn retarded_static_integer_delay_is_exact_and_tail_remains_explicit() {
        let samples: Vec<f64> = (0..128).map(|i| (f64::from(i) * 0.1).sin()).collect();
        let cfg = retarded_config();
        // c/fs = 1/32 m per frame, so this is an exact 16-frame delay.
        let positions = SourcePositionTrack::Static([0.0, 0.0, 0.5]);
        let listener = ListenerPoseTrack::Static(listener_at([0.0; 3]));
        let output = render_single(cfg, &samples, positions, listener, None).unwrap();
        assert_eq!(output.samples().len(), samples.len() + 16);
        for (i, value) in samples.iter().copied().enumerate() {
            assert_eq!(
                output.samples()[i + 16].left_fs.to_bits(),
                (value * det::sqrt(0.5)).to_bits()
            );
        }
        let clamped = render_single(
            SpatialAudioConfig {
                output_horizon: SpatialOutputHorizon::ClampToInputFrames,
                ..cfg
            },
            &samples,
            positions,
            listener,
            None,
        )
        .unwrap();
        assert_eq!(clamped.samples(), &output.samples()[..samples.len()]);
        assert_eq!(clamped.diagnostics().discarded_tail_frames, 16);
        assert_ne!(clamped.config_identity(), output.config_identity());
    }

    #[test]
    fn retarded_microphone_axes_are_read_at_reception_and_room_ir_keeps_its_clock() {
        let signal = [1.0; 64];
        let listeners: Vec<ListenerPose> = (0..128)
            .map(|i| {
                let mut pose = listener_at([0.0; 3]);
                if i >= 32 {
                    pose.forward_unit = [0.0, 0.0, -1.0];
                }
                pose
            })
            .collect();
        let mut cfg = retarded_config();
        cfg.microphone_directivity = MicrophoneDirectivity::Cardioid {
            rear_floor_gain: 0.0,
        };
        let positions = SourcePositionTrack::Static([0.0, 0.0, 0.5]);
        let listener = ListenerPoseTrack::PerFrame(&listeners);
        let dry = render_single(cfg, &signal, positions, listener, None).unwrap();
        assert!(dry.samples()[31].left_fs > 0.7);
        assert_eq!(dry.samples()[32].left_fs, 0.0);
        with_cx(false, |cx| {
            let room = StereoRoomImpulseResponse::try_new(
                8_192,
                vec![1.0, 0.25],
                vec![1.0, -0.25],
                SpatialAudioAuthority::Artistic,
                cx,
            )
            .unwrap();
            let wet = render_single(cfg, &signal, positions, listener, Some(&room)).unwrap();
            assert_eq!(wet.samples()[32].left_fs, 0.25 * dry.samples()[31].left_fs);
            assert_eq!(
                wet.samples()[32].right_fs,
                -0.25 * dry.samples()[31].right_fs
            );
            assert_eq!(wet.samples().len(), dry.samples().len() + 1);
        });
    }

    #[test]
    fn retarded_admission_refuses_aliasing_sonic_motion_and_missing_receiver_tail() {
        let signal = [0.0; 64];
        let positions: Vec<[f64; 3]> = (0..64)
            .map(|i| [0.0, 0.0, f64::from(i) * 0.03125])
            .collect();
        let listener = ListenerPoseTrack::Static(listener_at([0.0, 0.0, 64.0]));
        assert!(
            render_single(
                retarded_config(),
                &signal,
                SourcePositionTrack::PerFrame(&positions),
                listener,
                None
            )
            .is_err()
        );
        let approaching: Vec<[f64; 3]> = positions.iter().map(|p| p.map(|v| 0.75 * v)).collect();
        let narrow = SpatialAudioConfig {
            delay_policy: SpatialDelayPolicy::RetardedSinc {
                source_bandwidth_hz: 1_800.0,
                receiver_bandwidth_hz: 3_500.0,
            },
            ..retarded_config()
        };
        assert!(
            render_single(
                narrow,
                &signal,
                SourcePositionTrack::PerFrame(&approaching),
                listener,
                None
            )
            .is_err()
        );
        let short = [listener_at([0.0, 0.0, 1.0]); 64];
        assert!(
            render_single(
                retarded_config(),
                &signal,
                SourcePositionTrack::Static([0.0; 3]),
                ListenerPoseTrack::PerFrame(&short),
                None
            )
            .is_err()
        );
        with_cx(false, |cx| {
            let unsupported_band = SpatialAudioConfig {
                delay_policy: SpatialDelayPolicy::RetardedSinc {
                    source_bandwidth_hz: 3_000.0,
                    receiver_bandwidth_hz: 3_500.0,
                },
                ..retarded_config()
            };
            assert!(OfflineSpatializer::try_new(unsupported_band, cx).is_err());
        });
    }

    #[test]
    fn retarded_work_budget_cancellation_and_replay_use_the_original_transaction() {
        let signal = [0.25; 128];
        let source = SpatialAudioSource {
            source_identity: identity("retarded source"),
            signal: SpatialMonoSignal::Samples(&signal),
            positions: SourcePositionTrack::Static([0.0, 0.0, 0.5]),
            gain_linear: 1.0,
            authority: SpatialAudioAuthority::PhysicallyParameterized,
        };
        let input = SpatialAudioRenderInput {
            sources: &[source],
            listener: ListenerPoseTrack::Static(listener_at([0.0; 3])),
            room_ir: None,
        };
        with_cx(false, |cx| {
            let spatializer = OfflineSpatializer::try_new(retarded_config(), cx).unwrap();
            let expected = spatializer.spatialize(input, cx).unwrap();
            let mut calls = 0;
            assert!(matches!(
                spatialize_with_checkpoint(spatializer, input, &mut || {
                    calls += 1;
                    if calls == 8 {
                        Err(SpatialAudioError::Cancelled)
                    } else {
                        Ok(())
                    }
                }),
                Err(SpatialAudioError::Cancelled)
            ));
            assert_eq!(spatializer.spatialize(input, cx).unwrap(), expected);
            let mut limited = retarded_config();
            limited.budget.maximum_work_units = 2_000;
            assert!(matches!(
                OfflineSpatializer::try_new(limited, cx)
                    .unwrap()
                    .spatialize(input, cx),
                Err(SpatialAudioError::ResourceLimit {
                    resource: "work units",
                    ..
                })
            ));
        });
    }

    #[test]
    fn cancellation_before_and_during_work_publishes_no_output() {
        let cfg = config(48_000, 343.0, SpatialDelayPolicy::LinearFloorCeil);
        with_cx(true, |cx| {
            assert_eq!(
                OfflineSpatializer::try_new(cfg, cx),
                Err(SpatialAudioError::Cancelled)
            );
        });

        with_cx(false, |cx| {
            let spatializer = OfflineSpatializer::try_new(cfg, cx).unwrap();
            let samples = vec![0.1; 2_048];
            let source = SpatialAudioSource {
                source_identity: identity("cancel-source"),
                signal: SpatialMonoSignal::Samples(&samples),
                positions: SourcePositionTrack::Static([0.0, 0.0, 1.0]),
                gain_linear: 1.0,
                authority: SpatialAudioAuthority::PhysicallyParameterized,
            };
            let input = SpatialAudioRenderInput {
                sources: &[source],
                listener: ListenerPoseTrack::Static(listener_at([0.0; 3])),
                room_ir: None,
            };
            let mut polls = 0_usize;
            let result = spatialize_with_checkpoint(spatializer, input, &mut || {
                polls += 1;
                if polls == 5 {
                    Err(SpatialAudioError::Cancelled)
                } else {
                    Ok(())
                }
            });
            assert_eq!(result, Err(SpatialAudioError::Cancelled));
            assert_eq!(polls, 5);
        });
    }
}
