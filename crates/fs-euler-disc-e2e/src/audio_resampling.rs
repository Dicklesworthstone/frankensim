//! Deterministic, measure-first reconstruction of Euler-disc modal drive at 48 kHz.
//!
//! Source intervals are finite-volume cells: their retained force-time measures,
//! rather than point-sampled means, are conservatively rasterized onto the exact
//! audio clock. Localized contact/rolling drive is projected through each mode's
//! time-varying contact factor before anti-alias filtering; distributed base/gas
//! drive remains in component coordinates. A centered Blackman-Harris FIR uses
//! global-horizon even reflection and explicit offline group-delay compensation,
//! so output chunk boundaries never reset filter state.
//!
//! This is signal reconstruction for a model-derived soundtrack. A cadence-based
//! Nyquist ceiling is only an admission limit; it is not evidence that mechanics
//! channels contain all physical acoustic bandwidth, nor is the output calibrated
//! pressure or a validated radiation prediction.

use core::{f64::consts::PI, fmt, num::NonZeroUsize};

use fs_blake3::{ContentHash, DomainHasher};
use fs_evidence::{
    cinematic::{CinematicClock, CinematicClockDomain},
    cinematic_sound::{
        SOUND_MASTER_SAMPLE_RATE_HZ, SOUND_MASTER_VIDEO_RATE_HZ, SoundMode, SoundModeParticipation,
        SoundSynthesisConfig,
    },
};
use fs_exec::Cx;
use fs_math::{STRICT_CORE_GOLDEN_HASH, STRICT_CORE_SEMANTICS_VERSION, det};

use crate::audio_excitation::{AudioExcitationEvent, AudioExcitationReconstructionStatus};
use crate::control_stream::ContactEventMeasure;
use crate::coupled_runner::ContactTransitionKind;
use crate::{
    AUDIO_EXCITATION_ALGORITHM_VERSION, AudioExcitationGrid, AudioExcitationInterval,
    AudioExcitationMapper, MAX_MODAL_SPATIAL_PARTICIPATION, MODAL_SYNTHESIS_ALGORITHM_VERSION,
    ModalComponentValues, ModalDriveFrame, ModalSpatialParticipation, ModalSynthesisCheckpoint,
    ModalSynthesisChunk, ModalSynthesisError, ModalSynthesisModel, procedural_texture_unit_sample,
};

/// Version of the interval-to-audio reconstruction, checkpoint, and identity semantics.
pub const AUDIO_RESAMPLING_ALGORITHM_VERSION: u32 = 1;
/// Version of the exact Blackman-Harris coefficient and response-audit semantics.
pub const AUDIO_RECONSTRUCTION_FILTER_VERSION: u32 = 1;
/// Maximum supported odd FIR length.
pub const MAX_AUDIO_RECONSTRUCTION_FILTER_TAPS: usize = 4_097;
/// Maximum audio frames between cancellation polls in raster/convolution loops.
pub const AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES: usize = 64;
/// Maximum canonical modes traversed between cancellation polls in inner loops.
pub const AUDIO_RESAMPLING_CANCELLATION_POLL_MODES: usize = 8;
/// Source-clock endpoint mismatch admitted by the binary64 trajectory boundary [frames].
pub const MAX_SOURCE_CLOCK_ALIGNMENT_ERROR_FRAMES: f64 = 1.0e-6;
/// Near-integer event coordinate tolerance used only to remove binary64 roundoff [frames].
pub const EVENT_SAMPLE_SNAP_TOLERANCE_FRAMES: f64 = 1.0e-9;
/// Weakest accepted stopband contract for the production reconstruction path.
pub const MIN_AUDIO_FILTER_STOPBAND_ATTENUATION_DB: f64 = 80.0;
/// Largest accepted passband-ripple contract for the production path.
pub const MAX_AUDIO_FILTER_PASSBAND_RIPPLE_DB: f64 = 0.1;

const RESAMPLER_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-resampler.v1";
const FILTER_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-cinematic.audio-reconstruction-filter.v1";
const SOURCE_PAYLOAD_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-cinematic.audio-excitation-payload.v1";
const CHUNK_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-drive-chunk.v1";
const LN_10: f64 = core::f64::consts::LN_10;

/// Boundary extension used by the centered FIR.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioResamplingBoundaryPolicy {
    /// Reflect around each half-sample horizon edge. The global horizon, never
    /// an output chunk, owns the reflection boundary.
    HalfSampleEvenReflectionV1,
}

/// Fractional event placement rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioEventFractionalDelay {
    /// Split an impulse between its neighboring sample boundaries. This
    /// preserves impulse and first moment exactly for interior binary64 times.
    LinearTwoBoundaryV1,
}

/// Explicit deterministic low-pass design and acceptance thresholds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioReconstructionFilterSpec {
    /// Highest declared passband frequency [Hz].
    pub passband_edge_hz: f64,
    /// First declared stopband frequency [Hz].
    pub stopband_edge_hz: f64,
    /// Number of taps on either side of the centered coefficient.
    pub half_length: u16,
    /// Maximum measured passband magnitude ripple [dB].
    pub maximum_passband_ripple_db: f64,
    /// Minimum measured stopband attenuation [dB].
    pub minimum_stopband_attenuation_db: f64,
    /// Number of deterministic intervals in the admission response sweep.
    pub response_grid_intervals: u32,
}

/// Measured response and latency of the admitted physical-control filter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioReconstructionFilterDiagnostics {
    /// Exact odd tap count.
    pub tap_count: usize,
    /// Measured worst passband ripple [dB].
    pub measured_passband_ripple_db: f64,
    /// Measured worst stopband attenuation [dB].
    pub measured_stopband_attenuation_db: f64,
    /// Causal realization delay before compensation [frames].
    pub intrinsic_group_delay_frames: u32,
    /// Offline compensation applied by centered evaluation [frames].
    pub group_delay_compensation_frames: u32,
    /// Future frames needed to evaluate one published frame.
    pub required_lookahead_frames: u32,
    /// Published sample offset after compensation; always zero in v1.
    pub published_alignment_offset_frames: i32,
}

/// Explicit memory/work ceilings. Violations refuse before result allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AudioResamplingBudget {
    /// Maximum complete source intervals bound into one model.
    pub maximum_source_intervals: usize,
    /// Maximum exact master-clock audio frames.
    pub maximum_total_audio_frames: u64,
    /// Maximum frames in one atomically published output chunk.
    pub maximum_chunk_audio_frames: usize,
    /// Maximum admitted FIR taps.
    pub maximum_filter_taps: usize,
    /// Maximum values in each row-major output or FIR-halo mode array.
    pub maximum_chunk_mode_values: usize,
    /// Maximum retained timing/event records.
    pub maximum_events: usize,
    /// Maximum exact video-boundary synchronization markers.
    pub maximum_sync_markers: usize,
    /// Maximum multiply-adds estimated for one chunk, including its halo.
    pub maximum_chunk_multiply_adds: u64,
}

impl AudioResamplingBudget {
    /// Bounded reference-film budget: up to 12 seconds at 48 kHz, 256 modes,
    /// 4,097 FIR taps, and 65,536-frame transactional chunks.
    #[must_use]
    pub const fn reference_film() -> Self {
        Self {
            maximum_source_intervals: 1_048_576,
            maximum_total_audio_frames: 12 * SOUND_MASTER_SAMPLE_RATE_HZ as u64,
            maximum_chunk_audio_frames: 65_536,
            maximum_filter_taps: MAX_AUDIO_RECONSTRUCTION_FILTER_TAPS,
            maximum_chunk_mode_values: (65_536 + MAX_AUDIO_RECONSTRUCTION_FILTER_TAPS - 1) * 256,
            maximum_events: 65_536,
            maximum_sync_markers: 8_193,
            maximum_chunk_multiply_adds: 100_000_000_000,
        }
    }
}

/// Complete resampling input. Every field is bound into model identity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioResamplingModelInput {
    /// Exact reference-video clock.
    pub video_clock: CinematicClock,
    /// Exact reference-audio clock.
    pub audio_clock: CinematicClock,
    /// Declared highest mechanics-derived source frequency [Hz].
    pub declared_source_bandwidth_hz: f64,
    /// Physical continuous-control reconstruction filter.
    pub filter: AudioReconstructionFilterSpec,
    /// Global-horizon padding rule.
    pub boundary_policy: AudioResamplingBoundaryPolicy,
    /// Interior event-placement rule.
    pub event_fractional_delay: AudioEventFractionalDelay,
    /// Explicit resource ceilings.
    pub budget: AudioResamplingBudget,
}

/// Exact video-frame boundary to audio-sample boundary correspondence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AudioVideoSyncMarker {
    /// Absolute video clock tick.
    pub video_tick: i64,
    /// Absolute aligned audio clock tick.
    pub audio_tick: i64,
    /// Zero-based audio frame offset from the admitted audio-clock start.
    pub audio_frame_offset: u64,
}

/// Exact master-clock alignment, including the exclusive endpoint marker.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AudioVideoAlignment {
    /// Exact number of audio frames per video frame for the frozen master.
    pub audio_frames_per_video_frame: u32,
    /// Ordered frame-boundary markers, including the exclusive endpoint.
    pub markers: Vec<AudioVideoSyncMarker>,
    /// Accumulated drift by integer construction; always zero in v1.
    pub endpoint_drift_audio_frames: i64,
}

/// Placement receipt for one timing-only source event.
#[derive(Debug, Clone, PartialEq)]
pub struct ResampledAudioEvent {
    /// Original source event.
    pub source: AudioExcitationEvent,
    /// Exact direct binary64 conversion to the master sample coordinate.
    pub requested_sample_position: f64,
    /// Left emitted frame, if the event lies before the exclusive endpoint.
    pub left_frame_offset: Option<u64>,
    /// Right emitted frame for a non-integral interior event.
    pub right_frame_offset: Option<u64>,
    /// Fraction of an artistic impulse placed on the left boundary.
    pub left_weight: f64,
    /// Fraction placed on the right boundary.
    pub right_weight: f64,
    /// Emitted impulse centroid minus requested position [frames].
    pub centroid_error_frames: f64,
    /// Source localization bracket start in audio-frame coordinates.
    pub bracket_start_sample_position: f64,
    /// Source localization bracket end in audio-frame coordinates.
    pub bracket_end_sample_position: f64,
}

/// Immutable restart point. Absolute indexing reproduces all FIR halo and RNG state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AudioResamplingCheckpoint {
    model_identity: ContentHash,
    next_audio_frame_offset: u64,
}

impl AudioResamplingCheckpoint {
    /// Bound resampling-model identity.
    #[must_use]
    pub const fn model_identity(self) -> ContentHash {
        self.model_identity
    }

    /// First audio frame not yet published.
    #[must_use]
    pub const fn next_audio_frame_offset(self) -> u64 {
        self.next_audio_frame_offset
    }
}

/// Chunk-local diagnostics suitable for deterministic run logging.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioResamplingDiagnostics {
    /// Inclusive zero-based first output frame.
    pub start_audio_frame_offset: u64,
    /// Exclusive zero-based successor frame.
    pub end_audio_frame_offset: u64,
    /// Exact filter halo used on either side [frames].
    pub filter_half_length_frames: u32,
    /// Largest absolute distributed component force in the chunk [N].
    pub maximum_abs_distributed_force_n: f64,
    /// Largest absolute already-participated localized modal force [N].
    pub maximum_abs_localized_mode_force_n: f64,
    /// Number of uniquely owned event receipts in this chunk.
    pub owned_event_count: usize,
}

/// Atomically published 48 kHz drive and successor checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioResamplingChunk {
    modal_identity: ContentHash,
    start_audio_frame_offset: u64,
    /// Domain-separated identity of model plus exact output range.
    pub identity: ContentHash,
    /// Distributed component drive. Localized component fields are exactly zero
    /// because the row-major arrays below already contain spatial participation.
    pub drive_frames: Vec<ModalDriveFrame>,
    /// Row-major `(frame, canonical_mode)` localized modal force [N].
    pub preparticipated_localized_force_n: Vec<f64>,
    /// Row-major localized modal boundary impulse [N s].
    pub preparticipated_localized_impulse_n_s: Vec<f64>,
    /// Event receipts uniquely owned by this chunk.
    pub events: Vec<ResampledAudioEvent>,
    /// Video boundaries falling in this half-open chunk, plus the final endpoint
    /// when this is the terminal chunk.
    pub sync_markers: Vec<AudioVideoSyncMarker>,
    /// Chunk-local diagnostics.
    pub diagnostics: AudioResamplingDiagnostics,
    /// Immutable restart point for the following chunk.
    pub successor: AudioResamplingCheckpoint,
}

impl AudioResamplingChunk {
    /// Drive the exact modal model and successor position bound by this chunk,
    /// using the only valid localized-drive representation. This prevents both
    /// out-of-order synthesis and accidentally selecting `Declared`
    /// participation, which would silently discard contact/rolling drive.
    pub fn synthesize_modal(
        &self,
        modal: &ModalSynthesisModel,
        checkpoint: &ModalSynthesisCheckpoint,
        cx: &Cx<'_>,
    ) -> Result<ModalSynthesisChunk, ModalSynthesisError> {
        if modal.identity() != self.modal_identity {
            return Err(ModalSynthesisError::SoundConfigurationMismatch(
                "resampled drive modal identity",
            ));
        }
        if checkpoint.next_sample_frame() != self.start_audio_frame_offset {
            return Err(ModalSynthesisError::InvalidCheckpoint);
        }
        modal.synthesize_chunk(
            checkpoint,
            &self.drive_frames,
            ModalSpatialParticipation::PreparticipatedLocalizedDrive {
                generalized_force_n: &self.preparticipated_localized_force_n,
                boundary_impulse_n_s: &self.preparticipated_localized_impulse_n_s,
            },
            cx,
        )
    }
}

/// Typed admission or transactional reconstruction refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum AudioResamplingError {
    /// Cancellation was observed before atomic publication.
    Cancelled,
    /// Source intervals are empty, incomplete, noncontiguous, or malformed.
    InvalidSourceTimeline {
        interval: usize,
        field: &'static str,
    },
    /// The excitation mapper was admitted against a different modal model.
    ExcitationModalMismatch,
    /// Audio/video clock domain, rate, range, or exact correspondence is invalid.
    InvalidMasterClock(&'static str),
    /// Source endpoints disagree with the exact audio clock beyond tolerance.
    SourceClockMismatch {
        endpoint: &'static str,
        error_frames: f64,
    },
    /// Filter parameters or measured response do not meet their declaration.
    InvalidFilter(&'static str),
    /// Declared source content exceeds the conservative source-cadence ceiling.
    UnsupportedSourceBandwidth {
        requested_hz: f64,
        nominal_ceiling_hz: f64,
    },
    /// One named budget is zero or above a hard ceiling.
    InvalidBudget(&'static str),
    /// A requested resource exceeds its explicit budget.
    BudgetExceeded {
        artifact: &'static str,
        requested: u64,
        limit: u64,
    },
    /// Allocation refused after bounded preflight.
    Capacity {
        artifact: &'static str,
        requested: usize,
    },
    /// Event timing or physical/artistic authority separation is malformed.
    InvalidEvent { event: usize, field: &'static str },
    /// A nonzero event impulse cannot be represented before the exclusive endpoint.
    EventOutsideRepresentableRange { event: usize },
    /// A checkpoint belongs to another model or is outside this timeline.
    InvalidCheckpoint,
    /// Every exact audio frame has been published.
    Complete,
    /// The admitted high-level sound configuration disagrees with this model.
    SoundConfigurationMismatch(&'static str),
    /// A deterministic arithmetic result became non-finite.
    NonFiniteResult { frame: u64, field: &'static str },
}

impl fmt::Display for AudioResamplingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for AudioResamplingError {}

#[derive(Debug, Clone)]
struct FirKernel {
    coefficients: Vec<f64>,
    identity: ContentHash,
    diagnostics: AudioReconstructionFilterDiagnostics,
}

#[derive(Debug, Clone)]
struct PlannedEvent {
    receipt: ResampledAudioEvent,
    owner_frame_offset: u64,
    localized_mode_impulse_n_s: Vec<f64>,
}

#[derive(Debug, Clone)]
struct RawAudioSpan {
    first_frame_offset: u64,
    mode_count: usize,
    distributed_force_n: Vec<ModalComponentValues>,
    physical_localized_mode_force_n: Vec<f64>,
    artistic_localized_mode_force_n: Vec<f64>,
}

/// Admitted offline reconstruction model. Complete source intervals are retained
/// so every output chunk can recompute the same global FIR halo without mutable
/// filter history or chunk-boundary artifacts.
pub struct AudioResampler {
    identity: ContentHash,
    filter: FirKernel,
    artistic_filter: Option<FirKernel>,
    source_payload_identity: ContentHash,
    excitation_identity: ContentHash,
    modal_identity: ContentHash,
    modes: Vec<SoundMode>,
    grid: AudioExcitationGrid,
    intervals: Vec<AudioExcitationInterval>,
    events: Vec<PlannedEvent>,
    input: AudioResamplingModelInput,
    alignment: AudioVideoAlignment,
    total_audio_frames: u64,
    audio_frame_period_s: f64,
    source_start_offset_audio_frames: f64,
}

impl AudioResampler {
    /// Admit source payload, clocks, bandwidth, exact coefficients, identities,
    /// and budgets. No waveform or modal state is produced by construction.
    pub fn try_new(
        mapper: &AudioExcitationMapper<'_, '_>,
        modal: &ModalSynthesisModel,
        intervals: Vec<AudioExcitationInterval>,
        input: AudioResamplingModelInput,
        cx: &Cx<'_>,
    ) -> Result<Self, AudioResamplingError> {
        Self::try_new_with_checkpoint(mapper, modal, intervals, input, &mut || checkpoint(cx))
    }

    fn try_new_with_checkpoint(
        mapper: &AudioExcitationMapper<'_, '_>,
        modal: &ModalSynthesisModel,
        intervals: Vec<AudioExcitationInterval>,
        input: AudioResamplingModelInput,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<Self, AudioResamplingError> {
        checkpoint_fn()?;
        validate_budget(input.budget)?;
        if mapper.modal_identity() != modal.identity() {
            return Err(AudioResamplingError::ExcitationModalMismatch);
        }
        let (total_audio_frames, alignment, audio_frame_period_s) = validate_clocks(
            input.video_clock,
            input.audio_clock,
            input.budget,
            checkpoint_fn,
        )?;
        validate_source_timeline(
            mapper.grid(),
            modal.modes(),
            &intervals,
            input.budget,
            checkpoint_fn,
        )?;
        let source_start_offset_audio_frames =
            validate_source_clock(mapper.grid(), input.audio_clock)?;
        if !input.declared_source_bandwidth_hz.is_finite()
            || input.declared_source_bandwidth_hz <= 0.0
        {
            return Err(AudioResamplingError::InvalidFilter(
                "declared source bandwidth",
            ));
        }
        if input.declared_source_bandwidth_hz > mapper.grid().nominal_source_nyquist_ceiling_hz {
            return Err(AudioResamplingError::UnsupportedSourceBandwidth {
                requested_hz: input.declared_source_bandwidth_hz,
                nominal_ceiling_hz: mapper.grid().nominal_source_nyquist_ceiling_hz,
            });
        }
        if input.declared_source_bandwidth_hz > input.filter.passband_edge_hz {
            return Err(AudioResamplingError::InvalidFilter(
                "passband does not contain declared source bandwidth",
            ));
        }
        let filter = design_physical_filter(
            input.filter,
            mapper.grid().nominal_source_nyquist_ceiling_hz,
            input.budget,
            checkpoint_fn,
        )?;
        let texture_band = common_texture_band(&intervals, checkpoint_fn)?;
        let artistic_filter = match texture_band {
            None => None,
            Some((low_hz, high_hz)) => Some(design_artistic_filter(
                low_hz,
                high_hz,
                input.filter.half_length,
                input.budget,
                checkpoint_fn,
            )?),
        };
        let source_payload_identity = source_payload_identity(&intervals, checkpoint_fn)?;
        let events = plan_events(
            &intervals,
            modal.modes(),
            mapper.grid(),
            source_start_offset_audio_frames,
            total_audio_frames,
            checkpoint_fn,
        )?;
        if events.len() > input.budget.maximum_events {
            return Err(AudioResamplingError::BudgetExceeded {
                artifact: "resampled events",
                requested: events.len() as u64,
                limit: input.budget.maximum_events as u64,
            });
        }
        let identity = resampler_identity(
            mapper.identity(),
            modal.identity(),
            source_payload_identity,
            filter.identity,
            artistic_filter.as_ref().map(|kernel| kernel.identity),
            input,
            total_audio_frames,
        );
        let mut modes = Vec::new();
        reserve_exact(&mut modes, modal.modes().len(), "canonical sound modes")?;
        modes.extend_from_slice(modal.modes());
        checkpoint_fn()?;
        Ok(Self {
            identity,
            filter,
            artistic_filter,
            source_payload_identity,
            excitation_identity: mapper.identity(),
            modal_identity: modal.identity(),
            modes,
            grid: mapper.grid(),
            intervals,
            events,
            input,
            alignment,
            total_audio_frames,
            audio_frame_period_s,
            source_start_offset_audio_frames,
        })
    }

    /// Complete reconstruction-model identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact physical-control filter identity required by the sound config.
    #[must_use]
    pub const fn filter_identity(&self) -> ContentHash {
        self.filter.identity
    }

    /// Exact bound source-payload identity.
    #[must_use]
    pub const fn source_payload_identity(&self) -> ContentHash {
        self.source_payload_identity
    }

    /// Measured filter response and explicit compensated latency.
    #[must_use]
    pub const fn filter_diagnostics(&self) -> AudioReconstructionFilterDiagnostics {
        self.filter.diagnostics
    }

    /// Exact Blackman-Harris coefficients in centered order.
    #[must_use]
    pub fn filter_coefficients(&self) -> &[f64] {
        &self.filter.coefficients
    }

    /// Exact A/V boundary map, including the exclusive endpoint.
    #[must_use]
    pub const fn alignment(&self) -> &AudioVideoAlignment {
        &self.alignment
    }

    /// Exact number of output audio frames.
    #[must_use]
    pub const fn total_audio_frames(&self) -> u64 {
        self.total_audio_frames
    }

    /// Selected mechanics grid whose cadence bounds source admission.
    #[must_use]
    pub const fn source_grid(&self) -> AudioExcitationGrid {
        self.grid
    }

    /// Single fixed modal-drive frame period [s].
    #[must_use]
    pub const fn audio_frame_period_s(&self) -> f64 {
        self.audio_frame_period_s
    }

    /// Immutable zero-progress restart point.
    pub fn initial_checkpoint(
        &self,
        cx: &Cx<'_>,
    ) -> Result<AudioResamplingCheckpoint, AudioResamplingError> {
        checkpoint(cx)?;
        Ok(AudioResamplingCheckpoint {
            model_identity: self.identity,
            next_audio_frame_offset: 0,
        })
    }

    /// Validate the admitted high-level sound config against exact clocks and
    /// algorithm identities. Spatialization/room fields remain later-stage inputs.
    pub fn validate_sound_configuration(
        &self,
        sound: &SoundSynthesisConfig,
    ) -> Result<(), AudioResamplingError> {
        let sound_input = sound.input();
        if sound_input.excitation.identity() != self.excitation_identity
            || sound_input.excitation.version() != AUDIO_EXCITATION_ALGORITHM_VERSION
        {
            return Err(AudioResamplingError::SoundConfigurationMismatch(
                "excitation identity or version",
            ));
        }
        if sound_input.sound_model.identity() != self.modal_identity
            || sound_input.sound_model.version() != MODAL_SYNTHESIS_ALGORITHM_VERSION
            || sound_input.modes.as_slice() != self.modes.as_slice()
        {
            return Err(AudioResamplingError::SoundConfigurationMismatch(
                "modal identity, version, or modes",
            ));
        }
        if sound_input.resampler_identity != self.identity
            || sound_input.resampler_version != AUDIO_RESAMPLING_ALGORITHM_VERSION
            || sound_input.filter_identity != self.filter.identity
            || sound_input.filter_version != AUDIO_RECONSTRUCTION_FILTER_VERSION
        {
            return Err(AudioResamplingError::SoundConfigurationMismatch(
                "resampler or filter identity/version",
            ));
        }
        if sound_input.audio_clock != self.input.audio_clock
            || sound_input.video_clock != self.input.video_clock
        {
            return Err(AudioResamplingError::SoundConfigurationMismatch(
                "audio/video clocks",
            ));
        }
        Ok(())
    }

    /// Publish the next bounded output range transactionally. The predecessor
    /// checkpoint remains valid after any refusal or cancellation.
    pub fn resample_next_chunk(
        &self,
        sound: &SoundSynthesisConfig,
        prior: &AudioResamplingCheckpoint,
        maximum_frames: NonZeroUsize,
        cx: &Cx<'_>,
    ) -> Result<AudioResamplingChunk, AudioResamplingError> {
        self.resample_next_chunk_with_checkpoint(sound, prior, maximum_frames, &mut || {
            checkpoint(cx)
        })
    }

    fn resample_next_chunk_with_checkpoint(
        &self,
        sound: &SoundSynthesisConfig,
        prior: &AudioResamplingCheckpoint,
        maximum_frames: NonZeroUsize,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<AudioResamplingChunk, AudioResamplingError> {
        checkpoint_fn()?;
        self.validate_sound_configuration(sound)?;
        self.validate_checkpoint(prior)?;
        if prior.next_audio_frame_offset == self.total_audio_frames {
            return Err(AudioResamplingError::Complete);
        }
        if maximum_frames.get() > self.input.budget.maximum_chunk_audio_frames {
            return Err(AudioResamplingError::BudgetExceeded {
                artifact: "chunk audio frames",
                requested: maximum_frames.get() as u64,
                limit: self.input.budget.maximum_chunk_audio_frames as u64,
            });
        }
        let start = prior.next_audio_frame_offset;
        let requested = u64::try_from(maximum_frames.get()).map_err(|_| {
            AudioResamplingError::BudgetExceeded {
                artifact: "chunk audio frames",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_audio_frames as u64,
            }
        })?;
        let end = start
            .checked_add(requested)
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk audio frame successor",
                requested: u64::MAX,
                limit: self.total_audio_frames,
            })?
            .min(self.total_audio_frames);
        let count =
            usize::try_from(end - start).map_err(|_| AudioResamplingError::BudgetExceeded {
                artifact: "chunk audio frames",
                requested: end - start,
                limit: usize::MAX as u64,
            })?;
        let mode_values =
            count
                .checked_mul(self.modes.len())
                .ok_or(AudioResamplingError::BudgetExceeded {
                    artifact: "chunk mode values",
                    requested: u64::MAX,
                    limit: self.input.budget.maximum_chunk_mode_values as u64,
                })?;
        if mode_values > self.input.budget.maximum_chunk_mode_values {
            return Err(AudioResamplingError::BudgetExceeded {
                artifact: "chunk mode values",
                requested: mode_values as u64,
                limit: self.input.budget.maximum_chunk_mode_values as u64,
            });
        }
        let radius = i64::from(self.input.filter.half_length);
        let (raw_first, raw_last) =
            reflected_source_span(start, end, radius, self.total_audio_frames);
        let raw_count_u64 = raw_last
            .checked_sub(raw_first)
            .and_then(|value| value.checked_add(1))
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "raw FIR halo frames",
                requested: u64::MAX,
                limit: self.input.budget.maximum_total_audio_frames,
            })?;
        let raw_count =
            usize::try_from(raw_count_u64).map_err(|_| AudioResamplingError::BudgetExceeded {
                artifact: "raw FIR halo frames",
                requested: raw_count_u64,
                limit: usize::MAX as u64,
            })?;
        let raw_mode_values = raw_count.checked_mul(self.modes.len()).ok_or(
            AudioResamplingError::BudgetExceeded {
                artifact: "raw FIR halo mode values",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_mode_values as u64,
            },
        )?;
        if raw_mode_values > self.input.budget.maximum_chunk_mode_values {
            return Err(AudioResamplingError::BudgetExceeded {
                artifact: "raw FIR halo mode values",
                requested: raw_mode_values as u64,
                limit: self.input.budget.maximum_chunk_mode_values as u64,
            });
        }

        let tap_count = self.filter.coefficients.len() as u64;
        let mode_count = self.modes.len() as u64;
        let filtered_mode_banks = 1_u64 + u64::from(self.artistic_filter.is_some());
        let filtered_coordinates = mode_count
            .checked_mul(filtered_mode_banks)
            .and_then(|value| value.checked_add(3))
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            })?;
        let filter_work = (count as u64)
            .checked_mul(tap_count)
            .and_then(|value| value.checked_mul(filtered_coordinates))
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            })?;
        // Two contiguous partitions have at most N + M - 1 nonempty
        // intersections. Using every retained source interval is conservative
        // for a chunk-local raw span and avoids an unbounded preflight scan.
        let overlap_bound = raw_count_u64
            .checked_add(self.intervals.len() as u64)
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            })?;
        let localized_raster_coordinates = mode_count
            .checked_mul(4)
            .and_then(|value| value.checked_add(3))
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            })?;
        let raster_work = overlap_bound
            .checked_mul(localized_raster_coordinates)
            .and_then(|value| {
                self.artistic_filter.as_ref().map_or(Some(value), |_| {
                    raw_count_u64
                        .checked_mul(mode_count)
                        .and_then(|texture| texture.checked_mul(4))
                        .and_then(|texture| value.checked_add(texture))
                })
            })
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            })?;
        let event_work = (self.events.len() as u64)
            .checked_mul(mode_count)
            .and_then(|value| value.checked_mul(2))
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            })?;
        let work = filter_work
            .checked_add(raster_work)
            .and_then(|value| value.checked_add(event_work))
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: u64::MAX,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            })?;
        if work > self.input.budget.maximum_chunk_multiply_adds {
            return Err(AudioResamplingError::BudgetExceeded {
                artifact: "chunk multiply-adds",
                requested: work,
                limit: self.input.budget.maximum_chunk_multiply_adds,
            });
        }
        let raw = self.build_raw_span(raw_first, raw_last, checkpoint_fn)?;
        checkpoint_fn()?;

        let mut drive_frames = Vec::new();
        let mut localized_force = Vec::new();
        let mut localized_impulse = Vec::new();
        reserve_exact(&mut drive_frames, count, "distributed drive frames")?;
        reserve_exact(&mut localized_force, mode_values, "localized mode force")?;
        reserve_exact(
            &mut localized_impulse,
            mode_values,
            "localized mode impulse",
        )?;
        localized_impulse.resize(mode_values, 0.0);

        let mut maximum_abs_distributed_force_n = 0.0_f64;
        let mut maximum_abs_localized_mode_force_n = 0.0_f64;
        for local_frame in 0..count {
            if local_frame % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let global = start + local_frame as u64;
            let distributed = convolve_components(
                global,
                &raw,
                &self.filter.coefficients,
                self.total_audio_frames,
            );
            if !components_finite(distributed) {
                return Err(AudioResamplingError::NonFiniteResult {
                    frame: global,
                    field: "distributed filtered force",
                });
            }
            maximum_abs_distributed_force_n =
                maximum_abs_distributed_force_n.max(max_abs_components(distributed));
            drive_frames.push(ModalDriveFrame {
                localized_generalized_force_n: ModalComponentValues::ZERO,
                distributed_generalized_force_n: distributed,
                localized_boundary_impulse_n_s: ModalComponentValues::ZERO,
                distributed_boundary_impulse_n_s: ModalComponentValues::ZERO,
            });
            for mode_index in 0..self.modes.len() {
                if mode_index % AUDIO_RESAMPLING_CANCELLATION_POLL_MODES == 0 {
                    checkpoint_fn()?;
                }
                let physical = convolve_mode(
                    global,
                    mode_index,
                    &raw,
                    &self.filter.coefficients,
                    self.total_audio_frames,
                    false,
                );
                let artistic = self.artistic_filter.as_ref().map_or(0.0, |kernel| {
                    convolve_mode(
                        global,
                        mode_index,
                        &raw,
                        &kernel.coefficients,
                        self.total_audio_frames,
                        true,
                    )
                });
                let total = physical + artistic;
                if !total.is_finite() {
                    return Err(AudioResamplingError::NonFiniteResult {
                        frame: global,
                        field: "localized filtered modal force",
                    });
                }
                maximum_abs_localized_mode_force_n =
                    maximum_abs_localized_mode_force_n.max(total.abs());
                localized_force.push(total);
            }
        }
        self.add_event_impulses(start, end, &mut localized_impulse, checkpoint_fn)?;
        let events = self.chunk_events(start, end, checkpoint_fn)?;
        let sync_markers = self.chunk_markers(start, end, checkpoint_fn)?;
        checkpoint_fn()?;
        let successor = AudioResamplingCheckpoint {
            model_identity: self.identity,
            next_audio_frame_offset: end,
        };
        let diagnostics = AudioResamplingDiagnostics {
            start_audio_frame_offset: start,
            end_audio_frame_offset: end,
            filter_half_length_frames: u32::from(self.input.filter.half_length),
            maximum_abs_distributed_force_n,
            maximum_abs_localized_mode_force_n,
            owned_event_count: events.len(),
        };
        Ok(AudioResamplingChunk {
            modal_identity: self.modal_identity,
            start_audio_frame_offset: start,
            identity: chunk_identity(self.identity, start, end),
            drive_frames,
            preparticipated_localized_force_n: localized_force,
            preparticipated_localized_impulse_n_s: localized_impulse,
            events,
            sync_markers,
            diagnostics,
            successor,
        })
    }

    fn validate_checkpoint(
        &self,
        checkpoint: &AudioResamplingCheckpoint,
    ) -> Result<(), AudioResamplingError> {
        if checkpoint.model_identity != self.identity
            || checkpoint.next_audio_frame_offset > self.total_audio_frames
        {
            Err(AudioResamplingError::InvalidCheckpoint)
        } else {
            Ok(())
        }
    }

    fn build_raw_span(
        &self,
        first: u64,
        last: u64,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<RawAudioSpan, AudioResamplingError> {
        let count_u64 = last
            .checked_sub(first)
            .and_then(|value| value.checked_add(1))
            .ok_or(AudioResamplingError::Capacity {
                artifact: "raw FIR halo",
                requested: usize::MAX,
            })?;
        let count = usize::try_from(count_u64).map_err(|_| AudioResamplingError::Capacity {
            artifact: "raw FIR halo",
            requested: usize::MAX,
        })?;
        let mode_values =
            count
                .checked_mul(self.modes.len())
                .ok_or(AudioResamplingError::BudgetExceeded {
                    artifact: "raw FIR halo mode values",
                    requested: u64::MAX,
                    limit: self.input.budget.maximum_chunk_mode_values as u64,
                })?;
        if mode_values > self.input.budget.maximum_chunk_mode_values {
            return Err(AudioResamplingError::BudgetExceeded {
                artifact: "raw FIR halo mode values",
                requested: mode_values as u64,
                limit: self.input.budget.maximum_chunk_mode_values as u64,
            });
        }
        let mut distributed_force_n = Vec::new();
        let mut physical_localized_mode_force_n = Vec::new();
        let mut artistic_localized_mode_force_n = Vec::new();
        reserve_exact(&mut distributed_force_n, count, "raw distributed FIR halo")?;
        reserve_exact(
            &mut physical_localized_mode_force_n,
            mode_values,
            "raw physical localized FIR halo",
        )?;
        reserve_exact(
            &mut artistic_localized_mode_force_n,
            mode_values,
            "raw artistic localized FIR halo",
        )?;
        let mut localized_sums = Vec::new();
        reserve_exact(
            &mut localized_sums,
            self.modes.len(),
            "raw localized compensated sums",
        )?;
        localized_sums.resize(self.modes.len(), CompensatedSum::new());
        for local in 0..count {
            if local % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            localized_sums.fill(CompensatedSum::new());
            distributed_force_n.push(self.raster_frame(
                first + local as u64,
                &mut localized_sums,
                &mut physical_localized_mode_force_n,
                &mut artistic_localized_mode_force_n,
                checkpoint_fn,
            )?);
        }
        Ok(RawAudioSpan {
            first_frame_offset: first,
            mode_count: self.modes.len(),
            distributed_force_n,
            physical_localized_mode_force_n,
            artistic_localized_mode_force_n,
        })
    }

    fn raster_frame(
        &self,
        frame_offset: u64,
        localized: &mut [CompensatedSum],
        physical_output: &mut Vec<f64>,
        artistic_output: &mut Vec<f64>,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<ModalComponentValues, AudioResamplingError> {
        let dt = self.audio_frame_period_s;
        let t0 = (frame_offset as f64 - self.source_start_offset_audio_frames) * dt;
        let t1 = ((frame_offset + 1) as f64 - self.source_start_offset_audio_frames) * dt;
        if !t0.is_finite() || !t1.is_finite() || t1 <= t0 {
            return Err(AudioResamplingError::NonFiniteResult {
                frame: frame_offset,
                field: "relative audio cell",
            });
        }
        let mut distributed = [CompensatedSum::new(); 3];
        let source_start = self.grid.start_time_s;
        let mut index = self
            .intervals
            .partition_point(|interval| interval.end_time_s - source_start <= t0);
        let mut traversed_intervals = 0_usize;
        while let Some(interval) = self.intervals.get(index) {
            if traversed_intervals % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let interval_start = interval.start_time_s - source_start;
            let interval_end = interval.end_time_s - source_start;
            if interval_start >= t1 {
                break;
            }
            let overlap_start = interval_start.max(t0);
            let overlap_end = interval_end.min(t1);
            if overlap_end > overlap_start {
                let overlap = overlap_end - overlap_start;
                let distributed_measure = interval.distributed_force_time_measure_n_s();
                distributed[0].add(distributed_measure.disc * overlap / interval.duration_s);
                distributed[1].add(distributed_measure.glass_plate * overlap / interval.duration_s);
                distributed[2]
                    .add(distributed_measure.base_assembly * overlap / interval.duration_s);
                let localized_measure = interval.localized_force_time_measure_n_s();
                for (mode_index, (mode, envelope)) in self
                    .modes
                    .iter()
                    .zip(&interval.spatial_envelopes)
                    .enumerate()
                {
                    if mode_index % AUDIO_RESAMPLING_CANCELLATION_POLL_MODES == 0 {
                        checkpoint_fn()?;
                    }
                    let factor_integral = linear_factor_integral(
                        interval_start,
                        interval.duration_s,
                        overlap_start,
                        overlap_end,
                        envelope.start_factor,
                        envelope.end_factor,
                    );
                    let participated_measure =
                        participation_dot(mode.source_participation, localized_measure);
                    localized[mode_index]
                        .add(participated_measure * factor_integral / interval.duration_s);
                }
            }
            index += 1;
            traversed_intervals += 1;
        }
        let distributed_force_n = ModalComponentValues {
            disc: distributed[0].total() / dt,
            glass_plate: distributed[1].total() / dt,
            base_assembly: distributed[2].total() / dt,
        };
        physical_output.extend(localized.iter().copied().map(|sum| sum.total() / dt));
        self.texture_frame(frame_offset, 0.5, artistic_output, checkpoint_fn)?;
        Ok(distributed_force_n)
    }

    fn texture_frame(
        &self,
        frame_offset: u64,
        alpha_in_frame: f64,
        output: &mut Vec<f64>,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<(), AudioResamplingError> {
        if self.artistic_filter.is_none() {
            output.resize(output.len() + self.modes.len(), 0.0);
            return Ok(());
        }
        let relative_time_s = (frame_offset as f64 + alpha_in_frame
            - self.source_start_offset_audio_frames)
            * self.audio_frame_period_s;
        let source_start = self.grid.start_time_s;
        let index = self
            .intervals
            .partition_point(|interval| interval.end_time_s - source_start <= relative_time_s)
            .min(self.intervals.len() - 1);
        let interval = &self.intervals[index];
        let Some(texture) = interval.artistic_texture else {
            output.resize(output.len() + self.modes.len(), 0.0);
            return Ok(());
        };
        let interval_start = interval.start_time_s - source_start;
        let interval_alpha =
            ((relative_time_s - interval_start) / interval.duration_s).clamp(0.0, 1.0);
        let noise = procedural_texture_unit_sample(texture.stream_identity, frame_offset);
        for (mode_index, (mode, envelope)) in self
            .modes
            .iter()
            .zip(&interval.spatial_envelopes)
            .enumerate()
        {
            if mode_index % AUDIO_RESAMPLING_CANCELLATION_POLL_MODES == 0 {
                checkpoint_fn()?;
            }
            let factor = (envelope.end_factor - envelope.start_factor)
                .mul_add(interval_alpha, envelope.start_factor);
            output.push(
                participation_dot(mode.source_participation, texture.peak_force_envelope_n)
                    * factor
                    * noise,
            );
        }
        Ok(())
    }

    fn add_event_impulses(
        &self,
        start: u64,
        end: u64,
        output: &mut [f64],
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<(), AudioResamplingError> {
        for (event_index, event) in self.events.iter().enumerate() {
            if event_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            for (target, weight) in [
                (event.receipt.left_frame_offset, event.receipt.left_weight),
                (event.receipt.right_frame_offset, event.receipt.right_weight),
            ] {
                let Some(target) = target else { continue };
                if target < start || target >= end || weight == 0.0 {
                    continue;
                }
                let local_frame = usize::try_from(target - start).map_err(|_| {
                    AudioResamplingError::InvalidEvent {
                        event: event_index,
                        field: "planned event frame is not addressable in its output chunk",
                    }
                })?;
                let row = local_frame.checked_mul(self.modes.len()).ok_or(
                    AudioResamplingError::BudgetExceeded {
                        artifact: "localized event impulse row",
                        requested: u64::MAX,
                        limit: self.input.budget.maximum_chunk_mode_values as u64,
                    },
                )?;
                for (mode_index, impulse) in
                    event.localized_mode_impulse_n_s.iter().copied().enumerate()
                {
                    if mode_index % AUDIO_RESAMPLING_CANCELLATION_POLL_MODES == 0 {
                        checkpoint_fn()?;
                    }
                    let output_index =
                        row.checked_add(mode_index)
                            .ok_or(AudioResamplingError::InvalidEvent {
                                event: event_index,
                                field: "planned event exceeds its output impulse buffer",
                            })?;
                    let slot =
                        output
                            .get_mut(output_index)
                            .ok_or(AudioResamplingError::InvalidEvent {
                                event: event_index,
                                field: "planned event exceeds its output impulse buffer",
                            })?;
                    let next = *slot + impulse * weight;
                    if !next.is_finite() {
                        return Err(AudioResamplingError::NonFiniteResult {
                            frame: target,
                            field: "localized event impulse",
                        });
                    }
                    *slot = next;
                }
            }
        }
        Ok(())
    }

    fn chunk_events(
        &self,
        start: u64,
        end: u64,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<Vec<ResampledAudioEvent>, AudioResamplingError> {
        let mut owned_count = 0_usize;
        for (event_index, event) in self.events.iter().enumerate() {
            if event_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            if (start..end).contains(&event.owner_frame_offset) {
                owned_count += 1;
            }
        }
        let mut result = Vec::new();
        reserve_exact(&mut result, owned_count, "chunk event receipts")?;
        for (event_index, event) in self.events.iter().enumerate() {
            if event_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            if (start..end).contains(&event.owner_frame_offset) {
                result.push(event.receipt.clone());
            }
        }
        Ok(result)
    }

    fn chunk_markers(
        &self,
        start: u64,
        end: u64,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
    ) -> Result<Vec<AudioVideoSyncMarker>, AudioResamplingError> {
        let terminal = end == self.total_audio_frames;
        let mut owned_count = 0_usize;
        for (marker_index, marker) in self.alignment.markers.iter().enumerate() {
            if marker_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            if (start..end).contains(&marker.audio_frame_offset)
                || (terminal && marker.audio_frame_offset == end)
            {
                owned_count += 1;
            }
        }
        let mut result = Vec::new();
        reserve_exact(&mut result, owned_count, "chunk synchronization markers")?;
        for (marker_index, marker) in self.alignment.markers.iter().copied().enumerate() {
            if marker_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            if (start..end).contains(&marker.audio_frame_offset)
                || (terminal && marker.audio_frame_offset == end)
            {
                result.push(marker);
            }
        }
        Ok(result)
    }
}

fn validate_budget(budget: AudioResamplingBudget) -> Result<(), AudioResamplingError> {
    for (name, value) in [
        ("maximum_source_intervals", budget.maximum_source_intervals),
        (
            "maximum_chunk_audio_frames",
            budget.maximum_chunk_audio_frames,
        ),
        ("maximum_filter_taps", budget.maximum_filter_taps),
        (
            "maximum_chunk_mode_values",
            budget.maximum_chunk_mode_values,
        ),
        ("maximum_events", budget.maximum_events),
        ("maximum_sync_markers", budget.maximum_sync_markers),
    ] {
        if value == 0 {
            return Err(AudioResamplingError::InvalidBudget(name));
        }
    }
    if budget.maximum_total_audio_frames == 0 {
        return Err(AudioResamplingError::InvalidBudget(
            "maximum_total_audio_frames",
        ));
    }
    if budget.maximum_chunk_multiply_adds == 0 {
        return Err(AudioResamplingError::InvalidBudget(
            "maximum_chunk_multiply_adds",
        ));
    }
    if budget.maximum_filter_taps > MAX_AUDIO_RECONSTRUCTION_FILTER_TAPS {
        return Err(AudioResamplingError::InvalidBudget(
            "maximum_filter_taps hard ceiling",
        ));
    }
    Ok(())
}

fn validate_clocks(
    video: CinematicClock,
    audio: CinematicClock,
    budget: AudioResamplingBudget,
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<(u64, AudioVideoAlignment, f64), AudioResamplingError> {
    if video.domain() != CinematicClockDomain::Video
        || video.ticks_per_second_numerator() != SOUND_MASTER_VIDEO_RATE_HZ
        || video.ticks_per_second_denominator() != 1
    {
        return Err(AudioResamplingError::InvalidMasterClock(
            "video clock must be exact 24/1 Hz",
        ));
    }
    if audio.domain() != CinematicClockDomain::Audio
        || audio.ticks_per_second_numerator() != SOUND_MASTER_SAMPLE_RATE_HZ
        || audio.ticks_per_second_denominator() != 1
    {
        return Err(AudioResamplingError::InvalidMasterClock(
            "audio clock must be exact 48000/1 Hz",
        ));
    }
    if !same_rational_instant(video, video.start_tick(), audio, audio.start_tick())
        || !same_rational_instant(
            video,
            video.end_tick_exclusive(),
            audio,
            audio.end_tick_exclusive(),
        )
    {
        return Err(AudioResamplingError::InvalidMasterClock(
            "audio/video endpoints differ",
        ));
    }
    let audio_count_i128 = i128::from(audio.end_tick_exclusive()) - i128::from(audio.start_tick());
    if audio_count_i128 <= 0 {
        return Err(AudioResamplingError::InvalidMasterClock(
            "empty audio range",
        ));
    }
    let total_audio_frames =
        u64::try_from(audio_count_i128).map_err(|_| AudioResamplingError::BudgetExceeded {
            artifact: "total audio frames",
            requested: u64::MAX,
            limit: budget.maximum_total_audio_frames,
        })?;
    if total_audio_frames > budget.maximum_total_audio_frames {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "total audio frames",
            requested: total_audio_frames,
            limit: budget.maximum_total_audio_frames,
        });
    }
    let video_count_i128 = i128::from(video.end_tick_exclusive()) - i128::from(video.start_tick());
    if video_count_i128 <= 0 {
        return Err(AudioResamplingError::InvalidMasterClock(
            "empty video range",
        ));
    }
    let marker_count_i128 = video_count_i128 + 1;
    let marker_count =
        usize::try_from(marker_count_i128).map_err(|_| AudioResamplingError::Capacity {
            artifact: "A/V synchronization markers",
            requested: usize::MAX,
        })?;
    if marker_count > budget.maximum_sync_markers {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "A/V synchronization markers",
            requested: marker_count as u64,
            limit: budget.maximum_sync_markers as u64,
        });
    }
    let ratio_numerator = u64::from(audio.ticks_per_second_numerator())
        * u64::from(video.ticks_per_second_denominator());
    let ratio_denominator = u64::from(video.ticks_per_second_numerator())
        * u64::from(audio.ticks_per_second_denominator());
    if ratio_numerator % ratio_denominator != 0 {
        return Err(AudioResamplingError::InvalidMasterClock(
            "video boundaries are not integral audio boundaries",
        ));
    }
    let audio_frames_per_video_frame = u32::try_from(ratio_numerator / ratio_denominator)
        .map_err(|_| AudioResamplingError::InvalidMasterClock("clock ratio overflow"))?;
    let mut markers = Vec::new();
    reserve_exact(&mut markers, marker_count, "A/V synchronization markers")?;
    for (marker_index, video_tick) in (video.start_tick()..=video.end_tick_exclusive()).enumerate()
    {
        if marker_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        let numerator = i128::from(video_tick)
            * i128::from(video.ticks_per_second_denominator())
            * i128::from(audio.ticks_per_second_numerator());
        let denominator = i128::from(video.ticks_per_second_numerator())
            * i128::from(audio.ticks_per_second_denominator());
        if numerator % denominator != 0 {
            return Err(AudioResamplingError::InvalidMasterClock(
                "nonintegral interior video boundary",
            ));
        }
        let audio_tick = i64::try_from(numerator / denominator)
            .map_err(|_| AudioResamplingError::InvalidMasterClock("audio tick overflow"))?;
        let offset = u64::try_from(i128::from(audio_tick) - i128::from(audio.start_tick()))
            .map_err(|_| AudioResamplingError::InvalidMasterClock("negative audio marker"))?;
        markers.push(AudioVideoSyncMarker {
            video_tick,
            audio_tick,
            audio_frame_offset: offset,
        });
    }
    if markers.last().map(|marker| marker.audio_frame_offset) != Some(total_audio_frames) {
        return Err(AudioResamplingError::InvalidMasterClock(
            "exclusive endpoint marker mismatch",
        ));
    }
    Ok((
        total_audio_frames,
        AudioVideoAlignment {
            audio_frames_per_video_frame,
            markers,
            endpoint_drift_audio_frames: 0,
        },
        f64::from(audio.ticks_per_second_denominator())
            / f64::from(audio.ticks_per_second_numerator()),
    ))
}

fn validate_source_clock(
    grid: AudioExcitationGrid,
    audio: CinematicClock,
) -> Result<f64, AudioResamplingError> {
    let clock_start = clock_tick_time_s(audio, audio.start_tick());
    let clock_end = clock_tick_time_s(audio, audio.end_tick_exclusive());
    let rate = f64::from(audio.ticks_per_second_numerator())
        / f64::from(audio.ticks_per_second_denominator());
    for (endpoint, source, expected) in [
        ("start", grid.start_time_s, clock_start),
        ("end", grid.end_time_s, clock_end),
    ] {
        let error_frames = (source - expected).abs() * rate;
        if !error_frames.is_finite() || error_frames > MAX_SOURCE_CLOCK_ALIGNMENT_ERROR_FRAMES {
            return Err(AudioResamplingError::SourceClockMismatch {
                endpoint,
                error_frames,
            });
        }
    }
    let source_start_offset_audio_frames = (grid.start_time_s - clock_start) * rate;
    if !source_start_offset_audio_frames.is_finite() {
        return Err(AudioResamplingError::SourceClockMismatch {
            endpoint: "start",
            error_frames: source_start_offset_audio_frames,
        });
    }
    Ok(source_start_offset_audio_frames)
}

fn validate_source_timeline(
    grid: AudioExcitationGrid,
    modes: &[SoundMode],
    intervals: &[AudioExcitationInterval],
    budget: AudioResamplingBudget,
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<(), AudioResamplingError> {
    if intervals.is_empty() || intervals.len() != grid.interval_count {
        return Err(AudioResamplingError::InvalidSourceTimeline {
            interval: intervals.len(),
            field: "complete selected interval count",
        });
    }
    if intervals.len() > budget.maximum_source_intervals {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "source intervals",
            requested: intervals.len() as u64,
            limit: budget.maximum_source_intervals as u64,
        });
    }
    if !grid.minimum_interval_duration_s.is_finite()
        || !grid.maximum_interval_duration_s.is_finite()
        || !grid.nominal_mechanics_timestep_s.is_finite()
        || !grid.nominal_source_nyquist_ceiling_hz.is_finite()
        || grid.minimum_interval_duration_s <= 0.0
        || grid.maximum_interval_duration_s < grid.minimum_interval_duration_s
        || grid.nominal_mechanics_timestep_s <= 0.0
        || grid.reconstruction != AudioExcitationReconstructionStatus::RequiresBandLimitedResampling
    {
        return Err(AudioResamplingError::InvalidSourceTimeline {
            interval: intervals.len(),
            field: "source grid metadata",
        });
    }
    let mut event_count = 0_usize;
    let mut previous_end_bits = None;
    let mut previous_source_index: Option<usize> = None;
    let mut previous_event_time: Option<f64> = None;
    let mut minimum_duration_s = f64::INFINITY;
    let mut maximum_duration_s = 0.0_f64;
    for (interval_index, interval) in intervals.iter().enumerate() {
        if interval_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        if !interval.start_time_s.is_finite()
            || !interval.end_time_s.is_finite()
            || !interval.duration_s.is_finite()
            || interval.duration_s <= 0.0
            || interval.end_time_s <= interval.start_time_s
            || interval.end_time_s - interval.start_time_s != interval.duration_s
        {
            return Err(AudioResamplingError::InvalidSourceTimeline {
                interval: interval_index,
                field: "finite positive interval geometry",
            });
        }
        minimum_duration_s = minimum_duration_s.min(interval.duration_s);
        maximum_duration_s = maximum_duration_s.max(interval.duration_s);
        if previous_end_bits.is_some_and(|bits| bits != interval.start_time_s.to_bits()) {
            return Err(AudioResamplingError::InvalidSourceTimeline {
                interval: interval_index,
                field: "bit-exact interval continuity",
            });
        }
        previous_end_bits = Some(interval.end_time_s.to_bits());
        if interval.first_source_sample_index > interval.last_source_sample_index
            || previous_source_index.is_some_and(|previous| {
                previous.checked_add(1) != Some(interval.first_source_sample_index)
            })
        {
            return Err(AudioResamplingError::InvalidSourceTimeline {
                interval: interval_index,
                field: "contiguous source sample coverage",
            });
        }
        previous_source_index = Some(interval.last_source_sample_index);
        if interval.spatial_envelopes.len() != modes.len() {
            return Err(AudioResamplingError::InvalidSourceTimeline {
                interval: interval_index,
                field: "spatial envelope mode count",
            });
        }
        for (mode_index, (mode, envelope)) in
            modes.iter().zip(&interval.spatial_envelopes).enumerate()
        {
            if mode_index % AUDIO_RESAMPLING_CANCELLATION_POLL_MODES == 0 {
                checkpoint_fn()?;
            }
            if envelope.mode_id != mode.mode_id
                || !envelope.start_factor.is_finite()
                || !envelope.end_factor.is_finite()
                || envelope.start_factor.abs() > MAX_MODAL_SPATIAL_PARTICIPATION
                || envelope.end_factor.abs() > MAX_MODAL_SPATIAL_PARTICIPATION
            {
                return Err(AudioResamplingError::InvalidSourceTimeline {
                    interval: interval_index,
                    field: "canonical finite spatial envelope",
                });
            }
        }
        for values in [
            interval.mean_force_stems_n.contact,
            interval.mean_force_stems_n.rolling,
            interval.mean_force_stems_n.base,
            interval.mean_force_stems_n.gas,
            interval.force_time_stems_n_s.contact,
            interval.force_time_stems_n_s.rolling,
            interval.force_time_stems_n_s.base,
            interval.force_time_stems_n_s.gas,
            interval.mean_generalized_force_n,
            interval.generalized_force_time_n_s,
            interval.measure_residual_stems_n_s.contact,
            interval.measure_residual_stems_n_s.rolling,
            interval.measure_residual_stems_n_s.base,
            interval.measure_residual_stems_n_s.gas,
            interval.measure_residual_n_s,
            interval.localized_force_time_measure_n_s(),
            interval.distributed_force_time_measure_n_s(),
        ] {
            if !components_finite(values) {
                return Err(AudioResamplingError::InvalidSourceTimeline {
                    interval: interval_index,
                    field: "finite force-time measures",
                });
            }
        }
        if interval.mean_generalized_force_n != interval.mean_force_stems_n.sum()
            || interval.generalized_force_time_n_s != interval.force_time_stems_n_s.sum()
            || interval.measure_residual_n_s != interval.measure_residual_stems_n_s.sum()
        {
            return Err(AudioResamplingError::InvalidSourceTimeline {
                interval: interval_index,
                field: "source stem aggregate consistency",
            });
        }
        if interval.event_barrier != !interval.events.is_empty() {
            return Err(AudioResamplingError::InvalidSourceTimeline {
                interval: interval_index,
                field: "event barrier consistency",
            });
        }
        if let Some(texture) = interval.artistic_texture {
            if !components_finite(texture.peak_force_envelope_n)
                || !texture.band_low_hz.is_finite()
                || !texture.band_high_hz.is_finite()
            {
                return Err(AudioResamplingError::InvalidSourceTimeline {
                    interval: interval_index,
                    field: "finite artistic texture envelope",
                });
            }
        }
        let mut local_event_time = None;
        for (local_event_index, event) in interval.events.iter().enumerate() {
            if local_event_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let artistic_finite = event
                .artistic
                .is_none_or(|artistic| components_finite(artistic.impulse_n_s));
            if event.measure != ContactEventMeasure::TimingOnly
                || event.source_sample_index < interval.first_source_sample_index
                || event.source_sample_index > interval.last_source_sample_index
                || !event.time_s.is_finite()
                || !event.bracket_start_s.is_finite()
                || !event.bracket_end_s.is_finite()
                || event.bracket_start_s > event.time_s
                || event.time_s > event.bracket_end_s
                || event.time_s < interval.start_time_s
                || event.time_s > interval.end_time_s
                || local_event_time.is_some_and(|time| event.time_s <= time)
                || previous_event_time.is_some_and(|time| event.time_s <= time)
                || event.physical_impulse_n_s != ModalComponentValues::ZERO
                || (event.kind == ContactTransitionKind::Opening && event.artistic.is_some())
                || !artistic_finite
            {
                return Err(AudioResamplingError::InvalidEvent {
                    event: event_count.saturating_add(local_event_index),
                    field: "canonical timing-only source event",
                });
            }
            local_event_time = Some(event.time_s);
            previous_event_time = Some(event.time_s);
        }
        event_count = event_count.checked_add(interval.events.len()).ok_or(
            AudioResamplingError::BudgetExceeded {
                artifact: "source events",
                requested: u64::MAX,
                limit: budget.maximum_events as u64,
            },
        )?;
    }
    if intervals[0].start_time_s.to_bits() != grid.start_time_s.to_bits()
        || intervals
            .last()
            .is_none_or(|interval| interval.end_time_s.to_bits() != grid.end_time_s.to_bits())
    {
        return Err(AudioResamplingError::InvalidSourceTimeline {
            interval: intervals.len(),
            field: "grid endpoint binding",
        });
    }
    if minimum_duration_s.to_bits() != grid.minimum_interval_duration_s.to_bits()
        || maximum_duration_s.to_bits() != grid.maximum_interval_duration_s.to_bits()
        || (0.5 / maximum_duration_s).to_bits() != grid.nominal_source_nyquist_ceiling_hz.to_bits()
    {
        return Err(AudioResamplingError::InvalidSourceTimeline {
            interval: intervals.len(),
            field: "grid duration and cadence binding",
        });
    }
    if event_count > budget.maximum_events {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "source events",
            requested: event_count as u64,
            limit: budget.maximum_events as u64,
        });
    }
    let event_mode_values =
        event_count
            .checked_mul(modes.len())
            .ok_or(AudioResamplingError::BudgetExceeded {
                artifact: "event mode values",
                requested: u64::MAX,
                limit: budget.maximum_chunk_mode_values as u64,
            })?;
    if event_mode_values > budget.maximum_chunk_mode_values {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "event mode values",
            requested: event_mode_values as u64,
            limit: budget.maximum_chunk_mode_values as u64,
        });
    }
    checkpoint_fn()
}

fn design_physical_filter(
    spec: AudioReconstructionFilterSpec,
    nominal_source_nyquist_hz: f64,
    budget: AudioResamplingBudget,
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<FirKernel, AudioResamplingError> {
    let target_nyquist_hz = 0.5 * f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    if !spec.passband_edge_hz.is_finite()
        || !spec.stopband_edge_hz.is_finite()
        || spec.passband_edge_hz <= 0.0
        || spec.stopband_edge_hz <= spec.passband_edge_hz
        || spec.stopband_edge_hz > nominal_source_nyquist_hz
        || spec.stopband_edge_hz >= target_nyquist_hz
    {
        return Err(AudioResamplingError::InvalidFilter(
            "ordered pass/stop edges within source and target Nyquist",
        ));
    }
    if !spec.maximum_passband_ripple_db.is_finite()
        || spec.maximum_passband_ripple_db <= 0.0
        || spec.maximum_passband_ripple_db > MAX_AUDIO_FILTER_PASSBAND_RIPPLE_DB
        || !spec.minimum_stopband_attenuation_db.is_finite()
        || spec.minimum_stopband_attenuation_db < MIN_AUDIO_FILTER_STOPBAND_ATTENUATION_DB
    {
        return Err(AudioResamplingError::InvalidFilter(
            "production ripple/attenuation contract",
        ));
    }
    if !(8_192..=32_768).contains(&spec.response_grid_intervals) {
        return Err(AudioResamplingError::InvalidFilter(
            "response grid interval count",
        ));
    }
    let tap_count = usize::from(spec.half_length)
        .checked_mul(2)
        .and_then(|value| value.checked_add(1))
        .ok_or(AudioResamplingError::InvalidFilter("tap count overflow"))?;
    if spec.half_length < 4 || tap_count > budget.maximum_filter_taps {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "filter taps",
            requested: tap_count as u64,
            limit: budget.maximum_filter_taps as u64,
        });
    }
    let minimum_response_grid = u32::from(spec.half_length)
        .checked_mul(8)
        .map_or(u32::MAX, |value| value.max(8_192));
    if spec.response_grid_intervals < minimum_response_grid {
        return Err(AudioResamplingError::InvalidFilter(
            "response grid too coarse for filter length",
        ));
    }
    let response_work = u64::from(spec.response_grid_intervals)
        .checked_add(1)
        .and_then(|value| value.checked_mul(tap_count as u64))
        .and_then(|value| value.checked_mul(2))
        .ok_or(AudioResamplingError::BudgetExceeded {
            artifact: "filter response audit multiply-adds",
            requested: u64::MAX,
            limit: budget.maximum_chunk_multiply_adds,
        })?;
    if response_work > budget.maximum_chunk_multiply_adds {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "filter response audit multiply-adds",
            requested: response_work,
            limit: budget.maximum_chunk_multiply_adds,
        });
    }
    let cutoff_hz = 0.5 * (spec.passband_edge_hz + spec.stopband_edge_hz);
    let coefficients = design_windowed_lowpass(cutoff_hz, spec.half_length, checkpoint_fn)?;
    let (ripple_db, attenuation_db) = measure_lowpass_response(
        &coefficients,
        spec.passband_edge_hz,
        spec.stopband_edge_hz,
        spec.response_grid_intervals,
        checkpoint_fn,
    )?;
    if ripple_db > spec.maximum_passband_ripple_db {
        return Err(AudioResamplingError::InvalidFilter(
            "measured passband ripple",
        ));
    }
    if attenuation_db < spec.minimum_stopband_attenuation_db {
        return Err(AudioResamplingError::InvalidFilter(
            "measured stopband attenuation",
        ));
    }
    let diagnostics = AudioReconstructionFilterDiagnostics {
        tap_count,
        measured_passband_ripple_db: ripple_db,
        measured_stopband_attenuation_db: attenuation_db,
        intrinsic_group_delay_frames: u32::from(spec.half_length),
        group_delay_compensation_frames: u32::from(spec.half_length),
        required_lookahead_frames: u32::from(spec.half_length),
        published_alignment_offset_frames: 0,
    };
    let identity = filter_identity(
        b"physical-lowpass-blackman-harris-4-v1",
        &coefficients,
        spec.passband_edge_hz,
        spec.stopband_edge_hz,
        diagnostics,
    );
    Ok(FirKernel {
        coefficients,
        identity,
        diagnostics,
    })
}

fn design_artistic_filter(
    low_hz: f64,
    high_hz: f64,
    half_length: u16,
    budget: AudioResamplingBudget,
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<FirKernel, AudioResamplingError> {
    let nyquist = 0.5 * f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    if !low_hz.is_finite()
        || !high_hz.is_finite()
        || low_hz < 0.0
        || low_hz >= high_hz
        || high_hz >= nyquist
    {
        return Err(AudioResamplingError::InvalidFilter("artistic texture band"));
    }
    let tap_count = usize::from(half_length) * 2 + 1;
    if tap_count > budget.maximum_filter_taps {
        return Err(AudioResamplingError::BudgetExceeded {
            artifact: "artistic filter taps",
            requested: tap_count as u64,
            limit: budget.maximum_filter_taps as u64,
        });
    }
    let mut coefficients = design_windowed_lowpass(high_hz, half_length, checkpoint_fn)?;
    if low_hz > 0.0 {
        let low = design_windowed_lowpass(low_hz, half_length, checkpoint_fn)?;
        for (high_coefficient, low_coefficient) in coefficients.iter_mut().zip(low) {
            *high_coefficient -= low_coefficient;
        }
        let mut sum = CompensatedSum::new();
        for value in &coefficients {
            sum.add(*value);
        }
        coefficients[usize::from(half_length)] -= sum.total();
        let center_hz = 0.5 * (low_hz + high_hz);
        let gain = centered_frequency_response(&coefficients, center_hz).abs();
        if !gain.is_finite() || gain <= f64::EPSILON {
            return Err(AudioResamplingError::InvalidFilter(
                "artistic bandpass normalization",
            ));
        }
        for value in &mut coefficients {
            *value /= gain;
        }
    }
    let diagnostics = AudioReconstructionFilterDiagnostics {
        tap_count,
        measured_passband_ripple_db: 0.0,
        measured_stopband_attenuation_db: 0.0,
        intrinsic_group_delay_frames: u32::from(half_length),
        group_delay_compensation_frames: u32::from(half_length),
        required_lookahead_frames: u32::from(half_length),
        published_alignment_offset_frames: 0,
    };
    let identity = filter_identity(
        b"artistic-texture-band-blackman-harris-4-v1",
        &coefficients,
        low_hz,
        high_hz,
        diagnostics,
    );
    Ok(FirKernel {
        coefficients,
        identity,
        diagnostics,
    })
}

fn design_windowed_lowpass(
    cutoff_hz: f64,
    half_length: u16,
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<Vec<f64>, AudioResamplingError> {
    let radius = usize::from(half_length);
    let tap_count = radius * 2 + 1;
    let normalized_cutoff = cutoff_hz / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    let mut coefficients = Vec::new();
    reserve_exact(&mut coefficients, tap_count, "filter coefficients")?;
    coefficients.resize(tap_count, 0.0);
    for offset in 0..=radius {
        if offset % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        let phase = 2.0 * normalized_cutoff * offset as f64;
        let ideal = 2.0 * normalized_cutoff * sinc_pi(phase);
        let tap = radius + offset;
        let angle = 2.0 * PI * tap as f64 / (tap_count - 1) as f64;
        let window = 0.35875 - 0.48829 * det::cos(angle) + 0.14128 * det::cos(2.0 * angle)
            - 0.01168 * det::cos(3.0 * angle);
        let coefficient = ideal * window;
        coefficients[radius + offset] = coefficient;
        coefficients[radius - offset] = coefficient;
    }
    let mut sum = CompensatedSum::new();
    for value in &coefficients {
        sum.add(*value);
    }
    let total = sum.total();
    if !total.is_finite() || total.abs() <= f64::EPSILON {
        return Err(AudioResamplingError::InvalidFilter(
            "low-pass normalization",
        ));
    }
    for value in &mut coefficients {
        *value /= total;
    }
    let mut normalized_sum = CompensatedSum::new();
    for value in &coefficients {
        normalized_sum.add(*value);
    }
    coefficients[usize::from(half_length)] += 1.0 - normalized_sum.total();
    Ok(coefficients)
}

fn measure_lowpass_response(
    coefficients: &[f64],
    passband_hz: f64,
    stopband_hz: f64,
    grid_intervals: u32,
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<(f64, f64), AudioResamplingError> {
    let mut maximum_ripple_db = 0.0_f64;
    let mut maximum_stop_amplitude = 0.0_f64;
    let nyquist = 0.5 * f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    for index in 0..=grid_intervals {
        checkpoint_fn()?;
        let alpha = f64::from(index) / f64::from(grid_intervals);
        let pass_frequency = passband_hz * alpha;
        let pass_amplitude = centered_frequency_response(coefficients, pass_frequency).abs();
        maximum_ripple_db = maximum_ripple_db.max(amplitude_db(pass_amplitude).abs());
        let stop_frequency = (nyquist - stopband_hz).mul_add(alpha, stopband_hz);
        maximum_stop_amplitude = maximum_stop_amplitude
            .max(centered_frequency_response(coefficients, stop_frequency).abs());
    }
    let attenuation_db = if maximum_stop_amplitude == 0.0 {
        f64::INFINITY
    } else {
        -amplitude_db(maximum_stop_amplitude)
    };
    Ok((maximum_ripple_db, attenuation_db))
}

fn centered_frequency_response(coefficients: &[f64], frequency_hz: f64) -> f64 {
    debug_assert!(!coefficients.is_empty() && coefficients.len() % 2 == 1);
    let radius = coefficients.len() / 2;
    let omega = 2.0 * PI * frequency_hz / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    let mut response = CompensatedSum::new();
    response.add(coefficients[radius]);
    for offset in 1..=radius {
        let symmetric_pair = coefficients[radius - offset] + coefficients[radius + offset];
        response.add(symmetric_pair * det::cos(omega * offset as f64));
    }
    response.total()
}

fn amplitude_db(amplitude: f64) -> f64 {
    20.0 * det::ln(amplitude.max(f64::MIN_POSITIVE)) / LN_10
}

fn common_texture_band(
    intervals: &[AudioExcitationInterval],
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<Option<(f64, f64)>, AudioResamplingError> {
    let mut band = None;
    for (interval_index, interval) in intervals.iter().enumerate() {
        if interval_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        let Some(texture) = interval.artistic_texture else {
            continue;
        };
        let candidate = (texture.band_low_hz, texture.band_high_hz);
        if band.is_some_and(|existing| existing != candidate) {
            return Err(AudioResamplingError::InvalidSourceTimeline {
                interval: interval_index,
                field: "single artistic texture band",
            });
        }
        band = Some(candidate);
    }
    Ok(band)
}

fn plan_events(
    intervals: &[AudioExcitationInterval],
    modes: &[SoundMode],
    grid: AudioExcitationGrid,
    source_start_offset_audio_frames: f64,
    total_audio_frames: u64,
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<Vec<PlannedEvent>, AudioResamplingError> {
    let count = intervals
        .iter()
        .try_fold(0_usize, |count, interval| {
            count.checked_add(interval.events.len())
        })
        .ok_or(AudioResamplingError::BudgetExceeded {
            artifact: "planned resampled events",
            requested: u64::MAX,
            limit: usize::MAX as u64,
        })?;
    let mut result = Vec::new();
    reserve_exact(&mut result, count, "planned resampled events")?;
    let rate = f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    let mut flat_event_index = 0_usize;
    for (interval_index, interval) in intervals.iter().enumerate() {
        if interval_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        for (local_event_index, event) in interval.events.iter().enumerate() {
            if local_event_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            if event.physical_impulse_n_s != ModalComponentValues::ZERO {
                return Err(AudioResamplingError::InvalidEvent {
                    event: flat_event_index,
                    field: "timing-only physical impulse must remain zero",
                });
            }
            if event.kind == ContactTransitionKind::Opening && event.artistic.is_some() {
                return Err(AudioResamplingError::InvalidEvent {
                    event: flat_event_index,
                    field: "opening cannot acquire artistic reimpact impulse",
                });
            }
            let requested = snap_exact_audio_frame(
                (event.time_s - grid.start_time_s).mul_add(rate, source_start_offset_audio_frames),
            );
            if !requested.is_finite() || requested < 0.0 || requested > total_audio_frames as f64 {
                return Err(AudioResamplingError::InvalidEvent {
                    event: flat_event_index,
                    field: "event outside master audio range",
                });
            }
            let artistic_component = event
                .artistic
                .map_or(ModalComponentValues::ZERO, |artistic| artistic.impulse_n_s);
            let has_artistic_impulse = max_abs_components(artistic_component) != 0.0;
            let floor = requested.floor();
            let fraction = requested - floor;
            let floor_offset = floor as u64;
            let (left, right, left_weight, right_weight, centroid_error) = if has_artistic_impulse {
                let right_candidate = floor_offset.checked_add(1);
                if floor_offset >= total_audio_frames
                    || (fraction > 0.0
                        && right_candidate.is_none_or(|right| right >= total_audio_frames))
                {
                    return Err(AudioResamplingError::EventOutsideRepresentableRange {
                        event: flat_event_index,
                    });
                }
                let right = if fraction > 0.0 {
                    Some(right_candidate.ok_or(
                        AudioResamplingError::EventOutsideRepresentableRange {
                            event: flat_event_index,
                        },
                    )?)
                } else {
                    None
                };
                let left_weight = 1.0 - fraction;
                let right_weight = fraction;
                let centroid = floor_offset as f64 * left_weight
                    + right.map_or(0.0, |offset| offset as f64 * right_weight);
                (
                    Some(floor_offset),
                    right,
                    left_weight,
                    right_weight,
                    centroid - requested,
                )
            } else {
                let left = (floor_offset < total_audio_frames).then_some(floor_offset);
                (left, None, 0.0, 0.0, 0.0)
            };
            let interval_alpha =
                ((event.time_s - interval.start_time_s) / interval.duration_s).clamp(0.0, 1.0);
            let mut localized_mode_impulse_n_s = Vec::new();
            reserve_exact(
                &mut localized_mode_impulse_n_s,
                modes.len(),
                "event localized modal impulse",
            )?;
            for (mode_index, (mode, envelope)) in
                modes.iter().zip(&interval.spatial_envelopes).enumerate()
            {
                if mode_index % AUDIO_RESAMPLING_CANCELLATION_POLL_MODES == 0 {
                    checkpoint_fn()?;
                }
                let factor = (envelope.end_factor - envelope.start_factor)
                    .mul_add(interval_alpha, envelope.start_factor);
                let impulse =
                    participation_dot(mode.source_participation, artistic_component) * factor;
                if !impulse.is_finite() {
                    return Err(AudioResamplingError::InvalidEvent {
                        event: flat_event_index,
                        field: "finite participated artistic impulse",
                    });
                }
                localized_mode_impulse_n_s.push(impulse);
            }
            // Event times and their admissible brackets share one coordinate
            // system, so exact-frame normalization must be identical for all
            // three.  Snapping only `requested` can otherwise make a valid
            // terminal bracket appear infinitesimally too short after two
            // algebraically equivalent floating-point evaluations.
            let bracket_start_sample_position = snap_exact_audio_frame(
                (event.bracket_start_s - grid.start_time_s)
                    .mul_add(rate, source_start_offset_audio_frames),
            );
            let bracket_end_sample_position = snap_exact_audio_frame(
                (event.bracket_end_s - grid.start_time_s)
                    .mul_add(rate, source_start_offset_audio_frames),
            );
            if !bracket_start_sample_position.is_finite()
                || !bracket_end_sample_position.is_finite()
                || bracket_start_sample_position > requested
                || requested > bracket_end_sample_position
            {
                return Err(AudioResamplingError::InvalidEvent {
                    event: flat_event_index,
                    field: "event bracket in sample coordinates",
                });
            }
            let owner_frame_offset = left.unwrap_or(total_audio_frames - 1);
            result.push(PlannedEvent {
                receipt: ResampledAudioEvent {
                    source: *event,
                    requested_sample_position: requested,
                    left_frame_offset: left,
                    right_frame_offset: right,
                    left_weight,
                    right_weight,
                    centroid_error_frames: centroid_error,
                    bracket_start_sample_position,
                    bracket_end_sample_position,
                },
                owner_frame_offset,
                localized_mode_impulse_n_s,
            });
            flat_event_index += 1;
        }
    }
    Ok(result)
}

fn snap_exact_audio_frame(sample_position: f64) -> f64 {
    let nearest = sample_position.round();
    if (sample_position - nearest).abs() <= EVENT_SAMPLE_SNAP_TOLERANCE_FRAMES {
        nearest
    } else {
        sample_position
    }
}

fn reflected_source_span(start: u64, end: u64, radius: i64, total: u64) -> (u64, u64) {
    let first_virtual = i128::from(start) - i128::from(radius);
    let last_virtual = i128::from(end - 1) + i128::from(radius);
    let mut minimum = u64::MAX;
    let mut maximum = 0_u64;
    for virtual_index in first_virtual..=last_virtual {
        let reflected = reflect_half_sample_even(virtual_index, total);
        minimum = minimum.min(reflected);
        maximum = maximum.max(reflected);
    }
    (minimum, maximum)
}

fn reflect_half_sample_even(index: i128, length: u64) -> u64 {
    debug_assert!(length > 0);
    let length = i128::from(length);
    let period = 2 * length;
    let reduced = index.rem_euclid(period);
    let reflected = if reduced < length {
        reduced
    } else {
        period - 1 - reduced
    };
    debug_assert!((0..length).contains(&reflected));
    reflected as u64
}

fn convolve_components(
    output_frame: u64,
    raw: &RawAudioSpan,
    coefficients: &[f64],
    total_frames: u64,
) -> ModalComponentValues {
    let radius = (coefficients.len() / 2) as i128;
    let mut sums = [CompensatedSum::new(); 3];
    for (tap, coefficient) in coefficients.iter().copied().enumerate() {
        let virtual_index = i128::from(output_frame) + tap as i128 - radius;
        let source_index = reflect_half_sample_even(virtual_index, total_frames);
        let values = raw.distributed_force_n[(source_index - raw.first_frame_offset) as usize];
        sums[0].add(coefficient * values.disc);
        sums[1].add(coefficient * values.glass_plate);
        sums[2].add(coefficient * values.base_assembly);
    }
    ModalComponentValues {
        disc: sums[0].total(),
        glass_plate: sums[1].total(),
        base_assembly: sums[2].total(),
    }
}

fn convolve_mode(
    output_frame: u64,
    mode_index: usize,
    raw: &RawAudioSpan,
    coefficients: &[f64],
    total_frames: u64,
    artistic: bool,
) -> f64 {
    let radius = (coefficients.len() / 2) as i128;
    let mut sum = CompensatedSum::new();
    for (tap, coefficient) in coefficients.iter().copied().enumerate() {
        let virtual_index = i128::from(output_frame) + tap as i128 - radius;
        let source_index = reflect_half_sample_even(virtual_index, total_frames);
        let row = (source_index - raw.first_frame_offset) as usize;
        let values = if artistic {
            &raw.artistic_localized_mode_force_n
        } else {
            &raw.physical_localized_mode_force_n
        };
        sum.add(coefficient * values[row * raw.mode_count + mode_index]);
    }
    sum.total()
}

fn linear_factor_integral(
    interval_start_s: f64,
    interval_duration_s: f64,
    overlap_start_s: f64,
    overlap_end_s: f64,
    start_factor: f64,
    end_factor: f64,
) -> f64 {
    let u0 = overlap_start_s - interval_start_s;
    let u1 = overlap_end_s - interval_start_s;
    let factor0 = (end_factor - start_factor).mul_add(u0 / interval_duration_s, start_factor);
    let factor1 = (end_factor - start_factor).mul_add(u1 / interval_duration_s, start_factor);
    0.5 * (factor0 + factor1) * (u1 - u0)
}

fn participation_dot(participation: SoundModeParticipation, values: ModalComponentValues) -> f64 {
    let mut sum = CompensatedSum::new();
    sum.add(participation.disc * values.disc);
    sum.add(participation.glass_plate * values.glass_plate);
    sum.add(participation.base_assembly * values.base_assembly);
    sum.total()
}

fn components_finite(values: ModalComponentValues) -> bool {
    values.disc.is_finite() && values.glass_plate.is_finite() && values.base_assembly.is_finite()
}

fn max_abs_components(values: ModalComponentValues) -> f64 {
    values
        .disc
        .abs()
        .max(values.glass_plate.abs())
        .max(values.base_assembly.abs())
}

fn sinc_pi(value: f64) -> f64 {
    if value.abs() <= 1.0e-8 {
        let x = PI * value;
        let squared = x * x;
        1.0 + squared * (-1.0 / 6.0 + squared * (1.0 / 120.0 - squared / 5_040.0))
    } else {
        det::sin(PI * value) / (PI * value)
    }
}

fn clock_tick_time_s(clock: CinematicClock, tick: i64) -> f64 {
    tick as f64 * f64::from(clock.ticks_per_second_denominator())
        / f64::from(clock.ticks_per_second_numerator())
}

fn same_rational_instant(
    left: CinematicClock,
    left_tick: i64,
    right: CinematicClock,
    right_tick: i64,
) -> bool {
    i128::from(left_tick)
        * i128::from(left.ticks_per_second_denominator())
        * i128::from(right.ticks_per_second_numerator())
        == i128::from(right_tick)
            * i128::from(right.ticks_per_second_denominator())
            * i128::from(left.ticks_per_second_numerator())
}

fn source_payload_identity(
    intervals: &[AudioExcitationInterval],
    checkpoint_fn: &mut impl FnMut() -> Result<(), AudioResamplingError>,
) -> Result<ContentHash, AudioResamplingError> {
    let mut hasher = DomainHasher::new(SOURCE_PAYLOAD_IDENTITY_DOMAIN);
    hash_u64(&mut hasher, intervals.len() as u64);
    for (interval_index, interval) in intervals.iter().enumerate() {
        if interval_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        hash_u64(&mut hasher, interval.first_source_sample_index as u64);
        hash_u64(&mut hasher, interval.last_source_sample_index as u64);
        for value in [
            interval.start_time_s,
            interval.end_time_s,
            interval.duration_s,
        ] {
            hash_f64(&mut hasher, value);
        }
        match interval.visual_coverage.start_visualization_index {
            None => hasher.update(&[0]),
            Some(index) => {
                hasher.update(&[1]);
                hash_u64(&mut hasher, index as u64);
            }
        }
        hash_u64(
            &mut hasher,
            interval.visual_coverage.end_visualization_index as u64,
        );
        for stems in [
            interval.mean_force_stems_n,
            interval.force_time_stems_n_s,
            interval.measure_residual_stems_n_s,
        ] {
            for values in [stems.contact, stems.rolling, stems.base, stems.gas] {
                hash_components(&mut hasher, values);
            }
        }
        hash_components(&mut hasher, interval.mean_generalized_force_n);
        hash_components(&mut hasher, interval.generalized_force_time_n_s);
        hash_components(&mut hasher, interval.measure_residual_n_s);
        for availability in [
            interval.availability.contact,
            interval.availability.rolling,
            interval.availability.base,
            interval.availability.gas,
        ] {
            hasher.update(&[availability_tag(availability)]);
        }
        hasher.update(&[u8::from(interval.event_barrier)]);
        hash_u64(&mut hasher, interval.spatial_envelopes.len() as u64);
        for (mode_index, envelope) in interval.spatial_envelopes.iter().enumerate() {
            if mode_index % AUDIO_RESAMPLING_CANCELLATION_POLL_MODES == 0 {
                checkpoint_fn()?;
            }
            hash_u32(&mut hasher, envelope.mode_id);
            hash_f64(&mut hasher, envelope.start_factor);
            hash_f64(&mut hasher, envelope.end_factor);
            hasher.update(&[match envelope.source {
                crate::SpatialEnvelopeSource::DeclaredStatic => 1,
                crate::SpatialEnvelopeSource::ExactEndpointInterpolation => 2,
                crate::SpatialEnvelopeSource::HeldStartEndpoint => 3,
                crate::SpatialEnvelopeSource::HeldEndEndpoint => 4,
                crate::SpatialEnvelopeSource::MissingContactStatic => 5,
            }]);
        }
        match interval.artistic_texture {
            None => hasher.update(&[0]),
            Some(texture) => {
                hasher.update(&[1]);
                hasher.update(texture.stream_identity.as_bytes());
                hash_components(&mut hasher, texture.peak_force_envelope_n);
                hash_f64(&mut hasher, texture.band_low_hz);
                hash_f64(&mut hasher, texture.band_high_hz);
                hasher.update(&[availability_tag(texture.rolling_availability)]);
            }
        }
        hash_u64(&mut hasher, interval.events.len() as u64);
        for (event_index, event) in interval.events.iter().enumerate() {
            if event_index % AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            hash_u64(&mut hasher, event.source_sample_index as u64);
            hasher.update(&[match event.kind {
                ContactTransitionKind::Opening => 1,
                ContactTransitionKind::Reimpact => 2,
            }]);
            hash_f64(&mut hasher, event.time_s);
            hash_f64(&mut hasher, event.bracket_start_s);
            hash_f64(&mut hasher, event.bracket_end_s);
            hasher.update(&[match event.measure {
                ContactEventMeasure::TimingOnly => 1,
            }]);
            hash_components(&mut hasher, event.physical_impulse_n_s);
            match event.artistic {
                None => hasher.update(&[0]),
                Some(artistic) => {
                    hasher.update(&[1]);
                    hasher.update(artistic.stream_identity.as_bytes());
                    hash_components(&mut hasher, artistic.impulse_n_s);
                }
            }
        }
    }
    Ok(hasher.finalize())
}

fn availability_tag(availability: crate::ExcitationSourceAvailability) -> u8 {
    match availability {
        crate::ExcitationSourceAvailability::NotMapped => 1,
        crate::ExcitationSourceAvailability::Available => 2,
        crate::ExcitationSourceAvailability::Unavailable => 3,
    }
}

fn hash_components(hasher: &mut DomainHasher, values: ModalComponentValues) {
    hash_f64(hasher, values.disc);
    hash_f64(hasher, values.glass_plate);
    hash_f64(hasher, values.base_assembly);
}

fn hash_u32(hasher: &mut DomainHasher, value: u32) {
    hasher.update(&value.to_le_bytes());
}

fn hash_u64(hasher: &mut DomainHasher, value: u64) {
    hasher.update(&value.to_le_bytes());
}

fn hash_f64(hasher: &mut DomainHasher, value: f64) {
    hasher.update(&value.to_bits().to_le_bytes());
}

fn filter_identity(
    label: &[u8],
    coefficients: &[f64],
    lower_edge_hz: f64,
    upper_edge_hz: f64,
    diagnostics: AudioReconstructionFilterDiagnostics,
) -> ContentHash {
    let mut hasher = DomainHasher::new(FILTER_IDENTITY_DOMAIN);
    hash_u32(&mut hasher, AUDIO_RECONSTRUCTION_FILTER_VERSION);
    hash_u32(&mut hasher, STRICT_CORE_SEMANTICS_VERSION);
    hash_u64(&mut hasher, STRICT_CORE_GOLDEN_HASH);
    hasher.update(label);
    hash_f64(&mut hasher, lower_edge_hz);
    hash_f64(&mut hasher, upper_edge_hz);
    hash_u64(&mut hasher, diagnostics.tap_count as u64);
    hash_f64(&mut hasher, diagnostics.measured_passband_ripple_db);
    hash_f64(&mut hasher, diagnostics.measured_stopband_attenuation_db);
    hash_u32(&mut hasher, diagnostics.intrinsic_group_delay_frames);
    hash_u32(&mut hasher, diagnostics.group_delay_compensation_frames);
    hash_u32(&mut hasher, diagnostics.required_lookahead_frames);
    hasher.update(&diagnostics.published_alignment_offset_frames.to_le_bytes());
    for coefficient in coefficients {
        hash_f64(&mut hasher, *coefficient);
    }
    hasher.finalize()
}

fn resampler_identity(
    excitation_identity: ContentHash,
    modal_identity: ContentHash,
    source_payload_identity: ContentHash,
    filter_identity: ContentHash,
    artistic_filter_identity: Option<ContentHash>,
    input: AudioResamplingModelInput,
    total_audio_frames: u64,
) -> ContentHash {
    let mut hasher = DomainHasher::new(RESAMPLER_IDENTITY_DOMAIN);
    hash_u32(&mut hasher, AUDIO_RESAMPLING_ALGORITHM_VERSION);
    hash_u32(&mut hasher, AUDIO_RECONSTRUCTION_FILTER_VERSION);
    hash_u32(&mut hasher, STRICT_CORE_SEMANTICS_VERSION);
    hash_u64(&mut hasher, STRICT_CORE_GOLDEN_HASH);
    hasher.update(excitation_identity.as_bytes());
    hasher.update(modal_identity.as_bytes());
    hasher.update(source_payload_identity.as_bytes());
    hasher.update(filter_identity.as_bytes());
    match artistic_filter_identity {
        None => hasher.update(&[0]),
        Some(identity) => {
            hasher.update(&[1]);
            hasher.update(identity.as_bytes());
        }
    }
    hash_clock(&mut hasher, input.video_clock);
    hash_clock(&mut hasher, input.audio_clock);
    hash_f64(&mut hasher, input.declared_source_bandwidth_hz);
    hash_f64(&mut hasher, input.filter.passband_edge_hz);
    hash_f64(&mut hasher, input.filter.stopband_edge_hz);
    hasher.update(&input.filter.half_length.to_le_bytes());
    hash_f64(&mut hasher, input.filter.maximum_passband_ripple_db);
    hash_f64(&mut hasher, input.filter.minimum_stopband_attenuation_db);
    hash_u32(&mut hasher, input.filter.response_grid_intervals);
    hasher.update(&[match input.boundary_policy {
        AudioResamplingBoundaryPolicy::HalfSampleEvenReflectionV1 => 1,
    }]);
    hasher.update(&[match input.event_fractional_delay {
        AudioEventFractionalDelay::LinearTwoBoundaryV1 => 1,
    }]);
    for value in [
        input.budget.maximum_source_intervals as u64,
        input.budget.maximum_total_audio_frames,
        input.budget.maximum_chunk_audio_frames as u64,
        input.budget.maximum_filter_taps as u64,
        input.budget.maximum_chunk_mode_values as u64,
        input.budget.maximum_events as u64,
        input.budget.maximum_sync_markers as u64,
        input.budget.maximum_chunk_multiply_adds,
        total_audio_frames,
    ] {
        hash_u64(&mut hasher, value);
    }
    hasher.finalize()
}

fn chunk_identity(model_identity: ContentHash, start: u64, end: u64) -> ContentHash {
    let mut hasher = DomainHasher::new(CHUNK_IDENTITY_DOMAIN);
    hasher.update(model_identity.as_bytes());
    hash_u64(&mut hasher, start);
    hash_u64(&mut hasher, end);
    hasher.finalize()
}

fn hash_clock(hasher: &mut DomainHasher, clock: CinematicClock) {
    hasher.update(&[match clock.domain() {
        CinematicClockDomain::Simulation => 1,
        CinematicClockDomain::Video => 2,
        CinematicClockDomain::Audio => 3,
        CinematicClockDomain::Composition => 4,
        CinematicClockDomain::Timeless => 5,
    }]);
    hash_u32(hasher, clock.ticks_per_second_numerator());
    hash_u32(hasher, clock.ticks_per_second_denominator());
    hasher.update(&clock.start_tick().to_le_bytes());
    hasher.update(&clock.end_tick_exclusive().to_le_bytes());
}

fn reserve_exact<T>(
    values: &mut Vec<T>,
    requested: usize,
    artifact: &'static str,
) -> Result<(), AudioResamplingError> {
    values
        .try_reserve_exact(requested)
        .map_err(|_| AudioResamplingError::Capacity {
            artifact,
            requested,
        })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), AudioResamplingError> {
    cx.checkpoint().map_err(|_| AudioResamplingError::Cancelled)
}

#[derive(Debug, Clone, Copy)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    const fn new() -> Self {
        Self {
            sum: 0.0,
            correction: 0.0,
        }
    }

    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn total(self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn g0_half_sample_even_reflection_handles_one_and_two_frame_horizons() {
        for index in -32..=32 {
            assert_eq!(reflect_half_sample_even(index, 1), 0);
        }

        let expected = [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0];
        for (index, expected) in (-6..=7).zip(expected) {
            assert_eq!(reflect_half_sample_even(index, 2), expected);
        }
        assert_eq!(reflected_source_span(0, 1, 8, 1), (0, 0));
        assert_eq!(reflected_source_span(0, 2, 5, 2), (0, 1));

        let raw = RawAudioSpan {
            first_frame_offset: 0,
            mode_count: 1,
            distributed_force_n: vec![
                ModalComponentValues {
                    disc: 10.0,
                    glass_plate: 0.0,
                    base_assembly: 0.0,
                },
                ModalComponentValues {
                    disc: 20.0,
                    glass_plate: 0.0,
                    base_assembly: 0.0,
                },
            ],
            physical_localized_mode_force_n: vec![1.0, 3.0],
            artistic_localized_mode_force_n: vec![0.0, 0.0],
        };
        let coefficients = [0.25, 0.5, 0.25];
        assert_close(
            convolve_components(0, &raw, &coefficients, 2).disc,
            12.5,
            f64::EPSILON,
        );
        assert_close(
            convolve_components(1, &raw, &coefficients, 2).disc,
            17.5,
            f64::EPSILON,
        );
        assert_close(
            convolve_mode(0, 0, &raw, &coefficients, 2, false),
            1.5,
            f64::EPSILON,
        );
        assert_close(
            convolve_mode(1, 0, &raw, &coefficients, 2, false),
            2.5,
            f64::EPSILON,
        );
    }

    #[test]
    fn g0_linear_factor_integral_is_exactly_additive_for_varying_factors() {
        let whole = linear_factor_integral(100.0, 4.0, 100.0, 104.0, -2.0, 6.0);
        let left = linear_factor_integral(100.0, 4.0, 100.0, 101.0, -2.0, 6.0);
        let middle = linear_factor_integral(100.0, 4.0, 101.0, 103.0, -2.0, 6.0);
        let right = linear_factor_integral(100.0, 4.0, 103.0, 104.0, -2.0, 6.0);

        assert_eq!(whole, 8.0);
        assert_eq!(middle, 4.0);
        assert_eq!(left + middle + right, whole);
    }

    #[test]
    fn g5_exact_long_master_clocks_have_zero_marker_drift_at_large_epochs() {
        const VIDEO_FRAMES: i64 = 8_192;
        const VIDEO_START: i64 = 1_000_000_000;
        const AUDIO_PER_VIDEO: i64 = 2_000;

        let video = CinematicClock::try_new(
            CinematicClockDomain::Video,
            SOUND_MASTER_VIDEO_RATE_HZ,
            1,
            VIDEO_START,
            VIDEO_START + VIDEO_FRAMES,
        )
        .unwrap();
        let audio = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            VIDEO_START * AUDIO_PER_VIDEO,
            (VIDEO_START + VIDEO_FRAMES) * AUDIO_PER_VIDEO,
        )
        .unwrap();
        let mut budget = AudioResamplingBudget::reference_film();
        budget.maximum_total_audio_frames = (VIDEO_FRAMES * AUDIO_PER_VIDEO) as u64;
        budget.maximum_sync_markers = VIDEO_FRAMES as usize + 1;

        let mut polls = 0_usize;
        let (total, alignment, period_s) = validate_clocks(video, audio, budget, &mut || {
            polls += 1;
            Ok(())
        })
        .unwrap();

        assert_eq!(total, (VIDEO_FRAMES * AUDIO_PER_VIDEO) as u64);
        assert_eq!(alignment.audio_frames_per_video_frame, 2_000);
        assert_eq!(alignment.endpoint_drift_audio_frames, 0);
        assert_eq!(alignment.markers.len(), VIDEO_FRAMES as usize + 1);
        for (index, marker) in alignment.markers.iter().enumerate() {
            assert_eq!(marker.video_tick, VIDEO_START + index as i64);
            assert_eq!(
                marker.audio_tick,
                (VIDEO_START + index as i64) * AUDIO_PER_VIDEO
            );
            assert_eq!(marker.audio_frame_offset, index as u64 * 2_000);
        }
        assert_eq!(period_s, 1.0 / 48_000.0);
        assert!(
            polls > 1,
            "long marker generation must poll within the loop"
        );
    }

    #[test]
    fn g4_windowed_filter_design_observes_in_work_cancellation() {
        let mut polls = 0_usize;
        let error = design_windowed_lowpass(2_000.0, 200, &mut || {
            polls += 1;
            if polls == 3 {
                Err(AudioResamplingError::Cancelled)
            } else {
                Ok(())
            }
        })
        .unwrap_err();

        assert_eq!(error, AudioResamplingError::Cancelled);
        assert_eq!(polls, 3);
    }

    #[test]
    fn g0_centered_response_matches_independent_direct_dft_spots() {
        let coefficients = [0.03, 0.17, 0.60, 0.17, 0.03];
        let radius = (coefficients.len() / 2) as isize;
        for frequency_hz in [0.0, 1_234.5, 12_000.0, 23_999.0] {
            let mut direct_real = 0.0;
            let mut direct_imaginary = 0.0;
            for (tap, coefficient) in coefficients.iter().copied().enumerate() {
                let centered_offset = tap as isize - radius;
                let phase = -2.0 * PI * frequency_hz * centered_offset as f64
                    / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
                direct_real += coefficient * det::cos(phase);
                direct_imaginary += coefficient * det::sin(phase);
            }

            assert_close(direct_imaginary, 0.0, 2.0e-15);
            assert_close(
                centered_frequency_response(&coefficients, frequency_hz),
                direct_real,
                2.0e-15,
            );
        }
    }
}
