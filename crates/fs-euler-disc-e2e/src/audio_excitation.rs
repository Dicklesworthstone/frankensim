//! Source-clock mechanics-to-modal excitation mapping for the Euler-disc film.
//!
//! This module converts admitted interval measures into component-attributed
//! generalized-force measures. It deliberately does **not** sample the result at
//! 48 kHz: rational-clock resampling, anti-alias FIR design, filter latency, and
//! fractional-delay event placement belong to the following multirate stage.
//! Current contact events are timing-only, so their physical impulse is always
//! exactly zero. Optional procedural controls are separately typed as artistic.

use core::{fmt, num::NonZeroUsize};

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::cinematic_sound::{
    SOUND_MASTER_SAMPLE_RATE_HZ, SoundExcitationChannel, SoundExcitationControl,
    SoundModalComponent, SoundSynthesisConfig, SoundTrajectoryDisposition,
};
use fs_exec::Cx;
use fs_math::{STRICT_CORE_GOLDEN_HASH, STRICT_CORE_SEMANTICS_VERSION, det};
use fs_mbd::Vec3;
use fs_rand::{STREAM_SEMANTICS_VERSION, Stream, StreamKey as RandomStreamKey};

use crate::{
    AudioControlFilter, AudioVisualCoverage, ChannelControl, ChannelControlSet,
    CoarsenedAudioControls, ContactEventMeasure, ControlContactEvent,
    EULER_CONTROL_STREAM_SCHEMA_VERSION, EulerControlStream, EulerRenderTrajectoryArtifact,
    MODAL_SYNTHESIS_ALGORITHM_VERSION, ModalComponentValues, ModalSynthesisModel,
    RenderSampleDisposition, coupled_runner::ContactTransitionKind,
};

/// Version of the source-clock mapping, checkpoint, and identity semantics.
pub const AUDIO_EXCITATION_ALGORITHM_VERSION: u32 = 2;
/// Largest number of source intervals admitted by one transactional mapping call.
pub const MAX_AUDIO_EXCITATION_CHUNK_INTERVALS: usize = 65_536;
/// Largest supported azimuthal harmonic in the compact contact-shape model.
pub const MAX_AUDIO_EXCITATION_AZIMUTHAL_HARMONIC: u16 = 64;
/// Maximum delay between cancellation polls while mapping source intervals.
pub const AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS: usize = 64;

const MAPPER_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-excitation-mapper.v2";
const CHUNK_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-excitation-chunk.v2";
const TEXTURE_STREAM_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-texture-stream.v1";
const EVENT_STREAM_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-event-stream.v1";

/// Source-cadence reduction performed before excitation mapping.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioExcitationReduction {
    /// Preserve every exact accepted source interval.
    RawIntervals,
    /// Reuse the control stream's measure-first, event-barrier boxcar.
    WholeIntervalBoxcarV1 {
        /// Maximum event-free source intervals combined into one bin.
        intervals_per_bin: NonZeroUsize,
    },
}

/// Honest reconstruction status attached to every v2 mapped timeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioExcitationReconstructionStatus {
    /// Interval averages and measures are available, but the result is not an
    /// audio-rate or mathematically band-limited signal.
    RequiresBandLimitedResampling,
}

/// Whether one mechanics-derived source family was requested and retained.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExcitationSourceAvailability {
    /// No mapping selected this source family.
    NotMapped,
    /// The upstream control channel was explicitly present, including when zero.
    Available,
    /// A mapping requested the family, but the upstream channel was unavailable.
    Unavailable,
}

/// Availability of all source-attributed mechanics-derived stems.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AudioExcitationAvailability {
    /// Aggregate contact controls.
    pub contact: ExcitationSourceAvailability,
    /// Reduced rolling-resistance controls.
    pub rolling: ExcitationSourceAvailability,
    /// Reduced-base damping controls.
    pub base: ExcitationSourceAvailability,
    /// Exterior-gas body-work controls.
    pub gas: ExcitationSourceAvailability,
}

/// Component coordinates separated by their mechanics source family.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AudioExcitationStems {
    /// Aggregate contact contribution.
    pub contact: ModalComponentValues,
    /// Reduced rolling-resistance contribution.
    pub rolling: ModalComponentValues,
    /// Reduced-base damping contribution.
    pub base: ModalComponentValues,
    /// Exterior-gas body-work contribution.
    pub gas: ModalComponentValues,
}

impl AudioExcitationStems {
    /// Exact zero across every mechanics-derived source stem.
    pub const ZERO: Self = Self {
        contact: ModalComponentValues::ZERO,
        rolling: ModalComponentValues::ZERO,
        base: ModalComponentValues::ZERO,
        gas: ModalComponentValues::ZERO,
    };

    /// Deterministic contact, rolling, base, then gas sum.
    #[must_use]
    pub fn sum(self) -> ModalComponentValues {
        let mut sums = [CompensatedSum::new(); 3];
        for values in [self.contact, self.rolling, self.base, self.gas] {
            sums[0].add(values.disc);
            sums[1].add(values.glass_plate);
            sums[2].add(values.base_assembly);
        }
        ModalComponentValues {
            disc: sums[0].total(),
            glass_plate: sums[1].total(),
            base_assembly: sums[2].total(),
        }
    }

    /// Contact-location-dependent contribution. Only this drive may be
    /// multiplied by [`ModalSpatialEnvelope`] factors downstream.
    #[must_use]
    pub fn localized_sum(self) -> ModalComponentValues {
        component_add(self.contact, self.rolling)
    }

    /// Source-location-independent base/gas contribution. This drive retains
    /// the modal model's declared static participation downstream.
    #[must_use]
    pub fn distributed_sum(self) -> ModalComponentValues {
        component_add(self.base, self.gas)
    }
}

/// Explicit resource and headroom limits. A violation refuses; values are never
/// silently saturated or normalized.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioExcitationBudget {
    /// Maximum raw source intervals admitted before optional reduction.
    pub maximum_total_intervals: usize,
    /// Maximum intervals accepted by one transactional mapping call.
    pub maximum_chunk_intervals: usize,
    /// Maximum retained event records in one chunk.
    pub maximum_chunk_events: usize,
    /// Maximum interval-by-mode spatial envelopes materialized in one chunk.
    pub maximum_chunk_spatial_envelopes: usize,
    /// Maximum absolute force-valued source scalar [N].
    pub maximum_abs_source_force_n: f64,
    /// Maximum absolute work-rate-valued source scalar [W].
    pub maximum_abs_source_work_rate_w: f64,
    /// Maximum absolute mapped component force, per stem or summed [N].
    pub maximum_abs_generalized_force_n: f64,
    /// Maximum absolute mapped component force-time measure [N s].
    pub maximum_abs_force_time_measure_n_s: f64,
    /// Maximum admitted arithmetic `mean * duration - measure` residual [N s].
    pub maximum_measure_residual_n_s: f64,
    /// Maximum artistic rolling-texture force envelope [N].
    pub maximum_artistic_texture_envelope_n: f64,
    /// Maximum artistic event impulse magnitude [N s].
    pub maximum_artistic_event_impulse_n_s: f64,
}

impl AudioExcitationBudget {
    /// Bounded defaults suitable for one admitted film trajectory.
    #[must_use]
    pub const fn reference_film(maximum_total_intervals: usize) -> Self {
        let maximum_chunk_intervals =
            if maximum_total_intervals < MAX_AUDIO_EXCITATION_CHUNK_INTERVALS {
                maximum_total_intervals
            } else {
                MAX_AUDIO_EXCITATION_CHUNK_INTERVALS
            };
        Self {
            maximum_total_intervals,
            maximum_chunk_intervals,
            maximum_chunk_events: MAX_AUDIO_EXCITATION_CHUNK_INTERVALS,
            maximum_chunk_spatial_envelopes: 1_048_576,
            maximum_abs_source_force_n: 1.0e9,
            maximum_abs_source_work_rate_w: 1.0e12,
            maximum_abs_generalized_force_n: 1.0e9,
            maximum_abs_force_time_measure_n_s: 1.0e12,
            maximum_measure_residual_n_s: 1.0e-6,
            maximum_artistic_texture_envelope_n: 1.0e6,
            maximum_artistic_event_impulse_n_s: 1.0e3,
        }
    }
}

/// Compact declared source-location shape for one modal coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContactModeShape {
    /// Use the mode's declared static participation.
    Uniform,
    /// Signed cosine around the local disc/body or base azimuth.
    AzimuthalCosine {
        /// Nonzero integer nodal-diameter harmonic.
        harmonic: u16,
        /// Declared angular phase [rad].
        phase_rad: f64,
    },
}

/// Contact-location rule keyed by stable modal ID.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModeContactParticipationRule {
    /// Stable mode identifier from the bound modal model.
    pub mode_id: u32,
    /// Declared compact shape.
    pub shape: ContactModeShape,
}

/// Whether source-location participation is static or explicitly reconstructed.
#[derive(Debug, Clone, PartialEq)]
pub enum ContactParticipationPolicy {
    /// Emit exact factor one for every mode.
    DeclaredStatic,
    /// Require exactly one canonicalizable rule for every bound mode.
    ContactCoordinates {
        /// Rules may be supplied in any order and are canonicalized by mode ID.
        rules: Vec<ModeContactParticipationRule>,
    },
}

/// Provenance of one interval's spatial envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpatialEnvelopeSource {
    /// Static participation was explicitly selected.
    DeclaredStatic,
    /// Both exact closed endpoint contact coordinates were evaluated.
    ExactEndpointInterpolation,
    /// Only the exact opening endpoint contact coordinate was available.
    HeldStartEndpoint,
    /// Only the exact closing endpoint contact coordinate was available.
    HeldEndEndpoint,
    /// No exact closed endpoint contact coordinate was available.
    MissingContactStatic,
}

/// Start/end participation for one canonical mode over one source interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalSpatialEnvelope {
    /// Stable canonical mode ID.
    pub mode_id: u32,
    /// Participation multiplier at the interval start.
    pub start_factor: f64,
    /// Participation multiplier at the interval end.
    pub end_factor: f64,
    /// How the factors were obtained.
    pub source: SpatialEnvelopeSource,
}

/// Explicitly artistic procedural texture controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArtisticTextureConfig {
    /// User-selected deterministic seed.
    pub seed: u64,
    /// Peak rolling-noise envelope force per absolute rolling work rate [N/W].
    pub rolling_force_gain_n_per_w: f64,
    /// Component receiving rolling microtexture.
    pub rolling_target_component: SoundModalComponent,
    /// Lower declared texture band edge for the later resampler [Hz].
    pub band_low_hz: f64,
    /// Upper declared texture band edge for the later resampler [Hz].
    pub band_high_hz: f64,
    /// Maximum signed artistic impulse at a reimpact marker [N s].
    pub reimpact_impulse_n_s: f64,
    /// Component receiving the artistic reimpact impulse.
    pub reimpact_target_component: SoundModalComponent,
}

/// Per-interval envelope and random-access stream identity for later procedural
/// band-limited noise generation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArtisticTextureEnvelope {
    /// Domain-separated stream identity bound to source, seed, and interval.
    pub stream_identity: ContentHash,
    /// Peak component force envelope [N]. This is not a physical force sample.
    pub peak_force_envelope_n: ModalComponentValues,
    /// Lower declared texture band edge [Hz].
    pub band_low_hz: f64,
    /// Upper declared texture band edge [Hz].
    pub band_high_hz: f64,
    /// Whether the rolling channel was retained rather than merely numerical zero.
    pub rolling_availability: ExcitationSourceAvailability,
}

/// Optional artistic contribution attached to one timing-only reimpact marker.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArtisticEventExcitation {
    /// Domain-separated stream identity bound to exact event metadata.
    pub stream_identity: ContentHash,
    /// Declared artistic impulse [N s], kept separate from physical impulse.
    pub impulse_n_s: ModalComponentValues,
}

/// Complete mapper input. Every field is bound into mapper identity.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioExcitationModelInput {
    /// Canonical sound-config source mappings.
    pub mappings: Vec<SoundExcitationControl>,
    /// Optional source-cadence measure-first reduction.
    pub reduction: AudioExcitationReduction,
    /// Explicit location-participation policy.
    pub spatial_policy: ContactParticipationPolicy,
    /// Optional artistic texture controls.
    pub artistic_texture: Option<ArtisticTextureConfig>,
    /// Resource and headroom limits.
    pub budget: AudioExcitationBudget,
}

/// Exact selected source grid and its deliberately limited interpretation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioExcitationGrid {
    /// Number of selected raw intervals or coarsened bins.
    pub interval_count: usize,
    /// Exact first selected start [s].
    pub start_time_s: f64,
    /// Exact final selected end [s].
    pub end_time_s: f64,
    /// Smallest selected positive duration [s].
    pub minimum_interval_duration_s: f64,
    /// Largest selected positive duration [s].
    pub maximum_interval_duration_s: f64,
    /// Producer-declared nominal mechanics timestep [s].
    pub nominal_mechanics_timestep_s: f64,
    /// `0.5 / maximum_interval_duration_s`, exposed only as a conservative
    /// nominal ceiling and not as a stable-reconstruction proof.
    pub nominal_source_nyquist_ceiling_hz: f64,
    /// Required downstream treatment.
    pub reconstruction: AudioExcitationReconstructionStatus,
}

/// One retained contact event with an explicit zero physical impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AudioExcitationEvent {
    /// Source interval that owns this event.
    pub source_sample_index: usize,
    /// Opening or reimpact.
    pub kind: ContactTransitionKind,
    /// Exact localized time [s].
    pub time_s: f64,
    /// Inclusive bracket start [s].
    pub bracket_start_s: f64,
    /// Inclusive bracket end [s].
    pub bracket_end_s: f64,
    /// Explicit timing-only upstream measure.
    pub measure: ContactEventMeasure,
    /// Always zero in v2 because the source retains no event impulse.
    pub physical_impulse_n_s: ModalComponentValues,
    /// Optional separately declared artistic rendering contribution.
    pub artistic: Option<ArtisticEventExcitation>,
}

/// One mapped source-clock interval. Mean generalized forces and separately
/// retained force-time measures are both exposed so downstream resampling starts
/// from measures. The residual is an arithmetic self-consistency check, not
/// independent physical evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioExcitationInterval {
    /// First source interval represented by this output.
    pub first_source_sample_index: usize,
    /// Last source interval represented by this output, inclusive.
    pub last_source_sample_index: usize,
    /// Exact interval start [s].
    pub start_time_s: f64,
    /// Exact interval end [s].
    pub end_time_s: f64,
    /// Exact positive duration [s].
    pub duration_s: f64,
    /// Exact visualization endpoint coverage.
    pub visual_coverage: AudioVisualCoverage,
    /// Per-source duration-mean generalized force [N].
    pub mean_force_stems_n: AudioExcitationStems,
    /// Per-source generalized force-time measure [N s].
    pub force_time_stems_n_s: AudioExcitationStems,
    /// Deterministic sum of the transfer-mapped mean-force stems [N].
    pub mean_generalized_force_n: ModalComponentValues,
    /// Deterministic sum of the transfer-mapped force-time stems [N s].
    pub generalized_force_time_n_s: ModalComponentValues,
    /// Requested/retained state of each mechanics-derived source family.
    pub availability: AudioExcitationAvailability,
    /// Ordered timing-only event markers.
    pub events: Vec<AudioExcitationEvent>,
    /// True when the source reduction isolated an event barrier.
    pub event_barrier: bool,
    /// Canonical per-mode source-location envelope. These factors apply only
    /// to the localized contact/rolling stems, never the base/gas stems.
    pub spatial_envelopes: Vec<ModalSpatialEnvelope>,
    /// Optional separately labeled artistic rolling envelope.
    pub artistic_texture: Option<ArtisticTextureEnvelope>,
    /// Per-source `mean * duration - retained measure` [N s].
    pub measure_residual_stems_n_s: AudioExcitationStems,
    /// Sum of per-source arithmetic residuals for each component [N s].
    pub measure_residual_n_s: ModalComponentValues,
}

impl AudioExcitationInterval {
    /// Contact-location-dependent source-cadence mean generalized force [N].
    /// Downstream spatial reconstruction applies modal factors to this value
    /// only; it must preserve the corresponding interval measure.
    #[must_use]
    pub fn localized_mean_generalized_force_n(&self) -> ModalComponentValues {
        self.mean_force_stems_n.localized_sum()
    }

    /// Base/gas source-cadence mean generalized force [N]. It must retain
    /// declared static modal participation and must not be multiplied by
    /// contact-location factors.
    #[must_use]
    pub fn distributed_mean_generalized_force_n(&self) -> ModalComponentValues {
        self.mean_force_stems_n.distributed_sum()
    }

    /// Contact-location-dependent generalized force-time measure [N s].
    /// This is the authority-bearing quantity for downstream resampling.
    #[must_use]
    pub fn localized_force_time_measure_n_s(&self) -> ModalComponentValues {
        self.force_time_stems_n_s.localized_sum()
    }

    /// Base/gas generalized force-time measure [N s]. This is reconstructed
    /// without contact-location modulation downstream.
    #[must_use]
    pub fn distributed_force_time_measure_n_s(&self) -> ModalComponentValues {
        self.force_time_stems_n_s.distributed_sum()
    }

    /// Linearly evaluate the already-derived per-mode envelope at `alpha` in
    /// `[0,1]`. This is not audio-rate resampling.
    pub fn spatial_factors_at(&self, alpha: f64) -> Result<Vec<f64>, AudioExcitationError> {
        if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) {
            return Err(AudioExcitationError::InvalidSpatialEvaluation);
        }
        let mut factors = Vec::new();
        factors
            .try_reserve_exact(self.spatial_envelopes.len())
            .map_err(|_| AudioExcitationError::Capacity {
                artifact: "evaluated spatial factors",
                requested: self.spatial_envelopes.len(),
            })?;
        for envelope in &self.spatial_envelopes {
            factors.push(
                (envelope.end_factor - envelope.start_factor).mul_add(alpha, envelope.start_factor),
            );
        }
        Ok(factors)
    }
}

/// Chunk-local and cumulative mapping diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioExcitationDiagnostics {
    /// First selected interval index in this chunk.
    pub start_interval_index: usize,
    /// Exclusive selected interval index after this chunk.
    pub end_interval_index: usize,
    /// Timing-only events retained in this chunk.
    pub event_count: usize,
    /// Largest absolute interval measure-reconciliation residual [N s].
    pub maximum_abs_measure_residual_n_s: f64,
    /// Chunk-local summed transfer-mapped force-time stems [N s].
    pub chunk_force_time_stems_n_s: AudioExcitationStems,
    /// Cumulative summed transfer-mapped force-time stems [N s].
    pub cumulative_force_time_stems_n_s: AudioExcitationStems,
}

/// Immutable restart point for source-interval mapping.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioExcitationCheckpoint {
    mapper_identity: ContentHash,
    source_identity: ContentHash,
    next_interval_index: usize,
    last_end_time_bits: Option<u64>,
    cumulative_force_time: StemAccumulator,
}

impl AudioExcitationCheckpoint {
    /// Mapper identity to which the checkpoint is bound.
    #[must_use]
    pub const fn mapper_identity(&self) -> ContentHash {
        self.mapper_identity
    }

    /// Exact durable source-trajectory identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Next selected interval index.
    #[must_use]
    pub const fn next_interval_index(&self) -> usize {
        self.next_interval_index
    }

    /// Cumulative force-time stems represented before the next interval.
    #[must_use]
    pub fn cumulative_force_time_stems_n_s(&self) -> AudioExcitationStems {
        self.cumulative_force_time.total()
    }
}

/// Atomically published mapped source intervals and successor checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct AudioExcitationChunk {
    /// Domain-separated identity of mapper plus exact selected range.
    pub identity: ContentHash,
    /// Mapped intervals in strict source-clock order.
    pub intervals: Vec<AudioExcitationInterval>,
    /// Mapping diagnostics.
    pub diagnostics: AudioExcitationDiagnostics,
    /// Restart point for the following interval.
    pub successor: AudioExcitationCheckpoint,
}

/// Typed refusal from source-clock excitation admission or mapping.
#[derive(Debug, Clone, PartialEq)]
pub enum AudioExcitationError {
    /// Execution cancellation was observed before atomic publication.
    Cancelled,
    /// Controls are not pointer-bound to the trajectory inside the named artifact.
    SourceBindingMismatch,
    /// The durable source ended in numerical refusal and cannot mint excitation.
    NumericalRefusalSource,
    /// No positive-duration source interval is available.
    EmptySourceIntervals,
    /// Mapping controls were absent, unordered, duplicated, or malformed.
    InvalidMappings,
    /// The current control schema cannot honestly supply this selector.
    UnsupportedMapping(SoundExcitationChannel),
    /// Upstream measure-first source reduction refused with its exact reason.
    SourceReduction(crate::ControlStreamError),
    /// A named budget is zero, non-finite, negative, or above a hard ceiling.
    InvalidBudget(&'static str),
    /// Artistic texture configuration is malformed.
    InvalidArtisticTexture(&'static str),
    /// Spatial rules are incomplete, duplicated, or malformed.
    InvalidSpatialPolicy(&'static str),
    /// A checkpoint does not belong to this exact mapper/source/range.
    InvalidCheckpoint,
    /// Every selected source interval has already been consumed.
    Complete,
    /// The requested source chunk exceeds its explicit interval budget.
    ChunkIntervalBudgetExceeded { requested: usize, limit: usize },
    /// The selected chunk contains more events than admitted.
    ChunkEventBudgetExceeded { requested: usize, limit: usize },
    /// The selected chunk would materialize more modal envelopes than admitted.
    ChunkSpatialEnvelopeBudgetExceeded { requested: usize, limit: usize },
    /// Allocator refusal for an explicitly preflighted result collection.
    Capacity {
        artifact: &'static str,
        requested: usize,
    },
    /// An upstream or mapped scalar was non-finite.
    NonFinite {
        interval: usize,
        field: &'static str,
    },
    /// A source scalar or mapped result exceeded an explicit limit.
    LimitExceeded {
        interval: usize,
        field: &'static str,
        magnitude: f64,
        limit: f64,
    },
    /// Event records were not strict, bracketed, timing-only, and in range.
    InvalidEventOrder { interval: usize, event: usize },
    /// Selected source intervals were not finite, positive, or contiguous.
    InvalidSourceGrid {
        /// Selected interval index.
        interval: usize,
        /// Stable semantic field.
        field: &'static str,
    },
    /// A requested spatial evaluation coordinate was outside `[0,1]`.
    InvalidSpatialEvaluation,
    /// The admitted high-level sound configuration disagrees with the mapper.
    SoundConfigurationMismatch(&'static str),
}

impl fmt::Display for AudioExcitationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for AudioExcitationError {}

#[derive(Debug, Clone, Copy)]
struct PreparedModeRule {
    mode_id: u32,
    component: SoundModalComponent,
    shape: ContactModeShape,
}

enum SelectedIntervals<'trajectory> {
    Raw,
    Coarsened(CoarsenedAudioControls<'trajectory>),
}

/// Admitted, content-identified source-clock mapper.
pub struct AudioExcitationMapper<'controls, 'trajectory> {
    source_identity: ContentHash,
    identity: ContentHash,
    controls: &'controls EulerControlStream<'trajectory>,
    selected: SelectedIntervals<'trajectory>,
    mappings: Vec<SoundExcitationControl>,
    reduction: AudioExcitationReduction,
    mode_rules: Vec<PreparedModeRule>,
    modal_identity: ContentHash,
    modal_modes: Vec<fs_evidence::cinematic_sound::SoundMode>,
    base_normal_world: Vec3,
    artistic_texture: Option<ArtisticTextureConfig>,
    texture_root_identity: Option<ContentHash>,
    budget: AudioExcitationBudget,
    grid: AudioExcitationGrid,
}

impl<'controls, 'trajectory> AudioExcitationMapper<'controls, 'trajectory> {
    /// Admit exact source, mapping, modal, spatial, and artistic identities.
    pub fn try_new(
        artifact: &'trajectory EulerRenderTrajectoryArtifact,
        controls: &'controls EulerControlStream<'trajectory>,
        modal: &ModalSynthesisModel,
        mut input: AudioExcitationModelInput,
        cx: &Cx<'_>,
    ) -> Result<Self, AudioExcitationError> {
        checkpoint(cx)?;
        if !controls.is_bound_to(artifact.trajectory()) {
            return Err(AudioExcitationError::SourceBindingMismatch);
        }
        if matches!(
            artifact
                .trajectory()
                .samples()
                .last()
                .map(|sample| sample.input().disposition),
            Some(RenderSampleDisposition::NumericalRefusal(_))
        ) {
            return Err(AudioExcitationError::NumericalRefusalSource);
        }
        validate_budget(input.budget)?;
        validate_mappings(&input.mappings)?;
        validate_texture(input.artistic_texture, input.budget)?;
        let raw_interval_count = controls.audio().len();
        if raw_interval_count == 0 {
            return Err(AudioExcitationError::EmptySourceIntervals);
        }
        if raw_interval_count > input.budget.maximum_total_intervals {
            return Err(AudioExcitationError::LimitExceeded {
                interval: raw_interval_count,
                field: "raw source intervals",
                magnitude: raw_interval_count as f64,
                limit: input.budget.maximum_total_intervals as f64,
            });
        }
        let has_time_varying_spatial_rules = match &input.spatial_policy {
            ContactParticipationPolicy::DeclaredStatic => false,
            ContactParticipationPolicy::ContactCoordinates { rules } => rules
                .iter()
                .any(|rule| matches!(rule.shape, ContactModeShape::AzimuthalCosine { .. })),
        };
        if matches!(
            input.reduction,
            AudioExcitationReduction::WholeIntervalBoxcarV1 { .. }
        ) && has_time_varying_spatial_rules
        {
            return Err(AudioExcitationError::InvalidSpatialPolicy(
                "source coarsening cannot preserve time-varying contact participation",
            ));
        }

        let selected = match input.reduction {
            AudioExcitationReduction::RawIntervals => SelectedIntervals::Raw,
            AudioExcitationReduction::WholeIntervalBoxcarV1 { intervals_per_bin } => {
                SelectedIntervals::Coarsened(
                    controls.boxcar_coarsen(intervals_per_bin, cx).map_err(
                        |error| match error {
                            crate::ControlStreamError::Cancelled => AudioExcitationError::Cancelled,
                            other => AudioExcitationError::SourceReduction(other),
                        },
                    )?,
                )
            }
        };
        let interval_count = selected_len(controls, &selected);
        debug_assert!(interval_count > 0);
        if interval_count > input.budget.maximum_total_intervals {
            return Err(AudioExcitationError::LimitExceeded {
                interval: interval_count,
                field: "total selected intervals",
                magnitude: interval_count as f64,
                limit: input.budget.maximum_total_intervals as f64,
            });
        }

        let mut mode_rules = prepare_mode_rules(modal, &mut input.spatial_policy)?;
        mode_rules.sort_by_key(|rule| rule.mode_id);
        let grid = derive_grid(controls, &selected, cx)?;
        let source_identity = artifact.receipt().artifact_identity();
        let base_normal_world = artifact
            .trajectory()
            .metadata()
            .base_frame
            .orientation_base_to_world
            .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        let texture_root_identity = input
            .artistic_texture
            .map(|texture| texture_root_identity(source_identity, texture));
        let identity = mapper_identity(
            source_identity,
            modal.identity(),
            &input.mappings,
            input.reduction,
            &mode_rules,
            input.artistic_texture,
            input.budget,
            grid,
        );
        checkpoint(cx)?;
        Ok(Self {
            source_identity,
            identity,
            controls,
            selected,
            mappings: input.mappings,
            reduction: input.reduction,
            mode_rules,
            modal_identity: modal.identity(),
            modal_modes: modal.modes().to_vec(),
            base_normal_world,
            artistic_texture: input.artistic_texture,
            texture_root_identity,
            budget: input.budget,
            grid,
        })
    }

    /// Complete mapper identity, including source, modal model, mappings, RNG,
    /// math semantics, reduction, spatial rules, artistic controls, and limits.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Modal model identity against which mappings and spatial rules were admitted.
    #[must_use]
    pub const fn modal_identity(&self) -> ContentHash {
        self.modal_identity
    }

    /// Durable canonical source-trajectory identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Selected exact source grid and reconstruction no-claim.
    #[must_use]
    pub const fn grid(&self) -> AudioExcitationGrid {
        self.grid
    }

    /// Selected reduction rule.
    #[must_use]
    pub const fn reduction(&self) -> AudioExcitationReduction {
        self.reduction
    }

    /// Verify the L6 sound config names this exact source, excitation mapper,
    /// modal model, algorithm versions, mappings, modes, and terminal class.
    pub fn validate_sound_configuration(
        &self,
        sound: &SoundSynthesisConfig,
    ) -> Result<(), AudioExcitationError> {
        if sound.input().trajectory.identity() != self.source_identity {
            return Err(AudioExcitationError::SoundConfigurationMismatch(
                "trajectory identity",
            ));
        }
        if sound.input().excitation.identity() != self.identity
            || sound.input().excitation.version() != AUDIO_EXCITATION_ALGORITHM_VERSION
        {
            return Err(AudioExcitationError::SoundConfigurationMismatch(
                "excitation identity or version",
            ));
        }
        if sound.input().sound_model.identity() != self.modal_identity
            || sound.input().sound_model.version() != MODAL_SYNTHESIS_ALGORITHM_VERSION
            || sound.input().audio_clock.ticks_per_second_numerator() != SOUND_MASTER_SAMPLE_RATE_HZ
            || sound.input().audio_clock.ticks_per_second_denominator() != 1
        {
            return Err(AudioExcitationError::SoundConfigurationMismatch(
                "modal model identity, version, or audio clock",
            ));
        }
        if sound.input().excitation_controls != self.mappings
            || sound.input().modes != self.modal_modes
        {
            return Err(AudioExcitationError::SoundConfigurationMismatch(
                "mappings or modes",
            ));
        }
        let final_disposition = self
            .controls
            .source()
            .samples()
            .last()
            .map(|sample| sample.input().disposition)
            .ok_or(AudioExcitationError::EmptySourceIntervals)?;
        let expected = match final_disposition {
            RenderSampleDisposition::TerminalInclination => {
                SoundTrajectoryDisposition::PhysicalTerminal
            }
            RenderSampleDisposition::HorizonCensored => SoundTrajectoryDisposition::HorizonCensored,
            RenderSampleDisposition::NumericalRefusal(_) => {
                return Err(AudioExcitationError::NumericalRefusalSource);
            }
            RenderSampleDisposition::Continue => {
                return Err(AudioExcitationError::SoundConfigurationMismatch(
                    "source terminal disposition",
                ));
            }
        };
        if sound.input().trajectory_disposition != expected {
            return Err(AudioExcitationError::SoundConfigurationMismatch(
                "trajectory disposition",
            ));
        }
        Ok(())
    }

    /// Construct the zero-progress checkpoint for this exact mapper.
    pub fn initial_checkpoint(
        &self,
        cx: &Cx<'_>,
    ) -> Result<AudioExcitationCheckpoint, AudioExcitationError> {
        checkpoint(cx)?;
        Ok(AudioExcitationCheckpoint {
            mapper_identity: self.identity,
            source_identity: self.source_identity,
            next_interval_index: 0,
            last_end_time_bits: None,
            cumulative_force_time: StemAccumulator::new(),
        })
    }

    /// Map the next bounded group of complete source intervals transactionally.
    /// The predecessor checkpoint is immutable and no partial chunk is returned.
    pub fn map_next_chunk(
        &self,
        prior: &AudioExcitationCheckpoint,
        maximum_intervals: NonZeroUsize,
        cx: &Cx<'_>,
    ) -> Result<AudioExcitationChunk, AudioExcitationError> {
        self.map_next_chunk_with_checkpoint(prior, maximum_intervals, &mut || checkpoint(cx))
    }

    fn map_next_chunk_with_checkpoint(
        &self,
        prior: &AudioExcitationCheckpoint,
        maximum_intervals: NonZeroUsize,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioExcitationError>,
    ) -> Result<AudioExcitationChunk, AudioExcitationError> {
        checkpoint_fn()?;
        self.validate_checkpoint(prior)?;
        if maximum_intervals.get() > self.budget.maximum_chunk_intervals {
            return Err(AudioExcitationError::ChunkIntervalBudgetExceeded {
                requested: maximum_intervals.get(),
                limit: self.budget.maximum_chunk_intervals,
            });
        }
        let total = selected_len(self.controls, &self.selected);
        if prior.next_interval_index == total {
            return Err(AudioExcitationError::Complete);
        }
        let end = prior
            .next_interval_index
            .saturating_add(maximum_intervals.get())
            .min(total);
        let count = end - prior.next_interval_index;
        let mut event_count = 0_usize;
        for (local, index) in (prior.next_interval_index..end).enumerate() {
            if local % AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS == 0 {
                checkpoint_fn()?;
            }
            let view = selected_view(self.controls, &self.selected, index);
            event_count = event_count.checked_add(view.events.len()).ok_or(
                AudioExcitationError::ChunkEventBudgetExceeded {
                    requested: usize::MAX,
                    limit: self.budget.maximum_chunk_events,
                },
            )?;
        }
        if event_count > self.budget.maximum_chunk_events {
            return Err(AudioExcitationError::ChunkEventBudgetExceeded {
                requested: event_count,
                limit: self.budget.maximum_chunk_events,
            });
        }
        let spatial_envelope_count = count.checked_mul(self.mode_rules.len()).ok_or(
            AudioExcitationError::ChunkSpatialEnvelopeBudgetExceeded {
                requested: usize::MAX,
                limit: self.budget.maximum_chunk_spatial_envelopes,
            },
        )?;
        if spatial_envelope_count > self.budget.maximum_chunk_spatial_envelopes {
            return Err(AudioExcitationError::ChunkSpatialEnvelopeBudgetExceeded {
                requested: spatial_envelope_count,
                limit: self.budget.maximum_chunk_spatial_envelopes,
            });
        }

        let mut intervals = Vec::new();
        intervals
            .try_reserve_exact(count)
            .map_err(|_| AudioExcitationError::Capacity {
                artifact: "mapped excitation intervals",
                requested: count,
            })?;
        let mut cumulative = prior.cumulative_force_time;
        let mut chunk_accumulator = StemAccumulator::new();
        let mut maximum_abs_residual = 0.0_f64;
        for (local, index) in (prior.next_interval_index..end).enumerate() {
            if local % AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS == 0 {
                checkpoint_fn()?;
            }
            let mapped = self.map_interval(index, checkpoint_fn)?;
            cumulative.add(mapped.force_time_stems_n_s);
            chunk_accumulator.add(mapped.force_time_stems_n_s);
            maximum_abs_residual =
                maximum_abs_residual.max(max_abs_stems(mapped.measure_residual_stems_n_s));
            intervals.push(mapped);
        }
        checkpoint_fn()?;
        let last_end_time_bits = intervals
            .last()
            .map(|interval| interval.end_time_s.to_bits());
        let successor = AudioExcitationCheckpoint {
            mapper_identity: self.identity,
            source_identity: self.source_identity,
            next_interval_index: end,
            last_end_time_bits,
            cumulative_force_time: cumulative,
        };
        let diagnostics = AudioExcitationDiagnostics {
            start_interval_index: prior.next_interval_index,
            end_interval_index: end,
            event_count,
            maximum_abs_measure_residual_n_s: maximum_abs_residual,
            chunk_force_time_stems_n_s: chunk_accumulator.total(),
            cumulative_force_time_stems_n_s: cumulative.total(),
        };
        let identity = chunk_identity(self.identity, prior.next_interval_index, end);
        checkpoint_fn()?;
        Ok(AudioExcitationChunk {
            identity,
            intervals,
            diagnostics,
            successor,
        })
    }

    fn validate_checkpoint(
        &self,
        checkpoint: &AudioExcitationCheckpoint,
    ) -> Result<(), AudioExcitationError> {
        let total = selected_len(self.controls, &self.selected);
        if checkpoint.mapper_identity != self.identity
            || checkpoint.source_identity != self.source_identity
            || checkpoint.next_interval_index > total
        {
            return Err(AudioExcitationError::InvalidCheckpoint);
        }
        let expected_last = checkpoint.next_interval_index.checked_sub(1).map(|index| {
            selected_view(self.controls, &self.selected, index)
                .end_time_s
                .to_bits()
        });
        if checkpoint.last_end_time_bits != expected_last
            || !checkpoint.cumulative_force_time.is_finite()
        {
            return Err(AudioExcitationError::InvalidCheckpoint);
        }
        Ok(())
    }

    fn map_interval(
        &self,
        interval_index: usize,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioExcitationError>,
    ) -> Result<AudioExcitationInterval, AudioExcitationError> {
        checkpoint_fn()?;
        let view = selected_view(self.controls, &self.selected, interval_index);
        validate_events(view, interval_index)?;
        let mut mean_stems = AudioExcitationStems::ZERO;
        let mut measure_stems = AudioExcitationStems::ZERO;
        let mut availability = initial_availability(&self.mappings);
        for (mapping_index, mapping) in self.mappings.iter().enumerate() {
            if mapping_index % AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS == 0 {
                checkpoint_fn()?;
            }
            let scalar = source_scalar(
                view,
                mapping.channel,
                self.base_normal_world,
                interval_index,
            )?;
            mark_availability(&mut availability, scalar.class, scalar.available);
            if !scalar.available {
                continue;
            }
            let source_limit = match scalar.unit {
                SourceUnit::Newton => self.budget.maximum_abs_source_force_n,
                SourceUnit::Watt => self.budget.maximum_abs_source_work_rate_w,
            };
            check_limit(
                scalar.mean.abs(),
                source_limit,
                interval_index,
                scalar.unit.limit_field(),
            )?;
            let mapped_mean = scalar.mean * mapping.source_scale;
            let mapped_measure = scalar.measure * mapping.source_scale;
            if !mapped_mean.is_finite() || !mapped_measure.is_finite() {
                return Err(AudioExcitationError::NonFinite {
                    interval: interval_index,
                    field: "mapped generalized force or measure",
                });
            }
            add_to_stem(
                &mut mean_stems,
                scalar.class,
                mapping.target_component,
                mapped_mean,
            );
            add_to_stem(
                &mut measure_stems,
                scalar.class,
                mapping.target_component,
                mapped_measure,
            );
        }
        validate_stem_limits(
            mean_stems,
            self.budget.maximum_abs_generalized_force_n,
            interval_index,
            "mapped stem force",
        )?;
        validate_stem_limits(
            measure_stems,
            self.budget.maximum_abs_force_time_measure_n_s,
            interval_index,
            "mapped stem force-time measure",
        )?;
        let mean_generalized_force_n = mean_stems.sum();
        let generalized_force_time_n_s = measure_stems.sum();
        check_components_limit(
            mean_generalized_force_n,
            self.budget.maximum_abs_generalized_force_n,
            interval_index,
            "summed generalized force",
        )?;
        check_components_limit(
            generalized_force_time_n_s,
            self.budget.maximum_abs_force_time_measure_n_s,
            interval_index,
            "summed generalized force-time measure",
        )?;
        let measure_residual_stems_n_s =
            stem_sub(stem_scale(mean_stems, view.duration_s), measure_stems);
        validate_stem_limits(
            measure_residual_stems_n_s,
            self.budget.maximum_measure_residual_n_s,
            interval_index,
            "per-source mean-measure reconciliation residual",
        )?;
        let measure_residual_n_s = measure_residual_stems_n_s.sum();
        check_components_limit(
            measure_residual_n_s,
            self.budget.maximum_measure_residual_n_s,
            interval_index,
            "mean-measure reconciliation residual",
        )?;

        let spatial_envelopes = self.spatial_envelopes(view, interval_index, checkpoint_fn)?;
        let artistic_texture = self.texture_envelope(view, interval_index)?;
        let mut events = Vec::new();
        events.try_reserve_exact(view.events.len()).map_err(|_| {
            AudioExcitationError::Capacity {
                artifact: "mapped excitation events",
                requested: view.events.len(),
            }
        })?;
        for (event_index, event) in view.events.iter().copied().enumerate() {
            if event_index % AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS == 0 {
                checkpoint_fn()?;
            }
            events.push(self.map_event(event, event_index, interval_index)?);
        }
        Ok(AudioExcitationInterval {
            first_source_sample_index: view.first_source_sample_index,
            last_source_sample_index: view.last_source_sample_index,
            start_time_s: view.start_time_s,
            end_time_s: view.end_time_s,
            duration_s: view.duration_s,
            visual_coverage: view.visual_coverage,
            mean_force_stems_n: mean_stems,
            force_time_stems_n_s: measure_stems,
            mean_generalized_force_n,
            generalized_force_time_n_s,
            availability,
            events,
            event_barrier: view.event_barrier,
            spatial_envelopes,
            artistic_texture,
            measure_residual_stems_n_s,
            measure_residual_n_s,
        })
    }

    fn spatial_envelopes(
        &self,
        view: SourceIntervalView<'_>,
        interval_index: usize,
        checkpoint_fn: &mut impl FnMut() -> Result<(), AudioExcitationError>,
    ) -> Result<Vec<ModalSpatialEnvelope>, AudioExcitationError> {
        let mut envelopes = Vec::new();
        envelopes
            .try_reserve_exact(self.mode_rules.len())
            .map_err(|_| AudioExcitationError::Capacity {
                artifact: "modal spatial envelopes",
                requested: self.mode_rules.len(),
            })?;
        let start_contact = view
            .visual_coverage
            .start_visualization_index
            .and_then(|index| self.controls.visualization().get(index))
            .and_then(|point| point.contact);
        let end_contact = self
            .controls
            .visualization()
            .get(view.visual_coverage.end_visualization_index)
            .and_then(|point| point.contact);
        for (rule_index, rule) in self.mode_rules.iter().enumerate() {
            if rule_index % AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS == 0 {
                checkpoint_fn()?;
            }
            let envelope = match rule.shape {
                ContactModeShape::Uniform => ModalSpatialEnvelope {
                    mode_id: rule.mode_id,
                    start_factor: 1.0,
                    end_factor: 1.0,
                    source: SpatialEnvelopeSource::DeclaredStatic,
                },
                ContactModeShape::AzimuthalCosine {
                    harmonic,
                    phase_rad,
                } => {
                    if view.event_barrier || !view.events.is_empty() {
                        let [event] = view.events else {
                            return Err(AudioExcitationError::InvalidSpatialPolicy(
                                "nonuniform participation cannot reconstruct multiple events",
                            ));
                        };
                        let (contact, source) = match event.kind {
                            ContactTransitionKind::Opening => {
                                (start_contact, SpatialEnvelopeSource::HeldStartEndpoint)
                            }
                            ContactTransitionKind::Reimpact => {
                                (end_contact, SpatialEnvelopeSource::HeldEndEndpoint)
                            }
                        };
                        let factor = contact
                            .ok_or(AudioExcitationError::InvalidSpatialPolicy(
                                "event-side contact coordinate unavailable",
                            ))
                            .and_then(|contact| {
                                spatial_factor(contact, rule, harmonic, phase_rad)
                            })?;
                        ModalSpatialEnvelope {
                            mode_id: rule.mode_id,
                            start_factor: factor,
                            end_factor: factor,
                            source,
                        }
                    } else {
                        let start = start_contact
                            .map(|contact| spatial_factor(contact, rule, harmonic, phase_rad));
                        let end = end_contact
                            .map(|contact| spatial_factor(contact, rule, harmonic, phase_rad));
                        match (start.transpose()?, end.transpose()?) {
                            (Some(start_factor), Some(end_factor)) => ModalSpatialEnvelope {
                                mode_id: rule.mode_id,
                                start_factor,
                                end_factor,
                                source: SpatialEnvelopeSource::ExactEndpointInterpolation,
                            },
                            (Some(factor), None) => ModalSpatialEnvelope {
                                mode_id: rule.mode_id,
                                start_factor: factor,
                                end_factor: factor,
                                source: SpatialEnvelopeSource::HeldStartEndpoint,
                            },
                            (None, Some(factor)) => ModalSpatialEnvelope {
                                mode_id: rule.mode_id,
                                start_factor: factor,
                                end_factor: factor,
                                source: SpatialEnvelopeSource::HeldEndEndpoint,
                            },
                            (None, None) => ModalSpatialEnvelope {
                                mode_id: rule.mode_id,
                                start_factor: 1.0,
                                end_factor: 1.0,
                                source: SpatialEnvelopeSource::MissingContactStatic,
                            },
                        }
                    }
                }
            };
            if !envelope.start_factor.is_finite() || !envelope.end_factor.is_finite() {
                return Err(AudioExcitationError::NonFinite {
                    interval: interval_index,
                    field: "modal spatial envelope",
                });
            }
            envelopes.push(envelope);
        }
        Ok(envelopes)
    }

    fn texture_envelope(
        &self,
        view: SourceIntervalView<'_>,
        interval_index: usize,
    ) -> Result<Option<ArtisticTextureEnvelope>, AudioExcitationError> {
        let Some(texture) = self.artistic_texture else {
            return Ok(None);
        };
        let rolling = source_channel(view.channels.rolling);
        let (availability, magnitude_w) = match rolling {
            Some(channel) => (
                ExcitationSourceAvailability::Available,
                channel.signed_mean_work_rate_w.abs(),
            ),
            None => (ExcitationSourceAvailability::Unavailable, 0.0),
        };
        check_limit(
            magnitude_w,
            self.budget.maximum_abs_source_work_rate_w,
            interval_index,
            "artistic rolling source work rate",
        )?;
        let envelope_n = magnitude_w * texture.rolling_force_gain_n_per_w;
        if !envelope_n.is_finite() {
            return Err(AudioExcitationError::NonFinite {
                interval: interval_index,
                field: "artistic rolling texture envelope",
            });
        }
        check_limit(
            envelope_n,
            self.budget.maximum_artistic_texture_envelope_n,
            interval_index,
            "artistic rolling texture envelope",
        )?;
        let root =
            self.texture_root_identity
                .ok_or(AudioExcitationError::InvalidArtisticTexture(
                    "missing root identity",
                ))?;
        let stream_identity = interval_texture_identity(root, view);
        Ok(Some(ArtisticTextureEnvelope {
            stream_identity,
            peak_force_envelope_n: component_value(texture.rolling_target_component, envelope_n),
            band_low_hz: texture.band_low_hz,
            band_high_hz: texture.band_high_hz,
            rolling_availability: availability,
        }))
    }

    fn map_event(
        &self,
        event: ControlContactEvent,
        event_index: usize,
        interval_index: usize,
    ) -> Result<AudioExcitationEvent, AudioExcitationError> {
        let artistic = match (self.artistic_texture, event.kind) {
            (Some(texture), ContactTransitionKind::Reimpact)
                if texture.reimpact_impulse_n_s != 0.0 =>
            {
                let root = self.texture_root_identity.ok_or(
                    AudioExcitationError::InvalidArtisticTexture("missing root identity"),
                )?;
                let stream_identity = event_texture_identity(root, event, event_index);
                let unit = procedural_texture_unit_sample(stream_identity, 0);
                let impulse = texture.reimpact_impulse_n_s * unit;
                check_limit(
                    impulse.abs(),
                    self.budget.maximum_artistic_event_impulse_n_s,
                    interval_index,
                    "artistic event impulse",
                )?;
                Some(ArtisticEventExcitation {
                    stream_identity,
                    impulse_n_s: component_value(texture.reimpact_target_component, impulse),
                })
            }
            _ => None,
        };
        Ok(AudioExcitationEvent {
            source_sample_index: event.source_sample_index,
            kind: event.kind,
            time_s: event.time_s,
            bracket_start_s: event.bracket_start_s,
            bracket_end_s: event.bracket_end_s,
            measure: event.measure,
            physical_impulse_n_s: ModalComponentValues::ZERO,
            artistic,
        })
    }
}

/// Random-access deterministic unit noise in `[-1,1)`. The caller still owns
/// band-limiting at its admitted output clock.
#[must_use]
pub fn procedural_texture_unit_sample(
    stream_identity: ContentHash,
    absolute_sample_index: u64,
) -> f64 {
    let bytes = stream_identity.as_bytes();
    let seed = u64::from_le_bytes(bytes[0..8].try_into().expect("fixed hash width"));
    let kernel = u32::from_le_bytes(bytes[8..12].try_into().expect("fixed hash width"));
    let tile = u32::from_le_bytes(bytes[12..16].try_into().expect("fixed hash width"));
    let block = Stream::at(
        RandomStreamKey { seed, kernel, tile },
        absolute_sample_index,
    );
    let bits = (u64::from(block[1]) << 32) | u64::from(block[0]);
    let uniform = (bits >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0);
    uniform.mul_add(2.0, -1.0)
}

#[derive(Clone, Copy)]
struct SourceIntervalView<'a> {
    first_source_sample_index: usize,
    last_source_sample_index: usize,
    start_time_s: f64,
    end_time_s: f64,
    duration_s: f64,
    visual_coverage: AudioVisualCoverage,
    mean_base_normal_contact_force_n: Option<f64>,
    channels: ChannelControlSet,
    events: &'a [ControlContactEvent],
    event_barrier: bool,
}

fn selected_len(controls: &EulerControlStream<'_>, selected: &SelectedIntervals<'_>) -> usize {
    match selected {
        SelectedIntervals::Raw => controls.audio().len(),
        SelectedIntervals::Coarsened(coarsened) => coarsened.bins().len(),
    }
}

fn selected_view<'a>(
    controls: &'a EulerControlStream<'_>,
    selected: &'a SelectedIntervals<'_>,
    index: usize,
) -> SourceIntervalView<'a> {
    match selected {
        SelectedIntervals::Raw => {
            let interval = &controls.audio()[index];
            SourceIntervalView {
                first_source_sample_index: interval.source_sample_index,
                last_source_sample_index: interval.source_sample_index,
                start_time_s: interval.start_time_s,
                end_time_s: interval.end_time_s,
                duration_s: interval.duration_s,
                visual_coverage: interval.visual_coverage,
                mean_base_normal_contact_force_n: interval.mean_base_normal_contact_force_n,
                channels: interval.channels,
                events: &interval.events,
                event_barrier: !interval.events.is_empty(),
            }
        }
        SelectedIntervals::Coarsened(coarsened) => {
            let bin = &coarsened.bins()[index];
            SourceIntervalView {
                first_source_sample_index: bin.first_source_sample_index,
                last_source_sample_index: bin.last_source_sample_index,
                start_time_s: bin.start_time_s,
                end_time_s: bin.end_time_s,
                duration_s: bin.duration_s,
                visual_coverage: bin.visual_coverage,
                mean_base_normal_contact_force_n: bin.mean_base_normal_contact_force_n,
                channels: bin.channels,
                events: &bin.events,
                event_barrier: bin.event_barrier,
            }
        }
    }
}

fn derive_grid(
    controls: &EulerControlStream<'_>,
    selected: &SelectedIntervals<'_>,
    cx: &Cx<'_>,
) -> Result<AudioExcitationGrid, AudioExcitationError> {
    let count = selected_len(controls, selected);
    let first = selected_view(controls, selected, 0);
    let mut minimum = first.duration_s;
    let mut maximum = first.duration_s;
    let mut previous_end = first.start_time_s;
    for index in 0..count {
        if index % AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS == 0 {
            checkpoint(cx)?;
        }
        let view = selected_view(controls, selected, index);
        if !view.start_time_s.is_finite()
            || !view.end_time_s.is_finite()
            || !view.duration_s.is_finite()
            || view.duration_s <= 0.0
        {
            return Err(AudioExcitationError::NonFinite {
                interval: index,
                field: "selected source grid",
            });
        }
        if (view.end_time_s - view.start_time_s).to_bits() != view.duration_s.to_bits() {
            return Err(AudioExcitationError::InvalidSourceGrid {
                interval: index,
                field: "duration",
            });
        }
        if index > 0 && view.start_time_s.to_bits() != previous_end.to_bits() {
            return Err(AudioExcitationError::InvalidSourceGrid {
                interval: index,
                field: "contiguity",
            });
        }
        validate_events(view, index)?;
        minimum = minimum.min(view.duration_s);
        maximum = maximum.max(view.duration_s);
        previous_end = view.end_time_s;
    }
    let nominal = controls.source().metadata().timestep_s;
    let nyquist = 0.5 / maximum;
    if !nominal.is_finite() || nominal <= 0.0 || !nyquist.is_finite() {
        return Err(AudioExcitationError::NonFinite {
            interval: 0,
            field: "source cadence metadata",
        });
    }
    Ok(AudioExcitationGrid {
        interval_count: count,
        start_time_s: first.start_time_s,
        end_time_s: previous_end,
        minimum_interval_duration_s: minimum,
        maximum_interval_duration_s: maximum,
        nominal_mechanics_timestep_s: nominal,
        nominal_source_nyquist_ceiling_hz: nyquist,
        reconstruction: AudioExcitationReconstructionStatus::RequiresBandLimitedResampling,
    })
}

fn validate_mappings(mappings: &[SoundExcitationControl]) -> Result<(), AudioExcitationError> {
    if mappings.is_empty()
        || !mappings.windows(2).all(|pair| {
            (pair[0].channel, pair[0].target_component)
                < (pair[1].channel, pair[1].target_component)
        })
        || mappings
            .iter()
            .any(|mapping| !mapping.source_scale.is_finite() || mapping.source_scale == 0.0)
    {
        return Err(AudioExcitationError::InvalidMappings);
    }
    for mapping in mappings {
        match mapping.channel {
            SoundExcitationChannel::ContactNormalForce
            | SoundExcitationChannel::ContactSignedWorkRate
            | SoundExcitationChannel::RollingSignedWorkRate
            | SoundExcitationChannel::BaseDampingSignedWorkRate
            | SoundExcitationChannel::ExteriorGasBodySignedWorkRate => {}
            unsupported => return Err(AudioExcitationError::UnsupportedMapping(unsupported)),
        }
    }
    Ok(())
}

fn validate_budget(budget: AudioExcitationBudget) -> Result<(), AudioExcitationError> {
    if budget.maximum_total_intervals == 0 {
        return Err(AudioExcitationError::InvalidBudget(
            "maximum_total_intervals",
        ));
    }
    if budget.maximum_chunk_intervals == 0
        || budget.maximum_chunk_intervals > MAX_AUDIO_EXCITATION_CHUNK_INTERVALS
        || budget.maximum_chunk_intervals > budget.maximum_total_intervals
    {
        return Err(AudioExcitationError::InvalidBudget(
            "maximum_chunk_intervals",
        ));
    }
    if budget.maximum_chunk_events == 0 {
        return Err(AudioExcitationError::InvalidBudget("maximum_chunk_events"));
    }
    if budget.maximum_chunk_spatial_envelopes == 0 {
        return Err(AudioExcitationError::InvalidBudget(
            "maximum_chunk_spatial_envelopes",
        ));
    }
    for (field, value, strictly_positive) in [
        (
            "maximum_abs_source_force_n",
            budget.maximum_abs_source_force_n,
            true,
        ),
        (
            "maximum_abs_source_work_rate_w",
            budget.maximum_abs_source_work_rate_w,
            true,
        ),
        (
            "maximum_abs_generalized_force_n",
            budget.maximum_abs_generalized_force_n,
            true,
        ),
        (
            "maximum_abs_force_time_measure_n_s",
            budget.maximum_abs_force_time_measure_n_s,
            true,
        ),
        (
            "maximum_measure_residual_n_s",
            budget.maximum_measure_residual_n_s,
            false,
        ),
        (
            "maximum_artistic_texture_envelope_n",
            budget.maximum_artistic_texture_envelope_n,
            false,
        ),
        (
            "maximum_artistic_event_impulse_n_s",
            budget.maximum_artistic_event_impulse_n_s,
            false,
        ),
    ] {
        if !value.is_finite() || value < 0.0 || (strictly_positive && value == 0.0) {
            return Err(AudioExcitationError::InvalidBudget(field));
        }
    }
    Ok(())
}

fn validate_texture(
    texture: Option<ArtisticTextureConfig>,
    budget: AudioExcitationBudget,
) -> Result<(), AudioExcitationError> {
    let Some(texture) = texture else {
        return Ok(());
    };
    if !texture.rolling_force_gain_n_per_w.is_finite() || texture.rolling_force_gain_n_per_w < 0.0 {
        return Err(AudioExcitationError::InvalidArtisticTexture(
            "rolling_force_gain_n_per_w",
        ));
    }
    if !texture.band_low_hz.is_finite()
        || !texture.band_high_hz.is_finite()
        || texture.band_low_hz < 0.0
        || texture.band_low_hz >= texture.band_high_hz
        || texture.band_high_hz >= f64::from(SOUND_MASTER_SAMPLE_RATE_HZ) * 0.5
    {
        return Err(AudioExcitationError::InvalidArtisticTexture("texture band"));
    }
    if !texture.reimpact_impulse_n_s.is_finite()
        || texture.reimpact_impulse_n_s.abs() > budget.maximum_artistic_event_impulse_n_s
    {
        return Err(AudioExcitationError::InvalidArtisticTexture(
            "reimpact_impulse_n_s",
        ));
    }
    if texture.rolling_force_gain_n_per_w == 0.0 && texture.reimpact_impulse_n_s == 0.0 {
        return Err(AudioExcitationError::InvalidArtisticTexture(
            "empty artistic texture",
        ));
    }
    Ok(())
}

fn prepare_mode_rules(
    modal: &ModalSynthesisModel,
    policy: &mut ContactParticipationPolicy,
) -> Result<Vec<PreparedModeRule>, AudioExcitationError> {
    match policy {
        ContactParticipationPolicy::DeclaredStatic => Ok(modal
            .modes()
            .iter()
            .map(|mode| PreparedModeRule {
                mode_id: mode.mode_id,
                component: mode.component,
                shape: ContactModeShape::Uniform,
            })
            .collect()),
        ContactParticipationPolicy::ContactCoordinates { rules } => {
            rules.sort_by_key(|rule| rule.mode_id);
            if rules.len() != modal.modes().len()
                || rules
                    .windows(2)
                    .any(|pair| pair[0].mode_id == pair[1].mode_id)
            {
                return Err(AudioExcitationError::InvalidSpatialPolicy(
                    "mode rule coverage",
                ));
            }
            let mut prepared = Vec::new();
            prepared.try_reserve_exact(rules.len()).map_err(|_| {
                AudioExcitationError::Capacity {
                    artifact: "prepared spatial rules",
                    requested: rules.len(),
                }
            })?;
            for (rule, mode) in rules.iter().zip(modal.modes()) {
                if rule.mode_id != mode.mode_id {
                    return Err(AudioExcitationError::InvalidSpatialPolicy(
                        "mode rule identity",
                    ));
                }
                if let ContactModeShape::AzimuthalCosine {
                    harmonic,
                    phase_rad,
                } = rule.shape
                    && (harmonic == 0
                        || harmonic > MAX_AUDIO_EXCITATION_AZIMUTHAL_HARMONIC
                        || !phase_rad.is_finite())
                {
                    return Err(AudioExcitationError::InvalidSpatialPolicy(
                        "azimuthal shape",
                    ));
                }
                prepared.push(PreparedModeRule {
                    mode_id: mode.mode_id,
                    component: mode.component,
                    shape: rule.shape,
                });
            }
            Ok(prepared)
        }
    }
}

fn validate_events(
    view: SourceIntervalView<'_>,
    interval_index: usize,
) -> Result<(), AudioExcitationError> {
    let mut previous = None;
    for (event_index, event) in view.events.iter().enumerate() {
        if event.measure != ContactEventMeasure::TimingOnly
            || event.source_sample_index < view.first_source_sample_index
            || event.source_sample_index > view.last_source_sample_index
            || !event.time_s.is_finite()
            || !event.bracket_start_s.is_finite()
            || !event.bracket_end_s.is_finite()
            || event.bracket_start_s > event.time_s
            || event.time_s > event.bracket_end_s
            || event.time_s < view.start_time_s
            || event.time_s > view.end_time_s
            || previous.is_some_and(|time| event.time_s <= time)
        {
            return Err(AudioExcitationError::InvalidEventOrder {
                interval: interval_index,
                event: event_index,
            });
        }
        previous = Some(event.time_s);
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum SourceClass {
    Contact,
    Rolling,
    Base,
    Gas,
}

#[derive(Clone, Copy)]
enum SourceUnit {
    Newton,
    Watt,
}

impl SourceUnit {
    const fn limit_field(self) -> &'static str {
        match self {
            Self::Newton => "source force",
            Self::Watt => "source signed work rate",
        }
    }
}

#[derive(Clone, Copy)]
struct SourceScalar {
    class: SourceClass,
    unit: SourceUnit,
    available: bool,
    mean: f64,
    measure: f64,
}

fn source_scalar(
    view: SourceIntervalView<'_>,
    channel: SoundExcitationChannel,
    base_normal_world: Vec3,
    interval: usize,
) -> Result<SourceScalar, AudioExcitationError> {
    let (class, unit, available, mean, measure) = match channel {
        SoundExcitationChannel::ContactNormalForce => {
            let retained = source_channel(view.channels.contact);
            let available = view.mean_base_normal_contact_force_n.is_some();
            let mean = view.mean_base_normal_contact_force_n.unwrap_or(0.0);
            let measure = retained.map_or_else(
                || {
                    if available {
                        mean * view.duration_s
                    } else {
                        0.0
                    }
                },
                |value| value.force_time_measure_world_n_s.dot(base_normal_world),
            );
            (
                SourceClass::Contact,
                SourceUnit::Newton,
                available,
                mean,
                measure,
            )
        }
        SoundExcitationChannel::ContactSignedWorkRate => {
            work_scalar(SourceClass::Contact, view.channels.contact)
        }
        SoundExcitationChannel::RollingSignedWorkRate => {
            work_scalar(SourceClass::Rolling, view.channels.rolling)
        }
        SoundExcitationChannel::BaseDampingSignedWorkRate => {
            work_scalar(SourceClass::Base, view.channels.base)
        }
        SoundExcitationChannel::ExteriorGasBodySignedWorkRate => {
            work_scalar(SourceClass::Gas, view.channels.gas)
        }
        unsupported => return Err(AudioExcitationError::UnsupportedMapping(unsupported)),
    };
    if available && (!mean.is_finite() || !measure.is_finite()) {
        return Err(AudioExcitationError::NonFinite {
            interval,
            field: "source scalar or measure",
        });
    }
    Ok(SourceScalar {
        class,
        unit,
        available,
        mean,
        measure,
    })
}

fn work_scalar(
    class: SourceClass,
    channel: ChannelControl,
) -> (SourceClass, SourceUnit, bool, f64, f64) {
    match source_channel(channel) {
        Some(value) => (
            class,
            SourceUnit::Watt,
            true,
            value.signed_mean_work_rate_w,
            value.signed_work_j,
        ),
        None => (class, SourceUnit::Watt, false, 0.0, 0.0),
    }
}

fn source_channel(channel: ChannelControl) -> Option<crate::AvailableChannelControl> {
    channel.available()
}

fn initial_availability(mappings: &[SoundExcitationControl]) -> AudioExcitationAvailability {
    let mut result = AudioExcitationAvailability {
        contact: ExcitationSourceAvailability::NotMapped,
        rolling: ExcitationSourceAvailability::NotMapped,
        base: ExcitationSourceAvailability::NotMapped,
        gas: ExcitationSourceAvailability::NotMapped,
    };
    for mapping in mappings {
        let target = match mapping.channel {
            SoundExcitationChannel::ContactNormalForce
            | SoundExcitationChannel::ContactSignedWorkRate => &mut result.contact,
            SoundExcitationChannel::RollingSignedWorkRate => &mut result.rolling,
            SoundExcitationChannel::BaseDampingSignedWorkRate => &mut result.base,
            SoundExcitationChannel::ExteriorGasBodySignedWorkRate => &mut result.gas,
            _ => continue,
        };
        *target = ExcitationSourceAvailability::Unavailable;
    }
    result
}

fn mark_availability(
    availability: &mut AudioExcitationAvailability,
    class: SourceClass,
    available: bool,
) {
    let target = match class {
        SourceClass::Contact => &mut availability.contact,
        SourceClass::Rolling => &mut availability.rolling,
        SourceClass::Base => &mut availability.base,
        SourceClass::Gas => &mut availability.gas,
    };
    if available {
        *target = ExcitationSourceAvailability::Available;
    }
}

fn add_to_stem(
    stems: &mut AudioExcitationStems,
    class: SourceClass,
    component: SoundModalComponent,
    value: f64,
) {
    let stem = match class {
        SourceClass::Contact => &mut stems.contact,
        SourceClass::Rolling => &mut stems.rolling,
        SourceClass::Base => &mut stems.base,
        SourceClass::Gas => &mut stems.gas,
    };
    add_component(stem, component, value);
}

fn add_component(values: &mut ModalComponentValues, component: SoundModalComponent, value: f64) {
    match component {
        SoundModalComponent::Disc => values.disc += value,
        SoundModalComponent::GlassPlate => values.glass_plate += value,
        SoundModalComponent::BaseAssembly => values.base_assembly += value,
    }
}

fn component_value(component: SoundModalComponent, value: f64) -> ModalComponentValues {
    let mut result = ModalComponentValues::ZERO;
    add_component(&mut result, component, value);
    result
}

fn component_scale(values: ModalComponentValues, scale: f64) -> ModalComponentValues {
    ModalComponentValues {
        disc: values.disc * scale,
        glass_plate: values.glass_plate * scale,
        base_assembly: values.base_assembly * scale,
    }
}

fn component_add(left: ModalComponentValues, right: ModalComponentValues) -> ModalComponentValues {
    ModalComponentValues {
        disc: left.disc + right.disc,
        glass_plate: left.glass_plate + right.glass_plate,
        base_assembly: left.base_assembly + right.base_assembly,
    }
}

fn component_sub(left: ModalComponentValues, right: ModalComponentValues) -> ModalComponentValues {
    ModalComponentValues {
        disc: left.disc - right.disc,
        glass_plate: left.glass_plate - right.glass_plate,
        base_assembly: left.base_assembly - right.base_assembly,
    }
}

fn stem_scale(stems: AudioExcitationStems, scale: f64) -> AudioExcitationStems {
    AudioExcitationStems {
        contact: component_scale(stems.contact, scale),
        rolling: component_scale(stems.rolling, scale),
        base: component_scale(stems.base, scale),
        gas: component_scale(stems.gas, scale),
    }
}

fn stem_sub(left: AudioExcitationStems, right: AudioExcitationStems) -> AudioExcitationStems {
    AudioExcitationStems {
        contact: component_sub(left.contact, right.contact),
        rolling: component_sub(left.rolling, right.rolling),
        base: component_sub(left.base, right.base),
        gas: component_sub(left.gas, right.gas),
    }
}

fn max_abs_components(values: ModalComponentValues) -> f64 {
    values
        .disc
        .abs()
        .max(values.glass_plate.abs())
        .max(values.base_assembly.abs())
}

fn max_abs_stems(stems: AudioExcitationStems) -> f64 {
    [stems.contact, stems.rolling, stems.base, stems.gas]
        .into_iter()
        .map(max_abs_components)
        .fold(0.0, f64::max)
}

fn check_components_limit(
    values: ModalComponentValues,
    limit: f64,
    interval: usize,
    field: &'static str,
) -> Result<(), AudioExcitationError> {
    if !values.disc.is_finite()
        || !values.glass_plate.is_finite()
        || !values.base_assembly.is_finite()
    {
        return Err(AudioExcitationError::NonFinite { interval, field });
    }
    check_limit(max_abs_components(values), limit, interval, field)
}

fn validate_stem_limits(
    stems: AudioExcitationStems,
    limit: f64,
    interval: usize,
    field: &'static str,
) -> Result<(), AudioExcitationError> {
    for values in [stems.contact, stems.rolling, stems.base, stems.gas] {
        check_components_limit(values, limit, interval, field)?;
    }
    Ok(())
}

fn check_limit(
    magnitude: f64,
    limit: f64,
    interval: usize,
    field: &'static str,
) -> Result<(), AudioExcitationError> {
    if !magnitude.is_finite() {
        return Err(AudioExcitationError::NonFinite { interval, field });
    }
    if magnitude > limit {
        return Err(AudioExcitationError::LimitExceeded {
            interval,
            field,
            magnitude,
            limit,
        });
    }
    Ok(())
}

fn spatial_factor(
    contact: crate::ContactFrameCoordinates,
    rule: &PreparedModeRule,
    harmonic: u16,
    phase_rad: f64,
) -> Result<f64, AudioExcitationError> {
    let point = match rule.component {
        SoundModalComponent::Disc => contact.point_body_m,
        SoundModalComponent::GlassPlate | SoundModalComponent::BaseAssembly => contact.point_base_m,
    };
    if !point.x.is_finite() || !point.y.is_finite() || (point.x == 0.0 && point.y == 0.0) {
        return Err(AudioExcitationError::InvalidSpatialPolicy(
            "degenerate local contact azimuth",
        ));
    }
    let angle = det::atan2(point.y, point.x);
    let factor = det::cos(f64::from(harmonic).mul_add(angle, phase_rad));
    if factor.is_finite() {
        Ok(factor)
    } else {
        Err(AudioExcitationError::InvalidSpatialPolicy(
            "non-finite contact shape",
        ))
    }
}

fn mapper_identity(
    source_identity: ContentHash,
    modal_identity: ContentHash,
    mappings: &[SoundExcitationControl],
    reduction: AudioExcitationReduction,
    mode_rules: &[PreparedModeRule],
    texture: Option<ArtisticTextureConfig>,
    budget: AudioExcitationBudget,
    grid: AudioExcitationGrid,
) -> ContentHash {
    let mut bytes = Vec::with_capacity(512 + mappings.len() * 24 + mode_rules.len() * 32);
    push_u32(&mut bytes, AUDIO_EXCITATION_ALGORITHM_VERSION);
    bytes.extend_from_slice(&EULER_CONTROL_STREAM_SCHEMA_VERSION.to_le_bytes());
    push_u32(&mut bytes, STRICT_CORE_SEMANTICS_VERSION);
    push_u64(&mut bytes, STRICT_CORE_GOLDEN_HASH);
    push_u32(&mut bytes, STREAM_SEMANTICS_VERSION);
    bytes.extend_from_slice(source_identity.as_bytes());
    bytes.extend_from_slice(modal_identity.as_bytes());
    push_u32(&mut bytes, mappings.len() as u32);
    for mapping in mappings {
        bytes.push(mapping.channel as u8);
        bytes.push(mapping.target_component as u8);
        push_f64(&mut bytes, mapping.source_scale);
    }
    match reduction {
        AudioExcitationReduction::RawIntervals => bytes.push(1),
        AudioExcitationReduction::WholeIntervalBoxcarV1 { intervals_per_bin } => {
            bytes.push(2);
            push_u64(&mut bytes, intervals_per_bin.get() as u64);
            bytes.push(AudioControlFilter::WholeIntervalBoxcarV1 as u8);
        }
    }
    push_u32(&mut bytes, mode_rules.len() as u32);
    for rule in mode_rules {
        push_u32(&mut bytes, rule.mode_id);
        bytes.push(rule.component as u8);
        match rule.shape {
            ContactModeShape::Uniform => bytes.push(1),
            ContactModeShape::AzimuthalCosine {
                harmonic,
                phase_rad,
            } => {
                bytes.push(2);
                bytes.extend_from_slice(&harmonic.to_le_bytes());
                push_f64(&mut bytes, phase_rad);
            }
        }
    }
    match texture {
        None => bytes.push(0),
        Some(texture) => {
            bytes.push(1);
            push_u64(&mut bytes, texture.seed);
            push_f64(&mut bytes, texture.rolling_force_gain_n_per_w);
            bytes.push(texture.rolling_target_component as u8);
            push_f64(&mut bytes, texture.band_low_hz);
            push_f64(&mut bytes, texture.band_high_hz);
            push_f64(&mut bytes, texture.reimpact_impulse_n_s);
            bytes.push(texture.reimpact_target_component as u8);
        }
    }
    for value in [
        budget.maximum_total_intervals as u64,
        budget.maximum_chunk_intervals as u64,
        budget.maximum_chunk_events as u64,
        budget.maximum_chunk_spatial_envelopes as u64,
    ] {
        push_u64(&mut bytes, value);
    }
    for value in [
        budget.maximum_abs_source_force_n,
        budget.maximum_abs_source_work_rate_w,
        budget.maximum_abs_generalized_force_n,
        budget.maximum_abs_force_time_measure_n_s,
        budget.maximum_measure_residual_n_s,
        budget.maximum_artistic_texture_envelope_n,
        budget.maximum_artistic_event_impulse_n_s,
        grid.start_time_s,
        grid.end_time_s,
        grid.minimum_interval_duration_s,
        grid.maximum_interval_duration_s,
        grid.nominal_mechanics_timestep_s,
        grid.nominal_source_nyquist_ceiling_hz,
    ] {
        push_f64(&mut bytes, value);
    }
    push_u64(&mut bytes, grid.interval_count as u64);
    bytes.extend_from_slice(b"interval-measures-require-band-limited-resampling-v1");
    hash_domain(MAPPER_IDENTITY_DOMAIN, &bytes)
}

fn chunk_identity(mapper: ContentHash, start: usize, end: usize) -> ContentHash {
    let mut bytes = Vec::with_capacity(48);
    bytes.extend_from_slice(mapper.as_bytes());
    push_u64(&mut bytes, start as u64);
    push_u64(&mut bytes, end as u64);
    hash_domain(CHUNK_IDENTITY_DOMAIN, &bytes)
}

fn texture_root_identity(
    source_identity: ContentHash,
    texture: ArtisticTextureConfig,
) -> ContentHash {
    let mut bytes = Vec::with_capacity(112);
    bytes.extend_from_slice(source_identity.as_bytes());
    push_u32(&mut bytes, AUDIO_EXCITATION_ALGORITHM_VERSION);
    push_u32(&mut bytes, STREAM_SEMANTICS_VERSION);
    push_u64(&mut bytes, texture.seed);
    push_f64(&mut bytes, texture.rolling_force_gain_n_per_w);
    bytes.push(texture.rolling_target_component as u8);
    push_f64(&mut bytes, texture.band_low_hz);
    push_f64(&mut bytes, texture.band_high_hz);
    push_f64(&mut bytes, texture.reimpact_impulse_n_s);
    bytes.push(texture.reimpact_target_component as u8);
    hash_domain(TEXTURE_STREAM_DOMAIN, &bytes)
}

fn interval_texture_identity(root: ContentHash, view: SourceIntervalView<'_>) -> ContentHash {
    let mut bytes = Vec::with_capacity(80);
    bytes.extend_from_slice(root.as_bytes());
    push_u64(&mut bytes, view.first_source_sample_index as u64);
    push_u64(&mut bytes, view.last_source_sample_index as u64);
    push_f64(&mut bytes, view.start_time_s);
    push_f64(&mut bytes, view.end_time_s);
    hash_domain(TEXTURE_STREAM_DOMAIN, &bytes)
}

fn event_texture_identity(
    root: ContentHash,
    event: ControlContactEvent,
    event_index: usize,
) -> ContentHash {
    let mut bytes = Vec::with_capacity(104);
    bytes.extend_from_slice(root.as_bytes());
    push_u64(&mut bytes, event.source_sample_index as u64);
    push_u64(&mut bytes, event_index as u64);
    bytes.push(match event.kind {
        ContactTransitionKind::Opening => 1,
        ContactTransitionKind::Reimpact => 2,
    });
    push_f64(&mut bytes, event.time_s);
    push_f64(&mut bytes, event.bracket_start_s);
    push_f64(&mut bytes, event.bracket_end_s);
    hash_domain(EVENT_STREAM_DOMAIN, &bytes)
}

fn push_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_f64(bytes: &mut Vec<u8>, value: f64) {
    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), AudioExcitationError> {
    cx.checkpoint().map_err(|_| AudioExcitationError::Cancelled)
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
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

    fn is_finite(self) -> bool {
        self.sum.is_finite() && self.correction.is_finite() && self.total().is_finite()
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct ComponentAccumulator {
    disc: CompensatedSum,
    glass_plate: CompensatedSum,
    base_assembly: CompensatedSum,
}

impl ComponentAccumulator {
    const fn new() -> Self {
        Self {
            disc: CompensatedSum::new(),
            glass_plate: CompensatedSum::new(),
            base_assembly: CompensatedSum::new(),
        }
    }

    fn add(&mut self, value: ModalComponentValues) {
        self.disc.add(value.disc);
        self.glass_plate.add(value.glass_plate);
        self.base_assembly.add(value.base_assembly);
    }

    fn total(self) -> ModalComponentValues {
        ModalComponentValues {
            disc: self.disc.total(),
            glass_plate: self.glass_plate.total(),
            base_assembly: self.base_assembly.total(),
        }
    }

    fn is_finite(self) -> bool {
        self.disc.is_finite() && self.glass_plate.is_finite() && self.base_assembly.is_finite()
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct StemAccumulator {
    contact: ComponentAccumulator,
    rolling: ComponentAccumulator,
    base: ComponentAccumulator,
    gas: ComponentAccumulator,
}

impl StemAccumulator {
    const fn new() -> Self {
        Self {
            contact: ComponentAccumulator::new(),
            rolling: ComponentAccumulator::new(),
            base: ComponentAccumulator::new(),
            gas: ComponentAccumulator::new(),
        }
    }

    fn add(&mut self, stems: AudioExcitationStems) {
        self.contact.add(stems.contact);
        self.rolling.add(stems.rolling);
        self.base.add(stems.base);
        self.gas.add(stems.gas);
    }

    fn total(self) -> AudioExcitationStems {
        AudioExcitationStems {
            contact: self.contact.total(),
            rolling: self.rolling.total(),
            base: self.base.total(),
            gas: self.gas.total(),
        }
    }

    fn is_finite(self) -> bool {
        self.contact.is_finite()
            && self.rolling.is_finite()
            && self.base.is_finite()
            && self.gas.is_finite()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_evidence::cinematic_sound::{SoundMode, SoundModeParticipation};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion};

    use crate::{
        AudioExcitationModelInput, DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
        ModalSynthesisBudget, ModalSynthesisModelInput, RenderBaseFrame, RenderBaseModeState,
        RenderChannelAvailability, RenderContactBranch, RenderContactGeometry,
        RenderMassProperties, RenderSupportFeature, RenderTrajectory, RenderTrajectoryAuthority,
        RenderTrajectoryCodecBudget, RenderTrajectoryMetadata, RenderTrajectorySampleInput,
        RenderUnitSystem, RenderWorldFrame, coupled_runner::ChannelOwnership,
    };

    fn with_test_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x4558_4349_5441_5445,
                    kernel_id: 0x4155_4449,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn test_identity(label: &str) -> ContentHash {
        hash_domain(
            "org.frankensim.test.audio-excitation-cancellation.v1",
            label.as_bytes(),
        )
    }

    fn cancellation_fixture(cx: &Cx<'_>) -> EulerRenderTrajectoryArtifact {
        let mass = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.1, 0.1, 0.2)).unwrap();
        // Keep the fixture away from the symmetry-axis coordinate pole: the
        // trajectory contract intentionally refuses an undefined precession
        // angle at zero inclination before the cancellation seam is reached.
        let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.2).unwrap();
        let symmetry_axis_world = orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        let state = RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, 0.05), orientation).unwrap(),
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 1.0),
        )
        .unwrap();
        let input = RenderTrajectorySampleInput {
            interval_start_time_s: 0.0,
            time_s: 1.0,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            center_of_mass_world_m: state.pose().position_world(),
            orientation_body_to_world: orientation.components(),
            linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
            angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
            symmetry_axis_world,
            contact_branch: RenderContactBranch::Closed,
            contact_geometry: Some(RenderContactGeometry {
                point_world_m: Vec3::new(0.01, 0.0, 0.0),
                normal_world: Vec3::new(0.0, 0.0, 1.0),
                support_feature: RenderSupportFeature::ProfileFeature(1),
            }),
            signed_gap_m: 0.0,
            interval_contact_active: true,
            interval_normal_force_n: 0.0,
            contact_transitions: Vec::new(),
            base_mode: Some(RenderBaseModeState {
                displacement_m: 0.0,
                velocity_m_per_s: 0.0,
            }),
            channels: ChannelOwnership::default(),
            mechanical_energy_j: 1.0,
            energy_defect_j: 0.0,
            qois: DerivedEulerQois::from_state(state, mass, 0.0).unwrap(),
            disposition: RenderSampleDisposition::HorizonCensored,
            terminal_event: None,
        };
        let metadata = RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: test_identity("profile"),
            specimen_chart_identity: test_identity("chart"),
            mass_properties: RenderMassProperties {
                identity: test_identity("mass"),
                properties: mass,
            },
            initial_state: state,
            initial_base_mode: input.base_mode.unwrap(),
            base_model_identity: test_identity("base"),
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity: test_identity("model"),
            channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
            configuration_identity: test_identity("configuration"),
            configuration_fingerprint: 0x4155_4449_4f5f_4734,
            timestep_s: 1.0,
            producer_version: "audio-cancellation-test-v1".into(),
            applicability: "transactional cancellation fixture only".into(),
            no_claims: vec!["does not represent calibrated acoustics".into()],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        };
        let trajectory = RenderTrajectory::try_new(metadata, vec![input]).unwrap();
        EulerRenderTrajectoryArtifact::try_from_trajectory(
            test_identity("campaign"),
            trajectory,
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap()
    }

    fn cancellation_modal(cx: &Cx<'_>) -> ModalSynthesisModel {
        ModalSynthesisModel::try_new(
            ModalSynthesisModelInput {
                sample_rate_hz: SOUND_MASTER_SAMPLE_RATE_HZ,
                modes: vec![SoundMode {
                    mode_id: 1,
                    component: SoundModalComponent::Disc,
                    frequency_hz: 800.0,
                    damping_ratio: 0.02,
                    modal_mass_kg: 0.2,
                    source_participation: SoundModeParticipation {
                        disc: 1.0,
                        glass_plate: 0.0,
                        base_assembly: 0.0,
                    },
                    radiation_gain_fs_s_per_m: 0.1,
                    material_identity: test_identity("material"),
                    base_identity: test_identity("modal-base"),
                }],
                budget: ModalSynthesisBudget::reference_film(48_000),
            },
            cx,
        )
        .unwrap()
    }

    #[test]
    fn g4_injected_cancellation_does_not_publish_a_partial_mapping() {
        with_test_cx(|cx| {
            let artifact = cancellation_fixture(cx);
            let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
            let modal = cancellation_modal(cx);
            let mapper = AudioExcitationMapper::try_new(
                &artifact,
                &controls,
                &modal,
                AudioExcitationModelInput {
                    mappings: vec![SoundExcitationControl {
                        channel: SoundExcitationChannel::ContactNormalForce,
                        target_component: SoundModalComponent::Disc,
                        source_scale: 1.0,
                    }],
                    reduction: AudioExcitationReduction::RawIntervals,
                    spatial_policy: ContactParticipationPolicy::DeclaredStatic,
                    artistic_texture: None,
                    budget: AudioExcitationBudget::reference_film(1),
                },
                cx,
            )
            .unwrap();
            let initial = mapper.initial_checkpoint(cx).unwrap();
            let original = initial.clone();
            let mut polls = 0_usize;
            let result = mapper.map_next_chunk_with_checkpoint(
                &initial,
                NonZeroUsize::new(1).unwrap(),
                &mut || {
                    polls += 1;
                    if polls == 7 {
                        Err(AudioExcitationError::Cancelled)
                    } else {
                        Ok(())
                    }
                },
            );
            assert_eq!(
                polls, 7,
                "cancellation must follow completed interval mapping"
            );
            assert_eq!(result, Err(AudioExcitationError::Cancelled));
            assert_eq!(initial, original, "predecessor checkpoint is immutable");

            let replay = mapper
                .map_next_chunk_with_checkpoint(
                    &initial,
                    NonZeroUsize::new(1).unwrap(),
                    &mut || Ok(()),
                )
                .unwrap();
            assert_eq!(replay.intervals.len(), 1);
            assert_eq!(replay.successor.next_interval_index(), 1);
        });
    }

    #[test]
    fn procedural_noise_is_random_access_deterministic_and_bounded() {
        let identity = hash_domain(TEXTURE_STREAM_DOMAIN, b"unit-test");
        for index in [0, 1, 17, u32::MAX as u64 + 1] {
            let left = procedural_texture_unit_sample(identity, index);
            let right = procedural_texture_unit_sample(identity, index);
            assert_eq!(left.to_bits(), right.to_bits());
            assert!((-1.0..1.0).contains(&left));
        }
        assert_ne!(
            procedural_texture_unit_sample(identity, 0).to_bits(),
            procedural_texture_unit_sample(identity, 1).to_bits()
        );
    }
}
