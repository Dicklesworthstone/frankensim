//! Canonical crash-recovery checkpoints for opaque spectral render jobs.
//!
//! This module is a child of `tracer` so it can validate and restore the
//! private row-prefix accumulators without widening their visibility.  It
//! deliberately owns no filesystem or ledger policy: L6 streams the canonical
//! bytes into its transactional artifact store.

use super::*;
use crate::motion::{ShutterConvention, ShutterDistribution};
use fs_blake3::{ContentHash, DomainHasher};

/// Canonical checkpoint schema for uniform and adaptive pending films.
pub const RENDER_CHECKPOINT_SCHEMA_VERSION: u16 = 1;
/// Domain-separated digest over every checkpoint byte before the seal.
pub const RENDER_CHECKPOINT_CONTENT_DOMAIN: &str =
    "org.frankensim.fs-render.progress-checkpoint.v1";
/// Domain-separated identity of every renderer-owned job input that must be
/// identical across checkpoint write and restore.
pub const RENDER_CHECKPOINT_JOB_DOMAIN: &str =
    "org.frankensim.fs-render.progress-checkpoint-job.v1";
/// Domain for the runtime ISA and sorted feature set that constrain bitwise
/// deterministic resume.
pub const RENDER_CHECKPOINT_EXECUTION_ENVIRONMENT_DOMAIN: &str =
    "org.frankensim.fs-render.progress-checkpoint-execution-environment.v1";

const MAGIC: &[u8; 8] = b"FSRCP001";
const SEAL_MAGIC: &[u8; 8] = b"FSRSEAL1";
const SEAL_BYTES: u64 = 8 + 32 + 8;
const STREAM_CHUNK_BYTES: usize = 64 * 1024;
const RESTORE_POLL_PIXELS: u64 = 1_024;
const TILE_RECORD_BYTES: u64 = 8 + 5 * 4;
const UNIFORM_PIXEL_BYTES: u64 = 3 * 8;
const ADAPTIVE_PIXEL_BYTES: u64 = 9 * 8 + 4 + 1;

/// Exact estimator/AOV payload stored in a render checkpoint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum RenderCheckpointKind {
    /// Raw sequential XYZ sums for a uniform-SPP pending render.
    Uniform = 0,
    /// Raw sums plus Welford moments, counts, and stopping decisions.
    Adaptive = 1,
}

/// External identities that close the renderer's otherwise-borrowed assets.
///
/// The renderer validates its own settings, time mode, execution policy, and
/// numeric state.  L6 supplies content identities for the source artifact,
/// configuration, composed scene, frame, job, executable build, and producer
/// claim.  None may use the all-zero sentinel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RenderCheckpointBinding {
    source_artifact_identity: ContentHash,
    source_configuration_identity: ContentHash,
    scene_identity: ContentHash,
    frame_identity: ContentHash,
    render_job_identity: ContentHash,
    producer_build_identity: ContentHash,
    producer_claim_identity: ContentHash,
    generation: u64,
    predecessor_checkpoint: Option<ContentHash>,
}

impl RenderCheckpointBinding {
    /// Construct a complete durable binding.  Generation zero is the only
    /// root; later generations must name a nonzero predecessor.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        source_artifact_identity: ContentHash,
        source_configuration_identity: ContentHash,
        scene_identity: ContentHash,
        frame_identity: ContentHash,
        render_job_identity: ContentHash,
        producer_build_identity: ContentHash,
        producer_claim_identity: ContentHash,
        generation: u64,
        predecessor_checkpoint: Option<ContentHash>,
    ) -> Result<Self, RenderCheckpointError> {
        for (field, identity) in [
            ("source_artifact_identity", source_artifact_identity),
            (
                "source_configuration_identity",
                source_configuration_identity,
            ),
            ("scene_identity", scene_identity),
            ("frame_identity", frame_identity),
            ("render_job_identity", render_job_identity),
            ("producer_build_identity", producer_build_identity),
            ("producer_claim_identity", producer_claim_identity),
        ] {
            if identity.as_bytes().iter().all(|byte| *byte == 0) {
                return Err(RenderCheckpointError::InvalidBinding { field });
            }
        }
        if predecessor_checkpoint
            .is_some_and(|identity| identity.as_bytes().iter().all(|byte| *byte == 0))
        {
            return Err(RenderCheckpointError::InvalidBinding {
                field: "predecessor_checkpoint",
            });
        }
        if (generation == 0) != predecessor_checkpoint.is_none() {
            return Err(RenderCheckpointError::InvalidGenerationChain);
        }
        Ok(Self {
            source_artifact_identity,
            source_configuration_identity,
            scene_identity,
            frame_identity,
            render_job_identity,
            producer_build_identity,
            producer_claim_identity,
            generation,
            predecessor_checkpoint,
        })
    }

    /// Source trajectory/field artifact identity.
    #[must_use]
    pub const fn source_artifact_identity(self) -> ContentHash {
        self.source_artifact_identity
    }

    /// Complete source run/configuration identity.
    #[must_use]
    pub const fn source_configuration_identity(self) -> ContentHash {
        self.source_configuration_identity
    }

    /// Composed render-scene identity.
    #[must_use]
    pub const fn scene_identity(self) -> ContentHash {
        self.scene_identity
    }

    /// Exact frame/shot/segment identity.
    #[must_use]
    pub const fn frame_identity(self) -> ContentHash {
        self.frame_identity
    }

    /// Complete logical render-job identity.
    #[must_use]
    pub const fn render_job_identity(self) -> ContentHash {
        self.render_job_identity
    }

    /// Producer executable/build identity.
    #[must_use]
    pub const fn producer_build_identity(self) -> ContentHash {
        self.producer_build_identity
    }

    /// Producer claim/lease identity retained for later scheduler CAS.
    #[must_use]
    pub const fn producer_claim_identity(self) -> ContentHash {
        self.producer_claim_identity
    }

    /// Monotone checkpoint generation.
    #[must_use]
    pub const fn generation(self) -> u64 {
        self.generation
    }

    /// Immutable prior checkpoint in the generation chain.
    #[must_use]
    pub const fn predecessor_checkpoint(self) -> Option<ContentHash> {
        self.predecessor_checkpoint
    }
}

/// Canonical checkpoint evidence returned only after the final seal is emitted
/// or validated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RenderCheckpointReceipt {
    kind: RenderCheckpointKind,
    content_hash: ContentHash,
    byte_len: u64,
    progress: RenderProgress,
    binding: RenderCheckpointBinding,
}

impl RenderCheckpointReceipt {
    /// Stored estimator/AOV kind.
    #[must_use]
    pub const fn kind(self) -> RenderCheckpointKind {
        self.kind
    }

    /// Domain-separated digest of the canonical body (the footer is excluded).
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.content_hash
    }

    /// Exact sealed artifact length.
    #[must_use]
    pub const fn byte_len(self) -> u64 {
        self.byte_len
    }

    /// Complete-row progress represented by the payload.
    #[must_use]
    pub const fn progress(self) -> RenderProgress {
        self.progress
    }

    /// External identity closure stored in the artifact.
    #[must_use]
    pub const fn binding(self) -> RenderCheckpointBinding {
        self.binding
    }
}

/// Fail-closed canonical-codec or state-validation refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RenderCheckpointError {
    /// One required external identity used the absent all-zero sentinel.
    InvalidBinding {
        /// Absent binding field.
        field: &'static str,
    },
    /// Root/predecessor presence disagreed with the declared generation.
    InvalidGenerationChain,
    /// The exact encoded artifact exceeds the caller's byte budget.
    ByteLimitExceeded {
        /// Exact bytes required by the operation.
        required: u64,
        /// Caller-declared byte ceiling.
        limit: u64,
    },
    /// A fixed-size encoder counter overflowed.
    LengthOverflow,
    /// The bounded 64 KiB streaming scratch allocation failed.
    Allocation,
    /// Cancellation was observed before the next bounded chunk/tile.
    Cancelled,
    /// Input ended before a complete canonical field or seal.
    Truncated,
    /// Magic or reserved bytes did not name this codec.
    InvalidEnvelope,
    /// A different checkpoint schema was supplied.
    UnsupportedSchema {
        /// Schema version found in the envelope.
        found: u16,
    },
    /// Uniform/adaptive payload kind did not match the fresh pending job.
    KindMismatch,
    /// One external identity or generation-chain field differed.
    BindingMismatch {
        /// First mismatched binding field.
        field: &'static str,
    },
    /// One bit-affecting renderer semantic version differed.
    SemanticsMismatch {
        /// First mismatched bit-semantics component.
        field: &'static str,
    },
    /// Settings, time, execution, layout, or policy differed from the freshly
    /// admitted job.
    JobMismatch {
        /// First mismatched renderer job field.
        field: &'static str,
    },
    /// Tile records were reordered, duplicated, malformed, or out of bounds.
    InvalidTileState,
    /// A numeric accumulator, count, or decision could not have been produced
    /// by the pending-render transaction.
    InvalidPixelState,
    /// The domain-separated body digest or stored byte length disagreed.
    IntegrityMismatch,
    /// Canonical payload bytes remained after the derived final record.
    TrailingBytes,
}

impl core::fmt::Display for RenderCheckpointError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidBinding { field } => {
                write!(formatter, "render checkpoint binding {field} is absent")
            }
            Self::InvalidGenerationChain => {
                formatter.write_str("render checkpoint generation/predecessor chain is invalid")
            }
            Self::ByteLimitExceeded { required, limit } => write!(
                formatter,
                "render checkpoint needs {required} bytes, exceeding limit {limit}"
            ),
            Self::LengthOverflow => formatter.write_str("render checkpoint length overflow"),
            Self::Allocation => formatter.write_str("render checkpoint scratch allocation failed"),
            Self::Cancelled => formatter.write_str("render checkpoint operation cancelled"),
            Self::Truncated => formatter.write_str("render checkpoint is truncated"),
            Self::InvalidEnvelope => formatter.write_str("invalid render checkpoint envelope"),
            Self::UnsupportedSchema { found } => {
                write!(formatter, "unsupported render checkpoint schema {found}")
            }
            Self::KindMismatch => formatter.write_str("render checkpoint AOV kind mismatch"),
            Self::BindingMismatch { field } => {
                write!(formatter, "render checkpoint binding mismatch: {field}")
            }
            Self::SemanticsMismatch { field } => {
                write!(formatter, "render checkpoint semantics mismatch: {field}")
            }
            Self::JobMismatch { field } => {
                write!(formatter, "render checkpoint job mismatch: {field}")
            }
            Self::InvalidTileState => formatter.write_str("invalid render checkpoint tile state"),
            Self::InvalidPixelState => formatter.write_str("invalid render checkpoint pixel state"),
            Self::IntegrityMismatch => {
                formatter.write_str("render checkpoint seal or content digest mismatch")
            }
            Self::TrailingBytes => formatter.write_str("render checkpoint has trailing bytes"),
        }
    }
}

impl core::error::Error for RenderCheckpointError {}

/// Streaming emission failure that preserves the sink's native error.
#[derive(Debug)]
pub enum RenderCheckpointWriteError<E> {
    /// Codec, budget, state, or cancellation refusal.
    Checkpoint(RenderCheckpointError),
    /// Caller-provided artifact sink refused one bounded chunk.
    Sink(E),
}

impl<E> From<RenderCheckpointError> for RenderCheckpointWriteError<E> {
    fn from(error: RenderCheckpointError) -> Self {
        Self::Checkpoint(error)
    }
}

impl<E: core::fmt::Display> core::fmt::Display for RenderCheckpointWriteError<E> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Checkpoint(error) => error.fmt(formatter),
            Self::Sink(error) => write!(formatter, "render checkpoint sink refused: {error}"),
        }
    }
}

impl<E: core::fmt::Debug + core::fmt::Display> core::error::Error
    for RenderCheckpointWriteError<E>
{
}

trait BodySink {
    type Error;
    fn write(&mut self, bytes: &[u8]) -> Result<(), Self::Error>;
}

struct CountSink {
    len: u64,
}

impl BodySink for CountSink {
    type Error = RenderCheckpointError;

    fn write(&mut self, bytes: &[u8]) -> Result<(), Self::Error> {
        self.len = self
            .len
            .checked_add(
                u64::try_from(bytes.len()).map_err(|_| RenderCheckpointError::LengthOverflow)?,
            )
            .ok_or(RenderCheckpointError::LengthOverflow)?;
        Ok(())
    }
}

struct StreamSink<'cx, 'scope, 'emit, E, F>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    cx: &'cx Cx<'scope>,
    emit: &'emit mut F,
    buffer: Vec<u8>,
    hasher: DomainHasher,
    len: u64,
}

impl<'cx, 'scope, 'emit, E, F> StreamSink<'cx, 'scope, 'emit, E, F>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    fn try_new(cx: &'cx Cx<'scope>, emit: &'emit mut F) -> Result<Self, RenderCheckpointError> {
        let mut buffer = Vec::new();
        buffer
            .try_reserve_exact(STREAM_CHUNK_BYTES)
            .map_err(|_| RenderCheckpointError::Allocation)?;
        Ok(Self {
            cx,
            emit,
            buffer,
            hasher: DomainHasher::new(RENDER_CHECKPOINT_CONTENT_DOMAIN),
            len: 0,
        })
    }

    fn flush(&mut self) -> Result<(), RenderCheckpointWriteError<E>> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.cx
            .checkpoint()
            .map_err(|_| RenderCheckpointError::Cancelled)?;
        (self.emit)(&self.buffer).map_err(RenderCheckpointWriteError::Sink)?;
        self.buffer.clear();
        Ok(())
    }

    fn finish_body(mut self) -> Result<(ContentHash, u64), RenderCheckpointWriteError<E>> {
        self.flush()?;
        Ok((self.hasher.finalize(), self.len))
    }
}

impl<E, F> BodySink for StreamSink<'_, '_, '_, E, F>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    type Error = RenderCheckpointWriteError<E>;

    fn write(&mut self, mut bytes: &[u8]) -> Result<(), Self::Error> {
        self.hasher.update(bytes);
        self.len = self
            .len
            .checked_add(
                u64::try_from(bytes.len()).map_err(|_| RenderCheckpointError::LengthOverflow)?,
            )
            .ok_or(RenderCheckpointError::LengthOverflow)?;
        while !bytes.is_empty() {
            let available = STREAM_CHUNK_BYTES - self.buffer.len();
            let take = available.min(bytes.len());
            self.buffer.extend_from_slice(&bytes[..take]);
            bytes = &bytes[take..];
            if self.buffer.len() == STREAM_CHUNK_BYTES {
                self.flush()?;
            }
        }
        Ok(())
    }
}

fn put_u8<S: BodySink>(sink: &mut S, value: u8) -> Result<(), S::Error> {
    sink.write(&[value])
}

fn put_u16<S: BodySink>(sink: &mut S, value: u16) -> Result<(), S::Error> {
    sink.write(&value.to_le_bytes())
}

fn put_u32<S: BodySink>(sink: &mut S, value: u32) -> Result<(), S::Error> {
    sink.write(&value.to_le_bytes())
}

fn put_u64<S: BodySink>(sink: &mut S, value: u64) -> Result<(), S::Error> {
    sink.write(&value.to_le_bytes())
}

fn put_f64<S: BodySink>(sink: &mut S, value: f64) -> Result<(), S::Error> {
    put_u64(sink, value.to_bits())
}

fn put_hash<S: BodySink>(sink: &mut S, value: ContentHash) -> Result<(), S::Error> {
    sink.write(value.as_bytes())
}

fn encode_binding<S: BodySink>(
    sink: &mut S,
    binding: RenderCheckpointBinding,
) -> Result<(), S::Error> {
    for identity in [
        binding.source_artifact_identity,
        binding.source_configuration_identity,
        binding.scene_identity,
        binding.frame_identity,
        binding.render_job_identity,
        binding.producer_build_identity,
        binding.producer_claim_identity,
    ] {
        put_hash(sink, identity)?;
    }
    put_u64(sink, binding.generation)?;
    put_u8(sink, u8::from(binding.predecessor_checkpoint.is_some()))?;
    put_hash(
        sink,
        binding
            .predecessor_checkpoint
            .unwrap_or(ContentHash([0; 32])),
    )
}

fn sampler_tag(sampler: Sampler) -> u8 {
    match sampler {
        Sampler::Iid => 0,
        Sampler::OwenSobol => 1,
    }
}

fn strategy_tag(strategy: DirectStrategy) -> u8 {
    match strategy {
        DirectStrategy::NeeOnly => 0,
        DirectStrategy::BsdfOnly => 1,
        DirectStrategy::Mis => 2,
    }
}

fn mode_tag(mode: ExecMode) -> u8 {
    match mode {
        ExecMode::Deterministic => 0,
        ExecMode::Fast => 1,
    }
}

fn convention_tag(convention: ShutterConvention) -> u8 {
    match convention {
        ShutterConvention::Centered => 0,
        ShutterConvention::FrontLoaded => 1,
        ShutterConvention::BackLoaded => 2,
    }
}

fn time_fields(mode: FilmTimeMode) -> (u8, f64, f64, u8, u8, u32, u64, u64) {
    let shutter_fields = |tag, shutter: ShutterInterval, stream, shot| {
        let (distribution, strata) = match shutter.distribution() {
            ShutterDistribution::UniformCounterV1 => (0, 0),
            ShutterDistribution::StratifiedCounterV1 { strata } => (1, strata),
        };
        (
            tag,
            shutter.open_s(),
            shutter.close_s(),
            convention_tag(shutter.convention()),
            distribution,
            strata,
            stream,
            shot,
        )
    };
    match mode {
        FilmTimeMode::Uninitialized => (0, 0.0, 0.0, 0, 0, 0, 0, 0),
        FilmTimeMode::Static => (1, 0.0, 0.0, 0, 0, 0, 0, 0),
        FilmTimeMode::Motion {
            shutter,
            stream_identity,
        } => shutter_fields(2, shutter, stream_identity, 0),
        FilmTimeMode::Cinematic {
            shutter,
            stream_identity,
            shot_id,
        } => shutter_fields(3, shutter, stream_identity, shot_id),
    }
}

fn hash_budget(hasher: &mut DomainHasher, budget: Budget) {
    match budget.deadline {
        None => {
            hasher.update(&[0]);
            hasher.update(&0_u64.to_le_bytes());
        }
        Some(deadline) => {
            hasher.update(&[1]);
            hasher.update(&deadline.as_nanos().to_le_bytes());
        }
    }
    hasher.update(&budget.poll_quota.to_le_bytes());
    match budget.cost_quota {
        None => {
            hasher.update(&[0]);
            hasher.update(&0_u64.to_le_bytes());
        }
        Some(cost) => {
            hasher.update(&[1]);
            hasher.update(&cost.to_le_bytes());
        }
    }
    hasher.update(&[budget.priority]);
}

fn execution_environment_identity() -> ContentHash {
    let probe = fs_substrate::CapabilityProbe::topology_only();
    let mut hasher = DomainHasher::new(RENDER_CHECKPOINT_EXECUTION_ENVIRONMENT_DOMAIN);
    let isa = match probe.isa {
        fs_substrate::Isa::Aarch64Apple => 0,
        fs_substrate::Isa::Aarch64Other => 1,
        fs_substrate::Isa::X86_64 => 2,
        fs_substrate::Isa::Other => 3,
    };
    hasher.update(&[isa]);
    hasher.update(&(probe.features.len() as u64).to_le_bytes());
    for feature in probe.features {
        let feature = feature.as_bytes();
        hasher.update(&(feature.len() as u64).to_le_bytes());
        hasher.update(feature);
    }
    hasher.finalize()
}

#[allow(clippy::too_many_arguments)]
fn checkpoint_job_identity(
    kind: RenderCheckpointKind,
    settings: &Settings,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: &RenderExecutionConfig,
    policy: Option<AdaptiveSamplingConfig>,
) -> ContentHash {
    let mut hasher = DomainHasher::new(RENDER_CHECKPOINT_JOB_DOMAIN);
    hasher.update(&RENDER_CHECKPOINT_SCHEMA_VERSION.to_le_bytes());
    hasher.update(&[kind as u8]);
    for version in [
        TRACER_BIT_SEMANTICS_VERSION,
        MOTION_TRACER_BIT_SEMANTICS_VERSION,
        CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION,
        DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION,
        LIGHTING_TRACER_BIT_SEMANTICS_VERSION,
        policy.map_or(0, |_| ADAPTIVE_SAMPLING_SEMANTICS_VERSION),
    ] {
        hasher.update(&version.to_le_bytes());
    }
    for value in [
        settings.width,
        settings.height,
        settings.spp,
        settings.max_depth,
    ] {
        hasher.update(&value.to_le_bytes());
    }
    hasher.update(&[sampler_tag(settings.sampler)]);
    hasher.update(&[strategy_tag(settings.strategy)]);
    hasher.update(&settings.seed.to_le_bytes());
    let time = time_fields(requested_mode);
    hasher.update(&[time.0]);
    hasher.update(&time.1.to_bits().to_le_bytes());
    hasher.update(&time.2.to_bits().to_le_bytes());
    hasher.update(&[time.3, time.4]);
    hasher.update(&time.5.to_le_bytes());
    hasher.update(&time.6.to_le_bytes());
    hasher.update(&time.7.to_le_bytes());
    hasher.update(&[mode_tag(execution_mode)]);
    hash_budget(&mut hasher, execution_budget);
    hasher.update(execution_environment_identity().as_bytes());
    hasher.update(&execution.tile_width.to_le_bytes());
    hasher.update(&execution.tile_height.to_le_bytes());
    hasher.update(&(execution.workers as u64).to_le_bytes());
    hasher.update(&execution.memory_limit_bytes.to_le_bytes());
    hasher.update(&execution.run_id.0.to_le_bytes());
    hasher.update(&(execution.quantum_weights.len() as u64).to_le_bytes());
    for weight in &execution.quantum_weights {
        hasher.update(&weight.to_le_bytes());
    }
    match policy {
        None => hasher.update(&[0]),
        Some(policy) => {
            hasher.update(&[1]);
            hasher.update(&policy.minimum_samples().to_le_bytes());
            hasher.update(&policy.batch_samples().to_le_bytes());
            hasher.update(&policy.absolute_error().to_bits().to_le_bytes());
            hasher.update(&policy.relative_error().to_bits().to_le_bytes());
            hasher.update(&policy.dark_floor().to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

/// Deterministic identity of every renderer-owned uniform job input.
#[must_use]
pub fn uniform_checkpoint_job_identity(
    settings: &Settings,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: &RenderExecutionConfig,
) -> ContentHash {
    checkpoint_job_identity(
        RenderCheckpointKind::Uniform,
        settings,
        requested_mode,
        execution_mode,
        execution_budget,
        execution,
        None,
    )
}

/// Deterministic identity of every renderer-owned adaptive job input.
#[must_use]
pub fn adaptive_checkpoint_job_identity(
    settings: &Settings,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: &RenderExecutionConfig,
    policy: AdaptiveSamplingConfig,
) -> ContentHash {
    checkpoint_job_identity(
        RenderCheckpointKind::Adaptive,
        settings,
        requested_mode,
        execution_mode,
        execution_budget,
        execution,
        Some(policy),
    )
}

fn encode_budget<S: BodySink>(sink: &mut S, budget: Budget) -> Result<(), S::Error> {
    match budget.deadline {
        None => {
            put_u8(sink, 0)?;
            put_u64(sink, 0)?;
        }
        Some(deadline) => {
            put_u8(sink, 1)?;
            put_u64(sink, deadline.as_nanos())?;
        }
    }
    put_u32(sink, budget.poll_quota)?;
    match budget.cost_quota {
        None => {
            put_u8(sink, 0)?;
            put_u64(sink, 0)?;
        }
        Some(cost) => {
            put_u8(sink, 1)?;
            put_u64(sink, cost)?;
        }
    }
    put_u8(sink, budget.priority)
}

#[allow(clippy::too_many_arguments)]
fn encode_header<S: BodySink>(
    sink: &mut S,
    kind: RenderCheckpointKind,
    binding: RenderCheckpointBinding,
    settings: &Settings,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: &RenderExecutionConfig,
    layout: RenderTileLayout,
    policy: Option<AdaptiveSamplingConfig>,
    attempts: u64,
) -> Result<(), S::Error> {
    sink.write(MAGIC)?;
    put_u16(sink, RENDER_CHECKPOINT_SCHEMA_VERSION)?;
    put_u8(sink, kind as u8)?;
    put_u8(sink, 0)?;
    for version in [
        TRACER_BIT_SEMANTICS_VERSION,
        MOTION_TRACER_BIT_SEMANTICS_VERSION,
        CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION,
        DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION,
        LIGHTING_TRACER_BIT_SEMANTICS_VERSION,
        policy.map_or(0, |_| ADAPTIVE_SAMPLING_SEMANTICS_VERSION),
    ] {
        put_u32(sink, version)?;
    }
    encode_binding(sink, binding)?;
    for value in [
        settings.width,
        settings.height,
        settings.spp,
        settings.max_depth,
    ] {
        put_u32(sink, value)?;
    }
    put_u8(sink, sampler_tag(settings.sampler))?;
    put_u8(sink, strategy_tag(settings.strategy))?;
    put_u64(sink, settings.seed)?;
    let (time_tag, open, close, convention, distribution, strata, stream, shot) =
        time_fields(requested_mode);
    put_u8(sink, time_tag)?;
    put_f64(sink, open)?;
    put_f64(sink, close)?;
    put_u8(sink, convention)?;
    put_u8(sink, distribution)?;
    put_u32(sink, strata)?;
    put_u64(sink, stream)?;
    put_u64(sink, shot)?;
    put_u8(sink, mode_tag(execution_mode))?;
    encode_budget(sink, execution_budget)?;
    put_hash(sink, execution_environment_identity())?;
    put_u32(sink, execution.tile_width)?;
    put_u32(sink, execution.tile_height)?;
    put_u64(sink, execution.workers as u64)?;
    put_u64(sink, execution.memory_limit_bytes)?;
    put_u64(sink, execution.run_id.0)?;
    put_u32(sink, execution.quantum_weights.len() as u32)?;
    for weight in &execution.quantum_weights {
        put_u32(sink, *weight)?;
    }
    put_u32(sink, layout.tiles_x)?;
    put_u32(sink, layout.tiles_y)?;
    put_u64(sink, layout.tile_count)?;
    put_u64(sink, attempts)?;
    put_u64(sink, u64::from(settings.width) * u64::from(settings.height))?;
    match policy {
        None => {
            put_u8(sink, 0)?;
            put_u32(sink, 0)?;
            put_u32(sink, 0)?;
            put_f64(sink, 0.0)?;
            put_f64(sink, 0.0)?;
            put_f64(sink, 0.0)
        }
        Some(policy) => {
            put_u8(sink, 1)?;
            put_u32(sink, policy.minimum_samples())?;
            put_u32(sink, policy.batch_samples())?;
            put_f64(sink, policy.absolute_error())?;
            put_f64(sink, policy.relative_error())?;
            put_f64(sink, policy.dark_floor())
        }
    }
}

fn checkpoint_byte_len(
    kind: RenderCheckpointKind,
    binding: RenderCheckpointBinding,
    settings: &Settings,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: &RenderExecutionConfig,
    layout: RenderTileLayout,
    policy: Option<AdaptiveSamplingConfig>,
    attempts: u64,
) -> Result<u64, RenderCheckpointError> {
    let mut counter = CountSink { len: 0 };
    encode_header(
        &mut counter,
        kind,
        binding,
        settings,
        requested_mode,
        execution_mode,
        execution_budget,
        execution,
        layout,
        policy,
        attempts,
    )?;
    let pixels = u64::from(settings.width) * u64::from(settings.height);
    let pixel_bytes = match kind {
        RenderCheckpointKind::Uniform => UNIFORM_PIXEL_BYTES,
        RenderCheckpointKind::Adaptive => ADAPTIVE_PIXEL_BYTES,
    };
    counter
        .len
        .checked_add(
            layout
                .tile_count
                .checked_mul(TILE_RECORD_BYTES)
                .ok_or(RenderCheckpointError::LengthOverflow)?,
        )
        .and_then(|len| len.checked_add(pixels.checked_mul(pixel_bytes)?))
        .and_then(|len| len.checked_add(SEAL_BYTES))
        .ok_or(RenderCheckpointError::LengthOverflow)
}

fn encode_tile_header<S: BodySink>(
    sink: &mut S,
    tile: u64,
    bounds: RenderTileBounds,
    next_row: u32,
) -> Result<(), S::Error> {
    put_u64(sink, tile)?;
    for value in [bounds.x, bounds.y, bounds.width, bounds.height, next_row] {
        put_u32(sink, value)?;
    }
    Ok(())
}

fn pixel_index(width: u32, x: u32, y: u32) -> usize {
    y as usize * width as usize + x as usize
}

fn encode_uniform_state<E, F>(
    sink: &mut StreamSink<'_, '_, '_, E, F>,
    settings: &Settings,
    layout: RenderTileLayout,
    state: &Mutex<PendingRenderState>,
    next_rows: &[u32],
) -> Result<(), RenderCheckpointWriteError<E>>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    let row_capacity = u64::from(layout.tile_width.min(settings.width))
        .checked_mul(UNIFORM_PIXEL_BYTES)
        .and_then(|bytes| usize::try_from(bytes).ok())
        .ok_or(RenderCheckpointError::LengthOverflow)?;
    let mut row = Vec::new();
    row.try_reserve_exact(row_capacity)
        .map_err(|_| RenderCheckpointError::Allocation)?;
    for tile in 0..layout.tile_count {
        let bounds = layout
            .bounds(tile)
            .expect("validated pending tile remains inside its layout");
        encode_tile_header(sink, tile, bounds, next_rows[tile as usize])?;
        for y in bounds.y..bounds.y + bounds.height {
            row.clear();
            {
                let state = state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                for x in bounds.x..bounds.x + bounds.width {
                    for value in state.xyz[pixel_index(settings.width, x, y)] {
                        row.extend_from_slice(&value.to_bits().to_le_bytes());
                    }
                }
            }
            sink.write(&row)?;
        }
    }
    Ok(())
}

fn encode_adaptive_state<E, F>(
    sink: &mut StreamSink<'_, '_, '_, E, F>,
    settings: &Settings,
    layout: RenderTileLayout,
    state: &Mutex<PendingAdaptiveRenderState>,
    next_rows: &[u32],
) -> Result<(), RenderCheckpointWriteError<E>>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    let row_capacity = u64::from(layout.tile_width.min(settings.width))
        .checked_mul(ADAPTIVE_PIXEL_BYTES)
        .and_then(|bytes| usize::try_from(bytes).ok())
        .ok_or(RenderCheckpointError::LengthOverflow)?;
    let mut row = Vec::new();
    row.try_reserve_exact(row_capacity)
        .map_err(|_| RenderCheckpointError::Allocation)?;
    for tile in 0..layout.tile_count {
        let bounds = layout
            .bounds(tile)
            .expect("validated pending adaptive tile remains inside its layout");
        encode_tile_header(sink, tile, bounds, next_rows[tile as usize])?;
        for y in bounds.y..bounds.y + bounds.height {
            row.clear();
            {
                let state = state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                for x in bounds.x..bounds.x + bounds.width {
                    let pixel = pixel_index(settings.width, x, y);
                    for values in [
                        state.film.xyz[pixel],
                        state.film.mean_xyz[pixel],
                        state.film.m2_xyz[pixel],
                    ] {
                        for value in values {
                            row.extend_from_slice(&value.to_bits().to_le_bytes());
                        }
                    }
                    row.extend_from_slice(&state.film.sample_counts[pixel].to_le_bytes());
                    row.push(state.film.decisions[pixel] as u8);
                }
            }
            sink.write(&row)?;
        }
    }
    Ok(())
}

fn clone_next_rows(rows: &[u32]) -> Result<Vec<u32>, RenderCheckpointError> {
    let mut snapshot = Vec::new();
    snapshot
        .try_reserve_exact(rows.len())
        .map_err(|_| RenderCheckpointError::Allocation)?;
    snapshot.extend_from_slice(rows);
    Ok(snapshot)
}

fn finish_stream<E>(
    cx: &Cx<'_>,
    emit: &mut impl FnMut(&[u8]) -> Result<(), E>,
    digest: ContentHash,
    body_len: u64,
) -> Result<u64, RenderCheckpointWriteError<E>> {
    cx.checkpoint()
        .map_err(|_| RenderCheckpointError::Cancelled)?;
    let total = body_len
        .checked_add(SEAL_BYTES)
        .ok_or(RenderCheckpointError::LengthOverflow)?;
    let mut seal = [0_u8; SEAL_BYTES as usize];
    seal[..8].copy_from_slice(SEAL_MAGIC);
    seal[8..40].copy_from_slice(digest.as_bytes());
    seal[40..].copy_from_slice(&total.to_le_bytes());
    emit(&seal).map_err(RenderCheckpointWriteError::Sink)?;
    // A sink may request cancellation while performing the final durable
    // write. Observe that request before issuing a receipt so L6 drops its
    // uncommitted artifact writer instead of publishing a cancelled seal.
    cx.checkpoint()
        .map_err(|_| RenderCheckpointError::Cancelled)?;
    Ok(total)
}

impl PendingRender<'_> {
    /// Deterministic identity of every renderer-owned input that must remain
    /// identical across a uniform checkpoint write and restore.
    #[must_use]
    pub fn checkpoint_job_identity(&self) -> ContentHash {
        uniform_checkpoint_job_identity(
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
        )
    }

    /// Whether this job borrows the exact cinematic scene and camera objects.
    /// This process-local check lets an L6 adapter prevent cross-wiring an
    /// externally named scene to a pending render from another scene.
    #[must_use]
    pub fn checkpoint_uses_cinematic_sources(
        &self,
        scene: &Scene,
        camera: &AnimatedCamera,
    ) -> bool {
        std::ptr::eq(self.scene, scene)
            && matches!(
                self.camera_path,
                CameraPath::Cinematic {
                    camera: actual,
                    ..
                } if std::ptr::eq(actual, camera)
            )
    }

    /// Stream one canonical, self-sealed uniform checkpoint in bounded chunks.
    /// A sink error leaves publication policy to L6; `ArtifactWriter` rollback
    /// makes that path crash-safe without a filesystem dependency here.
    pub fn write_checkpoint<E>(
        &self,
        binding: RenderCheckpointBinding,
        max_bytes: u64,
        cx: &Cx<'_>,
        mut emit: impl FnMut(&[u8]) -> Result<(), E>,
    ) -> Result<RenderCheckpointReceipt, RenderCheckpointWriteError<E>> {
        cx.checkpoint()
            .map_err(|_| RenderCheckpointError::Cancelled)?;
        if binding.render_job_identity != self.checkpoint_job_identity() {
            return Err(RenderCheckpointError::JobMismatch {
                field: "render_job_identity",
            }
            .into());
        }
        let required = checkpoint_byte_len(
            RenderCheckpointKind::Uniform,
            binding,
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
            self.layout,
            None,
            self.attempts,
        )?;
        if required > max_bytes {
            return Err(RenderCheckpointError::ByteLimitExceeded {
                required,
                limit: max_bytes,
            }
            .into());
        }
        let next_rows = {
            let state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            clone_next_rows(&state.next_row)?
        };
        let progress = pending_progress(self.layout, &next_rows, self.attempts);
        let mut sink = StreamSink::try_new(cx, &mut emit)?;
        encode_header(
            &mut sink,
            RenderCheckpointKind::Uniform,
            binding,
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
            self.layout,
            None,
            self.attempts,
        )?;
        encode_uniform_state(
            &mut sink,
            &self.settings,
            self.layout,
            &self.state,
            &next_rows,
        )?;
        let (digest, body_len) = sink.finish_body()?;
        let byte_len = finish_stream(cx, &mut emit, digest, body_len)?;
        debug_assert_eq!(byte_len, required);
        Ok(RenderCheckpointReceipt {
            kind: RenderCheckpointKind::Uniform,
            content_hash: digest,
            byte_len,
            progress,
            binding,
        })
    }
}

impl PendingAdaptiveRender<'_> {
    /// Deterministic identity of every renderer-owned input that must remain
    /// identical across an adaptive checkpoint write and restore.
    #[must_use]
    pub fn checkpoint_job_identity(&self) -> ContentHash {
        adaptive_checkpoint_job_identity(
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
            self.policy,
        )
    }

    /// Whether this job borrows the exact cinematic scene and camera objects.
    #[must_use]
    pub fn checkpoint_uses_cinematic_sources(
        &self,
        scene: &Scene,
        camera: &AnimatedCamera,
    ) -> bool {
        std::ptr::eq(self.scene, scene)
            && matches!(
                self.camera_path,
                CameraPath::Cinematic {
                    camera: actual,
                    ..
                } if std::ptr::eq(actual, camera)
            )
    }

    /// Stream one canonical adaptive checkpoint including exact Welford AOVs.
    pub fn write_checkpoint<E>(
        &self,
        binding: RenderCheckpointBinding,
        max_bytes: u64,
        cx: &Cx<'_>,
        mut emit: impl FnMut(&[u8]) -> Result<(), E>,
    ) -> Result<RenderCheckpointReceipt, RenderCheckpointWriteError<E>> {
        cx.checkpoint()
            .map_err(|_| RenderCheckpointError::Cancelled)?;
        if binding.render_job_identity != self.checkpoint_job_identity() {
            return Err(RenderCheckpointError::JobMismatch {
                field: "render_job_identity",
            }
            .into());
        }
        let required = checkpoint_byte_len(
            RenderCheckpointKind::Adaptive,
            binding,
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
            self.layout,
            Some(self.policy),
            self.attempts,
        )?;
        if required > max_bytes {
            return Err(RenderCheckpointError::ByteLimitExceeded {
                required,
                limit: max_bytes,
            }
            .into());
        }
        let next_rows = {
            let state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            clone_next_rows(&state.next_row)?
        };
        let progress = pending_progress(self.layout, &next_rows, self.attempts);
        let mut sink = StreamSink::try_new(cx, &mut emit)?;
        encode_header(
            &mut sink,
            RenderCheckpointKind::Adaptive,
            binding,
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
            self.layout,
            Some(self.policy),
            self.attempts,
        )?;
        encode_adaptive_state(
            &mut sink,
            &self.settings,
            self.layout,
            &self.state,
            &next_rows,
        )?;
        let (digest, body_len) = sink.finish_body()?;
        let byte_len = finish_stream(cx, &mut emit, digest, body_len)?;
        debug_assert_eq!(byte_len, required);
        Ok(RenderCheckpointReceipt {
            kind: RenderCheckpointKind::Adaptive,
            content_hash: digest,
            byte_len,
            progress,
            binding,
        })
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, at: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], RenderCheckpointError> {
        let end = self
            .at
            .checked_add(count)
            .ok_or(RenderCheckpointError::Truncated)?;
        let out = self
            .bytes
            .get(self.at..end)
            .ok_or(RenderCheckpointError::Truncated)?;
        self.at = end;
        Ok(out)
    }

    fn u8(&mut self) -> Result<u8, RenderCheckpointError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, RenderCheckpointError> {
        Ok(u16::from_le_bytes(
            self.take(2)?.try_into().expect("exact reader width"),
        ))
    }

    fn u32(&mut self) -> Result<u32, RenderCheckpointError> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("exact reader width"),
        ))
    }

    fn u64(&mut self) -> Result<u64, RenderCheckpointError> {
        Ok(u64::from_le_bytes(
            self.take(8)?.try_into().expect("exact reader width"),
        ))
    }

    fn f64(&mut self) -> Result<f64, RenderCheckpointError> {
        Ok(f64::from_bits(self.u64()?))
    }

    fn hash(&mut self) -> Result<ContentHash, RenderCheckpointError> {
        let mut bytes = [0_u8; 32];
        bytes.copy_from_slice(self.take(32)?);
        Ok(ContentHash(bytes))
    }

    fn finish(self) -> Result<(), RenderCheckpointError> {
        if self.at == self.bytes.len() {
            Ok(())
        } else {
            Err(RenderCheckpointError::TrailingBytes)
        }
    }
}

fn decode_binding(
    reader: &mut Reader<'_>,
) -> Result<RenderCheckpointBinding, RenderCheckpointError> {
    let source_artifact_identity = reader.hash()?;
    let source_configuration_identity = reader.hash()?;
    let scene_identity = reader.hash()?;
    let frame_identity = reader.hash()?;
    let render_job_identity = reader.hash()?;
    let producer_build_identity = reader.hash()?;
    let producer_claim_identity = reader.hash()?;
    let generation = reader.u64()?;
    let predecessor_tag = reader.u8()?;
    let predecessor = reader.hash()?;
    let predecessor_checkpoint = match predecessor_tag {
        0 if predecessor.as_bytes().iter().all(|byte| *byte == 0) => None,
        1 => Some(predecessor),
        _ => return Err(RenderCheckpointError::InvalidGenerationChain),
    };
    RenderCheckpointBinding::try_new(
        source_artifact_identity,
        source_configuration_identity,
        scene_identity,
        frame_identity,
        render_job_identity,
        producer_build_identity,
        producer_claim_identity,
        generation,
        predecessor_checkpoint,
    )
}

fn require_binding(
    actual: RenderCheckpointBinding,
    expected: RenderCheckpointBinding,
) -> Result<(), RenderCheckpointError> {
    for (field, equal) in [
        (
            "source_artifact_identity",
            actual.source_artifact_identity == expected.source_artifact_identity,
        ),
        (
            "source_configuration_identity",
            actual.source_configuration_identity == expected.source_configuration_identity,
        ),
        (
            "scene_identity",
            actual.scene_identity == expected.scene_identity,
        ),
        (
            "frame_identity",
            actual.frame_identity == expected.frame_identity,
        ),
        (
            "render_job_identity",
            actual.render_job_identity == expected.render_job_identity,
        ),
        (
            "producer_build_identity",
            actual.producer_build_identity == expected.producer_build_identity,
        ),
        (
            "producer_claim_identity",
            actual.producer_claim_identity == expected.producer_claim_identity,
        ),
        ("generation", actual.generation == expected.generation),
        (
            "predecessor_checkpoint",
            actual.predecessor_checkpoint == expected.predecessor_checkpoint,
        ),
    ] {
        if !equal {
            return Err(RenderCheckpointError::BindingMismatch { field });
        }
    }
    Ok(())
}

fn require_u32(
    reader: &mut Reader<'_>,
    expected: u32,
    field: &'static str,
) -> Result<(), RenderCheckpointError> {
    if reader.u32()? == expected {
        Ok(())
    } else {
        Err(RenderCheckpointError::JobMismatch { field })
    }
}

fn require_u64(
    reader: &mut Reader<'_>,
    expected: u64,
    field: &'static str,
) -> Result<(), RenderCheckpointError> {
    if reader.u64()? == expected {
        Ok(())
    } else {
        Err(RenderCheckpointError::JobMismatch { field })
    }
}

fn require_budget(reader: &mut Reader<'_>, expected: Budget) -> Result<(), RenderCheckpointError> {
    let deadline_tag = reader.u8()?;
    let deadline_nanos = reader.u64()?;
    let poll_quota = reader.u32()?;
    let cost_tag = reader.u8()?;
    let cost_quota = reader.u64()?;
    let priority = reader.u8()?;
    if !matches!(deadline_tag, 0 | 1)
        || !matches!(cost_tag, 0 | 1)
        || (deadline_tag == 0 && deadline_nanos != 0)
        || (cost_tag == 0 && cost_quota != 0)
    {
        return Err(RenderCheckpointError::InvalidEnvelope);
    }
    let expected_deadline = expected
        .deadline
        .map_or((0, 0), |deadline| (1, deadline.as_nanos()));
    let expected_cost = expected.cost_quota.map_or((0, 0), |cost| (1, cost));
    if (deadline_tag, deadline_nanos) != expected_deadline
        || poll_quota != expected.poll_quota
        || (cost_tag, cost_quota) != expected_cost
        || priority != expected.priority
    {
        return Err(RenderCheckpointError::JobMismatch {
            field: "execution_budget",
        });
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn decode_header(
    reader: &mut Reader<'_>,
    expected_kind: RenderCheckpointKind,
    expected_binding: RenderCheckpointBinding,
    settings: &Settings,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: &RenderExecutionConfig,
    layout: RenderTileLayout,
    policy: Option<AdaptiveSamplingConfig>,
) -> Result<u64, RenderCheckpointError> {
    if reader.take(8)? != MAGIC.as_slice() {
        return Err(RenderCheckpointError::InvalidEnvelope);
    }
    let schema = reader.u16()?;
    if schema != RENDER_CHECKPOINT_SCHEMA_VERSION {
        return Err(RenderCheckpointError::UnsupportedSchema { found: schema });
    }
    if reader.u8()? != expected_kind as u8 {
        return Err(RenderCheckpointError::KindMismatch);
    }
    if reader.u8()? != 0 {
        return Err(RenderCheckpointError::InvalidEnvelope);
    }
    for (field, expected) in [
        ("tracer", TRACER_BIT_SEMANTICS_VERSION),
        ("motion", MOTION_TRACER_BIT_SEMANTICS_VERSION),
        (
            "cinematic_camera",
            CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION,
        ),
        ("dielectric", DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION),
        ("lighting", LIGHTING_TRACER_BIT_SEMANTICS_VERSION),
        (
            "adaptive",
            policy.map_or(0, |_| ADAPTIVE_SAMPLING_SEMANTICS_VERSION),
        ),
    ] {
        if reader.u32()? != expected {
            return Err(RenderCheckpointError::SemanticsMismatch { field });
        }
    }
    require_binding(decode_binding(reader)?, expected_binding)?;
    for (field, expected) in [
        ("width", settings.width),
        ("height", settings.height),
        ("spp", settings.spp),
        ("max_depth", settings.max_depth),
    ] {
        require_u32(reader, expected, field)?;
    }
    if reader.u8()? != sampler_tag(settings.sampler) {
        return Err(RenderCheckpointError::JobMismatch { field: "sampler" });
    }
    if reader.u8()? != strategy_tag(settings.strategy) {
        return Err(RenderCheckpointError::JobMismatch { field: "strategy" });
    }
    require_u64(reader, settings.seed, "seed")?;
    let expected_time = time_fields(requested_mode);
    let actual_time = (
        reader.u8()?,
        reader.f64()?,
        reader.f64()?,
        reader.u8()?,
        reader.u8()?,
        reader.u32()?,
        reader.u64()?,
        reader.u64()?,
    );
    if actual_time.0 != expected_time.0
        || actual_time.1.to_bits() != expected_time.1.to_bits()
        || actual_time.2.to_bits() != expected_time.2.to_bits()
        || actual_time.3 != expected_time.3
        || actual_time.4 != expected_time.4
        || actual_time.5 != expected_time.5
        || actual_time.6 != expected_time.6
        || actual_time.7 != expected_time.7
    {
        return Err(RenderCheckpointError::JobMismatch { field: "time_mode" });
    }
    if reader.u8()? != mode_tag(execution_mode) {
        return Err(RenderCheckpointError::JobMismatch {
            field: "execution_mode",
        });
    }
    require_budget(reader, execution_budget)?;
    if reader.hash()? != execution_environment_identity() {
        return Err(RenderCheckpointError::JobMismatch {
            field: "execution_environment",
        });
    }
    require_u32(reader, execution.tile_width, "tile_width")?;
    require_u32(reader, execution.tile_height, "tile_height")?;
    require_u64(reader, execution.workers as u64, "workers")?;
    require_u64(reader, execution.memory_limit_bytes, "memory_limit_bytes")?;
    require_u64(reader, execution.run_id.0, "run_id")?;
    let weight_count = reader.u32()? as usize;
    if weight_count != execution.quantum_weights.len() {
        return Err(RenderCheckpointError::JobMismatch {
            field: "quantum_weights",
        });
    }
    for expected in &execution.quantum_weights {
        require_u32(reader, *expected, "quantum_weights")?;
    }
    require_u32(reader, layout.tiles_x, "tiles_x")?;
    require_u32(reader, layout.tiles_y, "tiles_y")?;
    require_u64(reader, layout.tile_count, "tile_count")?;
    let attempts = reader.u64()?;
    require_u64(
        reader,
        u64::from(settings.width) * u64::from(settings.height),
        "pixel_count",
    )?;
    let policy_tag = reader.u8()?;
    let decoded_policy = (
        reader.u32()?,
        reader.u32()?,
        reader.f64()?,
        reader.f64()?,
        reader.f64()?,
    );
    match policy {
        None if policy_tag == 0
            && decoded_policy.0 == 0
            && decoded_policy.1 == 0
            && decoded_policy.2.to_bits() == 0
            && decoded_policy.3.to_bits() == 0
            && decoded_policy.4.to_bits() == 0 =>
        {
            Ok(attempts)
        }
        Some(expected)
            if policy_tag == 1
                && decoded_policy.0 == expected.minimum_samples()
                && decoded_policy.1 == expected.batch_samples()
                && decoded_policy.2.to_bits() == expected.absolute_error().to_bits()
                && decoded_policy.3.to_bits() == expected.relative_error().to_bits()
                && decoded_policy.4.to_bits() == expected.dark_floor().to_bits() =>
        {
            Ok(attempts)
        }
        _ => Err(RenderCheckpointError::JobMismatch {
            field: "adaptive_policy",
        }),
    }
}

fn open_checkpoint<'a>(
    bytes: &'a [u8],
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<(Reader<'a>, ContentHash), RenderCheckpointError> {
    let byte_len = u64::try_from(bytes.len()).map_err(|_| RenderCheckpointError::LengthOverflow)?;
    if byte_len > max_bytes {
        return Err(RenderCheckpointError::ByteLimitExceeded {
            required: byte_len,
            limit: max_bytes,
        });
    }
    let seal_start = bytes
        .len()
        .checked_sub(SEAL_BYTES as usize)
        .ok_or(RenderCheckpointError::Truncated)?;
    let (body, seal) = bytes.split_at(seal_start);
    if &seal[..8] != SEAL_MAGIC {
        return Err(RenderCheckpointError::IntegrityMismatch);
    }
    let stored_hash =
        ContentHash::from_slice(&seal[8..40]).ok_or(RenderCheckpointError::IntegrityMismatch)?;
    let stored_len = u64::from_le_bytes(
        seal[40..48]
            .try_into()
            .map_err(|_| RenderCheckpointError::Truncated)?,
    );
    if stored_len != byte_len {
        return Err(RenderCheckpointError::IntegrityMismatch);
    }
    let mut hasher = DomainHasher::new(RENDER_CHECKPOINT_CONTENT_DOMAIN);
    for chunk in body.chunks(STREAM_CHUNK_BYTES) {
        cx.checkpoint()
            .map_err(|_| RenderCheckpointError::Cancelled)?;
        hasher.update(chunk);
    }
    if hasher.finalize() != stored_hash {
        return Err(RenderCheckpointError::IntegrityMismatch);
    }
    Ok((Reader::new(body), stored_hash))
}

fn decode_tile_header(
    reader: &mut Reader<'_>,
    tile: u64,
    expected: RenderTileBounds,
) -> Result<u32, RenderCheckpointError> {
    let actual_tile = reader.u64()?;
    let actual = RenderTileBounds {
        x: reader.u32()?,
        y: reader.u32()?,
        width: reader.u32()?,
        height: reader.u32()?,
    };
    let next_row = reader.u32()?;
    if actual_tile != tile || actual != expected || next_row > expected.height {
        return Err(RenderCheckpointError::InvalidTileState);
    }
    Ok(next_row)
}

fn canonical_zero(values: impl IntoIterator<Item = f64>) -> bool {
    values.into_iter().all(|value| value.to_bits() == 0)
}

fn canonical_nonnegative_moment(value: f64) -> bool {
    value.is_finite() && value >= 0.0 && (value != 0.0 || value.to_bits() == 0)
}

fn adaptive_sum_matches_mean(sum: f64, mean: f64, samples: u32) -> bool {
    let samples = f64::from(samples);
    let expected = mean * samples;
    if !sum.is_finite() || !mean.is_finite() || !expected.is_finite() {
        return false;
    }
    // The sequential sum and Welford mean have different rounding paths.
    // Render radiance is finite, so a conservative O(n*eps) envelope accepts
    // reachable rounding drift while rejecting independently forged AOVs.
    let scale = sum.abs().max(expected.abs()).max(1.0);
    let tolerance = scale * f64::EPSILON * 64.0 * samples.max(1.0);
    (sum - expected).abs() <= tolerance
}

impl<'assets> PendingRender<'assets> {
    /// Restore a canonical uniform checkpoint into this freshly re-admitted
    /// scene/camera/settings job.  The input is whole and caller-bounded in v1;
    /// only emission is streaming.  On any failure `self` is consumed and no
    /// partially restored job can escape.
    pub fn restore_checkpoint(
        mut self,
        expected: RenderCheckpointBinding,
        bytes: &[u8],
        max_bytes: u64,
        cx: &Cx<'_>,
    ) -> Result<(Self, RenderCheckpointReceipt), RenderCheckpointError> {
        cx.checkpoint()
            .map_err(|_| RenderCheckpointError::Cancelled)?;
        let (mut reader, digest) = open_checkpoint(bytes, max_bytes, cx)?;
        let attempts = decode_header(
            &mut reader,
            RenderCheckpointKind::Uniform,
            expected,
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
            self.layout,
            None,
        )?;
        if expected.render_job_identity != self.checkpoint_job_identity() {
            return Err(RenderCheckpointError::JobMismatch {
                field: "render_job_identity",
            });
        }
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let mut ordinal = 0_u64;
        for tile in 0..self.layout.tile_count {
            cx.checkpoint()
                .map_err(|_| RenderCheckpointError::Cancelled)?;
            let bounds = self
                .layout
                .bounds(tile)
                .ok_or(RenderCheckpointError::InvalidTileState)?;
            let next_row = decode_tile_header(&mut reader, tile, bounds)?;
            state.next_row[tile as usize] = next_row;
            for local_y in 0..bounds.height {
                for x in bounds.x..bounds.x + bounds.width {
                    if ordinal.is_multiple_of(RESTORE_POLL_PIXELS) {
                        cx.checkpoint()
                            .map_err(|_| RenderCheckpointError::Cancelled)?;
                    }
                    ordinal += 1;
                    let xyz = [reader.f64()?, reader.f64()?, reader.f64()?];
                    if xyz.iter().any(|value| !value.is_finite())
                        || (local_y >= next_row && !canonical_zero(xyz))
                    {
                        return Err(RenderCheckpointError::InvalidPixelState);
                    }
                    let y = bounds.y + local_y;
                    state.xyz[pixel_index(self.settings.width, x, y)] = xyz;
                }
            }
        }
        reader.finish()?;
        self.attempts = attempts;
        let progress = pending_progress(self.layout, &state.next_row, attempts);
        drop(state);
        let receipt = RenderCheckpointReceipt {
            kind: RenderCheckpointKind::Uniform,
            content_hash: digest,
            byte_len: u64::try_from(bytes.len())
                .map_err(|_| RenderCheckpointError::LengthOverflow)?,
            progress,
            binding: expected,
        };
        Ok((self, receipt))
    }
}

impl<'assets> PendingAdaptiveRender<'assets> {
    /// Restore exact adaptive sums/moments/counts/decisions into a freshly
    /// re-admitted job.  Invalid or partially decoded state is never returned.
    pub fn restore_checkpoint(
        mut self,
        expected: RenderCheckpointBinding,
        bytes: &[u8],
        max_bytes: u64,
        cx: &Cx<'_>,
    ) -> Result<(Self, RenderCheckpointReceipt), RenderCheckpointError> {
        cx.checkpoint()
            .map_err(|_| RenderCheckpointError::Cancelled)?;
        let (mut reader, digest) = open_checkpoint(bytes, max_bytes, cx)?;
        let attempts = decode_header(
            &mut reader,
            RenderCheckpointKind::Adaptive,
            expected,
            &self.settings,
            self.requested_mode,
            self.execution_mode,
            self.execution_budget,
            &self.execution,
            self.layout,
            Some(self.policy),
        )?;
        if expected.render_job_identity != self.checkpoint_job_identity() {
            return Err(RenderCheckpointError::JobMismatch {
                field: "render_job_identity",
            });
        }
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let mut ordinal = 0_u64;
        for tile in 0..self.layout.tile_count {
            cx.checkpoint()
                .map_err(|_| RenderCheckpointError::Cancelled)?;
            let bounds = self
                .layout
                .bounds(tile)
                .ok_or(RenderCheckpointError::InvalidTileState)?;
            let next_row = decode_tile_header(&mut reader, tile, bounds)?;
            state.next_row[tile as usize] = next_row;
            for local_y in 0..bounds.height {
                for x in bounds.x..bounds.x + bounds.width {
                    if ordinal.is_multiple_of(RESTORE_POLL_PIXELS) {
                        cx.checkpoint()
                            .map_err(|_| RenderCheckpointError::Cancelled)?;
                    }
                    ordinal += 1;
                    let mut vectors = [[0.0_f64; 3]; 3];
                    for vector in &mut vectors {
                        for value in vector {
                            *value = reader.f64()?;
                        }
                    }
                    let samples = reader.u32()?;
                    let decision = match reader.u8()? {
                        0 => AdaptiveDecision::ErrorThreshold,
                        1 => AdaptiveDecision::MaximumSamples,
                        _ => return Err(RenderCheckpointError::InvalidPixelState),
                    };
                    let committed = local_y < next_row;
                    if vectors.iter().flatten().any(|value| !value.is_finite())
                        || vectors[2]
                            .iter()
                            .any(|value| !canonical_nonnegative_moment(*value))
                    {
                        return Err(RenderCheckpointError::InvalidPixelState);
                    }
                    if committed {
                        let accumulator = AdaptivePixelAccumulator {
                            sum_xyz: vectors[0],
                            mean_xyz: vectors[1],
                            m2_xyz: vectors[2],
                            samples,
                            decision: Some(decision),
                        };
                        if samples < self.policy.minimum_samples()
                            || samples > self.settings.spp
                            || (0..3).any(|channel| {
                                !adaptive_sum_matches_mean(
                                    vectors[0][channel],
                                    vectors[1][channel],
                                    samples,
                                )
                            })
                            || accumulator.decision(self.policy, self.settings.spp)
                                != Some(decision)
                        {
                            return Err(RenderCheckpointError::InvalidPixelState);
                        }
                    } else if samples != 0
                        || decision != AdaptiveDecision::MaximumSamples
                        || !vectors.iter().copied().all(canonical_zero)
                    {
                        return Err(RenderCheckpointError::InvalidPixelState);
                    }
                    let y = bounds.y + local_y;
                    let pixel = pixel_index(self.settings.width, x, y);
                    state.film.xyz[pixel] = vectors[0];
                    state.film.mean_xyz[pixel] = vectors[1];
                    state.film.m2_xyz[pixel] = vectors[2];
                    state.film.sample_counts[pixel] = samples;
                    state.film.decisions[pixel] = decision;
                }
            }
        }
        reader.finish()?;
        self.attempts = attempts;
        let progress = pending_progress(self.layout, &state.next_row, attempts);
        drop(state);
        let receipt = RenderCheckpointReceipt {
            kind: RenderCheckpointKind::Adaptive,
            content_hash: digest,
            byte_len: u64::try_from(bytes.len())
                .map_err(|_| RenderCheckpointError::LengthOverflow)?,
            progress,
            binding: expected,
        };
        Ok((self, receipt))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lighting::EnvironmentMap;
    use asupersync::types::Budget;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_blake3::hash_domain;
    use fs_exec::{CancelGate, StreamKey};
    use std::convert::Infallible;

    const MAX_CHECKPOINT_BYTES: u64 = 1 << 20;

    fn with_gate_cx<R>(gate: &CancelGate, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let arenas = ArenaPool::new(ArenaConfig::default());
        arenas.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 0x6368_6563_6b70_6f69,
                    kernel_id: 0x7273_746f_7265_5f76,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn scene() -> Scene {
        Scene {
            primitives: Vec::new(),
            lights: Vec::new(),
            environment: Some(
                EnvironmentMap::try_from_linear_srgb(
                    4,
                    2,
                    vec![
                        [0.3, 0.4, 0.5],
                        [0.7, 0.2, 0.1],
                        [0.2, 0.8, 0.4],
                        [0.9, 0.7, 0.3],
                        [0.4, 0.1, 0.8],
                        [0.6, 0.6, 0.6],
                        [0.1, 0.3, 0.9],
                        [0.8, 0.5, 0.2],
                    ],
                    0.37,
                )
                .expect("valid deterministic checkpoint environment"),
            ),
            camera: Camera {
                eye: Point3::new(0.0, 0.0, 0.0),
                forward: Vec3::new(1.0, 0.0, 0.0),
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 0.7,
            },
        }
    }

    fn test_settings(seed: u64) -> Settings {
        Settings {
            width: 3,
            height: 3,
            spp: 2,
            max_depth: 2,
            sampler: Sampler::Iid,
            strategy: DirectStrategy::Mis,
            seed,
        }
    }

    fn execution() -> RenderExecutionConfig {
        RenderExecutionConfig::try_new(2, 2, 1, 8 << 20, RunId(0x6368_6563_6b70_0001))
            .expect("valid checkpoint execution policy")
    }

    fn identity(label: &str) -> ContentHash {
        hash_domain(
            "org.frankensim.fs-render.checkpoint-test.v1",
            label.as_bytes(),
        )
    }

    fn binding_with_job(render_job_identity: ContentHash) -> RenderCheckpointBinding {
        RenderCheckpointBinding::try_new(
            identity("source"),
            identity("configuration"),
            identity("scene"),
            identity("frame"),
            render_job_identity,
            identity("build"),
            identity("claim"),
            0,
            None,
        )
        .expect("valid root binding")
    }

    fn binding() -> RenderCheckpointBinding {
        binding_with_job(identity("job"))
    }

    fn assert_film_bits_eq(actual: &Film, expected: &Film, context: &str) {
        assert_eq!(
            (
                actual.width,
                actual.height,
                actual.spp_done,
                actual.time_mode
            ),
            (
                expected.width,
                expected.height,
                expected.spp_done,
                expected.time_mode
            ),
            "{context}: film metadata"
        );
        for (pixel, (actual, expected)) in actual.xyz.iter().zip(&expected.xyz).enumerate() {
            for channel in 0..3 {
                assert_eq!(
                    actual[channel].to_bits(),
                    expected[channel].to_bits(),
                    "{context}: pixel={pixel} channel={channel}"
                );
            }
        }
    }

    fn assert_adaptive_bits_eq(actual: &AdaptiveFilm, expected: &AdaptiveFilm, context: &str) {
        assert_eq!(
            actual.sample_counts(),
            expected.sample_counts(),
            "{context}"
        );
        assert_eq!(actual.decisions(), expected.decisions(), "{context}");
        for pixel in 0..actual.xyz_sums().len() {
            for (label, actual, expected) in [
                ("sum", actual.xyz_sums()[pixel], expected.xyz_sums()[pixel]),
                (
                    "mean",
                    actual.running_means_xyz()[pixel],
                    expected.running_means_xyz()[pixel],
                ),
                ("m2", actual.m2_xyz()[pixel], expected.m2_xyz()[pixel]),
            ] {
                for channel in 0..3 {
                    assert_eq!(
                        actual[channel].to_bits(),
                        expected[channel].to_bits(),
                        "{context}: {label} pixel={pixel} channel={channel}"
                    );
                }
            }
        }
    }

    fn reseal_body(mut body: Vec<u8>) -> Vec<u8> {
        let mut hasher = DomainHasher::new(RENDER_CHECKPOINT_CONTENT_DOMAIN);
        hasher.update(&body);
        let digest = hasher.finalize();
        let total = u64::try_from(body.len())
            .expect("test body length fits u64")
            .checked_add(SEAL_BYTES)
            .expect("test sealed length fits u64");
        body.extend_from_slice(SEAL_MAGIC);
        body.extend_from_slice(digest.as_bytes());
        body.extend_from_slice(&total.to_le_bytes());
        body
    }

    fn uniform_state_offset(
        pending: &PendingRender<'_>,
        binding: RenderCheckpointBinding,
        body: &[u8],
    ) -> usize {
        let mut reader = Reader::new(body);
        decode_header(
            &mut reader,
            RenderCheckpointKind::Uniform,
            binding,
            &pending.settings,
            pending.requested_mode,
            pending.execution_mode,
            pending.execution_budget,
            &pending.execution,
            pending.layout,
            None,
        )
        .expect("decode valid uniform test header");
        reader.at
    }

    fn adaptive_state_offset(
        pending: &PendingAdaptiveRender<'_>,
        binding: RenderCheckpointBinding,
        body: &[u8],
    ) -> usize {
        let mut reader = Reader::new(body);
        decode_header(
            &mut reader,
            RenderCheckpointKind::Adaptive,
            binding,
            &pending.settings,
            pending.requested_mode,
            pending.execution_mode,
            pending.execution_budget,
            &pending.execution,
            pending.layout,
            Some(pending.policy),
        )
        .expect("decode valid adaptive test header");
        reader.at
    }

    #[test]
    fn g0_binding_generation_chain_is_fail_closed() {
        let zero = ContentHash([0; 32]);
        let valid = binding();
        assert_eq!(valid.generation(), 0);
        assert_eq!(valid.predecessor_checkpoint(), None);
        assert!(matches!(
            RenderCheckpointBinding::try_new(
                zero,
                identity("configuration"),
                identity("scene"),
                identity("frame"),
                identity("job"),
                identity("build"),
                identity("claim"),
                0,
                None,
            ),
            Err(RenderCheckpointError::InvalidBinding {
                field: "source_artifact_identity"
            })
        ));
        assert!(matches!(
            RenderCheckpointBinding::try_new(
                identity("source"),
                identity("configuration"),
                identity("scene"),
                identity("frame"),
                identity("job"),
                identity("build"),
                identity("claim"),
                0,
                Some(identity("prior")),
            ),
            Err(RenderCheckpointError::InvalidGenerationChain)
        ));
        assert!(matches!(
            RenderCheckpointBinding::try_new(
                identity("source"),
                identity("configuration"),
                identity("scene"),
                identity("frame"),
                identity("job"),
                identity("build"),
                identity("claim"),
                1,
                None,
            ),
            Err(RenderCheckpointError::InvalidGenerationChain)
        ));
    }

    #[test]
    fn g3_uniform_codec_round_trips_and_refuses_damage_stale_jobs_and_short_budgets() {
        let gate = CancelGate::new_clock_free();
        with_gate_cx(&gate, |cx| {
            let scene = scene();
            let settings = test_settings(0x6368_6563_6b70_0101);
            let execution = execution();
            let pending = PendingRender::begin_static(&scene, cx, settings, execution.clone())
                .expect("admit uniform checkpoint job");
            let binding = binding_with_job(pending.checkpoint_job_identity());
            let mut bytes = Vec::new();
            let receipt = pending
                .write_checkpoint::<Infallible>(binding, MAX_CHECKPOINT_BYTES, cx, |chunk| {
                    bytes.extend_from_slice(chunk);
                    Ok(())
                })
                .expect("encode uniform checkpoint");
            assert_eq!(receipt.kind(), RenderCheckpointKind::Uniform);
            assert_eq!(receipt.byte_len(), bytes.len() as u64);
            assert_eq!(receipt.progress().committed_tile_rows, 0);

            for prefix in 0..bytes.len() {
                assert!(
                    open_checkpoint(&bytes[..prefix], MAX_CHECKPOINT_BYTES, cx).is_err(),
                    "truncated prefix {prefix}/{} unexpectedly opened",
                    bytes.len()
                );
            }
            for offset in [0, bytes.len() / 2, bytes.len() - 1] {
                let mut damaged = bytes.clone();
                damaged[offset] ^= 0x01;
                assert!(
                    matches!(
                        open_checkpoint(&damaged, MAX_CHECKPOINT_BYTES, cx),
                        Err(RenderCheckpointError::IntegrityMismatch)
                    ),
                    "bit corruption at byte {offset} escaped the seal"
                );
            }

            let restored_seed =
                PendingRender::begin_static(&scene, cx, settings, execution.clone())
                    .expect("admit uniform restore target");
            let (restored, restored_receipt) = restored_seed
                .restore_checkpoint(binding, &bytes, MAX_CHECKPOINT_BYTES, cx)
                .expect("restore uniform checkpoint");
            assert_eq!(restored_receipt, receipt);
            let restored = restored.resume(cx).expect("finish restored uniform");
            let reference = PendingRender::begin_static(&scene, cx, settings, execution.clone())
                .expect("admit uniform reference")
                .resume(cx)
                .expect("finish uniform reference");
            assert_film_bits_eq(&restored.film, &reference.film, "uniform checkpoint");

            let mut changed_binding = binding;
            changed_binding.producer_claim_identity = identity("different-claim");
            let changed_target =
                PendingRender::begin_static(&scene, cx, settings, execution.clone())
                    .expect("admit binding mismatch target");
            assert!(matches!(
                changed_target.restore_checkpoint(
                    changed_binding,
                    &bytes,
                    MAX_CHECKPOINT_BYTES,
                    cx,
                ),
                Err(RenderCheckpointError::BindingMismatch {
                    field: "producer_claim_identity"
                })
            ));

            let changed_target = PendingRender::begin_static(
                &scene,
                cx,
                test_settings(settings.seed ^ 1),
                execution,
            )
            .expect("admit stale job target");
            assert!(matches!(
                changed_target.restore_checkpoint(binding, &bytes, MAX_CHECKPOINT_BYTES, cx),
                Err(RenderCheckpointError::JobMismatch { field: "seed" })
            ));

            let mut emitted = false;
            assert!(matches!(
                pending.write_checkpoint::<Infallible>(binding, receipt.byte_len() - 1, cx, |_| {
                    emitted = true;
                    Ok(())
                },),
                Err(RenderCheckpointWriteError::Checkpoint(
                    RenderCheckpointError::ByteLimitExceeded { .. }
                ))
            ));
            assert!(!emitted, "one-short admission emitted partial bytes");
        });
    }

    #[test]
    fn g3_adaptive_codec_round_trips_all_aovs_exactly() {
        let gate = CancelGate::new_clock_free();
        with_gate_cx(&gate, |cx| {
            let scene = scene();
            let settings = test_settings(0x6368_6563_6b70_0202);
            let execution = execution();
            let policy = AdaptiveSamplingConfig::try_new(2, 1, 0.0, 0.0, 0.0)
                .expect("valid adaptive checkpoint policy");
            let pending = PendingAdaptiveRender::begin_static(
                &scene,
                cx,
                settings,
                policy,
                execution.clone(),
            )
            .expect("admit adaptive checkpoint job");
            let binding = binding_with_job(pending.checkpoint_job_identity());
            let mut bytes = Vec::new();
            let receipt = pending
                .write_checkpoint::<Infallible>(binding, MAX_CHECKPOINT_BYTES, cx, |chunk| {
                    bytes.extend_from_slice(chunk);
                    Ok(())
                })
                .expect("encode adaptive checkpoint");
            assert_eq!(receipt.kind(), RenderCheckpointKind::Adaptive);
            let target = PendingAdaptiveRender::begin_static(
                &scene,
                cx,
                settings,
                policy,
                execution.clone(),
            )
            .expect("admit adaptive restore target");
            let (restored, restored_receipt) = target
                .restore_checkpoint(binding, &bytes, MAX_CHECKPOINT_BYTES, cx)
                .expect("restore adaptive checkpoint");
            assert_eq!(restored_receipt, receipt);
            let restored = restored.resume(cx).expect("finish restored adaptive");
            let reference =
                PendingAdaptiveRender::begin_static(&scene, cx, settings, policy, execution)
                    .expect("admit adaptive reference")
                    .resume(cx)
                    .expect("finish adaptive reference");
            assert_adaptive_bits_eq(&restored.film, &reference.film, "adaptive checkpoint");
        });
    }

    #[test]
    fn g0_well_sealed_semantically_invalid_state_is_refused() {
        let gate = CancelGate::new_clock_free();
        with_gate_cx(&gate, |cx| {
            let scene = scene();
            let settings = test_settings(0x6368_6563_6b70_0404);
            let execution = execution();
            let uniform = PendingRender::begin_static(&scene, cx, settings, execution.clone())
                .expect("admit semantic-corruption uniform job");
            let uniform_binding = binding_with_job(uniform.checkpoint_job_identity());
            let mut uniform_bytes = Vec::new();
            uniform
                .write_checkpoint::<Infallible>(
                    uniform_binding,
                    MAX_CHECKPOINT_BYTES,
                    cx,
                    |chunk| {
                        uniform_bytes.extend_from_slice(chunk);
                        Ok(())
                    },
                )
                .expect("encode semantic-corruption uniform fixture");
            let uniform_body = uniform_bytes[..uniform_bytes.len() - SEAL_BYTES as usize].to_vec();
            let uniform_state = uniform_state_offset(&uniform, uniform_binding, &uniform_body);

            let mut malformed_uniform = Vec::new();
            let mut reordered = uniform_body.clone();
            reordered[uniform_state..uniform_state + 8].copy_from_slice(&1_u64.to_le_bytes());
            malformed_uniform.push(("reordered tile", reseal_body(reordered)));
            let mut bad_row = uniform_body.clone();
            bad_row[uniform_state + 24..uniform_state + 28]
                .copy_from_slice(&(uniform.layout.tile_height + 1).to_le_bytes());
            malformed_uniform.push(("row prefix outside tile", reseal_body(bad_row)));
            let mut nonzero_uncommitted = uniform_body.clone();
            nonzero_uncommitted[uniform_state + TILE_RECORD_BYTES as usize
                ..uniform_state + TILE_RECORD_BYTES as usize + 8]
                .copy_from_slice(&1.0_f64.to_bits().to_le_bytes());
            malformed_uniform.push((
                "nonzero uncommitted pixel",
                reseal_body(nonzero_uncommitted),
            ));
            let mut trailing = uniform_body.clone();
            trailing.push(0xA5);
            malformed_uniform.push(("trailing field", reseal_body(trailing)));
            let mut missing = uniform_body.clone();
            missing.pop().expect("uniform body is nonempty");
            malformed_uniform.push(("missing final pixel byte", reseal_body(missing)));

            for (label, bytes) in malformed_uniform {
                let target = PendingRender::begin_static(&scene, cx, settings, execution.clone())
                    .expect("re-admit malformed uniform target");
                assert!(
                    target
                        .restore_checkpoint(uniform_binding, &bytes, MAX_CHECKPOINT_BYTES, cx,)
                        .is_err(),
                    "well-sealed {label} payload was accepted"
                );
            }

            let policy = AdaptiveSamplingConfig::try_new(2, 1, 0.0, 0.0, 0.0)
                .expect("admit semantic-corruption adaptive policy");
            let adaptive = PendingAdaptiveRender::begin_static(
                &scene,
                cx,
                settings,
                policy,
                execution.clone(),
            )
            .expect("admit semantic-corruption adaptive job")
            .advance_to_safe_point(cx, NonZeroU32::MIN)
            .expect("commit one row per adaptive tile")
            .into_pending();
            let adaptive_binding = binding_with_job(adaptive.checkpoint_job_identity());
            let mut adaptive_bytes = Vec::new();
            adaptive
                .write_checkpoint::<Infallible>(
                    adaptive_binding,
                    MAX_CHECKPOINT_BYTES,
                    cx,
                    |chunk| {
                        adaptive_bytes.extend_from_slice(chunk);
                        Ok(())
                    },
                )
                .expect("encode semantic-corruption adaptive fixture");
            let adaptive_body =
                adaptive_bytes[..adaptive_bytes.len() - SEAL_BYTES as usize].to_vec();
            let adaptive_state = adaptive_state_offset(&adaptive, adaptive_binding, &adaptive_body);
            let pixel = adaptive_state + TILE_RECORD_BYTES as usize;

            let mut contradictory_sum = adaptive_body.clone();
            let old_sum = f64::from_bits(u64::from_le_bytes(
                contradictory_sum[pixel..pixel + 8]
                    .try_into()
                    .expect("exact sum width"),
            ));
            contradictory_sum[pixel..pixel + 8]
                .copy_from_slice(&(old_sum + 1.0).to_bits().to_le_bytes());
            let mut negative_zero_m2 = adaptive_body.clone();
            negative_zero_m2[pixel + 48..pixel + 56]
                .copy_from_slice(&(-0.0_f64).to_bits().to_le_bytes());
            let mut missing_aov = adaptive_body.clone();
            missing_aov.pop().expect("adaptive body is nonempty");

            for (label, bytes) in [
                (
                    "contradictory sum/mean/count",
                    reseal_body(contradictory_sum),
                ),
                ("negative-zero M2", reseal_body(negative_zero_m2)),
                ("missing adaptive AOV byte", reseal_body(missing_aov)),
            ] {
                let target = PendingAdaptiveRender::begin_static(
                    &scene,
                    cx,
                    settings,
                    policy,
                    execution.clone(),
                )
                .expect("re-admit malformed adaptive target");
                assert!(
                    matches!(
                        target.restore_checkpoint(
                            adaptive_binding,
                            &bytes,
                            MAX_CHECKPOINT_BYTES,
                            cx,
                        ),
                        Err(RenderCheckpointError::InvalidPixelState)
                            | Err(RenderCheckpointError::Truncated)
                    ),
                    "well-sealed {label} payload was accepted"
                );
            }
        });
    }

    #[test]
    fn g4_sink_failure_and_final_seal_cancellation_issue_no_receipt() {
        let gate = CancelGate::new_clock_free();
        with_gate_cx(&gate, |cx| {
            let scene = scene();
            let pending = PendingRender::begin_static(
                &scene,
                cx,
                test_settings(0x6368_6563_6b70_0303),
                execution(),
            )
            .expect("admit cancellable checkpoint job");
            let binding = binding_with_job(pending.checkpoint_job_identity());
            let mut reentered_progress = None;
            let sink_error = pending.write_checkpoint(binding, MAX_CHECKPOINT_BYTES, cx, |_| {
                reentered_progress = Some(pending.progress());
                Err::<(), _>("injected checkpoint sink failure")
            });
            assert_eq!(reentered_progress, Some(pending.progress()));
            assert!(matches!(
                sink_error,
                Err(RenderCheckpointWriteError::Sink(
                    "injected checkpoint sink failure"
                ))
            ));

            let mut saw_seal = false;
            let cancelled = pending.write_checkpoint::<Infallible>(
                binding,
                MAX_CHECKPOINT_BYTES,
                cx,
                |chunk| {
                    if chunk.starts_with(SEAL_MAGIC) {
                        saw_seal = true;
                        gate.request();
                    }
                    Ok(())
                },
            );
            assert!(saw_seal, "fixture never reached the final checkpoint seal");
            assert!(matches!(
                cancelled,
                Err(RenderCheckpointWriteError::Checkpoint(
                    RenderCheckpointError::Cancelled
                ))
            ));
        });
    }
}
