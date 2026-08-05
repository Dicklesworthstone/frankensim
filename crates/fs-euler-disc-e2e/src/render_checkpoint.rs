//! Durable Euler cinematic-render checkpoints in the Design Ledger.
//!
//! This L6 adapter owns the identity boundary between an event-delimited
//! Euler frame and fs-render's versioned checkpoint codec. Checkpoint bytes
//! stream directly into an [`fs_ledger::ArtifactWriter`]; an interrupted or
//! failed write is rolled back by that writer and cannot replace an earlier
//! content-addressed checkpoint.
//!
//! Every read is bounded by a caller-supplied byte ceiling. This module makes
//! no blanket claim that a particular memory limit is sufficient for 4K:
//! resolution, adaptive state, and codec overhead all affect the required
//! budget, and fs-render remains the admission authority for a concrete job.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_ledger::{ArtifactInfo, Ledger, LedgerError, PutReceipt};
use fs_render::camera::{AnimatedCamera, CameraError, CutSide};
use fs_render::motion::{ShutterConvention, ShutterDistribution};
use fs_render::tracer::{
    AdaptiveSamplingConfig, FilmTimeMode, PendingAdaptiveRender, PendingRender,
    RenderCheckpointBinding, RenderCheckpointError, RenderCheckpointReceipt,
    RenderCheckpointWriteError, RenderExecutionConfig, Scene, Settings,
    adaptive_checkpoint_job_identity, uniform_checkpoint_job_identity,
};

use crate::render_scene_bridge::{EulerCinematicScene, EulerPreparedFrame};

/// Canonical identity schema for the Euler-to-renderer checkpoint adapter.
pub const EULER_RENDER_CHECKPOINT_IDENTITY_VERSION: u16 = 1;
/// Domain for one resolved, event-delimited Euler frame segment.
pub const EULER_RENDER_CHECKPOINT_FRAME_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.render-checkpoint-frame.v1";
/// Ledger artifact kind used for fs-render checkpoint codec v1 bytes.
pub const EULER_RENDER_CHECKPOINT_ARTIFACT_KIND: &str = "euler-render-checkpoint-v1";

/// Explicit producer and lineage inputs attached to a checkpoint binding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderCheckpointProvenance {
    producer_build: ContentHash,
    producer_claim: ContentHash,
    generation: u64,
    predecessor: Option<ContentHash>,
    predecessor_binding: Option<RenderCheckpointBinding>,
}

impl EulerRenderCheckpointProvenance {
    /// Construct root provenance from non-placeholder producer identities.
    /// Non-root lineage can only be derived from [`Self::try_successor`].
    pub fn try_root(
        producer_build: ContentHash,
        producer_claim: ContentHash,
    ) -> Result<Self, EulerRenderCheckpointError> {
        if is_zero(producer_build) {
            return Err(EulerRenderCheckpointError::InvalidProvenance(
                "producer_build must be nonzero",
            ));
        }
        if is_zero(producer_claim) {
            return Err(EulerRenderCheckpointError::InvalidProvenance(
                "producer_claim must be nonzero",
            ));
        }
        Ok(Self {
            producer_build,
            producer_claim,
            generation: 0,
            predecessor: None,
            predecessor_binding: None,
        })
    }

    /// Derive the next generation from a checkpoint this adapter published.
    ///
    /// This is the preferred non-root constructor: it obtains both the next
    /// generation and predecessor renderer-content identity from a typed,
    /// successfully stored receipt rather than accepting unrelated scalars.
    pub fn try_successor(
        producer_build: ContentHash,
        producer_claim: ContentHash,
        predecessor: EulerStoredRenderCheckpoint,
    ) -> Result<Self, EulerRenderCheckpointError> {
        let prior = predecessor.checkpoint;
        let generation = prior.binding().generation().checked_add(1).ok_or(
            EulerRenderCheckpointError::InvalidProvenance("checkpoint generation overflow"),
        )?;
        if is_zero(producer_build) || is_zero(producer_claim) {
            return Err(EulerRenderCheckpointError::InvalidProvenance(
                "successor producer build and claim must be nonzero",
            ));
        }
        Ok(Self {
            producer_build,
            producer_claim,
            generation,
            predecessor: Some(prior.content_hash()),
            predecessor_binding: Some(prior.binding()),
        })
    }

    /// Exact producer build identity.
    #[must_use]
    pub const fn producer_build(self) -> ContentHash {
        self.producer_build
    }

    /// Scientific or product claim identity under which the render runs.
    #[must_use]
    pub const fn producer_claim(self) -> ContentHash {
        self.producer_claim
    }

    /// Checkpoint generation supplied by orchestration; zero denotes the root.
    #[must_use]
    pub const fn generation(self) -> u64 {
        self.generation
    }

    /// Exact renderer-content identity of the preceding checkpoint, absent at the root.
    #[must_use]
    pub const fn predecessor(self) -> Option<ContentHash> {
        self.predecessor
    }
}

/// Canonical renderer binding produced from an admitted Euler scene and frame.
///
/// The inner binding is private so persistence through this adapter cannot
/// accidentally omit the source artifact, configuration, scene, or canonical
/// frame/job identities.
#[derive(Clone, Copy)]
pub struct EulerRenderCheckpointBinding<'scene> {
    renderer: RenderCheckpointBinding,
    scene: &'scene Scene,
    camera: &'scene AnimatedCamera,
}

impl fmt::Debug for EulerRenderCheckpointBinding<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("EulerRenderCheckpointBinding")
            .field("renderer", &self.renderer)
            .finish_non_exhaustive()
    }
}

impl PartialEq for EulerRenderCheckpointBinding<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.renderer == other.renderer
    }
}

impl Eq for EulerRenderCheckpointBinding<'_> {}

impl EulerRenderCheckpointBinding<'_> {
    /// Frozen fs-render binding for lower-level orchestration and inspection.
    #[must_use]
    pub const fn renderer(self) -> RenderCheckpointBinding {
        self.renderer
    }
}

/// Combined receipt for renderer serialization and immutable ledger storage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerStoredRenderCheckpoint {
    artifact: PutReceipt,
    checkpoint: RenderCheckpointReceipt,
}

impl EulerStoredRenderCheckpoint {
    /// Content-addressed Design Ledger receipt.
    #[must_use]
    pub const fn artifact(self) -> PutReceipt {
        self.artifact
    }

    /// Fs-render codec and progress receipt.
    #[must_use]
    pub const fn checkpoint(self) -> RenderCheckpointReceipt {
        self.checkpoint
    }
}

/// Fail-closed adapter diagnostics.
#[derive(Debug)]
pub enum EulerRenderCheckpointError {
    /// A producer identity or lineage field was a placeholder.
    InvalidProvenance(&'static str),
    /// Prepared frame belongs to a different beauty scene.
    PreparedFrameMismatch,
    /// Pending render was admitted from different scene/camera/job inputs.
    PendingJobMismatch(&'static str),
    /// Segment index was outside the prepared frame.
    InvalidPreparedSegment {
        /// Requested segment.
        index: usize,
        /// Available segment count.
        segment_count: usize,
    },
    /// A platform-sized count could not enter the canonical u64 identity.
    IdentityRange(&'static str),
    /// Fs-render refused the binding, codec, budget, or restored state.
    Renderer(RenderCheckpointError),
    /// Camera/shutter admission needed to construct exact frame-job identity failed.
    Camera(CameraError),
    /// Design Ledger storage or bounded retrieval failed.
    Ledger(LedgerError),
    /// The addressed content identity was absent from the ledger.
    MissingArtifact(ContentHash),
    /// Addressed bytes existed under a different ledger artifact contract.
    ArtifactKindMismatch {
        /// Required typed artifact kind.
        expected: &'static str,
        /// Stored artifact kind.
        actual: String,
    },
    /// Renderer and ledger disagreed on the serialized byte count.
    ArtifactLengthMismatch {
        /// Byte count returned by the ledger.
        ledger: u64,
        /// Byte count returned by fs-render's codec.
        renderer: u64,
    },
}

impl fmt::Display for EulerRenderCheckpointError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for EulerRenderCheckpointError {}

impl From<RenderCheckpointError> for EulerRenderCheckpointError {
    fn from(error: RenderCheckpointError) -> Self {
        Self::Renderer(error)
    }
}

impl From<LedgerError> for EulerRenderCheckpointError {
    fn from(error: LedgerError) -> Self {
        Self::Ledger(error)
    }
}

impl From<CameraError> for EulerRenderCheckpointError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}

/// Canonical identity of one resolved prepared-frame segment.
pub fn euler_render_checkpoint_frame_identity(
    prepared: &EulerPreparedFrame,
    segment_index: usize,
) -> Result<ContentHash, EulerRenderCheckpointError> {
    let segment = prepared.segments().get(segment_index).ok_or(
        EulerRenderCheckpointError::InvalidPreparedSegment {
            index: segment_index,
            segment_count: prepared.segments().len(),
        },
    )?;
    let segment_index = u64::try_from(segment_index)
        .map_err(|_| EulerRenderCheckpointError::IdentityRange("segment_index"))?;
    let shutter = segment.shutter();
    let mut hasher = DomainHasher::new(EULER_RENDER_CHECKPOINT_FRAME_IDENTITY_DOMAIN);
    hasher.update(&EULER_RENDER_CHECKPOINT_IDENTITY_VERSION.to_le_bytes());
    hasher.update(&[cut_side_tag(prepared.cut_side())]);
    hasher.update(&segment_index.to_le_bytes());
    hasher.update(&shutter.open_s().to_bits().to_le_bytes());
    hasher.update(&shutter.close_s().to_bits().to_le_bytes());
    hasher.update(&[shutter_convention_tag(shutter.convention())]);
    hash_shutter_distribution(&mut hasher, shutter.distribution());
    hasher.update(&segment.duration_weight().to_bits().to_le_bytes());
    Ok(hasher.finalize())
}

/// Build the complete binding for a fixed-spp pending render.
pub fn try_uniform_render_checkpoint_binding<'scene>(
    scene: &'scene EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    settings: &Settings,
    execution: &RenderExecutionConfig,
    pending: &PendingRender<'_>,
    provenance: EulerRenderCheckpointProvenance,
    cx: &Cx<'_>,
) -> Result<EulerRenderCheckpointBinding<'scene>, EulerRenderCheckpointError> {
    try_binding(
        scene,
        prepared,
        segment_index,
        settings,
        execution,
        None,
        pending.checkpoint_job_identity(),
        pending.checkpoint_uses_cinematic_sources(scene.scene(), scene.camera()),
        provenance,
        cx,
    )
}

/// Build the complete binding for an adaptive pending render.
pub fn try_adaptive_render_checkpoint_binding<'scene>(
    scene: &'scene EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    settings: &Settings,
    execution: &RenderExecutionConfig,
    adaptive: AdaptiveSamplingConfig,
    pending: &PendingAdaptiveRender<'_>,
    provenance: EulerRenderCheckpointProvenance,
    cx: &Cx<'_>,
) -> Result<EulerRenderCheckpointBinding<'scene>, EulerRenderCheckpointError> {
    try_binding(
        scene,
        prepared,
        segment_index,
        settings,
        execution,
        Some(adaptive),
        pending.checkpoint_job_identity(),
        pending.checkpoint_uses_cinematic_sources(scene.scene(), scene.camera()),
        provenance,
        cx,
    )
}

/// Stream one fixed-spp checkpoint into an atomic ledger artifact write.
pub fn store_uniform_render_checkpoint(
    ledger: &mut Ledger,
    pending: &PendingRender<'_>,
    binding: EulerRenderCheckpointBinding<'_>,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<EulerStoredRenderCheckpoint, EulerRenderCheckpointError> {
    require_uniform_pending_binding(pending, binding)?;
    let mut writer = ledger.artifact_writer(EULER_RENDER_CHECKPOINT_ARTIFACT_KIND)?;
    let checkpoint = pending
        .write_checkpoint(binding.renderer, max_bytes, cx, |chunk| writer.write(chunk))
        .map_err(map_write_error)?;
    let artifact = writer.finish(None)?;
    require_checkpoint_artifact_kind(ledger, artifact.hash)?;
    reconcile_receipts(artifact, checkpoint)
}

/// Stream one adaptive checkpoint into an atomic ledger artifact write.
pub fn store_adaptive_render_checkpoint(
    ledger: &mut Ledger,
    pending: &PendingAdaptiveRender<'_>,
    binding: EulerRenderCheckpointBinding<'_>,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<EulerStoredRenderCheckpoint, EulerRenderCheckpointError> {
    require_adaptive_pending_binding(pending, binding)?;
    let mut writer = ledger.artifact_writer(EULER_RENDER_CHECKPOINT_ARTIFACT_KIND)?;
    let checkpoint = pending
        .write_checkpoint(binding.renderer, max_bytes, cx, |chunk| writer.write(chunk))
        .map_err(map_write_error)?;
    let artifact = writer.finish(None)?;
    require_checkpoint_artifact_kind(ledger, artifact.hash)?;
    reconcile_receipts(artifact, checkpoint)
}

/// Load and restore one fixed-spp pending render under an explicit byte limit.
pub fn restore_uniform_render_checkpoint<'assets>(
    ledger: &Ledger,
    artifact: ContentHash,
    pending: PendingRender<'assets>,
    binding: EulerRenderCheckpointBinding<'_>,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<(PendingRender<'assets>, EulerStoredRenderCheckpoint), EulerRenderCheckpointError> {
    require_uniform_pending_binding(&pending, binding)?;
    let (bytes, artifact_receipt) = load_checkpoint_bytes(ledger, artifact, max_bytes)?;
    let (pending, receipt) = pending.restore_checkpoint(binding.renderer, &bytes, max_bytes, cx)?;
    reconcile_loaded_receipt(bytes.len(), receipt)?;
    let stored = reconcile_receipts(artifact_receipt, receipt)?;
    Ok((pending, stored))
}

/// Load and restore one adaptive pending render under an explicit byte limit.
pub fn restore_adaptive_render_checkpoint<'assets>(
    ledger: &Ledger,
    artifact: ContentHash,
    pending: PendingAdaptiveRender<'assets>,
    binding: EulerRenderCheckpointBinding<'_>,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<(PendingAdaptiveRender<'assets>, EulerStoredRenderCheckpoint), EulerRenderCheckpointError>
{
    require_adaptive_pending_binding(&pending, binding)?;
    let (bytes, artifact_receipt) = load_checkpoint_bytes(ledger, artifact, max_bytes)?;
    let (pending, receipt) = pending.restore_checkpoint(binding.renderer, &bytes, max_bytes, cx)?;
    reconcile_loaded_receipt(bytes.len(), receipt)?;
    let stored = reconcile_receipts(artifact_receipt, receipt)?;
    Ok((pending, stored))
}

fn try_binding<'scene>(
    scene: &'scene EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    settings: &Settings,
    execution: &RenderExecutionConfig,
    adaptive: Option<AdaptiveSamplingConfig>,
    pending_job_identity: ContentHash,
    uses_expected_sources: bool,
    provenance: EulerRenderCheckpointProvenance,
    cx: &Cx<'_>,
) -> Result<EulerRenderCheckpointBinding<'scene>, EulerRenderCheckpointError> {
    if prepared.scene_identity() != scene.scene_identity() {
        return Err(EulerRenderCheckpointError::PreparedFrameMismatch);
    }
    let segment = prepared.segments().get(segment_index).ok_or(
        EulerRenderCheckpointError::InvalidPreparedSegment {
            index: segment_index,
            segment_count: prepared.segments().len(),
        },
    )?;
    let exposure = scene
        .camera()
        .admit_shutter(cx, segment.shutter(), prepared.cut_side())?;
    let requested_mode = FilmTimeMode::Cinematic {
        shutter: segment.shutter(),
        stream_identity: settings.seed,
        shot_id: exposure.shot_id(),
    };
    let render_job = match adaptive {
        None => uniform_checkpoint_job_identity(
            settings,
            requested_mode,
            cx.mode(),
            cx.budget(),
            execution,
        ),
        Some(policy) => adaptive_checkpoint_job_identity(
            settings,
            requested_mode,
            cx.mode(),
            cx.budget(),
            execution,
            policy,
        ),
    };
    if !uses_expected_sources || pending_job_identity != render_job {
        return Err(EulerRenderCheckpointError::PendingJobMismatch(
            "pending render does not match the named scene, frame, settings, mode, budget, or execution policy",
        ));
    }
    let frame = euler_render_checkpoint_frame_identity(prepared, segment_index)?;
    let renderer = RenderCheckpointBinding::try_new(
        scene.source_trajectory_identity(),
        scene.source_configuration_identity(),
        scene.scene_identity(),
        frame,
        render_job,
        provenance.producer_build,
        provenance.producer_claim,
        provenance.generation,
        provenance.predecessor,
    )?;
    require_lineage_continuity(provenance.predecessor_binding, renderer)?;
    Ok(EulerRenderCheckpointBinding {
        renderer,
        scene: scene.scene(),
        camera: scene.camera(),
    })
}

fn require_uniform_pending_binding(
    pending: &PendingRender<'_>,
    binding: EulerRenderCheckpointBinding<'_>,
) -> Result<(), EulerRenderCheckpointError> {
    if !pending.checkpoint_uses_cinematic_sources(binding.scene, binding.camera)
        || pending.checkpoint_job_identity() != binding.renderer.render_job_identity()
    {
        return Err(EulerRenderCheckpointError::PendingJobMismatch(
            "uniform pending render and checkpoint binding were constructed from different inputs",
        ));
    }
    Ok(())
}

fn require_adaptive_pending_binding(
    pending: &PendingAdaptiveRender<'_>,
    binding: EulerRenderCheckpointBinding<'_>,
) -> Result<(), EulerRenderCheckpointError> {
    if !pending.checkpoint_uses_cinematic_sources(binding.scene, binding.camera)
        || pending.checkpoint_job_identity() != binding.renderer.render_job_identity()
    {
        return Err(EulerRenderCheckpointError::PendingJobMismatch(
            "adaptive pending render and checkpoint binding were constructed from different inputs",
        ));
    }
    Ok(())
}

fn require_lineage_continuity(
    predecessor: Option<RenderCheckpointBinding>,
    successor: RenderCheckpointBinding,
) -> Result<(), EulerRenderCheckpointError> {
    let Some(predecessor) = predecessor else {
        return Ok(());
    };
    if predecessor.source_artifact_identity() != successor.source_artifact_identity()
        || predecessor.source_configuration_identity() != successor.source_configuration_identity()
        || predecessor.scene_identity() != successor.scene_identity()
        || predecessor.frame_identity() != successor.frame_identity()
        || predecessor.render_job_identity() != successor.render_job_identity()
    {
        return Err(EulerRenderCheckpointError::InvalidProvenance(
            "successor checkpoint must retain source, configuration, scene, frame, and render-job identity",
        ));
    }
    Ok(())
}

fn load_checkpoint_bytes(
    ledger: &Ledger,
    artifact: ContentHash,
    max_bytes: u64,
) -> Result<(Vec<u8>, PutReceipt), EulerRenderCheckpointError> {
    let info = checkpoint_artifact_info(ledger, artifact)?;
    let bytes = ledger
        .get_artifact_bounded(&artifact, max_bytes)?
        .ok_or(EulerRenderCheckpointError::MissingArtifact(artifact))?;
    let receipt = PutReceipt {
        hash: artifact,
        len: info.len,
        deduped: true,
        chunked: info.chunk_count != 0,
    };
    Ok((bytes, receipt))
}

fn require_checkpoint_artifact_kind(
    ledger: &Ledger,
    artifact: ContentHash,
) -> Result<(), EulerRenderCheckpointError> {
    checkpoint_artifact_info(ledger, artifact).map(|_| ())
}

fn checkpoint_artifact_info(
    ledger: &Ledger,
    artifact: ContentHash,
) -> Result<ArtifactInfo, EulerRenderCheckpointError> {
    let info = ledger
        .artifact_info(&artifact)?
        .ok_or(EulerRenderCheckpointError::MissingArtifact(artifact))?;
    if info.kind != EULER_RENDER_CHECKPOINT_ARTIFACT_KIND {
        return Err(EulerRenderCheckpointError::ArtifactKindMismatch {
            expected: EULER_RENDER_CHECKPOINT_ARTIFACT_KIND,
            actual: info.kind,
        });
    }
    Ok(info)
}

fn reconcile_receipts(
    artifact: PutReceipt,
    checkpoint: RenderCheckpointReceipt,
) -> Result<EulerStoredRenderCheckpoint, EulerRenderCheckpointError> {
    if artifact.len != checkpoint.byte_len() {
        return Err(EulerRenderCheckpointError::ArtifactLengthMismatch {
            ledger: artifact.len,
            renderer: checkpoint.byte_len(),
        });
    }
    Ok(EulerStoredRenderCheckpoint {
        artifact,
        checkpoint,
    })
}

fn reconcile_loaded_receipt(
    byte_len: usize,
    checkpoint: RenderCheckpointReceipt,
) -> Result<(), EulerRenderCheckpointError> {
    let byte_len = u64::try_from(byte_len)
        .map_err(|_| EulerRenderCheckpointError::IdentityRange("checkpoint_bytes"))?;
    if byte_len != checkpoint.byte_len() {
        return Err(EulerRenderCheckpointError::ArtifactLengthMismatch {
            ledger: byte_len,
            renderer: checkpoint.byte_len(),
        });
    }
    Ok(())
}

fn map_write_error(error: RenderCheckpointWriteError<LedgerError>) -> EulerRenderCheckpointError {
    match error {
        RenderCheckpointWriteError::Checkpoint(error) => {
            EulerRenderCheckpointError::Renderer(error)
        }
        RenderCheckpointWriteError::Sink(error) => EulerRenderCheckpointError::Ledger(error),
    }
}

fn is_zero(identity: ContentHash) -> bool {
    identity.as_bytes().iter().all(|byte| *byte == 0)
}

const fn cut_side_tag(side: CutSide) -> u8 {
    match side {
        CutSide::Before => 0,
        CutSide::After => 1,
    }
}

const fn shutter_convention_tag(convention: ShutterConvention) -> u8 {
    match convention {
        ShutterConvention::Centered => 0,
        ShutterConvention::FrontLoaded => 1,
        ShutterConvention::BackLoaded => 2,
    }
}

fn hash_shutter_distribution(hasher: &mut DomainHasher, distribution: ShutterDistribution) {
    match distribution {
        ShutterDistribution::UniformCounterV1 => hasher.update(&[0]),
        ShutterDistribution::StratifiedCounterV1 { strata } => {
            hasher.update(&[1]);
            hasher.update(&strata.to_le_bytes());
        }
    }
}
