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
use core::num::NonZeroU32;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_ledger::{ArtifactInfo, ArtifactWriter, Ledger, LedgerError, PutReceipt};
use fs_render::camera::CutSide;
use fs_render::motion::{ShutterConvention, ShutterDistribution};
use fs_render::tracer::{
    AdaptiveRenderCheckpointYield, AdaptiveRenderOutput, AdaptiveRenderSuspension,
    AdaptiveSamplingConfig, ParkedRenderScope, PendingAdaptiveRender, PendingRender,
    RenderCheckpointBinding, RenderCheckpointError, RenderCheckpointReceipt,
    RenderCheckpointWriteError, RenderCheckpointYield, RenderExecutionConfig, RenderExecutionError,
    RenderExecutionOutput, RenderExecutionReport, RenderProgress, RenderSuspension, Settings,
};

use crate::render_scene_bridge::{EulerCinematicScene, EulerPreparedFrame, EulerSceneError};

/// Canonical identity schema for the Euler-to-renderer checkpoint adapter.
pub const EULER_RENDER_CHECKPOINT_IDENTITY_VERSION: u16 = 1;
/// Domain for one resolved, event-delimited Euler frame segment.
pub const EULER_RENDER_CHECKPOINT_FRAME_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.render-checkpoint-frame.v1";
/// Ledger artifact kind used for fs-render checkpoint codec v1 bytes.
pub const EULER_RENDER_CHECKPOINT_ARTIFACT_KIND: &str = "euler-render-checkpoint-v1";

/// Explicit producer identities attached to one durable checkpoint generation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EulerRenderCheckpointProducer {
    producer_build: ContentHash,
    producer_claim: ContentHash,
}

impl EulerRenderCheckpointProducer {
    /// Validate the build and claim identities of a checkpoint producer.
    pub fn try_new(
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
}

/// Restore-only description of the exact durable generation expected.
///
/// Publication never accepts this type: a fresh job cannot use an expectation
/// to mint a successor. Successor publication is derived only from the private
/// head carried by an [`EulerUniformCheckpointJob`] or
/// [`EulerAdaptiveCheckpointJob`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EulerRenderCheckpointExpectation {
    /// A generation-zero checkpoint.
    Root(EulerRenderCheckpointProducer),
    /// A generation immediately following a typed predecessor.
    Successor {
        /// Producer of the expected generation.
        producer: EulerRenderCheckpointProducer,
        /// Verified predecessor generation.
        predecessor: EulerStoredRenderCheckpoint,
    },
}

impl EulerRenderCheckpointExpectation {
    /// Expect a root generation produced by `producer`.
    #[must_use]
    pub const fn root(producer: EulerRenderCheckpointProducer) -> Self {
        Self::Root(producer)
    }

    /// Expect the generation immediately after `predecessor`.
    #[must_use]
    pub const fn successor(
        producer: EulerRenderCheckpointProducer,
        predecessor: EulerStoredRenderCheckpoint,
    ) -> Self {
        Self::Successor {
            producer,
            predecessor,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FrozenEulerFrameJob {
    source_artifact: ContentHash,
    source_configuration: ContentHash,
    scene: ContentHash,
    frame: ContentHash,
    render_job: ContentHash,
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

/// Sealed fixed-SPP durable job. Its pending state and lineage head cannot be
/// separated or replaced through the safe adapter API.
#[must_use = "advance, checkpoint, or resume the durable render job"]
pub struct EulerUniformCheckpointJob<'scene> {
    pending: PendingRender<'scene>,
    frame_job: FrozenEulerFrameJob,
    head: Option<EulerStoredRenderCheckpoint>,
}

/// Sealed adaptive durable job with the same state/lineage ownership rule.
#[must_use = "advance, checkpoint, or resume the durable adaptive render job"]
pub struct EulerAdaptiveCheckpointJob<'scene> {
    pending: PendingAdaptiveRender<'scene>,
    frame_job: FrozenEulerFrameJob,
    head: Option<EulerStoredRenderCheckpoint>,
}

/// Successful bounded fixed-SPP advance that retains the sealed durable job.
#[must_use = "inspect the attempt or recover the durable render job"]
pub struct EulerUniformCheckpointYield<'scene> {
    yielded: RenderCheckpointYield<'scene>,
    frame_job: FrozenEulerFrameJob,
    head: Option<EulerStoredRenderCheckpoint>,
}

/// Failed/cancelled fixed-SPP attempt that retains the sealed durable job.
#[must_use = "inspect the refusal or recover the durable render job"]
pub struct EulerUniformCheckpointSuspension<'scene> {
    suspended: RenderSuspension<'scene>,
    frame_job: FrozenEulerFrameJob,
    head: Option<EulerStoredRenderCheckpoint>,
}

/// Successful bounded adaptive advance retaining all statistical AOV state.
#[must_use = "inspect the attempt or recover the durable adaptive render job"]
pub struct EulerAdaptiveCheckpointYield<'scene> {
    yielded: AdaptiveRenderCheckpointYield<'scene>,
    frame_job: FrozenEulerFrameJob,
    head: Option<EulerStoredRenderCheckpoint>,
}

/// Failed/cancelled adaptive attempt retaining the sealed durable job.
#[must_use = "inspect the refusal or recover the durable adaptive render job"]
pub struct EulerAdaptiveCheckpointSuspension<'scene> {
    suspended: AdaptiveRenderSuspension<'scene>,
    frame_job: FrozenEulerFrameJob,
    head: Option<EulerStoredRenderCheckpoint>,
}

impl fmt::Debug for EulerUniformCheckpointJob<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("EulerUniformCheckpointJob")
            .field("progress", &self.progress())
            .field("head", &self.head)
            .finish_non_exhaustive()
    }
}

impl fmt::Debug for EulerAdaptiveCheckpointJob<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("EulerAdaptiveCheckpointJob")
            .field("progress", &self.progress())
            .field("head", &self.head)
            .finish_non_exhaustive()
    }
}

impl fmt::Debug for EulerUniformCheckpointYield<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("EulerUniformCheckpointYield")
            .field("yielded", &self.yielded)
            .field("head", &self.head)
            .finish_non_exhaustive()
    }
}

impl fmt::Debug for EulerAdaptiveCheckpointYield<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("EulerAdaptiveCheckpointYield")
            .field("yielded", &self.yielded)
            .field("head", &self.head)
            .finish_non_exhaustive()
    }
}

impl fmt::Debug for EulerUniformCheckpointSuspension<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("EulerUniformCheckpointSuspension")
            .field("suspended", &self.suspended)
            .field("head", &self.head)
            .finish_non_exhaustive()
    }
}

impl fmt::Debug for EulerAdaptiveCheckpointSuspension<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("EulerAdaptiveCheckpointSuspension")
            .field("suspended", &self.suspended)
            .field("head", &self.head)
            .finish_non_exhaustive()
    }
}

/// Fail-closed adapter diagnostics.
#[derive(Debug)]
pub enum EulerRenderCheckpointError {
    /// A producer identity or lineage field was a placeholder.
    InvalidProvenance(&'static str),
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
    /// Euler scene/frame or render admission failed.
    Scene(EulerSceneError),
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
    /// A typed predecessor receipt did not describe the artifact/state loaded
    /// from the ledger, so successor authority could not be reconstructed.
    PredecessorReceiptMismatch {
        /// First receipt component that disagreed.
        field: &'static str,
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

impl From<EulerSceneError> for EulerRenderCheckpointError {
    fn from(error: EulerSceneError) -> Self {
        Self::Scene(error)
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

/// Atomically admit a fixed-SPP pending render and freeze the exact Euler frame
/// identity used by durable checkpointing.
pub fn begin_uniform_checkpoint_job<'scene>(
    scene: &'scene EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    settings: Settings,
    execution: RenderExecutionConfig,
    cx: &Cx<'_>,
) -> Result<EulerUniformCheckpointJob<'scene>, EulerRenderCheckpointError> {
    let pending = scene.begin_segment_render(prepared, segment_index, settings, execution, cx)?;
    let frame_job = freeze_frame_job(
        scene,
        prepared,
        segment_index,
        pending.checkpoint_job_identity(),
    )?;
    Ok(EulerUniformCheckpointJob {
        pending,
        frame_job,
        head: None,
    })
}

/// Atomically admit an adaptive pending render and freeze the exact Euler frame
/// identity used by durable checkpointing.
pub fn begin_adaptive_checkpoint_job<'scene>(
    scene: &'scene EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    settings: Settings,
    adaptive: AdaptiveSamplingConfig,
    execution: RenderExecutionConfig,
    cx: &Cx<'_>,
) -> Result<EulerAdaptiveCheckpointJob<'scene>, EulerRenderCheckpointError> {
    let pending = scene.begin_segment_adaptive_render(
        prepared,
        segment_index,
        settings,
        adaptive,
        execution,
        cx,
    )?;
    let frame_job = freeze_frame_job(
        scene,
        prepared,
        segment_index,
        pending.checkpoint_job_identity(),
    )?;
    Ok(EulerAdaptiveCheckpointJob {
        pending,
        frame_job,
        head: None,
    })
}

impl<'scene> EulerUniformCheckpointJob<'scene> {
    /// Current opaque row/tile progress.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.pending.progress()
    }

    /// Latest durable generation carried by this exact state object.
    #[must_use]
    pub const fn head(&self) -> Option<EulerStoredRenderCheckpoint> {
        self.head
    }

    /// Advance to a bounded row-atomic safe point without exposing pixels.
    pub fn advance_to_safe_point(
        self,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<EulerUniformCheckpointYield<'scene>, EulerUniformCheckpointSuspension<'scene>> {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        match pending.advance_to_safe_point(cx, rows_per_incomplete_tile) {
            Ok(yielded) => Ok(EulerUniformCheckpointYield {
                yielded,
                frame_job,
                head,
            }),
            Err(suspended) => Err(EulerUniformCheckpointSuspension {
                suspended,
                frame_job,
                head,
            }),
        }
    }

    /// Parked-crew form of [`Self::advance_to_safe_point`].
    pub fn advance_to_safe_point_on_parked(
        self,
        parked: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<EulerUniformCheckpointYield<'scene>, EulerUniformCheckpointSuspension<'scene>> {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        match pending.advance_to_safe_point_on_parked(parked, cx, rows_per_incomplete_tile) {
            Ok(yielded) => Ok(EulerUniformCheckpointYield {
                yielded,
                frame_job,
                head,
            }),
            Err(suspended) => Err(EulerUniformCheckpointSuspension {
                suspended,
                frame_job,
                head,
            }),
        }
    }

    /// Complete and publish the film on a one-shot worker lane.
    pub fn resume(
        self,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, EulerUniformCheckpointSuspension<'scene>> {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        pending
            .resume(cx)
            .map_err(|suspended| EulerUniformCheckpointSuspension {
                suspended,
                frame_job,
                head,
            })
    }

    /// Complete and publish the film on an already parked worker crew.
    pub fn resume_on_parked(
        self,
        parked: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, EulerUniformCheckpointSuspension<'scene>> {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        pending
            .resume_on_parked(parked, cx)
            .map_err(|suspended| EulerUniformCheckpointSuspension {
                suspended,
                frame_job,
                head,
            })
    }
}

impl<'scene> EulerAdaptiveCheckpointJob<'scene> {
    /// Current opaque row/tile progress.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.pending.progress()
    }

    /// Latest durable generation carried by this exact state object.
    #[must_use]
    pub const fn head(&self) -> Option<EulerStoredRenderCheckpoint> {
        self.head
    }

    /// Advance adaptive state to a bounded row-atomic safe point.
    pub fn advance_to_safe_point(
        self,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<EulerAdaptiveCheckpointYield<'scene>, EulerAdaptiveCheckpointSuspension<'scene>>
    {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        match pending.advance_to_safe_point(cx, rows_per_incomplete_tile) {
            Ok(yielded) => Ok(EulerAdaptiveCheckpointYield {
                yielded,
                frame_job,
                head,
            }),
            Err(suspended) => Err(EulerAdaptiveCheckpointSuspension {
                suspended,
                frame_job,
                head,
            }),
        }
    }

    /// Parked-crew form of [`Self::advance_to_safe_point`].
    pub fn advance_to_safe_point_on_parked(
        self,
        parked: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<EulerAdaptiveCheckpointYield<'scene>, EulerAdaptiveCheckpointSuspension<'scene>>
    {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        match pending.advance_to_safe_point_on_parked(parked, cx, rows_per_incomplete_tile) {
            Ok(yielded) => Ok(EulerAdaptiveCheckpointYield {
                yielded,
                frame_job,
                head,
            }),
            Err(suspended) => Err(EulerAdaptiveCheckpointSuspension {
                suspended,
                frame_job,
                head,
            }),
        }
    }

    /// Complete and publish the adaptive film on a one-shot worker lane.
    pub fn resume(
        self,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, EulerAdaptiveCheckpointSuspension<'scene>> {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        pending
            .resume(cx)
            .map_err(|suspended| EulerAdaptiveCheckpointSuspension {
                suspended,
                frame_job,
                head,
            })
    }

    /// Complete and publish the adaptive film on an already parked crew.
    pub fn resume_on_parked(
        self,
        parked: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, EulerAdaptiveCheckpointSuspension<'scene>> {
        let Self {
            pending,
            frame_job,
            head,
        } = self;
        pending.resume_on_parked(parked, cx).map_err(|suspended| {
            EulerAdaptiveCheckpointSuspension {
                suspended,
                frame_job,
                head,
            }
        })
    }
}

impl<'scene> EulerUniformCheckpointYield<'scene> {
    /// Report for the bounded worker attempt.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        self.yielded.attempt_report()
    }

    /// Current row-atomic progress.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.yielded.progress()
    }

    /// Recover the sealed durable job.
    #[must_use]
    pub fn into_job(self) -> EulerUniformCheckpointJob<'scene> {
        EulerUniformCheckpointJob {
            pending: self.yielded.into_pending(),
            frame_job: self.frame_job,
            head: self.head,
        }
    }
}

impl<'scene> EulerAdaptiveCheckpointYield<'scene> {
    /// Report for the bounded adaptive worker attempt.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        self.yielded.attempt_report()
    }

    /// Current row-atomic progress.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.yielded.progress()
    }

    /// Recover the sealed durable adaptive job.
    #[must_use]
    pub fn into_job(self) -> EulerAdaptiveCheckpointJob<'scene> {
        EulerAdaptiveCheckpointJob {
            pending: self.yielded.into_pending(),
            frame_job: self.frame_job,
            head: self.head,
        }
    }
}

impl<'scene> EulerUniformCheckpointSuspension<'scene> {
    /// Structured renderer refusal.
    #[must_use]
    pub const fn cause(&self) -> &RenderExecutionError {
        self.suspended.cause()
    }

    /// Report for the refused worker attempt.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        self.suspended.attempt_report()
    }

    /// Current committed progress retained after refusal.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.suspended.progress()
    }

    /// Recover the sealed durable job for retry.
    #[must_use]
    pub fn into_job(self) -> EulerUniformCheckpointJob<'scene> {
        EulerUniformCheckpointJob {
            pending: self.suspended.into_pending(),
            frame_job: self.frame_job,
            head: self.head,
        }
    }
}

impl<'scene> EulerAdaptiveCheckpointSuspension<'scene> {
    /// Structured adaptive renderer refusal.
    #[must_use]
    pub const fn cause(&self) -> &RenderExecutionError {
        self.suspended.cause()
    }

    /// Report for the refused adaptive worker attempt.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        self.suspended.attempt_report()
    }

    /// Current committed progress retained after refusal.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.suspended.progress()
    }

    /// Recover the sealed durable adaptive job for retry.
    #[must_use]
    pub fn into_job(self) -> EulerAdaptiveCheckpointJob<'scene> {
        EulerAdaptiveCheckpointJob {
            pending: self.suspended.into_pending(),
            frame_job: self.frame_job,
            head: self.head,
        }
    }
}

/// Stream the next fixed-SPP generation into an atomic ledger artifact write.
/// The generation and predecessor are derived only from `job.head`.
pub fn store_uniform_render_checkpoint(
    ledger: &mut Ledger,
    job: &mut EulerUniformCheckpointJob<'_>,
    producer: EulerRenderCheckpointProducer,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<EulerStoredRenderCheckpoint, EulerRenderCheckpointError> {
    let binding = next_binding(job.frame_job, job.head, producer)?;
    let mut writer = ledger.artifact_writer(EULER_RENDER_CHECKPOINT_ARTIFACT_KIND)?;
    let mut streamed = 0_u64;
    let checkpoint = job
        .pending
        .write_checkpoint(binding, max_bytes, cx, |chunk| {
            write_counted(&mut writer, &mut streamed, chunk)
        })
        .map_err(map_write_error)?;
    require_streamed_length(streamed, checkpoint)?;
    let artifact = writer.finish(None)?;
    debug_assert_eq!(artifact.len, checkpoint.byte_len());
    let stored = EulerStoredRenderCheckpoint {
        artifact,
        checkpoint,
    };
    job.head = Some(stored);
    Ok(stored)
}

/// Stream the next adaptive generation into an atomic ledger artifact write.
pub fn store_adaptive_render_checkpoint(
    ledger: &mut Ledger,
    job: &mut EulerAdaptiveCheckpointJob<'_>,
    producer: EulerRenderCheckpointProducer,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<EulerStoredRenderCheckpoint, EulerRenderCheckpointError> {
    let binding = next_binding(job.frame_job, job.head, producer)?;
    let mut writer = ledger.artifact_writer(EULER_RENDER_CHECKPOINT_ARTIFACT_KIND)?;
    let mut streamed = 0_u64;
    let checkpoint = job
        .pending
        .write_checkpoint(binding, max_bytes, cx, |chunk| {
            write_counted(&mut writer, &mut streamed, chunk)
        })
        .map_err(map_write_error)?;
    require_streamed_length(streamed, checkpoint)?;
    let artifact = writer.finish(None)?;
    debug_assert_eq!(artifact.len, checkpoint.byte_len());
    let stored = EulerStoredRenderCheckpoint {
        artifact,
        checkpoint,
    };
    job.head = Some(stored);
    Ok(stored)
}

/// Atomically admit and restore one fixed-SPP durable job.
///
/// `max_bytes` bounds the root artifact, or the aggregate predecessor and
/// successor bytes when `expectation` names a successor.
pub fn restore_uniform_render_checkpoint<'scene>(
    ledger: &Ledger,
    artifact: ContentHash,
    scene: &'scene EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    settings: Settings,
    execution: RenderExecutionConfig,
    expectation: EulerRenderCheckpointExpectation,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<EulerUniformCheckpointJob<'scene>, EulerRenderCheckpointError> {
    let mut job =
        begin_uniform_checkpoint_job(scene, prepared, segment_index, settings, execution, cx)?;
    let (pending, receipt, artifact_receipt) = match expectation {
        EulerRenderCheckpointExpectation::Root(producer) => {
            let binding = binding_for(job.frame_job, producer, 0, None)?;
            let (bytes, artifact_receipt) = load_checkpoint_bytes(ledger, artifact, max_bytes)?;
            let (pending, receipt) = job
                .pending
                .restore_checkpoint(binding, &bytes, max_bytes, cx)?;
            reconcile_loaded_receipt(bytes.len(), receipt)?;
            (pending, receipt, artifact_receipt)
        }
        EulerRenderCheckpointExpectation::Successor {
            producer,
            predecessor,
        } => {
            let binding = successor_binding(job.frame_job, producer, predecessor)?;
            let ((predecessor_bytes, predecessor_artifact), (bytes, artifact_receipt)) =
                load_checkpoint_pair(ledger, predecessor, artifact, max_bytes)?;
            let predecessor_binding = predecessor.checkpoint().binding();
            let (pending, decoded_predecessor) = job.pending.restore_checkpoint(
                predecessor_binding,
                &predecessor_bytes,
                max_bytes,
                cx,
            )?;
            reconcile_loaded_receipt(predecessor_bytes.len(), decoded_predecessor)?;
            require_predecessor_receipts(predecessor, predecessor_artifact, decoded_predecessor)?;
            let (pending, receipt) =
                pending.restore_successor_checkpoint(binding, &bytes, max_bytes, cx)?;
            reconcile_loaded_receipt(bytes.len(), receipt)?;
            (pending, receipt, artifact_receipt)
        }
    };
    let stored = reconcile_receipts(artifact_receipt, receipt)?;
    job.pending = pending;
    job.head = Some(stored);
    Ok(job)
}

/// Atomically admit and restore one adaptive durable job.
///
/// `max_bytes` bounds the root artifact, or the aggregate predecessor and
/// successor bytes when `expectation` names a successor.
pub fn restore_adaptive_render_checkpoint<'scene>(
    ledger: &Ledger,
    artifact: ContentHash,
    scene: &'scene EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    settings: Settings,
    adaptive: AdaptiveSamplingConfig,
    execution: RenderExecutionConfig,
    expectation: EulerRenderCheckpointExpectation,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<EulerAdaptiveCheckpointJob<'scene>, EulerRenderCheckpointError> {
    let mut job = begin_adaptive_checkpoint_job(
        scene,
        prepared,
        segment_index,
        settings,
        adaptive,
        execution,
        cx,
    )?;
    let (pending, receipt, artifact_receipt) = match expectation {
        EulerRenderCheckpointExpectation::Root(producer) => {
            let binding = binding_for(job.frame_job, producer, 0, None)?;
            let (bytes, artifact_receipt) = load_checkpoint_bytes(ledger, artifact, max_bytes)?;
            let (pending, receipt) = job
                .pending
                .restore_checkpoint(binding, &bytes, max_bytes, cx)?;
            reconcile_loaded_receipt(bytes.len(), receipt)?;
            (pending, receipt, artifact_receipt)
        }
        EulerRenderCheckpointExpectation::Successor {
            producer,
            predecessor,
        } => {
            let binding = successor_binding(job.frame_job, producer, predecessor)?;
            let ((predecessor_bytes, predecessor_artifact), (bytes, artifact_receipt)) =
                load_checkpoint_pair(ledger, predecessor, artifact, max_bytes)?;
            let predecessor_binding = predecessor.checkpoint().binding();
            let (pending, decoded_predecessor) = job.pending.restore_checkpoint(
                predecessor_binding,
                &predecessor_bytes,
                max_bytes,
                cx,
            )?;
            reconcile_loaded_receipt(predecessor_bytes.len(), decoded_predecessor)?;
            require_predecessor_receipts(predecessor, predecessor_artifact, decoded_predecessor)?;
            let (pending, receipt) =
                pending.restore_successor_checkpoint(binding, &bytes, max_bytes, cx)?;
            reconcile_loaded_receipt(bytes.len(), receipt)?;
            (pending, receipt, artifact_receipt)
        }
    };
    let stored = reconcile_receipts(artifact_receipt, receipt)?;
    job.pending = pending;
    job.head = Some(stored);
    Ok(job)
}

fn freeze_frame_job(
    scene: &EulerCinematicScene<'_>,
    prepared: &EulerPreparedFrame,
    segment_index: usize,
    render_job: ContentHash,
) -> Result<FrozenEulerFrameJob, EulerRenderCheckpointError> {
    Ok(FrozenEulerFrameJob {
        source_artifact: scene.source_trajectory_identity(),
        source_configuration: scene.source_configuration_identity(),
        scene: scene.scene_identity(),
        frame: euler_render_checkpoint_frame_identity(prepared, segment_index)?,
        render_job,
    })
}

fn next_binding(
    frame_job: FrozenEulerFrameJob,
    head: Option<EulerStoredRenderCheckpoint>,
    producer: EulerRenderCheckpointProducer,
) -> Result<RenderCheckpointBinding, EulerRenderCheckpointError> {
    match head {
        None => binding_for(frame_job, producer, 0, None),
        Some(predecessor) => successor_binding(frame_job, producer, predecessor),
    }
}

fn successor_binding(
    frame_job: FrozenEulerFrameJob,
    producer: EulerRenderCheckpointProducer,
    predecessor: EulerStoredRenderCheckpoint,
) -> Result<RenderCheckpointBinding, EulerRenderCheckpointError> {
    let prior = predecessor.checkpoint();
    let generation = prior.binding().generation().checked_add(1).ok_or(
        EulerRenderCheckpointError::InvalidProvenance("checkpoint generation overflow"),
    )?;
    let successor = binding_for(frame_job, producer, generation, Some(prior.content_hash()))?;
    let predecessor = prior.binding();
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
    Ok(successor)
}

fn binding_for(
    frame_job: FrozenEulerFrameJob,
    producer: EulerRenderCheckpointProducer,
    generation: u64,
    predecessor: Option<ContentHash>,
) -> Result<RenderCheckpointBinding, EulerRenderCheckpointError> {
    Ok(RenderCheckpointBinding::try_new(
        frame_job.source_artifact,
        frame_job.source_configuration,
        frame_job.scene,
        frame_job.frame,
        frame_job.render_job,
        producer.producer_build,
        producer.producer_claim,
        generation,
        predecessor,
    )?)
}

fn write_counted(
    writer: &mut ArtifactWriter<'_>,
    streamed: &mut u64,
    chunk: &[u8],
) -> Result<(), LedgerError> {
    let chunk_len = u64::try_from(chunk.len()).map_err(|_| LedgerError::Invalid {
        field: "render_checkpoint_chunk_len".to_string(),
        problem: "checkpoint chunk length does not fit u64".to_string(),
    })?;
    let next = streamed
        .checked_add(chunk_len)
        .ok_or_else(|| LedgerError::Invalid {
            field: "render_checkpoint_stream_len".to_string(),
            problem: "checkpoint stream length overflowed u64".to_string(),
        })?;
    writer.write(chunk)?;
    *streamed = next;
    Ok(())
}

fn require_streamed_length(
    streamed: u64,
    checkpoint: RenderCheckpointReceipt,
) -> Result<(), EulerRenderCheckpointError> {
    if streamed != checkpoint.byte_len() {
        return Err(EulerRenderCheckpointError::ArtifactLengthMismatch {
            ledger: streamed,
            renderer: checkpoint.byte_len(),
        });
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

type LoadedCheckpoint = (Vec<u8>, PutReceipt);

/// Materialize a typed predecessor and its proposed successor under one
/// aggregate byte ceiling. Successor verification needs both immutable states
/// resident at once, so treating `max_bytes` as two independent per-artifact
/// limits would silently admit twice the declared checkpoint payload.
fn load_checkpoint_pair(
    ledger: &Ledger,
    predecessor: EulerStoredRenderCheckpoint,
    successor_artifact: ContentHash,
    max_bytes: u64,
) -> Result<(LoadedCheckpoint, LoadedCheckpoint), EulerRenderCheckpointError> {
    let predecessor_hash = predecessor.artifact().hash;
    let predecessor_info = checkpoint_artifact_info(ledger, predecessor_hash)?;
    let successor_info = checkpoint_artifact_info(ledger, successor_artifact)?;
    let required = predecessor_info
        .len
        .checked_add(successor_info.len)
        .ok_or(RenderCheckpointError::LengthOverflow)?;
    if required > max_bytes {
        return Err(RenderCheckpointError::ByteLimitExceeded {
            required,
            limit: max_bytes,
        }
        .into());
    }
    let predecessor_bytes = ledger
        .get_artifact_bounded(&predecessor_hash, predecessor_info.len)?
        .ok_or(EulerRenderCheckpointError::MissingArtifact(
            predecessor_hash,
        ))?;
    let successor_bytes = ledger
        .get_artifact_bounded(&successor_artifact, successor_info.len)?
        .ok_or(EulerRenderCheckpointError::MissingArtifact(
            successor_artifact,
        ))?;
    Ok((
        (
            predecessor_bytes,
            receipt_from_info(predecessor_info, predecessor_hash),
        ),
        (
            successor_bytes,
            receipt_from_info(successor_info, successor_artifact),
        ),
    ))
}

fn receipt_from_info(info: ArtifactInfo, hash: ContentHash) -> PutReceipt {
    PutReceipt {
        hash,
        len: info.len,
        deduped: true,
        chunked: info.chunk_count != 0,
    }
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

fn require_predecessor_receipts(
    expected: EulerStoredRenderCheckpoint,
    artifact: PutReceipt,
    checkpoint: RenderCheckpointReceipt,
) -> Result<(), EulerRenderCheckpointError> {
    let expected_artifact = expected.artifact();
    if artifact.hash != expected_artifact.hash {
        return Err(EulerRenderCheckpointError::PredecessorReceiptMismatch {
            field: "artifact_hash",
        });
    }
    if artifact.len != expected_artifact.len {
        return Err(EulerRenderCheckpointError::PredecessorReceiptMismatch {
            field: "artifact_length",
        });
    }
    if checkpoint != expected.checkpoint() {
        return Err(EulerRenderCheckpointError::PredecessorReceiptMismatch {
            field: "renderer_receipt",
        });
    }
    Ok(())
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
