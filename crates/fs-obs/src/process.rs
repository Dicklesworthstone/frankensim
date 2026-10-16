//! Bounded process-stream capture policy (i94v.7.3.3).
//!
//! Process runners own pipes, durable stores, and cancellation tokens. This
//! I/O-free module owns the deterministic admission decision: enqueue,
//! durably spill, explicitly drop, backpressure, or fail proof closure.
//!
//! Proof-critical frames are never converted into lossy gaps. Diagnostic and
//! telemetry frames may be omitted only after a quantified [`ProcessGap`] is
//! retained. A full gap ledger therefore applies backpressure too.

use crate::{
    EvidenceCompleteness, EvidenceIntegrity, PromotionEffect, ScopedReceiptOutcome, Severity,
};
use core::fmt;
use std::collections::VecDeque;

/// Current process-capture policy semantics.
pub const PROCESS_CAPTURE_POLICY_VERSION: u32 = 1;

fn valid_bounded_text(value: &str, max_len: usize) -> bool {
    !value.is_empty() && value.len() <= max_len && value.chars().all(|ch| !ch.is_control())
}

/// Validated identity of one child process or process-like execution scope.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ProcessId(String);

impl ProcessId {
    /// Validate and adopt a process identity.
    ///
    /// # Errors
    /// [`ProcessIdError`] when the identity is empty, longer than 256 bytes,
    /// or contains control characters.
    pub fn new(raw: impl Into<String>) -> Result<Self, ProcessIdError> {
        let raw = raw.into();
        if valid_bounded_text(&raw, 256) {
            Ok(Self(raw))
        } else {
            Err(ProcessIdError { raw })
        }
    }

    /// Validated raw identity.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Refusal for a malformed [`ProcessId`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProcessIdError {
    /// Rejected text.
    pub raw: String,
}

impl fmt::Display for ProcessIdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "process identity must be 1..=256 bytes without control characters; got {:?}",
            self.raw
        )
    }
}

impl core::error::Error for ProcessIdError {}

/// Child stream from which a frame originated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProcessStream {
    /// Standard output, captured as data rather than printed by a core crate.
    Stdout,
    /// Standard error, captured separately from standard output.
    Stderr,
}

impl ProcessStream {
    /// Stable policy/wire name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Stdout => "stdout",
            Self::Stderr => "stderr",
        }
    }
}

/// Loss contract attached to one process frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProcessEventClass {
    /// Proof, promotion, lifecycle, identity, or control-critical evidence.
    Critical,
    /// Diagnostics that may be coalesced or truncated with an explicit gap.
    Diagnostic,
    /// Performance/operational telemetry governed by declared sampling.
    Telemetry,
}

impl ProcessEventClass {
    /// Stable policy/wire name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Critical => "critical",
            Self::Diagnostic => "diagnostic",
            Self::Telemetry => "telemetry",
        }
    }
}

/// Admission token for an already-committed full-payload artifact.
///
/// The type proves use of the committed-pointer constructor. Independently
/// authenticating the external store and hash remains the sink's job.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct DurableArtifactPointer {
    artifact: String,
    content_hash: String,
}

impl DurableArtifactPointer {
    /// Adopt an already-committed artifact identity and retained content hash.
    ///
    /// # Errors
    /// [`ArtifactPointerError`] when either component is malformed.
    pub fn committed(
        artifact: impl Into<String>,
        content_hash: impl Into<String>,
    ) -> Result<Self, ArtifactPointerError> {
        let artifact = artifact.into();
        let content_hash = content_hash.into();
        if !valid_bounded_text(&artifact, 512) {
            return Err(ArtifactPointerError::InvalidArtifact { artifact });
        }
        if !(16..=128).contains(&content_hash.len())
            || !content_hash.bytes().all(|byte| byte.is_ascii_hexdigit())
        {
            return Err(ArtifactPointerError::InvalidContentHash { content_hash });
        }
        Ok(Self {
            artifact,
            content_hash,
        })
    }

    /// Retained artifact identity.
    #[must_use]
    pub fn artifact(&self) -> &str {
        &self.artifact
    }

    /// Retained hexadecimal content hash.
    #[must_use]
    pub fn content_hash(&self) -> &str {
        &self.content_hash
    }
}

/// Refusal for an invalid committed-artifact pointer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArtifactPointerError {
    /// Artifact identity is empty, oversized, or control-bearing.
    InvalidArtifact {
        /// Rejected identity.
        artifact: String,
    },
    /// Content hash is not 16..=128 hexadecimal characters.
    InvalidContentHash {
        /// Rejected hash.
        content_hash: String,
    },
}

impl fmt::Display for ArtifactPointerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidArtifact { artifact } => write!(
                f,
                "artifact identity must be 1..=512 bytes without control characters; got {artifact:?}"
            ),
            Self::InvalidContentHash { content_hash } => write!(
                f,
                "content hash must be 16..=128 hexadecimal characters; got {content_hash:?}"
            ),
        }
    }
}

impl core::error::Error for ArtifactPointerError {}

/// One framed stdout/stderr chunk. Payload bytes are opaque and may be
/// non-UTF-8.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProcessFrame {
    process: ProcessId,
    stream: ProcessStream,
    ordinal: u64,
    class: ProcessEventClass,
    severity: Severity,
    payload: Vec<u8>,
    original_bytes: usize,
    artifact: Option<DurableArtifactPointer>,
}

impl ProcessFrame {
    /// Construct a frame. Ordinals are one-based and monotone per stream.
    ///
    /// # Errors
    /// [`ProcessFrameError::ZeroOrdinal`] when `ordinal == 0`.
    pub fn new(
        process: ProcessId,
        stream: ProcessStream,
        ordinal: u64,
        class: ProcessEventClass,
        severity: Severity,
        payload: impl Into<Vec<u8>>,
    ) -> Result<Self, ProcessFrameError> {
        if ordinal == 0 {
            return Err(ProcessFrameError::ZeroOrdinal);
        }
        let payload = payload.into();
        let original_bytes = payload.len();
        Ok(Self {
            process,
            stream,
            ordinal,
            class,
            severity,
            payload,
            original_bytes,
            artifact: None,
        })
    }

    /// Bind an already-committed full-payload artifact to this frame.
    #[must_use]
    pub fn with_committed_artifact(mut self, artifact: DurableArtifactPointer) -> Self {
        self.artifact = Some(artifact);
        self
    }

    /// Process identity.
    #[must_use]
    pub fn process(&self) -> &ProcessId {
        &self.process
    }

    /// Source stream.
    #[must_use]
    pub const fn stream(&self) -> ProcessStream {
        self.stream
    }

    /// One-based source ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> u64 {
        self.ordinal
    }

    /// Loss class.
    #[must_use]
    pub const fn class(&self) -> ProcessEventClass {
        self.class
    }

    /// Event severity.
    #[must_use]
    pub const fn severity(&self) -> Severity {
        self.severity
    }

    /// Inline bytes retained in the queue.
    #[must_use]
    pub fn payload(&self) -> &[u8] {
        &self.payload
    }

    /// Original payload byte count before bounded inline projection.
    #[must_use]
    pub const fn original_bytes(&self) -> usize {
        self.original_bytes
    }

    /// Durable full-payload pointer, when committed before admission.
    #[must_use]
    pub fn artifact(&self) -> Option<&DurableArtifactPointer> {
        self.artifact.as_ref()
    }

    /// Consume this captured frame into the canonical typed fs-obs event
    /// stream. The frame's severity becomes the event-envelope severity so a
    /// runtime owner cannot accidentally stamp a conflicting level.
    #[must_use]
    pub fn into_event(self, emitter: &mut crate::Emitter, wall_ns: Option<u64>) -> crate::Event {
        let Self {
            process,
            stream,
            ordinal,
            class,
            severity,
            payload,
            original_bytes,
            artifact,
        } = self;
        let original_bytes = u64::try_from(original_bytes)
            .expect("supported Rust targets represent usize within u64");
        emitter.emit(
            severity,
            crate::EventKind::ProcessFrame {
                process,
                stream,
                ordinal,
                class,
                payload,
                original_bytes,
                artifact,
            },
            wall_ns,
        )
    }
}

/// Refusal for an invalid process frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProcessFrameError {
    /// Source ordinals are one-based.
    ZeroOrdinal,
}

impl fmt::Display for ProcessFrameError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("process frame ordinal must be greater than zero")
    }
}

impl core::error::Error for ProcessFrameError {}

/// Bounded process-capture policy. Both payload and gap metadata have explicit
/// resident ceilings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProcessCapturePolicy {
    policy_version: u32,
    queue_capacity: usize,
    gap_capacity: usize,
    max_inline_bytes: usize,
    telemetry_sample_every: u64,
}

impl ProcessCapturePolicy {
    /// Validate one policy. A telemetry stride of `N` retains ordinals
    /// `1, 1 + N, ...`; every omission still produces a gap.
    ///
    /// # Errors
    /// [`CapturePolicyError`] when any value is zero.
    pub fn new(
        policy_version: u32,
        queue_capacity: usize,
        gap_capacity: usize,
        max_inline_bytes: usize,
        telemetry_sample_every: u64,
    ) -> Result<Self, CapturePolicyError> {
        let zero = if policy_version == 0 {
            Some("policy_version")
        } else if queue_capacity == 0 {
            Some("queue_capacity")
        } else if gap_capacity == 0 {
            Some("gap_capacity")
        } else if max_inline_bytes == 0 {
            Some("max_inline_bytes")
        } else if telemetry_sample_every == 0 {
            Some("telemetry_sample_every")
        } else {
            None
        };
        if let Some(field) = zero {
            return Err(CapturePolicyError { field });
        }
        Ok(Self {
            policy_version,
            queue_capacity,
            gap_capacity,
            max_inline_bytes,
            telemetry_sample_every,
        })
    }

    /// Version recorded in every gap.
    #[must_use]
    pub const fn policy_version(&self) -> u32 {
        self.policy_version
    }

    /// Maximum resident frames.
    #[must_use]
    pub const fn queue_capacity(&self) -> usize {
        self.queue_capacity
    }

    /// Maximum resident non-coalesced gap records.
    #[must_use]
    pub const fn gap_capacity(&self) -> usize {
        self.gap_capacity
    }

    /// Maximum inline bytes per queued frame.
    #[must_use]
    pub const fn max_inline_bytes(&self) -> usize {
        self.max_inline_bytes
    }

    /// Deterministic telemetry sampling stride.
    #[must_use]
    pub const fn telemetry_sample_every(&self) -> u64 {
        self.telemetry_sample_every
    }
}

impl Default for ProcessCapturePolicy {
    fn default() -> Self {
        Self {
            policy_version: PROCESS_CAPTURE_POLICY_VERSION,
            queue_capacity: 1_024,
            gap_capacity: 256,
            max_inline_bytes: 64 * 1_024,
            telemetry_sample_every: 1,
        }
    }
}

/// Refusal for an unbounded or versionless policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CapturePolicyError {
    /// Zero-valued field.
    pub field: &'static str,
}

impl fmt::Display for CapturePolicyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "process capture policy field {} must be non-zero",
            self.field
        )
    }
}

impl core::error::Error for CapturePolicyError {}

/// Cancellation observation supplied at admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CaptureCancellation {
    /// Ordinary bounded backpressure is permitted.
    Running,
    /// Critical backpressure cannot wait indefinitely and must fail closed.
    Requested,
}

/// External sink failure observed before frame admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SinkFailure {
    /// Consumer closed its pipe.
    BrokenPipe,
    /// Retained storage is full.
    DiskExhausted,
    /// Declared retention budget is exhausted.
    RetentionBudgetExhausted,
    /// Sink is otherwise unavailable.
    Unavailable,
}

impl SinkFailure {
    const fn loss_reason(self) -> LossReason {
        match self {
            Self::BrokenPipe => LossReason::BrokenPipe,
            Self::DiskExhausted => LossReason::DiskExhausted,
            Self::RetentionBudgetExhausted => LossReason::RetentionBudgetExhausted,
            Self::Unavailable => LossReason::SinkUnavailable,
        }
    }
}

/// Sink state supplied with one decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SinkAvailability {
    /// Sink can accept queued work.
    Ready,
    /// Sink failure is already known.
    Failed(SinkFailure),
}

/// Stable reason for an explicit omission or inline projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LossReason {
    /// Deterministic sampling omitted telemetry.
    PolicySampling,
    /// Resident frame capacity was exhausted.
    QueueCapacity,
    /// Inline limit truncated data without retained detail.
    InlineLimit,
    /// Omitted inline detail exists in a committed artifact.
    DurableSpill,
    /// Pipe consumer closed unexpectedly.
    BrokenPipe,
    /// Retained storage exhausted disk allocation.
    DiskExhausted,
    /// Declared retention budget was exhausted.
    RetentionBudgetExhausted,
    /// Sink was otherwise unavailable.
    SinkUnavailable,
}

impl LossReason {
    /// Stable policy/wire name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::PolicySampling => "policy_sampling",
            Self::QueueCapacity => "queue_capacity",
            Self::InlineLimit => "inline_limit",
            Self::DurableSpill => "durable_spill",
            Self::BrokenPipe => "broken_pipe",
            Self::DiskExhausted => "disk_exhausted",
            Self::RetentionBudgetExhausted => "retention_budget_exhausted",
            Self::SinkUnavailable => "sink_unavailable",
        }
    }
}

/// Quantified omission range. A pointer means omitted inline bytes remain in
/// an already-committed artifact; no pointer means explicit loss.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProcessGap {
    /// Process identity.
    pub process: ProcessId,
    /// Source stream.
    pub stream: ProcessStream,
    /// Loss class.
    pub class: ProcessEventClass,
    /// First affected ordinal.
    pub first_ordinal: u64,
    /// Last affected ordinal.
    pub last_ordinal: u64,
    /// Frames presented in the range.
    pub original_count: u128,
    /// Frames still represented in the queue.
    pub emitted_count: u128,
    /// Whole frames omitted from the queue.
    pub dropped_count: u128,
    /// Bytes omitted from inline payloads.
    pub omitted_inline_bytes: u128,
    /// Omission reason.
    pub reason: LossReason,
    /// Exact policy version.
    pub policy_version: u32,
    /// Durable retained detail, when committed before recording.
    pub artifact: Option<DurableArtifactPointer>,
}

impl ProcessGap {
    /// Bytes with no durable retained pointer.
    #[must_use]
    pub fn lost_bytes(&self) -> u128 {
        if self.artifact.is_some() {
            0
        } else {
            self.omitted_inline_bytes
        }
    }

    /// Consume this quantified omission into the canonical typed fs-obs event
    /// stream. Gaps are warnings: proof-critical loss takes the separate
    /// fail-closed [`CaptureDecision::CriticalFailure`] path and is never
    /// laundered into an ordinary gap.
    #[must_use]
    pub fn into_event(self, emitter: &mut crate::Emitter, wall_ns: Option<u64>) -> crate::Event {
        let Self {
            process,
            stream,
            class,
            first_ordinal,
            last_ordinal,
            original_count,
            emitted_count,
            dropped_count,
            omitted_inline_bytes,
            reason,
            policy_version,
            artifact,
        } = self;
        emitter.emit(
            Severity::Warn,
            crate::EventKind::ProcessGap {
                process,
                stream,
                class,
                first_ordinal,
                last_ordinal,
                original_count,
                emitted_count,
                dropped_count,
                omitted_inline_bytes,
                reason,
                policy_version,
                artifact,
            },
            wall_ns,
        )
    }

    fn can_coalesce(&self, next: &Self) -> bool {
        self.process == next.process
            && self.stream == next.stream
            && self.class == next.class
            && self.reason == next.reason
            && self.policy_version == next.policy_version
            && self.artifact == next.artifact
            && self.last_ordinal.checked_add(1) == Some(next.first_ordinal)
    }

    fn coalesce(&mut self, next: &Self) {
        self.last_ordinal = next.last_ordinal;
        self.original_count += next.original_count;
        self.emitted_count += next.emitted_count;
        self.dropped_count += next.dropped_count;
        self.omitted_inline_bytes += next.omitted_inline_bytes;
    }
}

/// Resident and cumulative capture counters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CaptureMetrics {
    /// Frames consumed by enqueue/drop decisions.
    pub original_frames: u128,
    /// Frames admitted to the resident queue.
    pub enqueued_frames: u128,
    /// Whole diagnostic/telemetry frames dropped.
    pub dropped_frames: u128,
    /// Inline bytes omitted, whether lost or spilled.
    pub omitted_inline_bytes: u128,
    /// Frames whose omitted detail is durably retained.
    pub durable_spill_frames: u128,
    /// Decisions requiring drain/spill/retry.
    pub backpressure_decisions: u128,
    /// Unretained proof-critical frames.
    pub critical_failures: u128,
    /// Peak resident frame count.
    pub max_queue_depth: usize,
    /// Peak resident non-coalesced gap count.
    pub max_gap_depth: usize,
}

/// Why a caller must drain, spill, or retry without consuming a frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackpressureReason {
    /// Frame queue is full.
    QueueCapacity,
    /// Gap ledger is full and cannot coalesce the next gap.
    GapLedgerCapacity,
    /// Oversized critical frame lacks a committed full-payload artifact.
    DurableSpillRequired,
}

/// Sink/cancellation failure that makes proof closure non-eligible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriticalFailureReason {
    /// Critical event could not reach its sink.
    Sink(SinkFailure),
    /// Cancellation prevented waiting for capture to recover.
    CancelledWhileBackpressured(BackpressureReason),
}

/// Fail-closed effect for an unretained critical frame or missing receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReceiptImpact {
    /// Required evidence is incomplete.
    pub completeness: EvidenceCompleteness,
    /// Missing bytes cannot be authenticated.
    pub integrity: EvidenceIntegrity,
    /// Promotion must be refused/demoted.
    pub promotion: PromotionEffect,
}

impl ReceiptImpact {
    /// Standard non-eligible process-capture impact.
    #[must_use]
    pub const fn non_eligible() -> Self {
        Self {
            completeness: EvidenceCompleteness::Incomplete,
            integrity: EvidenceIntegrity::Unchecked,
            promotion: PromotionEffect::Demoted,
        }
    }

    /// Apply only the orthogonal closure dimensions to a scoped outcome.
    pub fn apply_to(self, outcome: &mut ScopedReceiptOutcome) {
        outcome.completeness = self.completeness;
        outcome.integrity = self.integrity;
        outcome.promotion = self.promotion;
    }
}

/// Explicit critical capture failure and receipt effect.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriticalCaptureFailure {
    /// Failure mechanism.
    pub reason: CriticalFailureReason,
    /// Mandatory fail-closed effect.
    pub impact: ReceiptImpact,
}

/// Refusal for a non-monotone source ordinal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CaptureOrderError {
    /// Source stream.
    pub stream: ProcessStream,
    /// Last consumed ordinal.
    pub previous: u64,
    /// Rejected ordinal.
    pub received: u64,
}

/// Result of offering one frame. Non-consumption variants return the untouched
/// frame so the caller retains ownership.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CaptureDecision {
    /// Frame entered the bounded queue.
    Enqueued {
        /// Source ordinal.
        ordinal: u64,
        /// Inline bytes omitted.
        omitted_inline_bytes: usize,
        /// Whether committed storage retains omitted detail.
        durably_spilled: bool,
    },
    /// Noncritical frame was consumed after a gap was retained.
    Dropped {
        /// Source ordinal.
        ordinal: u64,
        /// Recorded loss reason.
        reason: LossReason,
    },
    /// Caller must drain/spill/retry.
    Backpressure {
        /// Untouched frame.
        frame: ProcessFrame,
        /// Blocking resource.
        reason: BackpressureReason,
    },
    /// Critical evidence was not retained.
    CriticalFailure {
        /// Untouched frame.
        frame: ProcessFrame,
        /// Failure and evidence effect.
        failure: CriticalCaptureFailure,
    },
    /// Ordinal would rewind/fork a stream.
    Rejected {
        /// Untouched frame.
        frame: ProcessFrame,
        /// Ordering refusal.
        error: CaptureOrderError,
    },
}

/// Deterministic, I/O-free bounded capture state.
#[derive(Debug, Clone)]
pub struct ProcessCapture {
    policy: ProcessCapturePolicy,
    frames: VecDeque<ProcessFrame>,
    gaps: VecDeque<ProcessGap>,
    metrics: CaptureMetrics,
    last_stdout_ordinal: Option<u64>,
    last_stderr_ordinal: Option<u64>,
}

impl ProcessCapture {
    /// Start empty under a validated policy.
    #[must_use]
    pub fn new(policy: ProcessCapturePolicy) -> Self {
        Self {
            policy,
            frames: VecDeque::new(),
            gaps: VecDeque::new(),
            metrics: CaptureMetrics::default(),
            last_stdout_ordinal: None,
            last_stderr_ordinal: None,
        }
    }

    /// Policy in force.
    #[must_use]
    pub const fn policy(&self) -> &ProcessCapturePolicy {
        &self.policy
    }

    /// Current and cumulative metrics.
    #[must_use]
    pub const fn metrics(&self) -> CaptureMetrics {
        self.metrics
    }

    /// Resident frame count.
    #[must_use]
    pub fn frame_len(&self) -> usize {
        self.frames.len()
    }

    /// Resident non-coalesced gap count.
    #[must_use]
    pub fn gap_len(&self) -> usize {
        self.gaps.len()
    }

    /// Drain the oldest queued frame.
    pub fn pop_frame(&mut self) -> Option<ProcessFrame> {
        self.frames.pop_front()
    }

    /// Drain the oldest gap.
    pub fn pop_gap(&mut self) -> Option<ProcessGap> {
        self.gaps.pop_front()
    }

    /// Offer one frame under current cancellation and sink state. No I/O or
    /// lifecycle/control lock is held by this method.
    pub fn offer(
        &mut self,
        frame: ProcessFrame,
        cancellation: CaptureCancellation,
        sink: SinkAvailability,
    ) -> CaptureDecision {
        if let Some(previous) = self.last_ordinal(frame.stream)
            && frame.ordinal <= previous
        {
            let error = CaptureOrderError {
                stream: frame.stream,
                previous,
                received: frame.ordinal,
            };
            return CaptureDecision::Rejected { frame, error };
        }

        if frame.class == ProcessEventClass::Telemetry
            && !(frame.ordinal - 1).is_multiple_of(self.policy.telemetry_sample_every)
        {
            return self.drop_with_gap(frame, LossReason::PolicySampling);
        }

        if let SinkAvailability::Failed(failure) = sink {
            if frame.class == ProcessEventClass::Critical {
                return self.critical_failure(frame, CriticalFailureReason::Sink(failure));
            }
            return self.drop_with_gap(frame, failure.loss_reason());
        }

        if self.frames.len() == self.policy.queue_capacity {
            if frame.class == ProcessEventClass::Critical {
                return self.critical_backpressure(
                    frame,
                    cancellation,
                    BackpressureReason::QueueCapacity,
                );
            }
            return self.drop_with_gap(frame, LossReason::QueueCapacity);
        }

        if frame.payload.len() > self.policy.max_inline_bytes {
            return self.offer_oversized(frame, cancellation);
        }

        self.enqueue(frame, 0, false)
    }

    fn offer_oversized(
        &mut self,
        mut frame: ProcessFrame,
        cancellation: CaptureCancellation,
    ) -> CaptureDecision {
        let original_len = frame.payload.len();
        let has_artifact = frame.artifact.is_some();
        if frame.class == ProcessEventClass::Critical && !has_artifact {
            return self.critical_backpressure(
                frame,
                cancellation,
                BackpressureReason::DurableSpillRequired,
            );
        }

        let omitted_inline_bytes = original_len - self.policy.max_inline_bytes;
        let reason = if has_artifact {
            LossReason::DurableSpill
        } else {
            LossReason::InlineLimit
        };
        let gap = self.gap_for(&frame, 1, 0, omitted_inline_bytes, reason);
        if !self.record_gap(gap) {
            self.metrics.backpressure_decisions += 1;
            return CaptureDecision::Backpressure {
                frame,
                reason: BackpressureReason::GapLedgerCapacity,
            };
        }
        frame.payload.truncate(self.policy.max_inline_bytes);
        self.metrics.omitted_inline_bytes += omitted_inline_bytes as u128;
        if has_artifact {
            self.metrics.durable_spill_frames += 1;
        }
        self.enqueue(frame, omitted_inline_bytes, has_artifact)
    }

    fn critical_backpressure(
        &mut self,
        frame: ProcessFrame,
        cancellation: CaptureCancellation,
        reason: BackpressureReason,
    ) -> CaptureDecision {
        match cancellation {
            CaptureCancellation::Running => {
                self.metrics.backpressure_decisions += 1;
                CaptureDecision::Backpressure { frame, reason }
            }
            CaptureCancellation::Requested => self.critical_failure(
                frame,
                CriticalFailureReason::CancelledWhileBackpressured(reason),
            ),
        }
    }

    fn critical_failure(
        &mut self,
        frame: ProcessFrame,
        reason: CriticalFailureReason,
    ) -> CaptureDecision {
        self.metrics.critical_failures += 1;
        CaptureDecision::CriticalFailure {
            frame,
            failure: CriticalCaptureFailure {
                reason,
                impact: ReceiptImpact::non_eligible(),
            },
        }
    }

    fn drop_with_gap(&mut self, frame: ProcessFrame, reason: LossReason) -> CaptureDecision {
        let gap = self.gap_for(&frame, 0, 1, frame.payload.len(), reason);
        if !self.record_gap(gap) {
            self.metrics.backpressure_decisions += 1;
            return CaptureDecision::Backpressure {
                frame,
                reason: BackpressureReason::GapLedgerCapacity,
            };
        }
        self.consume_ordinal(frame.stream, frame.ordinal);
        self.metrics.original_frames += 1;
        self.metrics.dropped_frames += 1;
        self.metrics.omitted_inline_bytes += frame.payload.len() as u128;
        CaptureDecision::Dropped {
            ordinal: frame.ordinal,
            reason,
        }
    }

    fn enqueue(
        &mut self,
        frame: ProcessFrame,
        omitted_inline_bytes: usize,
        durably_spilled: bool,
    ) -> CaptureDecision {
        let ordinal = frame.ordinal;
        self.consume_ordinal(frame.stream, ordinal);
        self.frames.push_back(frame);
        self.metrics.original_frames += 1;
        self.metrics.enqueued_frames += 1;
        self.metrics.max_queue_depth = self.metrics.max_queue_depth.max(self.frames.len());
        CaptureDecision::Enqueued {
            ordinal,
            omitted_inline_bytes,
            durably_spilled,
        }
    }

    fn gap_for(
        &self,
        frame: &ProcessFrame,
        emitted_count: u128,
        dropped_count: u128,
        omitted_inline_bytes: usize,
        reason: LossReason,
    ) -> ProcessGap {
        ProcessGap {
            process: frame.process.clone(),
            stream: frame.stream,
            class: frame.class,
            first_ordinal: frame.ordinal,
            last_ordinal: frame.ordinal,
            original_count: 1,
            emitted_count,
            dropped_count,
            omitted_inline_bytes: omitted_inline_bytes as u128,
            reason,
            policy_version: self.policy.policy_version,
            artifact: frame.artifact.clone(),
        }
    }

    fn record_gap(&mut self, gap: ProcessGap) -> bool {
        if let Some(last) = self.gaps.back_mut()
            && last.can_coalesce(&gap)
        {
            last.coalesce(&gap);
            return true;
        }
        if self.gaps.len() == self.policy.gap_capacity {
            return false;
        }
        self.gaps.push_back(gap);
        self.metrics.max_gap_depth = self.metrics.max_gap_depth.max(self.gaps.len());
        true
    }

    const fn last_ordinal(&self, stream: ProcessStream) -> Option<u64> {
        match stream {
            ProcessStream::Stdout => self.last_stdout_ordinal,
            ProcessStream::Stderr => self.last_stderr_ordinal,
        }
    }

    fn consume_ordinal(&mut self, stream: ProcessStream, ordinal: u64) {
        match stream {
            ProcessStream::Stdout => self.last_stdout_ordinal = Some(ordinal),
            ProcessStream::Stderr => self.last_stderr_ordinal = Some(ordinal),
        }
    }
}

/// Validated final operation-receipt identity produced by a child.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FinalReceiptId(String);

impl FinalReceiptId {
    /// Validate and adopt a final receipt identity.
    ///
    /// # Errors
    /// [`FinalReceiptIdError`] when empty, oversized, or control-bearing.
    pub fn new(raw: impl Into<String>) -> Result<Self, FinalReceiptIdError> {
        let raw = raw.into();
        if valid_bounded_text(&raw, 512) {
            Ok(Self(raw))
        } else {
            Err(FinalReceiptIdError { raw })
        }
    }

    /// Validated receipt identity.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Refusal for an invalid final receipt identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinalReceiptIdError {
    /// Rejected text.
    pub raw: String,
}

impl fmt::Display for FinalReceiptIdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "final receipt identity must be 1..=512 bytes without control characters; got {:?}",
            self.raw
        )
    }
}

impl core::error::Error for FinalReceiptIdError {}

/// Observed child termination. This never substitutes for a final receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProcessTermination {
    /// OS exit code was observed.
    Exited(i32),
    /// Process ended from an OS signal.
    Signaled(u32),
    /// Host disappeared before terminal process state was observed.
    HostLost,
}

/// Reconciliation of process termination with the child's typed receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProcessClosure {
    /// A final receipt exists; its typed outcome is authoritative.
    Reconciled {
        /// Process identity.
        process: ProcessId,
        /// Observed termination.
        termination: ProcessTermination,
        /// Retained final receipt.
        final_receipt: FinalReceiptId,
    },
    /// No final receipt exists, so even exit zero remains an explicit gap.
    ObservationGap {
        /// Process identity.
        process: ProcessId,
        /// Observed termination.
        termination: ProcessTermination,
        /// Stable diagnostic reason.
        reason: &'static str,
        /// Mandatory fail-closed effect.
        impact: ReceiptImpact,
    },
}

impl ProcessClosure {
    /// Reconcile observed termination with an optional final receipt.
    #[must_use]
    pub fn reconcile(
        process: ProcessId,
        termination: ProcessTermination,
        final_receipt: Option<FinalReceiptId>,
    ) -> Self {
        match final_receipt {
            Some(final_receipt) => Self::Reconciled {
                process,
                termination,
                final_receipt,
            },
            None => Self::ObservationGap {
                process,
                termination,
                reason: "process terminated without a final operation receipt",
                impact: ReceiptImpact::non_eligible(),
            },
        }
    }
}
