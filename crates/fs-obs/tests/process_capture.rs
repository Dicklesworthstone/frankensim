//! Process-capture conformance tests for bounded queues, loss accounting,
//! cancellation, durable closure, and retained failure evidence.

use fs_obs::process::{
    BackpressureReason, CaptureCancellation, CaptureDecision, CriticalFailureReason,
    DurableArtifactPointer, FinalReceiptId, LossReason, ProcessCapture, ProcessCapturePolicy,
    ProcessClosure, ProcessEventClass, ProcessFrame, ProcessId, ProcessStream, ProcessTermination,
    SinkAvailability, SinkFailure,
};
use fs_obs::{
    DomainApplicability, Emitter, EpistemicGrade, EventKind, EvidenceCompleteness,
    EvidenceIntegrity, ExecutionDisposition, OperationalSupport, PredicateOutcome, PromotionEffect,
    ReceiptScope, ScopedReceiptOutcome, Severity, lint_failure_record, validate_line,
};

fn process() -> ProcessId {
    ProcessId::new("worker-7").expect("fixture process identity")
}

fn policy(
    queue_capacity: usize,
    gap_capacity: usize,
    max_inline_bytes: usize,
    telemetry_sample_every: u64,
) -> ProcessCapturePolicy {
    ProcessCapturePolicy::new(
        7,
        queue_capacity,
        gap_capacity,
        max_inline_bytes,
        telemetry_sample_every,
    )
    .expect("bounded fixture policy")
}

fn frame(
    stream: ProcessStream,
    ordinal: u64,
    class: ProcessEventClass,
    payload: &[u8],
) -> ProcessFrame {
    ProcessFrame::new(process(), stream, ordinal, class, Severity::Info, payload)
        .expect("valid fixture frame")
}

#[test]
fn constructors_reject_unbounded_or_ambiguous_state() {
    for args in [
        (0, 1, 1, 1, 1),
        (1, 0, 1, 1, 1),
        (1, 1, 0, 1, 1),
        (1, 1, 1, 0, 1),
        (1, 1, 1, 1, 0),
    ] {
        assert!(ProcessCapturePolicy::new(args.0, args.1, args.2, args.3, args.4).is_err());
    }
    assert!(ProcessId::new("").is_err());
    assert!(ProcessId::new("bad\nidentity").is_err());
    assert!(
        ProcessFrame::new(
            process(),
            ProcessStream::Stdout,
            0,
            ProcessEventClass::Critical,
            Severity::Error,
            b"closure".to_vec(),
        )
        .is_err()
    );
    assert!(DurableArtifactPointer::committed("artifact", "not-a-hex-digest").is_err());
}

#[test]
fn critical_frames_backpressure_and_never_become_lossy_gaps() {
    let mut capture = ProcessCapture::new(policy(1, 2, 64, 1));
    assert!(matches!(
        capture.offer(
            frame(
                ProcessStream::Stdout,
                1,
                ProcessEventClass::Critical,
                b"first"
            ),
            CaptureCancellation::Running,
            SinkAvailability::Ready,
        ),
        CaptureDecision::Enqueued { ordinal: 1, .. }
    ));
    assert!(matches!(
        capture.offer(
            frame(
                ProcessStream::Stdout,
                2,
                ProcessEventClass::Critical,
                b"closure"
            ),
            CaptureCancellation::Running,
            SinkAvailability::Ready,
        ),
        CaptureDecision::Backpressure {
            reason: BackpressureReason::QueueCapacity,
            ..
        }
    ));
    assert_eq!(capture.frame_len(), 1);
    assert_eq!(capture.gap_len(), 0);
    assert_eq!(capture.metrics().dropped_frames, 0);
}

#[test]
fn cancelled_critical_backpressure_demotes_only_closure_dimensions() {
    let mut capture = ProcessCapture::new(policy(1, 2, 64, 1));
    let _ = capture.offer(
        frame(
            ProcessStream::Stdout,
            1,
            ProcessEventClass::Diagnostic,
            b"fills",
        ),
        CaptureCancellation::Running,
        SinkAvailability::Ready,
    );
    let decision = capture.offer(
        frame(
            ProcessStream::Stderr,
            1,
            ProcessEventClass::Critical,
            b"proof",
        ),
        CaptureCancellation::Requested,
        SinkAvailability::Ready,
    );
    let CaptureDecision::CriticalFailure { failure, .. } = decision else {
        panic!("cancelled critical pressure must fail closed");
    };
    assert_eq!(
        failure.reason,
        CriticalFailureReason::CancelledWhileBackpressured(BackpressureReason::QueueCapacity)
    );

    let mut outcome = ScopedReceiptOutcome {
        scope: ReceiptScope::Operation,
        receipt: "receipt-1".into(),
        disposition: ExecutionDisposition::Completed,
        predicate: PredicateOutcome::Satisfied,
        evidence_methods: "conformance".into(),
        grade: EpistemicGrade::Verified,
        applicability: DomainApplicability::InDomain,
        support: OperationalSupport::Supported,
        completeness: EvidenceCompleteness::Complete,
        integrity: EvidenceIntegrity::Intact,
        promotion: PromotionEffect::Promoted,
        detail: String::new(),
    };
    failure.impact.apply_to(&mut outcome);
    assert_eq!(outcome.completeness, EvidenceCompleteness::Incomplete);
    assert_eq!(outcome.integrity, EvidenceIntegrity::Unchecked);
    assert_eq!(outcome.promotion, PromotionEffect::Demoted);
    assert_eq!(outcome.disposition, ExecutionDisposition::Completed);
    assert_eq!(outcome.predicate, PredicateOutcome::Satisfied);
}

#[test]
fn diagnostic_queue_loss_is_quantified_and_contiguous_gaps_coalesce() {
    let mut capture = ProcessCapture::new(policy(1, 1, 64, 1));
    let _ = capture.offer(
        frame(
            ProcessStream::Stderr,
            1,
            ProcessEventClass::Diagnostic,
            b"kept",
        ),
        CaptureCancellation::Running,
        SinkAvailability::Ready,
    );
    for ordinal in 2..=4 {
        assert_eq!(
            capture.offer(
                frame(
                    ProcessStream::Stderr,
                    ordinal,
                    ProcessEventClass::Diagnostic,
                    b"dropped",
                ),
                CaptureCancellation::Running,
                SinkAvailability::Ready,
            ),
            CaptureDecision::Dropped {
                ordinal,
                reason: LossReason::QueueCapacity,
            }
        );
    }
    assert_eq!(capture.gap_len(), 1);
    let gap = capture.pop_gap().expect("coalesced gap");
    assert_eq!((gap.first_ordinal, gap.last_ordinal), (2, 4));
    assert_eq!(
        (gap.original_count, gap.emitted_count, gap.dropped_count),
        (3, 0, 3)
    );
    assert_eq!(gap.omitted_inline_bytes, 21);
    assert_eq!(gap.lost_bytes(), 21);
    assert_eq!(gap.policy_version, 7);
}

#[test]
fn telemetry_sampling_is_deterministic_and_explicit() {
    let mut capture = ProcessCapture::new(policy(8, 8, 64, 3));
    for ordinal in 1..=7 {
        let decision = capture.offer(
            frame(
                ProcessStream::Stdout,
                ordinal,
                ProcessEventClass::Telemetry,
                &[ordinal as u8],
            ),
            CaptureCancellation::Running,
            SinkAvailability::Ready,
        );
        if matches!(ordinal, 1 | 4 | 7) {
            assert!(matches!(decision, CaptureDecision::Enqueued { .. }));
        } else {
            assert_eq!(
                decision,
                CaptureDecision::Dropped {
                    ordinal,
                    reason: LossReason::PolicySampling,
                }
            );
        }
    }
    assert_eq!(capture.metrics().enqueued_frames, 3);
    assert_eq!(capture.metrics().dropped_frames, 4);
    assert_eq!(capture.gap_len(), 2, "2-3 and 5-6 form two ranges");
}

#[test]
fn full_gap_ledger_backpressures_instead_of_hiding_loss() {
    let mut capture = ProcessCapture::new(policy(1, 1, 64, 2));
    let _ = capture.offer(
        frame(
            ProcessStream::Stdout,
            1,
            ProcessEventClass::Critical,
            b"fills",
        ),
        CaptureCancellation::Running,
        SinkAvailability::Ready,
    );
    assert!(matches!(
        capture.offer(
            frame(
                ProcessStream::Stdout,
                2,
                ProcessEventClass::Telemetry,
                b"sampled"
            ),
            CaptureCancellation::Running,
            SinkAvailability::Ready,
        ),
        CaptureDecision::Dropped { .. }
    ));
    assert!(matches!(
        capture.offer(
            frame(
                ProcessStream::Stderr,
                1,
                ProcessEventClass::Diagnostic,
                b"different-gap",
            ),
            CaptureCancellation::Running,
            SinkAvailability::Failed(SinkFailure::BrokenPipe),
        ),
        CaptureDecision::Backpressure {
            reason: BackpressureReason::GapLedgerCapacity,
            ..
        }
    ));
    assert_eq!(capture.metrics().dropped_frames, 1);
}

#[test]
fn oversized_critical_payload_requires_committed_durable_detail() {
    let mut capture = ProcessCapture::new(policy(2, 2, 4, 1));
    assert!(matches!(
        capture.offer(
            frame(
                ProcessStream::Stdout,
                1,
                ProcessEventClass::Critical,
                b"proof-payload",
            ),
            CaptureCancellation::Running,
            SinkAvailability::Ready,
        ),
        CaptureDecision::Backpressure {
            reason: BackpressureReason::DurableSpillRequired,
            ..
        }
    ));

    let pointer = DurableArtifactPointer::committed(
        "artifact://process/worker-7/stdout/1",
        "0123456789abcdef0123456789abcdef",
    )
    .expect("committed fixture pointer");
    let durable = frame(
        ProcessStream::Stdout,
        1,
        ProcessEventClass::Critical,
        b"proof-payload",
    )
    .with_committed_artifact(pointer.clone());
    assert_eq!(
        capture.offer(
            durable,
            CaptureCancellation::Running,
            SinkAvailability::Ready,
        ),
        CaptureDecision::Enqueued {
            ordinal: 1,
            omitted_inline_bytes: 9,
            durably_spilled: true,
        }
    );
    let queued = capture.pop_frame().expect("durable queue record");
    assert_eq!(queued.payload(), b"proo");
    assert_eq!(queued.original_bytes(), 13);
    assert_eq!(queued.artifact(), Some(&pointer));
    let gap = capture.pop_gap().expect("inline spill accounting");
    assert_eq!(gap.reason, LossReason::DurableSpill);
    assert_eq!(gap.omitted_inline_bytes, 9);
    assert_eq!(gap.lost_bytes(), 0);
}

#[test]
fn diagnostic_truncation_without_artifact_declares_byte_loss() {
    let mut capture = ProcessCapture::new(policy(2, 2, 3, 1));
    assert!(matches!(
        capture.offer(
            frame(
                ProcessStream::Stderr,
                1,
                ProcessEventClass::Diagnostic,
                b"abcdef",
            ),
            CaptureCancellation::Running,
            SinkAvailability::Ready,
        ),
        CaptureDecision::Enqueued {
            omitted_inline_bytes: 3,
            durably_spilled: false,
            ..
        }
    ));
    assert_eq!(
        capture.pop_frame().expect("queued prefix").payload(),
        b"abc"
    );
    let gap = capture.pop_gap().expect("explicit truncation gap");
    assert_eq!(gap.reason, LossReason::InlineLimit);
    assert_eq!((gap.emitted_count, gap.dropped_count), (1, 0));
    assert_eq!(gap.lost_bytes(), 3);
}

#[test]
fn hostile_binary_payloads_are_opaque_and_stream_ordinals_are_independent() {
    let mut capture = ProcessCapture::new(policy(4, 4, 64, 1));
    let hostile = [0, 0xff, b'\n', b'\r', b'"', b'\\'];
    let _ = capture.offer(
        frame(
            ProcessStream::Stdout,
            1,
            ProcessEventClass::Diagnostic,
            &hostile,
        ),
        CaptureCancellation::Running,
        SinkAvailability::Ready,
    );
    let _ = capture.offer(
        frame(
            ProcessStream::Stderr,
            1,
            ProcessEventClass::Diagnostic,
            b"stderr",
        ),
        CaptureCancellation::Running,
        SinkAvailability::Ready,
    );
    assert_eq!(
        capture.pop_frame().expect("stdout frame").payload(),
        hostile
    );
    assert_eq!(
        capture.pop_frame().expect("stderr frame").stream(),
        ProcessStream::Stderr
    );

    let CaptureDecision::Rejected { error, .. } = capture.offer(
        frame(
            ProcessStream::Stdout,
            1,
            ProcessEventClass::Diagnostic,
            b"fork",
        ),
        CaptureCancellation::Running,
        SinkAvailability::Ready,
    ) else {
        panic!("duplicate ordinal must refuse");
    };
    assert_eq!((error.previous, error.received), (1, 1));
}

#[test]
fn captured_binary_frame_projects_to_one_typed_canonical_event() {
    let pointer = DurableArtifactPointer::committed(
        "artifact://process/worker-7/stderr/1",
        "0123456789abcdef0123456789abcdef",
    )
    .expect("committed fixture pointer");
    let frame = ProcessFrame::new(
        process(),
        ProcessStream::Stderr,
        1,
        ProcessEventClass::Critical,
        Severity::Error,
        vec![0, 0xff, b'\n'],
    )
    .expect("valid binary frame")
    .with_committed_artifact(pointer);
    let mut emitter = Emitter::new("capture-session", "worker-7");
    let event = frame.into_event(&mut emitter, Some(99));

    assert_eq!(event.severity, Severity::Error);
    assert!(matches!(&event.kind, EventKind::ProcessFrame { .. }));
    assert_eq!(event.kind.kind_name(), "process_frame");
    lint_failure_record(&event).expect("typed frame stays lint-valid");
    let line = event.to_jsonl();
    validate_line(&line).expect("typed frame stays schema-valid");
    assert_eq!(
        line,
        "{\"v\":1,\"session\":\"capture-session\",\"scope\":\"worker-7\",\"seq\":0,\
         \"severity\":\"error\",\"kind\":\"process_frame\",\"payload\":{\
         \"process\":\"worker-7\",\"stream\":\"stderr\",\"ordinal\":1,\"class\":\"critical\",\
         \"payload_hex\":\"00ff0a\",\"original_bytes\":3,\
         \"artifact\":\"artifact://process/worker-7/stderr/1\",\
         \"artifact_content_hash\":\"0123456789abcdef0123456789abcdef\"},\"wall_ns\":99}"
    );
}

#[test]
fn quantified_gap_projects_counts_policy_and_explicit_absent_artifact() {
    let mut capture = ProcessCapture::new(policy(2, 2, 3, 1));
    let decision = capture.offer(
        frame(
            ProcessStream::Stdout,
            1,
            ProcessEventClass::Diagnostic,
            b"abcdef",
        ),
        CaptureCancellation::Running,
        SinkAvailability::Ready,
    );
    assert!(matches!(decision, CaptureDecision::Enqueued { .. }));
    let gap = capture.pop_gap().expect("inline truncation gap");
    let mut emitter = Emitter::new("capture-session", "worker-7");
    let event = gap.into_event(&mut emitter, None);

    assert_eq!(event.severity, Severity::Warn);
    assert!(matches!(&event.kind, EventKind::ProcessGap { .. }));
    assert_eq!(event.kind.kind_name(), "process_gap");
    lint_failure_record(&event).expect("balanced gap stays lint-valid");
    let line = event.to_jsonl();
    validate_line(&line).expect("typed gap stays schema-valid");
    assert_eq!(
        line,
        "{\"v\":1,\"session\":\"capture-session\",\"scope\":\"worker-7\",\"seq\":0,\
         \"severity\":\"warn\",\"kind\":\"process_gap\",\"payload\":{\
         \"process\":\"worker-7\",\"stream\":\"stdout\",\"class\":\"diagnostic\",\
         \"first_ordinal\":1,\"last_ordinal\":1,\"original_count\":1,\
         \"emitted_count\":1,\"dropped_count\":0,\"omitted_inline_bytes\":3,\
         \"reason\":\"inline_limit\",\"policy_version\":7,\"artifact\":null,\
         \"artifact_content_hash\":null}}"
    );
}

#[test]
fn malformed_typed_process_records_fail_the_event_lint() {
    let mut emitter = Emitter::new("capture-session", "worker-7");
    let gap = |original_count, emitted_count, dropped_count, reason, policy_version| {
        EventKind::ProcessGap {
            process: process(),
            stream: ProcessStream::Stdout,
            class: ProcessEventClass::Diagnostic,
            first_ordinal: 1,
            last_ordinal: 2,
            original_count,
            emitted_count,
            dropped_count,
            omitted_inline_bytes: 3,
            reason,
            policy_version,
            artifact: None,
        }
    };
    for kind in [
        gap(1, 1, 0, LossReason::InlineLimit, 7),
        gap(2, 1, 1, LossReason::InlineLimit, 0),
        gap(2, 2, 0, LossReason::DurableSpill, 7),
    ] {
        let event = emitter.emit(Severity::Warn, kind, None);
        assert!(
            lint_failure_record(&event).is_err(),
            "malformed typed process evidence must fail closed: {event:?}"
        );
    }

    let frame = emitter.emit(
        Severity::Info,
        EventKind::ProcessFrame {
            process: process(),
            stream: ProcessStream::Stderr,
            ordinal: 1,
            class: ProcessEventClass::Diagnostic,
            payload: b"four".to_vec(),
            original_bytes: 3,
            artifact: None,
        },
        None,
    );
    assert!(lint_failure_record(&frame).is_err());
}

#[test]
fn sink_failures_are_lossy_only_for_noncritical_classes() {
    let mut capture = ProcessCapture::new(policy(2, 2, 64, 1));
    assert!(matches!(
        capture.offer(
            frame(
                ProcessStream::Stderr,
                1,
                ProcessEventClass::Critical,
                b"proof"
            ),
            CaptureCancellation::Running,
            SinkAvailability::Failed(SinkFailure::DiskExhausted),
        ),
        CaptureDecision::CriticalFailure {
            failure: fs_obs::process::CriticalCaptureFailure {
                reason: CriticalFailureReason::Sink(SinkFailure::DiskExhausted),
                ..
            },
            ..
        }
    ));
    assert_eq!(
        capture.offer(
            frame(
                ProcessStream::Stdout,
                1,
                ProcessEventClass::Diagnostic,
                b"detail"
            ),
            CaptureCancellation::Running,
            SinkAvailability::Failed(SinkFailure::BrokenPipe),
        ),
        CaptureDecision::Dropped {
            ordinal: 1,
            reason: LossReason::BrokenPipe,
        }
    );
    assert_eq!(
        capture.pop_gap().expect("sink gap").reason,
        LossReason::BrokenPipe
    );
}

#[test]
fn termination_never_substitutes_for_a_final_receipt() {
    let missing = ProcessClosure::reconcile(process(), ProcessTermination::Exited(0), None);
    let ProcessClosure::ObservationGap { impact, .. } = missing else {
        panic!("zero exit without receipt must remain a gap");
    };
    assert_eq!(impact.completeness, EvidenceCompleteness::Incomplete);
    assert_eq!(impact.integrity, EvidenceIntegrity::Unchecked);
    assert_eq!(impact.promotion, PromotionEffect::Demoted);

    let receipt = FinalReceiptId::new("operation-receipt-deadbeef").expect("receipt id");
    let reconciled = ProcessClosure::reconcile(
        process(),
        ProcessTermination::Signaled(15),
        Some(receipt.clone()),
    );
    assert!(matches!(
        reconciled,
        ProcessClosure::Reconciled {
            termination: ProcessTermination::Signaled(15),
            final_receipt,
            ..
        } if final_receipt == receipt
    ));
}
