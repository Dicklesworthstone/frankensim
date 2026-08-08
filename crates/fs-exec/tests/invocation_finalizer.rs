//! G0/G4/G5 coverage for the affine post-cancel finalizer authority.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_domain;
use fs_exec::{
    Budget, CancelGate, CostUnits, Cx, EvaluationUnits, ExecMode, FinalizationObservation,
    FinalizationPublication, FinalizationReport, FinalizationResources, InvocationAdmitter,
    InvocationDisposition, InvocationError, InvocationLimits, InvocationPlanBinding,
    InvocationPublicationScope, InvocationReceipt, InvocationResources, MemoryBytes, OutputBytes,
    PollUnits, StreamKey, Time, VirtualClock, WorkUnits,
};

const STREAM: StreamKey = StreamKey {
    seed: 0xF1A1_1E20_2607_28,
    kernel_id: 0xF1A1,
    tile: 0,
    iteration: 0,
};

fn resources(work: u128, polls: u32, output: u64) -> InvocationResources {
    InvocationResources::new(
        WorkUnits::new(work),
        PollUnits::new(polls),
        CostUnits::new(0),
        EvaluationUnits::new(0),
        MemoryBytes::new(0),
        OutputBytes::new(output),
    )
}

fn identities(
    case: &str,
) -> (
    fs_blake3::ContentHash,
    fs_blake3::ContentHash,
    fs_blake3::ContentHash,
) {
    (
        hash_domain(
            "frankensim.fs-exec.test.finalizer.invocation.v1",
            case.as_bytes(),
        ),
        hash_domain("frankensim.fs-exec.test.finalizer.accuracy.v1", b"exact"),
        hash_domain(
            "frankensim.fs-exec.test.finalizer.capability.v1",
            b"unit-test",
        ),
    )
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>, &VirtualClock) -> R) -> R {
    let arenas = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = arenas.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            STREAM,
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx, &clock)
    });
    assert!(
        arenas.stats().quiescent(),
        "finalizer test leaked arena state: {}",
        arenas.stats().to_json()
    );
    result
}

#[test]
fn checked_plan_binding_separates_child_report_and_invocation_evidence() {
    fn run(plan: Option<InvocationPlanBinding>) -> (FinalizationReport, InvocationReceipt) {
        let gate = CancelGate::new_clock_free();
        let required = resources(0, 0, 0);
        let (invocation, accuracy, capability) = identities("plan-binding");
        with_cx(&gate, |cx, clock| {
            let admitter = InvocationAdmitter::new();
            let admission = match plan {
                Some(binding) => admitter.admit_bound(
                    invocation,
                    binding,
                    InvocationLimits::new(required, None, accuracy, capability),
                    required,
                ),
                None => admitter.admit(
                    invocation,
                    InvocationLimits::new(required, None, accuracy, capability),
                    required,
                ),
            }
            .expect("the zero-resource fixture admits");
            let mut root = admission.begin(cx, clock).expect("invocation begins");
            let child = root
                .split_finalizable_child(
                    "bound-plan-operation",
                    required,
                    FinalizationResources::default(),
                )
                .expect("finalizable child admits");
            let mut finalizer = child.begin_finalization();
            finalizer
                .abort_child_local_publication()
                .expect("zero-output publication seals");
            let report = finalizer.finish().expect("child evidence seals");
            let receipt = root.finish().expect("root evidence seals");
            (report, receipt)
        })
    }

    let schema_root = hash_domain(
        "frankensim.fs-exec.test.finalizer.plan-schema.v1",
        b"checked-plan",
    );
    let left_binding = InvocationPlanBinding::new(
        schema_root,
        7,
        hash_domain(
            "frankensim.fs-exec.test.finalizer.plan.v1",
            b"left-operation",
        ),
    );
    let right_binding = InvocationPlanBinding::new(
        schema_root,
        7,
        hash_domain(
            "frankensim.fs-exec.test.finalizer.plan.v1",
            b"right-operation",
        ),
    );

    let (left_report, left_receipt) = run(Some(left_binding));
    let (left_replay_report, left_replay_receipt) = run(Some(left_binding));
    let (right_report, right_receipt) = run(Some(right_binding));
    let (_, unbound_receipt) = run(None);

    assert_eq!(left_report, left_replay_report);
    assert_eq!(left_receipt, left_replay_receipt);
    assert_eq!(left_report.plan_binding(), Some(left_binding));
    assert_eq!(left_receipt.plan_binding(), Some(left_binding));
    assert_eq!(right_report.plan_binding(), Some(right_binding));
    assert_eq!(right_receipt.plan_binding(), Some(right_binding));
    assert_eq!(unbound_receipt.plan_binding(), None);

    assert_ne!(left_report.child(), right_report.child());
    assert_ne!(left_report.root(), right_report.root());
    assert_ne!(left_receipt.root(), right_receipt.root());
    assert!(left_report.join(&left_receipt).is_ok());
    assert!(matches!(
        left_report.join(&right_receipt),
        Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "plan-binding"
        })
    ));
    assert!(matches!(
        left_report.join(&unbound_receipt),
        Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "plan-binding"
        })
    ));
}

#[test]
fn successful_child_can_publish_only_after_finalizer_poll() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(3, 2, 8);
    let finalization = FinalizationResources::new(WorkUnits::new(2), PollUnits::new(1));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .expect("fixture total is representable");
    let (invocation, accuracy, capability) = identities("success");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .expect("complete plan admits");
        let mut root = admission.begin(cx, clock).expect("invocation begins");
        let mut child = root
            .split_finalizable_child("sparse-success", scientific, finalization)
            .expect("scientific and cleanup resources transfer once");
        child.charge_work(WorkUnits::new(3)).unwrap();
        child.poll("sparse-success.work").unwrap();

        let mut finalizer = child.begin_finalization();
        finalizer.charge_cleanup_work(WorkUnits::new(2)).unwrap();
        let prepared = finalizer
            .prepare_publication()
            .expect("final poll permits preparation");
        let mut destination = 7_u64;
        let replaced = finalizer
            .commit_child_local_publication(prepared, OutputBytes::new(8), &mut destination, 11_u64)
            .expect("prepared output commits once");
        assert_eq!(replaced, 7);
        assert_eq!(destination, 11);
        let report = finalizer.finish().expect("finalizer seals");
        let receipt = root.finish().expect("root seals after child");
        (report, receipt)
    });

    assert_eq!(report.disposition(), InvocationDisposition::Completed);
    assert_eq!(
        report.publication_scope(),
        InvocationPublicationScope::ChildLocal
    );
    assert_eq!(
        report.publication(),
        FinalizationPublication::Committed {
            bytes: OutputBytes::new(8)
        }
    );
    assert_eq!(report.consumed().work(), WorkUnits::new(2));
    assert_eq!(report.consumed().polls(), PollUnits::new(1));
    assert!(report.verifies_integrity());
    assert!(receipt.verifies_integrity());
    assert_eq!(receipt.output_retained_bytes(), 8);
    let joined = report.join(&receipt).expect("exact child joins");
    assert!(joined.verifies_integrity());
}

#[test]
fn finalizable_transaction_rejects_nested_output_authority_before_mutation() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 1);
    let finalization = FinalizationResources::new(WorkUnits::new(0), PollUnits::new(0));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let (invocation, accuracy, capability) = identities("nested-output-bypass");

    let receipt = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let mut child = root
            .split_finalizable_child("outer-transaction", scientific, finalization)
            .unwrap();
        let expected_ancestor = child.id();
        assert!(matches!(
            child.split_child("nested-output", resources(0, 0, 1)),
            Err(InvocationError::TransactionalOutputScopeViolation {
                ancestor,
                phase: "nested-output",
                requested: 1,
            }) if ancestor == expected_ancestor
        ));

        let mut finalizer = child.begin_finalization();
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        assert_eq!(report.disposition(), InvocationDisposition::Refused);
        assert_eq!(report.publication(), FinalizationPublication::Aborted);
        root.finish().unwrap()
    });

    assert_eq!(receipt.children().len(), 1);
    assert_eq!(receipt.output_retained_bytes(), 0);
    assert!(receipt.verifies_integrity());
}

#[test]
fn publication_refuses_abandoned_descendants_before_poll_or_destination_mutation() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 1);
    let finalization = FinalizationResources::new(WorkUnits::new(0), PollUnits::new(1));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.abandoned-descendant.v1",
        b"nested-unwind",
    );
    let (invocation, accuracy, capability) = identities("abandoned-descendant-publication");

    let receipt = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let outer_id;
        {
            let mut outer = root
                .split_finalizable_child("outer", scientific, finalization)
                .unwrap();
            outer_id = outer.id();
            let nested = outer
                .split_child("abandoned-nested", resources(0, 0, 0))
                .unwrap();
            drop(nested);

            let mut finalizer = outer.begin_finalization();
            assert!(matches!(
                finalizer.prepare_publication(),
                Err(InvocationError::LiveNestedChildren { count: 1 })
            ));
            assert_eq!(finalizer.remaining().polls(), PollUnits::new(1));
        }

        let nested = root.next_unfinished_child().unwrap();
        assert_ne!(nested.id(), outer_id);
        let recovered = root
            .recover_child_budget(nested.id(), "abandoned-nested.recover", reason)
            .unwrap();
        assert_eq!(recovered.finish().unwrap(), InvocationDisposition::Refused);

        let mut finalizer = root
            .recover_child_finalizer(outer_id, "outer.recover", reason)
            .unwrap();
        finalizer.abort_publication().unwrap();
        assert_eq!(
            finalizer.finish().unwrap().disposition(),
            InvocationDisposition::Refused
        );
        root.finish().unwrap()
    });

    assert_eq!(receipt.output_retained_bytes(), 0);
    assert!(receipt.verifies_integrity());
}

#[test]
fn cancellation_keeps_cleanup_authority_live_and_preserves_first_phase() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(4, 2, 4);
    let finalization = FinalizationResources::new(WorkUnits::new(2), PollUnits::new(1));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let (invocation, accuracy, capability) = identities("cancel");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let mut child = root
            .split_finalizable_child("sparse-cancel", scientific, finalization)
            .unwrap();
        child.charge_work(WorkUnits::new(1)).unwrap();
        gate.request();
        assert!(matches!(
            child.poll("sparse-cancel.work-poll"),
            Err(InvocationError::Cancelled {
                phase: "sparse-cancel.work-poll"
            })
        ));

        let mut finalizer = child.begin_finalization();
        finalizer.charge_cleanup_work(WorkUnits::new(2)).unwrap();
        assert_eq!(
            finalizer.poll_cleanup("sparse-cancel.cleanup").unwrap(),
            FinalizationObservation::Terminal(InvocationDisposition::Cancelled)
        );
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.disposition(), InvocationDisposition::Cancelled);
    assert_eq!(report.publication(), FinalizationPublication::Aborted);
    assert!(matches!(
        report.failure(),
        Some(InvocationError::Cancelled {
            phase: "sparse-cancel.work-poll"
        })
    ));
    assert_eq!(receipt.output_retained_bytes(), 0);
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn request_observed_at_finalization_begin_cannot_become_success() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 1);
    let finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let (invocation, accuracy, capability) = identities("cancel-at-transition");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("transition", scientific, finalization)
            .unwrap();
        gate.request();
        let mut finalizer = child.begin_finalization();
        finalizer.charge_cleanup_work(WorkUnits::new(1)).unwrap();
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert!(matches!(
        report.failure(),
        Some(InvocationError::Cancelled {
            phase: "child-finalization-begin"
        })
    ));
    assert_eq!(report.disposition(), InvocationDisposition::Cancelled);
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn finalization_overrun_refuses_but_does_not_disable_abort_and_seal() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 0);
    let finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let required = finalization.as_invocation_resources();
    let (invocation, accuracy, capability) = identities("cleanup-overrun");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("cleanup-overrun", scientific, finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        assert!(matches!(
            finalizer.charge_cleanup_work(WorkUnits::new(2)),
            Err(InvocationError::ResourceExceeded {
                resource: "finalization-work",
                requested: 2,
                available: 1,
            })
        ));
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.disposition(), InvocationDisposition::Refused);
    assert_eq!(report.consumed().work(), WorkUnits::new(0));
    assert_eq!(report.returned().work(), WorkUnits::new(1));
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn cancellation_racing_after_prepare_forces_unchanged_publication() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 16);
    let finalization = FinalizationResources::new(WorkUnits::new(0), PollUnits::new(1));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let (invocation, accuracy, capability) = identities("prepare-race");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("prepare-race", scientific, finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        let prepared = finalizer.prepare_publication().unwrap();
        gate.request();
        let mut destination = 7_u64;
        let error = finalizer
            .commit_publication(prepared, OutputBytes::new(16), &mut destination, 11_u64)
            .unwrap_err();
        assert_eq!(error.error(), &InvocationError::PublicationForbidden);
        assert_eq!(error.into_parts().1, 11);
        assert_eq!(destination, 7);
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.publication(), FinalizationPublication::Aborted);
    assert_eq!(report.disposition(), InvocationDisposition::Cancelled);
    assert_eq!(receipt.output_retained_bytes(), 0);
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn zero_cleanup_reserve_is_legal_for_a_no_output_terminal_step() {
    let gate = CancelGate::new_clock_free();
    let required = resources(0, 0, 0);
    let finalization = FinalizationResources::default();
    let (invocation, accuracy, capability) = identities("zero");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("zero", required, finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.granted(), FinalizationResources::default());
    assert_eq!(report.consumed(), FinalizationResources::default());
    assert_eq!(report.returned(), FinalizationResources::default());
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn finalization_and_join_roots_replay_deterministically() {
    fn run() -> (
        fs_blake3::ContentHash,
        fs_blake3::ContentHash,
        fs_blake3::ContentHash,
    ) {
        let gate = CancelGate::new_clock_free();
        let scientific = resources(1, 0, 0);
        let finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
        let required = scientific
            .checked_add(finalization.as_invocation_resources())
            .unwrap();
        let (invocation, accuracy, capability) = identities("replay");
        with_cx(&gate, |cx, clock| {
            let admission = InvocationAdmitter::new()
                .admit(
                    invocation,
                    InvocationLimits::new(required, None, accuracy, capability),
                    required,
                )
                .unwrap();
            let mut root = admission.begin(cx, clock).unwrap();
            let mut child = root
                .split_finalizable_child("replay", scientific, finalization)
                .unwrap();
            child.charge_work(WorkUnits::new(1)).unwrap();
            let mut finalizer = child.begin_finalization();
            finalizer.charge_cleanup_work(WorkUnits::new(1)).unwrap();
            finalizer.abort_publication().unwrap();
            let report = finalizer.finish().unwrap();
            let receipt = root.finish().unwrap();
            let joined = report.join(&receipt).unwrap();
            (report.root(), receipt.root(), joined.root())
        })
    }

    assert_eq!(run(), run());
}

#[test]
fn finalization_report_refuses_join_to_another_invocation() {
    fn one(case: &str) -> (fs_exec::FinalizationReport, fs_exec::InvocationReceipt) {
        let gate = CancelGate::new_clock_free();
        let required = resources(0, 0, 0);
        let (invocation, accuracy, capability) = identities(case);
        with_cx(&gate, |cx, clock| {
            let admission = InvocationAdmitter::new()
                .admit(
                    invocation,
                    InvocationLimits::new(required, None, accuracy, capability),
                    required,
                )
                .unwrap();
            let mut root = admission.begin(cx, clock).unwrap();
            let child = root
                .split_finalizable_child("join", required, FinalizationResources::default())
                .unwrap();
            let mut finalizer = child.begin_finalization();
            finalizer.abort_publication().unwrap();
            let report = finalizer.finish().unwrap();
            let receipt = root.finish().unwrap();
            (report, receipt)
        })
    }

    let (left, _) = one("join-left");
    let (_, right) = one("join-right");
    assert!(matches!(
        left.join(&right),
        Err(InvocationError::FinalizationReceiptMismatch {
            invariant: "child-exists"
        })
    ));
}

#[test]
fn child_local_publication_explicitly_survives_later_root_cancellation() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 8);
    let finalization = FinalizationResources::new(WorkUnits::new(0), PollUnits::new(1));
    let finalizable_grant = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let legacy_grant = resources(1, 0, 0);
    let required = finalizable_grant.checked_add(legacy_grant).unwrap();
    let (invocation, accuracy, capability) = identities("root-cancel-after-commit");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("committed-child", scientific, finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        let prepared = finalizer.prepare_publication().unwrap();
        let mut destination = 1_u64;
        finalizer
            .commit_child_local_publication(prepared, OutputBytes::new(8), &mut destination, 2_u64)
            .unwrap();
        assert_eq!(destination, 2);
        let report = finalizer.finish().unwrap();

        let mut legacy = root.split_child("legacy-child", legacy_grant).unwrap();
        legacy.charge_work(WorkUnits::new(1)).unwrap();
        assert_eq!(legacy.finish().unwrap(), InvocationDisposition::Completed);
        gate.request();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.disposition(), InvocationDisposition::Completed);
    assert_eq!(
        report.publication_scope(),
        InvocationPublicationScope::ChildLocal
    );
    assert_eq!(
        report.publication(),
        FinalizationPublication::Committed {
            bytes: OutputBytes::new(8)
        }
    );
    assert_eq!(receipt.disposition(), InvocationDisposition::Cancelled);
    assert!(matches!(
        receipt.failure(),
        Some(InvocationError::Cancelled {
            phase: "invocation-finalize"
        })
    ));
    assert_eq!(receipt.output_retained_bytes(), 8);
    assert!(receipt.verifies_integrity());
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn root_observes_cancel_after_aborted_child_and_keeps_output_empty() {
    let gate = CancelGate::new_clock_free();
    let finalization = FinalizationResources::default();
    let required = resources(1, 0, 0);
    let (invocation, accuracy, capability) = identities("root-cancel-after-abort");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("aborted-child", resources(0, 0, 0), finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();

        let mut legacy = root.split_child("legacy-after-abort", required).unwrap();
        legacy.charge_work(WorkUnits::new(1)).unwrap();
        assert_eq!(legacy.finish().unwrap(), InvocationDisposition::Completed);
        gate.request();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.disposition(), InvocationDisposition::Completed);
    assert_eq!(report.publication(), FinalizationPublication::Aborted);
    assert_eq!(receipt.disposition(), InvocationDisposition::Cancelled);
    assert!(matches!(
        receipt.failure(),
        Some(InvocationError::Cancelled {
            phase: "invocation-finalize"
        })
    ));
    assert_eq!(receipt.output_retained_bytes(), 0);
    assert!(receipt.verifies_integrity());
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn deadline_crossing_between_prepare_and_atomic_swap_aborts_without_mutation() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 8);
    let finalization = FinalizationResources::new(WorkUnits::new(0), PollUnits::new(1));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let (invocation, accuracy, capability) = identities("deadline-at-commit");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, Some(Time::from_nanos(5)), accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("deadline-child", scientific, finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        let prepared = finalizer.prepare_publication().unwrap();
        clock.advance(5);
        let mut destination = 13_u64;
        let error = finalizer
            .commit_publication(prepared, OutputBytes::new(8), &mut destination, 21_u64)
            .unwrap_err();
        assert_eq!(error.error(), &InvocationError::PublicationForbidden);
        assert_eq!(error.into_parts().1, 21);
        assert_eq!(destination, 13);
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.publication(), FinalizationPublication::Aborted);
    assert_eq!(report.disposition(), InvocationDisposition::Cancelled);
    assert!(matches!(
        report.failure(),
        Some(InvocationError::DeadlineExpired {
            phase: "child-finalization-commit",
            deadline_ns: 5,
            observed_ns: 5,
        })
    ));
    assert_eq!(receipt.output_retained_bytes(), 0);
    assert!(receipt.verifies_integrity());
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn premature_finish_is_non_consuming_and_cleanup_can_continue() {
    let gate = CancelGate::new_clock_free();
    let finalization = FinalizationResources::new(WorkUnits::new(2), PollUnits::new(0));
    let required = finalization.as_invocation_resources();
    let (invocation, accuracy, capability) = identities("premature-finish");

    let (first, replay, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("premature", resources(0, 0, 0), finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        assert_eq!(
            finalizer.finish(),
            Err(InvocationError::FinalizationIncomplete {
                step: "publication-seal"
            })
        );
        finalizer.charge_cleanup_work(WorkUnits::new(2)).unwrap();
        finalizer.abort_publication().unwrap();
        let first = finalizer.finish().unwrap();
        let replay = finalizer.finish().unwrap();
        assert_eq!(first, replay, "replay cannot return resources twice");
        let replay_from_root = root.replay_child_finalization(first.child()).unwrap();
        assert_eq!(first, replay_from_root);
        let receipt = root.finish().unwrap();
        (first, replay, receipt)
    });

    assert_eq!(first, replay);
    assert_eq!(first.consumed().work(), WorkUnits::new(2));
    assert!(receipt.verifies_integrity());
    assert!(first.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn sealed_publication_rejects_later_cleanup_without_relabeling_child() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(0, 0, 8);
    let finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(2));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let (invocation, accuracy, capability) = identities("post-commit-misuse");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("sealed", scientific, finalization)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        finalizer.charge_cleanup_work(WorkUnits::new(1)).unwrap();
        let prepared = finalizer.prepare_publication().unwrap();
        let mut destination = 1_u64;
        assert_eq!(
            finalizer
                .commit_publication(prepared, OutputBytes::new(8), &mut destination, 2_u64)
                .unwrap(),
            1
        );
        gate.request();
        assert_eq!(
            finalizer.poll_cleanup("must-not-observe-post-commit"),
            Err(InvocationError::PublicationAlreadySealed)
        );
        assert_eq!(
            finalizer.charge_cleanup_work(WorkUnits::new(1)),
            Err(InvocationError::PublicationAlreadySealed)
        );
        assert_eq!(
            finalizer.abort_publication(),
            Err(InvocationError::PublicationAlreadySealed)
        );
        assert_eq!(destination, 2);
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert_eq!(report.disposition(), InvocationDisposition::Completed);
    assert!(report.failure().is_none());
    assert_eq!(
        report.publication(),
        FinalizationPublication::Committed {
            bytes: OutputBytes::new(8)
        }
    );
    assert_eq!(receipt.disposition(), InvocationDisposition::Cancelled);
    assert!(receipt.verifies_integrity());
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn maximum_cleanup_reserve_round_trips_and_combined_overflow_latches() {
    let gate = CancelGate::new_clock_free();
    let maximum = FinalizationResources::new(WorkUnits::new(u128::MAX), PollUnits::new(u32::MAX));
    let required = maximum.as_invocation_resources();
    let (invocation, accuracy, capability) = identities("maximum-cleanup");

    let report = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("maximum", resources(0, 0, 0), maximum)
            .unwrap();
        let mut finalizer = child.begin_finalization();
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        assert!(report.join(&receipt).unwrap().verifies_integrity());
        report
    });
    assert_eq!(report.granted(), maximum);
    assert_eq!(report.returned(), maximum);
    assert_eq!(report.consumed(), FinalizationResources::default());

    let gate = CancelGate::new_clock_free();
    let (invocation, accuracy, capability) = identities("combined-overflow");
    let receipt = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        assert!(matches!(
            root.split_finalizable_child("overflow", resources(1, 0, 0), maximum),
            Err(InvocationError::ArithmeticOverflow { resource: "work" })
        ));
        root.finish().unwrap()
    });
    assert!(matches!(
        receipt.failure(),
        Some(InvocationError::ArithmeticOverflow { resource: "work" })
    ));
    assert_eq!(receipt.disposition(), InvocationDisposition::Refused);
    assert!(receipt.verifies_integrity());
}

#[test]
fn caught_scientific_unwind_recovers_fail_closed_cleanup_authority() {
    let gate = CancelGate::new_clock_free();
    let scientific = resources(2, 0, 0);
    let finalization = FinalizationResources::new(WorkUnits::new(2), PollUnits::new(0));
    let required = scientific
        .checked_add(finalization.as_invocation_resources())
        .unwrap();
    let reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.unwind-reason.v1",
        b"scientific",
    );
    let (invocation, accuracy, capability) = identities("scientific-unwind");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("scientific-unwind", scientific, finalization)
            .unwrap();
        let child_id = child.id();
        let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut child = child;
            child.charge_work(WorkUnits::new(1)).unwrap();
            panic!("seeded scientific unwind");
        }));
        assert!(unwind.is_err());

        let mut finalizer = root
            .recover_child_finalizer(child_id, "scientific-unwind.recover", reason)
            .unwrap();
        finalizer.charge_cleanup_work(WorkUnits::new(2)).unwrap();
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert!(matches!(
        report.failure(),
        Some(InvocationError::ExplicitRefusal {
            phase: "scientific-unwind.recover",
            reason: observed,
        }) if *observed == reason
    ));
    assert_eq!(report.disposition(), InvocationDisposition::Refused);
    assert_eq!(report.consumed().work(), WorkUnits::new(2));
    assert!(receipt.verifies_integrity());
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn caught_cleanup_unwind_resumes_remaining_reserve_and_preserves_first_cause() {
    let gate = CancelGate::new_clock_free();
    let finalization = FinalizationResources::new(WorkUnits::new(3), PollUnits::new(0));
    let required = finalization.as_invocation_resources();
    let reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.unwind-reason.v1",
        b"cleanup",
    );
    let (invocation, accuracy, capability) = identities("cleanup-unwind");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child("cleanup-unwind", resources(0, 0, 0), finalization)
            .unwrap();
        let child_id = child.id();
        let mut finalizer = child.begin_finalization();
        let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
            finalizer.charge_cleanup_work(WorkUnits::new(1)).unwrap();
            panic!("seeded cleanup unwind");
        }));
        assert!(unwind.is_err());

        let mut recovered = root
            .recover_child_finalizer(child_id, "cleanup-unwind.recover", reason)
            .unwrap();
        assert_eq!(
            recovered.remaining().work(),
            WorkUnits::new(2),
            "already charged cleanup cannot be recreated"
        );
        recovered.charge_cleanup_work(WorkUnits::new(2)).unwrap();
        recovered.abort_publication().unwrap();
        let report = recovered.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert!(matches!(
        report.failure(),
        Some(InvocationError::ExplicitRefusal {
            phase: "cleanup-unwind.recover",
            reason: observed,
        }) if *observed == reason
    ));
    assert_eq!(report.consumed().work(), WorkUnits::new(3));
    assert_eq!(report.disposition(), InvocationDisposition::Refused);
    assert!(receipt.verifies_integrity());
    assert!(report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn abandoned_nested_authority_is_recovered_inside_out_before_parent_seals() {
    let gate = CancelGate::new_clock_free();
    let parent_finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let nested_finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let nested_scientific = resources(1, 0, 0);
    let nested_total = nested_scientific
        .checked_add(nested_finalization.as_invocation_resources())
        .unwrap();
    let parent_scientific = nested_total;
    let required = parent_scientific
        .checked_add(parent_finalization.as_invocation_resources())
        .unwrap();
    let reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.unwind-reason.v1",
        b"nested",
    );
    let (invocation, accuracy, capability) = identities("nested-recovery");

    let (nested_report, parent_report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let mut parent = root
            .split_finalizable_child("parent", parent_scientific, parent_finalization)
            .unwrap();
        let nested = parent
            .split_finalizable_child("nested", nested_scientific, nested_finalization)
            .unwrap();
        let nested_id = nested.id();
        std::mem::forget(nested);
        let parent_id = parent.id();
        let mut parent_finalizer = parent.begin_finalization();
        assert_eq!(
            parent_finalizer.finish(),
            Err(InvocationError::LiveNestedChildren { count: 1 })
        );
        drop(parent_finalizer);

        let mut nested_finalizer = root
            .recover_child_finalizer(nested_id, "nested.recover", reason)
            .unwrap();
        nested_finalizer
            .charge_cleanup_work(WorkUnits::new(1))
            .unwrap();
        nested_finalizer.abort_publication().unwrap();
        let nested_report = nested_finalizer.finish().unwrap();
        drop(nested_finalizer);

        let mut parent_finalizer = root
            .recover_child_finalizer(parent_id, "parent.recover", reason)
            .unwrap();
        parent_finalizer
            .charge_cleanup_work(WorkUnits::new(1))
            .unwrap();
        parent_finalizer.abort_publication().unwrap();
        let parent_report = parent_finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (nested_report, parent_report, receipt)
    });

    assert_eq!(nested_report.disposition(), InvocationDisposition::Refused);
    assert_eq!(parent_report.disposition(), InvocationDisposition::Refused);
    assert!(receipt.verifies_integrity());
    assert!(nested_report.join(&receipt).unwrap().verifies_integrity());
    assert!(parent_report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn finalization_begin_preserves_deadline_before_requested_gate_precedence() {
    let gate = CancelGate::new_clock_free();
    let required = resources(0, 0, 0);
    let (invocation, accuracy, capability) = identities("begin-precedence");

    let (report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, Some(Time::from_nanos(5)), accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let child = root
            .split_finalizable_child(
                "begin-precedence",
                required,
                FinalizationResources::default(),
            )
            .unwrap();
        clock.advance(5);
        gate.request();
        let mut finalizer = child.begin_finalization();
        finalizer.abort_publication().unwrap();
        let report = finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (report, receipt)
    });

    assert!(matches!(
        report.failure(),
        Some(InvocationError::DeadlineExpired {
            phase: "child-finalization-begin",
            deadline_ns: 5,
            observed_ns: 5,
        })
    ));
    assert_eq!(report.disposition(), InvocationDisposition::Cancelled);
    assert!(receipt.verifies_integrity());
}

#[test]
fn abandoned_ordinary_nested_child_can_be_recovered_fail_closed() {
    let gate = CancelGate::new_clock_free();
    let parent_finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let nested_grant = resources(1, 0, 0);
    let required = nested_grant
        .checked_add(parent_finalization.as_invocation_resources())
        .unwrap();
    let reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.unwind-reason.v1",
        b"ordinary-nested",
    );
    let (invocation, accuracy, capability) = identities("ordinary-nested-recovery");

    let (parent_report, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let mut parent = root
            .split_finalizable_child("parent", nested_grant, parent_finalization)
            .unwrap();
        let nested = parent.split_child("ordinary-nested", nested_grant).unwrap();
        let nested_id = nested.id();
        let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
            let mut nested = nested;
            nested.charge_work(WorkUnits::new(1)).unwrap();
            panic!("seeded ordinary nested unwind");
        }));
        assert!(unwind.is_err());

        let parent_id = parent.id();
        let mut parent_finalizer = parent.begin_finalization();
        assert_eq!(
            parent_finalizer.finish(),
            Err(InvocationError::LiveNestedChildren { count: 1 })
        );
        drop(parent_finalizer);

        let ordinary = root
            .recover_child_budget(nested_id, "ordinary-nested.recover", reason)
            .unwrap();
        assert_eq!(ordinary.finish().unwrap(), InvocationDisposition::Refused);

        let mut parent_finalizer = root
            .recover_child_finalizer(parent_id, "parent.recover", reason)
            .unwrap();
        parent_finalizer
            .charge_cleanup_work(WorkUnits::new(1))
            .unwrap();
        parent_finalizer.abort_publication().unwrap();
        let parent_report = parent_finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (parent_report, receipt)
    });

    assert_eq!(parent_report.disposition(), InvocationDisposition::Refused);
    assert!(matches!(
        parent_report.failure(),
        Some(InvocationError::ExplicitRefusal {
            phase: "ordinary-nested.recover",
            reason: observed,
        }) if *observed == reason
    ));
    assert!(receipt.verifies_integrity());
    assert!(parent_report.join(&receipt).unwrap().verifies_integrity());
}

#[test]
fn root_finish_reports_deepest_abandoned_child_without_consuming_recovery_authority() {
    let gate = CancelGate::new_clock_free();
    let leaf_finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let leaf_total = leaf_finalization.as_invocation_resources();
    let parent_finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let required = leaf_total
        .checked_add(parent_finalization.as_invocation_resources())
        .unwrap();
    let reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.unwind-reason.v1",
        b"discover-without-saved-id",
    );
    let (invocation, accuracy, capability) = identities("discover-unfinished");

    let (leaf_report, parent_report, first, replay) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut parent = root
                .split_finalizable_child("discover-parent", leaf_total, parent_finalization)
                .unwrap();
            let _leaf = parent
                .split_finalizable_child("discover-leaf", resources(0, 0, 0), leaf_finalization)
                .unwrap();
            panic!("panic immediately after nested split, before id capture");
        }));
        assert!(unwind.is_err());

        let leaf = root
            .next_unfinished_child()
            .expect("deepest unfinished child is discoverable");
        assert_eq!(leaf.phase(), "discover-leaf");
        assert_eq!(leaf.live_children(), 0);
        assert_eq!(
            root.finish(),
            Err(InvocationError::UnfinishedChild { child: leaf.id() })
        );
        let mut leaf_finalizer = root
            .recover_child_finalizer(leaf.id(), "discover-leaf.recover", reason)
            .unwrap();
        leaf_finalizer
            .charge_cleanup_work(WorkUnits::new(1))
            .unwrap();
        leaf_finalizer.abort_publication().unwrap();
        let leaf_report = leaf_finalizer.finish().unwrap();
        drop(leaf_finalizer);

        let parent = root.next_unfinished_child().unwrap();
        assert_eq!(parent.phase(), "discover-parent");
        let mut parent_finalizer = root
            .recover_child_finalizer(parent.id(), "discover-parent.recover", reason)
            .unwrap();
        parent_finalizer
            .charge_cleanup_work(WorkUnits::new(1))
            .unwrap();
        parent_finalizer.abort_publication().unwrap();
        let parent_report = parent_finalizer.finish().unwrap();
        drop(parent_finalizer);
        assert!(root.unfinished_children().is_empty());

        let first = root.finish().unwrap();
        let replay = root.finish().unwrap();
        assert_eq!(first, replay);
        assert!(matches!(
            root.split_child("after-seal", resources(0, 0, 0)),
            Err(InvocationError::InvocationAlreadyFinalized)
        ));
        (leaf_report, parent_report, first, replay)
    });

    assert_eq!(first, replay);
    assert!(first.verifies_integrity());
    assert!(leaf_report.join(&first).unwrap().verifies_integrity());
    assert!(parent_report.join(&first).unwrap().verifies_integrity());
}

#[test]
fn recovered_sibling_inherits_one_explicit_failure_origin_and_receipt_verifies() {
    let gate = CancelGate::new_clock_free();
    let a_finalization = FinalizationResources::new(WorkUnits::new(1), PollUnits::new(0));
    let a_total = a_finalization.as_invocation_resources();
    let b_grant = resources(1, 0, 0);
    let required = a_total.checked_add(b_grant).unwrap();
    let b_reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.failure-origin.v1",
        b"ordinary-b",
    );
    let unwind_reason = hash_domain(
        "frankensim.fs-exec.test.finalizer.unwind-reason.v1",
        b"finalizable-a",
    );
    let (invocation, accuracy, capability) = identities("inherited-sibling-failure");

    let (a_report, b_id, receipt) = with_cx(&gate, |cx, clock| {
        let admission = InvocationAdmitter::new()
            .admit(
                invocation,
                InvocationLimits::new(required, None, accuracy, capability),
                required,
            )
            .unwrap();
        let mut root = admission.begin(cx, clock).unwrap();
        let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _a = root
                .split_finalizable_child("abandoned-a", resources(0, 0, 0), a_finalization)
                .unwrap();
            panic!("abandon A before capturing its id");
        }));
        assert!(unwind.is_err());

        let mut b = root.split_child("origin-b", b_grant).unwrap();
        let b_id = b.id();
        assert_eq!(
            b.refuse("origin-b.refusal", b_reason),
            InvocationError::ExplicitRefusal {
                phase: "origin-b.refusal",
                reason: b_reason,
            }
        );
        assert_eq!(b.finish().unwrap(), InvocationDisposition::Refused);

        let a = root.next_unfinished_child().unwrap();
        let mut a_finalizer = root
            .recover_child_finalizer(a.id(), "abandoned-a.recover", unwind_reason)
            .unwrap();
        a_finalizer.charge_cleanup_work(WorkUnits::new(1)).unwrap();
        a_finalizer.abort_publication().unwrap();
        let a_report = a_finalizer.finish().unwrap();
        let receipt = root.finish().unwrap();
        (a_report, b_id, receipt)
    });

    assert_eq!(receipt.failure_origin(), Some(b_id));
    let failed = receipt
        .children()
        .iter()
        .filter(|child| child.failure().is_some())
        .collect::<Vec<_>>();
    assert_eq!(failed.len(), 2);
    assert!(
        !failed
            .iter()
            .find(|child| child.id() == b_id)
            .unwrap()
            .failure_inherited()
    );
    assert!(
        failed
            .iter()
            .find(|child| child.id() != b_id)
            .unwrap()
            .failure_inherited()
    );
    assert!(receipt.verifies_integrity());
    assert!(a_report.join(&receipt).unwrap().verifies_integrity());
}
