//! G0/G4/G5 source-level contract tests for the zero-dependency constellation
//! admission core (bead `frankensim-sj31i.50.1`).

#[path = "../src/constellation_admission.rs"]
mod admission;

use admission::{
    AdmissionBudgets, AdmissionContext, AdmissionContextSpec, AdmissionMachine, AdmissionRule,
    AdmissionState, AdmissionTransition, AdmittedPhase, AnchorObservation, AuthorityId,
    BudgetCharge, CancellationCause, CancellationPhase, CommandClass, ComputeBudget, CxBinding,
    DeadlineBudget, DiagnosticPhase, DrainObligations, ExecutableCapability, ExecutableSlot,
    FetchAuthority, IoBudget, MAX_ADMISSION_BYTES, MAX_ADMISSION_EVENTS, NetworkBudget,
    PathCapability, PathSlot, PublicationAuthority, StateKind, TerminalPhase, TransitionKind,
    TrustAnchorState, transition_kind_may_apply,
};

const DEADLINE: u64 = 100;

fn id(byte: u8) -> AuthorityId {
    AuthorityId::try_from_bytes([byte; 32]).expect("nonzero fixture identity")
}

fn online_budgets() -> AdmissionBudgets {
    AdmissionBudgets::new(
        DeadlineBudget::new(id(4), DEADLINE),
        ComputeBudget {
            work_units: 10,
            memory_bytes: 20,
        },
        IoBudget {
            processes: 2,
            files: 4,
            output_bytes: 100,
        },
        NetworkBudget {
            requests: 3,
            bytes: 200,
        },
        2,
    )
}

fn online_context_with(
    budgets: AdmissionBudgets,
    trust_anchor: TrustAnchorState,
    reverse_capabilities: bool,
) -> AdmissionContext {
    let mut paths = vec![
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::ConstellationLock, id(11)),
        PathCapability::new(PathSlot::DestinationRoot, id(12)),
        PathCapability::new(PathSlot::PublicationTarget, id(13)),
    ];
    let mut executables = vec![ExecutableCapability::new(ExecutableSlot::Git, id(20))];
    if reverse_capabilities {
        paths.reverse();
        executables.reverse();
    }
    AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::Bootstrap,
        fetch: FetchAuthority::PinnedTransport { capability: id(30) },
        publication: PublicationAuthority::BootstrapReceipt { capability: id(13) },
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets,
        trust_anchor,
        path_capabilities: &paths,
        executable_capabilities: &executables,
    })
    .expect("valid online context")
}

fn online_context() -> AdmissionContext {
    online_context_with(
        online_budgets(),
        TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        false,
    )
}

fn offline_verify_context() -> AdmissionContext {
    let paths = [
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::ConstellationLock, id(11)),
    ];
    let executables = [ExecutableCapability::new(ExecutableSlot::Git, id(20))];
    AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::VerifyOnly,
        fetch: FetchAuthority::Offline,
        publication: PublicationAuthority::Prohibited,
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: AdmissionBudgets::new(
            DeadlineBudget::new(id(4), DEADLINE),
            ComputeBudget {
                work_units: 10,
                memory_bytes: 20,
            },
            IoBudget {
                processes: 2,
                files: 4,
                output_bytes: 100,
            },
            NetworkBudget {
                requests: 0,
                bytes: 0,
            },
            2,
        ),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &paths,
        executable_capabilities: &executables,
    })
    .expect("valid offline context")
}

fn preflight(machine: &mut AdmissionMachine) {
    machine
        .apply(AdmissionTransition::Preflight {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            at_tick: 1,
        })
        .expect("preflight");
}

fn admit(machine: &mut AdmissionMachine) {
    preflight(machine);
    machine
        .apply(AdmissionTransition::StabilityRecheck {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            at_tick: 2,
        })
        .expect("stable recheck");
}

fn assert_rejected_without_mutation(
    machine: &mut AdmissionMachine,
    transition: AdmissionTransition,
    rule: AdmissionRule,
) {
    let before = machine.encode_canonical().expect("encode before refusal");
    let before_state = machine.state();
    let before_consumption = machine.consumption();
    let before_events = machine.events().len();
    let error = machine
        .apply(transition)
        .expect_err("transition must refuse");
    assert_eq!(error.rule(), rule);
    assert_eq!(machine.state(), before_state);
    assert_eq!(machine.consumption(), before_consumption);
    assert_eq!(machine.events().len(), before_events);
    assert_eq!(
        machine.encode_canonical().expect("encode after refusal"),
        before,
        "a refused transition must be transactionally invisible"
    );
}

#[test]
fn g0_legal_path_binds_preflight_work_recheck_and_one_publication() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    assert_eq!(
        machine.state(),
        AdmissionState::Diagnostic(DiagnosticPhase::Created)
    );

    admit(&mut machine);
    for charge in [
        BudgetCharge::Work(10),
        BudgetCharge::Memory(20),
        BudgetCharge::Processes(2),
        BudgetCharge::Files(4),
        BudgetCharge::Output(100),
        BudgetCharge::Network {
            requests: 3,
            bytes: 200,
        },
    ] {
        machine
            .apply(AdmissionTransition::Charge { charge, at_tick: 3 })
            .expect("exact-cap charge");
    }
    machine
        .apply(AdmissionTransition::BeginPublication {
            quiescence: id(70),
            at_tick: 4,
        })
        .expect("close work before publication");
    machine
        .apply(AdmissionTransition::PublicationRecheck {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            fence: id(71),
            at_tick: 5,
        })
        .expect("publication recheck");
    machine
        .apply(AdmissionTransition::AuthorizePublication {
            receipt: id(60),
            at_tick: 6,
        })
        .expect("publication authorization");
    machine
        .apply(AdmissionTransition::FinalizePublication)
        .expect("publication finalization");

    assert_eq!(
        machine.state(),
        AdmissionState::Admitted(AdmittedPhase::Published { receipt: id(60) })
    );
    assert!(machine.state().attempt_is_finalized());
    for (sequence, event) in machine.events().iter().enumerate() {
        assert_eq!(usize::from(event.sequence()), sequence);
        assert_eq!(event.attempt(), 0);
    }
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::FinalizePublication,
        AdmissionRule::IllegalTransition,
    );
}

#[test]
fn g0_every_budget_is_exact_cap_and_n_plus_one_is_transactional() {
    let cases = [
        (BudgetCharge::Work(10), BudgetCharge::Work(1)),
        (BudgetCharge::Memory(20), BudgetCharge::Memory(1)),
        (BudgetCharge::Processes(2), BudgetCharge::Processes(1)),
        (BudgetCharge::Files(4), BudgetCharge::Files(1)),
        (BudgetCharge::Output(100), BudgetCharge::Output(1)),
        (
            BudgetCharge::Network {
                requests: 3,
                bytes: 200,
            },
            BudgetCharge::Network {
                requests: 1,
                bytes: 1,
            },
        ),
    ];
    for (exact, over) in cases {
        let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
        admit(&mut machine);
        machine
            .apply(AdmissionTransition::Charge {
                charge: exact,
                at_tick: 3,
            })
            .expect("n charge");
        assert_rejected_without_mutation(
            &mut machine,
            AdmissionTransition::Charge {
                charge: over,
                at_tick: 3,
            },
            AdmissionRule::BudgetExceeded,
        );
    }

    let max_budgets = AdmissionBudgets::new(
        DeadlineBudget::new(id(4), DEADLINE),
        ComputeBudget {
            work_units: u64::MAX,
            memory_bytes: 1,
        },
        IoBudget {
            processes: 1,
            files: 1,
            output_bytes: 1,
        },
        NetworkBudget {
            requests: 1,
            bytes: 1,
        },
        1,
    );
    let mut overflow = AdmissionMachine::try_new(online_context_with(
        max_budgets,
        TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        false,
    ))
    .expect("machine");
    admit(&mut overflow);
    overflow
        .apply(AdmissionTransition::Charge {
            charge: BudgetCharge::Work(u64::MAX),
            at_tick: 3,
        })
        .expect("maximum exact charge");
    assert_rejected_without_mutation(
        &mut overflow,
        AdmissionTransition::Charge {
            charge: BudgetCharge::Work(1),
            at_tick: 3,
        },
        AdmissionRule::BudgetOverflow,
    );
}

#[test]
fn g0_deadline_is_inclusive_and_uses_only_the_bound_clock_tick() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    machine
        .apply(AdmissionTransition::Preflight {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            at_tick: DEADLINE,
        })
        .expect("deadline tick is admitted inclusively");
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::StabilityRecheck {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            at_tick: DEADLINE + 1,
        },
        AdmissionRule::DeadlineExceeded,
    );
}

#[test]
fn g0_clock_ticks_never_regress_and_poll_work_is_bounded() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    machine
        .apply(AdmissionTransition::Preflight {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            at_tick: 10,
        })
        .expect("preflight");
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::PollCancellation {
            work_since_last_poll: 1,
            at_tick: 9,
        },
        AdmissionRule::ClockRegression,
    );
    machine
        .apply(AdmissionTransition::PollCancellation {
            work_since_last_poll: 8,
            at_tick: 10,
        })
        .expect("exact poll interval");
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::PollCancellation {
            work_since_last_poll: 9,
            at_tick: 11,
        },
        AdmissionRule::PollIntervalExceeded,
    );
    assert_eq!(machine.last_observed_tick(), Some(10));
}

#[test]
fn g0_transition_table_is_total_and_exact_phase_guards_do_not_mutate() {
    let states = [
        StateKind::Diagnostic,
        StateKind::Unanchored,
        StateKind::Admitted,
        StateKind::Refused,
        StateKind::Cancelled,
        StateKind::Indeterminate,
    ];
    let actions = [
        TransitionKind::Preflight,
        TransitionKind::StabilityRecheck,
        TransitionKind::Charge,
        TransitionKind::Refuse,
        TransitionKind::RequestCancellation,
        TransitionKind::Drain,
        TransitionKind::FinalizeTerminal,
        TransitionKind::DeclareIndeterminate,
        TransitionKind::BeginPublication,
        TransitionKind::PublicationRecheck,
        TransitionKind::FinalizePublication,
        TransitionKind::Retry,
        TransitionKind::PollCancellation,
        TransitionKind::AuthorizePublication,
        TransitionKind::PublicationFailed,
    ];
    let mut decisions = 0usize;
    for state in states {
        for action in actions {
            let _decision = transition_kind_may_apply(state, action);
            decisions += 1;
        }
    }
    assert_eq!(decisions, 90);
    assert!(!transition_kind_may_apply(
        StateKind::Indeterminate,
        TransitionKind::Retry
    ));
    assert!(!transition_kind_may_apply(
        StateKind::Refused,
        TransitionKind::FinalizePublication
    ));

    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    for transition in [
        AdmissionTransition::StabilityRecheck {
            snapshot: id(50),
            anchor: AnchorObservation::Unavailable,
            at_tick: 1,
        },
        AdmissionTransition::BeginPublication {
            quiescence: id(70),
            at_tick: 1,
        },
        AdmissionTransition::FinalizeTerminal { receipt: id(80) },
        AdmissionTransition::FinalizePublication,
    ] {
        assert_rejected_without_mutation(
            &mut machine,
            transition,
            AdmissionRule::IllegalTransition,
        );
    }
}

#[test]
fn g4_cancellation_requires_request_drain_finalize_and_retry_preserves_spend() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut machine);
    machine
        .apply(AdmissionTransition::Charge {
            charge: BudgetCharge::Work(6),
            at_tick: 3,
        })
        .expect("work reservation");
    machine
        .apply(AdmissionTransition::RequestCancellation {
            cause: CancellationCause::Requested,
            obligations: DrainObligations {
                processes: 2,
                files: 1,
                outputs: 1,
            },
            observation: id(70),
        })
        .expect("cancel request");

    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::RequestCancellation {
            cause: CancellationCause::ParentScope,
            obligations: DrainObligations::default(),
            observation: id(71),
        },
        AdmissionRule::IllegalTransition,
    );
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::FinalizeTerminal { receipt: id(80) },
        AdmissionRule::DrainIncomplete,
    );
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::Drain {
            completed: DrainObligations {
                processes: 3,
                files: 0,
                outputs: 0,
            },
            observation: id(72),
        },
        AdmissionRule::DrainOverrun,
    );
    machine
        .apply(AdmissionTransition::Drain {
            completed: DrainObligations {
                processes: 1,
                files: 1,
                outputs: 0,
            },
            observation: id(73),
        })
        .expect("partial drain");
    machine
        .apply(AdmissionTransition::Drain {
            completed: DrainObligations {
                processes: 1,
                files: 0,
                outputs: 1,
            },
            observation: id(74),
        })
        .expect("complete drain");
    machine
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(80) })
        .expect("cancel finalization");
    machine
        .apply(AdmissionTransition::Retry {
            next_cx: CxBinding::try_new(id(5), id(6), id(4), 8).expect("fresh retry Cx"),
            at_tick: 4,
        })
        .expect("retry");
    machine
        .apply(AdmissionTransition::Refuse {
            rule: AdmissionRule::BudgetExceeded,
        })
        .expect("definitive pre-effect retry refusal");
    machine
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(81) })
        .expect("refusal finalization");
    assert!(matches!(
        machine.state(),
        AdmissionState::Refused {
            rule: AdmissionRule::BudgetExceeded,
            phase: TerminalPhase::Finalized,
        }
    ));

    let stale_cx = machine.current_cx();
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::Retry {
            next_cx: stale_cx,
            at_tick: 4,
        },
        AdmissionRule::RetryCxNotFresh,
    );
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::Retry {
            next_cx: CxBinding::try_new(id(6), id(5), id(4), 8)
                .expect("role-swapped stale retry Cx"),
            at_tick: 4,
        },
        AdmissionRule::RetryCxNotFresh,
    );
    machine
        .apply(AdmissionTransition::Retry {
            next_cx: CxBinding::try_new(id(7), id(8), id(4), 8).expect("fresh retry Cx"),
            at_tick: 4,
        })
        .expect("bounded retry");
    assert_eq!(machine.attempt(), 2);
    assert_eq!(machine.consumption().retries, 2);
    assert_eq!(machine.consumption().work_units, 6);
    assert_eq!(machine.context().request_identity(), id(1));
    assert_eq!(machine.current_cx().cx(), id(7));
    assert_eq!(machine.current_cx().cancellation(), id(8));
    assert_eq!(
        machine.state(),
        AdmissionState::Diagnostic(DiagnosticPhase::Created)
    );
}

#[test]
fn g4_indeterminate_never_reuses_the_old_authority() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut machine);
    machine
        .apply(AdmissionTransition::DeclareIndeterminate {
            rule: AdmissionRule::DrainIncomplete,
        })
        .expect("uncertainty declaration");
    machine
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(80) })
        .expect("indeterminate record finalization");
    assert!(matches!(
        machine.state(),
        AdmissionState::Indeterminate {
            phase: TerminalPhase::Finalized,
            ..
        }
    ));
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::Retry {
            next_cx: CxBinding::try_new(id(5), id(6), id(4), 8).expect("fresh retry Cx"),
            at_tick: 4,
        },
        AdmissionRule::IndeterminateRetryForbidden,
    );
}

#[test]
fn g4_cleanup_uncertainty_escapes_cancellation_only_as_indeterminate() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut machine);
    machine
        .apply(AdmissionTransition::RequestCancellation {
            cause: CancellationCause::InjectedFault,
            obligations: DrainObligations {
                processes: 1,
                files: 1,
                outputs: 1,
            },
            observation: id(70),
        })
        .expect("cancel request");
    machine
        .apply(AdmissionTransition::DeclareIndeterminate {
            rule: AdmissionRule::DrainIncomplete,
        })
        .expect("failed cleanup remains uncertain");
    machine
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(80) })
        .expect("finalize uncertainty record");
    assert_eq!(
        machine.state(),
        AdmissionState::Indeterminate {
            rule: AdmissionRule::DrainIncomplete,
            phase: TerminalPhase::Finalized,
        }
    );
}

#[test]
fn g4_post_work_instability_and_commit_failure_never_become_refused() {
    let mut unstable = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut unstable);
    unstable
        .apply(AdmissionTransition::Charge {
            charge: BudgetCharge::Work(1),
            at_tick: 3,
        })
        .expect("charged work");
    unstable
        .apply(AdmissionTransition::BeginPublication {
            quiescence: id(70),
            at_tick: 4,
        })
        .expect("publication preparation");
    unstable
        .apply(AdmissionTransition::PublicationRecheck {
            snapshot: id(51),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            fence: id(71),
            at_tick: 5,
        })
        .expect("changed snapshot is a legal indeterminate observation");
    assert_eq!(
        unstable.state(),
        AdmissionState::Indeterminate {
            rule: AdmissionRule::PublicationStabilityChanged,
            phase: TerminalPhase::Pending,
        }
    );

    let mut failed_commit = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut failed_commit);
    failed_commit
        .apply(AdmissionTransition::BeginPublication {
            quiescence: id(70),
            at_tick: 3,
        })
        .expect("publication preparation");
    failed_commit
        .apply(AdmissionTransition::PublicationRecheck {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            fence: id(71),
            at_tick: 4,
        })
        .expect("publication fence");
    failed_commit
        .apply(AdmissionTransition::AuthorizePublication {
            receipt: id(60),
            at_tick: 5,
        })
        .expect("single-use authorization");
    failed_commit
        .apply(AdmissionTransition::PublicationFailed {
            rule: AdmissionRule::PublicationStabilityChanged,
        })
        .expect("failed external commit is indeterminate");
    assert!(matches!(
        failed_commit.state(),
        AdmissionState::Indeterminate {
            phase: TerminalPhase::Pending,
            ..
        }
    ));
    let failed_bytes = failed_commit
        .encode_canonical()
        .expect("failed publication receipt");
    assert_eq!(
        AdmissionMachine::decode_recorded(&failed_bytes)
            .expect("failed publication replay")
            .state(),
        failed_commit.state()
    );
}

#[test]
fn g0_event_capacity_always_retains_a_terminal_path() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    for _ in 0..MAX_ADMISSION_EVENTS - 2 {
        machine
            .apply(AdmissionTransition::Charge {
                charge: BudgetCharge::Work(0),
                at_tick: 1,
            })
            .expect("headroom-preserving event");
    }
    assert_rejected_without_mutation(
        &mut machine,
        AdmissionTransition::Charge {
            charge: BudgetCharge::Work(0),
            at_tick: 1,
        },
        AdmissionRule::TerminalHeadroomRequired,
    );
    machine
        .apply(AdmissionTransition::Refuse {
            rule: AdmissionRule::BudgetExceeded,
        })
        .expect("terminal fact uses reserved headroom");
    machine
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(80) })
        .expect("terminal finalization uses final slot");
    assert_eq!(machine.events().len(), MAX_ADMISSION_EVENTS);
    assert!(machine.state().attempt_is_finalized());
}

#[test]
fn g0_unanchored_and_mismatched_anchors_never_admit_work() {
    let mut unanchored = AdmissionMachine::try_new(online_context_with(
        online_budgets(),
        TrustAnchorState::Unanchored,
        false,
    ))
    .expect("unanchored machine");
    unanchored
        .apply(AdmissionTransition::Preflight {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            at_tick: 1,
        })
        .expect("unanchored observation is a terminal no-authority fact");
    assert_eq!(
        unanchored.state(),
        AdmissionState::Unanchored(TerminalPhase::Pending)
    );
    assert_rejected_without_mutation(
        &mut unanchored,
        AdmissionTransition::Charge {
            charge: BudgetCharge::Work(1),
            at_tick: 2,
        },
        AdmissionRule::IllegalTransition,
    );

    let mut mismatch = AdmissionMachine::try_new(online_context()).expect("machine");
    mismatch
        .apply(AdmissionTransition::Preflight {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(41),
                generation: 7,
            },
            at_tick: 1,
        })
        .expect("legal observation deterministically refuses");
    assert_eq!(
        mismatch.state(),
        AdmissionState::Refused {
            rule: AdmissionRule::TrustAnchorMismatch,
            phase: TerminalPhase::Pending,
        }
    );
}

#[test]
fn g0_offline_authority_is_independent_of_numeric_network_budget() {
    let mut offline = AdmissionMachine::try_new(offline_verify_context()).expect("machine");
    assert_rejected_without_mutation(
        &mut offline,
        AdmissionTransition::Charge {
            charge: BudgetCharge::Network {
                requests: 0,
                bytes: 0,
            },
            at_tick: 1,
        },
        AdmissionRule::NetworkAuthorityDenied,
    );

    let paths = [
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::ConstellationLock, id(11)),
    ];
    let executables = [ExecutableCapability::new(ExecutableSlot::Git, id(20))];
    let error = AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::VerifyOnly,
        fetch: FetchAuthority::PinnedTransport { capability: id(30) },
        publication: PublicationAuthority::Prohibited,
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: online_budgets(),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &paths,
        executable_capabilities: &executables,
    })
    .expect_err("verify-only cannot inherit fetch authority");
    assert_eq!(error.rule(), AdmissionRule::NetworkForbiddenForCommand);
}

#[test]
fn g5_canonical_receipt_replays_without_recreating_authority() {
    let mut machine = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut machine);
    machine
        .apply(AdmissionTransition::RequestCancellation {
            cause: CancellationCause::InjectedFault,
            obligations: DrainObligations::default(),
            observation: id(70),
        })
        .expect("cancel request");
    machine
        .apply(AdmissionTransition::Drain {
            completed: DrainObligations::default(),
            observation: id(71),
        })
        .expect("explicit zero-obligation drain");
    machine
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(80) })
        .expect("cancel finalization");
    let bytes = machine.encode_canonical().expect("canonical receipt");
    let recorded = AdmissionMachine::decode_recorded(&bytes).expect("inert replay");
    assert_eq!(recorded.canonical_bytes(), bytes);
    assert_eq!(recorded.request_identity(), id(1));
    assert_eq!(recorded.command(), CommandClass::Bootstrap);
    assert_eq!(recorded.context(), machine.context());
    assert_eq!(recorded.state(), machine.state());
    assert_eq!(recorded.consumption(), machine.consumption());
    assert_eq!(recorded.attempt(), machine.attempt());
    assert_eq!(recorded.current_cx(), machine.current_cx());
    assert_eq!(recorded.last_observed_tick(), machine.last_observed_tick());
    assert_eq!(recorded.history(), machine.history());
    assert_eq!(recorded.events(), machine.events());

    let left = AdmissionMachine::try_new(online_context_with(
        online_budgets(),
        TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        false,
    ))
    .expect("left");
    let right = AdmissionMachine::try_new(online_context_with(
        online_budgets(),
        TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        true,
    ))
    .expect("right");
    assert_eq!(
        left.encode_canonical().expect("left bytes"),
        right.encode_canonical().expect("right bytes"),
        "capability insertion order cannot change request identity bytes"
    );

    let shell_context = |reverse: bool| {
        let paths = [
            PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
            PathCapability::new(PathSlot::ConstellationLock, id(11)),
        ];
        let mut executables = vec![
            ExecutableCapability::new(ExecutableSlot::Git, id(20)),
            ExecutableCapability::new(ExecutableSlot::Shell, id(21)),
        ];
        if reverse {
            executables.reverse();
        }
        AdmissionContext::try_new(AdmissionContextSpec {
            request_identity: id(1),
            command: CommandClass::ShellVerifyOnly,
            fetch: FetchAuthority::Offline,
            publication: PublicationAuthority::Prohibited,
            cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
            budgets: AdmissionBudgets::new(
                DeadlineBudget::new(id(4), DEADLINE),
                ComputeBudget {
                    work_units: 1,
                    memory_bytes: 1,
                },
                IoBudget {
                    processes: 1,
                    files: 1,
                    output_bytes: 1,
                },
                NetworkBudget {
                    requests: 0,
                    bytes: 0,
                },
                0,
            ),
            trust_anchor: TrustAnchorState::Anchored {
                identity: id(40),
                generation: 7,
            },
            path_capabilities: &paths,
            executable_capabilities: &executables,
        })
        .expect("shell context")
    };
    let shell_left = AdmissionMachine::try_new(shell_context(false)).expect("shell left");
    let shell_right = AdmissionMachine::try_new(shell_context(true)).expect("shell right");
    assert_eq!(
        shell_left.encode_canonical().expect("shell left bytes"),
        shell_right.encode_canonical().expect("shell right bytes"),
        "executable insertion order cannot change canonical bytes"
    );
}

#[test]
fn g5_codec_round_trips_every_payload_family() {
    let mut published = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut published);
    published
        .apply(AdmissionTransition::PollCancellation {
            work_since_last_poll: 8,
            at_tick: 3,
        })
        .expect("poll");
    for charge in [
        BudgetCharge::Work(1),
        BudgetCharge::Memory(2),
        BudgetCharge::Processes(1),
        BudgetCharge::Files(1),
        BudgetCharge::Output(3),
        BudgetCharge::Network {
            requests: 1,
            bytes: 4,
        },
    ] {
        published
            .apply(AdmissionTransition::Charge { charge, at_tick: 3 })
            .expect("charge payload");
    }
    published
        .apply(AdmissionTransition::BeginPublication {
            quiescence: id(70),
            at_tick: 4,
        })
        .expect("prepare");
    published
        .apply(AdmissionTransition::PublicationRecheck {
            snapshot: id(50),
            anchor: AnchorObservation::Observed {
                identity: id(40),
                generation: 7,
            },
            fence: id(71),
            at_tick: 5,
        })
        .expect("fence");
    published
        .apply(AdmissionTransition::AuthorizePublication {
            receipt: id(60),
            at_tick: 6,
        })
        .expect("authorize");
    published
        .apply(AdmissionTransition::FinalizePublication)
        .expect("finalize");
    let published_bytes = published.encode_canonical().expect("published bytes");
    assert_eq!(
        AdmissionMachine::decode_recorded(&published_bytes)
            .expect("published replay")
            .state(),
        published.state()
    );

    let mut unavailable = AdmissionMachine::try_new(online_context()).expect("machine");
    unavailable
        .apply(AdmissionTransition::Preflight {
            snapshot: id(50),
            anchor: AnchorObservation::Unavailable,
            at_tick: 1,
        })
        .expect("unavailable anchor refusal");
    unavailable
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(80) })
        .expect("refusal finalization");
    let unavailable_bytes = unavailable.encode_canonical().expect("unavailable bytes");
    assert_eq!(
        AdmissionMachine::decode_recorded(&unavailable_bytes)
            .expect("unavailable replay")
            .state(),
        unavailable.state()
    );

    let mut retried = AdmissionMachine::try_new(online_context()).expect("machine");
    retried
        .apply(AdmissionTransition::Refuse {
            rule: AdmissionRule::BudgetExceeded,
        })
        .expect("explicit refusal");
    retried
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(80) })
        .expect("refusal finalization");
    let retry_cx = CxBinding::try_new(id(5), id(6), id(4), 8).expect("retry Cx");
    retried
        .apply(AdmissionTransition::Retry {
            next_cx: retry_cx,
            at_tick: 1,
        })
        .expect("retry");
    let retried_bytes = retried.encode_canonical().expect("retried bytes");
    let retried_record = AdmissionMachine::decode_recorded(&retried_bytes).expect("retried replay");
    assert_eq!(retried_record.context().cx(), online_context().cx());
    assert_eq!(retried_record.current_cx(), retry_cx);
    assert_eq!(retried_record.attempt(), 1);
    assert_eq!(retried_record.history(), retried.history());

    let mut declared = AdmissionMachine::try_new(online_context()).expect("machine");
    admit(&mut declared);
    declared
        .apply(AdmissionTransition::DeclareIndeterminate {
            rule: AdmissionRule::DrainIncomplete,
        })
        .expect("uncertainty declaration");
    declared
        .apply(AdmissionTransition::FinalizeTerminal { receipt: id(81) })
        .expect("uncertainty finalization");
    let declared_bytes = declared.encode_canonical().expect("declared bytes");
    assert_eq!(
        AdmissionMachine::decode_recorded(&declared_bytes)
            .expect("declared replay")
            .state(),
        declared.state()
    );

    for cut in 0..published_bytes.len() {
        assert!(
            AdmissionMachine::decode_recorded(&published_bytes[..cut]).is_err(),
            "transition-rich truncation at byte {cut} must refuse"
        );
    }
}

#[test]
fn g5_decoder_rejects_schema_truncation_tags_history_and_trailing_bytes() {
    let machine = AdmissionMachine::try_new(online_context()).expect("machine");
    let bytes = machine.encode_canonical().expect("canonical bytes");

    let mut bad_magic = bytes.clone();
    bad_magic[0] ^= 0xff;
    assert_eq!(
        AdmissionMachine::decode_recorded(&bad_magic)
            .expect_err("unknown magic")
            .rule(),
        AdmissionRule::UnknownMagic
    );

    let mut bad_version = bytes.clone();
    bad_version[8] = bad_version[8].wrapping_add(1);
    assert_eq!(
        AdmissionMachine::decode_recorded(&bad_version)
            .expect_err("unknown version")
            .rule(),
        AdmissionRule::UnknownSchema
    );

    let mut bad_schema = bytes.clone();
    let first_schema_byte = 8 + 2 + 2;
    bad_schema[first_schema_byte] ^= 1;
    assert_eq!(
        AdmissionMachine::decode_recorded(&bad_schema)
            .expect_err("unknown schema")
            .rule(),
        AdmissionRule::UnknownSchema
    );

    let domain_length_at = 8 + 2 + 2 + admission::ADMISSION_SCHEMA.len();
    let first_domain_byte = domain_length_at + 2;
    let mut bad_domain = bytes.clone();
    bad_domain[first_domain_byte] ^= 1;
    assert_eq!(
        AdmissionMachine::decode_recorded(&bad_domain)
            .expect_err("unknown domain")
            .rule(),
        AdmissionRule::UnknownSchema
    );

    let mut bad_domain_length = bytes.clone();
    bad_domain_length[domain_length_at..domain_length_at + 2]
        .copy_from_slice(&u16::MAX.to_le_bytes());
    assert_eq!(
        AdmissionMachine::decode_recorded(&bad_domain_length)
            .expect_err("oversized declared domain")
            .rule(),
        AdmissionRule::TruncatedEncoding
    );

    let request_at = first_domain_byte + admission::ADMISSION_DOMAIN.len();
    let mut zero_request = bytes.clone();
    zero_request[request_at..request_at + 32].fill(0);
    assert_eq!(
        AdmissionMachine::decode_recorded(&zero_request)
            .expect_err("zero request identity")
            .rule(),
        AdmissionRule::ZeroIdentity
    );

    let mut unknown_command = bytes.clone();
    unknown_command[request_at + 32] = 0xff;
    assert_eq!(
        AdmissionMachine::decode_recorded(&unknown_command)
            .expect_err("unknown command tag")
            .rule(),
        AdmissionRule::UnknownTag
    );

    let oversized = vec![0; MAX_ADMISSION_BYTES + 1];
    assert_eq!(
        AdmissionMachine::decode_recorded(&oversized)
            .expect_err("oversized envelope")
            .rule(),
        AdmissionRule::EncodedInputTooLarge
    );

    for cut in 0..bytes.len() {
        assert!(
            AdmissionMachine::decode_recorded(&bytes[..cut]).is_err(),
            "truncation at byte {cut} must refuse"
        );
    }

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert_eq!(
        AdmissionMachine::decode_recorded(&trailing)
            .expect_err("trailing byte")
            .rule(),
        AdmissionRule::TrailingBytes
    );

    let mut unknown_transition = bytes.clone();
    let count_at = unknown_transition.len() - 2;
    unknown_transition[count_at..].copy_from_slice(&1u16.to_le_bytes());
    unknown_transition.push(0xff);
    assert_eq!(
        AdmissionMachine::decode_recorded(&unknown_transition)
            .expect_err("unknown transition")
            .rule(),
        AdmissionRule::UnknownTag
    );

    let mut excessive_events = bytes.clone();
    let excessive_event_count =
        u16::try_from(MAX_ADMISSION_EVENTS + 1).expect("event fixture count fits canonical u16");
    excessive_events[count_at..].copy_from_slice(&excessive_event_count.to_le_bytes());
    assert_eq!(
        AdmissionMachine::decode_recorded(&excessive_events)
            .expect_err("excessive event count")
            .rule(),
        AdmissionRule::EventLimitExceeded
    );

    let mut impossible_history = bytes;
    let count_at = impossible_history.len() - 2;
    impossible_history[count_at..].copy_from_slice(&1u16.to_le_bytes());
    impossible_history.push(TransitionKind::FinalizeTerminal as u8);
    impossible_history.extend_from_slice(id(80).as_bytes());
    assert_eq!(
        AdmissionMachine::decode_recorded(&impossible_history)
            .expect_err("impossible history")
            .rule(),
        AdmissionRule::ImpossibleHistory
    );

    let canonical = AdmissionMachine::try_new(online_context())
        .expect("machine")
        .encode_canonical()
        .expect("canonical context");
    let mut workspace_row = vec![PathSlot::WorkspaceRoot as u8];
    workspace_row.extend_from_slice(id(10).as_bytes());
    let mut lock_row = vec![PathSlot::ConstellationLock as u8];
    lock_row.extend_from_slice(id(11).as_bytes());
    let workspace_at = canonical
        .windows(workspace_row.len())
        .position(|window| window == workspace_row.as_slice())
        .expect("workspace capability row");
    let lock_at = canonical
        .windows(lock_row.len())
        .position(|window| window == lock_row.as_slice())
        .expect("lock capability row");
    let mut reordered = canonical;
    reordered[workspace_at..workspace_at + workspace_row.len()].copy_from_slice(&lock_row);
    reordered[lock_at..lock_at + lock_row.len()].copy_from_slice(&workspace_row);
    assert_eq!(
        AdmissionMachine::decode_recorded(&reordered)
            .expect_err("noncanonical capability ordering")
            .rule(),
        AdmissionRule::NonCanonicalEncoding
    );
}

#[test]
fn g0_context_rejects_duplicate_missing_and_mismatched_capabilities() {
    let duplicate_paths = [
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::WorkspaceRoot, id(11)),
        PathCapability::new(PathSlot::ConstellationLock, id(12)),
        PathCapability::new(PathSlot::DestinationRoot, id(13)),
        PathCapability::new(PathSlot::PublicationTarget, id(14)),
    ];
    let executables = [ExecutableCapability::new(ExecutableSlot::Git, id(20))];
    let error = AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::Bootstrap,
        fetch: FetchAuthority::PinnedTransport { capability: id(30) },
        publication: PublicationAuthority::BootstrapReceipt { capability: id(14) },
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: online_budgets(),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &duplicate_paths,
        executable_capabilities: &executables,
    })
    .expect_err("duplicate path slot");
    assert_eq!(error.rule(), AdmissionRule::DuplicatePathCapability);

    let missing_paths = [
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::ConstellationLock, id(11)),
        PathCapability::new(PathSlot::DestinationRoot, id(12)),
    ];
    let error = AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::Bootstrap,
        fetch: FetchAuthority::PinnedTransport { capability: id(30) },
        publication: PublicationAuthority::BootstrapReceipt { capability: id(13) },
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: online_budgets(),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &missing_paths,
        executable_capabilities: &executables,
    })
    .expect_err("missing publication target");
    assert_eq!(error.rule(), AdmissionRule::MissingPathCapability);

    let paths = [
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::ConstellationLock, id(11)),
        PathCapability::new(PathSlot::DestinationRoot, id(12)),
        PathCapability::new(PathSlot::PublicationTarget, id(13)),
    ];
    let error = AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::Bootstrap,
        fetch: FetchAuthority::PinnedTransport { capability: id(30) },
        publication: PublicationAuthority::BootstrapReceipt { capability: id(99) },
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: online_budgets(),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &paths,
        executable_capabilities: &executables,
    })
    .expect_err("publication capability mismatch");
    assert_eq!(error.rule(), AdmissionRule::PublicationCapabilityMismatch);
}

#[test]
fn g0_command_classes_reject_unexpected_authority_slots() {
    let paths = [
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::ConstellationLock, id(11)),
        PathCapability::new(PathSlot::PublicationTarget, id(12)),
    ];
    let executables = [ExecutableCapability::new(ExecutableSlot::Git, id(20))];
    let error = AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::VerifyOnly,
        fetch: FetchAuthority::Offline,
        publication: PublicationAuthority::Prohibited,
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: AdmissionBudgets::new(
            DeadlineBudget::new(id(4), DEADLINE),
            ComputeBudget {
                work_units: 1,
                memory_bytes: 1,
            },
            IoBudget {
                processes: 1,
                files: 1,
                output_bytes: 1,
            },
            NetworkBudget {
                requests: 0,
                bytes: 0,
            },
            0,
        ),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &paths,
        executable_capabilities: &executables,
    })
    .expect_err("verify-only cannot carry a publication path capability");
    assert_eq!(error.rule(), AdmissionRule::UnexpectedPathCapability);

    let paths = [
        PathCapability::new(PathSlot::WorkspaceRoot, id(10)),
        PathCapability::new(PathSlot::ConstellationLock, id(11)),
    ];
    let executables = [
        ExecutableCapability::new(ExecutableSlot::Git, id(20)),
        ExecutableCapability::new(ExecutableSlot::Rch, id(21)),
    ];
    let error = AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::VerifyOnly,
        fetch: FetchAuthority::Offline,
        publication: PublicationAuthority::Prohibited,
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: AdmissionBudgets::new(
            DeadlineBudget::new(id(4), DEADLINE),
            ComputeBudget {
                work_units: 1,
                memory_bytes: 1,
            },
            IoBudget {
                processes: 1,
                files: 1,
                output_bytes: 1,
            },
            NetworkBudget {
                requests: 0,
                bytes: 0,
            },
            0,
        ),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &paths,
        executable_capabilities: &executables,
    })
    .expect_err("verify-only cannot inherit an RCH executable");
    assert_eq!(error.rule(), AdmissionRule::UnexpectedExecutableCapability);

    let error = AdmissionContext::try_new(AdmissionContextSpec {
        request_identity: id(1),
        command: CommandClass::ShellVerifyOnly,
        fetch: FetchAuthority::Offline,
        publication: PublicationAuthority::Prohibited,
        cx: CxBinding::try_new(id(2), id(3), id(4), 8).expect("fixture Cx"),
        budgets: AdmissionBudgets::new(
            DeadlineBudget::new(id(4), DEADLINE),
            ComputeBudget {
                work_units: 1,
                memory_bytes: 1,
            },
            IoBudget {
                processes: 1,
                files: 1,
                output_bytes: 1,
            },
            NetworkBudget {
                requests: 0,
                bytes: 0,
            },
            0,
        ),
        trust_anchor: TrustAnchorState::Anchored {
            identity: id(40),
            generation: 7,
        },
        path_capabilities: &paths,
        executable_capabilities: &[ExecutableCapability::new(ExecutableSlot::Git, id(20))],
    })
    .expect_err("shell adapter requires an explicit shell capability");
    assert_eq!(error.rule(), AdmissionRule::MissingExecutableCapability);
}
