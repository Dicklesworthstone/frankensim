//! G0/G3 integration battery for expansion-program risk session snapshots.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_exec::{CancelGate, DrainTracker, RunId};
use fs_govern::program_risks::{
    AssessmentStatus, ProgramRiskId, ProgramRiskObservation, TriggerComparison, program_risk,
};
use fs_session::{
    CapabilityToken, DurableGovernorNonce, Governor, ProgramRiskReportDisposition, SessionError,
    SessionId,
};

fn token(session: SessionId, scope: &str) -> CapabilityToken {
    CapabilityToken {
        session,
        ops: vec!["governance.*".to_string()],
        core_s: 60.0,
        mem_bytes: 1024 * 1024,
        wall_s: 60.0,
        cores: 1,
        ledger_scope: scope.to_string(),
    }
}

fn all_clear_observations() -> Vec<ProgramRiskObservation<'static>> {
    ProgramRiskId::ALL
        .into_iter()
        .map(|id| {
            let risk = program_risk(id);
            let value = match risk.trigger.comparison {
                TriggerComparison::GreaterThanOrEqual => risk.trigger.threshold - 1.0,
                TriggerComparison::LessThan => risk.trigger.threshold,
            };
            ProgramRiskObservation {
                id,
                value,
                unit: risk.trigger.unit,
                samples: risk.trigger.min_samples,
            }
        })
        .collect()
}

fn observations_with_pr001_trigger() -> Vec<ProgramRiskObservation<'static>> {
    let mut observations = all_clear_observations();
    observations[0].value = program_risk(ProgramRiskId::Pr001).trigger.threshold;
    observations
}

fn durable_counts(ledger: &fs_ledger::Ledger) -> [u64; 13] {
    [
        ledger.table_count("artifacts").expect("artifact count"),
        ledger
            .table_count("artifact_chunks")
            .expect("artifact-chunk count"),
        ledger.table_count("ops").expect("op count"),
        ledger.table_count("edges").expect("edge count"),
        ledger
            .table_count("artifact_output_seals")
            .expect("artifact-output-seal count"),
        ledger
            .table_count("op_artifact_edge_seals")
            .expect("op-artifact-edge-seal count"),
        ledger.table_count("session_claims").expect("claim count"),
        ledger
            .table_count("session_claim_discovery")
            .expect("claim-discovery count"),
        ledger
            .table_count("session_terminals")
            .expect("terminal count"),
        ledger
            .table_count("session_terminal_events")
            .expect("terminal-event count"),
        ledger
            .table_count("session_flush_batches")
            .expect("batch count"),
        ledger
            .table_count("session_flush_batch_members")
            .expect("batch-member count"),
        ledger.table_count("events").expect("event count"),
    ]
}

fn durable_ledger_path(case: &str) -> String {
    static NEXT: AtomicU64 = AtomicU64::new(1);
    let ordinal = NEXT.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!(
            "fs-session-program-risk-{}-{ordinal}-{case}.ledger",
            std::process::id()
        ))
        .to_string_lossy()
        .into_owned()
}

fn solver_checkpoint(
    ledger: &fs_ledger::Ledger,
    request: fs_session::PauseRequestId,
    label: &str,
) -> fs_ledger::session_registry::SolverCheckpointReceipt {
    let authority = request.checkpoint_authority();
    let mut run_bytes = [0_u8; 8];
    run_bytes.copy_from_slice(&authority.as_bytes()[..8]);
    let run = RunId(u64::from_le_bytes(run_bytes));
    let gate = CancelGate::new_clock_free();
    let tracker = DrainTracker::new(run, &gate);
    let worker = tracker.register_worker().expect("fixture worker");
    gate.request();
    drop(worker);
    let report = tracker.finalize().expect("fixture run drained");
    let snapshot = fs_exec::solver::envelope::seal(0x4653_434b_5054, 1, run.0, label.as_bytes());
    let artifact = ledger
        .put_artifact(
            fs_ledger::session_registry::SOLVER_STATE_ARTIFACT_KIND,
            &snapshot,
            None,
        )
        .expect("fixture solver-state artifact");
    ledger
        .attest_solver_checkpoint(fs_ledger::session_registry::SolverCheckpointClaim {
            session: request.session().0,
            pause_authority: authority,
            gate_generation: request.gate_generation(),
            solver_state_artifact: artifact.hash,
            drain_report: &report,
        })
        .expect("fixture checkpoint receipt")
}

#[test]
#[allow(clippy::too_many_lines)] // One complete persistence/replay/GC proof over the same report.
fn g0_tripped_risk_is_automatically_artifacted_rooted_and_surfaced_once() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("program-risk ledger");
    let governor = Governor::new();
    let session = SessionId(70_001);
    let open_id = governor
        .session_open_id(session, "program-risk-open")
        .expect("open authority");
    let open = governor
        .open_session(open_id, token(session, "program-risk"))
        .expect("session open");
    let open_flush = governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("persist open prerequisite");
    assert_eq!(open_flush.committed_terminals, 1);
    assert_eq!(open_flush.appended_rows, 1);
    assert!(!open_flush.remaining_dirty);

    let observations = observations_with_pr001_trigger();
    let first = governor
        .write_program_risk_session_end_report(&ledger, &open, 42, &observations)
        .expect("publish session-end snapshot");
    assert_eq!(first.disposition, ProgramRiskReportDisposition::Committed);
    assert_eq!(first.receipt.alert_count(), 1);
    assert_eq!(
        first.receipt.alerts(),
        vec![fs_session::ProgramRiskAlert {
            id: ProgramRiskId::Pr001,
            status: AssessmentStatus::Triggered,
        }]
    );
    assert_eq!(ledger.table_count("events").unwrap(), 2);
    let terminal = ledger
        .session_terminal(&first.receipt.report_id().content_hash())
        .expect("terminal query")
        .expect("program-risk terminal");
    assert_eq!(terminal.event_count, 1, "one owned alert-summary event");

    let register_info = ledger
        .artifact_info(&first.receipt.register_artifact())
        .expect("register metadata")
        .expect("register artifact");
    assert_eq!(
        register_info.kind,
        fs_session::PROGRAM_RISK_REGISTER_ARTIFACT_KIND
    );
    let report_info = ledger
        .artifact_info(&first.receipt.report_artifact())
        .expect("report metadata")
        .expect("report artifact");
    assert_eq!(
        report_info.kind,
        fs_session::PROGRAM_RISK_SESSION_REPORT_ARTIFACT_KIND
    );
    let report_json = String::from_utf8(
        ledger
            .get_artifact(&first.receipt.report_artifact())
            .expect("report read")
            .expect("report bytes"),
    )
    .expect("report utf8");
    assert!(report_json.contains("\"alert_count\":1"));
    assert!(report_json.contains("\"generation\":0"));
    assert!(report_json.contains("\"id\":\"PR-001\",\"status\":\"triggered\""));
    let producers = ledger
        .artifact_producer_ops_bounded(&first.receipt.report_artifact(), 1)
        .expect("bounded lineage query");
    assert!(!producers.truncated);
    assert_eq!(producers.op_ids, vec![first.receipt.lineage_op()]);
    assert_eq!(
        ledger
            .artifact_output_seal(&first.receipt.report_artifact())
            .expect("report output seal"),
        Some(first.receipt.lineage_op())
    );
    assert_eq!(
        ledger
            .op_artifact_edge_seal(first.receipt.lineage_op())
            .expect("lineage edge-set seal"),
        Some(2)
    );

    let counts = durable_counts(&ledger);
    let replay = governor
        .write_program_risk_session_end_report(&ledger, &open, 42, &observations)
        .expect("exact report replay");
    assert_eq!(replay.disposition, ProgramRiskReportDisposition::Replayed);
    assert_eq!(replay.receipt, first.receipt);
    assert_eq!(durable_counts(&ledger), counts, "exact replay adds no row");

    let altered = all_clear_observations();
    assert!(matches!(
        governor.write_program_risk_session_end_report(&ledger, &open, 42, &altered),
        Err(SessionError::MutationConflict {
            kind: "program-risk-report",
            ..
        })
    ));
    assert_eq!(
        durable_counts(&ledger),
        counts,
        "altered singleton retry is atomic"
    );

    let gc = ledger
        .gc_unreferenced_artifacts(false)
        .expect("lineage-aware GC");
    assert!(
        !gc.candidates
            .contains(&first.receipt.register_artifact().to_hex())
    );
    assert!(
        !gc.candidates
            .contains(&first.receipt.report_artifact().to_hex())
    );
    assert!(
        ledger
            .get_artifact(&first.receipt.register_artifact())
            .expect("register query after GC")
            .is_some()
    );
    assert!(
        ledger
            .get_artifact(&first.receipt.report_artifact())
            .expect("report query after GC")
            .is_some()
    );
    assert!(ledger.lint().expect("ledger lint").is_clean());
}

#[test]
fn g3_sealed_report_refuses_a_second_producer_before_replay() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("program-risk ledger");
    let governor = Governor::new();
    let session = SessionId(70_101);
    let open_id = governor
        .session_open_id(session, "program-risk-multiple-producers")
        .expect("open authority");
    let open = governor
        .open_session(open_id, token(session, "program-risk-multiple-producers"))
        .expect("session open");
    governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("persist open prerequisite");
    let observations = observations_with_pr001_trigger();
    let first = governor
        .write_program_risk_session_end_report(&ledger, &open, 42, &observations)
        .expect("publish canonical report");

    let explicits = fs_ledger::FiveExplicits {
        seed: b"foreign-producer",
        versions: "{}",
        budget: "{}",
        capability: "{}",
    };
    let foreign = ledger
        .begin_op(None, "{}", &explicits, 0)
        .expect("foreign producer op");
    let error = ledger
        .link(
            foreign,
            &first.receipt.report_artifact(),
            fs_ledger::EdgeRole::Out,
        )
        .expect_err("sealed report must reject a foreign output edge");
    assert!(matches!(
        error,
        fs_ledger::LedgerError::Invalid { field, problem }
            if field == "edge" && problem.contains("exclusive output-producer seal")
    ));
    ledger
        .finish_op(foreign, fs_ledger::OpOutcome::Ok, None, 0)
        .expect("finish foreign producer");

    let replay = governor
        .write_program_risk_session_end_report(&ledger, &open, 42, &observations)
        .expect("sealed exact report replay");
    assert_eq!(replay.disposition, ProgramRiskReportDisposition::Replayed);
    assert_eq!(replay.receipt, first.receipt);
}

#[test]
fn g3_sealed_lineage_refuses_a_third_edge_before_replay() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("program-risk ledger");
    let governor = Governor::new();
    let session = SessionId(70_102);
    let open_id = governor
        .session_open_id(session, "program-risk-extra-edge")
        .expect("open authority");
    let open = governor
        .open_session(open_id, token(session, "program-risk-extra-edge"))
        .expect("session open");
    governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("persist open prerequisite");
    let observations = observations_with_pr001_trigger();
    let first = governor
        .write_program_risk_session_end_report(&ledger, &open, 42, &observations)
        .expect("publish canonical report");

    let extra = ledger
        .put_artifact("program-risk-adversarial-extra", b"extra-edge", None)
        .expect("extra artifact");
    let error = ledger
        .link(
            first.receipt.lineage_op(),
            &extra.hash,
            fs_ledger::EdgeRole::In,
        )
        .expect_err("sealed lineage op must reject a third edge");
    assert!(matches!(
        error,
        fs_ledger::LedgerError::Invalid { field, problem }
            if field == "edge" && problem.contains("artifact-edge-set seal")
    ));

    let replay = governor
        .write_program_risk_session_end_report(&ledger, &open, 42, &observations)
        .expect("sealed exact report replay");
    assert_eq!(replay.disposition, ProgramRiskReportDisposition::Replayed);
    assert_eq!(replay.receipt, first.receipt);
}

#[test]
fn g0_report_refuses_unflushed_foreign_and_wrong_sink_contexts() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("owning ledger");
    let governor = Governor::new();
    let session = SessionId(70_002);
    let open_id = governor
        .session_open_id(session, "context-open")
        .expect("open authority");
    let open = governor
        .open_session(open_id, token(session, "context"))
        .expect("session open");
    let observations = all_clear_observations();

    assert!(matches!(
        governor.write_program_risk_session_end_report(&ledger, &open, 1, &observations),
        Err(SessionError::Persistence { .. })
    ));
    assert_eq!(ledger.table_count("artifacts").unwrap(), 0);

    governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("bind owning sink");
    let wrong_sink = fs_ledger::Ledger::open(":memory:").expect("wrong sink");
    assert!(matches!(
        governor.write_program_risk_session_end_report(&wrong_sink, &open, 1, &observations),
        Err(SessionError::LedgerScopeSinkMismatch { .. })
    ));
    assert_eq!(wrong_sink.table_count("artifacts").unwrap(), 0);

    let foreign = Governor::new();
    let foreign_session = SessionId(70_003);
    let foreign_id = foreign
        .session_open_id(foreign_session, "foreign-open")
        .expect("foreign open authority");
    let foreign_open = foreign
        .open_session(foreign_id, token(foreign_session, "foreign"))
        .expect("foreign open");
    assert!(matches!(
        governor.write_program_risk_session_end_report(&ledger, &foreign_open, 1, &observations),
        Err(SessionError::MutationAuthorityMismatch {
            kind: "program-risk-report",
            ..
        })
    ));
}

#[test]
fn g3_exact_replay_retains_the_original_generation_after_lifecycle_progress() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("generation ledger");
    let governor = Governor::new();
    let session = SessionId(70_005);
    let gate = Arc::new(CancelGate::new());
    let open_id = governor
        .session_open_id(session, "generation-open")
        .expect("open authority");
    let open = governor
        .open_session_gated(open_id, token(session, "generation"), Arc::clone(&gate))
        .expect("gated session open");
    governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("persist open prerequisite");
    let observations = all_clear_observations();
    let first = governor
        .write_program_risk_session_end_report(&ledger, &open, 51, &observations)
        .expect("generation-zero report");
    assert_eq!(first.receipt.generation(), 0);

    let action_id = governor
        .pressure_action_id(session, "advance-after-report")
        .expect("pressure authority");
    let pressure = governor
        .apply_memory_pressure(action_id, 3)
        .expect("request pause");
    let request_id = pressure
        .events()
        .iter()
        .find_map(|event| event.pause_request_id)
        .expect("pause request authority");
    let checkpoint = solver_checkpoint(&ledger, request_id, "post-report-checkpoint");
    let acknowledgement = governor
        .acknowledge_pause(request_id, &ledger, &checkpoint)
        .expect("acknowledge drained generation");
    governor
        .activate_resume(&acknowledgement)
        .expect("advance to the next execution generation");

    let counts = durable_counts(&ledger);
    let replay = governor
        .write_program_risk_session_end_report(&ledger, &open, 51, &observations)
        .expect("old report remains exactly replayable");
    assert_eq!(replay.disposition, ProgramRiskReportDisposition::Replayed);
    assert_eq!(replay.receipt, first.receipt);
    assert_eq!(replay.receipt.generation(), 0);
    assert_eq!(durable_counts(&ledger), counts);
}

#[test]
fn g3_durable_reopen_requires_typed_report_recovery_and_changes_no_rows() {
    let path = durable_ledger_path("recovery");
    let nonce = DurableGovernorNonce::from_bytes([0xA7; 32]);
    let session = SessionId(70_004);
    let capability = token(session, "durable-program-risk");
    let observations = observations_with_pr001_trigger();

    let ledger = fs_ledger::Ledger::open(&path).expect("durable program-risk ledger");
    let governor = Governor::new_durable(&ledger, nonce).expect("durable governor");
    let open_id = governor
        .session_open_id(session, "durable-program-risk-open")
        .expect("open authority");
    let open = governor
        .open_session(open_id, capability.clone())
        .expect("durable open");
    governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("persist durable open");
    let published = governor
        .write_program_risk_session_end_report(&ledger, &open, 77, &observations)
        .expect("persist durable program-risk report");
    let counts = durable_counts(&ledger);
    drop(governor);
    drop(ledger);

    let ledger = fs_ledger::Ledger::open(&path).expect("reopened program-risk ledger");
    let governor = Governor::new_durable(&ledger, nonce).expect("reopened governor");
    let recovered_open = governor
        .recover_open(&ledger, open_id, capability, None)
        .expect("recover open first");
    assert!(matches!(
        governor.write_program_risk_session_end_report(&ledger, &recovered_open, 77, &observations),
        Err(SessionError::DurableRecoveryIncomplete {
            remaining_claims: 1
        })
    ));

    let recovered = governor
        .recover_program_risk_report(&ledger, &recovered_open)
        .expect("typed report recovery");
    assert_eq!(recovered, published.receipt);
    assert_eq!(durable_counts(&ledger), counts);

    let replay = governor
        .write_program_risk_session_end_report(&ledger, &recovered_open, 77, &observations)
        .expect("post-recovery exact replay");
    assert_eq!(replay.disposition, ProgramRiskReportDisposition::Replayed);
    assert_eq!(replay.receipt, recovered);
    assert_eq!(
        durable_counts(&ledger),
        counts,
        "typed recovery and exact replay add no durable row"
    );
    assert!(ledger.lint().expect("reopened ledger lint").is_clean());
}

#[test]
fn g3_report_recovery_waits_for_its_pause_resume_generation() {
    let path = durable_ledger_path("generation-recovery-order");
    let nonce = DurableGovernorNonce::from_bytes([0xB8; 32]);
    let session = SessionId(70_006);
    let capability = token(session, "durable-program-risk-generation");
    let observations = all_clear_observations();

    let ledger = fs_ledger::Ledger::open(&path).expect("durable generation ledger");
    let governor = Governor::new_durable(&ledger, nonce).expect("durable governor");
    let initial_gate = Arc::new(CancelGate::new());
    let open_id = governor
        .session_open_id(session, "durable-program-risk-generation-open")
        .expect("open authority");
    let open = governor
        .open_session_gated(open_id, capability.clone(), Arc::clone(&initial_gate))
        .expect("durable gated open");
    governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("persist open prerequisite");

    let action_id = governor
        .pressure_action_id(session, "advance-before-program-risk-report")
        .expect("pressure authority");
    let pressure = governor
        .apply_memory_pressure(action_id, 3)
        .expect("request pause");
    let request_id = pressure
        .events()
        .iter()
        .find_map(|event| event.pause_request_id)
        .expect("pause request authority");
    let checkpoint = solver_checkpoint(&ledger, request_id, "program-risk-generation-checkpoint");
    let acknowledgement = governor
        .acknowledge_pause(request_id, &ledger, &checkpoint)
        .expect("acknowledge generation zero");
    let activation = governor
        .activate_resume(&acknowledgement)
        .expect("activate generation one");
    governor
        .flush_scope_to_ledger(&open.flush_permit(), &ledger)
        .expect("persist pause/resume lifecycle");

    let published = governor
        .write_program_risk_session_end_report(&ledger, &open, 88, &observations)
        .expect("publish generation-one report");
    assert_eq!(published.receipt.generation(), 1);
    let counts = durable_counts(&ledger);
    drop(governor);
    drop(ledger);

    let ledger = fs_ledger::Ledger::open(&path).expect("reopened generation ledger");
    let governor = Governor::new_durable(&ledger, nonce).expect("reopened governor");
    let recovered_initial_gate = Arc::new(CancelGate::new());
    let recovered_open = governor
        .recover_open(
            &ledger,
            open_id,
            capability,
            Some(Arc::clone(&recovered_initial_gate)),
        )
        .expect("recover gated open");
    assert!(matches!(
        governor.recover_program_risk_report(&ledger, &recovered_open),
        Err(SessionError::ProgramRiskReportGenerationAhead {
            id: 70_006,
            report_generation: 1,
            recovered_generation: 0,
        })
    ));

    assert_eq!(
        governor
            .recover_pressure(&ledger, action_id, 3)
            .expect("recover pause request"),
        pressure
    );
    let recovered_resume_gate = Arc::new(CancelGate::new());
    let recovered_acknowledgement = governor
        .recover_pause_acknowledgement(&ledger, request_id, &checkpoint, recovered_resume_gate)
        .expect("recover pause acknowledgement");
    assert_eq!(
        governor
            .recover_resume_activation(&ledger, &recovered_acknowledgement)
            .expect("recover generation-one activation"),
        activation
    );
    assert_eq!(
        governor
            .recover_program_risk_report(&ledger, &recovered_open)
            .expect("recover report after lifecycle"),
        published.receipt
    );
    assert_eq!(durable_counts(&ledger), counts);
}
