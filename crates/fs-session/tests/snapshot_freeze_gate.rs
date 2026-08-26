//! Integration coverage for predecessor-bound snapshot freeze gating
//! (bead frankensim-sj31i.52.5.4.1, session lane).
//!
//! A pause that was explicitly declared snapshot-freeze-bound may not
//! activate resume until a committed fs-exec freeze receipt — validated
//! against THIS governor and THE completed pause — is published as an
//! immutable ledger terminal row. Pauses nobody armed keep their legacy
//! activation path.

use fs_exec::freeze::{
    FreezeBoundaryLabels, FreezeOwnerBinding, FreezeResumeInputs, SnapshotFreezeRegistry,
};
use fs_session::{
    CapabilityToken, Governor, SessionError, SessionId, SnapshotFreezePublicationDisposition,
};
use std::sync::Arc;

fn token(session: u64, ledger_scope: &str) -> CapabilityToken {
    CapabilityToken {
        session: SessionId(session),
        ops: vec!["flux.*".to_string()],
        core_s: 1.0e9,
        mem_bytes: u64::MAX,
        wall_s: 1.0e9,
        cores: 1,
        ledger_scope: ledger_scope.to_string(),
    }
}

fn solver_checkpoint(
    ledger: &fs_ledger::Ledger,
    request: fs_session::PauseRequestId,
    label: &str,
) -> fs_ledger::session_registry::SolverCheckpointReceipt {
    let authority = request.checkpoint_authority();
    let mut run_bytes = [0_u8; 8];
    run_bytes.copy_from_slice(&authority.as_bytes()[..8]);
    let run = fs_exec::RunId(u64::from_le_bytes(run_bytes));
    let gate = fs_exec::CancelGate::new_clock_free();
    let tracker = fs_exec::DrainTracker::new(run, &gate);
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

struct GatedFixture {
    governor: Arc<Governor>,
    gate: Arc<fs_exec::CancelGate>,
    ledger: fs_ledger::Ledger,
    acknowledgement: fs_session::PauseAcknowledgement,
}

fn gated_fixture(session: u64, scope: &'static str) -> GatedFixture {
    let governor = Arc::new(Governor::new());
    let gate = Arc::new(fs_exec::CancelGate::new());
    let open_id = governor
        .session_open_id(SessionId(session), scope)
        .expect("fixture open id");
    governor
        .open_session_declared_gated(open_id, token(session, scope), Arc::clone(&gate))
        .expect("fixture gated session");
    let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
    gate.request();
    let events = governor
        .apply_memory_pressure(SessionId(session), 3)
        .expect("pause request minted");
    let request_id = events
        .last()
        .and_then(|event| event.pause_request_id)
        .expect("pause request authority");
    let checkpoint = solver_checkpoint(&ledger, request_id, "freeze-gate-fixture");
    let acknowledgement = governor
        .acknowledge_pause(request_id, &ledger, &checkpoint)
        .expect("pause acknowledged");
    GatedFixture {
        governor,
        gate,
        ledger,
        acknowledgement,
    }
}

/// Drive one full fs-exec freeze transaction bound to `governor` + `session`
/// and return the committed typed receipt.
fn commit_freeze(
    governor: &Governor,
    session: u64,
    nonce_tag: u8,
) -> fs_exec::freeze::SnapshotFreezeReceipt {
    use fs_exec::solver::snapshot_v2 as sv2;

    struct FixtureState {
        marker: u64,
    }
    impl fs_exec::solver::SolverStateV2 for FixtureState {
        const STATE_TYPE_ID_V2: sv2::SnapshotStateTypeIdV2 =
            sv2::SnapshotStateTypeIdV2::from_bytes([0xE1; 32]);
        const STATE_SCHEMA_ID_V2: sv2::SnapshotStateSchemaIdV2 =
            sv2::SnapshotStateSchemaIdV2::from_bytes([0xE2; 32]);
        const STATE_CODEC_ID_V2: sv2::SnapshotStateCodecIdV2 =
            sv2::SnapshotStateCodecIdV2::from_bytes([0xE3; 32]);
        const STATE_CODEC_VERSION_V2: u32 = 1;
        fn encode_v2(
            &self,
            encoder: &mut sv2::SnapshotEncoderV2<'_>,
        ) -> Result<(), sv2::SnapshotV2Error> {
            encoder.put_u64(self.marker)
        }
        fn decode_v2(
            decoder: &mut sv2::SnapshotDecoderV2<'_, '_>,
        ) -> Result<Self, sv2::SnapshotV2Error> {
            Ok(Self {
                marker: decoder.get_u64()?,
            })
        }
    }

    let binding = FreezeOwnerBinding {
        owner: *governor.identity().as_bytes(),
        session,
        solver_instance_generation: 1,
        instance_nonce: [nonce_tag; 32],
    };
    let registry = SnapshotFreezeRegistry::new(binding).expect("fixture registry");
    let labels = FreezeBoundaryLabels {
        pause_request: sv2::SnapshotPauseRequestIdV2::from_bytes([0xC1; 32]),
        gate_generation: 1,
    };
    let inputs = FreezeResumeInputs {
        algorithm: sv2::SnapshotAlgorithmIdV2::from_bytes([0xD1; 32]),
        algorithm_version: 1,
        problem: sv2::SnapshotProblemIdV2::from_bytes([0xD2; 32]),
        rng_counter: sv2::SnapshotRngCounterIdV2::from_bytes([0xD3; 32]),
        determinism: sv2::SnapshotDeterminismV2::Deterministic,
        execution_fingerprint: sv2::SnapshotExecutionFingerprintIdV2::from_bytes([0xD4; 32]),
        budget_state: sv2::SnapshotBudgetStateIdV2::from_bytes([0xD5; 32]),
        provenance: sv2::SnapshotProvenanceIdV2::from_bytes([0xD6; 32]),
    };
    let limits = sv2::SnapshotLimitsV2::new(
        1 << 20,
        64,
        fs_blake3::identity::CanonicalLimits::new(16_384, 4_096, 32, 32, 64),
        4_096,
        1 << 20,
        64,
    );
    let run_gate = fs_exec::CancelGate::new();
    run_gate.request();
    let tracker = fs_exec::DrainTracker::new(fs_exec::RunId(9), &run_gate);
    drop(tracker.register_worker().expect("fixture freeze worker"));
    let request = registry.begin_freeze().expect("freeze admitted");
    let cancellation = || false;
    let permit = request
        .freeze(
            FixtureState { marker: 7 },
            &tracker,
            labels,
            inputs,
            limits,
            &mut cancellation,
        )
        .expect("state frozen");
    permit.seal(cancellation).expect("envelope sealed")
}

#[test]
fn gated_resume_refuses_until_receipt_is_recorded_then_activates() {
    let fixture = gated_fixture(60, "freeze-gate-happy");
    fixture
        .governor
        .require_snapshot_freeze(&fixture.acknowledgement)
        .expect("gate armed");

    // Activation before any receipt must refuse.
    assert!(matches!(
        fixture.governor.activate_resume(&fixture.acknowledgement),
        Err(SessionError::SnapshotFreezeReceiptRequired { id: 60 })
    ));

    let freeze = commit_freeze(&fixture.governor, 60, 0x5A);
    let write = fixture
        .governor
        .record_snapshot_freeze_receipt(&fixture.ledger, 1_000, &fixture.acknowledgement, &freeze)
        .expect("receipt recorded");
    assert_eq!(
        write.disposition,
        SnapshotFreezePublicationDisposition::Committed
    );

    // Idempotent republication replays without side effects.
    let replayed = fixture
        .governor
        .record_snapshot_freeze_receipt(&fixture.ledger, 1_001, &fixture.acknowledgement, &freeze)
        .expect("receipt replayed");
    assert_eq!(
        replayed.disposition,
        SnapshotFreezePublicationDisposition::Replayed
    );

    fixture
        .governor
        .activate_resume(&fixture.acknowledgement)
        .expect("satisfied gate activates");
}

#[test]
fn foreign_session_receipt_fails_closed() {
    let fixture = gated_fixture(61, "freeze-gate-foreign");
    fixture
        .governor
        .require_snapshot_freeze(&fixture.acknowledgement)
        .expect("gate armed");
    // Bound to a DIFFERENT session number than the paused one.
    let freeze = commit_freeze(&fixture.governor, 999, 0x5B);
    assert!(matches!(
        fixture.governor.record_snapshot_freeze_receipt(
            &fixture.ledger,
            1,
            &fixture.acknowledgement,
            &freeze
        ),
        Err(SessionError::SnapshotFreezeBindingMismatch { id: 61, .. })
    ));
    // The gate remains unsatisfied.
    assert!(matches!(
        fixture.governor.activate_resume(&fixture.acknowledgement),
        Err(SessionError::SnapshotFreezeReceiptRequired { id: 61 })
    ));
}

#[test]
fn un_armed_pause_refuses_receipt_publication() {
    let fixture = gated_fixture(62, "freeze-gate-unarmed");
    let freeze = commit_freeze(&fixture.governor, 62, 0x5C);
    assert!(matches!(
        fixture.governor.record_snapshot_freeze_receipt(
            &fixture.ledger,
            1,
            &fixture.acknowledgement,
            &freeze
        ),
        Err(SessionError::SnapshotFreezeBindingMismatch { id: 62, .. })
    ));
}

#[test]
fn legacy_pauses_activate_without_any_freeze_machinery() {
    // Nobody armed the gate: the pre-existing pause/ack/activate path must be
    // untouched so snapshotless resumes keep working.
    let fixture = gated_fixture(63, "freeze-gate-legacy");
    fixture
        .governor
        .activate_resume(&fixture.acknowledgement)
        .expect("unarmed legacy path activates");
}
