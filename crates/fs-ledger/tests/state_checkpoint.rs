//! G0/G3 tests for portable semantic state-checkpoint receipts.

use fs_blake3::ContentHash;
use fs_ledger::{
    KnownStateSemantics, Ledger, LedgerError, RUNTIME_STATE_ARTIFACT_KIND, StateCheckpointClaim,
    StateSlotId,
};

fn ledger() -> Ledger {
    Ledger::open(":memory:").expect("open in-memory ledger")
}

fn known(code: ContentHash) -> KnownStateSemantics<'static> {
    KnownStateSemantics {
        law_id: "j2-plasticity",
        law_version: 4,
        state_schema_version: 2,
        canonical_parameters_hash: ContentHash([0x18; 32]),
        contract_and_code_hash: code,
    }
}

#[test]
fn exact_retry_dedupes_and_one_slot_accepts_successive_states() {
    let ledger = ledger();
    let slot = StateSlotId::from_content_hash(ContentHash([0x21; 32]));
    let code = ContentHash([0x31; 32]);
    let first_state = ledger
        .put_artifact(RUNTIME_STATE_ARTIFACT_KIND, b"state-step-10", None)
        .expect("retain first state")
        .hash;
    let first_claim = StateCheckpointClaim {
        state_slot: slot,
        semantics: known(code),
        runtime_state_artifact: first_state,
    };
    let first = ledger
        .record_state_checkpoint(first_claim)
        .expect("record first checkpoint");
    assert_eq!(
        ledger
            .record_state_checkpoint(first_claim)
            .expect("exact response-loss retry"),
        first
    );
    assert_eq!(
        ledger
            .table_count("semantic_state_checkpoint_receipts")
            .unwrap(),
        1
    );

    let second_state = ledger
        .put_artifact(RUNTIME_STATE_ARTIFACT_KIND, b"state-step-11", None)
        .expect("retain successor state")
        .hash;
    let second = ledger
        .record_state_checkpoint(StateCheckpointClaim {
            state_slot: slot,
            semantics: known(code),
            runtime_state_artifact: second_state,
        })
        .expect("record successive checkpoint for same slot");
    assert_ne!(first.content_hash(), second.content_hash());
    assert_eq!(first.state_slot(), second.state_slot());
    assert_eq!(
        ledger
            .table_count("semantic_state_checkpoint_receipts")
            .unwrap(),
        2
    );

    let first_loaded = ledger
        .load_state_checkpoint(first.content_hash(), known(code))
        .expect("load first under exact semantics")
        .expect("first checkpoint exists");
    assert_eq!(first_loaded.receipt(), &first);
    assert_eq!(first_loaded.state_bytes(), b"state-step-10");
    let second_loaded = ledger
        .load_state_checkpoint(second.content_hash(), known(code))
        .expect("load successor under exact semantics")
        .expect("second checkpoint exists");
    assert_eq!(second_loaded.receipt(), &second);
    assert_eq!(second_loaded.state_bytes(), b"state-step-11");
    ledger
        .verify_state_checkpoint_receipt(&second, known(code))
        .expect("transport candidate re-earns membership and semantics");
    assert!(ledger.lint().expect("checkpoint hygiene scan").is_clean());
}

#[test]
fn semantic_mismatches_are_typed_and_precede_state_materialization() {
    let ledger = ledger();
    let slot = StateSlotId::from_content_hash(ContentHash([0x22; 32]));
    let code = ContentHash([0x32; 32]);
    let state = ledger
        .put_artifact(RUNTIME_STATE_ARTIFACT_KIND, b"canonical-state", None)
        .expect("retain state")
        .hash;
    let receipt = ledger
        .record_state_checkpoint(StateCheckpointClaim {
            state_slot: slot,
            semantics: known(code),
            runtime_state_artifact: state,
        })
        .expect("record semantic checkpoint");

    let mut bumped_schema = known(code);
    bumped_schema.state_schema_version += 1;
    let error = ledger
        .load_state_checkpoint(receipt.content_hash(), bumped_schema)
        .expect_err("bumped state schema must refuse replay");
    assert_eq!(error.code(), "LedgerUnknownStateSemantics");
    let rendered = error.to_string();
    assert!(rendered.contains("state bytes withheld"));
    assert!(rendered.contains("state_schema_version"));
    let LedgerError::UnknownStateSemantics {
        stored_law,
        stored_law_version,
        stored_state_schema_version,
        expected_state_schema_version,
        differences,
        ..
    } = error
    else {
        panic!("bumped schema must produce typed unknown-semantics context");
    };
    assert_eq!(stored_law, "j2-plasticity");
    assert_eq!(stored_law_version, 4);
    assert_eq!(stored_state_schema_version, 2);
    assert_eq!(expected_state_schema_version, 3);
    assert_eq!(differences, vec!["state_schema_version"]);

    let mut bumped_law = known(code);
    bumped_law.law_version += 1;
    assert!(matches!(
        ledger.load_state_checkpoint(receipt.content_hash(), bumped_law),
        Err(LedgerError::UnknownStateSemantics {
            stored_law_version: 4,
            expected_law_version: 5,
            differences,
            ..
        }) if differences == vec!["law_version"]
    ));

    let mut other_law = known(code);
    other_law.law_id = "foreign-plasticity";
    assert!(matches!(
        ledger.load_state_checkpoint(receipt.content_hash(), other_law),
        Err(LedgerError::UnknownStateSemantics {
            stored_law,
            expected_law,
            differences,
            ..
        }) if stored_law == "j2-plasticity"
            && expected_law == "foreign-plasticity"
            && differences == vec!["law_id"]
    ));

    let mut changed_parameters = known(code);
    changed_parameters.canonical_parameters_hash = ContentHash([0x19; 32]);
    assert!(matches!(
        ledger.load_state_checkpoint(receipt.content_hash(), changed_parameters),
        Err(LedgerError::UnknownStateSemantics { differences, .. })
            if differences == vec!["canonical_parameters_hash"]
    ));
    assert!(matches!(
        ledger.load_state_checkpoint(
            receipt.content_hash(),
            known(ContentHash([0x33; 32]))
        ),
        Err(LedgerError::UnknownStateSemantics { differences, .. })
            if differences == vec!["contract_and_code_hash"]
    ));

    ledger
        .corrupt_artifact_for_test(&state)
        .expect("inject stored-state corruption");
    assert!(matches!(
        ledger.load_state_checkpoint(receipt.content_hash(), bumped_schema),
        Err(LedgerError::UnknownStateSemantics { .. })
    ));
    assert!(matches!(
        ledger.verify_state_checkpoint_receipt(&receipt, bumped_schema),
        Err(LedgerError::UnknownStateSemantics { .. })
    ));
    assert!(matches!(
        ledger.load_state_checkpoint(receipt.content_hash(), known(code)),
        Err(LedgerError::Corrupt { .. })
    ));
}

#[test]
fn invalid_or_missing_state_inputs_publish_no_receipt() {
    let ledger = ledger();
    let slot = StateSlotId::from_content_hash(ContentHash([0x23; 32]));
    let code = ContentHash([0x34; 32]);

    assert!(matches!(
        ledger.record_state_checkpoint(StateCheckpointClaim {
            state_slot: slot,
            semantics: known(code),
            runtime_state_artifact: ContentHash([0x99; 32]),
        }),
        Err(LedgerError::Invalid { field, problem })
            if field == "state_checkpoint.runtime_state_artifact"
                && problem.contains("does not exist")
    ));
    let wrong_kind = ledger
        .put_artifact("generic-bytes", b"state", None)
        .expect("wrong-kind artifact")
        .hash;
    assert!(matches!(
        ledger.record_state_checkpoint(StateCheckpointClaim {
            state_slot: slot,
            semantics: known(code),
            runtime_state_artifact: wrong_kind,
        }),
        Err(LedgerError::Invalid { field, problem })
            if field == "state_checkpoint.runtime_state_artifact"
                && problem.contains("constitutive-runtime-state")
    ));
    assert!(matches!(
        ledger.record_state_checkpoint(StateCheckpointClaim {
            state_slot: StateSlotId::from_content_hash(ContentHash([0; 32])),
            semantics: known(code),
            runtime_state_artifact: wrong_kind,
        }),
        Err(LedgerError::Invalid { field, .. }) if field == "state_checkpoint.state_slot"
    ));
    let mut zero_parameters = known(code);
    zero_parameters.canonical_parameters_hash = ContentHash([0; 32]);
    assert!(matches!(
        ledger.record_state_checkpoint(StateCheckpointClaim {
            state_slot: slot,
            semantics: zero_parameters,
            runtime_state_artifact: wrong_kind,
        }),
        Err(LedgerError::Invalid { field, .. })
            if field == "state_checkpoint.canonical_parameters_hash"
    ));
    assert!(matches!(
        ledger.record_state_checkpoint(StateCheckpointClaim {
            state_slot: slot,
            semantics: known(ContentHash([0; 32])),
            runtime_state_artifact: wrong_kind,
        }),
        Err(LedgerError::Invalid { field, .. })
            if field == "state_checkpoint.contract_and_code_hash"
    ));
    assert_eq!(
        ledger
            .table_count("semantic_state_checkpoint_receipts")
            .unwrap(),
        0
    );
}

#[test]
fn runtime_state_artifacts_are_gc_roots() {
    let ledger = ledger();
    let state = ledger
        .put_artifact(RUNTIME_STATE_ARTIFACT_KIND, b"rooted state", None)
        .expect("state artifact")
        .hash;
    let receipt = ledger
        .record_state_checkpoint(StateCheckpointClaim {
            state_slot: StateSlotId::from_content_hash(ContentHash([0x24; 32])),
            semantics: known(ContentHash([0x35; 32])),
            runtime_state_artifact: state,
        })
        .expect("checkpoint root");
    let unrelated = ledger
        .put_artifact("scratch", b"collect me", None)
        .expect("unrelated artifact")
        .hash;
    let report = ledger
        .gc_unreferenced_artifacts(false)
        .expect("collect unreachable artifacts");
    assert!(report.candidates.contains(&unrelated.to_hex()));
    assert!(!report.candidates.contains(&state.to_hex()));
    assert!(ledger.get_artifact(&unrelated).unwrap().is_none());
    assert!(
        ledger
            .load_state_checkpoint(receipt.content_hash(), known(ContentHash([0x35; 32])))
            .unwrap()
            .is_some()
    );
}
