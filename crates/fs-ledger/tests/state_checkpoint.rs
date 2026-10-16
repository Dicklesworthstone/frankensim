//! G0/G3/G4 tests for portable semantic state-checkpoint receipts.

use std::sync::atomic::{AtomicU32, Ordering};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use fs_blake3::ContentHash;
use fs_ledger::{
    KnownStateSemantics, Ledger, LedgerError, RUNTIME_STATE_ARTIFACT_KIND, StateCheckpointClaim,
    StateCheckpointReceipt, StateSlotId,
};

static NEXT_DB: AtomicU32 = AtomicU32::new(0);
const MAX_CRASH_RECEIPT_INDEX: u64 = 64;

struct ChildGuard(Option<std::process::Child>);

impl ChildGuard {
    fn child_mut(&mut self) -> &mut std::process::Child {
        self.0.as_mut().expect("child remains armed")
    }

    fn kill_running_and_wait(&mut self) -> std::io::Result<()> {
        let Some(child) = self.0.as_mut() else {
            return Ok(());
        };
        if let Some(status) = child.try_wait()? {
            self.0.take();
            return Err(std::io::Error::other(format!(
                "child exited before requested kill: {status}"
            )));
        }
        child.kill()?;
        let _ = child.wait()?;
        self.0.take();
        Ok(())
    }
}

impl Drop for ChildGuard {
    fn drop(&mut self) {
        if let Some(mut child) = self.0.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

fn temp_db(tag: &str) -> String {
    let sequence = NEXT_DB.fetch_add(1, Ordering::Relaxed);
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_nanos());
    std::env::temp_dir()
        .join(format!(
            "fs-ledger-state-checkpoint-{tag}-{}-{nonce}-{sequence}.db",
            std::process::id(),
        ))
        .display()
        .to_string()
}

fn cleanup_db(path: &str) {
    for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
        let _ = std::fs::remove_file(format!("{path}{suffix}"));
    }
}

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

fn checkpoint_bytes(index: u64) -> Vec<u8> {
    format!("canonical-state-step-{index:020}").into_bytes()
}

fn record_indexed_checkpoint(
    ledger: &Ledger,
    index: u64,
) -> Result<StateCheckpointReceipt, LedgerError> {
    let state = ledger
        .put_artifact(
            RUNTIME_STATE_ARTIFACT_KIND,
            &checkpoint_bytes(index),
            Some(r#"{"codec":"test-u64-index-v1"}"#),
        )?
        .hash;
    ledger.record_state_checkpoint(StateCheckpointClaim {
        state_slot: StateSlotId::from_content_hash(ContentHash([0x51; 32])),
        semantics: known(ContentHash([0x61; 32])),
        runtime_state_artifact: state,
    })
}

#[test]
fn state_checkpoint_crash_child_writer() {
    let Ok(path) = std::env::var("FS_LEDGER_STATE_CHECKPOINT_CRASH_DB") else {
        return;
    };
    let ledger = Ledger::open(&path).expect("crash child opens pre-migrated ledger");
    let mut index = 1u64;
    loop {
        record_indexed_checkpoint(&ledger, index).expect("crash child commits checkpoint");
        if index < MAX_CRASH_RECEIPT_INDEX {
            index += 1;
        }
    }
}

#[test]
fn committed_checkpoint_prefix_survives_kill_and_real_file_reopen() {
    let path = temp_db("crash-reopen");
    let (instance, baseline, baseline_transport) = {
        let ledger = Ledger::open(&path).expect("create file-backed ledger");
        let receipt = record_indexed_checkpoint(&ledger, 0).expect("commit baseline checkpoint");
        (ledger.instance_id(), receipt.clone(), receipt.to_bytes())
    };

    {
        let reopened = Ledger::open(&path).expect("reopen baseline ledger");
        assert_eq!(reopened.instance_id(), instance);
        let decoded = StateCheckpointReceipt::from_bytes(&baseline_transport)
            .expect("portable baseline receipt decodes after reopen");
        assert_eq!(decoded, baseline);
        reopened
            .verify_state_checkpoint_receipt(&decoded, known(ContentHash([0x61; 32])))
            .expect("reopened ledger revalidates portable receipt membership");
        let loaded = reopened
            .load_state_checkpoint(decoded.content_hash(), known(ContentHash([0x61; 32])))
            .expect("baseline replay query")
            .expect("baseline checkpoint survives reopen");
        assert_eq!(loaded.receipt(), &baseline);
        assert_eq!(loaded.state_bytes(), checkpoint_bytes(0));
    }

    let observer = Ledger::open(&path).expect("open crash observer");
    let deadline = Instant::now() + Duration::from_secs(5);
    let executable = std::env::current_exe().expect("current integration-test executable");
    let mut child = ChildGuard(Some(
        std::process::Command::new(executable)
            .args([
                "--exact",
                "state_checkpoint_crash_child_writer",
                "--nocapture",
            ])
            .env("FS_LEDGER_STATE_CHECKPOINT_CRASH_DB", &path)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .spawn()
            .expect("spawn checkpoint crash child"),
    ));
    loop {
        let committed = observer
            .table_count("semantic_state_checkpoint_receipts")
            .expect("observe committed checkpoint prefix");
        if committed >= 2 {
            break;
        }
        if let Some(status) = child.child_mut().try_wait().expect("poll crash child") {
            cleanup_db(&path);
            panic!("checkpoint crash child exited before kill: {status}");
        }
        if Instant::now() >= deadline {
            child
                .kill_running_and_wait()
                .expect("stop and reap wedged crash child");
            cleanup_db(&path);
            panic!("checkpoint crash child published no post-baseline receipt");
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    // This deliberately proves recovery after abrupt writer death. Without a
    // production transaction hook it does not claim the signal landed inside
    // one particular SQLite transaction.
    child
        .kill_running_and_wait()
        .expect("kill and reap checkpoint writer during traffic");
    drop(observer);

    let recovered = Ledger::open(&path).expect("recover checkpoint ledger after kill");
    assert_eq!(recovered.instance_id(), instance);
    let receipt_count = recovered
        .table_count("semantic_state_checkpoint_receipts")
        .expect("count recovered checkpoint prefix");
    assert!((2..=MAX_CRASH_RECEIPT_INDEX + 1).contains(&receipt_count));

    let oracle = ledger();
    let mut last_expected = None;
    for index in 0..receipt_count {
        let expected =
            record_indexed_checkpoint(&oracle, index).expect("build independent receipt");
        let loaded = recovered
            .load_state_checkpoint(expected.content_hash(), known(ContentHash([0x61; 32])))
            .expect("query recovered checkpoint")
            .unwrap_or_else(|| panic!("committed checkpoint index {index} is missing"));
        assert_eq!(loaded.receipt(), &expected, "receipt index {index}");
        assert_eq!(
            loaded.state_bytes(),
            checkpoint_bytes(index),
            "state bytes index {index}"
        );
        recovered
            .verify_state_checkpoint_receipt(&expected, known(ContentHash([0x61; 32])))
            .unwrap_or_else(|error| panic!("receipt index {index} failed verification: {error}"));
        last_expected = Some(expected);
    }

    let last_index = receipt_count - 1;
    let last_expected = last_expected.expect("nonempty recovered receipt prefix");
    assert_eq!(
        record_indexed_checkpoint(&recovered, last_index).expect("retry last checkpoint"),
        last_expected
    );
    assert_eq!(
        recovered
            .table_count("semantic_state_checkpoint_receipts")
            .expect("count after exact retry"),
        receipt_count
    );

    let mut bumped = known(ContentHash([0x61; 32]));
    bumped.state_schema_version += 1;
    let refusal = recovered
        .load_state_checkpoint(last_expected.content_hash(), bumped)
        .expect_err("recovered checkpoint must refuse bumped state schema");
    let rendered_refusal = refusal.to_string();
    assert!(rendered_refusal.contains("state bytes were withheld"));
    let LedgerError::UnknownStateSemantics {
        stored_law,
        stored_law_version,
        stored_state_schema_version,
        expected_law,
        expected_law_version,
        expected_state_schema_version,
        stored_parameters_hash,
        expected_parameters_hash,
        stored_contract_and_code_hash,
        expected_contract_and_code_hash,
        differences,
        ..
    } = refusal
    else {
        panic!("bumped recovered state schema must return forensic semantic context");
    };
    assert_eq!(stored_law, "j2-plasticity");
    assert_eq!(stored_law, expected_law);
    assert_eq!(stored_law_version, 4);
    assert_eq!(stored_law_version, expected_law_version);
    assert_eq!(stored_state_schema_version, 2);
    assert_eq!(expected_state_schema_version, 3);
    assert_eq!(stored_parameters_hash, ContentHash([0x18; 32]).to_hex());
    assert_eq!(stored_parameters_hash, expected_parameters_hash);
    assert_eq!(
        stored_contract_and_code_hash,
        ContentHash([0x61; 32]).to_hex()
    );
    assert_eq!(
        stored_contract_and_code_hash,
        expected_contract_and_code_hash
    );
    assert_eq!(differences, vec!["state_schema_version"]);
    assert!(recovered.lint().expect("post-crash lint").is_clean());
    assert!(
        recovered
            .verify_artifact_integrity()
            .expect("post-crash artifact scan")
            .is_clean()
    );
    println!(
        "{{\"suite\":\"fs-ledger/state-checkpoint\",\"case\":\"crash-reopen\",\
         \"recovered_receipts\":{receipt_count},\"last_receipt\":\"{}\",\
         \"law\":\"j2-plasticity\",\"law_version\":4,\"state_schema_version\":2,\
         \"instance_stable\":true,\"semantic_refusal_revalidated\":true,\
         \"lint_clean\":true,\"integrity_clean\":true}}",
        last_expected.content_hash().to_hex()
    );
    drop(recovered);
    cleanup_db(&path);
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
    assert!(rendered.contains("state bytes were withheld"));
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
