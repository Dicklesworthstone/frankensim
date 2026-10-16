//! G0/G3 tests for immutable quantity-schema migration receipts.

use fs_blake3::ContentHash;
use fs_ledger::{Ledger, LedgerError};
use fs_qty::json::{DimensionCrosswalkReceipt, decode_json, to_json};

fn ledger() -> Ledger {
    Ledger::open(":memory:").expect("open in-memory ledger")
}

fn retain_crosswalk_artifacts(
    ledger: &Ledger,
    old_json: &str,
) -> (DimensionCrosswalkReceipt, ContentHash, ContentHash) {
    let decoded = decode_json(old_json).expect("canonical legacy quantity JSON");
    let receipt = decoded
        .migration()
        .cloned()
        .expect("legacy decode must carry migration evidence");
    let new_json = to_json(decoded.qty()).expect("canonical six-base quantity JSON");
    let old = ledger
        .put_artifact("quantity-json-v1", old_json.as_bytes(), None)
        .expect("retain exact historical bytes");
    let new = ledger
        .put_artifact("quantity-json-v2", new_json.as_bytes(), None)
        .expect("retain exact canonical bytes");
    assert_eq!(receipt.old_hash(), old.hash);
    assert_eq!(receipt.new_hash(), new.hash);
    (receipt, old.hash, new.hash)
}

#[test]
fn qty_crosswalk_round_trips_dedupes_and_preserves_distinct_history() {
    const IMPLICIT_V1: &str = r#"{"value":0.12,"dims":[-1,1,-1,0,0]}"#;
    const EXPLICIT_V1: &str = r#"{"schema_version":1,"value":0.12,"dims":[-1,1,-1,0,0]}"#;
    const CANONICAL_V2: &str = r#"{"schema_version":2,"value":0.12,"dims":[-1,1,-1,0,0,0]}"#;
    const IMPLICIT_HASH: &str = "b97ca96f12cf487bc90760adad7257311fed950f95ab834c9107e51bf5f31ef1";
    const EXPLICIT_HASH: &str = "14a9d9ab56a9aa25c3ada13180ae342886bde302f5848d704edb4d5e1253cf24";
    const CANONICAL_HASH: &str = "8353a2a85f0de4a46f8cb31cb1673198c9bae9526b848369be545031d495bbb5";
    let ledger = ledger();
    let (implicit, implicit_old, canonical_new) = retain_crosswalk_artifacts(&ledger, IMPLICIT_V1);
    let (explicit, explicit_old, explicit_new) = retain_crosswalk_artifacts(&ledger, EXPLICIT_V1);
    assert_eq!(
        to_json(decode_json(IMPLICIT_V1).unwrap().qty()).unwrap(),
        CANONICAL_V2
    );
    assert_eq!(implicit_old, ContentHash::from_hex(IMPLICIT_HASH).unwrap());
    assert_eq!(explicit_old, ContentHash::from_hex(EXPLICIT_HASH).unwrap());
    assert_eq!(
        canonical_new,
        ContentHash::from_hex(CANONICAL_HASH).unwrap()
    );
    assert_ne!(implicit_old, explicit_old, "historical bytes stay distinct");
    assert_eq!(
        canonical_new, explicit_new,
        "both histories converge exactly"
    );

    let first = ledger
        .record_qty_dimension_crosswalk(&implicit)
        .expect("record implicit-v1 migration");
    assert_eq!(first.old_hash(), implicit_old);
    assert_eq!(first.new_hash(), canonical_new);
    assert!(!first.deduped());
    let replay = ledger
        .record_qty_dimension_crosswalk(&implicit)
        .expect("exact retry after response loss");
    assert!(replay.deduped());

    ledger
        .record_qty_dimension_crosswalk(&explicit)
        .expect("record explicit-v1 migration");
    assert_eq!(ledger.table_count("qty_dimension_crosswalks").unwrap(), 2);
    assert_eq!(
        ledger
            .qty_dimension_crosswalk(implicit_old)
            .expect("verified lookup"),
        Some(implicit)
    );
    assert_eq!(
        ledger
            .qty_dimension_crosswalk(explicit_old)
            .expect("verified lookup"),
        Some(explicit)
    );
    assert!(ledger.lint().expect("crosswalk hygiene scan").is_clean());
    assert_eq!(
        ledger
            .qty_dimension_crosswalk(ContentHash([0xA5; 32]))
            .expect("missing lookup"),
        None
    );
}

#[test]
fn qty_crosswalk_refuses_missing_or_corrupt_retained_bytes() {
    const OLD: &str = r#"{"value":3.5,"dims":[1,0,-1,0,0]}"#;
    let ledger = ledger();
    let decoded = decode_json(OLD).expect("legacy quantity");
    let receipt = decoded.migration().cloned().expect("migration receipt");
    ledger
        .put_artifact("quantity-json-v1", OLD.as_bytes(), None)
        .expect("retain source only");
    let Err(LedgerError::Invalid { field, problem }) =
        ledger.record_qty_dimension_crosswalk(&receipt)
    else {
        panic!("missing target must be a typed invalid-input refusal");
    };
    assert_eq!(field, "qty_dimension_crosswalk.new_hash");
    assert!(problem.contains(&receipt.new_hash().to_hex()));
    assert!(problem.contains("does not exist"));
    assert_eq!(ledger.table_count("qty_dimension_crosswalks").unwrap(), 0);

    let new_json = to_json(decoded.qty()).expect("canonical target");
    ledger
        .put_artifact("quantity-json-v2", new_json.as_bytes(), None)
        .expect("retain target");
    ledger
        .record_qty_dimension_crosswalk(&receipt)
        .expect("record after both artifacts exist");
    ledger
        .corrupt_artifact_for_test(&receipt.old_hash())
        .expect("inject stored-byte corruption");
    let Err(LedgerError::Corrupt { hash_hex, detail }) =
        ledger.qty_dimension_crosswalk(receipt.old_hash())
    else {
        panic!("tampered retained source must be a typed corruption refusal");
    };
    assert_eq!(hash_hex, receipt.old_hash().to_hex());
    assert!(detail.contains("content hash mismatch"));
}

#[test]
fn qty_crosswalk_recording_is_transaction_owned_and_atomic() {
    const OLD: &str = r#"{"value":1,"dims":[0,0,0,0,0]}"#;
    let ledger = ledger();
    let (receipt, old_hash, _) = retain_crosswalk_artifacts(&ledger, OLD);

    ledger.begin().expect("caller transaction");
    let Err(LedgerError::Invalid { field, problem }) =
        ledger.record_qty_dimension_crosswalk(&receipt)
    else {
        panic!("caller transaction must receive a typed ownership refusal");
    };
    assert_eq!(field, "qty_dimension_crosswalk.transaction");
    assert!(problem.contains("must own its transaction"));
    ledger.rollback().expect("rollback caller transaction");
    assert_eq!(ledger.table_count("qty_dimension_crosswalks").unwrap(), 0);

    ledger
        .record_qty_dimension_crosswalk(&receipt)
        .expect("owned transaction succeeds");
    assert!(
        ledger
            .qty_dimension_crosswalk(old_hash)
            .expect("lookup")
            .is_some()
    );
}

#[test]
fn qty_crosswalk_artifacts_are_gc_roots() {
    const OLD: &str = r#"{"value":2,"dims":[2,0,-2,0,0]}"#;
    let ledger = ledger();
    let (receipt, old_hash, new_hash) = retain_crosswalk_artifacts(&ledger, OLD);
    ledger
        .record_qty_dimension_crosswalk(&receipt)
        .expect("record crosswalk");
    let unrelated = ledger
        .put_artifact("scratch", b"unreferenced", None)
        .expect("unreferenced fixture")
        .hash;

    let dry = ledger
        .gc_unreferenced_artifacts(true)
        .expect("dry-run garbage collection");
    assert!(dry.candidates.contains(&unrelated.to_hex()));
    assert!(!dry.candidates.contains(&old_hash.to_hex()));
    assert!(!dry.candidates.contains(&new_hash.to_hex()));

    let live = ledger
        .gc_unreferenced_artifacts(false)
        .expect("collect only unreferenced artifact");
    assert_eq!(live.deleted, 1);
    assert!(ledger.get_artifact(&unrelated).unwrap().is_none());
    assert!(ledger.get_artifact(&old_hash).unwrap().is_some());
    assert!(ledger.get_artifact(&new_hash).unwrap().is_some());
    assert_eq!(
        ledger
            .qty_dimension_crosswalk(old_hash)
            .expect("crosswalk remains replayable"),
        Some(receipt)
    );
}
