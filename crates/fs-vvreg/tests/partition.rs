//! G0/G3 adversarial battery for purpose-typed V&V corpus access and
//! transitive calibration-taint enforcement (EXTREAL f85xj.7.1).

use fs_blake3::{ContentHash, hash_domain};
use fs_ledger::{Ledger, LedgerError};
use fs_qty::QtyAny;
use fs_vvreg::corpus::{ContextValue, DatasetPartition, corpus};
use fs_vvreg::partition::{
    DatasetPurpose, MAX_PARTITION_RECORD_BYTES, PARTITION_RECEIPT_ARTIFACT_KIND, PartitionLedger,
    PartitionReceiptRecord, PartitionRecordError, PartitionRefusal,
};
use std::sync::atomic::{AtomicU64, Ordering};

const CHT: &str = "fs-benchmark-cht-query-v1";
const MARTIN_MOYCE: &str = "martin-moyce-1952-square-column";
const PARTITION_RECORD_META: &str = r#"{"schema":"fs-vvreg-partition-receipt-wire-v1"}"#;

static NEXT_LEDGER: AtomicU64 = AtomicU64::new(0);

fn artifact(label: &str) -> ContentHash {
    hash_domain("org.frankensim.fs-vvreg.test-model.v1", label.as_bytes())
}

fn cht_context() -> [ContextValue; 1] {
    [ContextValue {
        name: "reference_cost_work_units".to_string(),
        value: QtyAny::dimensionless(250.0),
    }]
}

fn martin_moyce_context() -> [ContextValue; 1] {
    martin_moyce_context_at(1.0)
}

fn martin_moyce_context_at(t_star: f64) -> [ContextValue; 1] {
    [ContextValue {
        name: "t_star".to_string(),
        value: QtyAny::dimensionless(t_star),
    }]
}

fn receipt_records() -> Vec<PartitionReceiptRecord> {
    let mut partitions = PartitionLedger::capture(corpus());
    let repartition = partitions
        .repartition(
            CHT,
            DatasetPartition::Calibration,
            "reserve the exact CHT row for model calibration",
        )
        .expect("repartition receipt");
    let calibration = corpus()
        .query(
            &partitions,
            CHT,
            DatasetPurpose::Calibration,
            &cht_context(),
        )
        .expect("calibration access receipt");
    let model = partitions
        .register_model(artifact("durable-model"), &[&calibration], &[])
        .expect("model taint receipt");
    let held_out = corpus()
        .query(
            &partitions,
            MARTIN_MOYCE,
            DatasetPurpose::Validation,
            &martin_moyce_context(),
        )
        .expect("held-out access");
    let validation = partitions
        .validate_model(&model, &[&held_out])
        .expect("disjoint validation receipt");
    partitions
        .repartition(
            MARTIN_MOYCE,
            DatasetPartition::BlindHoldout,
            "seal the validation row for a later blind drill",
        )
        .expect("blind repartition");
    let blind_release = partitions
        .release_blind(
            MARTIN_MOYCE,
            artifact("durable-preregistration"),
            artifact("durable-blind-manifest"),
            "release after the frozen model and protocol were committed",
        )
        .expect("blind release receipt");

    vec![
        calibration.receipt().clone().into(),
        repartition.into(),
        blind_release.into(),
        model.into(),
        validation.into(),
    ]
}

fn partition_ledger_path() -> String {
    let sequence = NEXT_LEDGER.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!(
            "fs-vvreg-partition-receipts-{}-{sequence}.db",
            std::process::id()
        ))
        .display()
        .to_string()
}

#[test]
fn purpose_is_checked_independently_of_the_declared_partition_name() {
    let mut partitions = PartitionLedger::capture(corpus());
    let validation = corpus()
        .query(&partitions, CHT, DatasetPurpose::Validation, &cht_context())
        .expect("seeded validation row admits validation use");
    assert_eq!(
        validation.receipt().partition(),
        DatasetPartition::Validation
    );
    assert_eq!(validation.receipt().purpose(), DatasetPurpose::Validation);

    for attempted in [DatasetPurpose::Calibration, DatasetPurpose::BlindEvaluation] {
        assert!(matches!(
            corpus().query(&partitions, CHT, attempted, &cht_context()),
            Err(PartitionRefusal::PurposeMismatch {
                dataset_id,
                declared: DatasetPartition::Validation,
                attempted: observed,
            }) if dataset_id == CHT && observed == attempted
        ));
    }

    partitions
        .repartition(CHT, DatasetPartition::Training, "fit-only training split")
        .expect("justified repartition");
    let training = corpus()
        .query(
            &partitions,
            CHT,
            DatasetPurpose::Calibration,
            &cht_context(),
        )
        .expect("training rows enter the calibration-taint purpose");
    assert_eq!(training.receipt().partition(), DatasetPartition::Training);
    assert!(matches!(
        corpus().query(&partitions, CHT, DatasetPurpose::Validation, &cht_context(),),
        Err(PartitionRefusal::PurposeMismatch { .. })
    ));
}

#[test]
fn direct_and_transitive_validation_laundering_are_refused_with_paths() {
    let mut partitions = PartitionLedger::capture(corpus());
    let to_calibration = partitions
        .repartition(
            CHT,
            DatasetPartition::Calibration,
            "reserve the exact CHT row for model calibration",
        )
        .expect("repartition to calibration");
    assert!(to_calibration.stales_validation_claims());
    let calibration = corpus()
        .query(
            &partitions,
            CHT,
            DatasetPurpose::Calibration,
            &cht_context(),
        )
        .expect("calibration access");

    let model_a_id = artifact("model-a");
    let model_b_id = artifact("model-b");
    let model_a = partitions
        .register_model(model_a_id, &[&calibration], &[])
        .expect("direct model taint");
    let model_b = partitions
        .register_model(model_b_id, &[], &[&model_a])
        .expect("transitive model taint");

    let to_validation = partitions
        .repartition(
            CHT,
            DatasetPartition::Validation,
            "freeze the model and move the row to held-out evaluation",
        )
        .expect("repartition to validation");
    assert!(!to_validation.stales_validation_claims());
    let validation = corpus()
        .query(&partitions, CHT, DatasetPurpose::Validation, &cht_context())
        .expect("validation access");

    let direct = partitions
        .validate_model(&model_a, &[&validation])
        .expect_err("direct reuse must refuse");
    assert!(matches!(
        direct,
        PartitionRefusal::TaintIntersection {
            model_artifact,
            dataset_id,
            model_path,
            ..
        } if model_artifact == model_a_id && dataset_id == CHT && model_path == vec![model_a_id]
    ));

    let transitive = partitions
        .validate_model(&model_b, &[&validation])
        .expect_err("transitive laundering must refuse");
    assert!(matches!(
        transitive,
        PartitionRefusal::TaintIntersection {
            model_artifact,
            dataset_id,
            model_path,
            ..
        } if model_artifact == model_b_id
            && dataset_id == CHT
            && model_path == vec![model_b_id, model_a_id]
    ));
}

#[test]
fn disjoint_held_out_data_mints_only_a_taint_check_receipt() {
    let mut partitions = PartitionLedger::capture(corpus());
    partitions
        .repartition(CHT, DatasetPartition::Calibration, "calibration fixture")
        .unwrap();
    let calibration = corpus()
        .query(
            &partitions,
            CHT,
            DatasetPurpose::Calibration,
            &cht_context(),
        )
        .unwrap();
    let model = partitions
        .register_model(artifact("clean-model"), &[&calibration], &[])
        .unwrap();
    let held_out = corpus()
        .query(
            &partitions,
            MARTIN_MOYCE,
            DatasetPurpose::Validation,
            &martin_moyce_context(),
        )
        .unwrap();
    let receipt = partitions
        .validate_model(&model, &[&held_out])
        .expect("disjoint held-out input");
    assert_eq!(receipt.model_taint(), model.identity());
    assert_eq!(
        receipt.evaluation_accesses(),
        &[held_out.receipt().identity()]
    );
    assert_ne!(receipt.identity(), model.identity());
}

#[test]
fn repartition_invalidates_access_and_records_claim_staleness() {
    let mut partitions = PartitionLedger::capture(corpus());
    let old_validation = corpus()
        .query(&partitions, CHT, DatasetPurpose::Validation, &cht_context())
        .unwrap();
    let first = partitions
        .repartition(
            CHT,
            DatasetPartition::Calibration,
            "new calibration campaign uses the former validation row",
        )
        .unwrap();
    assert_eq!(first.generation(), 1);
    assert!(first.stales_validation_claims());
    let calibration = corpus()
        .query(
            &partitions,
            CHT,
            DatasetPurpose::Calibration,
            &cht_context(),
        )
        .unwrap();
    let model = partitions
        .register_model(artifact("stale-check"), &[&calibration], &[])
        .unwrap();
    assert!(matches!(
        partitions.validate_model(&model, &[&old_validation]),
        Err(PartitionRefusal::StaleAccess {
            receipt_generation: 0,
            current_generation: 1,
            ..
        })
    ));
    assert_eq!(partitions.events(), &[first]);
}

#[test]
fn stale_calibration_access_cannot_register_a_model() {
    let mut partitions = PartitionLedger::capture(corpus());
    partitions
        .repartition(CHT, DatasetPartition::Calibration, "calibration fixture")
        .unwrap();
    let stale = corpus()
        .query(
            &partitions,
            CHT,
            DatasetPurpose::Calibration,
            &cht_context(),
        )
        .unwrap();
    partitions
        .repartition(
            CHT,
            DatasetPartition::Validation,
            "freeze the calibration artifact before evaluation",
        )
        .unwrap();

    assert!(matches!(
        partitions.register_model(artifact("stale-calibration"), &[&stale], &[]),
        Err(PartitionRefusal::StaleAccess {
            receipt_generation: 1,
            current_generation: 2,
            ..
        })
    ));
}

#[test]
fn access_identity_binds_context_and_validation_retains_each_access() {
    let mut partitions = PartitionLedger::capture(corpus());
    partitions
        .repartition(CHT, DatasetPartition::Calibration, "calibration fixture")
        .unwrap();
    let calibration = corpus()
        .query(
            &partitions,
            CHT,
            DatasetPurpose::Calibration,
            &cht_context(),
        )
        .unwrap();
    let model = partitions
        .register_model(artifact("context-model"), &[&calibration], &[])
        .unwrap();
    let early = corpus()
        .query(
            &partitions,
            MARTIN_MOYCE,
            DatasetPurpose::Validation,
            &martin_moyce_context_at(1.0),
        )
        .unwrap();
    let late = corpus()
        .query(
            &partitions,
            MARTIN_MOYCE,
            DatasetPurpose::Validation,
            &martin_moyce_context_at(2.0),
        )
        .unwrap();

    assert_ne!(early.receipt().context(), late.receipt().context());
    assert_ne!(early.receipt().identity(), late.receipt().identity());
    let receipt = partitions.validate_model(&model, &[&early, &late]).unwrap();
    assert_eq!(receipt.evaluation_accesses().len(), 2);
}

#[test]
fn blind_evaluation_requires_a_generation_bound_release() {
    let mut partitions = PartitionLedger::capture(corpus());
    partitions
        .repartition(
            MARTIN_MOYCE,
            DatasetPartition::BlindHoldout,
            "seal this exact dataset generation for a blind drill",
        )
        .unwrap();
    assert!(matches!(
        corpus().query(
            &partitions,
            MARTIN_MOYCE,
            DatasetPurpose::BlindEvaluation,
            &martin_moyce_context(),
        ),
        Err(PartitionRefusal::BlindReleaseRequired {
            dataset_id,
            generation: 1,
        }) if dataset_id == MARTIN_MOYCE
    ));

    let release = partitions
        .release_blind(
            MARTIN_MOYCE,
            artifact("preregistration"),
            artifact("blind-manifest"),
            "release after the frozen model and protocol were committed",
        )
        .expect("non-zero exact blind identities");
    let access = corpus()
        .query(
            &partitions,
            MARTIN_MOYCE,
            DatasetPurpose::BlindEvaluation,
            &martin_moyce_context(),
        )
        .expect("released blind access");
    assert_eq!(access.receipt().blind_release(), Some(release.identity()));
    assert_eq!(partitions.blind_releases(), &[release]);
}

#[test]
fn taint_and_repartition_identities_ignore_caller_order() {
    let mut left = PartitionLedger::capture(corpus());
    let mut right = PartitionLedger::capture(corpus());
    for ledger in [&mut left, &mut right] {
        ledger
            .repartition(CHT, DatasetPartition::Calibration, "same transition")
            .unwrap();
        ledger
            .repartition(
                MARTIN_MOYCE,
                DatasetPartition::Calibration,
                "same second transition",
            )
            .unwrap();
    }
    assert_eq!(left.events(), right.events());

    let cht = corpus()
        .query(&left, CHT, DatasetPurpose::Calibration, &cht_context())
        .unwrap();
    let martin = corpus()
        .query(
            &left,
            MARTIN_MOYCE,
            DatasetPurpose::Calibration,
            &martin_moyce_context(),
        )
        .unwrap();
    let model_id = artifact("order-independent-model");
    let forward = left
        .register_model(model_id, &[&cht, &martin], &[])
        .unwrap();
    let reverse = left
        .register_model(model_id, &[&martin, &cht], &[])
        .unwrap();
    assert_eq!(forward, reverse);
    assert_eq!(forward.sources().len(), 2);
}

#[test]
fn validation_access_cannot_be_smuggled_into_model_training() {
    let partitions = PartitionLedger::capture(corpus());
    let validation = corpus()
        .query(&partitions, CHT, DatasetPurpose::Validation, &cht_context())
        .unwrap();
    assert!(matches!(
        partitions.register_model(artifact("smuggled"), &[&validation], &[]),
        Err(PartitionRefusal::WrongModelInputPurpose {
            dataset_id,
            purpose: DatasetPurpose::Validation,
        }) if dataset_id == CHT
    ));
}

#[test]
fn every_partition_receipt_record_round_trips_and_refuses_wire_tampering() {
    let records = receipt_records();
    assert_eq!(
        records
            .iter()
            .map(PartitionReceiptRecord::record_type)
            .collect::<Vec<_>>(),
        [
            "dataset-access",
            "repartition",
            "blind-release",
            "model-taint",
            "validation",
        ]
    );
    for record in &records {
        let encoded = record.encode();
        assert!(encoded.len() <= MAX_PARTITION_RECORD_BYTES);
        let decoded =
            PartitionReceiptRecord::decode(&encoded).expect("canonical record round trip");
        assert_eq!(decoded, *record);
        assert_eq!(decoded.semantic_identity(), record.semantic_identity());
    }

    let mut future = records[0].encode();
    future[8..12].copy_from_slice(&2_u32.to_le_bytes());
    assert_eq!(
        PartitionReceiptRecord::decode(&future),
        Err(PartitionRecordError::UnsupportedVersion { observed: 2 })
    );

    let mut unknown = records[0].encode();
    unknown[12] = 0xff;
    assert_eq!(
        PartitionReceiptRecord::decode(&unknown),
        Err(PartitionRecordError::UnknownVariant { observed: 0xff })
    );

    let mut tampered = records[0].encode();
    let last = tampered.len() - 1;
    tampered[last] ^= 1;
    assert!(matches!(
        PartitionReceiptRecord::decode(&tampered),
        Err(PartitionRecordError::IdentityMismatch {
            record: "dataset-access"
        })
    ));

    let mut truncated = records[4].encode();
    truncated.pop();
    assert!(matches!(
        PartitionReceiptRecord::decode(&truncated),
        Err(PartitionRecordError::Truncated {
            field: "validation identity"
        })
    ));

    let mut extended = records[2].encode();
    extended.push(0);
    assert_eq!(
        PartitionReceiptRecord::decode(&extended),
        Err(PartitionRecordError::TrailingBytes { observed: 1 })
    );
}

#[test]
fn partition_receipt_bytes_survive_ledger_reopen_dedupe_and_corruption_checks() {
    let records = receipt_records();
    let encoded = records
        .iter()
        .map(PartitionReceiptRecord::encode)
        .collect::<Vec<_>>();
    let path = partition_ledger_path();

    let hashes = {
        let ledger = Ledger::open(&path).expect("open durable partition ledger");
        encoded
            .iter()
            .map(|bytes| {
                let first = ledger
                    .put_artifact(
                        PARTITION_RECEIPT_ARTIFACT_KIND,
                        bytes,
                        Some(PARTITION_RECORD_META),
                    )
                    .expect("persist partition receipt bytes");
                assert!(!first.deduped);
                let duplicate = ledger
                    .put_artifact(
                        PARTITION_RECEIPT_ARTIFACT_KIND,
                        bytes,
                        Some(PARTITION_RECORD_META),
                    )
                    .expect("dedupe exact partition receipt bytes");
                assert!(duplicate.deduped);
                assert_eq!(duplicate.hash, first.hash);
                assert_eq!(duplicate.len, first.len);
                first.hash
            })
            .collect::<Vec<_>>()
    };

    let ledger = Ledger::open(&path).expect("reopen durable partition ledger");
    let read_limit =
        u64::try_from(MAX_PARTITION_RECORD_BYTES).expect("record byte cap fits in u64");
    for ((expected, bytes), hash) in records.iter().zip(&encoded).zip(&hashes) {
        let persisted = ledger
            .get_artifact_bounded(hash, read_limit)
            .expect("bounded durable read")
            .expect("persisted receipt is present");
        assert_eq!(persisted, *bytes);
        assert_eq!(
            PartitionReceiptRecord::decode(&persisted).expect("reverify persisted record"),
            *expected
        );
    }
    let intact = ledger
        .verify_artifact_integrity()
        .expect("scan intact receipt artifacts");
    assert_eq!(
        intact.checked,
        u64::try_from(records.len()).expect("record count fits in u64")
    );
    assert!(intact.corrupted.is_empty());

    ledger
        .corrupt_artifact_for_test(&hashes[0])
        .expect("inject receipt corruption");
    assert!(matches!(
        ledger.get_artifact_bounded(&hashes[0], read_limit),
        Err(LedgerError::Corrupt { .. })
    ));
    let corrupted = ledger
        .verify_artifact_integrity()
        .expect("scan corrupted receipt artifacts");
    assert_eq!(corrupted.corrupted, vec![hashes[0].to_hex()]);
}
