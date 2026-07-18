use std::mem::size_of;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalLimits, CanonicalSchema, ContentId, Field, FieldSpec,
    IdentityReceipt, SemanticId, TrustState, WireType, legacy::LegacyProvenanceV1,
};
use fs_ledger::{
    ARTIFACT_CONTENT_IDENTITY_ROW_VERSION, EDGE_CONTENT_IDENTITY_ROW_VERSION, EdgeRole,
    ExtensionTable, FiveExplicits, IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION, IdentityMigrationClaim,
    IdentityMigrationReceipt, IdentityMigrationWireError, Ledger, LedgerError,
    MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES, MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES,
    OP_CONTENT_IDENTITY_ROW_VERSION,
};

const LIMITS: CanonicalLimits = CanonicalLimits::new(64 * 1024, 16 * 1024, 8, 16, 4096);

enum DemoSemanticSchemaV1 {}

impl CanonicalSchema for DemoSemanticSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.test.ledger-migration.demo.v1";
    const NAME: &'static str = "ledger-migration-demo";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G0/G3 exact-byte identity migration fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("meaning", WireType::Bytes)];
}

type DemoSemanticId = SemanticId<DemoSemanticSchemaV1>;

enum OtherSemanticSchemaV1 {}

impl CanonicalSchema for OtherSemanticSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.test.ledger-migration.other.v1";
    const NAME: &'static str = "ledger-migration-other";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "G3 wrong-schema refusal fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("meaning", WireType::Bytes)];
}

type OtherSemanticId = SemanticId<OtherSemanticSchemaV1>;

fn semantic_receipt(meaning: &[u8]) -> IdentityReceipt<DemoSemanticId> {
    CanonicalEncoder::<DemoSemanticId, _>::new(LIMITS, || false)
        .expect("valid static migration fixture schema")
        .bytes(Field::new(0, "meaning"), meaning)
        .expect("bounded semantic fixture")
        .finish()
        .expect("complete semantic fixture")
}

fn claim<'a>(
    receipt: IdentityReceipt<DemoSemanticId>,
    legacy_bytes: &'a [u8],
    canonical_bytes: &'a [u8],
    semantic_rule: &'a str,
) -> IdentityMigrationClaim<'a, DemoSemanticId> {
    IdentityMigrationClaim {
        legacy_bytes,
        legacy_fnv: LegacyProvenanceV1::new(0xcbf2_9ce4_8422_2325),
        canonical_bytes,
        semantic_rule,
        receipt,
        audit: receipt.audit_record(),
    }
}

#[test]
fn receipt_identity_binds_exact_bytes_schema_and_audit_state() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let legacy = br#"{"legacy":"shape-a","provenance":1}"#;
    let canonical = br#"{"schema":1,"shape":"a"}"#;
    let semantic = semantic_receipt(b"shape-a");

    let first = ledger
        .record_identity_migration(claim(semantic, legacy, canonical, "demo-json-v0-to-v1"))
        .expect("record exact migration");
    assert!(!first.deduped());
    assert_eq!(first.legacy_content_id(), ContentId::of_bytes(legacy));
    assert_eq!(first.canonical_content_id(), ContentId::of_bytes(canonical));

    let stored = ledger
        .identity_migration_receipt(first.receipt_id())
        .expect("reverify stored receipt")
        .expect("stored receipt exists");
    assert_eq!(stored.legacy_bytes(), legacy);
    assert_eq!(stored.canonical_bytes(), canonical);
    assert_eq!(stored.legacy_fnv().value(), 0xcbf2_9ce4_8422_2325);
    assert_eq!(stored.semantic_rule(), "demo-json-v0-to-v1");
    assert_eq!(stored.trust_state(), TrustState::Unanchored);
    assert_eq!(
        stored.typed_semantic_id::<DemoSemanticId>(),
        Some(semantic.id())
    );
    let wire = stored.to_wire_bytes().expect("encode complete v1 receipt");
    assert!(wire.len() <= MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES);
    assert_eq!(
        u32::from_le_bytes(wire[8..12].try_into().unwrap()),
        IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION
    );
    let transported = IdentityMigrationReceipt::from_wire_bytes(&wire)
        .expect("decode and independently reconstruct complete v1 receipt");
    assert_eq!(transported, stored);
    assert_eq!(transported.legacy_fnv().value(), 0xcbf2_9ce4_8422_2325);
    assert_eq!(
        transported.typed_semantic_id::<DemoSemanticId>(),
        Some(semantic.id())
    );

    let retry = ledger
        .record_identity_migration(claim(semantic, legacy, canonical, "demo-json-v0-to-v1"))
        .expect("exact response-loss retry");
    assert!(retry.deduped());
    assert_eq!(retry.receipt_id(), first.receipt_id());

    let changed_legacy = ledger
        .record_identity_migration(claim(
            semantic,
            b"different-legacy",
            canonical,
            "demo-json-v0-to-v1",
        ))
        .unwrap();
    let changed_canonical = ledger
        .record_identity_migration(claim(
            semantic,
            legacy,
            b"different-canonical",
            "demo-json-v0-to-v1",
        ))
        .unwrap();
    let changed_rule = ledger
        .record_identity_migration(claim(semantic, legacy, canonical, "different-rule"))
        .unwrap();
    let changed_semantic = ledger
        .record_identity_migration(claim(
            semantic_receipt(b"shape-b"),
            legacy,
            canonical,
            "demo-json-v0-to-v1",
        ))
        .unwrap();
    let mut fnv_claim = claim(semantic, legacy, canonical, "demo-json-v0-to-v1");
    fnv_claim.legacy_fnv = LegacyProvenanceV1::new(7);
    let changed_fnv = ledger.record_identity_migration(fnv_claim).unwrap();
    for changed in [
        changed_legacy,
        changed_canonical,
        changed_rule,
        changed_semantic,
        changed_fnv,
    ] {
        assert_ne!(changed.receipt_id(), first.receipt_id());
    }
}

#[test]
fn wire_transport_refuses_truncation_extension_future_version_and_forged_id() {
    const WIRE_VERSION_OFFSET: usize = 8;
    const RECEIPT_ID_OFFSET: usize = WIRE_VERSION_OFFSET + size_of::<u32>();
    const LEGACY_LENGTH_OFFSET: usize = RECEIPT_ID_OFFSET + 32;
    const LEGACY_BYTES_OFFSET: usize = LEGACY_LENGTH_OFFSET + size_of::<u32>();

    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let semantic = semantic_receipt(b"wire-subject");
    let write = ledger
        .record_identity_migration(claim(
            semantic,
            b"legacy-wire-payload",
            b"canonical-wire-payload",
            "demo-wire-v0-to-v1",
        ))
        .expect("record wire fixture");
    let stored = ledger
        .identity_migration_receipt(write.receipt_id())
        .expect("reverify stored wire fixture")
        .expect("stored wire fixture exists");
    let wire = stored
        .to_wire_bytes()
        .expect("encode complete wire fixture");

    let mut truncated = wire.clone();
    let _ = truncated.pop();
    assert!(matches!(
        IdentityMigrationReceipt::from_wire_bytes(&truncated),
        Err(IdentityMigrationWireError::Truncated { .. })
    ));

    let mut extended = wire.clone();
    extended.push(0);
    assert_eq!(
        IdentityMigrationReceipt::from_wire_bytes(&extended),
        Err(IdentityMigrationWireError::TrailingBytes { remaining: 1 })
    );

    let unsupported_version = IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION + 1;
    let mut future_version = wire.clone();
    future_version[WIRE_VERSION_OFFSET..RECEIPT_ID_OFFSET]
        .copy_from_slice(&unsupported_version.to_le_bytes());
    assert_eq!(
        IdentityMigrationReceipt::from_wire_bytes(&future_version),
        Err(IdentityMigrationWireError::UnsupportedVersion {
            found: unsupported_version,
        })
    );

    let mut forged_id = wire.clone();
    forged_id[RECEIPT_ID_OFFSET] ^= 0x80;
    assert!(matches!(
        IdentityMigrationReceipt::from_wire_bytes(&forged_id),
        Err(IdentityMigrationWireError::ReceiptIdMismatch { .. })
    ));

    let mut wrong_magic = wire.clone();
    wrong_magic[0] ^= 0x80;
    assert_eq!(
        IdentityMigrationReceipt::from_wire_bytes(&wrong_magic),
        Err(IdentityMigrationWireError::Magic)
    );

    let mut mismatched_content = wire.clone();
    mismatched_content[LEGACY_BYTES_OFFSET] ^= 0x40;
    assert!(matches!(
        IdentityMigrationReceipt::from_wire_bytes(&mismatched_content),
        Err(IdentityMigrationWireError::InvalidReceipt { .. })
    ));

    let mut oversized_legacy = wire;
    oversized_legacy[LEGACY_LENGTH_OFFSET..LEGACY_LENGTH_OFFSET + size_of::<u32>()]
        .copy_from_slice(&u32::MAX.to_le_bytes());
    assert_eq!(
        IdentityMigrationReceipt::from_wire_bytes(&oversized_legacy),
        Err(IdentityMigrationWireError::FieldLength {
            field: "legacy_bytes",
            found: usize::try_from(u32::MAX).unwrap(),
            max: MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,
        })
    );
}

#[test]
fn typed_projection_refuses_a_different_schema() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let semantic = semantic_receipt(b"typed-subject");
    let write = ledger
        .record_identity_migration(claim(semantic, b"legacy", b"canonical", "demo-v0-to-v1"))
        .expect("record typed migration");
    let stored = ledger
        .identity_migration_receipt(write.receipt_id())
        .unwrap()
        .unwrap();
    assert_eq!(stored.typed_semantic_id::<OtherSemanticId>(), None);
    assert_eq!(
        stored.typed_semantic_id::<DemoSemanticId>(),
        Some(semantic.id())
    );
}

#[test]
fn ambiguous_legacy_candidates_are_bounded_and_never_selected() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let legacy = b"same-legacy-source";
    let semantic = semantic_receipt(b"same-subject");
    let first = ledger
        .record_identity_migration(claim(semantic, legacy, b"canonical-a", "demo-rule-a"))
        .unwrap();
    let second = ledger
        .record_identity_migration(claim(semantic, legacy, b"canonical-b", "demo-rule-b"))
        .unwrap();
    assert_ne!(first.receipt_id(), second.receipt_id());

    let existence = ledger
        .identity_migration_candidates(ContentId::of_bytes(legacy), 0)
        .unwrap();
    assert!(existence.receipt_ids().is_empty());
    assert!(existence.truncated());

    let one = ledger
        .identity_migration_candidates(ContentId::of_bytes(legacy), 1)
        .unwrap();
    assert_eq!(one.receipt_ids().len(), 1);
    assert!(one.truncated());

    let all = ledger
        .identity_migration_candidates(ContentId::of_bytes(legacy), 2)
        .unwrap();
    assert_eq!(all.receipt_ids().len(), 2);
    assert!(!all.truncated());
}

#[test]
fn payload_limit_refuses_before_any_row_is_published() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let oversized = vec![0xA5; MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES + 1];
    let semantic = semantic_receipt(b"bounded-subject");
    assert!(matches!(
        ledger.record_identity_migration(claim(
            semantic,
            &oversized,
            b"canonical",
            "demo-v0-to-v1",
        )),
        Err(LedgerError::Invalid { .. })
    ));
    assert_eq!(
        ledger.table_count("identity_migration_receipts").unwrap(),
        0
    );
}

#[test]
fn artifact_writes_dual_write_an_exact_typed_content_identity() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let bytes = b"artifact identity dual-write fixture";
    let write = ledger
        .put_artifact("identity-fixture", bytes, None)
        .expect("store exact artifact");

    let identity = ledger
        .artifact_content_identity(&write.hash)
        .expect("verify artifact content identity")
        .expect("stored artifact has a sidecar");
    assert_eq!(identity.artifact_hash(), write.hash);
    assert_eq!(identity.content_id(), ContentId::of_bytes(bytes));
    assert_eq!(
        identity.row_schema_version(),
        ARTIFACT_CONTENT_IDENTITY_ROW_VERSION
    );
    assert_eq!(
        ledger.table_count("artifact_content_identities").unwrap(),
        1
    );

    let retry = ledger
        .put_artifact("identity-fixture", bytes, None)
        .expect("dedupe exact artifact");
    assert!(retry.deduped);
    assert_eq!(
        ledger.table_count("artifact_content_identities").unwrap(),
        1,
        "artifact dedupe must not duplicate typed identity rows"
    );
}

#[test]
fn lineage_edges_dual_write_the_linked_artifact_content_identity() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let artifact = ledger
        .put_artifact("edge-identity-fixture", b"lineage payload", None)
        .expect("store linked artifact");
    let explicits = FiveExplicits {
        seed: b"edge-seed",
        versions: "{}",
        budget: "{}",
        capability: "{}",
    };
    let op = ledger
        .begin_op(None, "{}", &explicits, 1)
        .expect("begin lineage operation");
    ledger
        .link(op, &artifact.hash, EdgeRole::Out)
        .expect("link typed artifact output");

    let identity = ledger
        .edge_content_identity(op, &artifact.hash, EdgeRole::Out)
        .expect("verify edge content identity")
        .expect("linked edge has a sidecar");
    assert_eq!(identity.op(), op);
    assert_eq!(identity.role(), EdgeRole::Out);
    assert_eq!(identity.artifact_hash(), artifact.hash);
    assert_eq!(
        identity.content_id(),
        ContentId::of_bytes(b"lineage payload")
    );
    assert_eq!(
        identity.row_schema_version(),
        EDGE_CONTENT_IDENTITY_ROW_VERSION
    );
    assert_eq!(ledger.table_count("edge_content_identities").unwrap(), 1);
    assert_eq!(
        ledger
            .edge_content_identity(op, &artifact.hash, EdgeRole::In)
            .unwrap(),
        None,
        "role remains a separate part of edge identity"
    );
}

#[test]
fn operation_fields_receive_separate_exact_typed_content_identities() {
    let ledger = Ledger::open(":memory:").expect("fresh v18 ledger");
    let session = b"operation-session";
    let ir = r#"{"op":"identity-fixture","units":"m"}"#;
    let explicits = FiveExplicits {
        seed: b"operation-seed",
        versions: r#"{"fs-ledger":"0.1.0"}"#,
        budget: r#"{"memory_bytes":4096}"#,
        capability: r#"{"cores":2}"#,
    };
    let op = ledger
        .begin_op(Some(session), ir, &explicits, 1)
        .expect("atomically record operation and typed content sidecar");

    let identity = ledger
        .op_content_identity(op)
        .expect("independently re-hash operation fields")
        .expect("operation identity sidecar exists");
    assert_eq!(identity.op(), op);
    assert_eq!(
        identity.session_content_id(),
        Some(ContentId::of_bytes(session))
    );
    assert_eq!(identity.ir_content_id(), ContentId::of_bytes(ir.as_bytes()));
    assert_eq!(
        identity.seed_content_id(),
        ContentId::of_bytes(explicits.seed)
    );
    assert_eq!(
        identity.versions_content_id(),
        ContentId::of_bytes(explicits.versions.as_bytes())
    );
    assert_eq!(
        identity.budget_content_id(),
        ContentId::of_bytes(explicits.budget.as_bytes())
    );
    assert_eq!(
        identity.capability_content_id(),
        ContentId::of_bytes(explicits.capability.as_bytes())
    );
    assert_eq!(
        identity.row_schema_version(),
        OP_CONTENT_IDENTITY_ROW_VERSION
    );

    ledger
        .finish_op(op, fs_ledger::OpOutcome::Ok, None, 2)
        .expect("terminal envelope remains outside frozen-field identities");
    assert_eq!(ledger.op_content_identity(op).unwrap(), Some(identity));
    assert_eq!(ledger.table_count("op_content_identities").unwrap(), 1);

    ledger.begin().expect("open caller-owned transaction");
    let rolled_back = ledger
        .begin_op(None, r#"{"op":"rollback-fixture"}"#, &explicits, 3)
        .expect("record operation and sidecar inside caller transaction");
    assert_eq!(
        ledger
            .op_content_identity(rolled_back)
            .unwrap()
            .expect("uncommitted sidecar is visible to its transaction")
            .session_content_id(),
        None
    );
    ledger.rollback().expect("roll back caller transaction");
    assert_eq!(ledger.op(rolled_back).unwrap(), None);
    assert_eq!(ledger.op_content_identity(rolled_back).unwrap(), None);
    assert_eq!(ledger.table_count("op_content_identities").unwrap(), 1);
}

#[test]
fn explicit_receipt_binding_projects_only_the_exact_nominal_schema_and_roots_gc() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let bytes = b"semantic artifact bytes";
    let artifact = ledger
        .put_artifact("semantic-fixture", bytes, None)
        .expect("store canonical artifact");
    let semantic = semantic_receipt(b"semantic-artifact");
    let migration = ledger
        .record_identity_migration(claim(
            semantic,
            b"legacy semantic artifact",
            bytes,
            "semantic-artifact-v0-to-v1",
        ))
        .expect("record exact semantic receipt");

    let first = ledger
        .bind_artifact_semantic_identity(migration.receipt_id())
        .expect("bind retained canonical artifact");
    assert!(!first.deduped());
    assert_eq!(first.artifact_hash(), artifact.hash);
    let retry = ledger
        .bind_artifact_semantic_identity(migration.receipt_id())
        .expect("dedupe exact artifact semantic binding");
    assert!(retry.deduped());

    let stored = ledger
        .artifact_semantic_binding(&artifact.hash, migration.receipt_id())
        .expect("reverify artifact semantic binding")
        .expect("binding exists");
    assert_eq!(stored.artifact_hash(), artifact.hash);
    assert_eq!(
        stored.typed_semantic_id::<DemoSemanticId>(),
        Some(semantic.id())
    );
    assert_eq!(stored.typed_semantic_id::<OtherSemanticId>(), None);

    let gc = ledger
        .gc_unreferenced_artifacts(false)
        .expect("semantic binding is a GC root");
    assert!(!gc.candidates.contains(&artifact.hash.to_hex()));
    assert!(ledger.get_artifact(&artifact.hash).unwrap().is_some());
}

#[test]
fn artifact_semantic_candidates_preserve_ambiguity_and_missing_artifacts_refuse() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let bytes = b"shared semantic artifact";
    let artifact = ledger
        .put_artifact("semantic-fixture", bytes, None)
        .expect("store shared artifact");
    let semantic = semantic_receipt(b"shared-meaning");
    let first = ledger
        .record_identity_migration(claim(semantic, b"legacy-a", bytes, "shared-rule-a"))
        .unwrap();
    let second = ledger
        .record_identity_migration(claim(semantic, b"legacy-b", bytes, "shared-rule-b"))
        .unwrap();
    ledger
        .bind_artifact_semantic_identity(first.receipt_id())
        .unwrap();
    ledger
        .bind_artifact_semantic_identity(second.receipt_id())
        .unwrap();

    let existence = ledger
        .artifact_semantic_binding_candidates(&artifact.hash, 0)
        .unwrap();
    assert!(existence.receipt_ids().is_empty());
    assert!(existence.truncated());
    let one = ledger
        .artifact_semantic_binding_candidates(&artifact.hash, 1)
        .unwrap();
    assert_eq!(one.receipt_ids().len(), 1);
    assert!(one.truncated());
    let both = ledger
        .artifact_semantic_binding_candidates(&artifact.hash, 2)
        .unwrap();
    assert_eq!(both.receipt_ids().len(), 2);
    assert!(!both.truncated());

    let absent = ledger
        .record_identity_migration(claim(
            semantic,
            b"legacy-missing",
            b"canonical bytes not retained as an artifact",
            "missing-artifact-rule",
        ))
        .unwrap();
    assert!(matches!(
        ledger.bind_artifact_semantic_identity(absent.receipt_id()),
        Err(LedgerError::NotFound { .. })
    ));
    assert_eq!(ledger.table_count("artifact_semantic_bindings").unwrap(), 2);
}

#[test]
fn exact_evidence_binding_freezes_bytes_and_preserves_receipt_ambiguity() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let evidence_name = "gauntlet/exact-evidence";
    let body = r#"{"claim":"exact","value":17}"#;
    ledger
        .put_extension(ExtensionTable::Evidence, evidence_name, body)
        .expect("store exact evidence JSON");
    let semantic = semantic_receipt(b"exact-evidence");
    let first_receipt = ledger
        .record_identity_migration(claim(
            semantic,
            b"legacy evidence v0",
            body.as_bytes(),
            "evidence-v0-to-v1",
        ))
        .expect("record exact evidence receipt");

    let first = ledger
        .bind_evidence_semantic_identity(evidence_name, first_receipt.receipt_id())
        .expect("bind exact evidence bytes");
    assert_eq!(first.evidence_name(), evidence_name);
    assert_eq!(first.content_id(), ContentId::of_bytes(body.as_bytes()));
    assert!(!first.deduped());
    let retry = ledger
        .bind_evidence_semantic_identity(evidence_name, first_receipt.receipt_id())
        .expect("dedupe exact evidence binding");
    assert!(retry.deduped());

    let stored = ledger
        .evidence_semantic_binding(evidence_name, first_receipt.receipt_id())
        .expect("reverify evidence binding")
        .expect("evidence binding exists");
    assert_eq!(stored.evidence_name(), evidence_name);
    assert_eq!(stored.receipt_id(), first_receipt.receipt_id());
    assert_eq!(stored.content_id(), ContentId::of_bytes(body.as_bytes()));
    assert_eq!(
        stored.typed_semantic_id::<DemoSemanticId>(),
        Some(semantic.id())
    );
    assert_eq!(stored.typed_semantic_id::<OtherSemanticId>(), None);

    ledger
        .put_extension(ExtensionTable::Evidence, evidence_name, body)
        .expect("exact response-loss source retry remains legal");
    assert!(
        ledger
            .put_extension(
                ExtensionTable::Evidence,
                evidence_name,
                r#"{"claim":"reinterpreted","value":17}"#,
            )
            .is_err(),
        "a bound evidence body must not be silently reinterpreted"
    );
    assert_eq!(
        ledger
            .get_extension(ExtensionTable::Evidence, evidence_name)
            .unwrap()
            .as_deref(),
        Some(body)
    );

    let second_receipt = ledger
        .record_identity_migration(claim(
            semantic,
            b"different legacy evidence",
            body.as_bytes(),
            "independent-evidence-rule",
        ))
        .expect("record a second exact receipt");
    ledger
        .bind_evidence_semantic_identity(evidence_name, second_receipt.receipt_id())
        .expect("preserve a second exact interpretation");
    let existence = ledger
        .evidence_semantic_binding_candidates(evidence_name, 0)
        .unwrap();
    assert!(existence.receipt_ids().is_empty());
    assert!(existence.truncated());
    let one = ledger
        .evidence_semantic_binding_candidates(evidence_name, 1)
        .unwrap();
    assert_eq!(one.receipt_ids().len(), 1);
    assert!(one.truncated());
    let both = ledger
        .evidence_semantic_binding_candidates(evidence_name, 2)
        .unwrap();
    assert_eq!(both.receipt_ids().len(), 2);
    assert!(!both.truncated());
}

#[test]
fn evidence_binding_refuses_json_reformatting_missing_rows_and_partial_publication() {
    let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
    let evidence_name = "gauntlet/byte-sensitive";
    let compact = r#"{"a":1}"#;
    ledger
        .put_extension(ExtensionTable::Evidence, evidence_name, compact)
        .expect("store compact evidence JSON");
    let semantic = semantic_receipt(b"byte-sensitive-evidence");
    let reformatted = ledger
        .record_identity_migration(claim(
            semantic,
            b"legacy reformatted evidence",
            br#"{"a": 1}"#,
            "json-reformat-is-not-identity",
        ))
        .expect("record independently valid reformatted receipt");
    assert!(matches!(
        ledger.bind_evidence_semantic_identity(evidence_name, reformatted.receipt_id()),
        Err(LedgerError::Invalid { .. })
    ));
    assert_eq!(ledger.table_count("evidence_semantic_bindings").unwrap(), 0);

    let missing = ledger
        .record_identity_migration(claim(
            semantic,
            b"legacy missing evidence",
            compact.as_bytes(),
            "missing-evidence-rule",
        ))
        .expect("record receipt whose evidence row is absent");
    assert!(matches!(
        ledger.bind_evidence_semantic_identity("gauntlet/missing", missing.receipt_id()),
        Err(LedgerError::NotFound { .. })
    ));
    assert_eq!(ledger.table_count("evidence_semantic_bindings").unwrap(), 0);
}
