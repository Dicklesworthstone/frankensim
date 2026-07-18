use fs_blake3::identity::{
    CanonicalEncoder, CanonicalLimits, CanonicalSchema, ContentId, Field, FieldSpec,
    IdentityReceipt, SemanticId, TrustState, WireType, legacy::LegacyProvenanceV1,
};
use fs_ledger::{
    ARTIFACT_CONTENT_IDENTITY_ROW_VERSION, EDGE_CONTENT_IDENTITY_ROW_VERSION, EdgeRole,
    FiveExplicits, IdentityMigrationClaim, Ledger, LedgerError,
    MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,
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
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
fn typed_projection_refuses_a_different_schema() {
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
fn explicit_receipt_binding_projects_only_the_exact_nominal_schema_and_roots_gc() {
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
    let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
