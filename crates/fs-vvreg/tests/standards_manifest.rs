//! G0/G3/G5 conformance for the exact standards source manifest.
//!
//! These fixtures exercise fail-closed edition selection, supersession,
//! access and hash gates, canonical migration framing, deterministic identity,
//! and the metadata-only protected-text boundary.

use fs_vvreg::ContentHash;
use fs_vvreg::standards::{
    ClauseLocator, EditionChange, EditionChangeKind, JurisdictionProfile, MAX_STANDARD_RECORDS,
    ManifestError, ManifestLimits, ProtectedTextReference, ReferenceState, RuleBindingError,
    RuleBindingRequest, RuleUsePolicy, STANDARD_MANIFEST_SCHEMA_VERSION, SourceAccess, SourcePin,
    StandardEditionKey, StandardLicense, StandardManifest, StandardSourceDraft, StandardStatus,
};

const fn hash(byte: u8) -> ContentHash {
    ContentHash([byte; 32])
}

fn key(edition: &str) -> StandardEditionKey {
    StandardEditionKey::try_new("iso-12100", Some("1"), edition).expect("fixture key must be valid")
}

fn draft(edition: &str, status: StandardStatus) -> StandardSourceDraft {
    StandardSourceDraft {
        key: key(edition),
        title: "Safety of machinery metadata fixture".to_string(),
        changes: vec![
            EditionChange::try_new(EditionChangeKind::Amendment, "Amd 1:2024")
                .expect("fixture amendment must be valid"),
            EditionChange::try_new(EditionChangeKind::Corrigendum, "Cor 1:2025")
                .expect("fixture corrigendum must be valid"),
        ],
        status,
        jurisdiction: JurisdictionProfile::try_new("global", "machinery-safety")
            .expect("fixture jurisdiction must be valid"),
        source: ProtectedTextReference::try_new(
            "publisher-catalog",
            "document/iso-12100/part-1/edition",
        )
        .expect("fixture source reference must be valid"),
        source_pin: SourcePin::Pinned(hash(0x11)),
        license: StandardLicense::Restricted {
            terms: "publisher-subscription-v1".to_string(),
        },
        access: SourceAccess::Available,
        no_claim: "Metadata and lineage only; no conformity claim".to_string(),
    }
}

fn manifest(drafts: Vec<StandardSourceDraft>) -> StandardManifest {
    StandardManifest::try_new(drafts, ManifestLimits::DEFAULT)
        .expect("fixture manifest must be admitted")
}

fn rule_request(source: StandardEditionKey) -> RuleBindingRequest {
    RuleBindingRequest {
        rule_id: "minimum-clearance-rule".to_string(),
        source,
        clause: ClauseLocator::try_new("Clause 5.4, Table 2")
            .expect("fixture clause must be valid"),
        observed_source_hash: hash(0x11),
        derived_rule_hash: hash(0x22),
        reference_state: ReferenceState::Derived,
        use_policy: RuleUsePolicy::CurrentOnly,
    }
}

#[test]
fn g0_exact_edition_and_source_gates_fail_closed() {
    let current_key = key("2010");
    let current = manifest(vec![draft("2010", StandardStatus::Current)]);

    let provenance = current
        .bind_rule(rule_request(current_key.clone()))
        .expect("exact pinned current source must bind");
    assert_eq!(provenance.source(), &current_key);
    assert_eq!(provenance.source_hash(), hash(0x11));
    assert_eq!(provenance.derived_rule_hash(), hash(0x22));
    assert_eq!(provenance.reference_state(), ReferenceState::Derived);
    assert!(!provenance.historical());
    assert_eq!(provenance.manifest_identity(), current.identity());
    assert_eq!(
        provenance.source_identity(),
        current.records()[0].identity()
    );

    let wrong_key = key("2010-amd-unknown");
    assert_eq!(
        current.bind_rule(rule_request(wrong_key.clone())),
        Err(RuleBindingError::UnknownEdition { key: wrong_key })
    );

    let mut unread = rule_request(current_key.clone());
    unread.reference_state = ReferenceState::Unread;
    assert_eq!(
        current.bind_rule(unread),
        Err(RuleBindingError::UnreadSource)
    );

    let mut mismatch = rule_request(current_key.clone());
    mismatch.observed_source_hash = hash(0x33);
    assert_eq!(
        current.bind_rule(mismatch),
        Err(RuleBindingError::SourceHashMismatch {
            key: current_key.clone(),
            expected: hash(0x11),
            observed: hash(0x33),
        })
    );

    let mut unpinned_draft = draft("2010", StandardStatus::Current);
    unpinned_draft.source_pin = SourcePin::Unpinned;
    let unpinned = manifest(vec![unpinned_draft]);
    assert_eq!(
        unpinned.bind_rule(rule_request(current_key.clone())),
        Err(RuleBindingError::UnpinnedSource {
            key: current_key.clone(),
        })
    );

    let mut revoked_draft = draft("2010", StandardStatus::Current);
    revoked_draft.access = SourceAccess::Revoked {
        reason_code: "license-expired".to_string(),
    };
    let revoked = manifest(vec![revoked_draft]);
    assert_eq!(
        revoked.bind_rule(rule_request(current_key.clone())),
        Err(RuleBindingError::AccessRevoked { key: current_key })
    );
}

#[test]
fn g0_historical_editions_require_explicit_replay_policy() {
    let old_key = key("2010");
    let new_key = key("2026");
    let standards = manifest(vec![
        draft(
            "2010",
            StandardStatus::Superseded {
                by: new_key.clone(),
            },
        ),
        draft("2026", StandardStatus::Current),
    ]);

    assert_eq!(
        standards.bind_rule(rule_request(old_key.clone())),
        Err(RuleBindingError::HistoricalEdition {
            key: old_key.clone(),
        })
    );

    let mut replay = rule_request(old_key);
    replay.use_policy = RuleUsePolicy::HistoricalReplay;
    let provenance = standards
        .bind_rule(replay)
        .expect("historical replay may bind the exact superseded source");
    assert!(provenance.historical());
    assert_eq!(provenance.use_policy(), RuleUsePolicy::HistoricalReplay);

    let current = standards
        .bind_rule(rule_request(new_key))
        .expect("declared successor remains current");
    assert!(!current.historical());

    let withdrawn = manifest(vec![draft(
        "2003",
        StandardStatus::Withdrawn {
            reason_code: "publisher-withdrawn".to_string(),
        },
    )]);
    let withdrawn_key = key("2003");
    assert_eq!(
        withdrawn.bind_rule(rule_request(withdrawn_key.clone())),
        Err(RuleBindingError::HistoricalEdition {
            key: withdrawn_key.clone(),
        })
    );
    let mut withdrawn_replay = rule_request(withdrawn_key);
    withdrawn_replay.use_policy = RuleUsePolicy::HistoricalReplay;
    assert!(
        withdrawn
            .bind_rule(withdrawn_replay)
            .expect("withdrawn exact source may bind for replay only")
            .historical()
    );
}

#[test]
fn g0_edition_collisions_and_malformed_change_chains_refuse() {
    let collision = StandardManifest::try_new(
        vec![
            draft("2010", StandardStatus::Current),
            draft("2010", StandardStatus::Current),
        ],
        ManifestLimits::DEFAULT,
    );
    assert_eq!(
        collision,
        Err(ManifestError::EditionCollision { key: key("2010") })
    );

    let mut duplicate_change = draft("2010", StandardStatus::Current);
    duplicate_change.changes.push(
        EditionChange::try_new(EditionChangeKind::Amendment, "Amd 1:2024")
            .expect("fixture amendment must be valid"),
    );
    assert_eq!(
        StandardManifest::try_new(vec![duplicate_change], ManifestLimits::DEFAULT),
        Err(ManifestError::DuplicateChange {
            key: key("2010"),
            designation: "Amd 1:2024".to_string(),
        })
    );
}

#[test]
fn g0_supersession_graph_requires_closed_acyclic_exact_keys() {
    let old_key = key("2010");
    let unknown_key = key("2030");

    let self_edge = draft(
        "2010",
        StandardStatus::Superseded {
            by: old_key.clone(),
        },
    );
    assert_eq!(
        StandardManifest::try_new(vec![self_edge], ManifestLimits::DEFAULT),
        Err(ManifestError::SelfSupersession {
            key: old_key.clone(),
        })
    );

    let unknown = draft(
        "2010",
        StandardStatus::Superseded {
            by: unknown_key.clone(),
        },
    );
    assert_eq!(
        StandardManifest::try_new(vec![unknown], ManifestLimits::DEFAULT),
        Err(ManifestError::UnknownSupersessionTarget {
            key: old_key.clone(),
            target: unknown_key,
        })
    );

    let newer_key = key("2026");
    let cycle = StandardManifest::try_new(
        vec![
            draft(
                "2010",
                StandardStatus::Superseded {
                    by: newer_key.clone(),
                },
            ),
            draft(
                "2026",
                StandardStatus::Superseded {
                    by: old_key.clone(),
                },
            ),
        ],
        ManifestLimits::DEFAULT,
    );
    assert_eq!(
        cycle,
        Err(ManifestError::SupersessionCycle { key: old_key })
    );
}

#[test]
fn g3_manifest_identity_is_order_invariant_but_change_order_is_semantic() {
    let older = draft(
        "2010",
        StandardStatus::Withdrawn {
            reason_code: "publisher-withdrawn".to_string(),
        },
    );
    let newer = draft("2026", StandardStatus::Current);
    let forward = manifest(vec![older.clone(), newer.clone()]);
    let reverse = manifest(vec![newer, older.clone()]);
    assert_eq!(forward.canonical_bytes(), reverse.canonical_bytes());
    assert_eq!(forward.identity(), reverse.identity());

    let mut reordered_changes = older;
    reordered_changes.changes.reverse();
    let reordered = manifest(vec![reordered_changes]);
    let original = manifest(vec![draft(
        "2010",
        StandardStatus::Withdrawn {
            reason_code: "publisher-withdrawn".to_string(),
        },
    )]);
    assert_ne!(reordered.identity(), original.identity());
    assert_ne!(
        reordered.records()[0].identity(),
        original.records()[0].identity()
    );
}

#[test]
fn g5_every_source_field_moves_the_sealed_identity() {
    let baseline_draft = draft("2010", StandardStatus::Current);
    let baseline = manifest(vec![baseline_draft.clone()]);
    let baseline_row = baseline.records()[0].identity();
    let baseline_manifest = baseline.identity();

    let assert_moved = |changed: StandardSourceDraft| {
        let changed = manifest(vec![changed]);
        assert_ne!(changed.records()[0].identity(), baseline_row);
        assert_ne!(changed.identity(), baseline_manifest);
    };

    let mut changed = baseline_draft.clone();
    changed.key = key("2011");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.key = StandardEditionKey::try_new("iso-13849", Some("1"), "2010")
        .expect("changed standard id must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.key = StandardEditionKey::try_new("iso-12100", Some("2"), "2010")
        .expect("changed part must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.title.push_str(" revised");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.changes[0] = EditionChange::try_new(EditionChangeKind::Amendment, "Amd 2:2024")
        .expect("changed amendment must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.changes[0] = EditionChange::try_new(EditionChangeKind::Corrigendum, "Amd 1:2024")
        .expect("changed change kind must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.status = StandardStatus::Withdrawn {
        reason_code: "publisher-withdrawn".to_string(),
    };
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.jurisdiction = JurisdictionProfile::try_new("us", "machinery-safety")
        .expect("changed jurisdiction must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.jurisdiction =
        JurisdictionProfile::try_new("global", "robotics").expect("changed profile must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.source = ProtectedTextReference::try_new("publisher-catalog", "document/revision")
        .expect("changed source must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.source = ProtectedTextReference::try_new(
        "alternate-publisher-catalog",
        "document/iso-12100/part-1/edition",
    )
    .expect("changed catalog must be valid");
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.source_pin = SourcePin::Pinned(hash(0x12));
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.license = StandardLicense::UserSupplied {
        policy: "authorized-local-copy-v2".to_string(),
    };
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.license = StandardLicense::Restricted {
        terms: "publisher-subscription-v2".to_string(),
    };
    assert_moved(changed);

    let mut changed = baseline_draft.clone();
    changed.access = SourceAccess::Revoked {
        reason_code: "license-expired".to_string(),
    };
    assert_moved(changed);

    let mut changed = baseline_draft;
    changed.no_claim.push_str("; revised boundary");
    assert_moved(changed);
}

#[test]
fn g5_rule_identity_covers_read_state_clause_and_derived_bytes() {
    let standards = manifest(vec![draft("2010", StandardStatus::Current)]);
    let baseline_request = rule_request(key("2010"));
    let baseline = standards
        .bind_rule(baseline_request.clone())
        .expect("baseline request must bind");
    assert_eq!(baseline.rule_id(), "minimum-clearance-rule");
    assert_eq!(baseline.clause().as_str(), "Clause 5.4, Table 2");
    assert_eq!(baseline.use_policy(), RuleUsePolicy::CurrentOnly);

    let mut changed = baseline_request.clone();
    changed.reference_state = ReferenceState::Read;
    let read = standards.bind_rule(changed).expect("read source must bind");
    assert_ne!(read.identity(), baseline.identity());

    let mut changed = baseline_request.clone();
    changed.reference_state = ReferenceState::Reproduced;
    let reproduced = standards
        .bind_rule(changed)
        .expect("reproduced source must bind");
    assert_ne!(reproduced.identity(), baseline.identity());

    let mut changed = baseline_request.clone();
    changed.clause = ClauseLocator::try_new("Clause 5.5").expect("changed clause must be valid");
    let clause = standards
        .bind_rule(changed)
        .expect("changed clause must bind");
    assert_ne!(clause.identity(), baseline.identity());

    let mut changed = baseline_request.clone();
    changed.rule_id = "alternate-clearance-rule".to_string();
    let rule = standards
        .bind_rule(changed)
        .expect("changed rule id must bind");
    assert_ne!(rule.identity(), baseline.identity());

    let mut changed = baseline_request.clone();
    changed.use_policy = RuleUsePolicy::HistoricalReplay;
    let replay_policy = standards
        .bind_rule(changed)
        .expect("historical policy may still bind a current source");
    assert_eq!(replay_policy.use_policy(), RuleUsePolicy::HistoricalReplay);
    assert_ne!(replay_policy.identity(), baseline.identity());

    let mut changed = baseline_request;
    changed.derived_rule_hash = hash(0x23);
    let derived = standards
        .bind_rule(changed)
        .expect("changed derived bytes must bind");
    assert_ne!(derived.identity(), baseline.identity());

    let mut changed_source = draft("2010", StandardStatus::Current);
    changed_source.title.push_str(" revised metadata");
    let changed_manifest = manifest(vec![changed_source]);
    let rebound = changed_manifest
        .bind_rule(rule_request(key("2010")))
        .expect("same exact source key in changed manifest must bind");
    assert_ne!(rebound.manifest_identity(), baseline.manifest_identity());
    assert_ne!(rebound.source_identity(), baseline.source_identity());
    assert_ne!(rebound.identity(), baseline.identity());
}

#[test]
fn g0_zero_hashes_refuse_before_provenance_can_be_minted() {
    let mut zero_source = draft("2010", StandardStatus::Current);
    zero_source.source_pin = SourcePin::Pinned(ContentHash([0; 32]));
    assert_eq!(
        StandardManifest::try_new(vec![zero_source], ManifestLimits::DEFAULT),
        Err(ManifestError::ZeroHash {
            field: "source_hash",
        })
    );

    let standards = manifest(vec![draft("2010", StandardStatus::Current)]);
    let mut zero_rule = rule_request(key("2010"));
    zero_rule.derived_rule_hash = ContentHash([0; 32]);
    assert_eq!(
        standards.bind_rule(zero_rule),
        Err(RuleBindingError::InvalidRequest(ManifestError::ZeroHash {
            field: "derived_rule_hash",
        }))
    );
}

#[test]
fn g5_canonical_round_trip_and_migration_header_are_frozen() {
    let empty = manifest(Vec::new());
    assert_eq!(
        empty.canonical_bytes(),
        &[b'F', b'S', b'M', b'F', 1, 0, 0, 0, 0, 0, 0, 0]
    );

    let standards = manifest(vec![draft("2010", StandardStatus::Current)]);
    let bytes = standards.canonical_bytes();
    assert_eq!(&bytes[..4], b"FSMF");
    assert_eq!(
        &bytes[4..8],
        &STANDARD_MANIFEST_SCHEMA_VERSION.to_le_bytes()
    );
    assert_eq!(&bytes[8..12], &1_u32.to_le_bytes());

    let decoded = StandardManifest::decode_canonical(bytes, ManifestLimits::DEFAULT)
        .expect("canonical v1 bytes must round-trip");
    assert_eq!(decoded, standards);

    let mut future = bytes.to_vec();
    future[4..8].copy_from_slice(&2_u32.to_le_bytes());
    assert_eq!(
        StandardManifest::decode_canonical(&future, ManifestLimits::DEFAULT),
        Err(ManifestError::UnknownSchema { version: 2 })
    );

    let truncated = &bytes[..bytes.len() - 1];
    assert!(matches!(
        StandardManifest::decode_canonical(truncated, ManifestLimits::DEFAULT),
        Err(ManifestError::UnexpectedEof { .. })
    ));

    let mut trailing = bytes.to_vec();
    trailing.push(0);
    assert_eq!(
        StandardManifest::decode_canonical(&trailing, ManifestLimits::DEFAULT),
        Err(ManifestError::TrailingBytes { remaining: 1 })
    );
}

#[test]
fn g5_semantically_valid_unsorted_wire_rows_are_noncanonical() {
    let older = manifest(vec![draft("2010", StandardStatus::Current)]);
    let newer = manifest(vec![draft("2026", StandardStatus::Current)]);
    let mut unsorted = Vec::new();
    unsorted.extend_from_slice(b"FSMF");
    unsorted.extend_from_slice(&STANDARD_MANIFEST_SCHEMA_VERSION.to_le_bytes());
    unsorted.extend_from_slice(&2_u32.to_le_bytes());
    unsorted.extend_from_slice(&newer.canonical_bytes()[12..]);
    unsorted.extend_from_slice(&older.canonical_bytes()[12..]);

    assert_eq!(
        StandardManifest::decode_canonical(&unsorted, ManifestLimits::DEFAULT),
        Err(ManifestError::NonCanonicalEncoding)
    );
}

#[test]
fn g0_explicit_resource_caps_refuse_before_unbounded_work() {
    let relaxed_hard_cap = ManifestLimits {
        records: MAX_STANDARD_RECORDS + 1,
        ..ManifestLimits::DEFAULT
    };
    assert_eq!(
        StandardManifest::try_new(Vec::new(), relaxed_hard_cap),
        Err(ManifestError::LimitExceedsHardMaximum {
            field: "records",
            offered: MAX_STANDARD_RECORDS + 1,
            maximum: MAX_STANDARD_RECORDS,
        })
    );

    let record_limited = ManifestLimits {
        records: 1,
        ..ManifestLimits::DEFAULT
    };
    assert_eq!(
        StandardManifest::try_new(
            vec![
                draft("2010", StandardStatus::Current),
                draft("2026", StandardStatus::Current),
            ],
            record_limited,
        ),
        Err(ManifestError::RecordLimit { count: 2, limit: 1 })
    );

    let change_limited = ManifestLimits {
        changes_per_record: 1,
        ..ManifestLimits::DEFAULT
    };
    assert_eq!(
        StandardManifest::try_new(vec![draft("2010", StandardStatus::Current)], change_limited,),
        Err(ManifestError::ChangeLimit {
            key: key("2010"),
            count: 2,
            limit: 1,
        })
    );

    let byte_limited = ManifestLimits {
        manifest_bytes: 12,
        ..ManifestLimits::DEFAULT
    };
    assert!(matches!(
        StandardManifest::try_new(vec![draft("2010", StandardStatus::Current)], byte_limited,),
        Err(ManifestError::ManifestByteLimit { limit: 12, .. })
    ));
}

#[test]
fn g3_redaction_retains_lineage_escapes_json_and_leaks_no_protected_text() {
    const PROTECTED_PARAGRAPH: &str =
        "COPYRIGHTED NORMATIVE PARAGRAPH MUST NEVER ENTER THE MANIFEST";

    let mut source = draft("2010", StandardStatus::Current);
    source.title = "Public bibliographic title".to_string();
    source.source =
        ProtectedTextReference::try_new("publisher-\"catalog\"", "document\\edition/2010")
            .expect("JSON-hostile metadata remains a valid external locator");
    source.license = StandardLicense::Restricted {
        terms: "confidential-contract-code".to_string(),
    };
    source.access = SourceAccess::Revoked {
        reason_code: "license-expired".to_string(),
    };
    source.no_claim = "internal no-claim prose".to_string();
    let standards = manifest(vec![source]);

    let canonical = String::from_utf8_lossy(standards.canonical_bytes());
    let rows = standards.redacted_rows();
    let row = &rows[0];
    assert!(!canonical.contains(PROTECTED_PARAGRAPH));
    assert!(!row.contains(PROTECTED_PARAGRAPH));
    assert!(!row.contains("confidential-contract-code"));
    assert!(!row.contains("license-expired"));
    assert!(!row.contains("internal no-claim prose"));
    assert!(!row.contains("Public bibliographic title"));
    assert!(row.contains("iso-12100"));
    assert!(row.contains("2010"));
    assert!(row.contains(&hash(0x11).to_hex()));
    assert!(row.contains(&standards.records()[0].identity().to_hex()));
    assert!(row.contains("publisher-\\\"catalog\\\""));
    assert!(row.contains("document\\\\edition/2010"));
    assert!(row.contains("\"access\":\"revoked\""));
}
