//! G0/G3 battery for bounded concrete traceability-source loading.

use fs_blake3::hash_bytes;
use fs_govern::{
    LoadedTraceabilityGenerationError, TRACEABILITY_FILESYSTEM_ADAPTER_VERSION,
    TRACEABILITY_OWNER_JOIN_VERSION, TraceabilityField, TraceabilityFileSpec,
    TraceabilityFilesystemField, TraceabilityFilesystemLimits, TraceabilityOwnerJoinField,
    TraceabilitySourceKind, audit_traceability_owner_join,
    generate_traceability_ledger_from_loaded_sources, generate_traceability_ledger_from_snapshot,
    load_traceability_source_snapshot, proof_obligations, requirements,
};
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const BEADS_PATH: &str = ".beads/issues.jsonl";
const CONTRACT_PATH: &str = "crates/demo/CONTRACT.md";
const REGISTRY_PATH: &str = "registry/traceability.rs";

const BEADS_BYTES: &[u8] = br#"{"title":"first","id":"bead-a","meta":{"id":"nested"}}
{"id":"bead-b","deps":[{"id":"ignored"}]}
"#;
const CONTRACT_BYTES: &[u8] = b"# Demo contract\n\nDeclaration-only fixture.\n";
const REGISTRY_BYTES: &[u8] = b"pub const SCHEMA: &str = \"traceability-v1\";\n";

static NEXT_ROOT: AtomicU64 = AtomicU64::new(0);

fn unique_root(label: &str) -> PathBuf {
    let serial = NEXT_ROOT.fetch_add(1, Ordering::Relaxed);
    let root = std::env::temp_dir().join(format!(
        "fs-govern-traceability-{label}-{}-{serial}",
        std::process::id()
    ));
    std::fs::create_dir(&root).expect("fixture root is unique within the test process");
    root
}

fn write_source(root: &Path, relative: &str, bytes: &[u8]) {
    let path = root.join(relative);
    std::fs::create_dir_all(path.parent().expect("fixture source has a parent"))
        .expect("fixture parent directory");
    std::fs::write(path, bytes).expect("fixture source bytes");
}

fn write_fixture(root: &Path, beads: &[u8], contract: &[u8], registry: &[u8]) {
    write_source(root, BEADS_PATH, beads);
    write_source(root, CONTRACT_PATH, contract);
    write_source(root, REGISTRY_PATH, registry);
}

fn write_complete_fixture(root: &Path, beads: &[u8], registry: &[u8]) {
    write_fixture(root, beads, CONTRACT_BYTES, registry);
}

fn complete_specs() -> [TraceabilityFileSpec<'static>; 3] {
    [
        TraceabilityFileSpec {
            kind: TraceabilitySourceKind::Registry,
            relative_path: Path::new(REGISTRY_PATH),
        },
        TraceabilityFileSpec {
            kind: TraceabilitySourceKind::Beads,
            relative_path: Path::new(BEADS_PATH),
        },
        TraceabilityFileSpec {
            kind: TraceabilitySourceKind::Contract,
            relative_path: Path::new(CONTRACT_PATH),
        },
    ]
}

fn audit_with_beads(beads: &[u8]) -> fs_govern::TraceabilityFilesystemAudit {
    let root = unique_root("bad-beads");
    write_complete_fixture(&root, beads, REGISTRY_BYTES);
    load_traceability_source_snapshot(
        &root,
        &complete_specs(),
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("invalid Beads input must fail closed")
}

fn audit_with_text(
    label: &str,
    contract: &[u8],
    registry: &[u8],
) -> fs_govern::TraceabilityFilesystemAudit {
    let root = unique_root(label);
    write_fixture(&root, BEADS_BYTES, contract, registry);
    load_traceability_source_snapshot(
        &root,
        &complete_specs(),
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("invalid text source must fail closed")
}

fn owner_beads_json(status: &str, omitted: Option<&str>, extra: &[&str]) -> String {
    let mut ids = BTreeSet::new();
    for obligation in proof_obligations() {
        for owner in obligation.owner_beads {
            if Some(*owner) != omitted {
                ids.insert(*owner);
            }
        }
    }
    ids.extend(extra.iter().copied());
    let mut jsonl = String::new();
    for id in ids {
        writeln!(
            jsonl,
            "{{\"id\":\"{id}\",\"status\":\"{status}\",\"dependencies\":[]}}"
        )
        .expect("writing to a String is infallible");
    }
    jsonl
}

fn loaded_owner_fixture(
    label: &str,
    status: &str,
    omitted: Option<&str>,
    extra: &[&str],
) -> fs_govern::LoadedTraceabilitySourceSnapshot {
    let root = unique_root(label);
    let beads = owner_beads_json(status, omitted, extra);
    write_complete_fixture(&root, beads.as_bytes(), REGISTRY_BYTES);
    load_traceability_source_snapshot(
        &root,
        &complete_specs(),
        TraceabilityFilesystemLimits::default(),
    )
    .expect("owner fixture is a valid concrete source snapshot")
}

#[test]
fn concrete_sources_bind_exact_bytes_in_canonical_order() {
    let root = unique_root("canonical");
    write_complete_fixture(&root, BEADS_BYTES, REGISTRY_BYTES);

    let loaded = load_traceability_source_snapshot(
        &root,
        &complete_specs(),
        TraceabilityFilesystemLimits::default(),
    )
    .expect("complete concrete source snapshot");
    let receipt = loaded.receipt();
    assert_eq!(
        receipt.adapter_version(),
        TRACEABILITY_FILESYSTEM_ADAPTER_VERSION
    );
    assert_eq!(receipt.limits(), TraceabilityFilesystemLimits::default());
    assert_eq!(
        receipt.total_bytes(),
        BEADS_BYTES.len() + CONTRACT_BYTES.len() + REGISTRY_BYTES.len()
    );
    assert_eq!(receipt.sources().len(), 3);
    assert_eq!(receipt.sources()[0].locator(), BEADS_PATH);
    assert_eq!(receipt.sources()[1].locator(), CONTRACT_PATH);
    assert_eq!(receipt.sources()[2].locator(), REGISTRY_PATH);

    let beads_receipt = receipt
        .sources()
        .iter()
        .find(|source| source.kind() == TraceabilitySourceKind::Beads)
        .expect("Beads receipt");
    assert_eq!(beads_receipt.content_identity(), hash_bytes(BEADS_BYTES));
    assert_eq!(beads_receipt.byte_count(), BEADS_BYTES.len());
    assert_eq!(beads_receipt.beads_record_count(), Some(2));
    assert_eq!(
        beads_receipt.beads_ids(),
        &["bead-a".to_string(), "bead-b".to_string()]
    );
    assert!(
        receipt
            .sources()
            .iter()
            .filter(|source| source.kind() != TraceabilitySourceKind::Beads)
            .all(|source| source.beads_record_count().is_none())
    );
    assert_eq!(
        receipt.source_snapshot_identity(),
        loaded.snapshot().identity()
    );

    let artifact = generate_traceability_ledger_from_snapshot(
        requirements(),
        proof_obligations(),
        loaded.snapshot(),
    )
    .expect("canonical declarations bind to concrete source bytes");
    assert_eq!(
        artifact.source_snapshot_identity(),
        receipt.source_snapshot_identity()
    );
    assert!(
        artifact
            .json()
            .contains("\"authority\":\"declaration-only\"")
    );
}

#[test]
fn source_order_is_nonsemantic_and_content_mutation_moves_the_root() {
    let first_root = unique_root("order-forward");
    write_complete_fixture(&first_root, BEADS_BYTES, REGISTRY_BYTES);
    let forward = load_traceability_source_snapshot(
        &first_root,
        &complete_specs(),
        TraceabilityFilesystemLimits::default(),
    )
    .expect("forward source order");

    let mut reversed_specs = complete_specs();
    reversed_specs.reverse();
    let reverse = load_traceability_source_snapshot(
        &first_root,
        &reversed_specs,
        TraceabilityFilesystemLimits::default(),
    )
    .expect("reverse source order");
    assert_eq!(forward, reverse);

    let changed_root = unique_root("content-changed");
    write_complete_fixture(
        &changed_root,
        BEADS_BYTES,
        b"pub const SCHEMA: &str = \"traceability-v2\";\n",
    );
    let changed = load_traceability_source_snapshot(
        &changed_root,
        &complete_specs(),
        TraceabilityFilesystemLimits::default(),
    )
    .expect("changed registry source");
    assert_ne!(forward.snapshot().identity(), changed.snapshot().identity());
}

#[test]
fn loaded_owner_join_emits_a_declaration_only_ledger_and_receipt() {
    let loaded = loaded_owner_fixture("owner-complete", "closed", None, &[]);
    let audit = audit_traceability_owner_join(proof_obligations(), &loaded);
    assert!(
        audit.ok(),
        "unexpected owner diagnostics: {:?}",
        audit.diagnostics
    );

    let expected_total = proof_obligations()
        .iter()
        .map(|obligation| obligation.owner_beads.len())
        .sum();
    let expected_distinct = proof_obligations()
        .iter()
        .flat_map(|obligation| obligation.owner_beads.iter().copied())
        .collect::<BTreeSet<_>>()
        .len();
    assert_eq!(audit.total_owner_references, expected_total);
    assert_eq!(audit.distinct_owner_beads, expected_distinct);
    assert_eq!(audit.indexed_beads, expected_distinct);
    assert_eq!(audit.beads_source_count, 1);

    let generated = generate_traceability_ledger_from_loaded_sources(
        requirements(),
        proof_obligations(),
        &loaded,
    )
    .expect("every declared owner is present");
    assert_eq!(
        generated.owner_join().version(),
        TRACEABILITY_OWNER_JOIN_VERSION
    );
    assert_eq!(
        generated.owner_join().total_owner_references(),
        expected_total
    );
    assert_eq!(
        generated.owner_join().distinct_owner_beads(),
        expected_distinct
    );
    assert_eq!(generated.owner_join().indexed_beads(), expected_distinct);
    assert_eq!(generated.owner_join().beads_source_count(), 1);
    assert_eq!(
        generated.owner_join().source_snapshot_identity(),
        loaded.snapshot().identity()
    );
    assert_eq!(
        generated.ledger().source_snapshot_identity(),
        loaded.snapshot().identity()
    );
    assert!(
        generated
            .ledger()
            .json()
            .contains("\"authority\":\"declaration-only\"")
    );
}

#[test]
fn owner_join_refuses_one_missing_id_with_exact_obligation_context() {
    const MISSING: &str = "frankensim-ext-couple-port-schema-3feh";
    let loaded = loaded_owner_fixture("owner-missing", "open", Some(MISSING), &[]);
    let audit = audit_traceability_owner_join(proof_obligations(), &loaded);
    assert!(!audit.ok());
    assert_eq!(
        audit
            .diagnostics
            .iter()
            .filter(|diagnostic| diagnostic.owner_bead == MISSING)
            .count(),
        1
    );
    assert!(audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.proof_obligation_id == "PO-2"
            && diagnostic.owner_bead == MISSING
            && diagnostic.field == TraceabilityOwnerJoinField::OwnerBead
            && diagnostic.reason.contains("absent")
    }));

    let error = generate_traceability_ledger_from_loaded_sources(
        requirements(),
        proof_obligations(),
        &loaded,
    )
    .expect_err("missing owner refuses without a partial ledger");
    assert!(matches!(
        error,
        LoadedTraceabilityGenerationError::OwnerJoin(owner_audit)
            if owner_audit.diagnostics == audit.diagnostics
    ));
}

#[test]
fn tracker_status_is_bound_but_never_interpreted_as_owner_authority() {
    let open = loaded_owner_fixture("owners-open", "open", None, &[]);
    let closed = loaded_owner_fixture("owners-closed", "closed", None, &[]);
    let open_ledger = generate_traceability_ledger_from_loaded_sources(
        requirements(),
        proof_obligations(),
        &open,
    )
    .expect("open lexical owners are present");
    let closed_ledger = generate_traceability_ledger_from_loaded_sources(
        requirements(),
        proof_obligations(),
        &closed,
    )
    .expect("closed lexical owners are present");

    assert_ne!(open.snapshot().identity(), closed.snapshot().identity());
    assert_eq!(
        open_ledger.owner_join().total_owner_references(),
        closed_ledger.owner_join().total_owner_references()
    );
    assert_eq!(
        open_ledger.owner_join().distinct_owner_beads(),
        closed_ledger.owner_join().distinct_owner_beads()
    );
}

#[test]
fn declaration_validation_precedes_the_loaded_owner_join() {
    let loaded = loaded_owner_fixture("invalid-declaration", "open", None, &[]);
    let mut rows = requirements().to_vec();
    rows[0].owner_artifact = "";
    let error =
        generate_traceability_ledger_from_loaded_sources(&rows, proof_obligations(), &loaded)
            .expect_err("orphaned declaration refuses before source joining");
    assert!(matches!(
        error,
        LoadedTraceabilityGenerationError::Declaration(audit)
            if audit.diagnostics.iter().any(|diagnostic| {
                diagnostic.requirement_id == "B1"
                    && diagnostic.field == TraceabilityField::OwnerArtifact
            })
    ));
}

#[test]
fn beads_envelope_and_canonical_id_index_fail_closed() {
    let cases: [(&[u8], &str); 6] = [
        (
            b"{\"id\":\"duplicate\"}\n{\"id\":\"duplicate\"}\n",
            "repeats canonical id",
        ),
        (b"{\"title\":\"missing id\"}\n", "no canonical top-level"),
        (
            br#"{"id":"bead-a","bad":"\q"}
"#,
            "invalid JSON string escape",
        ),
        (
            b"{\"id\":\"bead-a\",\"nested\":[1,2}\n",
            "delimiters are mismatched",
        ),
        (
            b"{\"id\":\"bead-a\",\"id\":\"bead-b\"}\n",
            "repeats the top-level id key",
        ),
        (b"{\"\\u0069d\":\"bead-a\"}\n", "no canonical top-level"),
    ];

    for (beads, expected_reason) in cases {
        let audit = audit_with_beads(beads);
        assert!(
            audit.diagnostics.iter().any(|diagnostic| {
                diagnostic.field == TraceabilityFilesystemField::BeadsJsonl
                    && diagnostic.reason.contains(expected_reason)
            }),
            "missing {expected_reason:?} in {:?}",
            audit.diagnostics
        );
    }
}

#[test]
fn contract_and_registry_text_sources_fail_closed() {
    let cases: [(&str, &[u8], &[u8], &str); 3] = [
        (
            "empty-contract",
            b"",
            REGISTRY_BYTES,
            "source file is empty",
        ),
        (
            "nul-contract",
            b"contract\0content",
            REGISTRY_BYTES,
            "contains a NUL byte",
        ),
        (
            "invalid-registry-utf8",
            CONTRACT_BYTES,
            &[0xff],
            "text source is not UTF-8",
        ),
    ];

    for (label, contract, registry, expected_reason) in cases {
        let audit = audit_with_text(label, contract, registry);
        assert!(
            audit.diagnostics.iter().any(|diagnostic| {
                diagnostic.field == TraceabilityFilesystemField::Content
                    && diagnostic.reason.contains(expected_reason)
            }),
            "missing {expected_reason:?} in {:?}",
            audit.diagnostics
        );
    }
}

#[test]
fn configured_record_line_file_and_total_caps_are_enforced() {
    let root = unique_root("caps");
    write_complete_fixture(&root, BEADS_BYTES, REGISTRY_BYTES);

    let mut record_limits = TraceabilityFilesystemLimits::default();
    record_limits.max_beads_records = 1;
    let record_audit = load_traceability_source_snapshot(&root, &complete_specs(), record_limits)
        .expect_err("two Beads records exceed the configured cap");
    assert!(record_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::BeadsJsonl
            && diagnostic.reason.contains("1-record maximum")
    }));

    let mut line_limits = TraceabilityFilesystemLimits::default();
    line_limits.max_beads_line_bytes = 8;
    let line_audit = load_traceability_source_snapshot(&root, &complete_specs(), line_limits)
        .expect_err("long Beads line exceeds the configured cap");
    assert!(line_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::BeadsJsonl
            && diagnostic.reason.contains("configured maximum is 8")
    }));

    let mut nesting_limits = TraceabilityFilesystemLimits::default();
    nesting_limits.max_beads_json_nesting = 1;
    let nesting_audit = load_traceability_source_snapshot(&root, &complete_specs(), nesting_limits)
        .expect_err("nested Beads object exceeds the configured cap");
    assert!(nesting_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::BeadsJsonl
            && diagnostic
                .reason
                .contains("nesting exceeds configured maximum 1")
    }));

    let mut file_limits = TraceabilityFilesystemLimits::default();
    file_limits.max_source_bytes = 16;
    file_limits.max_total_bytes = 64;
    let file_audit = load_traceability_source_snapshot(&root, &complete_specs(), file_limits)
        .expect_err("source file exceeds the configured cap");
    assert!(file_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::Content
            && diagnostic.reason.contains("configured maximum is 16")
    }));

    let largest = BEADS_BYTES
        .len()
        .max(CONTRACT_BYTES.len())
        .max(REGISTRY_BYTES.len());
    let total = BEADS_BYTES.len() + CONTRACT_BYTES.len() + REGISTRY_BYTES.len();
    let mut total_limits = TraceabilityFilesystemLimits::default();
    total_limits.max_source_bytes = largest;
    total_limits.max_total_bytes = total - 1;
    let total_audit = load_traceability_source_snapshot(&root, &complete_specs(), total_limits)
        .expect_err("aggregate source bytes exceed the configured cap");
    assert!(total_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::Content
            && diagnostic.reason.contains("total-byte maximum")
    }));
}

#[test]
fn paths_metadata_coverage_and_limit_configuration_refuse_ambiguity() {
    let root = unique_root("admission");
    write_complete_fixture(&root, BEADS_BYTES, REGISTRY_BYTES);

    let mut invalid_path_specs = complete_specs();
    invalid_path_specs[0].relative_path = Path::new("../outside-registry.rs");
    let path_audit = load_traceability_source_snapshot(
        &root,
        &invalid_path_specs,
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("parent traversal refuses before reads");
    assert!(path_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::RelativePath
            && diagnostic.reason.contains("strict root-relative")
    }));

    let mut current_dir_specs = complete_specs();
    current_dir_specs[0].relative_path = Path::new("./registry/traceability.rs");
    let current_dir_audit = load_traceability_source_snapshot(
        &root,
        &current_dir_specs,
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("current-directory spelling is not canonical");
    assert!(current_dir_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::RelativePath
            && (diagnostic.reason.contains("strict root-relative")
                || diagnostic.reason.contains("without '.'"))
    }));

    let mut missing_specs = complete_specs();
    missing_specs[0].relative_path = Path::new("registry/missing.rs");
    let missing_audit = load_traceability_source_snapshot(
        &root,
        &missing_specs,
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("missing regular file refuses");
    assert!(missing_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::Metadata
            && diagnostic.reason.contains("cannot inspect source file")
    }));

    let incomplete = [complete_specs()[0], complete_specs()[1]];
    let coverage_audit = load_traceability_source_snapshot(
        &root,
        &incomplete,
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("all three source classes are mandatory");
    assert!(coverage_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::Snapshot
            && diagnostic.reason.contains("no contract artifact")
    }));

    let invalid_limits = TraceabilityFilesystemLimits {
        max_sources: 0,
        ..TraceabilityFilesystemLimits::default()
    };
    let limit_audit = load_traceability_source_snapshot(&root, &complete_specs(), invalid_limits)
        .expect_err("zero source cap is invalid configuration");
    assert!(limit_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::Limits
            && diagnostic.reason.contains("max_sources must be in")
    }));
}

#[test]
fn non_regular_sources_and_duplicate_locators_refuse() {
    let directory_root = unique_root("non-regular");
    write_source(&directory_root, BEADS_PATH, BEADS_BYTES);
    write_source(&directory_root, CONTRACT_PATH, CONTRACT_BYTES);
    std::fs::create_dir_all(directory_root.join(REGISTRY_PATH)).expect("directory-shaped source");
    let directory_audit = load_traceability_source_snapshot(
        &directory_root,
        &complete_specs(),
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("directory cannot stand in for source bytes");
    assert!(directory_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::Metadata
            && diagnostic.reason.contains("not a regular file")
    }));

    let duplicate_root = unique_root("duplicate-locator");
    write_complete_fixture(&duplicate_root, BEADS_BYTES, REGISTRY_BYTES);
    let duplicate_specs = [
        complete_specs()[0],
        complete_specs()[1],
        TraceabilityFileSpec {
            kind: TraceabilitySourceKind::Contract,
            relative_path: Path::new(BEADS_PATH),
        },
    ];
    let duplicate_audit = load_traceability_source_snapshot(
        &duplicate_root,
        &duplicate_specs,
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("one locator cannot be relabeled");
    assert!(duplicate_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::RelativePath
            && diagnostic.reason.contains("duplicate source locator")
    }));

    let duplicate_id_root = unique_root("duplicate-beads-id");
    write_source(
        &duplicate_id_root,
        ".beads/first.jsonl",
        b"{\"id\":\"same-bead\"}\n",
    );
    write_source(
        &duplicate_id_root,
        ".beads/second.jsonl",
        b"{\"id\":\"same-bead\"}\n",
    );
    write_source(&duplicate_id_root, CONTRACT_PATH, CONTRACT_BYTES);
    write_source(&duplicate_id_root, REGISTRY_PATH, REGISTRY_BYTES);
    let duplicate_id_specs = [
        TraceabilityFileSpec {
            kind: TraceabilitySourceKind::Beads,
            relative_path: Path::new(".beads/first.jsonl"),
        },
        TraceabilityFileSpec {
            kind: TraceabilitySourceKind::Beads,
            relative_path: Path::new(".beads/second.jsonl"),
        },
        complete_specs()[0],
        complete_specs()[2],
    ];
    let duplicate_id_audit = load_traceability_source_snapshot(
        &duplicate_id_root,
        &duplicate_id_specs,
        TraceabilityFilesystemLimits::default(),
    )
    .expect_err("one Bead id cannot appear in two bound sources");
    assert!(duplicate_id_audit.diagnostics.iter().any(|diagnostic| {
        diagnostic.field == TraceabilityFilesystemField::BeadsJsonl
            && diagnostic.reason.contains("more than one Beads source")
    }));
}
